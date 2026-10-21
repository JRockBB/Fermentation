module CocoaBenchmarkDiagnostics

using FermentModels
using TOML
include("fit_cocoa_benchmark.jl")
include("cocoa_mle_diagnostics.jl")
using .CocoaBenchmarkFits
using .CocoaMLEDiagnostics
const B = CocoaBenchmarkFits.CocoaBenchmarkData

"""
Inspect saved benchmark candidates at their fitted tolerances and at 100x
tighter tolerances. Output is new and checkpointed per trial/tolerance pair.
No diagnostic in this report accepts a candidate as an optimum.
"""
function diagnose_benchmark(fit_directory; output, trials=nothing,
        tolerance_factors=(1.0,0.01), fd_step=1e-4)
    ispath(output) && throw(ArgumentError("diagnostic output already exists"))
    !isempty(tolerance_factors) && all(f->f isa Real && isfinite(f) && 0 < f <= 1,
        tolerance_factors) && length(unique(tolerance_factors)) == length(tolerance_factors) ||
        throw(ArgumentError("tolerance factors must be distinct and in (0,1]"))
    isfinite(fd_step) && fd_step > 0 || throw(ArgumentError("fd_step must be positive"))
    settings = TOML.parsefile(joinpath(fit_directory,"settings.toml"))
    available = B.benchmark_trials(;scope=Symbol(settings["scope"]))
    # New runs explicitly list their subset; the original pilot used its whole scope.
    saved_ids = Symbol.(get(settings,"trials",String.([t.trial for t in available])))
    available = CocoaBenchmarkFits.select_trials(available,saved_ids)
    selected = CocoaBenchmarkFits.select_trials(available,trials)
    candidates = Dict(t.trial=>CocoaBenchmarkFits.saved_candidate(fit_directory,t) for t in selected)
    names = (COCOA_PARAMETERS..., :sigma)
    lower = NamedTuple{names}(ntuple(_->settings["mle_bounds"][1],25))
    upper = NamedTuple{names}(ntuple(_->settings["mle_bounds"][2],25))
    initial = NamedTuple{names}(Tuple(sqrt(lower[n]*upper[n]) for n in names))
    solver = Symbol(settings["solver"])
    implementation = ("scripts/diagnose_cocoa_benchmark.jl","scripts/cocoa_mle_diagnostics.jl",
        "scripts/fit_cocoa_benchmark.jl","src/core/identifiability.jl",
        CocoaBenchmarkFits.LIKELIHOOD_IMPLEMENTATION...)
    metadata = Dict(
        "fit_directory"=>replace(relpath(abspath(fit_directory),B.ROOT),'\\'=>'/'),
        "fit_settings_sha256"=>B.digest(joinpath(fit_directory,"settings.toml")),
        "trials"=>String.([t.trial for t in selected]),"tolerance_factors"=>collect(tolerance_factors),
        "finite_difference_step"=>fd_step,"solver"=>string(solver),
        "diagnostics_accept_fit"=>false,
        "interpretation"=>"Local numerical derivatives and information rank do not establish an optimum or structural identifiability.",
        "projected_gradient_definition"=>"z - clamp(z - gradient(NLL), log(lower), log(upper)); z=log(parameters)",
        "information_definition"=>"SVD of sigma-standardized mean Jacobian in log parameters, augmented by sqrt(2n) in log(sigma).",
        "parameter_coordinates"=>"normalized by trial state maxima",
        "julia_version"=>string(VERSION),
        "implementation_sha256"=>Dict(p=>B.digest(joinpath(B.ROOT,p)) for p in implementation))
    mkpath(output)
    CocoaBenchmarkFits.write_toml(joinpath(output,"settings.toml"),metadata)
    rows = Tuple[]
    columns = ("trial","tolerance_factor","solve_succeeded","nll","projected_gradient_inf",
        "fd_scaled_error","fd_half_step_scaled_error","information_rank","information_condition",
        "active_bound_count","saved_sigma","conditional_sigma","observations")
    for trial in selected
        saved = candidates[trial.trial]
        directory = joinpath(output,string(trial.trial))
        mkpath(directory)
        CocoaBenchmarkFits.write_toml(joinpath(directory,"provenance.toml"),saved.provenance)
        for (index,factor) in enumerate(tolerance_factors)
            model = cocoa_fit_model(;initial,lower,upper,solver,
                abstol=factor*settings["abstol"],reltol=factor*settings["reltol"],
                maxiters=settings["solver_maxiters"])
            prediction = model.predict(saved.parameters,trial.data)
            report = Dict{String,Any}("trial"=>string(trial.trial),"tolerance_factor"=>factor,
                "abstol"=>factor*settings["abstol"],"reltol"=>factor*settings["reltol"],
                "solve_succeeded"=>!isnothing(prediction),"diagnostics_accept_fit"=>false)
            if !isnothing(prediction)
                d = CocoaMLEDiagnostics.point_diagnostics(model,trial.data,saved.parameters;fd_step)
                merge!(report,Dict("nll"=>d.nll,"projected_gradient_inf"=>d.projected_gradient_inf,
                    "observation_count"=>d.observation_count,
                    "finite_difference_step"=>fd_step,
                    "finite_difference_scaled_errors"=>[c.max_scaled_error for c in d.finite_differences],
                    "information_rank"=>d.information_rank,"information_parameter_count"=>length(names),
                    "information_condition_number"=>d.information_condition,
                    "information_rank_relative_tolerance"=>d.rank_rtol,
                    "information_singular_values"=>d.singular_values,
                    "bound_relative_tolerance"=>d.bound_rtol,
                    "active_lower"=>String.(collect(names)[findall(d.active_lower)]),
                    "active_upper"=>String.(collect(names)[findall(d.active_upper)]),
                    "conditional_sigma_mle"=>d.conditional_sigma,"saved_sigma"=>saved.parameters.sigma))
                B.write_csv(joinpath(directory,"parameters_$(index).csv"),
                    ("parameter","value","gradient_log","projected_gradient",
                     "fd_gradient","fd_half_step_gradient","lower_log_range_fraction",
                     "upper_log_range_fraction","active_lower","active_upper","weakest_information_direction"),
                    [(name,saved.parameters[name],d.gradient[i],d.projected[i],
                        d.finite_differences[1].gradient[i],d.finite_differences[2].gradient[i],
                        d.lower_distance[i],d.upper_distance[i],d.active_lower[i],d.active_upper[i],
                        d.weakest_direction[i]) for (i,name) in enumerate(names)])
                push!(rows,(trial.trial,factor,true,d.nll,d.projected_gradient_inf,
                    (c.max_scaled_error for c in d.finite_differences)...,
                    d.information_rank,d.information_condition,count(d.active_lower .| d.active_upper),
                    saved.parameters.sigma,d.conditional_sigma,d.observation_count))
            else
                push!(rows,(trial.trial,factor,false,ntuple(_->missing,10)...))
            end
            CocoaBenchmarkFits.write_toml(joinpath(directory,"report_$(index).toml"),report)
            B.write_csv(joinpath(output,"comparison.csv"),columns,rows)
            println("Diagnosed ",trial.trial," at tolerance factor ",factor,
                "; projected gradient=",get(report,"projected_gradient_inf","unavailable"))
            flush(stdout)
        end
    end
    return rows
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (2,3) || error("Usage: diagnose_cocoa_benchmark.jl FIT_DIRECTORY NEW_OUTPUT [TRIAL,TRIAL,...]")
    CocoaBenchmarkDiagnostics.diagnose_benchmark(ARGS[1];output=ARGS[2],
        trials=length(ARGS)==3 ? Symbol.(split(ARGS[3],',')) : nothing)
end
