module CocoaMLEDiagnostics

using FermentModels
using ForwardDiff
using LinearAlgebra
using TOML
using SHA
include("prepare_cocoa.jl")
using .CocoaSourceData: cocoa_trial
include("cocoa_fit_output.jl")

"""
Inspect a parameter point without accepting it as an optimum. Gradients are with
respect to log parameters. The projected gradient is z-clamp(z-gradient,lo,hi),
so outward gradients at a bound do not count as feasible descent directions.
Finite differences evaluate the unbounded likelihood, including near search
bounds. Expected-information sensitivities include the shared noise parameter.
"""
function point_diagnostics(model, data, parameters; bound_rtol=1e-5,
                           rank_rtol=sqrt(eps(Float64)), fd_step=1e-4)
    0 < bound_rtol < 0.5 || throw(ArgumentError("invalid bound tolerance"))
    0 < rank_rtol < 1 || throw(ArgumentError("invalid rank tolerance"))
    isfinite(fd_step) && fd_step > 0 || throw(ArgumentError("invalid finite-difference step"))
    data.sigma isa Symbol || throw(ArgumentError("requires a fitted shared noise scale"))
    names = keys(model.initial)
    Set(keys(parameters)) == Set(names) || throw(ArgumentError("parameter names differ"))
    p = NamedTuple{names}(Tuple(parameters[name] for name in names))
    lo, hi = log.(collect(values(model.lower))), log.(collect(values(model.upper)))
    z = log.(collect(values(p)))
    all(isfinite,z) && all(lo .<= z .<= hi) || throw(ArgumentError("point outside bounds"))
    decode = x -> NamedTuple{names}(Tuple(exp.(x)))
    objective = x -> gaussian_nll(model,data,decode(x))
    nll = objective(z)
    gradient = ForwardDiff.gradient(objective,z)
    all(isfinite,gradient) || error("Nonfinite likelihood gradient")
    projected = z - clamp.(z-gradient,lo,hi)
    lower_distance, upper_distance = (z-lo)./(hi-lo), (hi-z)./(hi-lo)
    active_lower, active_upper = lower_distance .<= bound_rtol, upper_distance .<= bound_rtol
    finite_differences = map((fd_step,fd_step/2)) do step
        derivative = map(eachindex(z)) do i
            a, b = copy(z), copy(z)
            a[i] += step
            b[i] -= step
            (objective(a)-objective(b))/(2step)
        end
        (; step, gradient=derivative,
           max_scaled_error=maximum(abs.(derivative-gradient)./max.(1.,abs.(gradient))))
    end
    sensitivity = local_sensitivities(model,data,p;rtol=rank_rtol)
    noise_index = findfirst(==(data.sigma),names)
    noise_row = zeros(1,length(names))
    noise_row[noise_index] = sqrt(2length(sensitivity.observed_indices))
    # A mean-only Jacobian has a zero sigma column even when sigma is estimable.
    information_factor = vcat(sensitivity.jacobian,noise_row)
    decomposition = svd(information_factor;full=size(information_factor,1)<length(names))
    cutoff = rank_rtol * maximum(decomposition.S)
    numerical_rank = count(>(cutoff),decomposition.S)
    prediction = model.predict(p,data)
    residual = data.values[sensitivity.observed_indices] - prediction[sensitivity.observed_indices]
    conditional_sigma = sqrt(sum(abs2,residual)/length(residual))
    return (; names, parameters=p, nll, gradient, projected,
        projected_gradient_inf=norm(projected,Inf), lower_distance, upper_distance,
        active_lower, active_upper, bound_rtol, rank_rtol, finite_differences,
        information_factor, singular_values=decomposition.S, information_rank=numerical_rank,
        information_condition=numerical_rank < length(names) ? Inf :
            maximum(decomposition.S)/minimum(decomposition.S),
        weakest_direction=decomposition.V[:,end], conditional_sigma,
        observation_count=length(residual))
end

file_hash(path) = bytes2hex(open(sha256,path))

"""Write a new, provenance-linked diagnostic directory for one saved MLE start."""
function diagnose_cocoa_mle(fit_directory; start=1, output)
    ispath(output) && throw(ArgumentError("diagnostic output already exists"))
    settings_path = joinpath(fit_directory,"settings.toml")
    checkpoint = joinpath(fit_directory,"mle_start_$(lpad(start,3,'0')).toml")
    settings = TOML.parsefile(settings_path)
    settings["method"] == "mle" || throw(ArgumentError("requires an MLE run"))
    saved = TOML.parsefile(checkpoint)
    root = normpath(joinpath(@__DIR__,".."))
    source = cocoa_trial(Symbol(settings["trial"]))
    raw = joinpath(root,source.source_file)
    file_hash(raw) == settings["raw_sha256"] || error("Raw input differs from fitted data")
    data = cocoa_observations(source.initial,source.times,source.observations)
    template = load_params(joinpath(root,"params","cocoa_$(settings["trial"])_scaled.toml")).values
    parameters = cocoa_fit_initial(template,checkpoint)
    lower = NamedTuple{keys(parameters)}(ntuple(_->settings["mle_lower"],length(parameters)))
    upper = NamedTuple{keys(parameters)}(ntuple(_->settings["mle_upper"],length(parameters)))
    # A valid saved point can sit exactly at a bound; only the constructor's
    # unused starting point must remain strictly interior.
    initial = NamedTuple{keys(parameters)}(Tuple(sqrt(lower[n]*upper[n]) for n in keys(parameters)))
    solver = settings["solver"] == "Tsit5" ? :tsit5 :
        settings["solver"] == "Rodas5P" ? :rodas5p : error("Unknown saved solver")
    model = cocoa_fit_model(;initial,lower,upper,solver,abstol=settings["abstol"],
        reltol=settings["reltol"],maxiters=settings["solver_maxiters"])
    diagnostic = point_diagnostics(model,data,parameters)
    reference_prediction = model.predict(parameters,data)
    # Compare solvers at the saved tolerances and the same solver at 10x tighter
    # tolerances. Failed verification solves remain failures in the report.
    solver_checks = map([(solver,0.1),(solver == :tsit5 ? :rodas5p : :tsit5,1.0)]) do (algorithm,factor)
        check_model = cocoa_fit_model(;initial,lower,upper,solver=algorithm,
            abstol=factor*settings["abstol"],reltol=factor*settings["reltol"],
            maxiters=settings["solver_maxiters"])
        prediction = check_model.predict(parameters,data)
        Dict("solver"=>String(algorithm),"tolerance_factor"=>factor,
            "solve_succeeded"=>!isnothing(prediction),
            "nll"=>isnothing(prediction) ? Inf : gaussian_nll(check_model,data,parameters),
            "maximum_normalized_prediction_difference"=>isnothing(prediction) ? Inf :
                maximum(abs.(prediction-reference_prediction)))
    end
    report = Dict(
        "trial"=>settings["trial"], "checkpoint_converged"=>saved["converged"],
        "checkpoint_retcode"=>saved["retcode"],
        "checkpoint_time_limit_reached"=>saved["time_limit_reached"],
        "checkpoint_objective"=>saved["objective"], "reevaluated_nll"=>diagnostic.nll,
        "diagnostics_accept_fit"=>false,
        "interpretation"=>"Local point diagnostics do not establish a local/global optimum or structural identifiability.",
        "projected_gradient_definition"=>"z - clamp(z - gradient(NLL), log(lower), log(upper)); z=log(parameters)",
        "projected_gradient_inf"=>diagnostic.projected_gradient_inf,
        "bound_relative_tolerance"=>diagnostic.bound_rtol,
        "active_lower"=>[String(diagnostic.names[i]) for i in findall(diagnostic.active_lower)],
        "active_upper"=>[String(diagnostic.names[i]) for i in findall(diagnostic.active_upper)],
        "information_rank"=>diagnostic.information_rank,
        "information_parameter_count"=>length(diagnostic.names),
        "information_rank_relative_tolerance"=>diagnostic.rank_rtol,
        "information_condition_number"=>diagnostic.information_condition,
        "information_singular_values"=>diagnostic.singular_values,
        "information_definition"=>"SVD of sigma-standardized mean Jacobian in log parameters, augmented by sqrt(2n) in log(sigma).",
        "conditional_sigma_mle"=>diagnostic.conditional_sigma,
        "saved_sigma"=>parameters[data.sigma],
        "observation_count"=>diagnostic.observation_count,
        "finite_difference_checks"=>[Dict("step"=>c.step,"max_scaled_error"=>c.max_scaled_error)
            for c in diagnostic.finite_differences],
        "solver_checks"=>solver_checks,
        "source"=>"doi:10.1098/rsos.180964; doi:10.5061/dryad.321d33v",
        "parameter_coordinates"=>"normalized by original trial state maxima",
        "time_unit"=>"h", "concentration_source_unit"=>"mg/g pulp",
        "julia_version"=>string(VERSION),
        "checkpoint"=>replace(relpath(abspath(checkpoint),root),'\\'=>'/'),
        "checkpoint_sha256"=>file_hash(checkpoint),
        "fit_settings_sha256"=>file_hash(settings_path),
        "raw_file"=>replace(source.source_file,'\\'=>'/'), "raw_sha256"=>file_hash(raw),
        "implementation_sha256"=>Dict(path=>file_hash(joinpath(root,path)) for path in
            ("scripts/cocoa_mle_diagnostics.jl","scripts/cocoa_fit_output.jl","scripts/prepare_cocoa.jl",
             "src/cocoa/baseline.jl","src/cocoa/inference.jl","src/core/estimation.jl",
             "src/core/identifiability.jl")))
    # Validate serialization before creating the output directory.
    report_text = sprint(io -> TOML.print(io,report;sorted=true))
    mkpath(output)
    open(joinpath(output,"report.toml"),"w") do io
        write(io,report_text)
    end
    open(joinpath(output,"parameters.csv"),"w") do io
        println(io,"parameter,value,gradient_log,projected_gradient,lower_log_range_fraction,upper_log_range_fraction,active_lower,active_upper,fd_gradient,fd_half_step_gradient,weakest_information_direction")
        for (i,name) in enumerate(diagnostic.names)
            println(io,join((name,parameters[name],diagnostic.gradient[i],diagnostic.projected[i],
                diagnostic.lower_distance[i],diagnostic.upper_distance[i],
                diagnostic.active_lower[i],diagnostic.active_upper[i],
                diagnostic.finite_differences[1].gradient[i],diagnostic.finite_differences[2].gradient[i],
                diagnostic.weakest_direction[i]),','))
        end
    end
    println(settings["trial"]," point diagnostics: NLL=",diagnostic.nll,
        "; projected gradient=",diagnostic.projected_gradient_inf,
        "; information rank=",diagnostic.information_rank,"/",length(diagnostic.names))
    flush(stdout)
    return report
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (2,3) || error("Usage: cocoa_mle_diagnostics.jl FIT_DIRECTORY NEW_OUTPUT [START]")
    CocoaMLEDiagnostics.diagnose_cocoa_mle(ARGS[1];output=ARGS[2],
        start=length(ARGS)==3 ? parse(Int,ARGS[3]) : 1)
end
