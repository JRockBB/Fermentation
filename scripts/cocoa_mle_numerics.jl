include("cocoa_mle_diagnostics.jl")
using .CocoaMLEDiagnostics
using FermentModels, ForwardDiff, TOML, SHA

"""
Check a saved MLE point across two integrators, two tolerance pairs, and four
finite-difference step sizes. This changes numerical settings only, not the
likelihood, bounds, or parameter point. Failed solves remain explicit.
"""
function check_cocoa_mle_numerics(fit_directory; start=1, output)
    ispath(output) && throw(ArgumentError("output directory already exists"))
    root = normpath(joinpath(@__DIR__,".."))
    settings_path = joinpath(fit_directory,"settings.toml")
    settings = TOML.parsefile(settings_path)
    settings["method"] == "mle" || throw(ArgumentError("requires an MLE fit"))
    checkpoint = joinpath(fit_directory,"mle_start_$(lpad(start,3,'0')).toml")
    saved = TOML.parsefile(checkpoint)
    # Earlier fit_cocoa.jl versions fixed gtol=1e-6 without saving that field.
    gradient_tolerance = get(settings,"mle_gradient_tolerance",1e-6)
    source = CocoaMLEDiagnostics.cocoa_trial(Symbol(settings["trial"]))
    raw_path = joinpath(root,source.source_file)
    CocoaMLEDiagnostics.file_hash(raw_path) == settings["raw_sha256"] || error("Raw input changed")
    data = cocoa_observations(source.initial,source.times,source.observations)
    names = (COCOA_PARAMETERS...,:sigma)
    template = NamedTuple{names}(ntuple(_->0.5,length(names)))
    p = CocoaMLEDiagnostics.cocoa_fit_initial(template,checkpoint)
    lower = NamedTuple{names}(ntuple(_->settings["mle_lower"],length(names)))
    upper = NamedTuple{names}(ntuple(_->settings["mle_upper"],length(names)))
    initial = NamedTuple{names}(Tuple(sqrt(lower[n]*upper[n]) for n in names))
    z = log.(collect(values(p)))
    checks = []
    for (algorithm,atol,rtol) in ((:rodas5p,1e-9,1e-8),(:rodas5p,1e-11,1e-10),
                                  (:tsit5,1e-9,1e-8),(:tsit5,1e-11,1e-10))
        model = cocoa_fit_model(;initial,lower,upper,solver=algorithm,
            abstol=atol,reltol=rtol,maxiters=settings["solver_maxiters"])
        loss = x -> gaussian_nll(model,data,NamedTuple{names}(Tuple(exp.(x))))
        check = Dict{String,Any}("solver"=>String(algorithm),"abstol"=>atol,"reltol"=>rtol,
            "solver_maxiters"=>settings["solver_maxiters"],"completed"=>false,
            "stationarity_passed"=>false,"gradient_tolerance"=>gradient_tolerance)
        checks_fd = []
        try
            check["nll"] = loss(z)
            gradient = ForwardDiff.gradient(loss,z)
            check["gradient"] = gradient
            check["projected_gradient_inf"] =
                maximum(abs.(z-clamp.(z-gradient,log.(collect(values(lower))),log.(collect(values(upper))))))
            check["stationarity_passed"] = check["projected_gradient_inf"] <= gradient_tolerance
            for step in (1e-3,1e-4,1e-5,1e-6)
                fd = map(eachindex(z)) do i
                    a,b = copy(z),copy(z)
                    a[i] += step
                    b[i] -= step
                    (loss(a)-loss(b))/(2step)
                end
                error = maximum(abs.(fd-gradient)./max.(1.,abs.(gradient)))
                push!(checks_fd,Dict("step"=>step,"max_scaled_error"=>error,"gradient"=>fd))
                println(settings["trial"]," ",algorithm," atol=",atol," step=",step," error=",error)
                flush(stdout)
            end
            check["completed"] = true
        catch err
            err isa FermentModels.CoreModels._PredictionFailure || rethrow()
            check["error"] = sprint(showerror,err)
        end
        check["finite_differences"] = checks_fd
        push!(checks,check)
    end
    report = Dict("trial"=>settings["trial"],"parameter_names"=>String.(collect(names)),
        "parameters"=>Dict(String(n)=>p[n] for n in names),
        "parameter_coordinates"=>"normalized by original trial state maxima; derivatives in log parameters",
        "checks"=>checks,"diagnostics_accept_fit"=>false,
        "checkpoint_converged"=>saved["converged"],"checkpoint_retcode"=>saved["retcode"],
        "gradient_tolerance_source"=>haskey(settings,"mle_gradient_tolerance") ?
            "fit settings" : "legacy scripts/fit_cocoa.jl fixed gtol=1e-6",
        "finite_difference_error_definition"=>"max(abs(FD-AD)/max(1,abs(AD))) over parameters",
        "checkpoint"=>replace(relpath(abspath(checkpoint),root),'\\'=>'/'),
        "checkpoint_sha256"=>CocoaMLEDiagnostics.file_hash(checkpoint),
        "fit_settings_sha256"=>CocoaMLEDiagnostics.file_hash(settings_path),
        "raw_file"=>replace(source.source_file,'\\'=>'/'),
        "raw_sha256"=>CocoaMLEDiagnostics.file_hash(raw_path),
        "julia_version"=>string(VERSION),
        "source"=>"doi:10.1098/rsos.180964; doi:10.5061/dryad.321d33v",
        "implementation_sha256"=>Dict(path=>CocoaMLEDiagnostics.file_hash(joinpath(root,path))
            for path in ("scripts/cocoa_mle_numerics.jl","scripts/cocoa_mle_diagnostics.jl",
                "scripts/cocoa_fit_output.jl","scripts/prepare_cocoa.jl",
                "src/cocoa/baseline.jl","src/cocoa/inference.jl","src/core/estimation.jl")))
    content = sprint(io->TOML.print(io,report;sorted=true))
    mkpath(output)
    write(joinpath(output,"report.toml"),content)
    return report
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (2,3) || error("Usage: cocoa_mle_numerics.jl FIT_DIRECTORY NEW_OUTPUT [START]")
    check_cocoa_mle_numerics(ARGS[1];output=ARGS[2],start=length(ARGS)==3 ? parse(Int,ARGS[3]) : 1)
end
