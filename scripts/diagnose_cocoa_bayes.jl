module CocoaBayesDiagnostics

using FermentModels
using ForwardDiff
using Distributions: logpdf
using TOML
include("fit_cocoa_benchmark_bayes.jl")
const C = CocoaBenchmarkBayes
const B = C.B
const NAMES = (COCOA_PARAMETERS...,:sigma)

"""Read complete checkpoints only when the data, likelihood, and priors still match."""
function saved_draws(directory,trial;root=B.ROOT)
    settings_path = joinpath(directory,"settings.toml")
    settings = TOML.parsefile(settings_path)
    get(settings,"method",nothing) == "bayes" || throw(ArgumentError("requires a Bayesian run"))
    string(trial.trial) in settings["trials"] || throw(ArgumentError("trial absent from saved run"))
    for file in (C.F.LIKELIHOOD_IMPLEMENTATION...,"src/core/bayes.jl")
        get(get(settings,"implementation_sha256",Dict()),file,nothing) == B.digest(joinpath(root,file)) ||
            throw(ArgumentError("saved scientific implementation differs: $file"))
    end
    path = joinpath(directory,string(trial.trial))
    input = TOML.parsefile(joinpath(path,"input.toml"))
    expected = Dict("trial"=>string(trial.trial),"source_file"=>trial.source_file,
        "source_sha256"=>trial.source_sha256,"source_branch"=>trial.branch,
        "hidden_lactate"=>trial.hidden_lactate,"initial_normalized"=>trial.data.initial,
        "times_h"=>trial.data.times,"state_order"=>String.(collect(COCOA_STATES)),
        "seed"=>settings["trial_seeds"][string(trial.trial)],
        "solver_maxiters"=>settings["trial_solver_maxiters"][string(trial.trial)])
    all(get(input,k,nothing)==v for (k,v) in expected) ||
        throw(ArgumentError("saved input differs from current trial"))
    for (key,default) in (("solver_maxiters_scope","solve"),("solver_max_step_attempts",500))
        get(input,key,default) == get(settings,key,default) ||
            throw(ArgumentError("saved input solver budget differs: $key"))
    end
    get(settings,"solver_maxiters_scope","solve") in ("solve","output_interval") ||
        throw(ArgumentError("invalid saved solver budget scope"))
    priors = C.benchmark_priors(;hidden_lactate=trial.hidden_lactate,root)
    prior_file = TOML.parsefile(joinpath(path,"priors.toml"))
    Set(keys(prior_file)) == Set(String.(NAMES)) || throw(ArgumentError("prior names differ"))
    all(all(get(prior_file[string(n)],string(k),nothing)==v for
        (k,v) in pairs(priors[n].provenance)) for n in NAMES) ||
        throw(ArgumentError("saved prior provenance differs"))
    report = TOML.parsefile(joinpath(path,"report.toml"))
    get(report,"trial",nothing) == string(trial.trial) || throw(ArgumentError("report trial differs"))
    get(report,"status",nothing) in ("posterior_diagnostics_passed","posterior_diagnostics_failed") ||
        throw(ArgumentError("requires a completed sampling report"))
    samples,chains = settings["samples_per_chain"],settings["chains"]
    all(v->v isa Integer && !(v isa Bool) && v>0,(samples,chains)) ||
        throw(ArgumentError("invalid saved sample dimensions"))
    rows = NamedTuple[]
    files = ["input.toml","priors.toml","report.toml"]
    for chain in 1:chains
        file = "chain_$(lpad(chain,3,'0')).csv"
        lines = readlines(joinpath(path,file))
        expected_header = ["chain","draw",String.(collect(NAMES))...,"divergent","tree_depth"]
        !isempty(lines) && split(first(lines),',') == expected_header ||
            throw(ArgumentError("checkpoint header differs"))
        length(lines) == samples+1 || throw(ArgumentError("incomplete checkpoint"))
        for (draw,line) in enumerate(lines[2:end])
            fields = split(line,',';keepempty=true)
            length(fields) == length(expected_header) || throw(ArgumentError("incomplete checkpoint row"))
            parse(Int,fields[1]) == chain && parse(Int,fields[2]) == draw ||
                throw(ArgumentError("checkpoint chain/draw order differs"))
            values = parse.(Float64,fields[3:end-2])
            all(v->isfinite(v) && v>0,values) || throw(ArgumentError("invalid parameter draw"))
            divergent,depth = parse(Bool,fields[end-1]),parse(Int,fields[end])
            0<=depth<=settings["max_depth"] || throw(ArgumentError("invalid tree depth"))
            push!(rows,(;chain,draw,parameters=NamedTuple{NAMES}(Tuple(values)),
                divergent,tree_depth=depth))
        end
        push!(files,file)
    end
    get(report,"retained_draws",nothing) == length(rows) ||
        throw(ArgumentError("report/checkpoint draw count differs"))
    get(report,"observations",nothing) == count(!ismissing,trial.data.values) ||
        throw(ArgumentError("report observation count differs"))
    get(report,"divergences",nothing) == count(r->r.divergent,rows) &&
        get(report,"depth_hits",nothing) == count(r->r.tree_depth>=settings["max_depth"],rows) ||
        throw(ArgumentError("report/checkpoint sampler diagnostics differ"))
    provenance = Dict("parent_directory"=>replace(relpath(abspath(directory),root),'\\'=>'/'),
        "parent_settings_sha256"=>B.digest(settings_path),
        "parent_files_sha256"=>Dict(file=>B.digest(joinpath(path,file)) for file in files),
        "source_sha256"=>trial.source_sha256,"parent_status"=>report["status"])
    return (;rows,settings,input,report,priors,provenance)
end

"""
Select an explicit diagnostic subset: an evenly spaced grid, parameter extrema
per chain, flagged divergent draws, and their preceding retained states.
A flagged retained state is not necessarily the leapfrog state that diverged.
"""
function select_draws(rows;per_chain=3,tail_parameters=[:yc11])
    per_chain isa Integer && !(per_chain isa Bool) && per_chain>0 ||
        throw(ArgumentError("per_chain must be a positive integer"))
    isempty(rows) && throw(ArgumentError("no draws"))
    tail_parameters isa AbstractVector && all(n->n isa Symbol,tail_parameters) &&
        length(unique(tail_parameters))==length(tail_parameters) &&
        all(n->n in keys(first(rows).parameters),tail_parameters) ||
        throw(ArgumentError("invalid tail parameter selection"))
    reasons = Dict{Int,Set{String}}()
    add(i,reason) = push!(get!(()->Set{String}(),reasons,i),reason)
    for chain in sort(unique(r.chain for r in rows))
        indices = findall(r->r.chain==chain,rows)
        n = min(per_chain,length(indices))
        grid = n==1 ? [cld(length(indices),2)] : round.(Int,range(1,length(indices);length=n))
        for i in indices[grid]
            add(i,"grid")
        end
        for name in tail_parameters
            values = [rows[i].parameters[name] for i in indices]
            add(indices[argmin(values)],"$(name)_minimum")
            add(indices[argmax(values)],"$(name)_maximum")
        end
        for (j,i) in enumerate(indices)
            if rows[i].divergent
                add(i,"divergent_retained")
                j>1 && add(indices[j-1],"preceding_divergent")
            end
        end
    end
    indices = sort(collect(keys(reasons));by=i->(rows[i].chain,rows[i].draw))
    return [(;rows[i]...,reasons=join(sort(collect(reasons[i])),';')) for i in indices]
end

"""Negative log density in z=log(parameters), including the log-Jacobian."""
function log_coordinate_objective(model,data,priors,z)
    names = keys(model.initial)
    parameters = NamedTuple{names}(Tuple(exp.(z)))
    return gaussian_nll(model,data,parameters) -
        sum(logpdf(priors[n].distribution,parameters[n]) for n in names) - sum(z)
end

"""Numerical gradient checks at a positive point; numerical MLE bounds are unused."""
function point_check(model,data,priors,parameters;fd_step=1e-4)
    isfinite(fd_step) && fd_step>0 || throw(ArgumentError("invalid finite-difference step"))
    names = keys(model.initial)
    Set(keys(parameters)) == Set(names) == Set(keys(priors)) ||
        throw(ArgumentError("parameter/prior names differ"))
    all(n->isfinite(parameters[n]) && parameters[n]>0,names) ||
        throw(ArgumentError("parameters must be finite and positive"))
    z = log.([parameters[n] for n in names])
    result = Dict{String,Any}("status"=>"failed","fd_step"=>fd_step,"diagnostics_accept_fit"=>false,
        "parameter_order"=>String.(collect(names)))
    try
        prediction = model.predict(parameters,data)
        isnothing(prediction) && throw(FermentModels.CoreModels._PredictionFailure("simulation failed"))
        length(prediction) == length(data.values) || throw(DimensionMismatch("prediction length differs"))
        all(isfinite,prediction) ||
            throw(FermentModels.CoreModels._PredictionFailure("nonfinite latent prediction"))
        objective = x->log_coordinate_objective(model,data,priors,x)
        value = objective(z)
        gradient = ForwardDiff.gradient(objective,z)
        finite_differences = map([fd_step,fd_step/2]) do step
            fd = map(eachindex(z)) do i
                a,b = copy(z),copy(z)
                a[i] += step
                b[i] -= step
                (objective(a)-objective(b))/(2step)
            end
            Dict("step"=>step,"gradient"=>fd,
                "max_scaled_error"=>maximum(abs.(fd-gradient)./max.(1.,abs.(gradient))))
        end
        isfinite(value) && all(isfinite,gradient) &&
            all(all(isfinite,c["gradient"]) for c in finite_differences) ||
            throw(FermentModels.CoreModels._PredictionFailure("nonfinite density or derivative"))
        merge!(result,Dict("status"=>"complete","nll"=>gaussian_nll(model,data,parameters),
            "negative_log_density_log_coordinates"=>value,"gradient_log"=>gradient,
            "prediction_normalized"=>prediction,"finite_differences"=>finite_differences))
    catch err
        err isa FermentModels.CoreModels._PredictionFailure || rethrow()
        result["error"] = sprint(showerror,err)
    end
    return result
end

"""Compare saved and tighter tolerances, including the alternate integrator."""
function diagnose_bayes(directory;output,trials=nothing,per_chain=3,
        tail_parameters=[:yc11],tolerance_factor=.01,fd_step=1e-4)
    ispath(output) && throw(ArgumentError("use a new diagnostic output directory"))
    isfinite(tolerance_factor) && 0<tolerance_factor<1 ||
        throw(ArgumentError("tolerance_factor must be in (0,1)"))
    isfinite(fd_step) && fd_step>0 || throw(ArgumentError("invalid finite-difference step"))
    settings = TOML.parsefile(joinpath(directory,"settings.toml"))
    available = filter(t->string(t.trial) in settings["trials"],
        B.benchmark_trials(;scope=Symbol(settings["scope"])))
    selected = C.F.select_trials(available,trials)
    isempty(selected) && throw(ArgumentError("no selected trials"))
    saved = Dict(t.trial=>saved_draws(directory,t) for t in selected)
    selections = Dict(t.trial=>select_draws(saved[t.trial].rows;per_chain,tail_parameters) for t in selected)
    initial = NamedTuple{NAMES}(ntuple(_->.5,length(NAMES)))
    lower = NamedTuple{NAMES}(ntuple(_->1e-5,length(NAMES)))
    upper = NamedTuple{NAMES}(ntuple(_->5.,length(NAMES)))
    solver = Symbol(settings["solver"])
    cases = [("saved",solver,1.),("tighter",solver,tolerance_factor),
        ("alternate_tighter",solver==:tsit5 ? :rodas5p : :tsit5,tolerance_factor)]
    models = Dict((t.trial,label)=>cocoa_fit_model(;initial,lower,upper,solver=algorithm,
        abstol=factor*settings["abstol"],reltol=factor*settings["reltol"],
        maxiters=saved[t.trial].input["solver_maxiters"],
        maxiters_scope=Symbol(get(settings,"solver_maxiters_scope","solve")),
        max_step_attempts=get(settings,"solver_max_step_attempts",500)) for
        t in selected for (label,algorithm,factor) in cases)
    implementation = ("scripts/diagnose_cocoa_bayes.jl","scripts/fit_cocoa_benchmark_bayes.jl",
        "scripts/fit_cocoa_benchmark.jl","src/core/bayes.jl",C.F.LIKELIHOOD_IMPLEMENTATION...)
    metadata = Dict("parent_directory"=>replace(relpath(abspath(directory),B.ROOT),'\\'=>'/'),
        "parent_settings_sha256"=>B.digest(joinpath(directory,"settings.toml")),
        "trials"=>String.([t.trial for t in selected]),"per_chain_grid"=>per_chain,
        "tail_parameters"=>String.(tail_parameters),"tolerance_factor"=>tolerance_factor,"fd_step"=>fd_step,
        "solver_maxiters_scope"=>get(settings,"solver_maxiters_scope","solve"),
        "solver_max_step_attempts"=>get(settings,"solver_max_step_attempts",500),
        "density_definition"=>"Gaussian NLL minus source log priors minus sum(z); parameters=exp(z).",
        "error_definition"=>"max(abs(FD-AD)/max(1,abs(AD))); cross-case gradient uses saved AD in denominator.",
        "diagnostics_accept_fit"=>false,"julia_version"=>string(VERSION),
        "limitations"=>"Selected retained points only; no global numerical certificate or MCMC convergence claim. A divergence flag does not identify its failed leapfrog state.",
        "implementation_sha256"=>Dict(f=>B.digest(joinpath(B.ROOT,f)) for f in implementation))
    mkpath(output)
    C.F.write_toml(joinpath(output,"settings.toml"),metadata)
    rows = Tuple[]
    columns = ("trial","chain","draw","reasons","case","solver","abstol","reltol","solver_maxiters",
        "status","nll","negative_log_density_log_coordinates","fd_scaled_error","fd_half_step_scaled_error",
        "delta_nll_from_saved","max_prediction_difference_from_saved","max_scaled_gradient_difference_from_saved")
    for trial in selected
        path = joinpath(output,string(trial.trial))
        mkpath(path)
        C.F.write_toml(joinpath(path,"provenance.toml"),saved[trial.trial].provenance)
        for draw in selections[trial.trial]
            reference = nothing
            for (label,algorithm,factor) in cases
                check = point_check(models[(trial.trial,label)],trial.data,
                    saved[trial.trial].priors,draw.parameters;fd_step)
                label=="saved" && (reference=check)
                if check["status"]=="complete" && reference["status"]=="complete"
                    check["delta_nll_from_saved"] = check["nll"]-reference["nll"]
                    check["max_prediction_difference_from_saved"] =
                        maximum(abs.(check["prediction_normalized"]-reference["prediction_normalized"]))
                    check["max_scaled_gradient_difference_from_saved"] =
                        maximum(abs.(check["gradient_log"]-reference["gradient_log"])./
                            max.(1.,abs.(reference["gradient_log"])))
                end
                merge!(check,Dict("chain"=>draw.chain,"draw"=>draw.draw,"reasons"=>draw.reasons,
                    "case"=>label,"solver"=>string(algorithm),"abstol"=>factor*settings["abstol"],
                    "reltol"=>factor*settings["reltol"],"solver_maxiters"=>saved[trial.trial].input["solver_maxiters"],
                    "parameters"=>Dict(string(n)=>draw.parameters[n] for n in NAMES)))
                C.F.write_toml(joinpath(path,"chain_$(lpad(draw.chain,3,'0'))_draw_$(lpad(draw.draw,6,'0'))_$label.toml"),check)
                fd = get(check,"finite_differences",nothing)
                push!(rows,(trial.trial,draw.chain,draw.draw,draw.reasons,label,algorithm,
                    check["abstol"],check["reltol"],check["solver_maxiters"],check["status"],
                    get(check,"nll",missing),get(check,"negative_log_density_log_coordinates",missing),
                    isnothing(fd) ? missing : fd[1]["max_scaled_error"],
                    isnothing(fd) ? missing : fd[2]["max_scaled_error"],
                    (get(check,k,missing) for k in ("delta_nll_from_saved",
                        "max_prediction_difference_from_saved","max_scaled_gradient_difference_from_saved"))...))
                B.write_csv(joinpath(output,"comparison.csv"),columns,rows)
            end
            println("Checked ",trial.trial," chain ",draw.chain," draw ",draw.draw)
            flush(stdout)
        end
    end
    return rows
end

end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS) in (2,3) || error("Usage: diagnose_cocoa_bayes.jl FIT_DIRECTORY NEW_OUTPUT [TRIAL,TRIAL,...]")
    CocoaBayesDiagnostics.diagnose_bayes(ARGS[1];output=ARGS[2],
        trials=length(ARGS)==3 ? Symbol.(split(ARGS[3],',')) : nothing)
end
