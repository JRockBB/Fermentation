module CocoaBenchmarkFits

using FermentModels
using Random
using Statistics
using TOML
include("prepare_cocoa_benchmark.jl")
using .CocoaBenchmarkData

const LIKELIHOOD_IMPLEMENTATION = ("scripts/prepare_cocoa_benchmark.jl",
    "src/cocoa/baseline.jl", "src/cocoa/inference.jl", "src/core/estimation.jl")

"""Select trials without changing their scope order or the corresponding RNG seeds."""
function select_trials(available, ids)
    isnothing(ids) && return available
    ids isa AbstractVector && !isempty(ids) && all(id->id isa Symbol,ids) ||
        throw(ArgumentError("trials must be a nonempty vector of trial Symbols"))
    length(unique(ids)) == length(ids) || throw(ArgumentError("duplicate trial identifiers"))
    issubset(Set(ids),Set(t.trial for t in available)) ||
        throw(ArgumentError("requested trial is absent from this scope"))
    return filter(t->t.trial in ids,available)
end

"""
Load a saved candidate only when its source, state order, data preparation, model,
and likelihood match the current trial. Report and start checkpoint must agree.
Previous acceptance is never transferred to a new fit.
"""
function saved_candidate(fit_directory, trial; root=CocoaBenchmarkData.ROOT)
    settings_path = joinpath(fit_directory,"settings.toml")
    input_path = joinpath(fit_directory,string(trial.trial),"input.toml")
    report_path = joinpath(fit_directory,string(trial.trial),"report.toml")
    settings = TOML.parsefile(settings_path)
    get(settings,"method",nothing) == "mle" || throw(ArgumentError("restart requires an MLE run"))
    for path in LIKELIHOOD_IMPLEMENTATION
        get(get(settings,"implementation_sha256",Dict()),path,nothing) ==
            CocoaBenchmarkData.digest(joinpath(root,path)) ||
            throw(ArgumentError("saved likelihood implementation differs: $path"))
    end
    input = TOML.parsefile(input_path)
    expected = Dict("source_file"=>trial.source_file,"source_sha256"=>trial.source_sha256,
        "source_branch"=>trial.branch,"hidden_lactate"=>trial.hidden_lactate,
        "initial_normalized"=>trial.data.initial,"times_h"=>trial.data.times,
        "state_order"=>String.(collect(COCOA_STATES)))
    all(get(input,k,nothing)==v for (k,v) in expected) ||
        throw(ArgumentError("saved input differs from trial $(trial.trial)"))
    report = TOML.parsefile(report_path)
    get(report,"trial",nothing) == string(trial.trial) || throw(ArgumentError("saved trial differs"))
    start = get(report,"candidate_start",0)
    start isa Integer && !(start isa Bool) && start > 0 ||
        throw(ArgumentError("saved trial has no candidate"))
    checkpoint_path = joinpath(fit_directory,string(trial.trial),"start_$(lpad(start,3,'0')).toml")
    checkpoint = TOML.parsefile(checkpoint_path)
    parameters = get(report,"candidate_parameters",Dict())
    names = (COCOA_PARAMETERS..., :sigma)
    Set(keys(parameters)) == Set(String.(names)) ||
        throw(ArgumentError("candidate must contain all 24 kinetics and sigma"))
    parameters == get(checkpoint,"parameters",nothing) ||
        throw(ArgumentError("candidate report and start checkpoint disagree"))
    all(v->v isa Real && !(v isa Bool) && isfinite(v) && v > 0,values(parameters)) ||
        throw(ArgumentError("candidate parameters must be finite and positive"))
    point = NamedTuple{names}(Tuple(Float64(parameters[string(n)]) for n in names))
    provenance = Dict("directory"=>replace(relpath(abspath(fit_directory),root),'\\'=>'/'),
        "candidate_start"=>start,"previous_status"=>get(report,"status","unknown"),
        "previous_accepted"=>get(report,"accepted",false),
        "files_sha256"=>Dict(replace(relpath(path,fit_directory),'\\'=>'/')=>
            CocoaBenchmarkData.digest(path) for path in
            (settings_path,input_path,report_path,checkpoint_path)))
    return (;parameters=point,settings,provenance)
end

"""Move boundary candidates a recorded, small distance into the log search box."""
function restart_initial(parameters, lower, upper; margin=1e-8)
    0 < margin < 0.5 || throw(ArgumentError("invalid restart margin"))
    names = keys(lower)
    Set(keys(parameters)) == Set(names) == Set(keys(upper)) ||
        throw(ArgumentError("restart parameter names differ"))
    all(n->lower[n] <= parameters[n] <= upper[n],names) ||
        throw(ArgumentError("restart candidate lies outside current bounds"))
    adjusted = Symbol[]
    initial = NamedTuple{names}(map(names) do name
        lo, hi = log(lower[name]), log(upper[name])
        interior_lo, interior_hi = exp(lo+margin*(hi-lo)), exp(hi-margin*(hi-lo))
        value = clamp(parameters[name],interior_lo,interior_hi)
        value == parameters[name] || push!(adjusted,name)
        value
    end)
    return (;initial,adjusted,margin)
end

"""Training scores at a specified parameter point; never substitutes for PSIS-LOO."""
function training_scores(model, data, parameters)
    prediction = model.predict(parameters, data)
    isnothing(prediction) && return nothing
    indices = findall(!ismissing, data.values)
    residual = prediction[indices] - data.values[indices]
    return (; n=length(indices), nll=gaussian_nll(model,data,parameters),
        normalized_rmse=sqrt(sum(abs2,residual)/length(indices)),
        prediction, indices)
end

function run_record(run)
    Dict("objective"=>run.objective, "accepted"=>run.converged, "retcode"=>run.retcode,
        "optimizer_converged"=>run.optimizer_converged, "stationarity_passed"=>run.stationarity_passed,
        "projected_gradient_inf"=>ismissing(run.projected_gradient_inf) ? "unknown" : run.projected_gradient_inf,
        "gradient_tolerance"=>run.gradient_tolerance, "elapsed_seconds"=>run.elapsed_seconds,
        "time_limit_reached"=>ismissing(run.time_limit_reached) ? "unknown" : run.time_limit_reached,
        "error"=>isnothing(run.error) ? "" : run.error,
        "initial"=>Dict(string(k)=>v for (k,v) in pairs(run.initial)),
        "parameters"=>isnothing(run.parameters) ? Dict{String,Float64}() :
            Dict(string(k)=>v for (k,v) in pairs(run.parameters)))
end

function write_toml(path, value)
    open(path,"w") do io
        TOML.print(io,value;sorted=true)
    end
end

"""
Fit the unchanged baseline sequentially to a benchmark or all source trials.
Output must be new. MLE bounds, initialization, budgets, and stationarity are
computational settings, not biological estimates. Every completed start is saved.
Unaccepted candidates retain training scores with an explicit status.
"""
function run_benchmark(; output, scope=:benchmark, starts=1, maxiters=300,
        outer_iterations=10, maxtime=120.0, seed=141085, solver=:rodas5p,
        abstol=1e-9, reltol=1e-8, solver_maxiters=10_000, gtol=1e-6,
        trials=nothing, restart_from=nothing)
    ispath(output) && throw(ArgumentError("Use a new output directory to preserve previous runs"))
    starts isa Integer && !(starts isa Bool) && starts > 0 || throw(ArgumentError("starts must be positive"))
    maxiters isa Integer && !(maxiters isa Bool) && maxiters > 0 || throw(ArgumentError("maxiters must be positive"))
    outer_iterations isa Integer && !(outer_iterations isa Bool) && outer_iterations > 0 ||
        throw(ArgumentError("outer_iterations must be positive"))
    all(v->v isa Real && !(v isa Bool) && isfinite(v) && v > 0,(maxtime,gtol)) ||
        throw(ArgumentError("maxtime and gtol must be finite and positive"))
    seed isa Integer && !(seed isa Bool) && seed >= 0 || throw(ArgumentError("seed must be nonnegative"))
    available = CocoaBenchmarkData.benchmark_trials(; scope)
    selected = select_trials(available,trials)
    trial_seeds = Dict(t.trial=>seed+i-1 for (i,t) in enumerate(available))
    published = Dict(r.trial=>r for r in CocoaBenchmarkData.published_baseline())
    names = (COCOA_PARAMETERS..., :sigma)
    # Prior location is a reproducible optimizer start, NOT a parameter estimate.
    initial = NamedTuple{names}(ntuple(_->0.5,25))
    lower = NamedTuple{names}(ntuple(_->1e-5,25))
    upper = NamedTuple{names}(ntuple(_->5.0,25))
    model = cocoa_fit_model(; initial,lower,upper,solver,abstol,reltol,maxiters=solver_maxiters)
    root = CocoaBenchmarkData.ROOT
    # Validate every restart before creating output or beginning a potentially long batch.
    restarts = isnothing(restart_from) ? Dict() :
        Dict(t.trial=>saved_candidate(restart_from,t) for t in selected)
    initials = Dict(t.trial=>(isnothing(restart_from) ?
        (;initial,adjusted=Symbol[],margin=0.0) :
        restart_initial(restarts[t.trial].parameters,lower,upper)) for t in selected)
    implementation = ("scripts/fit_cocoa_benchmark.jl",LIKELIHOOD_IMPLEMENTATION...)
    settings = Dict("scope"=>string(scope), "method"=>"mle", "starts"=>starts,
        "maxiters"=>maxiters, "outer_iterations"=>outer_iterations, "maxtime_seconds"=>maxtime,
        "seed"=>seed, "solver"=>string(solver), "abstol"=>abstol, "reltol"=>reltol,
        "solver_maxiters"=>solver_maxiters,"gradient_tolerance"=>gtol,
        "trials"=>String.([t.trial for t in selected]),
        "trial_seeds"=>Dict(string(k)=>trial_seeds[k] for k in (t.trial for t in selected)),
        "initialization"=>isnothing(restart_from) ?
            "0.5 for all parameters; numerical starting point based on source kinetic-prior location; not a biological estimate" :
            "verified saved candidate first; additional starts uniform in log bounds; acceptance recomputed",
        "mle_bounds"=>[1e-5,5.0], "julia_version"=>string(VERSION),
        "paper_metric"=>"PSIS-LOO and relative BMA weights; not directly comparable to MLE training scores",
        "source"=>"doi:10.1098/rsos.210274; doi:10.5281/zenodo.5510370",
        "implementation_sha256"=>Dict(p=>CocoaBenchmarkData.digest(joinpath(root,p)) for p in implementation))
    isnothing(restart_from) || (settings["restart_from"] =
        replace(relpath(abspath(restart_from),root),'\\'=>'/'))
    mkpath(output)
    write_toml(joinpath(output,"settings.toml"),settings)
    columns = ("trial","in_paper_benchmark","observations","status","candidate_nll",
        "candidate_normalized_rmse","projected_gradient_inf","elapsed_seconds",
        "published_baseline_success","published_max_rhat","published_min_bulk_ess",
        "published_min_tail_ess","published_bma_weight","published_psis_loo","comparison")
    comparison_rows = Tuple[]
    for (index, trial) in enumerate(selected)
        directory = joinpath(output,string(trial.trial))
        mkpath(directory)
        input = Dict{String,Any}(
            "source_file"=>trial.source_file,"source_sha256"=>trial.source_sha256,
            "source_branch"=>trial.branch,"hidden_lactate"=>trial.hidden_lactate,
            "initial_normalized"=>trial.data.initial,"times_h"=>trial.data.times,
            "state_order"=>String.(collect(COCOA_STATES)),
            "seed"=>trial_seeds[trial.trial],
            "optimizer_initial"=>Dict(string(k)=>v for (k,v) in pairs(initials[trial.trial].initial)))
        if !isnothing(restart_from)
            input["restart"] = merge(restarts[trial.trial].provenance,Dict(
                "log_range_margin"=>initials[trial.trial].margin,
                "adjusted_parameters"=>String.(initials[trial.trial].adjusted)))
        end
        write_toml(joinpath(directory,"input.toml"),input)
        model = cocoa_fit_model(;initial=initials[trial.trial].initial,lower,upper,
            solver,abstol,reltol,maxiters=solver_maxiters)
        println("Starting ",index,"/",length(selected)," ",trial.trial); flush(stdout)
        started = time()
        result = fit_mle(model,trial.data;starts,maxiters,outer_iterations,maxtime,gtol,
            rng=MersenneTwister(trial_seeds[trial.trial]),
            on_start_complete=(i,run)->write_toml(joinpath(directory,"start_$(lpad(i,3,'0')).toml"),run_record(run)))
        finite_runs = findall(r->!isnothing(r.parameters) && isfinite(r.objective),result.runs)
        candidate_index = result.converged ? result.best_start :
            isempty(finite_runs) ? nothing : finite_runs[argmin([result.runs[i].objective for i in finite_runs])]
        status = result.converged ? "accepted_mle" :
            isnothing(candidate_index) ? "prediction_failure" : "unverified_candidate"
        candidate = isnothing(candidate_index) ? nothing : result.runs[candidate_index]
        scores = isnothing(candidate) ? nothing : training_scores(model,trial.data,candidate.parameters)
        report = Dict("trial"=>string(trial.trial),"status"=>status,"accepted"=>result.converged,
            "elapsed_seconds"=>time()-started,"candidate_start"=>something(candidate_index,0),
            "uncertainty_status"=>string(result.uncertainty.status),
            "published_psis_loo"=>"missing: archived posterior outputs absent",
            "direct_metric_comparison_available"=>false)
        if !isnothing(scores)
            report["candidate_nll"] = scores.nll
            report["candidate_normalized_rmse"] = scores.normalized_rmse
            report["candidate_parameters"] = Dict(string(k)=>v for (k,v) in pairs(candidate.parameters))
            report["projected_gradient_inf"] = ismissing(candidate.projected_gradient_inf) ? "unknown" : candidate.projected_gradient_inf
            ntime = length(trial.data.times)
            CocoaBenchmarkData.write_csv(joinpath(directory,"predictions.csv"),
                ("time_h","variable","observed_normalized","candidate_prediction_normalized"),
                [(t,state,trial.data.values[(j-1)*ntime+i],scores.prediction[(j-1)*ntime+i])
                    for (j,state) in enumerate(COCOA_STATES) for (i,t) in enumerate(trial.data.times)])
        end
        write_toml(joinpath(directory,"report.toml"),report)
        reference = get(published,trial.trial,nothing)
        push!(comparison_rows,(trial.trial,trial.in_paper_benchmark,count(!ismissing,trial.data.values),
            status,get(report,"candidate_nll",missing),get(report,"candidate_normalized_rmse",missing),
            get(report,"projected_gradient_inf",missing),report["elapsed_seconds"],
            isnothing(reference) ? false : reference.success,
            (isnothing(reference) ? missing : getproperty(reference,k) for k in (:max_rhat,:min_bulk_ess,:min_tail_ess,:bma_weight))...,
            missing,"MLE training scores versus published Bayesian diagnostics; PSIS-LOO unavailable"))
        CocoaBenchmarkData.write_csv(joinpath(output,"comparison.csv"),columns,comparison_rows)
        println("Completed ",trial.trial,": ",status,"; NLL=",get(report,"candidate_nll","unavailable"))
        flush(stdout)
    end
    return comparison_rows
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    options = Dict{Symbol,Any}()
    integer_options = (:starts,:maxiters,:outer_iterations,:seed,:solver_maxiters)
    float_options = (:maxtime,:abstol,:reltol,:gtol)
    for argument in ARGS
        parts = split(replace(argument,r"^--"=>""),'=';limit=2)
        length(parts)==2 || error("Use --name=value")
        key = Symbol(parts[1])
        key in (:output,:scope,:solver,:trials,:restart_from,integer_options...,float_options...) || error("Unknown option $key")
        options[key] = key in integer_options ? parse(Int,parts[2]) :
            key in float_options ? parse(Float64,parts[2]) : key in (:scope,:solver) ? Symbol(parts[2]) :
            key == :trials ? Symbol.(split(parts[2],',')) : parts[2]
    end
    haskey(options,:output) || error("Provide --output=<new-directory>")
    CocoaBenchmarkFits.run_benchmark(;options...)
end
