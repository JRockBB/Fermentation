module CocoaRuntimeDiagnostics

using FermentModels, ForwardDiff, OrdinaryDiffEq, Statistics, TOML
include("diagnose_cocoa_bayes.jl")
const D = CocoaBayesDiagnostics
const C,B = D.C,D.B

require(condition,message) = condition === true ? nothing : throw(ArgumentError(message))

"""Validate complete schema-2 runtime logs, including retained checkpoint correspondence."""
function read_runtime(file,settings,draws;sampler_file=joinpath(dirname(file),"sampler_statistics.csv"))
    require(get(settings,"runtime_schema_version",0)==2,"requires runtime schema 2; older logs lack sampler work statistics")
    samples,warmup,chains = settings["samples_per_chain"],settings["warmup_per_chain"],settings["chains"]
    require(all(x->x isa Integer && !(x isa Bool) && x>0,(samples,chains)) &&
        warmup isa Integer && !(warmup isa Bool) && warmup>=0,"invalid sampler dimensions")
    lines = readlines(file)
    require(!isempty(lines) && split(first(lines),',')==collect(C.RUNTIME_HEADER),"runtime header differs")
    require(length(lines)-1==chains*(samples+warmup),"incomplete runtime log")
    checkpoints = Dict((r.chain,r.draw)=>r for r in draws)
    require(length(checkpoints)==length(draws)==samples*chains,"incomplete or duplicate checkpoints")
    counts,previous = zeros(Int,chains),zeros(chains)
    global_elapsed = 0.
    rows = NamedTuple[]
    for line in lines[2:end]
        f = split(line,',';keepempty=true)
        require(length(f)==length(C.RUNTIME_HEADER),"runtime row width differs")
        chain,iteration,draw = parse(Int,f[1]),parse(Int,f[2]),parse(Int,f[4])
        require(1<=chain<=chains,"runtime chain outside range")
        require(iteration==counts[chain]+1 && iteration<=samples+warmup,"runtime call order differs")
        phase = iteration==1 ? "initialization" : iteration<=warmup ? "warmup" : "sampling"
        require(f[3]==phase && draw==max(0,iteration-warmup),"runtime phase/draw differs")
        elapsed,duration = parse(Float64,f[5]),parse(Float64,f[6])
        require(isfinite(elapsed) && isfinite(duration) && duration>=0 && elapsed>=global_elapsed &&
            isapprox(elapsed-previous[chain],duration;atol=1e-9,rtol=1e-10),"runtime timing differs")
        if iteration==1
            require(all(isempty,f[7:end]),"initialization must have missing transition statistics")
            stats = (;divergent=missing,tree_depth=missing,
                NamedTuple{C.RUNTIME_STATISTICS}(ntuple(_->missing,length(C.RUNTIME_STATISTICS)))...)
        else
            divergent,depth,n_steps = parse(Bool,f[7]),parse(Int,f[8]),parse(Int,f[9])
            step_size,acceptance,energy,error,max_error = parse.(Float64,f[10:14])
            require(0<=depth<=settings["max_depth"] && n_steps>0,"invalid tree depth or leapfrog count")
            require(isfinite(step_size) && step_size>0,"invalid leapfrog step size")
            require(0<=acceptance<=1 || (divergent && !isfinite(acceptance)),"invalid acceptance rate")
            require(all(isfinite,(energy,error,max_error)) || divergent,"nonfinite energy without divergence")
            stats = (;divergent,tree_depth=depth,n_steps,step_size,acceptance_rate=acceptance,
                hamiltonian_energy=energy,hamiltonian_energy_error=error,
                max_hamiltonian_energy_error=max_error)
        end
        if draw>0
            require(haskey(checkpoints,(chain,draw)),"runtime draw missing from checkpoints")
            kept = checkpoints[(chain,draw)]
            require(kept.divergent==coalesce(stats.divergent,false) &&
                kept.tree_depth==coalesce(stats.tree_depth,0),"runtime/checkpoint diagnostics differ")
        end
        push!(rows,(;chain,iteration,phase,draw,elapsed,duration,stats...))
        counts[chain],previous[chain],global_elapsed = iteration,elapsed,elapsed
    end
    require(all(==(samples+warmup),counts),"incomplete runtime chains")
    require(isfile(sampler_file),"missing independently exported sampler statistics")
    raw = readlines(sampler_file)
    require(!isempty(raw) && split(first(raw),',')==collect(C.SAMPLER_HEADER) &&
        length(raw)==samples*chains+1,"sampler statistics header or count differs")
    events = Dict((r.chain,r.draw)=>r for r in rows if r.draw>0)
    keys = (:divergent,:tree_depth,C.RUNTIME_STATISTICS...)
    for (i,line) in enumerate(raw[2:end])
        f = split(line,',';keepempty=true)
        require(length(f)==length(C.SAMPLER_HEADER),"sampler statistics row width differs")
        chain,draw = cld(i,samples),mod1(i,samples)
        require(parse(Int,f[1])==chain && parse(Int,f[2])==draw,"sampler statistics order differs")
        event = events[(chain,draw)]
        for (key,field) in zip(keys,f[3:end])
            value = isempty(field) ? missing : key==:divergent ? parse(Bool,field) :
                key in (:tree_depth,:n_steps) ? parse(Int,field) : parse(Float64,field)
            require(isequal(value,getproperty(event,key)),"runtime/raw-chain statistic differs: $key")
        end
    end
    return rows
end

"""Descriptive timings by chain and phase; seconds per leapfrog are not ODE timings."""
function summarize_runtime(rows)
    summaries = Tuple[]
    for chain in sort(unique(r.chain for r in rows)), phase in ("initialization","warmup","sampling")
        selected = filter(r->r.chain==chain && r.phase==phase,rows)
        isempty(selected) && continue
        duration = [r.duration for r in selected]
        transitions = filter(r->!ismissing(r.n_steps),selected)
        steps = sum((r.n_steps for r in transitions);init=0)
        push!(summaries,(chain,phase,length(selected),median(duration),quantile(duration,.95),maximum(duration),
            steps,isempty(transitions) ? missing : sum(r.duration for r in transitions)/steps,
            isempty(transitions) ? missing : minimum(r.step_size for r in transitions),
            isempty(transitions) ? missing : maximum(r.step_size for r in transitions),
            isempty(transitions) ? missing : mean(r.acceptance_rate for r in transitions),
            count(r->r.divergent===true,transitions),
            count(r->!isfinite(r.max_hamiltonian_energy_error),transitions)))
    end
    return summaries
end

"""Select slow retained calls, their predecessors, and a median-duration control per chain."""
function select_runtime_draws(events,draws;per_chain=1)
    require(per_chain isa Integer && !(per_chain isa Bool) && per_chain>0,"per_chain must be positive")
    lookup = Dict((r.chain,r.draw)=>r for r in draws)
    reasons = Dict{Tuple{Int,Int},Set{String}}()
    add(key,reason) = push!(get!(()->Set{String}(),reasons,key),reason)
    for chain in sort(unique(r.chain for r in events))
        rows = filter(r->r.chain==chain && r.phase=="sampling",events)
        isempty(rows) && continue
        ordered = sort(rows;by=r->(-r.duration,r.draw))
        for row in first(ordered,min(per_chain,length(ordered)))
            add((chain,row.draw),"slow_retained_call")
            row.draw>1 && add((chain,row.draw-1),"preceding_slow_call")
        end
        target = median([r.duration for r in rows])
        control = first(sort(rows;by=r->(abs(r.duration-target),r.draw)))
        add((chain,control.draw),"median_duration_control")
    end
    return [(;lookup[key]...,reasons=join(sort(collect(reasons[key])),';'))
        for key in sort(collect(keys(reasons)))]
end

# Warm each callable before measuring. Preserve each repeat, including GC time;
# timings are descriptive and never influence a model or an acceptance decision.
function measure_repeated(f,repeats)
    require(repeats isa Integer && !(repeats isa Bool) && repeats>0,"repeats must be positive")
    f()
    samples = [(@timed f()) for _ in 1:repeats]
    return (value=last(samples).value,seconds=[s.time for s in samples],
        bytes=[s.bytes for s in samples],gc_seconds=[s.gctime for s in samples])
end

"""Replay a saved state, measuring primal ODE and log-coordinate gradient costs separately."""
function profile_point(trial,saved,parameters;repeats=3)
    settings = saved.settings
    names = keys(parameters)
    initial = NamedTuple{names}(ntuple(_->.5,length(names)))
    lower,upper = NamedTuple{names}(ntuple(_->1e-5,length(names))),NamedTuple{names}(ntuple(_->5.,length(names)))
    solver = Symbol(settings["solver"])
    require(solver in (:tsit5,:rodas5p),"unknown solver")
    options = (abstol=settings["abstol"],reltol=settings["reltol"],
        maxiters=settings["trial_solver_maxiters"][string(trial.trial)] )
    scope = Symbol(get(settings,"solver_maxiters_scope","solve"))
    attempts = get(settings,"solver_max_step_attempts",500)
    model = cocoa_fit_model(;initial,lower,upper,solver,options...,maxiters_scope=scope,max_step_attempts=attempts)
    algorithm = solver==:tsit5 ? Tsit5() : Rodas5P()
    function trajectory()
        problem = cocoa_problem(trial.data.initial,(trial.data.t0,last(trial.data.times)),parameters)
        solution = if scope==:output_interval
            FermentModels.Cocoa._solve_output_intervals(problem,algorithm,trial.data.times;
                options...,max_step_attempts=attempts)
        else
            solved = solve(problem,algorithm;options...,saveat=trial.data.times,save_start=false,dense=false)
            successful_retcode(solved) ? solved : nothing
        end
        isnothing(solution) && throw(FermentModels.CoreModels._PredictionFailure("saved state failed to integrate"))
        return solution
    end
    result = Dict{String,Any}("status"=>"failed","diagnostics_accept_fit"=>false,
        "parameter_order"=>String.(collect(names)),"parameters"=>Dict(string(k)=>v for (k,v) in pairs(parameters)))
    try
        ode = measure_repeated(trajectory,repeats)
        prediction = vec(permutedims(Array(ode.value)))
        require(prediction==model.predict(parameters,trial.data),"profile trajectory differs from likelihood prediction")
        z = log.(collect(values(parameters)))
        objective = x->D.log_coordinate_objective(model,trial.data,saved.priors,x)
        gradient = measure_repeated(()->ForwardDiff.gradient(objective,z),repeats)
        require(all(isfinite,gradient.value) && isfinite(objective(z)),"nonfinite posterior derivative")
        merge!(result,Dict("status"=>"complete","prediction_normalized"=>prediction,
            "nll"=>gaussian_nll(model,trial.data,parameters),"gradient_log"=>gradient.value,
            "accepted_ode_steps"=>ode.value.stats.naccept,"rejected_ode_steps"=>ode.value.stats.nreject,
            "rhs_evaluations"=>ode.value.stats.nf,
            "ode_seconds"=>ode.seconds,"ode_bytes"=>ode.bytes,"ode_gc_seconds"=>ode.gc_seconds,
            "gradient_seconds"=>gradient.seconds,"gradient_bytes"=>gradient.bytes,"gradient_gc_seconds"=>gradient.gc_seconds))
    catch err
        err isa FermentModels.CoreModels._PredictionFailure || rethrow()
        result["error"] = sprint(showerror,err)
    end
    return result
end

"""Audit completed runtime records and profile a provenance-checked retained-state subset."""
function diagnose_runtime(directory;output,trials=nothing,per_chain=1,repeats=3)
    require(!ispath(output),"use a new diagnostic output directory")
    require(repeats isa Integer && !(repeats isa Bool) && repeats>0,"repeats must be positive")
    require(per_chain isa Integer && !(per_chain isa Bool) && per_chain>0,"per_chain must be positive")
    settings_hash = B.digest(joinpath(directory,"settings.toml"))
    settings = TOML.parsefile(joinpath(directory,"settings.toml"))
    require(get(settings,"runtime_schema_version",0)==2,"requires runtime schema 2")
    runner = "scripts/fit_cocoa_benchmark_bayes.jl"
    require(get(settings["implementation_sha256"],runner,nothing)==B.digest(joinpath(B.ROOT,runner)),"runtime writer implementation differs")
    available = filter(t->string(t.trial) in settings["trials"],B.benchmark_trials(;scope=Symbol(settings["scope"])))
    selected = C.F.select_trials(available,trials)
    require(!isempty(selected),"no selected trials")
    saved = Dict(t.trial=>D.saved_draws(directory,t) for t in selected)
    require(all(s->s.provenance["parent_settings_sha256"]==settings_hash,values(saved)),"settings changed during diagnosis")
    runtime_hashes = Dict(string(t.trial)*"/"*file=>B.digest(joinpath(directory,string(t.trial),file))
        for t in selected for file in ("runtime.csv","sampler_statistics.csv"))
    events = Dict(t.trial=>read_runtime(joinpath(directory,string(t.trial),"runtime.csv"),settings,saved[t.trial].rows) for t in selected)
    selections = Dict(t.trial=>select_runtime_draws(events[t.trial],saved[t.trial].rows;per_chain) for t in selected)
    hashes = merge(Dict("settings.toml"=>settings_hash),runtime_hashes)
    for trial in selected
        for (file,hash) in saved[trial.trial].provenance["parent_files_sha256"]
            hashes[string(trial.trial)*"/"*file] = hash
        end
    end
    metadata = Dict("status"=>"profiling","runtime_integrity_passed"=>true,"diagnostics_accept_fit"=>false,
        "parent_directory"=>replace(relpath(abspath(directory),B.ROOT),'\\'=>'/'),"parent_sha256"=>hashes,
        "trials"=>string.([t.trial for t in selected]),"per_chain"=>per_chain,"repeats"=>repeats,
        "julia_version"=>string(VERSION),"julia_threads"=>Threads.nthreads(),
        "environment_sha256"=>Dict(f=>B.digest(joinpath(B.ROOT,f)) for f in
            ("Project.toml","Manifest-v$(VERSION.major).$(VERSION.minor).toml","Manifest.toml") if isfile(joinpath(B.ROOT,f))),
        "limitations"=>"Elapsed completion gaps include scheduling, callback work, compilation, and GC. Seconds per leapfrog are not ODE timings. Retained/predecessor states do not reconstruct the sampled Hamiltonian path or a divergent leapfrog state. Replayed costs are warmed descriptive measurements, not a causal attribution or acceptance test. Warmup parameter states are not retained. Older runtime logs cannot recover these statistics.",
        "implementation_sha256"=>Dict(f=>B.digest(joinpath(B.ROOT,f)) for f in
            ("scripts/diagnose_cocoa_runtime.jl","scripts/diagnose_cocoa_bayes.jl",runner,"src/core/bayes.jl",C.F.LIKELIHOOD_IMPLEMENTATION...)))
    mkpath(output)
    C.F.write_toml(joinpath(output,"settings.toml"),metadata)
    for trial in selected
        id = string(trial.trial)
        path = joinpath(output,id)
        mkpath(path)
        B.write_csv(joinpath(path,"runtime_summary.csv"),("chain","phase","calls","median_seconds","p95_seconds",
            "max_seconds","leapfrog_steps","seconds_per_leapfrog","min_step_size","max_step_size",
            "mean_acceptance_rate","divergences","nonfinite_max_energy_errors"),summarize_runtime(events[trial.trial]))
        B.write_csv(joinpath(path,"selection.csv"),("chain","draw","reasons"),
            [(r.chain,r.draw,r.reasons) for r in selections[trial.trial]])
        comparisons = Tuple[]
        for row in selections[trial.trial]
            report = profile_point(trial,saved[trial.trial],row.parameters;repeats)
            merge!(report,Dict("chain"=>row.chain,"draw"=>row.draw,"selection_reasons"=>row.reasons))
            C.F.write_toml(joinpath(path,"chain_$(row.chain)_draw_$(row.draw).toml"),report)
            complete = report["status"]=="complete"
            push!(comparisons,(row.chain,row.draw,report["status"],
                (get(report,k,missing) for k in ("accepted_ode_steps","rejected_ode_steps","rhs_evaluations"))...,
                complete ? median(report["ode_seconds"]) : missing,
                complete ? median(report["gradient_seconds"]) : missing))
            B.write_csv(joinpath(path,"costs.csv"),("chain","draw","status","accepted_ode_steps",
                "rejected_ode_steps","rhs_evaluations","median_ode_seconds","median_gradient_seconds"),comparisons)
        end
    end
    for (file,hash) in hashes
        require(B.digest(joinpath(directory,file))==hash,"parent changed during diagnosis")
    end
    metadata["status"] = "complete"
    metadata["runtime_events_checked"] = sum(length,values(events))
    metadata["retained_draws_checked"] = sum(s->length(s.rows),values(saved))
    metadata["profiled_states"] = sum(length,values(selections))
    C.F.write_toml(joinpath(output,"settings.toml"),metadata)
    return metadata
end

end

if abspath(PROGRAM_FILE)==@__FILE__
    2<=length(ARGS)<=3 || error("Usage: diagnose_cocoa_runtime.jl FIT_DIRECTORY NEW_OUTPUT [TRIAL,TRIAL,...]")
    result = CocoaRuntimeDiagnostics.diagnose_runtime(ARGS[1];output=ARGS[2],
        trials=length(ARGS)==3 ? Symbol.(split(ARGS[3],',')) : nothing)
    println("Runtime integrity passed; profiled ",result["profiled_states"]," retained states. Sampling acceptance unchanged.")
end
