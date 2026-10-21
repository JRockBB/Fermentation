using FermentModels
using Random
using Statistics
using TOML
using SHA
include("prepare_cocoa.jl")
using .CocoaSourceData: cocoa_trial, published_summaries, TRIALS
include("cocoa_fit_output.jl")

function mle_run_report(run)
    return Dict("objective" => run.objective, "converged" => run.converged,
        "retcode" => run.retcode, "error" => isnothing(run.error) ? "" : run.error,
        "initial" => Dict(string(k)=>v for (k,v) in pairs(run.initial)),
        "parameters" => isnothing(run.parameters) ? Dict{String,Float64}() :
            Dict(string(k)=>v for (k,v) in pairs(run.parameters)),
        "iterations" => ismissing(run.iterations) ? -1 : run.iterations,
        "elapsed_seconds" => run.elapsed_seconds,
        "optimizer_converged" => run.optimizer_converged,
        "stationarity_passed" => run.stationarity_passed,
        "projected_gradient_inf" => ismissing(run.projected_gradient_inf) ? "unknown" : run.projected_gradient_inf,
        "gradient_tolerance" => run.gradient_tolerance,
        "time_limit_reached" => ismissing(run.time_limit_reached) ? "unknown" : run.time_limit_reached)
end

"""
Run an explicitly configured source-model fit and retain diagnostics even when
convergence fails. The output directory must be new. Defaults reproduce the
source's retained sample count, warmup, chain count, target acceptance, and depth.
Julia/Turing uses prior initialization, not the original absolute-Normal starts.
"""
function run_cocoa_fit(trial::Symbol; method=:bayes, output,
                       samples=2000, warmup=1000, chains=4, starts=10,
                       maxiters=1000, outer_iterations=20, target_accept=0.985,
                       max_depth=12, seed=100785, abstol=1e-6, reltol=1e-6,
                       solver_maxiters=100_000_000, maxtime=nothing, solver=:tsit5,
                       chain_execution=:serial, initial_file=nothing)
    trial in TRIALS || throw(ArgumentError("choose camu, box1, or box2"))
    method in (:mle,:bayes) || throw(ArgumentError("choose mle or bayes"))
    method == :bayes && !isnothing(maxtime) && throw(ArgumentError("maxtime applies to MLE starts only"))
    method == :bayes && !isnothing(initial_file) && throw(ArgumentError("initial_file applies to MLE starts only"))
    chain_execution in (:serial,:threads) || throw(ArgumentError("choose serial or threads for chain_execution"))
    method == :mle && chain_execution != :serial && throw(ArgumentError("chain_execution applies to Bayesian fits only"))
    ispath(output) && throw(ArgumentError("output directory exists; preserve prior fit runs"))
    root = normpath(joinpath(@__DIR__, ".."))
    source = cocoa_trial(trial)
    data = cocoa_observations(source.initial, source.times, source.observations)
    parameter_file = joinpath(root, "params", "cocoa_$(trial)_scaled.toml")
    initial = load_params(parameter_file).values
    isnothing(initial_file) || (initial = cocoa_fit_initial(initial,initial_file))
    # Numerical MLE search bounds, not biological estimates or Bayesian support.
    lower = NamedTuple{keys(initial)}(ntuple(_ -> 1e-5,25))
    upper = NamedTuple{keys(initial)}(ntuple(_ -> 5.,25))
    model = cocoa_fit_model(; initial, lower, upper, abstol, reltol, maxiters=solver_maxiters, solver)
    rng = MersenneTwister(seed)
    mkpath(output)
    settings = Dict("trial" => String(trial), "method" => String(method), "seed" => seed,
        "samples" => samples, "warmup" => warmup, "chains" => chains, "starts" => starts,
        "maxiters" => maxiters, "outer_iterations" => outer_iterations,
        "target_accept" => target_accept, "max_depth" => max_depth,
        "abstol" => abstol, "reltol" => reltol, "solver" => solver == :tsit5 ? "Tsit5" : "Rodas5P",
        "solver_maxiters" => solver_maxiters,
        "mle_lower" => 1e-5, "mle_upper" => 5.0,
        "mle_gradient_tolerance" => 1e-6,
        "mle_initialization" => isnothing(initial_file) ?
            "published means first; remaining starts uniform in log bounds" :
            "parameters from initial_file first; remaining starts uniform in log bounds",
        "bayesian_bounds" => "positive and unbounded above; MLE bounds are not imposed",
        "bayesian_initialization" => "independent prior draws",
        "chain_execution" => String(chain_execution), "julia_threads" => Threads.nthreads(),
        "julia_version" => string(VERSION), "script" => "scripts/fit_cocoa.jl",
        "source" => "doi:10.1098/rsos.180964; doi:10.5061/dryad.321d33v",
        "raw_file" => replace(source.source_file, '\\'=>'/'),
        "raw_sha256" => bytes2hex(open(sha256, joinpath(root, source.source_file))),
        "parameter_file_sha256" => bytes2hex(open(sha256, parameter_file)))
    settings["implementation_sha256"] = Dict(path => bytes2hex(open(sha256,joinpath(root,path)))
        for path in ("src/cocoa/baseline.jl", "src/cocoa/inference.jl",
            "src/core/estimation.jl", "src/core/bayes.jl", "scripts/fit_cocoa.jl",
            "scripts/cocoa_fit_output.jl"))
    settings["mle_time_limit_seconds"] = isnothing(maxtime) ? "unlimited" : maxtime
    if !isnothing(initial_file)
        settings["initial_file"] = replace(relpath(abspath(initial_file),root),'\\'=>'/')
        settings["initial_file_sha256"] = bytes2hex(open(sha256,initial_file))
    end
    open(joinpath(output, "settings.toml"), "w") do io
        TOML.print(io, settings; sorted=true)
    end
    println("Starting ", method, " fit for ", trial, "; ", length(data.values), " likelihood values and 25 parameters.")
    flush(stdout)
    start_time = time()
    result = if method == :mle
        function save_start(index, run)
            path = joinpath(output, "mle_start_$(lpad(index,3,'0')).toml")
            open(path, "w") do io
                TOML.print(io, mle_run_report(run); sorted=true)
            end
            println("Completed start ", index, "/", starts, ": ", run.retcode,
                "; NLL=", run.objective, "; seconds=", round(run.elapsed_seconds; digits=2),
                "; time limit=", run.time_limit_reached,
                "; stationarity=", run.stationarity_passed,
                "; projected gradient=", run.projected_gradient_inf)
            flush(stdout)
        end
        fit_mle(model, data; starts, rng, maxiters, outer_iterations, gtol=1e-6, maxtime,
            on_start_complete=save_start)
    else
        column_names = keys(model.initial)
        for chain in 1:chains
            open(joinpath(output,"chain_$(lpad(chain,3,'0')).csv"),"w") do io
                write_cocoa_draw_header(io,column_names)
            end
        end
        function save_draw(chain,draw,parameters,diagnostics)
            open(joinpath(output,"chain_$(lpad(chain,3,'0')).csv"),"a") do io
                write_cocoa_draw(io,chain,draw,parameters,diagnostics,column_names)
            end
            if draw == 1 || draw % 100 == 0 || draw == samples
                println("Chain ",chain,": saved ",draw,"/",samples," retained draws.")
                flush(stdout)
            end
        end
        fit_bayes(model, data, cocoa_priors(); samples, warmup, chains, target_accept, max_depth,
            predictive_draws=samples*chains, rng, progress=true, chain_execution, on_draw=save_draw)
    end
    report = Dict{String,Any}("elapsed_seconds" => time()-start_time,
        "method" => String(method), "trial" => String(trial))
    if method == :mle
        report["converged"] = result.converged
        report["uncertainty_status"] = String(result.uncertainty.status)
        report["objective"] = result.objective
        report["best_start"] = isnothing(result.best_start) ? 0 : result.best_start
        report["runs"] = mle_run_report.(result.runs)
        estimates = result.parameters
    else
        report["converged"] = result.diagnostics.passed
        report["divergences"] = result.diagnostics.divergences
        report["depth_hits"] = result.diagnostics.depth_hits
        report["max_rhat"] = maximum(row.rhat for row in values(result.summary))
        report["min_ess_bulk"] = minimum(row.ess_bulk for row in values(result.summary))
        report["min_ess_tail"] = minimum(row.ess_tail for row in values(result.summary))
        report["rhat_threshold"] = result.diagnostics.rhat_threshold
        report["minimum_ess"] = result.diagnostics.minimum_ess
        estimates = NamedTuple{keys(result.summary)}(Tuple(row.mean for row in values(result.summary)))
        open(joinpath(output,"draws.csv"),"w") do io
            println(io,"chain,draw,",join(String.(keys(result.draws)),','))
            for chain in 1:chains, draw in 1:samples
                println(io,chain,',',draw,',',join((result.draws[name][draw,chain] for name in keys(result.draws)),','))
            end
        end
        open(joinpath(output,"predictive.csv"),"w") do io
            println(io,"time_h,variable,observed,latent_mean,latent_median,latent_lower,latent_upper,observation_lower,observation_upper,unit")
            for (i,state) in enumerate(COCOA_STATES), (j,t) in enumerate(data.times)
                row = (i-1)*length(data.times)+j
                scale = data.scales[i]
                values = (t,state,source.observations[j,i],
                    mean(result.predictive.latent_samples[row,:])*scale,
                    result.predictive.latent.median[row]*scale,
                    result.predictive.latent.lower[row]*scale,
                    result.predictive.latent.upper[row]*scale,
                    result.predictive.observed.lower[row]*scale,
                    result.predictive.observed.upper[row]*scale,"mg/g pulp")
                println(io,join(values,','))
            end
        end
    end
    if !isnothing(estimates)
        report["parameters"] = Dict(String(k)=>v for (k,v) in pairs(estimates))
        reference = filter(row -> row.experiment == trial, published_summaries())
        inside = [row.q025 <= estimates[row.parameter] <= row.q975 for row in reference]
        report["all_estimates_in_published_intervals"] = all(inside)
        report["estimates_outside_published_intervals"] = String.([reference[i].parameter for i in eachindex(inside) if !inside[i]])
        open(joinpath(output,"parameter_comparison.csv"),"w") do io
            println(io,"parameter,estimate,published_mean,published_q025,published_q975,inside_published_interval,rhat,ess_bulk,ess_tail")
            for (i,row) in enumerate(reference)
                diagnostics = method == :bayes ? (result.summary[row.parameter].rhat,
                    result.summary[row.parameter].ess_bulk,result.summary[row.parameter].ess_tail) : ("","","")
                println(io,join((row.parameter,estimates[row.parameter],row.mean,row.q025,row.q975,inside[i],diagnostics...),','))
            end
        end
    end
    report["reproduction_accepted"] = report["converged"] && get(report,"all_estimates_in_published_intervals",false)
    open(joinpath(output,"report.toml"),"w") do io
        TOML.print(io,report;sorted=true)
    end
    println("Convergence diagnostics passed: ", report["converged"], "; reproduction accepted: ", report["reproduction_accepted"])
    return (; result, report)
end

if abspath(PROGRAM_FILE) == @__FILE__
    options = Dict{String,String}()
    for argument in ARGS
        pieces = split(replace(argument,r"^--"=>""),'=';limit=2)
        length(pieces) == 2 || error("Use --name=value arguments")
        options[pieces[1]] = pieces[2]
    end
    allowed = Set(["trial","method","output","samples","warmup","chains","starts","maxiters","outer_iterations","seed","target_accept","max_depth","abstol","reltol","solver_maxiters","maxtime","solver","chain_execution","initial_file"])
    all(in(allowed),keys(options)) || error("Unknown fit option")
    haskey(options,"output") || error("Provide --output=<new-directory> to preserve previous runs")
    trial = Symbol(get(options,"trial","camu"))
    kwargs = Dict{Symbol,Any}(:method=>Symbol(get(options,"method","bayes")),:output=>options["output"])
    haskey(options,"solver") && (kwargs[:solver] = Symbol(options["solver"]))
    haskey(options,"chain_execution") && (kwargs[:chain_execution] = Symbol(options["chain_execution"]))
    haskey(options,"initial_file") && (kwargs[:initial_file] = options["initial_file"])
    for name in ("samples","warmup","chains","starts","maxiters","outer_iterations","seed","max_depth","solver_maxiters")
        haskey(options,name) && (kwargs[Symbol(name)]=parse(Int,options[name]))
    end
    for name in ("target_accept","abstol","reltol","maxtime")
        haskey(options,name) && (kwargs[Symbol(name)]=parse(Float64,options[name]))
    end
    if trial == :all
        haskey(kwargs,:initial_file) && error("initial_file is trial-specific; choose one trial")
        ispath(kwargs[:output]) && error("Batch output directory already exists")
        for dataset in TRIALS
            trial_options = merge(kwargs, Dict(:output=>joinpath(kwargs[:output],String(dataset))))
            run_cocoa_fit(dataset;trial_options...)
        end
    else
        run_cocoa_fit(trial;kwargs...)
    end
end
