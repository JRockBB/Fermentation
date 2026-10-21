module CocoaBenchmarkBayes

using FermentModels
using Random
using Statistics
using TOML
include("fit_cocoa_benchmark.jl")
include("cocoa_fit_output.jl")
const F = CocoaBenchmarkFits
const B = F.CocoaBenchmarkData

"""Same distributions as cocoa_priors(), with the independently inspected 2022 source."""
function benchmark_priors(; hidden_lactate=false, root=B.ROOT)
    path = joinpath(B.RAW,B.ARCHIVE_DIR,"Stan_models",
        hidden_lactate ? "camu2008/mmc0.stan" : "mm0.stan")
    source_hash = B.digest(joinpath(root,path))
    original = cocoa_priors()
    return NamedTuple{keys(original)}(map(keys(original)) do name
        distribution = original[name].distribution
        provenance = (source="doi:10.1098/rsos.210274",
            location="baseline Stan model block; positive parameter declarations",
            file=replace(path,'\\'=>'/'),sha256=source_hash,
            coordinates="trial-normalized states; time in hours",parameter=String(name),
            distribution=name == :sigma ? "half-Cauchy(0,1)" : "positive-truncated Normal(0.5,0.3)",
            support="positive; unbounded above")
        (;distribution,provenance)
    end)
end

"""Normalized pointwise Gaussian log likelihood; never drops a failed draw."""
function pointwise_loglikelihood(model, data, parameters)
    data.sigma isa Symbol && hasproperty(parameters,data.sigma) ||
        throw(ArgumentError("requires a named shared noise parameter"))
    sigma = parameters[data.sigma]
    sigma isa Real && !(sigma isa Bool) && isfinite(sigma) && sigma > 0 ||
        throw(ArgumentError("noise scale must be finite and positive"))
    indices = findall(!ismissing,data.values)
    isempty(indices) && throw(ArgumentError("no observed rows"))
    prediction = model.predict(parameters,data)
    isnothing(prediction) && error("Prediction failed while exporting log likelihood; no draws may be omitted")
    length(prediction) == length(data.values) || throw(DimensionMismatch("prediction length differs"))
    all(i->isfinite(prediction[i]) && isfinite(data.values[i]),indices) ||
        error("Nonfinite value at an observed row")
    normalizer = log(sigma)+log(2pi)/2
    values = [-normalizer-((data.values[i]-prediction[i])/sigma)^2/2 for i in indices]
    all(isfinite,values) || error("Nonfinite log likelihood; no draws may be omitted")
    return (;indices,values)
end

function observation_rows(trial)
    ntime = length(trial.data.times)
    indices = findall(!ismissing,trial.data.values)
    return [("obs_$i",i,trial.data.times[mod1(i,ntime)],COCOA_STATES[cld(i,ntime)],
        trial.data.values[i]) for i in indices]
end

"""Write all retained joint draws, retaining chain identity and an explicit observation map."""
function write_log_likelihood(output,model,trial,result)
    B.write_csv(joinpath(output,"likelihood_observations.csv"),
        ("column","flat_index","time_h","variable","observed_normalized"),observation_rows(trial))
    names = keys(model.initial)
    keys(result.draws) == names || throw(ArgumentError("posterior parameter order differs"))
    dimensions = size(first(values(result.draws)))
    length(dimensions) == 2 && all(size(v)==dimensions for v in values(result.draws)) ||
        throw(DimensionMismatch("draws must share samples-by-chains dimensions"))
    samples, chains = dimensions
    samples > 0 && chains > 0 || throw(ArgumentError("posterior contains no draws"))
    indices = findall(!ismissing,trial.data.values)
    open(joinpath(output,"log_likelihood.csv"),"w") do io
        println(io,join(("chain","draw",("obs_$i" for i in indices)...),','))
        for chain in 1:chains, draw in 1:samples
            parameters = NamedTuple{names}(Tuple(result.draws[n][draw,chain] for n in names))
            likelihood = pointwise_loglikelihood(model,trial.data,parameters)
            println(io,join((chain,draw,likelihood.values...),','))
        end
    end
    return samples*chains
end

"""Save normalized bands for every state, physical bands only where a scale exists."""
function write_predictive(output,trial,result)
    predictive = result.predictive
    ntime = length(trial.data.times)
    rows = map(eachindex(trial.data.values)) do i
        state_index, time_index = cld(i,ntime), mod1(i,ntime)
        scale = trial.data.scales[state_index]
        latent_mean = mean(predictive.latent_samples[i,:])
        normalized = (latent_mean,predictive.latent.median[i],predictive.latent.lower[i],
            predictive.latent.upper[i],predictive.observed.lower[i],predictive.observed.upper[i])
        physical = map(v->ismissing(scale) ? missing : v*scale,normalized)
        (trial.data.times[time_index],COCOA_STATES[state_index],trial.data.values[i],
            normalized...,scale,trial.observations[time_index,state_index],physical...)
    end
    B.write_csv(joinpath(output,"predictive.csv"),
        ("time_h","variable","observed_normalized","latent_mean_normalized","latent_median_normalized",
         "latent_lower_normalized","latent_upper_normalized","observation_lower_normalized",
         "observation_upper_normalized","source_scale_mg_g","observed_mg_g","latent_mean_mg_g",
         "latent_median_mg_g","latent_lower_mg_g","latent_upper_mg_g","observation_lower_mg_g",
         "observation_upper_mg_g"),rows)
    samples = size(first(values(result.draws)),1)
    B.write_csv(joinpath(output,"predictive_selection.csv"),("column","chain","draw"),
        [(column,cld(index,samples),mod1(index,samples)) for
            (column,index) in enumerate(predictive.draw_indices)])
end

"""
Sequential baseline Bayesian benchmark with independent chain RNGs, retained-draw
checkpoints, and likelihood exports. Defaults use the archived sample counts and
ODE caps. Julia/Turing initialization and integrators differ from the source.
Output must be new. Sampling diagnostics are necessary checks, not proof of
convergence or reproduction. Pointwise likelihoods are not a computed PSIS-LOO score.
"""
function run_benchmark_bayes(;output,scope=:benchmark,trials=nothing,samples=2000,
        warmup=1000,chains=4,target_accept=0.995,max_depth=10,seed=141085,
        solver=:tsit5,abstol=1e-6,reltol=1e-6,solver_maxiters=:source,
        predictive_draws=500,chain_execution=:serial)
    ispath(output) && throw(ArgumentError("Use a new output directory to preserve previous runs"))
    for (name,value,minimum) in ((:samples,samples,10),(:warmup,warmup,0),(:chains,chains,1),
            (:max_depth,max_depth,1),(:seed,seed,0),(:predictive_draws,predictive_draws,1))
        value isa Integer && !(value isa Bool) && value >= minimum ||
            throw(ArgumentError("$name must be an integer at least $minimum"))
    end
    target_accept isa Real && !(target_accept isa Bool) && isfinite(target_accept) &&
        0 < target_accept < 1 || throw(ArgumentError("target_accept must be in (0,1)"))
    chain_execution in (:serial,:threads) || throw(ArgumentError("choose serial or threads"))
    solver_maxiters == :source || (solver_maxiters isa Integer && !(solver_maxiters isa Bool) &&
        solver_maxiters > 0) || throw(ArgumentError("solver_maxiters must be :source or a positive integer"))
    available = B.benchmark_trials(;scope)
    selected = F.select_trials(available,trials)
    seeds = Dict(t.trial=>seed+i-1 for (i,t) in enumerate(available))
    caps = Dict(t.trial=>(solver_maxiters == :source ? (t.hidden_lactate ? 1_000_000 : 10_000) :
        solver_maxiters) for t in selected)
    names = (COCOA_PARAMETERS..., :sigma)
    # Required only by FitModel's constructor; these MLE bounds are not applied to NUTS.
    initial = NamedTuple{names}(ntuple(_->0.5,25))
    lower = NamedTuple{names}(ntuple(_->1e-5,25))
    upper = NamedTuple{names}(ntuple(_->5.,25))
    models = Dict(t.trial=>cocoa_fit_model(;initial,lower,upper,solver,abstol,reltol,
        maxiters=caps[t.trial]) for t in selected)
    implementation = ("scripts/fit_cocoa_benchmark_bayes.jl","scripts/fit_cocoa_benchmark.jl",
        "scripts/cocoa_fit_output.jl","src/core/bayes.jl",F.LIKELIHOOD_IMPLEMENTATION...)
    settings = Dict("method"=>"bayes","scope"=>string(scope),"trials"=>String.([t.trial for t in selected]),
        "samples_per_chain"=>samples,"warmup_per_chain"=>warmup,"chains"=>chains,
        "target_accept"=>target_accept,"max_depth"=>max_depth,"seed"=>seed,
        "trial_seeds"=>Dict(string(t.trial)=>seeds[t.trial] for t in selected),
        "solver"=>string(solver),"abstol"=>abstol,"reltol"=>reltol,
        "solver_maxiters"=>string(solver_maxiters),
        "trial_solver_maxiters"=>Dict(string(t.trial)=>caps[t.trial] for t in selected),
        "predictive_draws"=>min(predictive_draws,samples*chains),
        "chain_execution"=>string(chain_execution),"julia_threads"=>Threads.nthreads(),
        "julia_version"=>string(VERSION),"time_unit"=>"h","source_concentration_unit"=>"mg/g pulp",
        "parameter_coordinates"=>"trial-normalized ODE coefficients and normalized observation sigma",
        "bayesian_support"=>"positive and unbounded above; FitModel MLE bounds are not applied",
        "initialization"=>"independent prior draws; differs from source absolute-Normal(0.5,0.2) initialization",
        "initial_step_size"=>"automatic Turing selection; source specifies 0.1",
        "source_sampler"=>"RStan; 4 chains, 3000 total iterations with 1/3 warmup; adapt_delta 0.995-0.999; depth 10; seed 141085",
        "source_integrator"=>"Stan RK45; atol=rtol=1e-6; mm0 cap 10000; mmc0 cap 1000000",
        "source"=>"doi:10.1098/rsos.210274; doi:10.5281/zenodo.5510370",
        "likelihood_unit"=>"one observed state-time cell, normalized using fixed full-trial scales",
        "likelihood_export"=>"all retained joint draws; missing observations omitted; initial conditions excluded",
        "checkpoint_limitation"=>"retained draws only; no adaptation state, so checkpoints cannot resume NUTS",
        "implementation_sha256"=>Dict(p=>B.digest(joinpath(B.ROOT,p)) for p in implementation))
    mkpath(output)
    F.write_toml(joinpath(output,"settings.toml"),settings)
    references = Dict(r.trial=>r for r in B.published_baseline())
    comparison = Tuple[]
    for (number,trial) in enumerate(selected)
        directory = joinpath(output,string(trial.trial))
        mkpath(directory)
        F.write_toml(joinpath(directory,"input.toml"),Dict("trial"=>string(trial.trial),
            "source_file"=>trial.source_file,"source_sha256"=>trial.source_sha256,
            "source_branch"=>trial.branch,"hidden_lactate"=>trial.hidden_lactate,
            "initial_normalized"=>trial.data.initial,"times_h"=>trial.data.times,
            "state_order"=>String.(collect(COCOA_STATES)),"seed"=>seeds[trial.trial],
            "solver_maxiters"=>caps[trial.trial]))
        priors = benchmark_priors(;hidden_lactate=trial.hidden_lactate)
        F.write_toml(joinpath(directory,"priors.toml"),Dict(string(n)=>
            Dict(string(k)=>v for (k,v) in pairs(priors[n].provenance)) for n in names))
        for chain in 1:chains
            open(joinpath(directory,"chain_$(lpad(chain,3,'0')).csv"),"w") do io
                write_cocoa_draw_header(io,names)
            end
        end
        function checkpoint(chain,draw,parameters,diagnostics)
            open(joinpath(directory,"chain_$(lpad(chain,3,'0')).csv"),"a") do io
                write_cocoa_draw(io,chain,draw,parameters,diagnostics,names)
            end
            if draw == 1 || draw % 100 == 0 || draw == samples
                println(trial.trial," chain ",chain,": saved ",draw,"/",samples," retained draws")
                flush(stdout)
            end
        end
        report = Dict{String,Any}("trial"=>string(trial.trial),"status"=>"sampling",
            "diagnostics_passed"=>false,"direct_metric_comparison_available"=>false,
            "psis_loo"=>"not computed; pointwise log likelihood exported separately",
            "published_psis_loo"=>"missing: archived posterior outputs absent")
        F.write_toml(joinpath(directory,"report.toml"),report)
        println("Sampling ",number,"/",length(selected)," ",trial.trial); flush(stdout)
        started = time()
        try
            result = fit_bayes(models[trial.trial],trial.data,priors;samples,warmup,chains,
                target_accept,max_depth,predictive_draws,rng=MersenneTwister(seeds[trial.trial]),
                chain_execution,on_draw=checkpoint)
            merge!(report,Dict("status"=>"postprocessing","diagnostics_passed"=>result.diagnostics.passed,
                "divergences"=>result.diagnostics.divergences,"depth_hits"=>result.diagnostics.depth_hits,
                "max_rhat"=>maximum(r.rhat for r in values(result.summary)),
                "min_ess_bulk"=>minimum(r.ess_bulk for r in values(result.summary)),
                "min_ess_tail"=>minimum(r.ess_tail for r in values(result.summary)),
                "rhat_threshold"=>result.diagnostics.rhat_threshold,
                "minimum_ess"=>result.diagnostics.minimum_ess,
                "retained_draws"=>samples*chains,"observations"=>count(!ismissing,trial.data.values)))
            F.write_toml(joinpath(directory,"report.toml"),report)
            B.write_csv(joinpath(directory,"parameters.csv"),
                ("parameter","mean","median","lower_95","upper_95","rhat","ess_bulk","ess_tail"),
                [(n,(getproperty(result.summary[n],k) for k in
                    (:mean,:median,:lower,:upper,:rhat,:ess_bulk,:ess_tail))...) for n in names])
            write_predictive(directory,trial,result)
            report["likelihood_draws"] = write_log_likelihood(directory,models[trial.trial],trial,result)
            report["status"] = result.diagnostics.passed ? "posterior_diagnostics_passed" : "posterior_diagnostics_failed"
            report["elapsed_seconds"] = time()-started
            F.write_toml(joinpath(directory,"report.toml"),report)
        catch err
            report["status"] = "failed"
            report["diagnostics_passed"] = false
            report["error"] = sprint(showerror,err)
            report["elapsed_seconds"] = time()-started
            F.write_toml(joinpath(directory,"report.toml"),report)
            rethrow()
        end
        reference = get(references,trial.trial,nothing)
        push!(comparison,(trial.trial,trial.in_paper_benchmark,report["observations"],report["status"],
            report["max_rhat"],report["min_ess_bulk"],report["min_ess_tail"],
            report["divergences"],report["depth_hits"],report["retained_draws"],report["elapsed_seconds"],
            isnothing(reference) ? false : reference.success,
            (isnothing(reference) ? missing : getproperty(reference,k) for k in (:max_rhat,:min_bulk_ess,:min_tail_ess,:bma_weight))...,
            missing,missing))
        B.write_csv(joinpath(output,"comparison.csv"),
            ("trial","in_paper_benchmark","observations","status","max_rhat","min_ess_bulk","min_ess_tail",
             "divergences","depth_hits","retained_draws","elapsed_seconds","published_baseline_success",
             "published_max_rhat","published_min_bulk_ess","published_min_tail_ess","published_bma_weight",
             "psis_loo","published_psis_loo"),comparison)
        println("Completed ",trial.trial,": ",report["status"]); flush(stdout)
    end
    return comparison
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    options = Dict{Symbol,Any}()
    integers = (:samples,:warmup,:chains,:max_depth,:seed,:predictive_draws)
    floats = (:target_accept,:abstol,:reltol)
    symbols = (:scope,:solver,:chain_execution)
    for argument in ARGS
        parts = split(replace(argument,r"^--"=>""),'=';limit=2)
        length(parts)==2 || error("Use --name=value")
        key = Symbol(parts[1])
        key in (:output,:trials,:solver_maxiters,integers...,floats...,symbols...) ||
            error("Unknown option $key")
        options[key] = key in integers ? parse(Int,parts[2]) :
            key in floats ? parse(Float64,parts[2]) : key in symbols ? Symbol(parts[2]) :
            key == :trials ? Symbol.(split(parts[2],',')) :
            key == :solver_maxiters ? (parts[2]=="source" ? :source : parse(Int,parts[2])) : parts[2]
    end
    haskey(options,:output) || error("Provide --output=<new-directory>")
    CocoaBenchmarkBayes.run_benchmark_bayes(;options...)
end
