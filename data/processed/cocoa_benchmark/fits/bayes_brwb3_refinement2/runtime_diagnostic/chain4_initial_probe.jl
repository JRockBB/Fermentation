using FermentModels, Random, Turing, Optimization
include("../scripts/fit_cocoa_benchmark_bayes.jl")
const C = CocoaBenchmarkBayes
const K = FermentModels.CoreModels

# A bounded initialization diagnostic, not a fitted posterior or resumed chain.
function probe()
    trial = only(filter(t->t.trial==:brwb3,C.B.benchmark_trials(;scope=:all)))
    names = (COCOA_PARAMETERS...,:sigma)
    initial = NamedTuple{names}(ntuple(_->.5,25))
    lower = NamedTuple{names}(ntuple(_->1e-5,25))
    upper = NamedTuple{names}(ntuple(_->5.,25))
    model = cocoa_fit_model(;initial,lower,upper,solver=:tsit5,abstol=1e-8,reltol=1e-8,maxiters=10_000)
    observations = K._observations(trial.data)
    joint,_ = K._bayes_priors(model,C.benchmark_priors(;hidden_lactate=false))
    posterior = K._gaussian_posterior(model,trial.data,observations.indices,observations.y,observations,joint,names)
    parent_rng = MersenneTwister(141100)
    rng = copy(parent_rng)
    seed = rand(parent_rng,UInt,4)[4]
    Random.seed!(rng,seed)
    sampler = Turing.NUTS(1000,.995;max_depth=12,adtype=Optimization.AutoForwardDiff())
    println("Diagnostic only; chain=4; seed=",seed,"; requested warmup steps=10; adaptation schedule=1000")
    flush(stdout)
    started = time()
    try
        _,state = Turing.AbstractMCMC.step(rng,posterior,sampler;
            initial_params=Turing.InitFromPrior(),nadapts=1000,discard_sample=true,verbose=false)
        println("warmup_step=1; elapsed_seconds=",time()-started); flush(stdout)
        for step in 2:10
            _,state = Turing.AbstractMCMC.step(rng,posterior,sampler,state;
                nadapts=1000,discard_sample=true,verbose=false)
            println("warmup_step=",step,"; elapsed_seconds=",time()-started); flush(stdout)
        end
        println("Bounded initialization diagnostic completed; no posterior acceptance claimed.")
    catch err
        println("exception_class=",nameof(typeof(err))); flush(stdout)
        if hasproperty(err,:msg) && getproperty(err,:msg) isa AbstractString
            println("message=",first(getproperty(err,:msg),1000)); flush(stdout)
        end
        for frame in first(stacktrace(catch_backtrace()),24)
            println(frame.func," at ",frame.file,":",frame.line)
        end
        flush(stdout)
        exit(1)
    end
end
probe()
