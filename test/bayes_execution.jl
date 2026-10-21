using Random
using Distributions: LogNormal
using Turing: @varname
using FlexiChains: Extra

@testset "Callback diagnostics stay local while another chain waits" begin
    entered,release,starting = Channel{Nothing}(1),Channel{Nothing}(1),Channel{Nothing}(1)
    kept,events = [],[]
    function checkpoint(chain,draw,p,d)
        push!(kept,(;chain,draw,p,d))
        if chain == 1
            put!(entered,nothing)
            take!(release) # models checkpoint I/O yielding while holding the lock
        end
    end
    callback = FermentModels.CoreModels._bayes_callback((:rate,),checkpoint,
        (chain,iteration,phase,d)->push!(events,(;chain,iteration,phase,d)),1)
    first_transition = (params=Dict(@varname(theta)=>[.8]),stats=(numerical_error=false,tree_depth=2))
    second_transition = (params=Dict(@varname(theta)=>[1.2]),stats=(numerical_error=true,tree_depth=7))
    task1 = @async callback(nothing,nothing,nothing,first_transition,nothing,2;chain_number=1)
    take!(entered)
    task2 = @async begin
        put!(starting,nothing)
        callback(nothing,nothing,nothing,second_transition,nothing,2;chain_number=2)
    end
    take!(starting) # task2 runs into the locked callback before task1 is released
    put!(release,nothing)
    wait(task1); wait(task2)
    @test getproperty.(kept,:chain) == getproperty.(events,:chain) == [1,2]
    @test getproperty.(kept,:d) == [(divergent=e.d.divergent,tree_depth=e.d.tree_depth) for e in events] ==
        [(divergent=false,tree_depth=2),(divergent=true,tree_depth=7)]
    @test getproperty.(getproperty.(kept,:p),:rate) == [.8,1.2]
    @test all(e->e.phase==:sampling,events)
end

@testset "Bayesian chain execution and retained-draw checkpoints" begin
    model = FitModel((p,d)->fill(log(p.rate),length(d.values));initial=(rate=1.,),
        lower=(rate=0.1,),upper=(rate=10.,))
    data = (values=[0.2,0.3,0.4],sigma=0.2)
    prior = (rate=LogNormal(0.,0.7),)
    checkpoints = []
    serial = fit_bayes(model,data,prior;samples=40,warmup=30,chains=2,
        predictive_draws=20,rng=MersenneTwister(1403),
        on_draw=(chain,draw,p,d)->push!(checkpoints,(;chain,draw,p,d)))
    @test length(checkpoints) == 80 # all warmup excluded
    @test all(row.p.rate == serial.draws.rate[row.draw,row.chain] for row in checkpoints)
    @test [row.draw for row in checkpoints if row.chain == 1] == collect(1:40)
    @test count(row->row.d.divergent,checkpoints) == serial.diagnostics.divergences
    @test count(row->row.d.tree_depth >= 10,checkpoints) == serial.diagnostics.depth_hits
    threaded_checkpoints = []
    threaded = fit_bayes(model,data,prior;samples=40,warmup=30,chains=2,
        predictive_draws=20,rng=MersenneTwister(1403),chain_execution=:threads,
        on_draw=(chain,draw,p,d)->push!(threaded_checkpoints,(;chain,draw,p,d)))
    @test threaded.settings.chain_execution == :threads
    @test threaded.settings.julia_threads == Threads.nthreads()
    @test threaded.draws == serial.draws # chain RNG streams do not depend on scheduling
    @test threaded.diagnostics == serial.diagnostics
    @test length(threaded_checkpoints) == 80
    @test all(row.p.rate == threaded.draws.rate[row.draw,row.chain] for row in threaded_checkpoints)
    @test all([row.draw for row in threaded_checkpoints if row.chain == c] == collect(1:40) for c in 1:2)
    for execution in (:serial,:threads)
        reference = execution == :serial ? serial : threaded
        events,kept = [],[]
        monitored = fit_bayes(model,data,prior;samples=40,warmup=30,chains=2,
            predictive_draws=20,rng=MersenneTwister(1403),chain_execution=execution,
            on_draw=(chain,draw,p,d)->push!(kept,(;chain,draw,p,d)),
            on_iteration=(chain,iteration,phase,d)->push!(events,(;chain,iteration,phase,d)))
        @test monitored.draws == serial.draws
        @test monitored.diagnostics == serial.diagnostics
        @test monitored.predictive.latent_samples == reference.predictive.latent_samples
        @test monitored.predictive.observation_samples == reference.predictive.observation_samples
        @test length(kept) == 80 && length(events) == 140
        @test all(row.p.rate == monitored.draws.rate[row.draw,row.chain] for row in kept)
        for chain in 1:2
            rows = filter(e->e.chain==chain,events)
            @test getproperty.(rows,:iteration) == collect(1:70)
            @test first(rows).phase == :initialization && ismissing(first(rows).d.divergent)
            @test all(e->e.phase==:warmup,rows[2:30])
            @test all(e->e.phase==:sampling,rows[31:end])
            retained = filter(r->r.chain==chain,kept)
            @test getproperty.(retained,:d) ==
                [(divergent=e.d.divergent,tree_depth=e.d.tree_depth) for e in rows[31:end]]
            @test all(r.d == (divergent=monitored.chain[Extra(:numerical_error)][r.draw,chain],
                tree_depth=monitored.chain[Extra(:tree_depth)][r.draw,chain]) for r in retained)
            for key in (:n_steps,:step_size,:acceptance_rate,:hamiltonian_energy,
                    :hamiltonian_energy_error,:max_hamiltonian_energy_error)
                @test ismissing(getproperty(first(rows).d,key))
                @test all(isequal(getproperty(e.d,key),monitored.chain[Extra(key)][e.iteration-30,chain])
                    for e in rows[31:end])
            end
        end
    end
    events = []
    no_warmup = fit_bayes(model,data,prior;samples=10,warmup=0,chains=1,
        rng=MersenneTwister(1403),on_iteration=(args...)->push!(events,args),on_draw=(args...)->nothing)
    @test length(events) == 10 && size(no_warmup.draws.rate) == (10,1)
    @test_throws ErrorException fit_bayes(model,data,prior;samples=10,warmup=5,chains=1,
        on_iteration=(args...)->error("runtime log unavailable"))
    @test_throws ArgumentError fit_bayes(model,data,prior;chain_execution=:invalid)
    @test_throws ErrorException fit_bayes(model,data,prior;samples=10,warmup=0,chains=1,
        on_draw=(args...)->error("checkpoint unavailable"))
end
