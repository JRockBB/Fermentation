using OrdinaryDiffEq
using ForwardDiff

@testset "Output-interval solver budgets" begin
    integrate = FermentModels.Cocoa._solve_output_intervals
    problem = ODEProblem((u,p,t)->-p*u,[1.],(0.,1.),1.)
    times = [.25,.5,.75,1.]
    events = []
    solution = integrate(problem,Tsit5(),times;maxiters=1,adaptive=false,dt=.25,
        on_interval=e->push!(events,e))
    @test !isnothing(solution)
    @test getproperty.(events,:accepted) == ones(Int,4)
    @test getproperty.(events,:attempted) == ones(Int,4)
    @test all(e->e.failure==:none,events)
    reference = solve(problem,Tsit5();saveat=times,save_start=false,dense=false,
        adaptive=false,dt=.25,maxiters=4)
    @test Array(solution) == Array(reference)
    @test !successful_retcode(solve(problem,Tsit5();adaptive=false,dt=.25,maxiters=1,
        verbose=DEVerbosity(SciMLLogging.None())))
    # One accepted step can interpolate several outputs; no restart or forced stop.
    empty!(events)
    dense_times = [.1,.2,.3,.4,.5,1.]
    solution = integrate(problem,Tsit5(),dense_times;maxiters=1,adaptive=false,dt=.5,
        on_interval=e->push!(events,e))
    @test getproperty.(events,:accepted) == [1,0,0,0,0,1]
    @test solution.t == dense_times
    @test Array(solution) == Array(solve(problem,Tsit5();saveat=dense_times,
        save_start=false,dense=false,adaptive=false,dt=.5))
    # Exhaustion in a later interval rejects the entire trajectory.
    empty!(events)
    @test isnothing(integrate(problem,Tsit5(),[.25,1.];maxiters=1,adaptive=false,dt=.25,
        on_interval=e->push!(events,e)))
    @test last(events).interval == 2 && last(events).failure == :output_step_limit
    @test last(events).accepted == 1
    # Rejected attempts are distinct from accepted steps, and have a finite guard.
    empty!(events)
    @test isnothing(integrate(problem,Tsit5(),[.001,1.];maxiters=1,dt=1.,
        abstol=1e-12,reltol=1e-12,on_interval=e->push!(events,e)))
    @test first(events).failure == :none && first(events).accepted == 1
    @test first(events).attempted > 1
    empty!(events)
    @test isnothing(integrate(problem,Tsit5(),[1.];maxiters=100,max_step_attempts=1,
        dt=1.,abstol=1e-12,reltol=1e-12,verbose=DEVerbosity(SciMLLogging.None()),on_interval=e->push!(events,e)))
    @test only(events).failure == :step_attempt_limit
    @test only(events).attempted == 1
    # Failures other than the budget remain failures, and coding errors propagate.
    @test isnothing(integrate(problem,Tsit5(),[1.];maxiters=100,
        unstable_check=(dt,u,p,t)->true,verbose=DEVerbosity(SciMLLogging.None())))
    broken = ODEProblem((u,p,t)->error("broken RHS"),[1.],(0.,1.))
    @test_throws ErrorException integrate(broken,Tsit5(),[1.];maxiters=100)
    for algorithm in (Tsit5(),Rodas5P())
        prediction(p) = Array(integrate(ODEProblem((u,p,t)->-p*u,[1.],(0.,1.),p),
            algorithm,times;maxiters=1000,abstol=1e-11,reltol=1e-11))
        @test vec(prediction(1.)) ≈ exp.(-times) rtol=1e-8
        @test vec(ForwardDiff.derivative(prediction,1.)) ≈ -times.*exp.(-times) rtol=1e-7
    end
end

@testset "Cocoa budget mode trajectory and gradient invariance" begin
    # Synthetic positive coefficients exercise all eight states; not estimates.
    names = (COCOA_PARAMETERS...,:sigma)
    initial = NamedTuple{names}(ntuple(_->.5,25))
    lower,upper = NamedTuple{names}(ntuple(_->1e-5,25)),NamedTuple{names}(ntuple(_->5.,25))
    data = (initial=ones(8),t0=0.,times=[.1,1.,4.],values=zeros(24),sigma=:sigma)
    for solver in (:tsit5,:rodas5p)
        a = cocoa_fit_model(;initial,lower,upper,solver,abstol=1e-9,reltol=1e-9)
        b = cocoa_fit_model(;initial,lower,upper,solver,abstol=1e-9,reltol=1e-9,
            maxiters_scope=:output_interval)
        @test a.predict(initial,data) == b.predict(initial,data)
        z = log.(collect(values(initial)))
        loss(model,z) = gaussian_nll(model,data,NamedTuple{names}(Tuple(exp.(z))))
        @test ForwardDiff.gradient(z->loss(a,z),z) == ForwardDiff.gradient(z->loss(b,z),z)
        limited = cocoa_fit_model(;initial,lower,upper,solver,maxiters_scope=:output_interval,maxiters=1)
        @test isnothing(limited.predict(initial,data))
        @test_throws FermentModels.CoreModels._PredictionFailure gaussian_nll(limited,data,initial)
    end
    for options in ((maxiters_scope=:unknown,),(maxiters=true,),(max_step_attempts=0,),
            (max_step_attempts=true,),(maxiters=big(typemax(Int))+1,))
        @test_throws ArgumentError cocoa_fit_model(;initial,lower,upper,options...)
    end
end
