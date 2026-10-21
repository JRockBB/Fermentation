@testset "Optimization time budgets preserve non-convergence" begin
    # A delayed nonlinear prediction ensures the tiny budget is exercised;
    # assertions concern fit status, not flaky wall-clock upper bounds.
    function delayed_prediction(p,d)
        sleep(0.01)
        return [10(p.b-p.a^2),1-p.a]
    end
    model = FitModel(delayed_prediction;initial=(a=0.2,b=2.),
        lower=(a=0.01,b=0.01),upper=(a=5.,b=5.))
    data = (values=[0.,0.],sigma=1.)
    completed = []
    limited = fit_mle(model,data;starts=1,maxtime=1e-6,
        on_start_complete=(index,run)->push!(completed,(index,run)))
    @test !limited.converged
    @test isnothing(limited.parameters) && limited.uncertainty.status == :no_converged_fit
    @test length(limited.runs) == 1 && !only(limited.runs).converged
    @test only(limited.runs).time_limit_reached === true
    @test only(limited.runs).elapsed_seconds > 0
    @test completed == [(1,only(limited.runs))]
    for invalid in (0.,-1.,Inf,NaN,true)
        @test_throws ArgumentError fit_mle(model,data;starts=1,maxtime=invalid)
    end
end

@testset "Completed-start callbacks retain successes and failures" begin
    model = FitModel((p,d)->p.a > 3 ? nothing : fill(p.a,3);
        initial=(a=1.,),lower=(a=0.01,),upper=(a=5.,))
    data = (values=[1.5,2.,2.5],sigma=0.5)
    completed = []
    result = fit_mle(model,data;starts=[(a=4.,),(a=2.,)],
        on_start_complete=(index,run)->push!(completed,(index,run)))
    @test result.converged && result.best_start == 2
    @test length(completed) == 2 && first.(completed) == [1,2]
    @test isequal(last.(completed),result.runs)
    @test result.runs[1].retcode == "PredictionFailure"
    @test ismissing(result.runs[1].time_limit_reached)
    @test result.runs[2].time_limit_reached === false
    # Checkpoint failures must not silently discard results or proceed to more starts.
    calls = Ref(0)
    @test_throws ErrorException fit_mle(model,data;starts=2,
        on_start_complete=(i,r)->(calls[] += 1; error("checkpoint unavailable")))
    @test calls[] == 1
end
