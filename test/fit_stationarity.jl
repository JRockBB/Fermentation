@testset "MLE stationarity is separate from native termination" begin
    model = FitModel((p,d)->fill(p.a,3);initial=(a=1.,),
        lower=(a=0.1,),upper=(a=3.,))
    # The constrained optimum is at a=3: a large outward gradient is valid.
    result = fit_mle(model,(values=[3.,4.,5.],sigma=1.);starts=1)
    run = only(result.runs)
    @test result.converged && run.optimizer_converged && run.stationarity_passed
    @test run.projected_gradient_inf <= run.gradient_tolerance
    @test result.parameters.a ≈ 3. atol=1e-6

    # Finite precision can stop the optimizer without satisfying an arbitrarily
    # strict gradient tolerance. Such a point must not become an accepted MLE.
    nonlinear = FitModel((p,d)->[p.a,p.a^2];initial=(a=1.,),
        lower=(a=0.1,),upper=(a=3.,))
    strict = fit_mle(nonlinear,(values=[1.,3.],sigma=1.);starts=1,gtol=1e-30)
    stopped = only(strict.runs)
    @test !strict.converged && isnothing(strict.parameters)
    @test !stopped.stationarity_passed
    @test stopped.projected_gradient_inf > stopped.gradient_tolerance
    # Preserve native termination independently, even if the optimizer says success.
    @test stopped.optimizer_converged == (stopped.retcode == "Success")
    @test strict.uncertainty.status == :no_converged_fit

    failing = FitModel((p,d)->nothing;initial=(a=1.,),
        lower=(a=0.1,),upper=(a=3.,))
    failed = only(fit_mle(failing,(values=[1.],sigma=1.);starts=1).runs)
    @test !failed.optimizer_converged && !failed.stationarity_passed
    @test ismissing(failed.projected_gradient_inf)
end
