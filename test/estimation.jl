using Random
using LinearAlgebra

@testset "Bounded multistart Gaussian estimation" begin
    predict = (p,d) -> p.intercept .+ p.slope .* d.times
    model = FitModel(predict; initial=(intercept=1., slope=1.),
        lower=(slope=0.05, intercept=0.05), upper=(intercept=20., slope=10.))
    times = collect(range(0.,4.; length=30))
    sigma = 0.15 .+ 0.03 .* times
    truth = (intercept=2., slope=1.5)
    noise = randn(MersenneTwister(20261019), length(times))
    data = (; times, values=predict(truth, (; times)) + sigma .* noise, sigma)
    result = fit_mle(model, data; starts=4, rng=MersenneTwister(82))
    @test result.converged
    @test length(result.runs) == 4
    @test result.runs[1].initial == model.initial
    @test result.objective == minimum(r.objective for r in result.runs if r.converged)
    @test all(r -> !r.converged || (0.05 <= r.parameters.intercept <= 20 && 0.05 <= r.parameters.slope <= 10), result.runs)
    @test result.uncertainty.status == :ok
    @test all(result.uncertainty.intervals[name].lower < truth[name] < result.uncertainty.intervals[name].upper for name in keys(truth))

    # Independent weighted least-squares solution and covariance, in physical units.
    design = hcat(ones(length(times)), times) ./ sigma
    reference = design \ (data.values ./ sigma)
    @test collect(values(result.parameters)) ≈ reference rtol=1e-6
    transform = Diagonal(collect(values(result.parameters)))
    @test transform * result.uncertainty.covariance_log * transform ≈ inv(design' * design) rtol=1e-6
    @test result.residuals ≈ data.values-result.predictions
    @test result.standardized_residuals ≈ result.residuals ./ sigma
    @test result.objective ≈ sum(abs2, result.standardized_residuals)/2 + sum(log, sigma) + length(times)*log(2pi)/2
    @test result.objective ≈ gaussian_nll(model, data, result.parameters)
    repeated = fit_mle(model, data; starts=4, rng=MersenneTwister(82))
    @test getproperty.(result.runs, :initial) == getproperty.(repeated.runs, :initial)
    @test collect(values(repeated.parameters)) ≈ collect(values(result.parameters))

    missing_values = Union{Missing,Float64}[data.values...]
    missing_values[5] = missing
    incomplete = (; times, values=missing_values, sigma)
    explicit = [(slope=0.6, intercept=0.8), (intercept=3., slope=3.)]
    masked = fit_mle(model, incomplete; starts=explicit)
    @test masked.converged
    @test ismissing(masked.residuals[5]) && ismissing(masked.standardized_residuals[5])
    @test masked.observed_indices == setdiff(1:30, [5])
    reduced = fit_mle(model, (; times=times[masked.observed_indices],
        values=data.values[masked.observed_indices], sigma=sigma[masked.observed_indices]); starts=1)
    @test collect(values(masked.parameters)) ≈ collect(values(reduced.parameters)) rtol=1e-5

    constant = FitModel((p,d)->fill(p.level, length(d.values)); initial=(level=0.5,),
                       lower=(level=0.1,), upper=(level=1.,))
    boundary = fit_mle(constant, (values=fill(2., 5), sigma=0.1); starts=1)
    @test boundary.converged
    @test boundary.parameters.level ≈ 1. atol=1e-5
    @test boundary.uncertainty.status == :active_bounds
    @test isnothing(boundary.uncertainty.intervals)
    confounded = FitModel((p,d)->fill(p.a*p.b, length(d.values)); initial=(a=1., b=1.),
                         lower=(a=0.01,b=0.01), upper=(a=100.,b=100.))
    ambiguous = fit_mle(confounded, (values=fill(2., 10), sigma=0.1); starts=1)
    @test ambiguous.converged
    @test ambiguous.uncertainty.status == :rank_deficient
    @test isnothing(ambiguous.uncertainty.intervals)
    weak = fit_mle(constant, (values=[0.5], sigma=100.); starts=1)
    @test weak.uncertainty.status == :interval_crosses_bounds
    @test isnothing(weak.uncertainty.intervals)

    invalid = FitModel((p,d)->nothing; initial=(a=1.,), lower=(a=0.1,), upper=(a=10.,))
    failed = fit_mle(invalid, (values=[1.], sigma=0.1); starts=2, rng=MersenneTwister(1))
    @test !failed.converged && isnothing(failed.parameters) && isnothing(failed.best_start)
    @test length(failed.runs) == 2 && all(r -> r.retcode == "PredictionFailure", failed.runs)
    partial = FitModel((p,d)->p.a > 5 ? nothing : fill(p.a,length(d.values));
        initial=(a=1.,), lower=(a=0.1,), upper=(a=10.,))
    resumed = fit_mle(partial, (values=fill(2.,5),sigma=0.1); starts=[(a=8.,),(a=1.,)])
    @test resumed.converged && resumed.best_start == 2
    @test resumed.runs[1].retcode == "PredictionFailure"
    bug = FitModel((p,d)->throw(ArgumentError("broken model")); initial=(a=1.,), lower=(a=0.1,), upper=(a=10.,))
    @test_throws ArgumentError fit_mle(bug, (values=[1.], sigma=1.); starts=1)
    short = fit_mle(model, data; starts=[(intercept=10.,slope=5.)], maxiters=1, outer_iterations=1)
    @test !short.converged && isnothing(short.parameters)
    @test !short.runs[1].converged

    @test_throws DomainError FitModel(predict; initial=(a=1.,), lower=(a=0.,), upper=(a=2.,))
    @test_throws ArgumentError FitModel(predict; initial=(a=1.,), lower=(b=0.1,), upper=(a=2.,))
    @test_throws DomainError fit_mle(model, (; times, values=data.values, sigma=0.); starts=1)
    @test_throws ArgumentError fit_mle(model, (values=[missing], sigma=1.); starts=1)
    @test_throws DimensionMismatch fit_mle(model, (; times, values=data.values, sigma=[1.]); starts=1)
    @test_throws ArgumentError fit_mle(model, data; starts=0)
    @test_throws DomainError fit_mle(model, data; starts=[(intercept=30., slope=1.)])
end

@testset "Noisy Monod ODE parameter recovery" begin
    spec = ModelSpec(; organisms=(:microbe,), substrates=(:sugar,), products=(:acid,),
                     stoichiometry=reshape([1.,-2.,1.],3,1))
    function predict_monod(parameters, data)
        initial = [0.2, 4., 0.] .* one(parameters.mu + parameters.Ks)
        problem = fermentation_problem(spec, initial, (0.,8.),
            (mu_max=[parameters.mu], Ks=reshape([parameters.Ks],1,1)))
        trajectory = solve(problem, Tsit5(); abstol=1e-10, reltol=1e-9)
        successful_retcode(trajectory) || return nothing
        return [trajectory(t)[1] for t in data.times]
    end
    truth = (mu=0.8, Ks=0.7)
    times = collect(0.25:0.25:7.)
    clean = predict_monod(truth, (; times))
    sigma = 0.02
    data = (; times, values=clean + sigma .* randn(MersenneTwister(314), length(clean)), sigma)
    model = FitModel(predict_monod; initial=(mu=0.5,Ks=1.),
                     lower=(mu=0.1,Ks=0.05), upper=(mu=2.,Ks=5.))
    result = fit_mle(model, data; starts=3, rng=MersenneTwister(12))
    @test result.converged
    @test result.uncertainty.status == :ok
    @test all(result.uncertainty.intervals[name].lower < truth[name] < result.uncertainty.intervals[name].upper for name in keys(truth))
    @test maximum(abs.(result.predictions-clean)) < 2sigma
    @test result.parameters.mu ≈ truth.mu rtol=0.15
    @test result.parameters.Ks ≈ truth.Ks rtol=0.3
end
