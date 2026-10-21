using Distributions: LogNormal, Normal
using Statistics

@testset "Bayesian Gaussian inference and predictive coverage" begin
    prior = literature_prior(1.; log_sd=0.7, source="synthetic test", location="test/bayes.jl",
        width_source="synthetic test", width_location="test/bayes.jl")
    @test median(prior.distribution) ≈ 1.
    @test prior.provenance.log_sd == 0.7
    @test_throws ArgumentError literature_prior(1.; log_sd=0.7, source="", location="test",
        width_source="synthetic", width_location="test")
    @test_throws DomainError literature_prior(1.; log_sd=0., source="synthetic", location="test",
        width_source="synthetic", width_location="test")

    # Log observations and a log-normal prior give an exact Gaussian posterior
    # for log(rate), independent of Turing's numerical implementation.
    model = FitModel((p,d)->fill(log(p.rate),length(d.values)); initial=(rate=1.,),
        lower=(rate=0.1,), upper=(rate=10.,))
    truth, sigma = 1.5, 0.2
    observed = log(truth) .+ sigma .* randn(MersenneTwister(821),32)
    data = (values=observed, sigma=sigma)
    result = fit_bayes(model, data, (rate=prior,); samples=1000, warmup=500,
        chains=2, predictive_draws=1000, rng=MersenneTwister(159))
    @test size(result.draws.rate) == (1000,2)
    @test result.prior_provenance.rate.source == "synthetic test"
    @test result.settings.warmup == 500
    @test result.diagnostics.passed
    @test result.summary.rate.lower < truth < result.summary.rate.upper
    variance = inv(inv(0.7^2) + length(observed)/sigma^2)
    center = variance * sum(observed)/sigma^2
    logs = log.(result.draws.rate)
    @test abs(mean(logs)-center) < 4sqrt(variance/result.summary.rate.ess_bulk)
    @test var(vec(logs)) ≈ variance rtol=0.2
    @test result.summary.rate.mean ≈ exp(center+variance/2) rtol=0.02
    @test result.summary.rate.lower ≈ exp(center+sqrt(variance)*quantile(Normal(),0.025)) rtol=0.02
    @test result.summary.rate.upper ≈ exp(center+sqrt(variance)*quantile(Normal(),0.975)) rtol=0.02

    heldout = log(truth) .+ sigma .* randn(MersenneTwister(202),400)
    future = (values=fill(missing,length(heldout)), sigma=sigma)
    predictive = posterior_predictive(model, result, future; draws=1500, rng=MersenneTwister(203))
    coverage = mean(predictive.observed.lower .<= heldout .<= predictive.observed.upper)
    @test 0.90 < coverage < 0.99
    @test mean(predictive.observed.upper-predictive.observed.lower) > 4mean(predictive.latent.upper-predictive.latent.lower)
    @test size(predictive.latent_samples) == (400,1500)
    @test length(unique(predictive.draw_indices)) == 1500
    @test_throws DomainError posterior_predictive(model, result, (values=[missing],sigma=0.))
    @test_throws ArgumentError fit_bayes(model, data, (other=LogNormal(),); samples=10)
    @test_throws ArgumentError fit_bayes(model, data, (rate=Normal(),); samples=10)
    @test_throws ArgumentError fit_bayes(model, data, (rate=prior,); samples=1)
end

@testset "NUTS through a fermentation ODE" begin
    spec = ModelSpec(; organisms=(:microbe,), substrates=(:sugar,), products=(:acid,),
                     stoichiometry=reshape([1.,-2.,1.],3,1))
    function prediction(p, data)
        problem = fermentation_problem(spec, [0.2,4.,0.] .* one(p.mu), (0.,5.),
            (mu_max=[p.mu], Ks=reshape([0.7],1,1)))
        trajectory = solve(problem, Tsit5(); abstol=1e-10, reltol=1e-9)
        successful_retcode(trajectory) || return nothing
        return [trajectory(t)[1] for t in data.times]
    end
    truth = (mu=0.8,)
    times = collect(0.25:0.25:5.)
    sigma = 0.03
    data = (; times, values=prediction(truth,(; times)) + sigma .* randn(MersenneTwister(122),length(times)), sigma)
    model = FitModel(prediction; initial=(mu=0.5,), lower=(mu=0.1,), upper=(mu=2.,))
    result = fit_bayes(model, data, (mu=LogNormal(log(0.7),0.3),);
        samples=500, warmup=400, chains=2, predictive_draws=100, rng=MersenneTwister(123))
    @test result.summary.mu.lower < truth.mu < result.summary.mu.upper
    @test result.summary.mu.mean ≈ truth.mu rtol=0.05
    @test result.summary.mu.rhat < 1.05
    @test result.diagnostics.divergences == 0
    @test ismissing(result.prior_provenance.mu)
    @test size(result.predictive.observation_samples) == (length(times),100)
end
