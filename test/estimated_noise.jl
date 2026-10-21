using Distributions: Cauchy, truncated, logpdf
using Statistics
using Random

@testset "Estimated Gaussian scale: analytic MLE and profile" begin
    y = 2 .+ 0.3randn(MersenneTwister(429), 40)
    data = (values=y, sigma=:noise)
    model = FitModel((p,d)->fill(p.level, length(d.values)); initial=(level=1.5, noise=0.5),
        lower=(level=0.1, noise=0.01), upper=(level=5., noise=2.))
    result = fit_mle(model, data; starts=3, rng=MersenneTwister(18))
    level = mean(y)
    sigma = sqrt(mean(abs2, y .- level))
    @test result.converged
    @test result.parameters.level ≈ level rtol=1e-6
    @test result.parameters.noise ≈ sigma rtol=1e-6
    @test result.objective ≈ length(y)*(log(sigma)+log(2pi)/2+0.5) rtol=1e-8
    @test result.standardized_residuals ≈ (y .- level) ./ sigma rtol=1e-5
    @test result.uncertainty.status == :ok
    @test result.uncertainty.covariance_log ≈ Diagonal([sigma^2/(length(y)*level^2), 1/(2length(y))]) rtol=1e-6
    @test gaussian_nll(model, data, (level=level,noise=2sigma)) ≈
        length(y)*(log(2sigma)+log(2pi)/2+0.125)
    # A noise-only parameter has no mean sensitivity, despite its nonzero
    # variance information in the likelihood.
    sensitivity = local_sensitivities(model, data, result.parameters)
    @test sensitivity.column_norms.noise == 0
    @test sensitivity.rank == 1
    grid = sigma .* [0.5,0.75,1.,1.5,2.]
    profile = profile_likelihood(model, data, :noise, grid; starts=1)
    @test profile.reference_valid
    @test all(p -> p.converged, profile.points)
    expected = length(y) .* (2log.(grid ./ sigma) + sigma^2 ./ grid.^2 .- 1)
    @test getproperty.(profile.points, :deviance) ≈ expected atol=1e-7
    # Profiling the mean must refit sigma, not freeze it at the reference MLE.
    level_grid = level .+ sigma .* [-1.,0.,1.]
    mean_profile = profile_likelihood(model, data, :level, level_grid; starts=1)
    @test mean_profile.reference_valid && all(p -> p.converged, mean_profile.points)
    @test getproperty.(mean_profile.points, :deviance) ≈
        length(y) .* log.(1 .+ ((level_grid .- level) ./ sigma).^2) atol=1e-7
    @test [p.parameters.noise for p in mean_profile.points] ≈
        sqrt.(sigma^2 .+ (level_grid .- level).^2) rtol=1e-5
    missing_data = (values=[missing;y],sigma=:noise)
    masked = fit_mle(model, missing_data; starts=1)
    @test masked.parameters.noise ≈ sigma rtol=1e-6
    @test ismissing(first(masked.standardized_residuals))
    @test_throws ArgumentError fit_mle(model, (values=y,sigma=:absent); starts=1)
end

@testset "Estimated scale posterior against independent quadrature" begin
    y = 2 .+ 0.3randn(MersenneTwister(631), 24)
    prior = truncated(Cauchy(0.,1.); lower=0.)
    model = FitModel((p,d)->fill(2., length(d.values)); initial=(noise=0.5,),
        lower=(noise=0.01,), upper=(noise=3.,))
    data = (values=y,sigma=:noise)
    result = fit_bayes(model, data, (noise=prior,); samples=1000, warmup=500, chains=2,
        predictive_draws=1000, rng=MersenneTwister(632))
    # Numerical quadrature in log(sigma) of the explicit scalar posterior,
    # including its change-of-variable Jacobian; independent of the fit code.
    z = collect(range(-7.,4.; length=20001))
    squared_error = sum(abs2, y .- 2.)
    logweights = [logpdf(prior,exp(x)) - length(y)*x - squared_error/(2exp(2x)) + x for x in z]
    weights = exp.(logweights .- maximum(logweights))
    weights ./= sum(weights)
    reference_mean = sum(exp.(z) .* weights)
    cumulative = cumsum(weights)
    reference_bounds = exp.(z[searchsortedfirst.(Ref(cumulative),[0.025,0.975])])
    @test result.diagnostics.passed
    @test result.summary.noise.mean ≈ reference_mean rtol=0.03
    @test result.summary.noise.lower ≈ reference_bounds[1] rtol=0.04
    @test result.summary.noise.upper ≈ reference_bounds[2] rtol=0.04
    selected_sigma = vec(result.draws.noise)[result.predictive.draw_indices]
    @test result.predictive.observation_sigma_samples == repeat(selected_sigma', length(y), 1)
    heldout = 2 .+ 0.3randn(MersenneTwister(633), 400)
    predictive = posterior_predictive(model, result, (values=fill(missing,400),sigma=:noise);
        draws=1500, rng=MersenneTwister(634))
    @test 0.90 < mean(predictive.observed.lower .<= heldout .<= predictive.observed.upper) < 0.99
    @test predictive.latent.lower == predictive.latent.upper == fill(2.,400)
end
