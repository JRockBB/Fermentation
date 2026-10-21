"""
    literature_prior(median; log_sd, source, location, width_source, width_location)

Construct a log-normal prior with explicit provenance for both its median and
log-scale standard deviation. No numerical width is inferred from a confidence
label. Pass the returned object as a named entry in `fit_bayes`'s `priors`.
"""
function literature_prior(median::Real; log_sd::Real, source::AbstractString,
                          location::AbstractString, width_source::AbstractString,
                          width_location::AbstractString)
    _finite_real(median, log_sd)
    median > 0 && log_sd > 0 || throw(DomainError((median,log_sd), "prior median and log_sd must be positive"))
    all(x -> !isempty(strip(x)), (source,location,width_source,width_location)) ||
        throw(ArgumentError("both the prior center and width require source and location"))
    provenance = (; median, log_sd, source, location, width_source, width_location)
    return (; distribution=Distributions.LogNormal(log(median), log_sd), provenance)
end

function _bayes_priors(model, priors)
    Set(keys(priors)) == Set(keys(model.initial)) || throw(ArgumentError("one prior per fitted parameter is required"))
    names = keys(model.initial)
    distributions = map(names) do name
        entry = priors[name]
        distribution = hasproperty(entry, :distribution) ? entry.distribution : entry
        distribution isa Distributions.ContinuousUnivariateDistribution ||
            throw(ArgumentError("$name requires a continuous univariate prior"))
        minimum(distribution) >= 0 || throw(ArgumentError("$name prior must have nonnegative support"))
        return distribution
    end
    provenance = NamedTuple{names}(Tuple(hasproperty(priors[name], :provenance) ?
        priors[name].provenance : missing for name in names))
    return Distributions.product_distribution(collect(distributions)), provenance
end

Turing.@model function _gaussian_posterior(model, data, indices, y, observations, joint_prior, names)
    theta ~ joint_prior
    parameters = NamedTuple{names}(Tuple(theta))
    prediction, sigma = try
        (_fit_prediction(model, parameters, data, observations.indices), _noise_sigma(observations, parameters))
    catch err
        err isa _PredictionFailure || rethrow()
        Turing.@addlogprob! -Inf
        return nothing
    end
    for k in eachindex(indices)
        y[k] ~ Distributions.Normal(prediction[indices[k]], sigma[k])
    end
    return nothing
end

function _posterior_summary(draws)
    return NamedTuple{keys(draws)}(map(values(draws)) do samples
        interval = Statistics.quantile(vec(samples), [0.025,0.5,0.975])
        return (; mean=Statistics.mean(samples), median=interval[2], lower=interval[1], upper=interval[3],
            rhat=only(MCMCDiagnosticTools.rhat(samples; kind=:rank)),
            ess_bulk=only(MCMCDiagnosticTools.ess(samples; kind=:bulk)),
            ess_tail=only(MCMCDiagnosticTools.ess(samples; kind=:tail)))
    end)
end

function _predictive_bands(matrix)
    q = [Statistics.quantile(row, [0.025,0.5,0.975]) for row in eachrow(matrix)]
    return (lower=first.(q), median=getindex.(q,2), upper=last.(q))
end

function _predictive_sigma(data, parameters)
    hasproperty(data, :values) && hasproperty(data, :sigma) || throw(ArgumentError("prediction data requires values and sigma"))
    data.values isa AbstractVector && !isempty(data.values) || throw(ArgumentError("prediction values must be a nonempty vector"))
    noise = data.sigma isa Symbol ? (hasproperty(parameters, data.sigma) ? parameters[data.sigma] :
        throw(ArgumentError("unknown predictive noise parameter $(data.sigma)"))) : data.sigma
    noise isa Real || noise isa AbstractVector || throw(ArgumentError("prediction sigma must be numeric or a fitted parameter name"))
    sigma = noise isa Real ? fill(noise,length(data.values)) : collect(noise)
    length(sigma) == length(data.values) || throw(DimensionMismatch("prediction sigma must match values"))
    _finite_real(sigma...)
    all(>(0), sigma) || throw(DomainError(sigma, "prediction standard deviations must be positive"))
    return sigma
end

"""
    posterior_predictive(model, posterior, data; draws=500, rng=Random.default_rng())

Predict at new data rows using retained posterior samples, without refitting.
`data.values` determines row count (values may all be missing); `data.sigma` is
known future observation uncertainty or a Symbol naming the sampled shared SD.
In the latter case each selected joint draw supplies its own noise scale.
Return pointwise 95% bands and samples
for both latent model means and noisy observations. Selected posterior draws are
sampled without replacement. A failed prediction raises an error rather than
silently deleting posterior samples.
"""
function posterior_predictive(model::FitModel, posterior, data; draws=500, rng=Random.default_rng())
    draws isa Integer && draws > 0 || throw(ArgumentError("predictive draws must be positive"))
    _predictive_sigma(data, model.initial)
    n = length(data.values)
    names = keys(model.initial)
    keys(posterior.draws) == names || throw(ArgumentError("posterior parameter names differ from model"))
    total = length(first(values(posterior.draws)))
    all(length(x) == total for x in values(posterior.draws)) || throw(DimensionMismatch("posterior draw lengths differ"))
    total > 0 || throw(ArgumentError("posterior contains no draws"))
    indices = Random.randperm(rng, total)[1:min(draws,total)]
    columns = map(indices) do index
        parameters = NamedTuple{names}(Tuple(posterior.draws[name][index] for name in names))
        return (; mean=_fit_prediction(model, parameters, data, 1:n), sigma=_predictive_sigma(data, parameters))
    end
    latent = hcat(getproperty.(columns, :mean)...)
    sigma = hcat(getproperty.(columns, :sigma)...)
    observed = latent + sigma .* Random.randn(rng, n, length(indices))
    return (; latent=_predictive_bands(latent), observed=_predictive_bands(observed),
             latent_samples=latent, observation_samples=observed,
             observation_sigma_samples=sigma, draw_indices=indices)
end

# Keep transition-local values in a separate scope from fit_bayes's final
# summary. Assigning a closure variable also named `diagnostics` in fit_bayes
# would box/share it between chains, even while a callback holds the I/O lock.
function _bayes_callback(names,on_draw,on_iteration,warmup)
    isnothing(on_draw) && isnothing(on_iteration) && return nothing
    monitor = !isnothing(on_iteration)
    callback_lock = ReentrantLock()
    return function (rng,model,sampler,transition,state,iteration;chain_number,kwargs...)
        draw = monitor ? iteration-warmup : iteration
        iteration_diagnostics = (divergent=get(transition.stats,:numerical_error,missing),
            tree_depth=get(transition.stats,:tree_depth,missing),
            n_steps=get(transition.stats,:n_steps,missing),
            step_size=get(transition.stats,:step_size,missing),
            acceptance_rate=get(transition.stats,:acceptance_rate,missing),
            hamiltonian_energy=get(transition.stats,:hamiltonian_energy,missing),
            hamiltonian_energy_error=get(transition.stats,:hamiltonian_energy_error,missing),
            max_hamiltonian_energy_error=get(transition.stats,:max_hamiltonian_energy_error,missing))
        lock(callback_lock) do
            if draw > 0 && !isnothing(on_draw)
                parameters = NamedTuple{names}(Tuple(transition.params[Turing.@varname(theta)]))
                on_draw(chain_number,draw,parameters,
                    (divergent=coalesce(iteration_diagnostics.divergent,false),
                     tree_depth=coalesce(iteration_diagnostics.tree_depth,0)))
            end
            if monitor
                phase = iteration == 1 ? :initialization : iteration <= warmup ? :warmup : :sampling
                on_iteration(chain_number,iteration,phase,iteration_diagnostics)
            end
        end
        return nothing
    end
end

"""
    fit_bayes(model, data, priors; samples=1000, warmup=1000, chains=4,
              target_accept=0.8, max_depth=10, predictive_draws=500,
              predictive_data=data, rng=Random.default_rng(), progress=false,
              chain_execution=:serial, on_draw=nothing, on_iteration=nothing)

Sample positive parameters with Turing NUTS and ForwardDiff using the same
Gaussian observation model as `fit_mle`, including a shared estimated SD when
`data.sigma` names its parameter. Priors are a NamedTuple of
positive-support distributions or `literature_prior` results. FitModel's MLE
bounds are not imposed on priors: use explicitly truncated priors if needed.

Chains start independently from their priors and discard all warmup. They run
serially by default; `chain_execution=:threads` uses Turing's threaded ensemble
(start Julia with multiple threads). An optional
`on_draw(chain_index, draw_index, parameters, diagnostics)` receives every retained
draw, with `diagnostics.divergent` and `diagnostics.tree_depth`. Callback calls are
serialized across threads, exclude warmup, and propagate errors. A checkpoint of
draws alone does not contain sampler adaptation state and cannot resume sampling.
`on_iteration(chain_index, iteration, phase, diagnostics)` optionally observes
every completed sampler call, including initialization and discarded warmup.
Phases are `:initialization`, `:warmup`, and `:sampling`; initialization has missing
transition diagnostics. Turing counts initialization in its discard count, so
positive `warmup` discards that call and the next `warmup-1` transitions. Monitoring
temporarily stores these calls, then removes them before summaries/prediction.
It does not change adaptation or random streams. Both callbacks are serialized
and propagate errors. A stalled call produces no completion event.
Iteration diagnostics also include the sampler's `n_steps` (leapfrog count),
`step_size`, `acceptance_rate`, `hamiltonian_energy`, `hamiltonian_energy_error`,
and `max_hamiltonian_energy_error`. These are copied from the completed transition,
without additional density evaluations. They are missing at initialization;
nonfinite energy errors on divergent transitions are preserved. The retained-draw
callback keeps its two-field divergence/tree-depth interface.
Return the raw chain, named draws (samples × chains), marginal 95%
intervals, rank R-hat, bulk/tail ESS, divergences, maximum-depth hits, provenance,
and posterior predictive bands. `diagnostics.passed` requires multiple chains,
R-hat < 1.01, both ESS >= 100*chains, and no divergence or depth-limit hits. These
checks are necessary diagnostics, not proof of convergence or model validity.
"""
function fit_bayes(model::FitModel, data, priors::NamedTuple; samples=1000, warmup=1000,
                   chains=4, target_accept=0.8, max_depth=10, predictive_draws=500,
                   predictive_data=data, rng=Random.default_rng(), progress=false,
                   chain_execution=:serial, on_draw=nothing, on_iteration=nothing)
    samples isa Integer && samples >= 10 || throw(ArgumentError("samples must be at least 10 per chain"))
    warmup isa Integer && warmup >= 0 || throw(ArgumentError("warmup must be nonnegative"))
    chains isa Integer && chains > 0 || throw(ArgumentError("chains must be positive"))
    max_depth isa Integer && max_depth > 0 || throw(ArgumentError("max_depth must be positive"))
    _finite_real(target_accept)
    0 < target_accept < 1 || throw(ArgumentError("target_accept must lie strictly between 0 and 1"))
    predictive_draws isa Integer && predictive_draws > 0 || throw(ArgumentError("predictive_draws must be positive"))
    chain_execution in (:serial,:threads) || throw(ArgumentError("chain_execution must be :serial or :threads"))
    observations = _observations(data)
    _noise_sigma(observations, model.initial)
    _predictive_sigma(predictive_data, model.initial)
    joint_prior, prior_provenance = _bayes_priors(model, priors)
    names = keys(model.initial)
    posterior_model = _gaussian_posterior(model, data, observations.indices,
        observations.y, observations, joint_prior, names)
    sampler = Turing.NUTS(warmup, Float64(target_accept); max_depth, adtype=Optimization.AutoForwardDiff())
    monitor = !isnothing(on_iteration)
    callback = _bayes_callback(names,on_draw,on_iteration,warmup)
    ensemble = chain_execution == :serial ? Turing.MCMCSerial() : Turing.MCMCThreads()
    chain = Turing.sample(rng, posterior_model, sampler, ensemble,
        monitor ? samples+warmup : samples, chains;
        num_warmup=warmup, discard_initial=monitor ? 0 : warmup, progress, verbose=false,
        initial_params=[Turing.InitFromPrior() for _ in 1:chains], callback)
    monitor && (chain = getindex(chain;iter=(warmup+1):(warmup+samples)))
    draws = NamedTuple{names}(Tuple(Array(chain[Turing.@varname(theta[i])]) for i in eachindex(names)))
    all(size(draw) == (samples,chains) for draw in values(draws)) || error("Unexpected posterior draw dimensions")
    summary = _posterior_summary(draws)
    divergences = count(isequal(true), Array(chain[FlexiChains.Extra(:numerical_error)]))
    depth_hits = count(d->!ismissing(d) && d>=max_depth, Array(chain[FlexiChains.Extra(:tree_depth)]))
    minimum_ess = 100chains
    passed = chains >= 2 && divergences == 0 && depth_hits == 0 && all(values(summary)) do row
        isfinite(row.rhat) && row.rhat < 1.01 &&
            isfinite(row.ess_bulk) && row.ess_bulk >= minimum_ess &&
            isfinite(row.ess_tail) && row.ess_tail >= minimum_ess
    end
    diagnostics = (; passed, divergences, depth_hits, rhat_threshold=1.01, minimum_ess)
    settings = (; samples, warmup, chains, target_accept, max_depth, chain_execution,
                 julia_threads=Threads.nthreads())
    result = (; chain, draws, summary, diagnostics, prior_provenance, settings)
    predictive = posterior_predictive(model, result, predictive_data; draws=predictive_draws, rng)
    return (; result..., predictive)
end
