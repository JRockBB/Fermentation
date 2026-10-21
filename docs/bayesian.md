# Bayesian estimation (P1.6)

`fit_bayes(model, data, priors)` uses Turing NUTS with ForwardDiff, sharing the
`FitModel` and independent Gaussian observation likelihood from
[point estimation](estimation.md). Measurement SDs can be known inputs, or
`data.sigma=:noise` can name a shared positive parameter to estimate jointly.
That parameter requires its own prior. Missing observations contribute no likelihood.

## Explicit priors

Supply one continuous positive-support prior per parameter, keyed by name:

```julia
using FermentModels, Random
using Distributions: LogNormal
priors = (mu=LogNormal(log(0.8),0.3), Ks=LogNormal(log(0.7),0.5))
posterior = fit_bayes(model, data, priors; rng=MersenneTwister(42))
```

These values are synthetic examples. Bare distributions are permitted, but their
returned provenance is `missing`. For literature-based defaults use
`literature_prior(median; log_sd, source, location, width_source, width_location)`.
It builds a log-normal distribution whose **median** is the literature value and
records separate references for the center and log-scale width. A source giving
only a point estimate does not determine a prior width. The implementation never
maps `high`, `medium`, or `low` confidence to an invented numerical width.

`FitModel` bounds govern MLE only. Bayesian support comes from the supplied priors;
use explicitly truncated distributions if there are hard parameter limits. Invalid
simulations (`predict` returns `nothing` or nonfinite observed predictions) have
zero posterior density. Programming errors still propagate. A prior dominated by
invalid simulations needs revision based on the model and evidence, not silent
clipping of sampled values.

## Sampling and diagnostics

Defaults are 1,000 retained draws per chain, 1,000 warmup steps, four serial chains,
and target acceptance 0.8. Chains initialize independently from their priors.
All warmup draws are discarded. An explicit RNG controls sampling and subsequent
predictive simulation without changing global random state; exact sequences can
differ between Julia/package versions.

Set `chain_execution=:threads` and start Julia with `--threads=4` to run chains
concurrently. Settings record the execution mode and Julia thread count. Each
chain uses an independent seeded stream; serial and threaded modes reproduce
the same parameter draws for a fixed seed in the tested environment. Subsequent
predictive noise draws can differ between modes because their parent RNG states
differ after sampling.

For durable progress, pass
`on_draw=(chain_index, draw_index, parameters, diagnostics) -> ...`.
The callback receives each retained draw, excluding warmup, with divergence and
tree-depth diagnostics. Calls are serialized across threads, so callbacks can
write files or collect results without a separate lock. Errors propagate. Draw
files preserve diagnostic output but lack adaptation state needed to resume
NUTS; incomplete runs still need full convergence assessment.

To observe initialization and warmup as well, pass
`on_iteration=(chain_index, iteration, phase, diagnostics) -> ...`.
Phases are `:initialization`, `:warmup`, and `:sampling`. Initialization has
missing divergence and tree-depth fields. Turing includes initialization in
the discard count: a positive `warmup` discards that call and the next
`warmup-1` transitions. Monitoring preserves adaptation and seeded sampling;
it temporarily stores discarded calls, then removes them before constructing
the returned posterior and predictions. Both callbacks are serialized together,
with `on_draw` first when a draw is retained. Each event reports a completed
call, so a stalled call produces no new event. The
[benchmark runner](cocoa_benchmark.md) uses these events for durable runtime
records and optional stops that preserve partial checkpoints.

Iteration diagnostics also carry `n_steps`, `step_size`, `acceptance_rate`,
`hamiltonian_energy`, `hamiltonian_energy_error`, and
`max_hamiltonian_energy_error`, copied directly from the completed sampler
transition without extra density evaluations. `n_steps` counts leapfrog steps;
it does not count ODE integration steps. Initialization leaves these fields
missing, and nonfinite energy errors on divergent transitions remain visible.
The `on_draw` diagnostic tuple still contains just divergence and tree depth.

The result contains the raw Turing chain, named `draws` matrices (draws × chains),
`summary`, `diagnostics`, `settings`, `prior_provenance`, and `predictive`.
Parameter summaries give posterior means, medians, marginal 95% credible intervals,
rank-normalized split R̂, and bulk/tail effective sample sizes. Diagnostics pass
only with at least two chains, all R̂ < 1.01, both ESS ≥ 100 × chain count, no
divergences, and no maximum-depth hits. Failed diagnostics remain visible in the
result; the code does not relabel a completed sampling run as converged.

Those thresholds follow the [MCMCDiagnosticTools guidance](https://turinglang.org/MCMCDiagnosticTools.jl/)
and the [rank-normalization paper](https://arxiv.org/abs/1903.08008).
They cannot prove that all posterior modes were explored or that the model is
scientifically adequate. Inspect the chains, priors, and residuals as well.

## Predictive uncertainty

The default predictive report uses up to 500 retained posterior draws without
replacement. It supplies two pointwise 95% bands:

- `latent`: variability in the model's predicted observation mean.
- `observed`: variability in a future observation, adding Gaussian measurement
  noise with the supplied `sigma` or the selected posterior draw's fitted sigma.

Both sets of simulated values are returned alongside the bands. These are
pointwise intervals, not a simultaneous band for an entire trajectory. Their
scale follows the observation function; a band for log CFU remains on that scale.

For held-out inputs, call
`posterior_predictive(model, posterior, future_data; draws=500, rng=...)` without
refitting. `future_data.values` sets the row count and may consist entirely of
`missing`; future `sigma` can be known and positive or name a sampled scale.
`observation_sigma_samples` records the SD used for every row and selected draw,
preserving its joint association with the kinetic parameters. A failed prediction
raises an error rather than removing difficult draws and narrowing the band.

Validation compares NUTS with an exact conjugate posterior in log parameters,
checks held-out synthetic observation coverage, and propagates gradients through
a fermentation ODE for parameter recovery. The inferred-scale posterior is also
checked against independent one-dimensional quadrature with a half-Cauchy prior.
The [cocoa fit workflow](cocoa_fitting.md) supplies the original food-specific
priors and reports convergence separately from completed sampling.

API references: [Turing sampling options](https://turinglang.org/docs/usage/sampling-options/)
and [Bayesian differential equations](https://turinglang.org/docs/tutorials/bayesian-differential-equations/).
