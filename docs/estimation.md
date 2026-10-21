# Point estimation (P1.5)

`fit_mle(model, data; starts=50)` fits **positive scalar parameters** under an
independent Gaussian observation model with known standard deviations or one
jointly estimated shared standard deviation. It
uses Optimization.jl with Fminbox/LBFGS and ForwardDiff, optimizing natural logs
of parameters inside explicit finite bounds. It does not infer measurement units,
biological bounds or conversions from counts to concentrations.

For expensive models, `maxtime=60.0` supplies an optimizer time budget in seconds
per start; the default is unlimited. Optim checks it between iterations, so
compilation and a single model evaluation can exceed the budget. Bound ODE work
separately with solver step limits. A stopped, non-converged optimization remains
a failed run and does not supply accepted parameters or uncertainty intervals.

## Specify the observation model

```julia
predict = (p, d) -> p.intercept .+ p.slope .* d.times
model = FitModel(predict;
    initial=(intercept=1., slope=1.),
    lower=(intercept=0.05, slope=0.05),
    upper=(intercept=20., slope=10.))
data = (times=[0.,1.,2.], values=[2.1,3.4,5.1], sigma=0.2)
result = fit_mle(model, data; starts=10, rng=MersenneTwister(42))
```

These numbers are synthetic. Import `FermentModels` and `Random` before running
the snippet. Parameter names must match between initial values and bounds;
their order may differ. Initial values must lie strictly inside the bounds.
Parameters fixed by the experiment belong in the prediction function's closure.

`predict(parameters, data)` must return a vector in the same row order and units
as `data.values`. Any extra fields in `data`, such as time, experimental condition,
or measured variable, are available to that function. Thus one model can map an
ODE solution to population counts, log counts, acid concentration, or pH. These
observation transformations must be explicit. Gaussian noise on log10 CFU is
different from Gaussian noise on CFU; `sigma` must use the chosen scale.

For ODE predictions, promote initial states with `physical_u0 .* one(p.mu)` (or a
representative fitted scalar) so ForwardDiff can propagate through the solver.
Check `successful_retcode(solution)` and return `nothing` for expected simulation
failure. A failed prediction abandons that start and is recorded; the remaining
starts continue. Incorrect shapes and programming errors propagate to the caller.

## Likelihood and missing data

For observed rows, with residual `r_i = y_i - prediction_i`, the negative log
likelihood is:

```text
NLL = sum(0.5*(r_i/sigma_i)^2 + log(sigma_i) + 0.5*log(2*pi))
```

`sigma` is a positive scalar, a vector matching `values`, or a Symbol naming a
positive fitted parameter (for example `sigma=:noise` with `noise` included in
the model's initial values and bounds). The last option jointly estimates one
shared observation SD. Its `n*log(sigma)` normalization is retained during
optimization; dropping it would give the wrong noise estimate. Missing observations
are omitted from the likelihood and remain `missing` in returned residuals.
Uncertainties at missing rows are ignored. NaN and infinity are rejected in
observed data. `gaussian_nll(model, data, parameters)` exposes the same normalized
likelihood without parameter bounds or a prior. Correlated
errors, censored measurements, and count likelihoods require future extensions.

## Multistart results and convergence

The first start uses `model.initial`; later starts sample uniformly within the
log bounds. Supply a seeded RNG for reproducibility, or pass an explicit vector
of parameter NamedTuples as `starts`. No global random seed is changed.

Results include `parameters`, `objective` (normalized NLL), `predictions`,
`residuals`, `standardized_residuals`, `observed_indices`, `best_start`, and `runs`.
Each run records its initial values, final values, objective, convergence flag,
solver return code, iteration count, and any expected prediction failure.
It also records elapsed seconds and `time_limit_reached`, taken from the
optimizer's stop diagnostics. This flag is `missing` when a failed prediction
prevents the optimizer from returning. A time-limited run is not accepted even
if its generic return code suggests success.
`optimizer_converged` preserves the native success status. Acceptance additionally
requires `stationarity_passed`: an independently recomputed log-parameter gradient
must satisfy `maximum(abs(z-clamp(z-gradient,lo,hi))) <= gtol`. Here `lo` and `hi`
are log bounds. This removes outward descent at an active bound while retaining
feasible descent directions. The norm and tolerance are saved as
`projected_gradient_inf` and `gradient_tolerance`. Native success from an unchanged
step alone cannot bypass this check. If the derivative solve fails, its finite
parameter point is retained, the norm is `missing`, and stationarity fails.
Only successful, stationary, finite runs are eligible as the best fit. If none converges,
`converged=false` and estimates are `nothing`; diagnostics remain available.
Stationarity is a necessary local check, not a proof of a local or global minimum.
Multiple starts improve the search but do not certify a global minimum.

`maxiters` limits each inner LBFGS solve. `outer_iterations` limits Fminbox barrier
updates; `gtol` sets inner and outer gradient tolerances. Inspect convergence and
the spread of final objectives before interpreting a fit.

Pass `on_start_complete=(index, run) -> ...` to save each completed start before
the next one begins. It receives both successful starts and expected prediction
failures. Callback errors propagate so a failed checkpoint does not silently
lose diagnostics while the search continues.

## Local uncertainty

At the best fit, ForwardDiff computes the mean-prediction Jacobian divided by
the observation SDs with respect to log parameters, freezing the SD at the fit
for that derivative. With an estimated shared scale, expected Gaussian
information additionally includes `2n * grad(log(sigma))*grad(log(sigma))'`.
An SVD gives the inverse-information covariance approximation. For each
parameter the reported interval is `exp(log(estimate) ± 1.959964*standard_error)`.
Known variances are not rescaled from residuals. An inferred sigma is its MLE
inside the likelihood; this is an asymptotic interval, not a finite-sample
Student-t or chi-square interval.

`uncertainty.status` explains whether these **approximate local 95% intervals**
are available. Intervals are withheld for numerical rank deficiency, active
bounds, nonfinite derivatives, or intervals crossing the fit bounds. Rank uses
a relative singular-value threshold of `sqrt(eps)`; a bound is treated as active
within `1e-5 * max(1, log-bound width)`. These are numerical diagnostic thresholds,
not biological parameters. A rank failure reports missing uncertainty rather than
silently using a pseudoinverse that could suggest false precision.

These intervals rely on local linearity and an interior optimum. They do not
replace profile likelihoods (P1.7), Bayesian inference (P1.6), or checks on the
observation model. Passing a synthetic recovery test is not evidence of calibrated
food predictions or guaranteed repeated-sample coverage.

Validation includes agreement with an independently computed weighted linear
least-squares solution and covariance, noisy two-parameter Monod recovery, and
explicit boundary, rank-deficiency, failure, and missing-data checks. Estimated
noise is checked against the exact Gaussian mean/SD MLE, expected information,
and the analytic profile likelihood for sigma.

Library references: [OptimizationOptimJL](https://docs.sciml.ai/Optimization/stable/optimization_packages/optim/)
and [Optimization automatic differentiation](https://docs.sciml.ai/Optimization/stable/API/ad/).
