# Sensitivity and identifiability (P1.7)

For an estimated shared noise scale, standardized mean sensitivities freeze
sigma at the supplied parameter values. A noise-only parameter has a zero mean
sensitivity column, although its variance contribution can be identified by the
likelihood. The MLE information matrix includes that contribution; the mean
Jacobian alone does not. Profile likelihood supports both the noise parameter
and kinetic parameters with the noise scale refitted as a nuisance parameter.

The shared inference tools can diagnose parameters that the observations do not
constrain. They use the same `FitModel`, known or estimated Gaussian uncertainties,
and missing-observation handling as [point estimation](estimation.md).

## Local sensitivities

```julia
fit = fit_mle(model, data; starts=10)
sens = local_sensitivities(model, data, fit.parameters)
sens.jacobian
sens.parameter_names
sens.rank
sens.unresolved_directions
```

The default matrix has entries
`J[i,j] = (1/sigma[i]) * d prediction[i] / d log(parameter[j])`.
Rows correspond to `observed_indices`, columns to `parameter_names`. Relative
parameter changes and measurement noise are therefore on comparable scales.
Use `scale=:parameter` for derivatives in physical parameter units and
`standardized=false` to omit division by sigma. Column norms measure local
influence; influence alone does not establish that a parameter can be estimated
independently of other parameters.

The SVD reports singular values, numerical rank, and condition number. Singular
values below `rtol * maximum(singular_values)` are treated as unresolved, with
default `rtol=sqrt(eps(Float64))`. Rank deficiency produces an infinite condition
number and `status=:rank_deficient`. Columns of `unresolved_directions` are the
corresponding right singular vectors, in the chosen parameter coordinates. The
null space is retained even with fewer observations than parameters.

This uses [ForwardDiff's Jacobian API](https://juliadiff.org/ForwardDiff.jl/stable/user/api/),
including differentiation through the ODE solver. Check solver tolerances and
rank-tolerance dependence before interpreting a nearly singular matrix. A local
rank check is not a proof of global or structural identifiability.

## Profile likelihood

```julia
using Random
profile = profile_likelihood(model, data, :mu, range(0.4, 1.2; length=41);
                             starts=10, rng=MersenneTwister(42))
profile.status
profile.points
profile.regions
```

The example grid is synthetic. Choose bounds and resolution appropriate to the
actual model. The helper first fits all parameters, then fixes the selected
parameter at each grid value and optimizes the others. Each reduced fit starts
from the reference nuisance estimates when interior, or the model's initial
values at a bound, plus log-uniform random starts. Grid values must increase
strictly and remain within the original bounds. All starts and failures remain
available in `reference.runs` and each point's `runs`.

Each successful point reports `deviance = 2 * (profile_nll - reference_nll)`.
The default nominal 95% cutoff is the 0.95 quantile of chi-square with one degree
of freedom. `supported_values` lists sampled values below that cutoff. `regions`
keeps disconnected supported stretches separate and reports adjacent brackets
for crossings. It does not interpolate or claim that the grid supplies exact
interval endpoints. This approach follows the profile-likelihood method of
[Raue et al. (2009), doi:10.1093/bioinformatics/btp358](https://pubmed.ncbi.nlm.nih.gov/19505944/).

Statuses distinguish the numerical outcomes:

| Status | Meaning |
|---|---|
| `flat_on_grid` | Every point converged with absolute deviance at most `flat_tol` (default `1e-4`). |
| `crossings_bracketed` | Every supported stretch has successful excluded neighbors on both sides. |
| `grid_limited` | At least one supported stretch reaches an edge of the chosen grid. |
| `no_supported_grid_points` | The grid missed the supported region; refine or extend it. |
| `incomplete_grid` | At least one point failed. Failed points create gaps, never crossing brackets. |
| `better_profile_optimum` | A grid fit improved on the reference beyond `objective_tol * max(1,abs(reference_nll))`; support reporting is withheld. Increase starts or revisit the optimization. |
| `no_reference_fit` | No full-model start converged; its diagnostics are returned without a profile. |

The cutoff is an asymptotic, regular-model approximation. It need not provide
95% coverage with active bounds, rank deficiency, small samples, or multiple
modes. Inspect `reference.uncertainty`, all starts, and the profile shape.
Nuisance-parameter bounds can artificially close a profile. A flat finite grid
signals an unresolved direction over that grid, not an unlimited confidence
interval or a mathematical proof of structural non-identifiability.

## Known non-identifiable example

For `prediction = a*b`, data determine the product, while increasing `a` and
decreasing `b` by the same relative amount leaves predictions unchanged. The
synthetic test and `examples/identifiability.jl` verify both a rank-one Jacobian
for two parameters and a flat profile with `b = fitted_product/a`. Independent
checks compare an identifiable linear-model profile to its analytical Gaussian
solution and ODE sensitivities to central finite differences.
