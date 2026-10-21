# Cocoa baseline (P2.3)

The eight-state M4 equations are implemented in `src/cocoa/baseline.jl` with
five additive growth routes and three product-dependent mortality terms.
The [transcription](cocoa_equations.md) maps every term and parameter to the
paper and the archived Stan model. This is a source-model reproduction, not
a general-purpose recipe or a newly calibrated cocoa process.

## Coordinates and API

`COCOA_STATES` is `(G,F,E,L,A,Xy,Xl,Xa)`. Time is hours. Physical states use
mg/g pulp, including dry microbial biomass. The original normalized coordinates
divide each state by its own trial maximum, including the fixed initial row.
No conversion from a pulp-mass basis to mmol/L is inferred.

```julia
using FermentModels, OrdinaryDiffEq
include("scripts/prepare_cocoa.jl")
data = CocoaSourceData.cocoa_trial(:camu)
p = load_params("params/cocoa_camu_scaled.toml").values
problem = cocoa_problem(data.initial ./ data.scales, (0.0,144.0), p)
solution = solve(problem, Tsit5(); abstol=1e-11, reltol=1e-10)
physical_at_24h = solution(24.0) .* data.scales
```

`cocoa_physical_parameters(p, scales)` converts the 24 kinetic parameters to
physical coordinates; `cocoa_scaled_parameters` is its inverse. Both exclude
extra fields, including `sigma`. Physical observation SDs are
`p.sigma .* scales`, not a shared SD in mg/g pulp. `yc3` uses the fructose
maximum, following supplement S3 and correcting the documented source-script
rescaling error. Rate and mortality parameters in the normalized TOMLs have
units h^-1; the other coefficients and sigma are dimensionless.

`cocoa_problem` validates dimensions and signs, promotes initial states for
ForwardDiff, and rejects negative or nonfinite accepted solver states.
`cocoa_rhs!` evaluates negative trial states at zero. The Contois flux is
continuously extended to zero at lactate=AAB=0. Derivative tests concern the
smooth interior; this boundary extension is not a claim of differentiability
at the origin. The source omits CO2 and dead biomass, so an eight-state total
mass-conservation assertion would be inappropriate.

## Data and published parameter extraction

The three parameter TOMLs contain the printed normalized posterior means from
supplement Tables S4-S6, pages 11-13, with DOI, table, units, original fit script,
and dataset provenance for every value. The [raw transcription](../data/raw/transcribed/cocoa2018/README.md)
also retains posterior SDs and marginal 2.5%, 50%, and 97.5% quantiles.
The `medium` confidence label reflects rounded, trial-specific fitted values;
it does not set a prior width. These are published estimates, not fresh fits.

Run `scripts/prepare_cocoa.jl` to verify the raw inputs and generate the parameter
TOMLs and tidy observations in `data/processed/cocoa/moreno2018/`. Initial rows
are explicitly separated from the 336 likelihood observations. Raw files are
never overwritten. `examples/cocoa_baseline.jl` writes dense physical
trajectories and normalized residual summaries under `examples/output/`.

## Validation and the remaining reproduction gap

Tests cover analytical mortality and additive-growth limits, the nonnegative
boundary vector field, normalized/physical trajectory equivalence, all three
published parameter sets, and ForwardDiff derivatives through the solver.

For an external check, 168 central-curve points were digitized and visually
reviewed from article Figures 3-5 as reproduced in the author's thesis Figures
2.3-2.5. The [digitization record](../data/raw/digitized/moreno2018-figures/README.md)
includes source images, calibration, attribution, and a reproducible extraction
script. `scripts/compare_cocoa_figures.jl` preserves a preselected tolerance of
2% of each panel's calibrated vertical range.

The rounded-mean parameter simulations **do not pass that figure-reproduction
criterion**. Maximum discrepancies are about 6.0% (Camu), 7.1% (box 1), and
4.9% (box 2). The test suite retains this as an explicit expected-failure check.
The [comparison table](../data/processed/cocoa/moreno2018/figure_comparison.csv)
and TOML summary report individual points and panels; tolerances are not tuned.

The archive lacks posterior draws and the plotting script. Marginal parameter
means do not determine a nonlinear model's posterior mean trajectory, and the
printed coefficients are rounded. Those are possible contributors, not an
established explanation of every difference. P2.4 must reproduce the joint
posterior and identify the plotted summary before this acceptance gap can close.
Do not construct bands by independently sampling the reported marginal quantiles.

The [paper's figure captions](https://upload.wikimedia.org/wikipedia/commons/5/50/A_mathematical_model_of_cocoa_bean_fermentation.pdf)
(Figures 3-5, printed pages 9-10) identify simulations and 95% prediction
intervals, but do not specify whether the central curve is a posterior mean,
median, or trajectory at a point estimate. The fit comparison therefore reports
mean and median posterior trajectories separately without assigning either
interpretation to the published curve.

The latest [P2.4 fits](cocoa_fitting.md#recorded-results-and-remaining-work)
pass convergence and parameter-interval checks for all three datasets. Box 1's
joint-posterior mean and median curves both pass the 2% figure criterion. Camu's
median passes, while its mean misses one point (2.04%); box 2's mean passes,
while its median misses one point (2.07%). Neither central summary passes across
all three trials. These results do not identify the source's plotted summary
or change the failed regression at printed parameter means.

Evaluating the printed marginal medians also fails the fixed figure criterion:
maximum errors are 6.80% for Camu, 7.92% for box 1, and 6.32% for box 2.
`scripts/compare_cocoa_summaries.jl` reproduces this check alongside the printed
means, preserving input hashes, units, and all 336 comparisons in the
[point-summary report](../data/processed/cocoa/moreno2018/point_summary_comparison.toml)
and [point table](../data/processed/cocoa/moreno2018/point_summary_comparison.csv).
These marginal-parameter trajectories do not resolve the unknown source
plotting statistic.

The normalized observation RMSEs at the printed means are 0.1428, 0.1576, and
0.1600, respectively. They are descriptive fit errors, not a separate validation
dataset. The source's shared unknown sigma and exact priors are required for
fit reproduction.

## Figures

After running the baseline example, generate the figures with
`julia --project=scripts scripts/cocoa_figures.jl`. Red lines are simulations at
the printed parameter means, black dots are source observations, and blue
crosses are digitized published central curves. No new uncertainty band is
claimed.

- [Camu](figures/cocoa_camu.png)
- [Box 1](figures/cocoa_box1.png)
- [Box 2](figures/cocoa_box2.png)
