# Runnable examples

From the repository root:

```sh
julia --project examples/core_demo.jl
```

The demo combines Monod growth, temperature and pH response functions, acid
inhibition, buffer chemistry, logarithmic populations, and two dilution feeds.
It writes `examples/output/core_demo.csv`; optionally supply a different CSV path
as the script's sole argument. Re-running replaces that generated output. Duplicate
times preserve the states immediately before and after feeding. At a temperature
change the displayed temperature is the right-continuous scheduled value.

**All constants are synthetic.** Biomass uses synthetic mmol/L equivalents so
the example's conservation weights are explicit. It does not estimate CFU,
starter peak time, or food safety. A food model needs sourced kinetics, population
conversion, yields, substrate composition, and calibrated buffer chemistry.

The demo holds its effective buffer constant across feeds, corresponding to a
fresh mixture with the same buffer composition. Its fixed pKas and pKw are
synthetic simplifying assumptions, not a temperature-dependent chemistry model.

## Fit synthetic observations

```sh
julia --project examples/fit_monod.jl
```

This example adds reproducible Gaussian noise to a Monod trajectory and fits
growth rate and half-saturation concentration from ten starting points. Outputs
in `examples/output/monod_fit/` include observations and residuals, all starting
points and optimizer results, and a TOML report with local 95% intervals and the
known synthetic truth. A directory can be supplied as the sole command-line
argument; re-running replaces generated outputs. No fitted food parameters are
written to `params/`. See [the likelihood and uncertainty assumptions](../docs/estimation.md).

## Bayesian recovery and predictive intervals

```sh
julia --project examples/bayes_monod.jl
```

This synthetic example samples the growth-rate posterior with NUTS, holding the
half-saturation concentration and yields known. It saves retained draws, sampler
diagnostics, and separate latent and noisy-observation bands at held-out times in
`examples/output/monod_bayes/`. A directory can be supplied as its sole argument.
Both the observations and prior values are synthetic. Coverage of only 20 held-out
points is descriptive; the automated test uses a larger synthetic held-out sample.
See [Bayesian assumptions and diagnostics](../docs/bayesian.md).

## Detect confounded parameters

```sh
julia --project examples/identifiability.jl
```

For the synthetic observation `a*b`, the data constrain the product but cannot
separate its factors. This example saves the sensitivity matrix, likelihood
profile, and diagnosis to `examples/output/identifiability/`. It reports rank one
for two parameters and a flat profile. An optional sole argument changes the
output directory. See [interpretation and limits](../docs/identifiability.md).

## Published cocoa baseline

```sh
julia --project scripts/prepare_cocoa.jl
julia --project examples/cocoa_baseline.jl
julia --project scripts/compare_cocoa_figures.jl
julia --project=scripts scripts/cocoa_figures.jl
```

Unlike the synthetic engine examples above, this example uses the three source
datasets and printed posterior means from supplement Tables S4-S6. It writes
physical mg/g-pulp trajectories under `examples/output/cocoa_baseline/` and a
residual report. The comparison script retains the documented discrepancy
against digitized published curves. The example provides no reconstructed
posterior band and does not claim a fully reproduced fit. See
[baseline validation](../docs/cocoa_baseline.md) and [fit reproduction](../docs/cocoa_fitting.md).
