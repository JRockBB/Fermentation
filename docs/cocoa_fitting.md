# Cocoa fit reproduction (P2.4)

`scripts/fit_cocoa.jl` fits each archived dataset with the eight source equations,
24 normalized kinetic parameters, and one jointly estimated observation SD.
It preserves settings, input hashes, estimates, and convergence diagnostics in
a new output directory. Existing fit directories are never overwritten.
The [baseline comparison](cocoa_baseline.md) remains a separate unmet acceptance
check; implementing a sampler does not establish reproduction of the paper.

## Observation model and priors

`cocoa_observations(initial, times, observations)` takes physical states in
mg/g pulp and returns normalized, state-major likelihood data. Each state is
divided by its maximum across the initial row and subsequent observations.
The initial row fixes the ODE initial condition and is excluded from the
likelihood. The three datasets contribute 128, 104, and 104 observations.

`data.sigma=:sigma` asks the shared inference engine to estimate a positive
common SD. The likelihood includes both the residual sum of squares and
`n*log(sigma)`. In physical coordinates the SD for state i is
`sigma*scales[i]`. Predictive observation draws use the noise scale from the
same joint draw as the latent trajectory.

`cocoa_priors()` returns the exact source distributions, with provenance:
independent positive-truncated Normal(0.5,0.3) kinetic priors and a half-Cauchy
prior with scale 1 for sigma. Their support is unbounded above. No prior width
is inferred from the parameter files' confidence labels or rounded summaries.

`cocoa_fit_model` uses Tsit5 with absolute and relative tolerances 1e-6 and a
step limit of 100,000,000, matching the archived numerical tolerances and step
cap. The original integrator is Stan RK45. Turing chains initialize independently
from the priors; the original drivers initialize with absolute Normal(0.5,0.2)
draws. Turing uses its default heuristic search for the initial leapfrog step
size, whereas the source drivers specify 0.01. Subsequent step sizes adapt during
warmup. These algorithmic differences are retained explicitly.

`solver=:rodas5p` (CLI `--solver=rodas5p`) selects a stiff integrator for difficult
optimization proposals. Tsit5 remains the default. Both methods integrate the
same source equations with the requested tolerances and step cap; no parameters,
priors, or likelihood terms change. Solver choice is saved in `settings.toml`.

## Running fits

From the repository root, after installing the main project dependencies:

```sh
julia --startup-file=no --check-bounds=yes --project scripts/prepare_cocoa.jl
julia --startup-file=no --check-bounds=yes --project scripts/fit_cocoa.jl --trial=camu --method=mle --output=data/processed/cocoa/moreno2018/fits/camu_mle_new
julia --startup-file=no --check-bounds=yes --project scripts/fit_cocoa.jl --trial=camu --method=bayes --output=data/processed/cocoa/moreno2018/fits/camu_bayes_full
```

Choose `--trial=box1` or `--trial=box2` for the other original datasets, or
`--trial=all` to run all three serially in one Julia process, with a separate
subdirectory per trial. The
Bayesian defaults match the source's 2,000 retained samples per chain, 1,000
warmup iterations, four chains, target acceptance 0.985, maximum tree depth 12,
and seed 100785. Chains run serially by default. Use Julia's `--threads=4` and
the runner's `--chain_execution=threads` to sample four independent chains
concurrently. The full configuration is a substantial
ODE inference job; a smaller pilot can check the workflow:

```sh
julia --startup-file=no --check-bounds=yes --project scripts/fit_cocoa.jl --trial=camu --method=bayes --samples=200 --warmup=200 --chains=2 --output=data/processed/cocoa/moreno2018/fits/camu_pilot_new
julia --startup-file=no --check-bounds=yes --threads=4 --project scripts/fit_cocoa.jl --trial=camu --method=bayes --chain_execution=threads --output=data/processed/cocoa/moreno2018/fits/camu_threaded_full
```

Run one Julia process at a time on memory-constrained machines. Independent fit
processes each load and compile the inference stack; an earlier batch of fits
and tests exceeded this Windows machine's paging-file limit. Use `--trial=all`
to process datasets sequentially in one process, with two or four sampling
threads within that process as resources permit.

MLE uses ten starts, beginning at the printed posterior means, and explicit
computational bounds 1e-5 to 5 for all normalized parameters. The remaining
starts are uniform in log coordinates. These bounds are search choices, not
published biological constraints, and do not restrict the Bayesian posterior.
Active-bound or rank-deficient fits do not receive local confidence intervals.
`--starts`, `--maxiters`, and `--outer_iterations` configure the bounded optimizer.
`--solver_maxiters` separately caps each ODE solve. Its default matches the source
step limit; lowering it can prevent extreme MLE proposals from consuming a very
long solve. A proposal that hits this limit is an explicit failed start, not a
converged estimate or a change to the differential equations.
`--maxtime=60` additionally gives each MLE start a 60-second optimizer budget
(default unlimited). This is a soft limit checked between iterations, not a
deadline that interrupts Julia compilation or an individual ODE call. Use both
controls when exploring the broad MLE bounds.

Every completed MLE start is saved immediately, with its native return code,
elapsed seconds, and the optimizer's explicit time-limit flag. Every retained
Bayesian draw is appended to a separate file for its chain, together with
divergence and tree-depth diagnostics. This preserves completed work if a later
start, chain, or predictive calculation fails. Warmup is excluded. These files
are diagnostic checkpoints, not resumable sampler states or accepted partial fits.

For MLE, `--initial_file=<mle_start_NNN.toml>` starts a new optimization from a
saved start's `[parameters]`, keeping the same strict search bounds. The output
directory must still be new; the input path and SHA-256 are recorded. This
reuses a parameter point, not the old optimizer's state, and does not mark the
old run as converged. Choose the matching trial. The option is unavailable for
Bayesian fits and `--trial=all`; Bayesian initialization continues to use the
source priors.

MLE acceptance now also requires an independently recomputed projected gradient
at or below 1e-6. Reports preserve `optimizer_converged`, `stationarity_passed`,
`projected_gradient_inf`, and `gradient_tolerance` separately from the native
return code. This prevents small-step termination from being accepted solely
because the optimizer returned `Success`. Historical runs retain their original
flags; use their subsequent diagnostic reports when assessing stationarity.

To inspect a completed MLE start, including one that did not converge:

```sh
julia --startup-file=no --check-bounds=yes --project scripts/cocoa_mle_diagnostics.jl data/processed/cocoa/moreno2018/fits/box1_mle_restart data/processed/cocoa/moreno2018/fits/box1_diagnostics_new
```

The optional third argument selects the start number (default 1). The new output
directory contains `report.toml` and `parameters.csv`, linked to the checkpoint,
fit settings, raw data, and diagnostic implementation by SHA-256. Diagnostics
reevaluate the likelihood, check log-parameter derivatives against central
finite differences at two step sizes, and compare integrators and tighter
tolerances at the same saved point. Failed verification solves remain explicit.

The projected gradient is `z - clamp(z - gradient(NLL), log(lower), log(upper))`,
where `z=log(parameters)`. It retains feasible descent directions while removing
an outward gradient at an exact search bound. The report also lists parameters
within 1e-5 of either end of their log search range. This bound proximity is
descriptive and does not round a parameter onto a bound for the gradient check.

Local information uses the sigma-standardized prediction Jacobian in log
coordinates, augmented by the shared-noise contribution `sqrt(2n)` in the
log-sigma column. Omitting that row would incorrectly call sigma unidentifiable
because it does not affect the ODE means. Singular values, numerical rank,
condition number, and the weakest direction describe this point and numerical
tolerance; they are not a structural-identifiability proof or confidence
intervals. The conditional noise MLE `sqrt(RSS/n)` provides a separate analytical
check. These diagnostics do not change the original convergence flags or accept
an unfinished optimization.

For a broader check of derivative accuracy at the saved point:

```sh
julia --startup-file=no --check-bounds=yes --project scripts/cocoa_mle_numerics.jl data/processed/cocoa/moreno2018/fits/box1_mle_restart_long data/processed/cocoa/moreno2018/fits/box1_numerics_new
```

This compares Rodas5P and Tsit5 with absolute/relative tolerances 1e-9/1e-8 and
1e-11/1e-10, using central finite-difference log steps 1e-3 through 1e-6. The
original per-solve step cap remains in force and failures are retained. Bounds,
data, likelihood, and the parameter point remain fixed. Smaller differences
are not automatically more accurate because subtracting two approximate
likelihoods amplifies integration error. Solver tolerances govern local state
error; see the [SciML options documentation](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/).

## Reports and acceptance

- `settings.toml`: numerical choices, source identifiers, input and code hashes, and
  Julia version. A directory containing only settings is an incomplete run.
- `report.toml`: elapsed time, optimizer runs or sampler diagnostics, parameter
  estimates, and comparison with published marginal 95% credible intervals.
- `mle_start_NNN.toml`: each completed optimization start, written before the next.
- `chain_NNN.csv`: retained draws and sampling diagnostics, saved during sampling.
- `parameter_comparison.csv`: each estimate beside its source mean and interval;
  Bayesian runs also include rank R-hat and bulk/tail ESS.
- `draws.csv`: retained joint Bayesian draws with chain and iteration indices.
- `predictive.csv`: physical posterior mean/median, latent 95% intervals, and
  observation 95% intervals at the original observation times.

Bayesian diagnostic acceptance requires multiple chains, every rank R-hat below
1.01, bulk and tail ESS at least 100 times the chain count, and no retained
divergences or depth-limit hits. A short pilot will often fail these checks.
Its predictions and parameter summaries remain diagnostic output, not accepted
estimates. `reproduction_accepted` in the fit report additionally requires all
point estimates to fall inside the paper's reported intervals. For MLE that is
a comparison against a Bayesian reference, not a guarantee implied by the
likelihood optimum. This flag concerns the P2.4 parameter criterion only.

Compare the joint-posterior predictions with the external figure points using:

```sh
julia --project scripts/compare_cocoa_fit.jl data/processed/cocoa/moreno2018/fits/camu_pilot_new
julia --project=scripts scripts/cocoa_fit_figures.jl data/processed/cocoa/moreno2018/fits/camu_pilot_new
julia --project=scripts scripts/cocoa_trace_figures.jl data/processed/cocoa/moreno2018/fits/camu_pilot_new
```

This writes `curve_comparison.csv` and `.toml`, reports both posterior mean and
median errors using the unchanged 2% axis-range tolerance, and gates acceptance
on convergence. The archived plotting script is absent, so the paper's plotted
central summary remains unidentified. Agreement with one summary alone would
not establish the source's plotting method.

The plotting command writes physical-unit predictive PNG/SVG figures with
convergence status prominently displayed. It separates latent and observation
uncertainty, and does not clip negative Gaussian observation bounds.
The trace command writes `traces.png` and `.svg` for all 25 parameters, coloring
panel titles red when their recorded R-hat or ESS check fails. Every retained
draw is shown by chain, with warmup excluded. These plots supplement the
numerical convergence checks; they do not replace them.

To separately inspect trajectories at the printed marginal parameter means and
medians, run:

```sh
julia --startup-file=no --check-bounds=yes --project scripts/compare_cocoa_summaries.jl
```

This writes
`point_summary_comparison.csv` and `.toml` under the processed cocoa directory,
with input and implementation hashes, physical units, and the same 2% criterion.
It does not refit data or treat marginal summaries as joint posterior draws.

## Recorded results and remaining work

The October 19, 2026 runs are retained under
[`data/processed/cocoa/moreno2018/fits/`](../data/processed/cocoa/moreno2018/fits/README.md),
including unsuccessful and interrupted runs. The latest Bayesian fit for each
dataset passes convergence and places all 25 posterior means inside the paper's
reported intervals. The parameter-reproduction criterion is met for all three
datasets. Figure agreement remains dependent on the central summary, as shown
separately below.

The box fits use two concurrent chains, 500 warmup iterations and 1,000 retained
draws per chain. The accepted Camu run uses the source's four chains, 1,000 warmup
iterations and 2,000 retained draws per chain. Each fit runs in one Julia 1.12.6
process, with target acceptance 0.995, maximum tree depth 12, and seed 100785.
Tsit5 uses absolute/relative tolerances 1e-6 and a 100,000-step cap. The box sample
counts, target acceptance, integrator, step cap, and initialization differ from
the original driver; the source equations, observations, likelihood, priors,
and acceptance thresholds are unchanged. Each run saves its exact settings
and input/implementation hashes.

To repeat the two-chain configuration in a new batch directory, or the longer
Camu run in its own directory:

```sh
julia --startup-file=no --check-bounds=yes --threads=2 --project scripts/fit_cocoa.jl --trial=all --method=bayes --samples=1000 --warmup=500 --chains=2 --chain_execution=threads --target_accept=0.995 --solver_maxiters=100000 --output=data/processed/cocoa/moreno2018/fits/bayes_check_new
julia --startup-file=no --check-bounds=yes --threads=4 --project scripts/fit_cocoa.jl --trial=camu --method=bayes --samples=2000 --warmup=1000 --chains=4 --chain_execution=threads --target_accept=0.995 --solver_maxiters=100000 --output=data/processed/cocoa/moreno2018/fits/camu_full_counts_new
```

| Trial / report | Maximum R-hat | Minimum bulk / tail ESS | Divergences / depth hits | Means inside published intervals | Accepted |
|---|---:|---:|---:|---:|---|
| [Camu](../data/processed/cocoa/moreno2018/fits/camu_bayes_full_counts/report.toml) | 1.00202 | 3430.7 / 2266.4 | 0 / 0 | 25 / 25 | Yes |
| [Box 1](../data/processed/cocoa/moreno2018/fits/box1_bayes_extended/report.toml) | 1.00967 | 603.8 / 488.8 | 0 / 0 | 25 / 25 | Yes |
| [Box 2](../data/processed/cocoa/moreno2018/fits/box2_bayes_extended/report.toml) | 1.00410 | 576.5 / 473.5 | 0 / 0 | 25 / 25 | Yes |

The thresholds are R-hat < 1.01, both ESS >= 100 times the chain count (400 for
Camu, 200 for each box), and zero divergences/depth hits. The four-chain Camu
run took 6,189 seconds including predictive calculations on this machine.

Earlier Camu attempts remain separate: the two-chain run at target acceptance
0.985 had two retained divergences; the two-chain retry at 0.995 failed tail ESS
for mu4 (102.9) and mu5 (190.5). Both remain unaccepted. The subsequent four-chain
run passes all gates on its own; no draws or diagnostics were combined across
runs to claim success.

| Trial / comparison | Maximum mean error | Maximum median error | Points outside 2%, mean / median | Curve acceptance |
|---|---:|---:|---:|---|
| [Camu](../data/processed/cocoa/moreno2018/fits/camu_bayes_full_counts/curve_comparison.toml) | 2.04% | 1.63% | 1 / 0 | Median passes; mean fails |
| [Box 1](../data/processed/cocoa/moreno2018/fits/box1_bayes_extended/curve_comparison.toml) | 1.63% | 0.94% | 0 / 0 | Both summaries pass |
| [Box 2](../data/processed/cocoa/moreno2018/fits/box2_bayes_extended/curve_comparison.toml) | 1.97% | 2.07% | 0 / 1 | Mean passes; median fails |

Each comparison uses 56 digitized central-curve points and the preselected 2%
axis-range tolerance. Camu's mean misses at yeast biomass at 12 h (2.0422%);
box 2's median misses at fructose at 12 h (2.0690%). Neither mean nor median
passes the fixed criterion across all three datasets. Choosing a different
summary for each trial would not identify the source's plotting method.
These are summaries of joint posterior trajectories, not trajectories evaluated
at marginal parameter summaries. The P2.3 expected-failure regression remains
unchanged for simulations at printed parameter means.

The latest checkpoint audits compare 200,000 Camu parameter values and 50,000
values per box against final `draws.csv` by parameter name, chain, and draw
index. They match exactly, as do the divergence/depth counts in the reports.
Predictive plots and all-parameter trace plots are saved beside each report.

| Trial | Predictive intervals | Parameter traces |
|---|---|---|
| Camu | [PNG](../data/processed/cocoa/moreno2018/fits/camu_bayes_full_counts/predictive.png) | [PNG](../data/processed/cocoa/moreno2018/fits/camu_bayes_full_counts/traces.png) |
| Box 1 | [PNG](../data/processed/cocoa/moreno2018/fits/box1_bayes_extended/predictive.png) | [PNG](../data/processed/cocoa/moreno2018/fits/box1_bayes_extended/traces.png) |
| Box 2 | [PNG](../data/processed/cocoa/moreno2018/fits/box2_bayes_extended/predictive.png) | [PNG](../data/processed/cocoa/moreno2018/fits/box2_bayes_extended/traces.png) |

The earlier tighter-tolerance Camu MLE reported native optimizer convergence with normalized negative log likelihood
-75.2937 and sigma 0.13437. Several coefficients reached the explicit search
bounds, so local intervals were withheld. Twelve of 25 estimates lie outside
the published Bayesian intervals. MLE and the source's regularized posterior
are different estimation targets; neither the priors nor bounds were changed
to conceal this difference.

Box MLE attempts remain unaccepted. A Rodas5P box-1 restart from its saved point
improved tentative NLL from -65.7352 to -65.8189, with several coefficients
approaching numerical bounds. It returned `MaxIters` and an explicit time-limit
flag after 189 seconds. For box 2, the first stiff-solver start reached tentative
NLL -60.1060 and optimizer limits; the second hit the 10,000-step ODE cap.
Per-start files preserve those outcomes. A stiff solver covers difficult regimes
but does not guarantee optimization convergence across the broad search bounds.

Subsequent one-start restarts used the same bounds, Rodas5P, 1e-9/1e-8 tolerances,
10,000-step ODE cap, 1,000 inner iterations, and 20 outer iterations, with a
600-second soft budget. Both used the earlier acceptance rule without the
independent stationarity gate; their original fit flags are preserved.

| Trial / checkpoint | NLL before → after | Native status / elapsed | Projected gradient before → after | Stationary at 1e-6 |
|---|---:|---|---:|---|
| [Box 1](../data/processed/cocoa/moreno2018/fits/box1_mle_restart_long/mle_start_001.toml) | -65.818909 → -65.818944 | Success / 500.5 s | 0.12344 → 0.012231 | No |
| [Box 2](../data/processed/cocoa/moreno2018/fits/box2_mle_restart_long/mle_start_001.toml) | -60.106032 → -60.106332 | MaxIters, time limit / 827.5 s | 0.0037993 → 0.000080134 | No |

The new gate would reject both saved points. Box 1 illustrates why native
success alone is insufficient: its `report.toml` records `converged=true` under
the earlier rule, but its independently recomputed gradient exceeds the requested
tolerance. Its reproduction flag is already false, with 15 of 25 estimates
outside the published Bayesian intervals. No historical report was relabeled.
Box 2's time overrun reflects the soft optimizer budget, not a longer requested
limit.

The [Box 1](../data/processed/cocoa/moreno2018/fits/box1_mle_diagnostics_after/report.toml)
and [Box 2](../data/processed/cocoa/moreno2018/fits/box2_mle_diagnostics_after/report.toml)
local information factors have numerical rank 25 at relative cutoff
`sqrt(eps(Float64))`, with condition numbers 6.61e5 and 7.71e6. The smallest
singular values are 9.86e-5 and 8.45e-6. These measurements indicate weak local
directions, not structural non-identifiability. Six box-1 and eight box-2
parameters lie within the stated bound-proximity threshold. Local confidence
intervals are not justified by these diagnostics.

The [Box 1 numerical scan](../data/processed/cocoa/moreno2018/fits/box1_mle_numerics/report.toml)
and [Box 2 numerical scan](../data/processed/cocoa/moreno2018/fits/box2_mle_numerics/report.toml)
hold those parameter points fixed:

| Trial | Integrator / absolute, relative tolerance | Projected gradient | FD disagreement at log step 1e-6 |
|---|---|---:|---:|
| Box 1 | Rodas5P / 1e-9, 1e-8 | 0.01223 | 0.04155 |
| Box 1 | Rodas5P / 1e-11, 1e-10 | 0.06299 | 0.003357 |
| Box 1 | Tsit5 / 1e-9, 1e-8 | 1.25584 | 0.26567 |
| Box 1 | Tsit5 / 1e-11, 1e-10 | 0.04437 | 0.01077 |
| Box 2 | Rodas5P / 1e-9, 1e-8 | 0.00008013 | 0.0004812 |
| Box 2 | Rodas5P / 1e-11, 1e-10 | 0.0009626 | 0.00001563 |
| Box 2 | Tsit5 / both pairs | ODE step cap | Unavailable |

FD disagreement is `max(abs(FD-AD)/max(1,abs(AD)))` over all 25 log derivatives.
For box 2, the tighter Rodas5P run agrees best at step 1e-5 (1.42e-6 disagreement).
For box 1, even the tighter solver gradients differ appreciably despite almost
identical NLL values. All completed configurations fail stationarity. These
results support further numerical refinement before claiming an MLE; a larger
time budget alone does not resolve the derivative-accuracy issue. They do not
establish structural non-identifiability.

The Bayesian parameter criterion now passes on all datasets under the recorded
configurations. Verified MLE stationarity for the boxes, MLE agreement with Bayesian
reference intervals, and a common central-curve reproduction remain unresolved.
The older Camu MLE predates the independent stationarity gate and has not been
reclassified under it.
The original RStan numerical configuration has not been run. The source's
plotted central summary is still unidentified. No hosted-CI success is claimed.

## Validation

Generic estimated-noise tests match an analytical Gaussian MLE, expected Fisher
information, and profile likelihood. A half-Cauchy noise posterior is checked
against independent one-dimensional numerical quadrature and held-out predictive
coverage. Cocoa checks verify normalization, exclusion of initial observations,
state-major ordering, full likelihood normalization, and the exact source priors.
The baseline tests separately check equations, analytic limits, positivity,
physical/normalized coordinate equivalence, and derivatives through the ODE.
Stiff-solver checks compare predictions and all 25 likelihood derivatives
against Tsit5 on each original dataset, plus an analytic fast-mortality case.
Execution checks compare seeded serial/threaded draws, exclude warmup from
checkpoints, preserve failures and time-limit flags, and verify every CSV
column against named parameter values in both loader and equation order.
The full suites on Julia 1.10.11 and 1.12.6 each pass 591 checks with two threads
and normal compiled modules. Each suite retains one explicit expected failure for the
published-curve gap. The time-budget tests use a delayed nonlinear prediction
to verify that a stopped optimization returns no accepted estimate; they do not
assert a fragile wall-clock upper bound.
Saved-start checks also verify exact named-parameter recovery and rejection of
incomplete parameter tables or starting values outside the original bounds.
The 28 added stationarity/diagnostic checks cover analytical interior and
constrained optima, failure to meet the requested gradient tolerance,
prediction failures, shared-noise information, redundant parameters, missing
observations, and parameter-name ordering. Julia 1.10's `Pkg.test()` took
6m03s; Julia 1.12's full `test/runtests.jl` took 7m27s on this machine.
An additional 96-check artifact audit on Julia 1.12 verifies checkpoint and raw
hashes, all named parameter values, bound labels, gradient norms, numerical-scan
coverage, and the saved finite-difference error calculations. It also confirms
that the scan's implementation hashes match the current files.

A separate 14-check analysis run on Julia 1.12.6 verifies the marginal-summary
diagnostic: all six trial/summary combinations and 336 comparison rows are
present, the mean results reproduce the existing baseline report to 1e-12,
and the original tolerance and unidentified source-summary status are preserved.
The median-parameter trajectories also fail the 2% criterion in all three trials;
see the [baseline discussion](cocoa_baseline.md).

The observation solver saves only the requested times, preserving state-major
ordering without retaining every integration step. A local 50-gradient benchmark
on Julia 1.12.6 reduced allocations from 76.4 MB to 12.9 MB and elapsed time from
0.19 s to 0.08 s after compilation. Prediction and gradient differences from the
dense-interpolation path were below 5e-16. These are local measurements, not a
portable performance guarantee; regression tests cover output equivalence and
derivatives.
