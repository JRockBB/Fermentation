# P2.5 baseline fit attempts

`mle_pilot/` is the first sequential numerical survey of all 32 available
trials, including the 23 retained in the paper. It was run on Julia 1.10.11
on 2026-10-20:

```sh
julia --startup-file=no --check-bounds=yes --project scripts/fit_cocoa_benchmark.jl --output=data/processed/cocoa_benchmark/fits/mle_pilot --scope=all --starts=1 --maxiters=60 --outer_iterations=2 --maxtime=30
```

All 32 runs completed with finite candidate predictions. **None is an
accepted MLE.** Every optimizer returned `MaxIters` and every projected
log-gradient exceeded the unchanged `1e-6` stationarity tolerance. One run,
`ghhp1`, additionally reached the soft optimizer time limit. The short
budget was an initial survey, not a convergence claim.

The sum of reported per-trial elapsed times is 349.018 seconds; startup,
tests, and some compilation are outside those timers. Training normalized
RMSE spans 0.06547–0.35059 across these unverified candidates. This range
does not establish a ranking of model quality across differently sampled
datasets.

[comparison.csv](mle_pilot/comparison.csv) contains all trials and keeps
candidate training NLL/RMSE separate from published Bayesian convergence
statistics and rounded BMA weights. Published PSIS-LOO entries remain empty.
Each trial directory contains the source identity, a complete start
checkpoint, a report, and normalized predictions. Root settings include
implementation hashes and all numerical controls.

An independent artifact audit passed 360 checks: source and implementation
hashes, initial conditions, complete trial coverage, saved observations,
parameter identities, candidate labels, and NLL/RMSE recomputed from the
saved prediction CSVs. No intervals or biological parameter estimates are
accepted from this survey.

## Numerical diagnostics and restarts

`mle_pilot_diagnostics/` inspects all 32 pilot candidates at the fitted
Rodas5P tolerances (`1e-9` absolute, `1e-8` relative) and 100 times tighter
tolerances. All 64 cases completed, including centered finite differences
at log-parameter steps `1e-4` and `5e-5`. At the original tolerances,
14 candidates have local information rank below 25 and 15 have at least
one parameter near a numerical search bound. These are local diagnostics,
not evidence that a particular biological mechanism is absent.

The [diagnostic table](mle_pilot_diagnostics/comparison.csv) and per-parameter
CSVs retain the actual errors, singular values, and bound distances. The
largest half-step scaled gradient discrepancy at the pilot tolerances is
0.004647 (`ecpb2`); finite-difference error can reflect both integration
error and finite step size. These scans do not certify the `1e-6`
stationarity threshold for every parameter point.

`mle_refinement1/` restarts `ghhp6`, `brwb3`, and `ghhp2` from their pilot
candidates with 1,000 inner iterations, 20 barrier updates, a soft
180-second budget per start, and Rodas5P tolerances `1e-11`/`1e-10`.
The 10,000-step ODE cap and `1e-6` stationarity gate remain unchanged.
Parent file hashes and any small inward adjustment at a bound are recorded
in each new input file. The pilot files are preserved.

| Trial | Pilot NLL | Refined NLL | Pilot normalized RMSE | Refined normalized RMSE | Refined projected gradient |
|---|---:|---:|---:|---:|---:|
| `ghhp6` | 11.8657 | unavailable | 0.35059 | unavailable | unavailable |
| `brwb3` | -39.4259 | -56.2145 | 0.16039 | 0.13473 | 0.03211 |
| `ghhp2` | -50.8775 | -69.8610 | 0.14903 | 0.12440 | 1.73471 |

The [refinement table](mle_refinement1/comparison.csv) retains all three
outcomes. `ghhp6` encountered an ODE step-limit failure during its search;
the attempt produced no new retained candidate. Its original pilot point
is still available. Both other trials improved their training scores but
reached the optimizer time limit and failed the stationarity gate.
**None of these refinements is an accepted MLE.**

The finite refined candidates are inspected separately in
`mle_refinement1_diagnostics/`, at their fitted tolerances and another
100 times tighter. The failed attempt has no new parameter point to
inspect, so that diagnostic batch explicitly selects `brwb3,ghhp2`.
These training improvements remain distinct from the paper's posterior
predictive comparison and do not supply biological estimates or intervals.
The refined projected gradients change by less than `1.1e-8` under the
additional tolerance tightening; both remain well above `1e-6`.

An expanded audit passes **1,111 checks** across the pilot, refinement, and
both diagnostic batches. It verifies parent-file hashes, source identity,
saved observations and parameters, independently reconstructed NLL/RMSE,
reported gradient discrepancies, and the explicit failed-attempt records.
The historical pilot runner hash remains historical; its unchanged
scientific implementation hashes are checked against current files.

## Bayesian benchmark run

`bayes_brwb3_ghhp4_v2/` exercises both source families: `brwb3` has eight
measured series and 96 observed state-time cells; `ghhp4` has seven measured
series and 105 cells, with lactate represented only as a latent state.
The run uses the 2022 source priors, two threaded chains, 500 warmup and
1,000 retained draws per chain, target acceptance 0.995, tree depth 12,
and Tsit5 with absolute/relative tolerances `1e-6`. Per-trial source ODE caps
are 10,000 and 1,000,000 respectively. The scope is `all`; seeds are 141100
and 141115, preserving the archived trial order.

```sh
julia --startup-file=no --check-bounds=yes --threads=2 --project scripts/fit_cocoa_benchmark_bayes.jl --output=data/processed/cocoa_benchmark/fits/new_bayes_run --scope=all --trials=brwb3,ghhp4 --samples=1000 --warmup=500 --chains=2 --target_accept=0.995 --max_depth=12 --predictive_draws=500 --chain_execution=threads
```

The [settings](bayes_brwb3_ghhp4_v2/settings.toml) record source identities,
implementation hashes, prior initialization, and differences from the
archived RStan workflow. Each completed trial exports every retained draw,
all pointwise Gaussian log likelihoods, and predictive summaries from 500
selected joint draws. Initial conditions and unmeasured lactate are excluded
from the likelihood. Hidden lactate has only normalized latent bands; no
physical scale or measurement band is assigned to it.

Both trials completed on Julia 1.10.11 on 2026-10-20 and exported all 2,000
retained draws per trial. **Neither passes the sampling diagnostics.**

| Trial | Maximum R-hat | Minimum bulk ESS | Minimum tail ESS | Divergences | Depth hits | Reported elapsed seconds |
|---|---:|---:|---:|---:|---:|---:|
| `brwb3` | 1.01915 | 258.47 | 353.19 | 2 | 0 | 1315.28 |
| `ghhp4` | 1.00980 | 301.76 | 98.64 | 0 | 0 | 695.37 |

The unchanged gate requires R-hat below 1.01, both ESS measures at least 200
for these two-chain runs, and zero divergences/depth hits. `brwb3` fails
R-hat and divergence checks; `ghhp4` fails tail ESS for `yc11`. The largest
`brwb3` R-hat values occur in `mu5`, `yc11`, and `k3`. These results identify
where further sampling diagnostics are needed; the parameter summaries and
bands remain unaccepted. No biological conclusion is drawn from a failed
sampling check.

The [Bayesian comparison table](bayes_brwb3_ghhp4_v2/comparison.csv) retains
both outcomes alongside the published baseline diagnostics and BMA weights.
The source reports successful fits for both trials, using a different
sampling workflow. Our pointwise likelihood exports are complete, but
computed and published PSIS-LOO entries remain empty.

An independent audit passes **475 checks** over these files. It recomputes
all marginal summaries and sampling diagnostics from the checkpoints,
reconstructs the latent bands from all 500 recorded predictive selections
per trial, and replays selected pointwise likelihoods through the ODE and
an independent Gaussian log-density calculation. It also checks all 402,000
exported likelihood values for finiteness and complete observation/draw
mapping, physical units, implementation/source hashes, comparison-table
consistency, and the interrupted attempt's original source snapshot.

Predictive and trace figures are saved as PNG and SVG, with source hashes
in each trial's `figure_sources.toml`. All four PNGs were visually reviewed;
failed diagnostic status is prominent, and hidden lactate has no physical
or observation band.

| Trial | Predictive panels | Parameter traces |
|---|---|---|
| `brwb3` | [PNG](bayes_brwb3_ghhp4_v2/brwb3/predictive.png) | [PNG](bayes_brwb3_ghhp4_v2/brwb3/traces.png) |
| `ghhp4` | [PNG](bayes_brwb3_ghhp4_v2/ghhp4/predictive.png) | [PNG](bayes_brwb3_ghhp4_v2/ghhp4/traces.png) |

`bayes_brwb3_ghhp4/` preserves an earlier attempt interrupted during sampling
to correct hidden-lactate export semantics. Its
[interruption record](bayes_brwb3_ghhp4/aborted_run.toml) labels it unaccepted,
and its original runner snapshot matches the implementation hash in that
attempt's settings. Its partial checkpoints are not a completed fit.

Preserve this directory when running a new experiment. See
[the benchmark workflow](../../../../docs/cocoa_benchmark.md) for the
remaining convergence and predictive-comparison work.

## Four-chain `ghhp4` refinement

`bayes_ghhp4_refinement1/` increases sampling to four chains, 1,000 warmup
and 2,000 retained draws per chain. The
[experiment record](bayes_ghhp4_refinement1.experiment.toml) was written
before sampling; it records the failed parent fit, parent file hashes,
environment hashes, changes, and the exact command.

This is a fresh fit with new adaptation, not a checkpoint resume or pooling
of earlier draws. The equations, data, priors, trial seed, Tsit5 tolerances,
ODE cap, target acceptance 0.995, and depth limit 12 are unchanged. With
four chains, the existing ESS rule requires both bulk and tail ESS of at
least 400; R-hat must remain below 1.01, with zero divergences and depth
hits.

The run completed on Julia 1.10.11 on 2026-10-20 and **passes the unchanged
sampling diagnostic rule**. All 8,000 retained draws and all 840,000
pointwise likelihood values are exported.

| Run | Retained draws | Maximum R-hat | Minimum bulk ESS | Minimum tail ESS | Divergences | Depth hits |
|---|---:|---:|---:|---:|---:|---:|
| Two-chain pilot | 2,000 | 1.00980 | 301.76 | 98.64 | 0 | 0 |
| Four-chain refinement | 8,000 | 1.00226 | 3219.46 | 1920.08 | 0 | 0 |

Reported fitting and export time is 3,917.48 seconds. The
[comparison table](bayes_ghhp4_refinement1/comparison.csv) includes the
published, rounded baseline diagnostics. Passing these sampling checks
does not establish model validity or reproduce published predictive scores.
Computed and published PSIS-LOO entries remain empty.

An independent **282-check audit** verifies the experiment and environment
hashes, unchanged inputs and priors, all saved marginal diagnostics,
complete draw/observation mapping, and all 840,000 finite likelihood
values. It replays selected likelihoods through the ODE and reconstructs
latent bands from all 500 recorded predictive selections. The original
`yc11` tail-ESS failure is absent in this run: its tail ESS is 2261.76.

The [predictive panels](bayes_ghhp4_refinement1/ghhp4/predictive.png) and
[four-chain traces](bayes_ghhp4_refinement1/ghhp4/traces.png) were visually
reviewed. PNG and SVG versions are saved with input/script hashes;
15 figure integrity checks pass. Lactate remains a normalized latent
state without a physical scale or measurement band.

## Saved posterior numerical checks

`bayes_pilot_numerics1/` checks 14 selected `brwb3` points and 10 selected
`ghhp4` points from the original two-chain pilot. Selection includes a
three-point grid and `yc11` extrema per chain, both flagged divergent
retained draws, and their preceding retained states. All **72 cases**
complete with the saved Tsit5 tolerances, 100 times tighter Tsit5
tolerances, and Rodas5P at the tighter tolerances. Per-trial source ODE caps
remain in force.

The [numerical comparison](bayes_pilot_numerics1/comparison.csv) records
every point and case. Maximum half-step finite-difference discrepancies
in the log-coordinate posterior gradient are:

| Trial | Tsit5, `1e-6` | Tsit5, `1e-8` | Rodas5P, `1e-8` |
|---|---:|---:|---:|
| `brwb3` | 0.004091 | 0.0002966 | 0.00002724 |
| `ghhp4` | 0.0004399 | 0.00001014 | 0.00008786 |

Errors use `max(abs(FD-AD)/max(1,abs(AD)))`, with a log-coordinate
finite-difference step of `5e-5`; the `1e-4` results and full vectors are
also saved. Across the selected points, tightening or switching integrators
changes NLL by less than `2.88e-5` and normalized predictions by less than
`2.96e-6`. The two tighter integrators' NLLs differ by at most `8.92e-6`
for `brwb3` and `5.06e-7` for `ghhp4`.

These results motivated tighter integration controls in the subsequent
`brwb3` refinements. Rodas5P has the smaller worst gradient discrepancy for
that subset, while the tighter Tsit5 result is better for the selected
`ghhp4` points. This does not identify the cause of a divergence: the failed
leapfrog state is not retained, and selected points do not certify the
entire posterior. Neither pilot's acceptance label changes.

An independent **1,709-check audit** verifies parent and implementation
hashes, selection coverage and reasons, parameter identities and ordering,
the prior/Jacobian contribution to each density, finite-difference errors,
and every reported cross-case discrepancy.

## `brwb3` solver cost and interrupted experiment

The [first refinement experiment](bayes_brwb3_refinement1.experiment.toml)
specified four chains, 1,000 warmup and 2,000 retained draws per chain,
Rodas5P at `abstol = reltol = 1e-8`, the source ODE cap of 10,000,
target acceptance 0.995, and depth limit 12. It was stopped for computational
cost before any retained draw, after **2,272.25 seconds of process elapsed
time** and 6,643.05 CPU seconds across threads. These times include startup.
All four chain files contain only their headers. The
[interruption record](bayes_brwb3_refinement1/aborted_run.toml) marks the
attempt unaccepted; the last sampling report is an unfinished checkpoint.
No adaptation state was saved, so it cannot be resumed or treated as a fit.

A [solver timing comparison](bayes_brwb3_solver_timing1/summary.csv) then
evaluated posterior gradients at 10 points from the completed pilot: one
middle draw per chain, `yc11` extrema, divergent retained draws, and their
predecessors, with duplicates removed. It uses the same log-coordinate
posterior as the numerical checks, including priors and the Jacobian.
Five evaluations per point and case give **150 timed gradients** after
compilation warmup. The three cases alternate within each repetition.

| Solver and tolerance | Median gradient time |
|---|---:|
| Tsit5, `1e-6` | 3.8300 ms |
| Tsit5, `1e-8` | 4.6412 ms |
| Rodas5P, `1e-8` | 38.4006 ms |

Tighter Tsit5 is **8.27 times faster** than tighter Rodas5P on this gradient
benchmark. These are measurements under shared machine load with one Julia
thread, not end-to-end NUTS timings or a numerical accuracy certificate.
The directory retains every timing, selected parameter point and gradient,
parent/implementation/environment hashes, and the exact timing script.
The results support using tighter Tsit5 for the fresh
[second refinement experiment](bayes_brwb3_refinement2.experiment.toml),
with the same source priors and unchanged diagnostic rule.

The second experiment requested four chains with 1,000 warmup and 2,000
retained draws per chain, using Tsit5 at `abstol = reltol = 1e-8`. It was
**interrupted and remains unaccepted** on 2026-10-20 after two retained
divergences had already violated the unchanged zero-divergence rule.
The remaining computation was stopped after checkpoint progress ceased.
The [interruption record](bayes_brwb3_refinement2/aborted_run.toml) preserves
the reason, observed counts, runtime, and diagnostic limitations.

| Chain | Retained draws saved | Divergent retained draw |
|---|---:|---:|
| 1 | 2,000 | 1,343 |
| 2 | 1,478 | None observed |
| 3 | 2,000 | 1,391 |
| 4 | 0 | None observed |

There are **5,478 partial retained draws**, with two observed divergences
and no observed depth-limit hits. Process elapsed time was 7,268.89 seconds,
with 17,034.83 CPU seconds across threads. These timings include startup,
shared machine load, and brief native thread inspection. The inspection
did not establish the cause of the delay. An isolated diagnostic using
the fourth chain's seed and the same numerical controls completed its
first ten warmup steps in 197.46 seconds, including step compilation;
its [script snapshot](bayes_brwb3_refinement2/runtime_diagnostic/chain4_initial_probe.jl)
and [log](bayes_brwb3_refinement2/runtime_diagnostic/chain4_initial_probe.log)
are retained. Those steps are not pooled with the experiment.

The original `report.toml` remains an unfinished sampling checkpoint;
`aborted_run.toml` supplies the terminal status. No completed joint R-hat/ESS,
posterior summary, likelihood export, or predictive bands are claimed.
The completed-fit auditor rejects this directory. Tighter tolerances and
the longer requested run did not eliminate retained divergences, so further
work must address the sampling and numerical behavior before accepting a fit.

A separate **379-check audit** verifies both interrupted experiments,
their preregistrations, solver timing records, and the partial checkpoint
identities and divergence flags. The saved
[interruption audit](bayes_brwb3_refinement2.interruption-audit.toml) records
file hashes and explicitly marks the experiment incomplete and unaccepted.

In the experiments above, "source cap" means the source's integer value
passed to Julia's whole-trajectory `maxiters`. Stan counts its step limit
between output times. This difference is documented in the
[solver-control comparison](../../../../docs/cocoa_benchmark.md#bayesian-baseline-workflow)
and remains a limitation of those historical fits; saved controls and acceptance
labels are unchanged. The later interval-budget experiments below use an
explicitly recorded scope.

## Complete saved-fit replay audits

The reusable [artifact auditor](../../../../scripts/audit_cocoa_benchmark_bayes.jl)
was run with `replay_per_chain=2000`, so every retained draw in these completed
experiments was replayed through its saved ODE controls. Gaussian log densities
were recomputed independently of the likelihood export function. Both audit
outputs pass artifact integrity:

| Experiment and audit | Draws replayed | Likelihood values checked | Sampling diagnostics |
|---|---:|---:|---|
| [Two-trial pilot](bayes_brwb3_ghhp4_v2.audit.toml) | 4,000 | 402,000 | Both fail |
| [`ghhp4` refinement](bayes_ghhp4_refinement1.audit.toml) | 8,000 | 840,000 | Pass |

The audits also reconstruct latent predictive bands from all 500 recorded
selections per trial: 1,500 selected predictions in total. They recompute
marginal diagnostics, verify source identities and physical conversions,
and confirm that parent files remain unchanged. The pilot and refinement
audits record 27,591 and 49,857 integrity checks respectively. These are
artifact checks, separate from the package's unit-test count.

Measurement-noise draws were not retained, so those bands receive only
structural and unit checks. Integrity does not establish convergence or
biological validity. The failed pilots remain unaccepted, and no PSIS-LOO
comparison is inferred from these replays.

## Output-interval solver budgets

The [budget comparison](bayes_brwb3_solver_budget1/settings.toml) reuses the ten
`brwb3` points selected in the earlier gradient timing experiment, including
divergent retained states, their predecessors, and `yc11` extrema. Its
[script](bayes_brwb3_solver_budget1/experiment.jl) verifies source/prior identity,
each point against its original checkpoint, parent hashes, the previous
implementation snapshots, and the environment hashes before comparing modes.
It does not bypass the completed-fit reader's implementation checks.

All 30 cases give bit-identical predictions, Gaussian NLLs, and log-coordinate
posterior gradients between whole-solve and output-interval budgeting. The
60 gradient evaluations also agree with the saved timing gradients at the
declared numerical tolerance. Each case checks its exact accepted-step boundary:
the maximum interval count succeeds, and one fewer fails. The comparison passes
302 case/integrity checks, in addition to its provenance preflight.

| Solver/tolerances | Points | Largest interval accepted-step count | Largest whole-trajectory attempted-step count |
|---|---:|---:|---:|
| Tsit5, `1e-6` | 10 | 373 | 1,022 |
| Tsit5, `1e-8` | 10 | 615 | 1,450 |
| Rodas5P, `1e-8` | 10 | 203 | 1,695 |

Both modes therefore remain below the source integer cap of 10,000 at these
points. This corrects budget accounting without demonstrating that the mismatch
caused earlier divergences. Retained states do not identify the failed leapfrog
states, and this subset is not a global numerical certificate. The two historical
failed pilots and accepted `ghhp4` refinement retain their original status.
Their saved replay audits predate this implementation change and remain tied to
their recorded code hashes.

## Paired runtime diagnostic

The [declared pilot](bayes_brwb3_budget_pilot1/experiment.toml) compares both budget
modes using the same four chain seeds, Tsit5 tolerances `1e-8`, target acceptance
`0.995`, 50 discarded calls, and 25 retained draws per chain. Depth 8 and a
600-second soft deadline per case bound this diagnostic. The process uses two
Julia threads. These short adaptation/sample counts and the depth limit are
computational choices, not the source's full fitting protocol.

The [driver](bayes_brwb3_budget_pilot1/experiment.jl) records its controls before
sampling, writes initialization/warmup/sampling completion records, audits every
likelihood row of any completed case, and compares corresponding checkpoint
rows. Each case starts fresh; no earlier draws are pooled. Compilation and
scheduling differ with case order, so the elapsed times are not a speed benchmark.

Both original cases completed with 100 retained draws, identical paired
checkpoints, zero retained divergences, and 80 depth hits. Maximum R-hat was
1.59931 and minimum bulk/tail ESS 9.84/14.44; both remain unaccepted. Their
scientific export audits each pass 1,449 checks, replaying every likelihood
row and all 50 selected latent predictions per case.

An additional independent runtime check then rejected the logs: in both cases,
chain 3 draw 1 has checkpoint depth 7 but runtime depth 8. See the
[discrepancy record](bayes_brwb3_budget_pilot1.runtime-discrepancy.toml) and
[reproducing check](bayes_brwb3_budget_pilot1.runtime-check.py). The completed-fit
audits checked scientific exports and aggregate flags, not correspondence of
each runtime row. These files are preserved, and their runtime logs are not
certified. The [exact callback source](bayes_brwb3_budget_pilot1.callback-source.jl)
is saved with its hash verified against the original settings.

The callback shared a boxed `diagnostics` variable with the final fit summary.
Another chain could overwrite it while checkpoint I/O yielded under the lock.
The fix gives each callback its own transition diagnostics in a separate scope;
a channel-gated task test forces this overlap even on a single Julia thread.
The earlier callback had the same shared-variable pattern, so historical
threaded per-draw diagnostic fields cannot be certified against unavailable raw
sampler transitions. Parameter draws and diagnostics calculated directly from
the raw chain are independent of that callback variable. This finding does not
explain the previous long-run delay or change historical fit acceptance labels.

The [corrected repeat](bayes_brwb3_budget_pilot2/experiment.toml) uses the same
declared controls and fresh initialization. Both cases again complete, with all
100 paired checkpoint rows identical. The independent
[runtime audit](bayes_brwb3_budget_pilot2.runtime-audit.toml), using a saved
[standard-library Python checker](bayes_brwb3_budget_pilot2.runtime-audit.py),
verifies all 600 initialization/warmup/sampling events, all 200 retained draws,
per-row diagnostic correspondence, phase boundaries, monotonic timings, matching
exports, implementation hashes, and unchanged parent files. No warmup row enters
the retained checkpoints.
All eight corrected checkpoint files also match the pre-fix pilot byte for byte;
the observed discrepancy was in its runtime log, and the corrected repeat does
not change those saved draws.

| Corrected case | Retained draws | Divergences | Depth hits | Maximum R-hat | Minimum bulk/tail ESS | Sampling accepted |
|---|---:|---:|---:|---:|---:|---|
| Whole solve | 100 | 0 | 80 | 1.59931 | 9.84 / 14.44 | No |
| Output interval | 100 | 0 | 80 | 1.59931 | 9.84 / 14.44 | No |

The corrected [whole-solve audit](bayes_brwb3_budget_pilot2/solve.audit.toml) and
[interval-budget audit](bayes_brwb3_budget_pilot2/output_interval.audit.toml)
each pass 1,449 scientific-export checks: 200 draws and 19,200 pointwise
likelihood values replayed across the pair, plus all 100 selected latent
predictions. Another 24 checks compare checkpoint counts and runtime ordering.
These are separate from the package test count.

In the corrected interval case's [runtime records](bayes_brwb3_budget_pilot2/output_interval/brwb3/runtime.csv),
the per-chain median gap between warmup completions is 1.07–1.15 seconds
(initialization excluded); the longest gap is 7.41 seconds. Retained-call
medians are 0.96–1.11 seconds, with a longest gap of 8.02 seconds. Slow calls
therefore occur during retained sampling too. These are elapsed completion
gaps including scheduling and compilation, not isolated integration timings.

Thus the accounting change leaves this short sampling run unchanged. It does
not resolve the earlier `brwb3` divergence/convergence problem. The next sampling
investigation should use the now reliable runtime records to examine costly
warmup and posterior exploration before another long fit.

## Descriptive parameter correlations

A [provenance-checked scan](bayes_brwb3_geometry1/settings.json) summarizes all
2,000 retained `brwb3` draws from the original unaccepted pilot. It verifies
checkpoint hashes against the earlier full replay audit, source/prior hashes,
and draw order. It performs no model replay or refit. The saved
[script](bayes_brwb3_geometry1/experiment.py) computes all 300 log-coordinate
parameter-pair correlations within each chain and pooled (900 values), plus
sample log-coordinate scales. The [ranking](bayes_brwb3_geometry1/ranked_pairs.csv)
uses the smaller absolute within-chain correlation, with signs retained.

| Pair | Chain 1 Pearson correlation of log values | Chain 2 |
|---|---:|---:|
| `mu5`, `yc11` | -0.673 | -0.661 |
| `mu3`, `ks3` | 0.655 | 0.660 |
| `mu5`, `yc8` | -0.577 | -0.609 |

These are descriptive features of an unaccepted pilot, useful for choosing
further posterior diagnostics. They are not population posterior estimates,
an identifiability result, or evidence that correlation caused the divergences
or slow warmup. No parameter file, prior, equation, or acceptance label changes.

## Extended runtime and saved-state cost profiles

The [telemetry repeat](bayes_brwb3_runtime_pilot1/experiment.toml) uses the same
four seeds and computational controls as the corrected interval-budget pilot:
50 discarded calls and 25 retained draws per chain, depth 8, target acceptance
0.995, Tsit5 tolerances `1e-8`, the source cap of 10,000 accepted steps between
outputs, and two Julia threads. The [driver](bayes_brwb3_runtime_pilot1/experiment.jl)
records the experiment before sampling and preserves the exact previous callback
and runner sources after checking their hashes against the reference settings.
No prior, equation, observation, or acceptance gate changes.

The new runtime schema copies leapfrog counts, step sizes, acceptance rates,
and Hamiltonian energy statistics from each completed sampler transition.
Initialization has missing statistics. Retained statistics are also exported
independently from the raw sampler chain. All nine scientific CSV exports,
including the four checkpoints, are byte-identical to the reference pilot.
The [scientific audit](bayes_brwb3_runtime_pilot1/fit.audit.toml) passes 1,449
checks, replaying all 100 draws, 9,600 likelihood values, and 50 selected latent
predictions. The [runtime audit and profile record](bayes_brwb3_runtime_pilot1/runtime_diagnostics/settings.toml)
validates all 300 completion events and exact agreement of every retained
statistic with the raw-chain export. Reference and parent hashes remain unchanged.

The fit remains unaccepted: maximum R-hat 1.59931, minimum bulk/tail ESS
9.84/14.44, zero retained divergences, and 80 depth hits. Six divergences occur
in discarded warmup. The short adaptation and depth limit are diagnostic
controls, not the published fitting protocol. Chain 4 uses 255 leapfrog steps
on every retained call. Per-chain median retained completion gaps range from
0.855 to 1.156 seconds; the longest is chain 1 draw 2 at 10.140 seconds.
These gaps include scheduling, compilation, callback work, and GC. Initial
completion times also include time waiting for a thread to start a chain.

The deterministic [selection](bayes_brwb3_runtime_pilot1/runtime_diagnostics/brwb3/selection.csv)
contains each chain's slowest retained call, its predecessor, and a
median-duration control: 12 distinct states. Each state receives one unmeasured
warm call followed by three separately recorded primal solves and three
log-coordinate posterior gradient evaluations. All profiles complete and
reproduce the likelihood trajectory. The [cost table](bayes_brwb3_runtime_pilot1/runtime_diagnostics/brwb3/costs.csv)
shows substantial variation in integration work:

| Selected state | Selection reason | Accepted / rejected ODE steps | Median gradient time (ms) |
|---|---|---:|---:|
| Chain 1, draw 2 | Slowest call in chain 1 | 147 / 1 | 2.70 |
| Chain 2, draw 14 | Slowest call in chain 2 | 329 / 117 | 5.57 |
| Chain 3, draw 10 | Slowest call in chain 3 | 141 / 0 | 3.97 |
| Chain 4, draw 14 | Median-duration control | 134 / 0 | 4.16 |
| Chain 4, draw 24 | Predecessor of slowest call | 1,269 / 1,262 | 37.96 |
| Chain 4, draw 25 | Slowest call in chain 4 | 1,201 / 1,180 | 27.31 |

The final two chain-4 states therefore provide concrete expensive integration
cases. The other ten selected states have median gradient times of 2.46–5.57 ms.
However, chain 1's longest completion gap leaves an inexpensive retained state.
A retained state and its predecessor do not reconstruct the Hamiltonian path;
these measurements cannot establish the cause of a slow transition, a divergence,
or the earlier long-run stall. Warmup parameter states were not retained.
These expensive states motivate the solver comparison below.

## Solver comparison at the costly states

The [declared comparison](bayes_brwb3_runtime_solvers1/experiment.toml) evaluates
all 12 previously selected states with Tsit5 and Rodas5P, each at absolute and
relative tolerances `1e-8` and `1e-10`: 48 cases. The source interval cap of
10,000 accepted steps and separate 500-attempt retry guard stay fixed. Source,
prior, checkpoint, implementation, and environment hashes are checked before
evaluation. The [driver](bayes_brwb3_runtime_solvers1/experiment.jl) and exact
selection are preserved. No states are refitted or pooled into a posterior.

Before evaluation, the experiment records numerical screening thresholds:
maximum normalized prediction difference `1e-6`, absolute NLL difference `1e-5`,
and maximum scaled log-coordinate gradient difference `1e-4` against both
`1e-10` solver references. Centered finite differences at log-coordinate steps
`1e-4` and `5e-5` must each agree with automatic differentiation within `1e-4`
on the same scaled measure, `max(abs(a-b)/max(1,abs(b)))`. The two fine references
must also pass these checks. These are exploratory numerical criteria, not
published parameter values or changes to MCMC acceptance.

Profiles retain primal integration work and three warmed timing repetitions.
Gradient costs are then measured in five separate rounds, rotating solver order
by point and round, with all callables warmed first. Every measurement retains
time, allocation bytes, GC time, and execution position. A Rodas5P candidate
can advance to a bounded pilot only if each expensive state has at least a 20%
median speed improvement and is faster in at least four of five matched rounds,
while the median paired time ratio across all 12 states is no greater than 1.25.
These performance thresholds were also recorded before evaluation.

All 48 cases complete and pass the numerical screen. All Tsit5 `1e-8` baseline
predictions, NLLs, gradients, and integration counts exactly reproduce the parent
profiles. The experiment passes 504 Julia integrity checks. An independent
[Python audit](bayes_brwb3_runtime_solvers1.audit.py) recomputes the comparison
metrics, finite-difference errors, timing summaries, and screening decisions;
its [report](bayes_brwb3_runtime_solvers1.audit.json) passes 100 hash checks and
verifies all 48 cases and 240 timed gradients. Parent artifacts remain unchanged.
The required Julia 1.10.11 `Pkg.test()` run also passes all 1,454 checks, with
the one existing expected cocoa figure failure.

| Solver / tolerances | Median paired gradient-time ratio to Tsit5 `1e-8` | Chain 4 draw 24 (ms) | Chain 4 draw 25 (ms) | Numerical screen |
|---|---:|---:|---:|---|
| Tsit5 / `1e-8` | 1.00 | 22.07 | 18.38 | Pass |
| Tsit5 / `1e-10` | 1.53 | 17.89 | 16.03 | Pass |
| Rodas5P / `1e-8` | 8.36 | 123.83 | 93.75 | Pass |
| Rodas5P / `1e-10` | 18.24 | 147.97 | 138.04 | Pass |

The [full comparison](bayes_brwb3_runtime_solvers1/comparison.csv) shows that
Tsit5 `1e-8` differs from the fine references by at most `3.80e-8` in normalized
prediction, `2.07e-7` in NLL, and `4.55e-6` in scaled gradient on this subset.
The fine references differ by at most `1.39e-6` in scaled gradient. Across all
cases, the largest finite-difference discrepancy is `2.53e-5`.

Rodas5P reduces accepted and rejected step counts at the two expensive states,
but its gradient evaluations are 5.10–5.61 times slower at `1e-8`. Neither
Rodas5P configuration passes the predeclared performance screen. Tighter Tsit5
is modestly faster at those two states, but its median paired cost across the
subset is higher and the existing tolerance passes the numerical checks.
The evidence supports retaining Tsit5 `1e-8`; no alternate-solver pilot is
launched, and no fitting default or acceptance label changes.

These are warmed local costs under shared machine load, not end-to-end sampler
timings. Numerical agreement at selected states is not a global guarantee or
an explanation of earlier divergences and stalls. These results motivate the
bounded longer-adaptation pilot below, using the current solver and runtime records.

## Bounded longer-adaptation pilot

The [declared experiment](bayes_brwb3_adaptation_pilot1/experiment.toml) increases
discarded calls from 50 to 250 and retained draws from 25 to 100 per chain.
Initialization occupies the first discarded call, leaving 249 warmup transitions
per chain. The four seeds, two Julia threads, depth 8, target acceptance 0.995,
Tsit5 tolerances `1e-8`, source interval cap of 10,000 accepted ODE steps, and
500-attempt retry guard remain fixed. Source equations, priors, observations,
units, and acceptance gates are unchanged. Each chain adapts afresh; no draws
are resumed or pooled with earlier fits.

The [driver](bayes_brwb3_adaptation_pilot1/experiment.jl) records controls and
follow-up criteria before sampling. It enables a 1,800-second soft sampling
limit and a stop on any retained divergence. The
[external watchdog](bayes_brwb3_adaptation_pilot1/watchdog.ps1) limits the launched
Julia process to 2,700 seconds, including loading and postprocessing. All four
chains finish their 350 calls and 100 retained draws. No sampling stop or
watchdog termination occurs. The last sampling callback is at 1,139.5 seconds;
the driver's full sampling, export, and audit test set finishes in 20m17.1s.

| Comparison | Retained draws per chain | Maximum R-hat | Minimum bulk / tail ESS | Retained divergences | Depth hits |
|---|---:|---:|---:|---:|---:|
| Reference: 50 discarded calls | 25 | 1.59931 | 9.84 / 14.44 | 0 | 80 / 100 |
| Longer adaptation: full run | 100 | 1.07215 | 47.49 / 39.09 | 0 | 387 / 400 |
| Longer adaptation: first block | 25 | 1.13617 | 42.20 / 27.11 | 0 | 98 / 100 |
| Longer adaptation: last block | 25 | 1.30703 | 14.57 / 31.26 | 0 | 98 / 100 |

The [comparison table](bayes_brwb3_adaptation_pilot1/mixing_comparison.csv) also
records minimum ESS per draw. Bulk ESS per draw rises from 0.0984 to 0.1187,
while tail ESS per draw falls from 0.1444 to 0.0977. Full-run comparisons change
both adaptation and retained length; the R-hat improvement cannot be assigned
to adaptation alone. The short matched blocks are descriptive and do not
establish stationary sampling or reliable posterior estimates.

The [fit report](bayes_brwb3_adaptation_pilot1/fit/brwb3/report.toml) remains
`posterior_diagnostics_failed`: maximum R-hat exceeds 1.01, minimum bulk/tail
ESS is below 400, and the depth-hit fraction increases from 80% to 96.75%.
Zero retained divergences satisfies only one part of the acceptance rule.
Ten divergences occur in discarded warmup. The predeclared follow-up rule
requires improved mixing and a lower depth-hit fraction before advancing to
a longer run; this result does not satisfy it. Next, investigate sensitivity
to the diagnostic depth cap with a bounded comparison and unchanged acceptance
criteria. This result does not establish the cause of earlier divergences or
the interrupted long run.

The [scientific audit](bayes_brwb3_adaptation_pilot1/fit.audit.toml) passes 3,249
checks, replaying all 400 retained draws, 38,400 pointwise likelihood values,
and all 50 selected latent predictions. Measurement-noise draws are unavailable,
so observation bands receive structural and unit checks only.
The [runtime validation](bayes_brwb3_adaptation_pilot1/runtime_diagnostics/settings.toml)
checks all 1,400 completion events and exact agreement of every retained
sampler statistic with the independently exported raw-chain values. Twelve
deterministically selected cost profiles complete. The driver passes another
19 integrity assertions, separately from the scientific and package test counts.
The required Julia 1.10.11 `Pkg.test()` run passes 1,454 checks with the one
existing expected cocoa figure failure (18m44.3s in the test summary), using
two threads and normal compiled modules.

Per-chain median retained completion gaps range from 1.008 to 1.256 seconds;
the longest is 11.406 seconds. All four chains have a median of 255 leapfrog
steps per retained call. Warmup's longest completion gap is 17.622 seconds in
chain 4. Initial completion times include thread-queue delay and must not be
interpreted as isolated initialization costs. The
[saved-state profiles](bayes_brwb3_adaptation_pilot1/runtime_diagnostics/brwb3/costs.csv)
include chain 4 draw 8 with 943 accepted and 766 rejected ODE steps and a
23.24 ms median gradient time. The other selected states range from 2.30 to
11.19 ms. These local replays do not reconstruct a Hamiltonian path or the
unsaved warmup states; elapsed completion gaps also include scheduling and GC.

The [independent audit](bayes_brwb3_adaptation_pilot1.audit.json) verifies 165
hashes, all runtime/checkpoint/raw-statistic correspondences, and reported
comparison ratios. Parent artifacts and implementation/environment identities
remain unchanged. It explicitly records a launcher bookkeeping limitation:
the [process receipt](bayes_brwb3_adaptation_pilot1.process.json) says the child
exited but contains a null native exit code. The original
[checker](bayes_brwb3_adaptation_pilot1.audit.py) stopped at its strict zero-exit
assertion without writing a report. That checker, receipt, and logs are preserved.
The [revised checker](bayes_brwb3_adaptation_pilot1.audit-v2.py) verifies completed
scientific outputs, the final 19-pass driver summary, empty stderr, and log
hashes, while leaving native process success unknown. It does not infer exit 0
or change sampling acceptance. A separate
[launcher probe](bayes_brwb3_adaptation_pilot1.launcher-probe.ps1) retains the
native process handle before waiting and fails if the exit code is missing;
its [report](bayes_brwb3_adaptation_pilot1.launcher-probe.json) verifies both
exit 0 and exit 7 for future launchers. It cannot recover the pilot's exit code.

## Bounded tree-depth comparison

The [declared experiment](bayes_brwb3_depth_pilot1/experiment.toml) compares
depth 9 with the completed depth-8 longer-adaptation pilot above. It keeps
250 discarded calls and 100 retained draws per chain, the four seeds, two
threads, target acceptance 0.995, Tsit5 `1e-8`, source ODE budget, retry guard,
data, equations, and priors fixed. Chains initialize and adapt afresh. The
[driver](bayes_brwb3_depth_pilot1/experiment.jl) and
[corrected watchdog](bayes_brwb3_depth_pilot1/watchdog.ps1) are saved with hashes.
The soft sampling limit increases from 1,800 to 2,700 seconds and the external
process limit from 2,700 to 3,600 seconds to accommodate more sampler work.
The retained-divergence stop remains enabled.

Before sampling, the experiment declares that a longer run can be considered
only after a completed audit, zero retained divergences, a lower cap-hit
fraction, no worse maximum R-hat, and higher minimum bulk and tail ESS per
retained leapfrog. This exploratory comparison does not relax the existing
sampling-acceptance rule. Equal retained lengths separate the depth comparison
from the previous change in sample count. Raising the cap can mechanically
reduce cap hits, so that metric alone is insufficient. ESS per retained
leapfrog excludes warmup and variation in gradient cost; elapsed sampling
time is also reported. One seed set does not establish reproducible efficiency.

All four chains complete their 350 calls and 100 retained draws. No retained
divergence, soft timeout, or watchdog termination occurs. The last sampling
callback is at 2,144.1 seconds; the full driver test set, including export and
audits, takes 38m08.7s. The [process receipt](bayes_brwb3_depth_pilot1.process.json)
records native exit code 0 using the corrected handle-retention launcher.
The reference's null exit code and its independent completion audit remain
preserved and explicitly identified in the new metadata.

| Metric | Depth 8 | Depth 9 |
|---|---:|---:|
| Retained draws | 400 | 400 |
| Maximum R-hat | 1.07215 | 1.05479 |
| Minimum bulk / tail ESS | 47.49 / 39.09 | 96.72 / 74.07 |
| Retained divergences | 0 | 0 |
| Hits at each run's depth cap | 387 (96.75%) | 133 (33.25%) |
| Retained leapfrog steps | 100,996 | 152,048 |
| Minimum bulk ESS per 1,000 retained leapfrogs | 0.4702 | 0.6361 |
| Minimum tail ESS per 1,000 retained leapfrogs | 0.3871 | 0.4871 |
| Elapsed sampling seconds | 1,139.5 | 2,144.1 |
| Minimum bulk / tail ESS per elapsed second | 0.04168 / 0.03431 | 0.04511 / 0.03454 |

The [full comparison](bayes_brwb3_depth_pilot1/depth_comparison.csv) passes every
condition of the predeclared longer-run screen. Bulk and tail ESS per retained
leapfrog improve by about 35% and 26%, respectively, with 51% more retained
leapfrog work. Measured sampling time increases by 88%; bulk/tail ESS per
elapsed second improves by only 8.2%/0.7%. Thus the work-normalized gains do
not establish a corresponding wall-time speedup. Timing includes warmup,
compilation, scheduling, and GC, while gradients vary in cost across states.

The [fit report](bayes_brwb3_depth_pilot1/fit/brwb3/report.toml) remains
`posterior_diagnostics_failed`. Maximum R-hat exceeds 1.01, minimum bulk/tail
ESS is below 400, and depth hits persist. Per-chain depth hits are 1, 15, 30,
and 87, so chain 4 remains the most affected. Seven divergences occur in
discarded warmup, with none in retained sampling. Passing the exploratory
screen supports a separately bounded longer four-chain fit at depth 9; it
does not accept this fit or establish that earlier divergences are resolved.

The [scientific audit](bayes_brwb3_depth_pilot1/fit.audit.toml) passes 3,249 checks,
replaying all 400 draws, 38,400 pointwise likelihood values, and 50 selected
latent predictions. Measurement-noise draws are unavailable, so observation
bands receive structural and unit checks only. The
[runtime validation](bayes_brwb3_depth_pilot1/runtime_diagnostics/settings.toml)
checks all 1,400 events and exact agreement of every retained statistic with
the raw-chain export. Twelve selected state-cost profiles complete. The driver
passes another 19 integrity assertions. The independent
[Python checker](bayes_brwb3_depth_pilot1.audit.py) verifies 168 hashes and
recomputes work, timing, diagnostic ratios, and the screen decision in its
[report](bayes_brwb3_depth_pilot1.audit.json). All parent artifacts, implementation
hashes, and environment hashes remain unchanged.

Chain 4's longest discarded-call completion gap is 21.45 seconds; its median
retained gap is 3.85 seconds, versus 1.20-2.39 seconds for the other chains.
Initialization completion times include waiting for a worker thread. The
[selected profiles](bayes_brwb3_depth_pilot1/runtime_diagnostics/brwb3/costs.csv)
have median gradient times of 4.84-18.44 ms; the most costly selected state,
chain 2 draw 90, requires 484 accepted and 240 rejected ODE steps. These are
different states from the depth-8 profiles. They do not replay discarded
warmup or reconstruct Hamiltonian paths and cannot explain the earlier delays.

This experiment changes computational controls and adds artifacts only. The
package's latest full Julia 1.10.11 suite remains the preceding 1,454-pass run
with one expected cocoa figure failure; it was not repeated for unchanged
production source. The new driver, replay audits, independent comparison
checks, and corrected native exit receipt provide this run's validation.

## Longer depth-9 convergence attempt

The local refinement `bayes_brwb3_depth_refinement1/` follows the passing
exploratory screen above. Four fresh chains each discard 1,000
calls and target 1,000 retained draws, with initialization included in the
discard count. Depth 9, target acceptance 0.995, the same four seeds, two
threads, Tsit5 `1e-8`, the source interval ODE budget, and retry guard remain
fixed. Priors, equations, observations, and units are unchanged. No earlier
draws or adaptation state are resumed or pooled. Predictive selection increases
to 500 draws if the fit completes.

The driver `experiment.jl` records a 10,800-second sampling limit and a stop on
any retained divergence. Its watchdog `watchdog.ps1` allows 12,600 seconds for
the launched process, including loading and postprocessing, and retains the
native process handle to capture its exit code. The declared acceptance rule
remains R-hat below 1.01, bulk and tail ESS at least 400 for every parameter,
zero retained divergences, and zero depth hits across four complete chains.

On completion, the plan replays all 4,000 likelihood draws and 500 selected
latent predictions, validates all 8,000 runtime events against raw retained
statistics, profiles selected states, and saves full-run and matched 100-draw
block comparisons. A stopped run receives only partial-checkpoint and completed
runtime-event audits, without joint convergence or posterior-export claims.
Changes in both adaptation and retained length prevent attributing improvements
to either alone. The comparison code passes its eight-check preflight on the
audited pilot.

Sampling is in progress; no convergence result is claimed yet. The active run
and its sidecars are temporarily ignored by Git so commits cannot capture
partial outputs. Remove those two scoped ignore rules and publish the artifacts
when the run and its audits are finalized.
