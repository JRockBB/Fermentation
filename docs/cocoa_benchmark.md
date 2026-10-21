# P2.5: 2022 cocoa benchmark

The [2022 paper](https://doi.org/10.1098/rsos.210274) screened 32 trials with
32 model variants. Its main analysis retained 23 trials with at least one
successful variant. All 32 observation datasets are obtainable from the
[archived v1.0 repository](https://doi.org/10.5281/zenodo.5510370), so this
implementation preserves both the retained set and the nine unsuccessful
trials. `datasets.csv` identifies membership explicitly.

## Acquisition and provenance

The archive, repository metadata, full-text XML, and
[supplement version 1](https://doi.org/10.6084/m9.figshare.19149492.v1)
are preserved in `data/raw/cocoa_benchmark/`, with acquisition sidecars.
The archive and PDF match the repository-supplied MD5 checksums, and their
SHA-256 values are pinned in the fetch and verification scripts.
[The member inventory](cocoa_benchmark_archive_files.csv) records all 78
extracted files. Repeated extraction verifies existing bytes.

The 2022 archive's Zenodo license is `other-open`, with no explicit license
file inside. This differs from the 2018 archive's CC0 declaration. We retain
it as attributed source material without asserting CC0 or relicensing its
code. The article and publisher supplement are CC BY 4.0. Our parser never
executes the downloaded R/Stan programs; the existing Julia baseline is used.

```sh
julia --startup-file=no --project scripts/fetch/cocoa_benchmark.jl
powershell -File scripts/fetch/extract_cocoa_benchmark.ps1
julia --startup-file=no --project scripts/prepare_cocoa_benchmark.jl
```

The extraction helper uses Windows/.NET; its output is platform-independent.
For other systems, extract the zip into `data/raw/cocoa_benchmark/` without
altering member bytes, then run the Julia preparation script. Tests and
preparation use the checked-in sources and require no network access.

## Observation representation

`posteriors/post_sims1.R` contains 28 complete eight-state datasets.
`posteriors/post_sims2.R` contains the four Camu-2008 datasets `ghhp2` through
`ghhp5`. The parser reads only numeric literals in their dataset branches,
checks dimensions and time ordering, and records the source file, branch,
and SHA-256 for each trial.

Times are hours. Observed metabolites and microbial biomass retain the
archive's mg/g-pulp units. The archive has already converted microbial
counts to biomass; we do not repeat or reverse that conversion. The canonical
state order is `(G,F,E,L,A,Xy,Xl,Xa)`. The four incomplete datasets use
`(G,F,E,A,Xy,Xl,Xa)` in R and append unobserved `L` last in Stan, so the loader
reorders them explicitly.

Observed series are divided by their own maxima, including the fixed
initial row. Initial values are excluded from the likelihood. For the four
incomplete trials, every lactic-acid observation remains `missing`. Their
physical initial lactic acid and scale also remain `missing`. The source
sets its **normalized hidden initial state to zero** (`post_sims2.R:139`);
the loader preserves that assumption without labelling it a measurement or
assigning a physical scale. All eight states are simulated, but only the
seven observed series enter the Gaussian likelihood.

Generated trial CSVs include physical values, roles, normalization scales,
and normalized values. Empty physical/scale fields are intentional.
`unobserved_initial` distinguishes the hidden initial condition from a
measured zero. The three overlapping datasets match the original 2018
archive's observations exactly.

```julia
include("scripts/prepare_cocoa_benchmark.jl")
using .CocoaBenchmarkData

trials = CocoaBenchmarkData.benchmark_trials()             # 23 retained trials
all_trials = CocoaBenchmarkData.benchmark_trials(scope=:all) # all 32
trial = only(filter(t -> t.trial == :ghhp2, trials))
# trial.data is ready for cocoa_fit_model / fit_mle / fit_bayes.
# Missing observations are omitted by the existing likelihood.
```

## Mechanism inventory

`mechanisms.csv` lists the five additions and their supplement equations.
`mechanism_variants.csv` enumerates the baseline and all 31 combinations,
including paths to both source Stan families and independently checked
parameter counts. There are 24 kinetic parameters in the baseline and
43 in the full combination, plus one shared Gaussian noise SD in each.
The variants are catalogued; this task does not implement their extra ODEs.

| Mechanism | Additional kinetics | Equation reference |
|---|---:|---|
| M1: first-order ethanol, lactate, acetate losses | 3 | Supplement 5–8 |
| M2: LAB fructose uptake and associated products | 7 | Supplement 9–14 |
| M3: acetate production by yeast | 2 | Supplement 15 |
| M4: lactate uptake by yeast and ethanol production | 4 | Supplement 16–19 |
| M5: acetate uptake and growth of AAB | 3 | Supplement 20–22 |

Two source inconsistencies are relevant to future extension work. The
supplement's M4 heading says “production” of lactic acid, but equations
16–19 and `mm4.stan` describe consumption. Supplement equation 13 omits the
baseline acetate-from-lactate term; `mm2.stan` retains it. These are recorded
here rather than silently resolved in a new extension implementation.

## Fit workflow and comparisons

```sh
julia --startup-file=no --project scripts/fit_cocoa_benchmark.jl --output=data/processed/cocoa_benchmark/fits/new_run --scope=all --starts=1 --maxiters=300 --outer_iterations=10 --maxtime=120
```

Runs are sequential, require a new output directory, and checkpoint each
completed start before moving on. The runner records source and implementation
hashes, numerical settings, candidate predictions, normalized training RMSE,
and Gaussian negative log likelihood. Its default stiff solver uses absolute
tolerance `1e-9`, relative tolerance `1e-8`, and a 10,000-step cap.

The numerical search bounds `[1e-5,5]` match the prior P2.4 MLE experiments.
Starting every parameter at `0.5` is a reproducible computational choice
based on the published kinetic-prior location, not an estimated biological
parameter vector. Additional starts sample the log bounds. MLE acceptance
requires native convergence, no optimizer time-limit flag, and the existing
independent projected log-gradient check at `1e-6`. The time budget is soft:
an inner optimization step, compilation, and diagnostic evaluation can overrun
it. An unaccepted finite point is labelled `unverified_candidate`; failed
simulation proposals are retained as `prediction_failure`.

To refine selected candidates, use a new output directory and identify the
previous run explicitly:

```sh
julia --startup-file=no --check-bounds=yes --project scripts/fit_cocoa_benchmark.jl --output=data/processed/cocoa_benchmark/fits/new_refinement --scope=all --trials=ghhp6,brwb3,ghhp2 --restart_from=data/processed/cocoa_benchmark/fits/new_run --maxiters=1000 --outer_iterations=20 --maxtime=180 --abstol=1e-11 --reltol=1e-10
```

Trial selection preserves the scope's original order and per-trial RNG seed,
so selecting a subset does not change its random starts. The first start uses
the saved candidate; any additional starts still sample the log search bounds.
Every selected restart is validated before output creation: raw source identity,
branch, initial conditions, observation times, state order, and hashes of the
data preparation, ODE, inference, and likelihood implementation must match.
The candidate report must agree with its per-start checkpoint. A mismatch
stops the batch; changed scientific code requires a fresh run or a separately
reviewed migration.

Each new input records hashes of its parent settings, input, report, and
checkpoint. An initial point within `1e-8` of either end of the log search
range is moved that fraction inside the box, with every affected parameter
recorded. No previous acceptance flag is carried forward. Historical files
retain their original implementation hashes even when the runner changes.
In particular, the archived `mle_pilot` predates the solver-budget implementation;
the current reader rejects it. Its historical diagnostic reports remain available,
but a new replay or restart requires the recorded implementation. The examples
here use `new_run` generated by the current code. Software tests create disposable
synthetic fixtures with freshly evaluated scores and separately verify rejection
of the historical run.

Inspect candidates before spending more optimization time:

```sh
julia --startup-file=no --check-bounds=yes --project scripts/diagnose_cocoa_benchmark.jl data/processed/cocoa_benchmark/fits/new_run data/processed/cocoa_benchmark/fits/new_diagnostics
# Optionally append a comma-separated trial subset, e.g. ghhp6,brwb3,ghhp2.
```

This produces per-trial reports and a comparison CSV at the saved ODE tolerances
and at 100 times tighter tolerances. It saves AD gradients, two centered
finite-difference checks, distances to numerical bounds, and singular values
of the local Gaussian information factor, including the noise parameter.
It also reports the conditional best noise SD at the fixed kinetic point.
The finite-difference errors divide each component's absolute difference by
`max(1,abs(AD gradient))` before taking the maximum. These are numerical
diagnostics, not an optimum acceptance rule or proof of structural
non-identifiability. Original source and implementation hashes accompany
every diagnostic batch; each completed tolerance check is saved.
Every selected trial must have a saved finite candidate. For a batch containing
`prediction_failure` trials, inspect their start checkpoints and pass the
remaining trial IDs explicitly; a failed attempt has no new point to diagnose.

The source's baseline diagnostics and rounded BMA weights are transcribed
from supplement Tables S3–S7 in `published_baseline.csv`. Eighteen of the
23 retained trials have a reported successful baseline fit. Missing values
represent printed dashes. [Transcription provenance and corrections](../data/raw/transcribed/cocoa2022/README.md)
are preserved separately from processed data.

The comparison CSV keeps our MLE training scores separate from the paper's
Bayesian diagnostics and BMA weights. BMA weights depend on the competing
variants and cannot be reproduced from a baseline MLE. The paper uses
PSIS-LOO for predictive fit quality; its plotting script reads posterior
`.Rsave` files that are absent from the archive. `GENERAL_INFO.txt` says
original outputs are available on request. Consequently, exact published
PSIS-LOO values remain missing and direct metric comparison is unavailable.
No contact with the authors has been initiated.

## Bayesian baseline workflow

The benchmark Bayesian runner uses the existing hand-written baseline and
Gaussian likelihood, with the source's positive-truncated `Normal(0.5,0.3)`
kinetic priors and half-`Cauchy(0,1)` noise prior. These distributions match
the 2018 implementation; their new provenance points to the independently
inspected 2022 `Stan_models/mm0.stan` and `Stan_models/camu2008/mmc0.stan`
model blocks. Numerical MLE bounds are not applied to Bayesian draws.

```sh
julia --startup-file=no --project scripts/fit_cocoa_benchmark_bayes.jl --output=data/processed/cocoa_benchmark/fits/new_bayes --scope=all --trials=brwb3,ghhp4 --samples=2000 --warmup=1000 --chains=4
```

Trials run sequentially. `--chain_execution=threads` runs a trial's chains
concurrently within one Julia process; start Julia with `--threads=2` or
more. Selecting a subset retains each trial's scope-specific seed.
Every output directory must be new. Retained draws are checkpointed per
chain, including divergence and tree-depth indicators. Warmup is excluded.
Draw checkpoints lack sampler adaptation state and cannot resume NUTS.

The archived driver comments specify four chains, 3,000 total iterations
with one-third warmup, target acceptance in `0.995–0.999`, tree depth 10,
initial step size 0.1, and seed 141085. The source initializes parameters
with absolute `Normal(0.5,0.2)` draws. The new runner defaults to 2,000
retained draws, 1,000 warmup, four chains, target acceptance 0.995, depth 10,
and base seed 141085; it uses Turing prior initialization and automatic
step-size selection. These differences are explicit in saved settings.
The archived driver accepts configurable arguments; these recorded defaults
are not a claim that every published fit used identical tuning.

Both source baseline families use Stan RK45 with absolute and relative
tolerances `1e-6`, but their step caps differ: `mm0.stan` uses 10,000 and
`mmc0.stan` uses 1,000,000. The runner defaults to Julia's Tsit5 at the same
tolerances and selects those caps per trial with `--solver_maxiters=source`.
A positive integer overrides both caps. Rodas5P is available via
`--solver=rodas5p`. The integrators are different implementations; these
settings do not establish an exact RStan reproduction.

New Bayesian benchmark runs default to `--solver_maxiters_scope=output_interval`.
This limits **accepted steps between requested output times**, consistent with
the [Stan ODE controls](https://mc-stan.org/docs/2_24/stan-users-guide/control-parameters-for-ode-solving.html)
and its [dense-output RK45 call](https://mc-stan.org/math/ode__rk45_8hpp_source.html).
The implementation keeps one OrdinaryDiffEq integrator and its adaptive history.
It interpolates observations without forcing time stops. A step crossing several
output times is charged to the first pending output; subsequent interpolated
outputs require no additional steps.

`--solver_max_step_attempts=500` separately bounds attempted steps, including
rejections, within each accepted step. This follows the distinct retry guard in
[Boost 1.72's dense-output stepper](https://github.com/boostorg/odeint/blob/boost-1.72.0/include/boost/numeric/odeint/stepper/dense_output_runge_kutta.hpp)
and its [step/observer adapters](https://github.com/boostorg/odeint/blob/boost-1.72.0/include/boost/numeric/odeint/integrate/check_adapter.hpp).
Julia's cumulative iteration counter is preserved; only its attempt ceiling is
advanced between successful steps. Failure of either allowance rejects the
entire prediction. Solver algorithms, sensitivity/error control, and initial
step selection still differ from Stan.

`--solver_maxiters_scope=solve` retains the previous whole-trajectory limit on
attempted steps. The package-level `cocoa_fit_model` still defaults to this mode
for compatibility, as do the existing MLE and 2018 workflows. Previous saved
runs with no explicit scope field mean `solve`. Their artifacts and acceptance
labels remain unchanged. The completed-fit reader continues to require matching
implementation hashes: audits of older code remain historical records, and
replaying those fits requires their recorded implementation. No hash guard is
disabled to treat them as new interval-budget runs.

Every new benchmark run also writes `runtime.csv` per trial. Each row records a
completed sampler call, chain, phase, retained-draw index, monotonic elapsed
time, time since that chain's previous completion, and divergence/tree depth
when available. Initialization has no transition diagnostics. In the pinned
Turing interface it occupies the first discarded slot, followed by `warmup-1`
discarded transitions; the existing discard/adaptation convention is preserved.
Monitoring temporarily retains warmup internally and removes it before all
posterior summaries, checkpoints, likelihoods, and predictive exports.

Runtime schema 2 also records leapfrog counts (`n_steps`), leapfrog step size,
mean acceptance probability, Hamiltonian energy, energy error, and maximum
energy error from each completed transition. These are
[AdvancedHMC sampler statistics](https://github.com/TuringLang/AdvancedHMC.jl/blob/main/src/trajectory.jl),
not ODE step counts. The callback adds no density evaluations. Initialization
fields remain empty; divergent transitions can retain nonfinite energy errors.
Completed runs separately export `sampler_statistics.csv` from the raw retained
chain, allowing the callback records to be checked independently. Historical
logs without these fields cannot reconstruct them.

`--max_sampling_seconds=600` gives each trial a soft sampling deadline, checked
after completed calls; compilation and scheduling count toward it. It cannot
interrupt a stalled call. `--stop_on_divergence=true` stops after the first
retained divergence, after saving that draw. Both options preserve partial
checkpoints and write a `sampling_stopped` report with completed-call and
retained-draw counts. Warmup divergences do not trigger this stopping rule.
The default is unlimited time and no early divergence stop. Neither setting
changes the fit acceptance criteria. A partial fit has no joint acceptance,
likelihood export, or predictive comparison.

Callback diagnostics are local to each sampler call, and both hooks execute
under one lock. A forced task-interleaving test guards against a previous shared
variable that could mix diagnostic fields between chains while checkpoint I/O
yielded. The initial diagnostic pilot's discrepant runtime logs are preserved;
its corrected repeat has an independent per-row runtime/checkpoint audit.
Historical threaded per-draw diagnostic fields have this documented limitation;
parameter draws and diagnostics computed directly from the raw sampler chain
are independent of that callback variable. See the
[fit history](../data/processed/cocoa_benchmark/fits/README.md#paired-runtime-diagnostic).

For a short runtime diagnostic with an explicit budget (not a converged fit):

```sh
julia --startup-file=no --check-bounds=yes --threads=2 --project scripts/fit_cocoa_benchmark_bayes.jl --output=data/processed/cocoa_benchmark/fits/new_runtime_check --scope=all --trials=brwb3 --samples=25 --warmup=50 --chains=4 --max_depth=8 --abstol=1e-8 --reltol=1e-8 --chain_execution=threads --solver_maxiters_scope=output_interval --max_sampling_seconds=600
```

Each completed trial saves:

- `chain_###.csv`: all retained joint draws, with chain/draw identifiers and
  per-draw diagnostics.
- `sampler_statistics.csv`: retained transition diagnostics exported directly
  from the raw chain, independently of the runtime callback.
- `parameters.csv`: marginal summaries, rank R-hat, bulk ESS, and tail ESS.
- `predictive.csv`: separate latent and observation bands, based on up to
  `--predictive_draws=500` draws; `predictive_selection.csv` records which
  joint draws were selected. Physical columns use mg/g pulp. Hidden lactate
  has only normalized latent bands: its measurement bands, physical values,
  and physical scale remain empty.
- `log_likelihood.csv`: a row for every retained joint draw and a column
  for each observed state-time cell. `likelihood_observations.csv` maps
  these columns to canonical state order, time, and normalized observation.
  Missing states and the fixed initial row are excluded. Failed evaluations
  stop export rather than silently discarding draws.
- Input, prior, settings, and report TOMLs with source and implementation
  hashes. Reports distinguish sampling, postprocessing, completed diagnostic
  results, and errors; failed runs retain their checkpoints.

`posterior_diagnostics_passed` requires at least two chains, rank R-hat
below 1.01, bulk and tail ESS at least `100 × chains`, and no divergences
or tree-depth-limit hits. Unconverged draws and their summaries remain
explicitly labelled `posterior_diagnostics_failed`; they are not accepted
biological estimates. Passing these checks is necessary evidence, not proof
of convergence or of a valid biological model.

Pointwise likelihoods use the source's fixed full-trial normalization and
initial states. Their export prepares the inputs for a future PSIS-LOO
calculation; it is not itself a PSIS-LOO score. Both computed and published
PSIS-LOO columns therefore remain empty in the comparison table. Baseline
fits alone also cannot reproduce weights across all 32 competing variants.

Plot a completed trial from its saved files using the separate plotting
environment:

```sh
julia --startup-file=no --project=scripts scripts/cocoa_benchmark_figures.jl data/processed/cocoa_benchmark/fits/new_bayes/brwb3
```

This writes PNG/SVG predictive panels and all-parameter trace plots with
diagnostic status visible. Measured series use mg/g pulp; hidden lactate
uses a clearly labelled normalized latent axis. A figure provenance TOML
records the plotting script and input file hashes.

The main text describes bulk/tail ESS above 100, while supplement Section 4
states 400. We retain the printed success labels, including baseline entries
with ESS below 400, without changing our own convergence requirements.

## Auditing saved Bayesian fits

```sh
julia --startup-file=no --check-bounds=yes --project scripts/audit_cocoa_benchmark_bayes.jl data/processed/cocoa_benchmark/fits/new_bayes data/processed/cocoa_benchmark/fits/new_bayes.audit.toml
```

An optional third argument selects comma-separated trial names. The audit
requires completed reports and matching source, prior, and implementation
hashes. An `aborted_run.toml` record explicitly disqualifies an interrupted
experiment, even if completed-looking files are present. It recalculates
marginal summaries, rank R-hat, bulk/tail ESS, and
the unchanged acceptance rule from the ordered chain checkpoints. It checks
the comparison table against those results and the published transcription.

Every pointwise likelihood row is checked for finiteness, width, and exact
chain/draw identity. The complete observation map must match the source
state/time cells, excluding initial conditions and missing observations.
By default, the first, middle, and last retained draw of each chain are
replayed through the saved ODE solver to verify their likelihoods. The Julia
function's `replay_per_chain` keyword controls the evenly spaced replay
grid; a value of one selects the middle draw.
Set it to at least the saved samples per chain to replay every likelihood
row. Replay evaluates Gaussian log densities independently of the export
function.

All recorded predictive selections are replayed to reconstruct the latent
bands. Normalized and physical observations, scales, and unit conversions
are checked. Hidden lactate must have only normalized latent bands. The
Gaussian observation-noise draws were not saved, so measurement bands can
only be checked for finite ordered endpoints and consistent physical
conversion, not reproduced exactly.

The audit writes a new TOML only after every check succeeds and confirms
that the checked parent files have not changed. It records input and audit
implementation hashes, replay counts, and separate integrity and sampling
statuses. It never edits the fit. A correctly saved fit that fails MCMC
diagnostics can pass the integrity audit; this does not accept its parameter
estimates or establish model validity. Likelihood replay checks a subset by
default; the recorded replay count identifies whether every draw was checked.
The audit does not compute PSIS-LOO.

## Numerical checks of saved posterior draws

```sh
julia --startup-file=no --check-bounds=yes --project scripts/diagnose_cocoa_bayes.jl data/processed/cocoa_benchmark/fits/new_bayes data/processed/cocoa_benchmark/fits/new_bayes_numerics
```

An optional third argument selects comma-separated trial names. The reader
requires a completed sampling report, complete ordered chain checkpoints,
matching source inputs and prior provenance, and unchanged scientific
implementation hashes. It checks retained counts and divergence/depth
flags against the report before creating an output directory. Use a new
directory for each diagnostic run.

The default subset contains three evenly spaced retained draws per chain,
each chain's minimum and maximum `yc11`, every divergent retained draw,
and its preceding retained state when present. Duplicate selections are
merged, with all selection reasons recorded. A flagged retained draw does
not necessarily identify the leapfrog state that caused the divergence.

Each selected point is evaluated with the saved integrator/tolerances,
with 100 times tighter absolute and relative tolerances, and with the
alternate integrator (Tsit5 or Rodas5P) at those tighter tolerances. The
saved per-trial ODE step cap is preserved. These checks evaluate the
negative log density in explicit log-parameter coordinates:
Gaussian NLL minus the source log priors minus the log-Jacobian `sum(z)`,
where `parameters = exp(z)`. Numerical MLE search bounds are not applied.

ForwardDiff gradients are compared with centered finite differences at
steps `1e-4` and `5e-5`. The scaled error is the largest component of
`abs(FD-AD)/max(1,abs(AD))`. Cross-case reports retain changes in NLL,
normalized latent predictions, and gradients; gradient differences use
the saved case's gradient in the denominator. A failed solve remains an
explicit failed case. Programming errors propagate.

`comparison.csv` is checkpointed after every case; per-case TOMLs retain
the parameter point, full gradient, finite differences, latent predictions,
and numerical controls. Parent settings, inputs, priors, reports, chains,
and diagnostic implementations are hashed. These selected-point checks
do not certify numerical accuracy over the whole posterior or establish
MCMC convergence. They never change the fit's acceptance status.

## Validation and remaining work

To inspect a completed run that records runtime schema 2:

```sh
julia --startup-file=no --check-bounds=yes --project scripts/diagnose_cocoa_runtime.jl data/processed/cocoa_benchmark/fits/new_runtime_check data/processed/cocoa_benchmark/fits/new_runtime_diagnostics brwb3
```

This requires matching source, priors, inputs, and implementation hashes. It
checks all runtime rows for call order, phase boundaries, timestamp consistency,
retained checkpoint correspondence, and exact agreement with the independently
exported raw-chain statistics. Its summaries report time, leapfrog work, step
sizes, and acceptance by chain and phase. Elapsed seconds per leapfrog include
scheduling, callback work, compilation, and GC; they are not ODE timings.

For each chain it selects the slowest retained call, its preceding retained
state when available, and a median-duration control. Selection is deterministic,
including tie handling. Separate warmed replays record accepted/rejected ODE
steps, RHS evaluations, repeated primal solve and log-coordinate gradient
timings, allocation bytes, and GC seconds. Per-point TOMLs retain parameters,
predictions, gradients, selection reasons, and failures. Nothing is refitted.
The retained state is not the complete Hamiltonian trajectory; these replays
cannot identify a failed leapfrog state or establish why a transition was slow.
Warmup parameter states are not retained. Parent hashes are checked again at
completion, and no diagnostic changes a fit's sampling acceptance.

The first [telemetry repeat and cost profiles](../data/processed/cocoa_benchmark/fits/README.md#extended-runtime-and-saved-state-cost-profiles)
preserve all nine scientific exports from the preceding bounded pilot byte for
byte. Its scientific audit passes 1,449 checks; runtime validation checks 300
events and all 100 retained draws against the independent raw-chain export.
All 12 selected profiles complete. Two chain-4 states have over 1,100 rejected
ODE steps and 27–38 ms median gradient costs, compared with 2–6 ms for the other
selected states. These provide specific numerical cases for solver comparisons;
they do not establish the cause of earlier divergences or delays.
The pilot still fails the unchanged convergence/depth criteria.

The subsequent [48-case solver comparison](../data/processed/cocoa_benchmark/fits/README.md#solver-comparison-at-the-costly-states)
evaluates both existing solvers at `1e-8` and `1e-10`, with numerical and
performance screens recorded in advance. All cases pass the numerical screen,
including finite-difference checks at two step sizes. Five rotating timing
rounds per case yield 240 gradient measurements. Rodas5P `1e-8` is 5.10–5.61
times slower at the two expensive states and has a median paired cost ratio of
8.36 across all 12 states. Neither alternate configuration passes the performance
screen. The current Tsit5 `1e-8` settings are retained, with no alternate-solver
pilot launched. The comparison passes 504 Julia integrity checks; an independent
audit recomputes its metrics and decisions and verifies 100 hashes. These local
measurements do not establish global numerical accuracy or sampling convergence.

The [longer-adaptation pilot](../data/processed/cocoa_benchmark/fits/README.md#bounded-longer-adaptation-pilot)
then retains 100 draws per chain after 250 discarded calls, using the same four
seeds, depth 8, and Tsit5 settings. All 400 draws complete with zero retained
divergences. Maximum R-hat is 1.07215 and minimum bulk/tail ESS is 47.49/39.09;
387 draws hit the depth cap. The fit fails the unchanged diagnostic rule.
The full-run R-hat improvement cannot be attributed to adaptation alone because
retained length also increased. Matched 25-draw blocks and ESS per draw are
saved as descriptive comparisons. Persistent depth hits and poor mixing prevent
advancement under the experiment's declared follow-up rule.

The scientific audit passes 3,249 checks, replaying all 38,400 likelihood values
and 50 selected latent predictions. Runtime validation checks 1,400 events and
all 400 retained statistics against raw-chain exports. Twelve selected cost
profiles complete. An independent audit verifies 165 hashes and records a
launcher limitation: its native child exit code is unavailable. The receipt and
initial failed checker are preserved; a revised checker verifies complete
scientific outputs and the driver's 19 passing assertions without inferring a
zero exit code. A separate probe confirms a handle-retention fix for future
launchers.

The matched [depth-9 comparison](../data/processed/cocoa_benchmark/fits/README.md#bounded-tree-depth-comparison)
keeps 250 discarded calls and 100 retained draws per chain. All 400 draws
complete without divergences. Maximum R-hat improves to 1.05479, minimum
bulk/tail ESS to 96.72/74.07, and depth hits fall from 387 to 133. Retained
leapfrog work increases from 100,996 to 152,048; minimum bulk/tail ESS per
leapfrog improves by 35%/26%. These results pass the prospectively recorded
screen for considering a longer run, while the fit remains unaccepted.
Elapsed sampling time rises from 1,139.5 to 2,144.1 seconds; measured bulk/tail
ESS per elapsed second improves by only 8.2%/0.7%. Wall times include warmup,
compilation, scheduling, and GC, and a single seed set does not establish a
reproducible speed improvement.

All 3,249 scientific replay checks and 19 driver integrity assertions pass.
The runtime audit checks 1,400 events and all 400 retained statistics;
12 selected profiles complete. An independent audit verifies 168 hashes,
recomputes the comparison decision, and confirms native exit code 0 from the
corrected launcher. No production source, model, prior, parameter file, or
dependency changes accompany this experiment.

The [longer depth-9 refinement](../data/processed/cocoa_benchmark/fits/README.md#longer-depth-9-convergence-attempt)
subsequently completes 4,000 retained draws across four fresh chains, each
with 1,000 discarded calls. Sampling finishes in 8,885.8 seconds within its
10,800-second limit, with no retained divergence. Maximum R-hat 1.00684 and
minimum bulk ESS 727.28 pass, but minimum tail ESS 390.54 for `ks4` is below
400. There are 1,904 depth-limit hits (47.6%), concentrated in chains 2 and 4.
The fit remains unaccepted under the unchanged rule.

The longer run's minimum bulk/tail ESS per retained leapfrog falls by 30%/51%
relative to the pilot. Its first and last 100-draw blocks have maximum R-hat
1.12843 and 1.20979, respectively. These short blocks are descriptive, and
changes in both warmup and retained length prevent attributing differences to
either alone. The improved full-run R-hat does not resolve the failed gate.
The next step is diagnosis of `ks4` tail mixing and persistent depth saturation
before another bounded sampler experiment, rather than extending this run.

All 25,749 scientific checks pass, including replay of 384,000 pointwise
likelihood values and 500 selected latent predictions. The runtime audit
checks 8,000 events and exact correspondence with all 4,000 raw retained
statistics; 12 state-cost profiles complete. The driver passes 23 integrity
assertions, and an independent audit verifies 169 hashes and native exit code
0. Source data, equations, priors, numerical controls other than the declared
run lengths, and implementation hashes remain unchanged. Observation-noise
draws were not saved, so measurement bands receive structural and unit checks
only. No PSIS-LOO score or accepted predictive comparison is claimed.

Tests check immutable archive contents, all trial dimensions and scales,
the three overlapping 2018 datasets, seven-series state reordering, strict
literal parsing, missing-lactate likelihood exclusion, all 32 variant
parameter counts, and reproducible processed outputs. MLE training metrics
are checked against the Gaussian likelihood directly. Restart checks cover
source and model mismatches, report/checkpoint agreement, boundary handling,
stable subset seeds, and a real hidden-lactate restart. Diagnostic report
checks validate likelihoods, observation counts, parameter order, gradients,
and parent-file hashes.
Bayesian checks independently match both source prior blocks, replay
Gaussian likelihoods, preserve every retained draw and its chain identity,
and verify that hidden lactate has no exported measurement or physical bands.
An intentionally short one-chain integration run must fail the diagnostic
gate while still producing complete checkpoints and likelihood exports.
Another 94 checks cover analytic posterior gradients including the
log-Jacobian, missing observations, unbounded prior support, deterministic
draw selection, corrupt checkpoints, provenance mismatches, failure
handling, and real Tsit5/Rodas5P report generation.
These 94 checks pass on Julia 1.10.11 and 1.12.6.
Another 136 checks validate the reusable artifact auditor on both Julia
versions. Synthetic checkpoints with real source metadata exercise full
likelihood replay, latent-band reconstruction, missing-lactate fields,
published-reference failures, unit labels, interrupted-run rejection, and
corruption of diagnostics, draw identities, likelihoods, and physical units.
The fixtures are software tests, not fitted biological estimates.

The telemetry changes pass 377 focused checks on Julia 1.10.11 and 1.12.6:
104 callback/sampling checks, 112 Bayesian benchmark checks, 94 posterior
numerical checks, and 67 runtime/profile checks. Runtime tests reject corrupted
ordering, timestamps, checkpoint fields, and each exported raw-chain statistic;
they also cover initialization, zero warmup, divergent nonfinite energy errors,
deterministic selection, explicit integration failures, and parent provenance.
Solver-budget checks in the full suite cover
accepted versus rejected steps, exact allowance boundaries, interpolated
outputs, failure propagation, and unchanged trajectories and gradients.
Callback tests force overlapping tasks and match every recorded diagnostic
to its raw sampler transition. Monitoring preserves seeded parameter and
predictive draws within each execution mode. The first
[bounded fit survey](../data/processed/cocoa_benchmark/fits/README.md)
completed all 32 trials with finite candidates; none passed optimizer
convergence or stationarity. The saved comparison table covers all trials,
and an expanded 1,111-check audit validates source identities, parent files,
saved observations, parameter values, scores, diagnostic reports, and
acceptance/failure labels across the pilot and refinements. The 64-case pilot
scan covers all 32 datasets. Two of the three longer restarts improve their
training scores but fail stationarity and reach the time limit; the third
hits an ODE step limit. No refined point is accepted as an MLE.
The first Bayesian pilot completes `brwb3` and hidden-lactate `ghhp4`, saving
2,000 retained draws per trial and complete pointwise likelihoods. `brwb3`
fails R-hat and divergence checks; `ghhp4` fails tail ESS. Both remain
explicitly unaccepted. A separate 475-check artifact audit reproduces the
saved marginal diagnostics, selected-draw latent bands, and likelihoods,
and checks source/implementation hashes, units, and complete draw mapping.
A fresh four-chain `ghhp4` refinement retains 8,000 draws after 1,000 warmup
per chain. It passes the unchanged diagnostic rule: maximum R-hat 1.00226,
minimum bulk/tail ESS 3219.46/1920.08, and zero divergences or depth hits.
Its source priors, equations, data, solver controls, and target acceptance
are unchanged; its experiment record precedes sampling. It exports all
840,000 pointwise likelihood values. The earlier failed pilot is preserved.
The refinement passes an independent 282-check artifact audit. A separate
72-case numerical scan of 24 points from the original pilots completes
with both integrators. Tighter tolerances reduce the worst gradient
discrepancy for both trials; NLL changes stay below `2.88e-5` on this subset.
Its reports pass 1,709 artifact checks. Full results and limitations are
recorded in the [fit history](../data/processed/cocoa_benchmark/fits/README.md).
The reusable auditor subsequently replays all 12,000 completed pilot and
`ghhp4` refinement draws, verifying 1,242,000 likelihood values and all
1,500 selected latent predictions. Its saved audit reports preserve the
two failed pilot labels and the refinement's passing sampling status.
A tighter-tolerance four-chain `brwb3` attempt was interrupted after saving
5,478 partial draws with two retained divergences. The unchanged
zero-divergence rule already precluded acceptance. Its interruption record
preserves all partial chains; no completed joint diagnostics, likelihood
export, or predictive bands are claimed. A separate ten-step initialization
probe completes without an exception, leaving the runtime delay unexplained.
The required `Pkg.test()` run on Julia 1.10.11 passes 1,454 checks with the
one existing expected cocoa figure failure (18 minutes 44.3 seconds in the
test summary). Both Julia validations use two threads and normal compiled
modules. Historical artifacts were not rewritten.
That full-suite result predates the depth comparison and longer refinement;
both new experiments received their own replay and integrity checks without
repeating unchanged package tests.

P2.5 remains open until baseline fits are adequately converged and a
like-for-like predictive comparison is available. A bounded initial survey
is useful for finding simulation and optimization failures; it does not
establish parameter estimates or reproduce the paper's Bayesian fits.
