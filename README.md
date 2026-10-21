# FermentModels

Mechanistic models of food fermentations in Julia: a shared modeling engine,
a reproduction of a published cocoa fermentation model, and practical sourdough
and vegetable fermentation calculators.

## Implementation status

The scaffold and data conventions (**P0.1–P0.3**) and shared simulation engine
(**P1.1–P1.7**) are implemented:

- TOML parameters with validated source, units, and confidence metadata.
- Temperature, pH, acid, and salt response functions.
- Acid speciation and an effective-buffer charge-balance solver.
- Monod/Contois ODE assembly, explicit stoichiometry, and optional log populations.
- Scheduled dilution feeds and step or linear temperature schedules.
- Bounded multistart Gaussian fitting with known or jointly estimated uncertainties,
  independent stationarity checks, convergence diagnostics, and local confidence intervals.
- Bayesian NUTS inference, sourced prior metadata, convergence diagnostics,
  and separate latent and observation predictive intervals.
- Local sensitivities, rank diagnostics, and profile likelihoods.
- A checksum-verified, inspected cocoa archive with its original three datasets
  and R/Stan source files (**P2.1**).
- A source-linked cocoa equation and scaling transcription, checked against the
  publisher's supplement (**P2.2**).
- The eight-state cocoa baseline, three sourced parameter sets, tidy source
  observations, and comparisons with digitized published curves (**P2.3**).
- A cocoa MLE/Bayesian workflow using the original likelihood and priors, with
  saved draws, predictive bands, and reproduction diagnostics (**P2.4**).
- Optional stiff integration, concurrent Bayesian chains within one process,
  per-start/draw checkpoints, and trace plots for longer cocoa fits.
- Saved-point MLE diagnostics for gradients, search bounds, local information,
  and sensitivity to integration tolerances.
- A verified 2022 cocoa benchmark archive, all 32 observation datasets with
  the published 23-trial subset marked, missing-lactate handling, a 32-variant
  mechanism inventory, and a sequential baseline fitting runner with
  provenance-checked restarts, trial subsets, and numerical diagnostics (**P2.5**).
  The Bayesian benchmark runner records the 2022 priors, checkpoints every
  retained draw, and exports predictive bands and complete pointwise likelihoods.
  Saved-draw checks compare posterior gradients, integration tolerances, and
  alternate solvers while preserving checkpoint provenance and failure labels.
  A reusable artifact audit recalculates sampling diagnostics, replays likelihoods
  and latent predictive bands, and checks source identities and physical units.
  New Bayesian benchmark runs count accepted ODE steps between output times,
  with a separate retry guard, and record initialization/warmup progress.
  Optional time and divergence stops preserve partial checkpoints.
  Extended runtime logs expose leapfrog work and energy diagnostics, with
  independent raw-chain exports and saved-state ODE/gradient cost profiles.
  A 48-case comparison at selected costly states supports retaining the current
  Tsit5 solver: both solvers pass numerical checks, but Rodas5P costs more there.
  A longer-adaptation pilot completes 400 retained draws without divergences,
  but remains unaccepted: maximum R-hat 1.072 and 387 depth-limit hits.
  A matched depth-9 pilot reduces depth hits to 133 and improves minimum
  bulk/tail ESS per leapfrog by 35%/26%. It passes the declared screen for
  considering a longer run, while still failing convergence criteria.
  The longer four-chain fit completes 4,000 retained draws without divergences:
  maximum R-hat 1.00684 and minimum bulk/tail ESS 727.28/390.54. It remains
  unaccepted because tail ESS is below 400 and 1,904 draws hit the depth cap.
  All 25,749 scientific replay checks pass, with 8,000 runtime events and
  4,000 raw retained statistics verified. Native process exit code is 0.
  See [benchmark provenance and comparison limits](docs/cocoa_benchmark.md).

The cocoa baseline uses extracted published estimates. Its figure reproduction
still has a documented discrepancy. Fresh Bayesian fits for all three datasets
pass parameter and convergence checks. Figure agreement remains incomplete:
Camu's mean and box 2's median each have one point outside the fixed 2% tolerance.
Sourdough and vegetable models remain planned. See
[the implementation plan](docs/PLAN.md) for task status and acceptance checks.

## Quickstart

Use Julia 1.10.8 or later in the 1.x series. From this directory:

```sh
julia --startup-file=no --project -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project -e 'using Pkg; Pkg.test()'
julia --startup-file=no --project examples/core_demo.jl
julia --startup-file=no --project examples/fit_monod.jl
julia --startup-file=no --project examples/bayes_monod.jl
julia --startup-file=no --project examples/identifiability.jl
julia --startup-file=no --project examples/cocoa_baseline.jl
```

The package name is `FermentModels`, regardless of the checkout directory name.
The [core demo](examples/core_demo.jl) simulates growth, acid production, a cool
interval, and two feeds, then writes a CSV to `examples/output/`. Read the
[example assumptions](examples/README.md) before adapting it.

The package uses OrdinaryDiffEq, DiffEqCallbacks, Optimization, ForwardDiff, and
Turing, with Julia's TOML standard library for parameter files.
Tests cover analytical growth limits, chemistry, conservation, positivity,
automatic differentiation, event timing, repeated feeding, and data provenance.
Tests include analytical likelihood checks, synthetic ODE recovery, held-out
predictive coverage, known non-identifiability, cocoa coordinate conversions,
and estimated-noise inference. One explicit expected failure records the unmet
2% published-curve tolerance; it is not treated as a successful reproduction.
The cocoa source verifier also checks the original archive and supplement.
The full suite passes 1,454 checks on Julia 1.10.11, using two threads and normal
compiled modules, with the one documented expected failure. Tests cover
stiff-solver agreement, seeded serial/threaded
sampling, checkpoint integrity, bounded restarts from saved MLE parameters,
independent stationarity and local-information diagnostics, and benchmark
provenance, scaling, state ordering, missing observations, and complete Bayesian
likelihood exports. The saved MLE benchmark artifacts pass 1,111 integrity
checks; the two Bayesian pilot fits pass a separate 475-check artifact audit.
The four-chain refinement passes 282 artifact checks, and the 72 saved-point
solver comparisons pass another 1,709 checks.
Both initial Bayesian pilot fits fail sampling diagnostics. A larger
four-chain `ghhp4` refinement passes the unchanged criteria with 8,000
retained draws. A tighter-tolerance four-chain `brwb3` attempt was interrupted
after two divergences, preserving 5,478 partial draws; it remains unaccepted.
All 377 focused telemetry checks pass on Julia 1.12.6: 104 callback/sampling
checks, 112 Bayesian benchmark checks, 94 posterior numerical checks, and
67 runtime/profile checks. Coverage includes task interleaving, exact agreement
with raw-chain statistics, unchanged seeded sampling, missing-lactate exports,
corruption rejection, and saved-state ODE/gradient cost profiles.
Disposable synthetic MLE fixtures carry fresh scores and implementation hashes;
historical fit artifacts retain their original provenance.
The reusable auditor replays all 12,000 completed pilot/refinement draws,
checking 1,242,000 likelihood values and reconstructing all 1,500 selected
latent predictions without changing sampling acceptance labels.

CI is configured to run the same instantiate/test command on pushes and pull
requests, with Julia 1.10 and 1.12 on Linux and Windows. A hosted run requires
pushing the repository to GitHub. Local manifests and Julia depot caches are
ignored; package compatibility is declared in `Project.toml`.
Initial compilation of the inference stack can take tens of minutes on a fresh
environment. CI caches dependencies and uses `JULIA_NUM_PRECOMPILE_TASKS=1` to
limit compiler memory use. When switching Julia
minor versions locally, resolve a separate manifest for each version instead of
reusing another version's standard-library dependency graph.

## Load parameters

```julia
using FermentModels

result = load_params("params/template.toml")
result.values        # NamedTuple of scalar values; unknowns are Julia missing
result.provenance    # Named rows: name, units, source, location, notes, confidence
```

The template contains an unresolved placeholder, not measured parameter values.
See [the parameter format](params/README.md) for required fields and how to record
missing values. The loader preserves source units; conversion to internal model
units belongs at the model-loading boundary.

## Ground rules

- Parameters come from published data or documented fits. Every value needs
  provenance; unavailable values remain missing.
- Raw data is preserved as acquired. Cleaning scripts write to `data/processed/`.
- Models use small, interpretable ODE systems and explicit observation models.
- Predictions must carry uncertainty bands because real ingredients and microbes
  differ from published systems.
- Vegetable safety estimates will be educational, not food-safety certification;
  they do not replace tested recipes or a pH meter.

See [agent instructions](AGENTS.md), [the tooling assessment](docs/TOOLING.md),
and [the source catalog](docs/SOURCES.md).

## Engine documentation

- [Secondary response equations and source checks](docs/secondary_models.md)
- [Acid speciation and buffer chemistry](docs/chemistry.md)
- [Model specification, kinetics, units, and conservation](docs/core_engine.md)
- [Feeding and temperature schedules](docs/events.md)
- [Point estimation and local uncertainty](docs/estimation.md)
- [Bayesian inference and predictive uncertainty](docs/bayesian.md)
- [Sensitivity and profile likelihood](docs/identifiability.md)
- [Original cocoa code, datasets, and reproduction issues](docs/cocoa_original.md)
- [Cocoa equation transcription and source discrepancies](docs/cocoa_equations.md)
- [Cocoa baseline, published parameters, and figure comparison](docs/cocoa_baseline.md)
- [Cocoa fitting workflow and reproduction diagnostics](docs/cocoa_fitting.md)
- [Response-function plots](docs/figures/secondary_models.png)

To fetch or verify the original cocoa archive:

```sh
julia --project scripts/fetch/cocoa.jl
julia --project scripts/fetch/cocoa_supplement.jl
julia --project scripts/verify_cocoa.jl
```

This reuses the existing archive only when its checksum and provenance agree;
it does not overwrite raw data. The extracted files and original workflow are
documented in P2.1, with equation transcription in P2.2. Generate the P2.3
parameter files and processed data, then run the baseline comparison:

```sh
julia --project scripts/prepare_cocoa.jl
julia --project examples/cocoa_baseline.jl
julia --project scripts/compare_cocoa_figures.jl
```

P2.4 commands and outstanding acceptance checks are in the fitting documentation.

Plotting dependencies live in a separate environment:

```sh
julia --project=scripts -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
julia --project=scripts scripts/secondary_figures.jl
julia --project=scripts scripts/cocoa_figures.jl
```

## Layout

```text
Project.toml              Package identity and compatibility bounds
src/FermentModels.jl      Package entry point
src/core/                Parameter loader and shared simulation engine
src/cocoa/               Cocoa module
src/sourdough/           Sourdough module
src/vegetables/          Vegetable module
params/                  Sourced cocoa parameters and placeholder template
data/raw/                Original downloads with provenance sidecars
data/processed/          Tidy observations, provenance, and fit reports
scripts/                 Reproducible fetching and response-function plots
examples/                Core simulation/inference and sourced cocoa examples
notebooks/               Pluto exploration and calculators (planned)
test/runtests.jl          Test entry point
docs/                    Plan, conventions, tooling, and literature sources
.github/workflows/ci.yml  GitHub Actions tests
```
