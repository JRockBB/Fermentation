# Plan

Phases are ordered by dependency. Each task is sized for one coding-agent session and
has an ID, inputs, a deliverable and a "done when" check. Tasks marked 🔎 need literature
lookup or data digitizing as well as coding; flag anything you could not find rather
than filling gaps.

---

## Shared modeling framework

All three foods use the same structure.

**State** `u = [X₁…Xₙ, S₁…Sₘ, P₁…Pₖ]`: microbial populations, substrates (sugars) and
products (lactic acid, acetic acid, ethanol, CO₂).

**Growth rate of organism i**

```
μᵢ = μmaxᵢ · γ_T(T) · γ_pH(pH) · γ_HA(undissociated acids) · γ_NaCl(salt) · f(S)
```

- **Primary dynamics.** `f(S)` is Monod or Contois uptake. Use an optional Baranyi lag
  term for slow starts.
- **Secondary models (the γ terms).** Each is a dimensionless factor between 0 and 1.
  - Temperature: cardinal temperature model with inflection (Rosso) and the Gänzle 1998 form.
  - pH: cardinal pH.
  - Undissociated acid: linear or power-law inhibition up to a minimum inhibitory concentration (MIC).
  - Salt: cardinal or linear inhibition.
- **Chemistry.**
  - pH is computed from acid concentrations using a buffer model.
  - The undissociated acid fraction is `1 / (1 + 10^(pH − pKa))`, with pKa about 3.86 for lactic and 4.76 for acetic.
- **Events.**
  - Feeding: a dilution event (sourdough).
  - Temperature schedules: piecewise or interpolated T(t), e.g. fridge retards.
- **Observation model.** States map to what was measured (log CFU, g/L, pH, rise
  height), with a noise model.

---

## Phase 0 — Scaffold

### P0.1 Repository scaffold
- **Status:** Scaffold implemented. Local `Pkg.test()` passes on Julia 1.10.11 and
  1.12.6 (the 1.12 run used `--compiled-modules=no`). GitHub Actions is configured
  for pushes and pull requests. The first hosted CI run is pending a push to GitHub.
- **Output:**
  - Julia package `FermentModels`, using the `src/` layout from the README
  - `[compat]` pins
  - empty module files
  - `test/runtests.jl`
- **Done when:** `Pkg.test()` passes on an empty test set, and CI (GitHub Actions)
  runs tests on push.

### P0.2 Parameter and provenance loader
- **Status:** Implemented as exported `load_params(path)`, with a placeholder
  template and format documentation in `params/`. All 91 tests pass via `Pkg.test()`
  on Julia 1.10.11 and via `test/runtests.jl` on Julia 1.12.6 with
  `--compiled-modules=no` at that scaffold stage. Sourced numerical cocoa
  parameter files were subsequently added in P2.3; current full-suite validation
  uses normal compiled modules.
- **Output:**
  - `load_params(path)` that reads the TOML format in `AGENTS.md` and returns a
    NamedTuple of values plus a provenance table
  - Rejects entries with no `source` unless `confidence = "placeholder"`
- **Done when:** tests cover a valid file, a missing source, and unit-string presence.

### P0.3 Data conventions
- **Status:** Implemented. The cocoa archive was fetched with URL, UTC timestamp,
  SHA-256 and license metadata; its MD5 matches the Zenodo record. Download tests
  cover reuse, corruption, incomplete pairs, and checksum failures (109 tests total).
- **Output:**
  - `data/raw/README.md` and `data/processed/README.md` describing the layout
  - A `scripts/fetch/` stub that downloads a file and records URL, date and checksum
- **Done when:** one real file is fetched end-to-end (use the cocoa Zenodo zip, P2.1).

---

## Phase 1 — Core engine

### P1.1 Secondary models
- **Status:** Implemented with ForwardDiff checks and source-verified equations.
  PNG and SVG response plots are saved in `docs/figures/` and visually reviewed.
  The Gänzle exponent sign is corrected against the original paper; its lower
  temperature tail has no finite cutoff. See `docs/secondary_models.md`.
- **Output:** `src/core/secondary.jl` with:
  - `γ_T_ctmi(T; Tmin, Topt, Tmax)` (Rosso cardinal temperature model)
  - `γ_T_ganzle(T; a, b, c, Tref)` (Gänzle 1998 form, `a·x^b·e^(-c·x)` with positive `c` and `x = Tref − T`, normalized to 1 at the optimum; sign corrected against Table 1)
  - `γ_pH_cardinal`
  - `γ_acid(HA; MIC, α)`
  - `γ_salt(NaCl; max, α)`
- **Done when:**
  - Each function is in [0, 1], 0 outside its limits and 1 at its optimum
  - Plots are saved to `docs/figures/` for visual review
  - Functions are AD-compatible (ForwardDiff test)

### P1.2 Acid and pH chemistry
- **Status:** Implemented with explicit buffer parameters, a bracketed charge-balance
  solve, and implicit derivative checks. Titration and monotonicity tests pass;
  food-specific buffer calibration remains P4.2.
- **Output:** `src/core/chemistry.jl` with:
  - undissociated-fraction helpers
  - a buffer model `pH(acids, buffer_params)` that solves the charge balance for a
    substrate with given buffer capacity
  - a simple Henderson–Hasselbalch fallback
- **Done when:**
  - Titration-curve tests match hand calculations
  - Monotonic: more acid gives lower pH

### P1.3 Model assembly
- **Status:** Implemented with Monod/Contois uptake, explicit stoichiometry,
  organism-specific modifiers, and optional natural-log populations. Tests match
  the analytic Monod solution and check positivity, conservation, and derivatives
  through the ODE solver. See `docs/core_engine.md`.
- **Output:**
  - A small spec type (organisms, substrates, products, yield/stoichiometry matrix, which γ terms apply to whom)
  - `build_rhs(spec)` returning `f!(du, u, p, t)`
  - Optional log-transformed populations, for numerical stability near zero
- **Done when:**
  - A 1-organism, 1-substrate, 1-product system reproduces the analytic logistic/Monod limit
  - Positivity and mass balance tests pass

### P1.4 Events and environment
- **Status:** Implemented. Exact feed timing, endpoint feeds, callback reuse,
  repeated 1:1:1 cycle convergence, and analytic temperature-step growth checks
  pass. All 243 tests pass via `Pkg.test()` on Julia 1.10.11 and 1.12.6 (the latter
  with `--compiled-modules=no`). The combined `examples/core_demo.jl` also runs
  successfully and writes a 36-hour CSV; no calibrated food predictions are implied.
- **Output:**
  - Dilution callback (multiply all states by a factor, then add fresh substrate)
  - Time-varying temperature input `T(t)` from a schedule
- **Done when:**
  - Repeated 1:N:N feedings reach a periodic steady state in a test
  - A step change in T changes the growth rate as expected

### P1.5 Point estimation
- **Status:** Implemented with named positive parameters, finite bounds, Gaussian
  measurement uncertainties, reproducible multistart diagnostics, and local 95%
  intervals. Weighted-regression and noisy Monod recovery tests pass (291 tests
  total on Julia 1.10.11). See `docs/estimation.md` and `examples/fit_monod.jl`.
- **Output:** `fit_mle(model, data; starts=50)`, multistart with Optimization.jl
  (LBFGS + ForwardDiff). Parameters fitted on log scale with bounds. Returns the best
  fit, all starts, and residuals.
- **Done when:** synthetic recovery test passes (AGENTS.md → Testing).

### P1.6 Bayesian estimation
- **Status:** Implemented with explicit prior provenance, NUTS/ForwardDiff,
  rank R-hat and bulk/tail ESS, and separate latent/observation predictive bands.
  Exact-posterior, held-out coverage, and Monod ODE recovery checks pass (319 tests
  total on Julia 1.10.11). See `docs/bayesian.md` and `examples/bayes_monod.jl`.
- **Output:** `fit_bayes(model, data, priors)` using Turing.jl NUTS. Priors default
  to log-normal centred on literature values with source-based widths. Includes
  posterior predictive bands and a convergence report (R̂, ESS).
- **Done when:** synthetic recovery test passes, and posterior predictive coverage is
  about 95 % on held-out synthetic data.

### P1.7 Sensitivity and identifiability
- **Status:** Implemented with scaled ForwardDiff Jacobians, SVD diagnostics,
  multistart nuisance fits, and likelihood crossing brackets that preserve failed
  grid points. All 373 tests pass on Julia 1.10.11 and 1.12.6 with normal compiled
  modules, including 54 sensitivity and profile checks and a known
  non-identifiable toy model. The runnable example also passes.
  See `docs/identifiability.md`.
- **Output:** local sensitivities (ForwardDiff) and a profile-likelihood helper. These
  identify which parameters the data can and cannot constrain.
- **Done when:** a known non-identifiable toy model is flagged.

---

## Phase 2 — Cocoa (academic)

Source: Moreno-Zambrano, Grimbs, Ullrich & Hütt (2018), *R. Soc. Open Sci.*
doi:10.1098/rsos.180964. Data and code: doi:10.5061/dryad.321d33v (Zenodo record
4935403, CC0).

### P2.1 Fetch and inspect original materials 🔎
- **Status:** Complete. Seven original files are extracted byte-for-byte and
  verified by `scripts/verify_cocoa.jl`. The three R-embedded datasets contain
  336 likelihood observations. `docs/cocoa_original.md` records the workflow,
  scaling, priors, inferred noise scale, and a `yc3` rescaling inconsistency
  subsequently confirmed against supplement S3. Source verification passes on
  Julia 1.10.11 and 1.12.6. Original RStan fits have not been run.
- **Output:**
  - `data/raw/cocoa/CocoaModel.zip`, unpacked
  - `docs/cocoa_original.md` listing every file, the language and software used, and
    the datasets included
- **Done when:** we know exactly what the original code does and which data it fits.

### P2.2 Transcribe the model
- **Status:** Transcribed in `docs/cocoa_equations.md`: all eight equations, the 24
  parameter names and dimensional units, source initial conditions, normalized
  priors, observation likelihood, and scaling conversions. Stan and extracted
  paper text agree on the dynamics; rendered supplement S3 pages 7–10 confirm
  the scaling and the `yc3` error in `table.R`. Dimensional units are explicit
  where printed main-table labels disagree. The supplement is preserved with
  DOI/license/hash provenance. Fitted numerical parameters remain P2.3.
- **Output:** `docs/cocoa_equations.md` with all equations, state variables, parameters,
  units, initial conditions and priors from the paper and code, each cross-referenced
  to the paper.
  - States: glucose, fructose, ethanol, lactic acid, acetic acid, yeast, LAB, AAB
  - Uptake kinetics: Monod and Contois
- **Done when:** a human reviewer can check each line against the paper.

### P2.3 Implement the baseline
- **Status:** Equations, coordinate conversion, three sourced parameter sets,
  processed observations, simulation example, and figure comparisons are
  implemented. Analytical, positivity, coordinate-invariance, and derivative
  tests pass. The 2% published-curve criterion remains an explicit expected
  failure: printed-mean simulations differ by up to 7.1% of the calibrated
  vertical range. Original joint draws/plotting code are absent; P2.4 must
  resolve the posterior-prediction comparison. See `docs/cocoa_baseline.md`.
- **Output:** `src/cocoa/baseline.jl`, built on the core engine.
- **Done when:** simulating with the paper's reported parameter estimates reproduces the
  published trajectories within plotting tolerance, as a regression test.

### P2.4 Reproduce the fits
- **Status:** Source likelihood and priors, jointly inferred Gaussian noise,
  reproducible MLE/NUTS runner, retained draws, predictive intervals, and
  parameter/curve comparison reports are implemented. Analytical MLE/profile
  and independent posterior-quadrature checks validate the noise extension.
  The earlier Camu MLE reports native convergence at search bounds and fails the
  published-interval comparison; it has not been reassessed under the new
  stationarity gate. Bayesian fits for all three datasets pass convergence and all 25
  source intervals. Camu uses four chains with 1,000 warmup and 2,000 retained
  draws each; the boxes use two chains with 500 warmup and 1,000 retained draws
  each. All use target acceptance 0.995. A common figure reproduction remains
  unresolved: Camu's mean misses one point at 2.04%, and box 2's median misses
  one at 2.07%. Longer box MLE restarts remain unverified: box 1 returns native
  success but fails an independently recomputed gradient check; box 2 reaches
  its time limit. New MLE acceptance requires projected-gradient stationarity
  in addition to native success. Saved-point checks record bounds, local
  information, finite-difference agreement, and solver/tolerance sensitivity;
  the derivative-accuracy limits remain explicit. Per-start optimizer
  time budgets bound expensive searches alongside ODE step caps. Optional stiff solves,
  concurrent chains within one process, per-start/draw checkpoints, MLE restarts
  from saved parameter points, and trace plots are implemented. The full suites
  on Julia 1.10.11 and 1.12.6 each pass 591 checks with two threads; one documented
  expected failure remains in both. Multi-process fitting exceeded the local
  paging-file limit, so follow-up jobs run sequentially in one Julia process.
  See `docs/cocoa_fitting.md` for results, commands, and diagnostic gates.
- **Output:** a MLE and Bayesian fit to the original datasets.
- **Done when:** estimates fall within the paper's reported credible intervals.
  Discrepancies are documented, not tuned away.

### P2.5 Benchmark set 🔎
- **Status:** Archived source and supplement acquired with checksums; all 32
  screened datasets assembled, with the 23-trial published subset marked.
  Four trials have unmeasured lactic acid, preserved as missing observations
  with the source's fixed normalized hidden initial condition. All 32
  mechanism combinations and parameter counts are catalogued. Published
  baseline diagnostics/BMA weights and a sequential bounded MLE runner are
  implemented, with provenance-checked restarts, stable trial selection, and
  saved numerical diagnostics. The initial bounded survey completes all 32
  trials with finite candidates; none meets convergence/stationarity requirements.
  A 64-case numerical scan covers every pilot point at two tolerance settings.
  Three longer restarts produce two improved, still unaccepted candidates and
  one ODE step-limit failure; no acceptance threshold was relaxed. All 406
  focused benchmark checks pass on Julia 1.10.11 and 1.12.6. Saved-artifact
  audits pass 1,111 MLE checks and 475 Bayesian checks on Julia 1.10.11.
  The full Julia 1.10.11 package suite passes 1,454 checks with one existing
  expected failure.
  A Bayesian runner now records the 2022 priors, retained-draw checkpoints,
  predictive bands, and complete pointwise likelihoods, with saved-file plots.
  Two 2,000-draw pilot fits complete: `brwb3` fails R-hat/divergence checks,
  and hidden-lactate `ghhp4` fails tail ESS. Neither is accepted; thresholds
  remain unchanged.
  A fresh four-chain `ghhp4` refinement with 1,000 warmup and 2,000 retained
  draws per chain passes the same diagnostic rule: maximum R-hat 1.00226,
  minimum bulk/tail ESS 3219.46/1920.08, and no divergences or depth hits.
  Saved-draw numerical checks compare log-coordinate posterior gradients
  across solver tolerances and integrators without changing fit acceptance.
  All 72 cases at 24 pilot points complete; tighter tolerances improve the
  worst gradient discrepancies with NLL changes below `2.88e-5`. Independent
  audits pass 282 refinement checks and 1,709 numerical-report checks.
  All 94 new posterior numerical-diagnostic tests pass on Julia 1.10.11
  and 1.12.6.
  A reusable completed-fit auditor now independently recomputes marginal
  diagnostics, replays Gaussian likelihoods and latent predictive bands,
  and checks source identities, units, checkpoint order, and parent hashes.
  It rejects interrupted experiments and keeps artifact integrity separate
  from sampling acceptance. Measurement-noise draws were not retained,
  so their exported bands receive structural and unit checks only.
  All 136 audit tests pass on Julia 1.10.11 and 1.12.6.
  Full replay of the completed pilots and `ghhp4` refinement verifies all
  12,000 draws, 1,242,000 likelihood values, and 1,500 selected latent
  predictions; their sampling acceptance labels remain unchanged.
  A four-chain `brwb3` attempt at tighter Tsit5 tolerances was interrupted
  with 5,478 partial draws and two retained divergences. It remains
  unaccepted; the incomplete run supplies no joint R-hat/ESS or predictive
  export. A ten-step initialization diagnostic for its fourth chain
  completes without an exception, but does not explain the runtime delay.
  Another 379 checks validate the interrupted experiments, timing records,
  and partial checkpoint integrity without accepting the fit.
  New Bayesian benchmark runs apply source caps to accepted steps between
  output times, preserving the continuous integrator and its adaptive history.
  A separate per-step attempt guard bounds retries. The historical whole-solve
  mode remains explicit, and saved artifacts retain strict implementation hashes.
  Runtime records now include initialization and discarded warmup completions;
  optional soft time limits and retained-divergence stops preserve partial draws.
  At ten archived `brwb3` points, 30 paired solver/tolerance cases give identical
  predictions, likelihoods, and gradients in both budget modes. All cases stay
  below the source cap; 302 comparison checks pass. This does not identify the
  cause of the previous divergences.
  A paired short four-chain pilot exposed a shared callback variable that could
  mix diagnostic fields between chains during checkpoint I/O. The initial logs
  are preserved with a failed runtime check. The corrected callback has a forced
  task-interleaving regression test. A fresh paired pilot has identical draws in
  both modes; 2,898 scientific-export audit checks pass, and an independent audit
  verifies all 600 runtime events against 200 retained draws. Both short cases
  fail the unchanged convergence/depth criteria; this is not a converged refit.
  The preceding budget/callback validation passed 756 focused checks on Julia
  1.12.6. MLE restart and diagnostic tests
  now use disposable synthetic fixtures with fresh scores; separate tests retain
  rejection of historical implementation hashes. A descriptive log-correlation
  scan of the earlier unaccepted `brwb3` pilot is saved for further diagnostics;
  it supplies no convergence claim or explanation of the earlier slow warmup.
  Extended runtime records now expose leapfrog work, step sizes, acceptance,
  and energy errors; an independent raw-chain export validates retained fields.
  A fresh bounded four-chain repeat preserves all nine scientific exports
  byte for byte, passes 1,449 scientific audit checks, and verifies 300 runtime
  events against 100 retained draws. All 12 selected state-cost profiles complete.
  The final two chain-4 states require over 1,100 rejected ODE steps each and
  27–38 ms median gradient evaluations, compared with 2–6 ms at the other selected
  states. Retained states do not reconstruct the sampler path or explain the
  earlier stall. This short pilot remains explicitly unaccepted.
  The telemetry changes pass 377 focused checks on Julia 1.12.6, covering
  callbacks, Bayesian exports, posterior numerics, runtime corruption, and
  cost-profile provenance. The full Julia 1.10.11 suite passes 1,454 checks
  with the one existing expected cocoa figure failure.
  A subsequent 48-case comparison evaluates all 12 selected states with Tsit5
  and Rodas5P at `1e-8` and `1e-10`. All numerical screens and 504 Julia integrity
  checks pass; an independent audit verifies 100 hashes and all 240 timing rows.
  Rodas5P `1e-8` is 5.10–5.61 times slower at the two expensive states and has
  a median paired gradient-time ratio of 8.36 across the subset. Neither alternate
  configuration passes the predeclared performance screen, so no alternate-solver
  pilot is launched. Tsit5 `1e-8` remains unchanged.
  A subsequent bounded pilot increases discarded calls from 50 to 250 and
  retained draws from 25 to 100 per chain, with the same four seeds, depth 8,
  target acceptance, solver, data, and priors. All 400 retained draws complete
  without divergences, within the 1,800-second sampling limit. Maximum R-hat
  improves to 1.07215, but minimum bulk/tail ESS is only 47.49/39.09 and 387
  draws hit the depth limit (96.75%, up from 80%). The fit remains unaccepted.
  All 3,249 scientific-export checks pass, replaying 38,400 likelihood values
  and 50 latent predictions; runtime validation checks 1,400 events against
  all 400 raw retained statistics. An independent audit verifies 165 hashes.
  The launcher's native exit code was not captured; the original receipt is
  preserved, and the revised audit explicitly distinguishes verified scientific
  completion from an unknown process exit code. A separate launcher probe
  validates handle retention for future runs with both zero and nonzero exits.
  The matched depth-9 comparison subsequently completes all 400 retained draws
  without divergences. Depth hits fall from 387 to 133, maximum R-hat improves
  to 1.05479, and minimum bulk/tail ESS rises to 96.72/74.07. Retained leapfrog
  work increases by 51%; minimum bulk/tail ESS per leapfrog improves by 35%/26%.
  All conditions of the predeclared longer-run screen pass, but the fit still
  fails the unchanged convergence and depth criteria. Sampling takes 2,144.1
  seconds within its 2,700-second limit. The driver passes 19 integrity checks,
  the scientific replay passes 3,249 checks, and an independent audit verifies
  168 hashes, 1,400 runtime events, all 400 raw retained statistics, and the
  screen decision. The corrected launcher captures native exit code 0.
  Next, use depth 9 in a separately bounded longer four-chain fit; preserve
  the divergence stop, runtime records, and unchanged acceptance criteria.
  Converged benchmark fits and a like-for-like PSIS-LOO comparison
  remain open; original posterior outputs are absent from the source archive.
  See `docs/cocoa_benchmark.md`.
- **Source:** Moreno-Zambrano, Ullrich & Hütt (2022), "Exploring cocoa bean
  fermentation mechanisms by kinetic modelling", PMC8847890: 32 model
  variants screened on 32 datasets; 23 retained after at least one successful fit.
- **Output:** the 23 datasets assembled into `data/processed/cocoa_benchmark/` (from
  supplements, or digitized with provenance), and their mechanism variants listed.
- **Done when:** the baseline is fitted to all datasets that can be obtained, with a
  table of fit quality vs. the 2022 paper.

### P2.6 Extensions (pick later)
- Arrhenius or cardinal temperature dependence (the author's thesis explores Arrhenius terms)
- Pulp/bean temperature as a state, rather than an input
- Hierarchical model across trials: shared kinetics, trial-specific environment

---

## Phase 3 — Sourdough (kitchen)

**Questions:** (a) When will my starter peak after a 1:N:N feed at temperature T?
(b) How do T and feed ratio shift yeast vs. bacteria, and so sourness?

Anchor source: Gänzle, Ehmann & Hammes (1998), AEM 64:2616, PMC106434.

### P3.1 Extract Gänzle parameters 🔎
- **Output:** `params/sourdough_ganzle1998.toml` with growth-rate dependence on
  temperature, pH, NaCl, lactate, acetate and ethanol for *L. sanfranciscensis* and
  *C. milleri* (now *Kazachstania humilis*).
- **Known qualitative values** (use as test expectations):
  - LAB grow at pH 3.9–6.7, are fully inhibited by 4 % NaCl, and tolerate more than
    160 mmol/L undissociated acetic acid
  - Yeast optimum is 27 °C, growth is unaffected by pH 3.5–7, it tolerates up to 8 %
    NaCl, and is fully inhibited at 150 mmol/L undissociated acetic acid
- **Done when:** every parameter has a table/equation reference, or is marked missing.

### P3.2 Two-organism starter model
- **States:**
  - yeast `Y` and LAB `L`
  - fermentable sugar `S` (maltose/glucose, lumped)
  - lactic `La`, acetic `Ac`, ethanol `E`, CO₂ `C`
- **Rise model:** a simple gas-retention model, `height ∝ retained CO₂`, with leakage
  that grows once production slows. Hydration fixed at 100 % initially.
- **Done when:**
  - Qualitative checks pass: warmer is faster; above about 29 °C yeast lags and LAB
    dominate; larger N gives a later peak
  - Sugar exhaustion or acid inhibition ends the rise

### P3.3 Feeding cycles
- **Output:** dilution events (P1.4) for repeated 1:N:N feeds; periodic steady state
  per (N, T).
- **Done when:** `peak_time(N, T)` and `acidity_at_peak(N, T)` return values over a
  grid of N ∈ {1…10} and T ∈ {16…30 °C}.

### P3.4 Calibration against baker data 🔎
- **Data:**
  - **Community rise-time tables (low trust, approximate).** Example: a 10 % inoculation
    takes about 12.5 h mix-to-bake at 65 °F and under 8 h at 80 °F (The Fresh Loaf,
    comment 27403).
  - **Landis et al. 2021 rise-rate data** (Dryad doi:10.5061/dryad.0p2ngf1z1, about 40
    starters) for between-starter spread.
- **Output:** a fitted scale factor or a small set of free parameters (e.g. initial
  population, gas retention), with the Gänzle kinetics as priors.
- **Done when:** predictions vs. tables are plotted, with residuals and honest error bars.

### P3.5 Sensitivity analysis
- **Output:** which inputs move peak time most (T, N, hydration, flour sugar).
- **Done when:** a figure plus a short write-up in `docs/` are committed.

### P3.6 Kitchen calculator v1
- **Output:** a Pluto notebook with sliders (T, N, hydration) showing a predicted rise
  curve with an uncertainty band and peak time.
- **Done when:** it runs from a clean environment.

---

## Phase 4 — Vegetables (kitchen)

**Questions:** given salt % and temperature, (a) when is pH ≤ 4.6? (b) when does it reach
target acidity / tartness? (c) how many days for a 5-log reduction of *E. coli*?

### P4.1 Literature data harvest 🔎
Collect pH, titratable acidity and organism time series under different salt and
temperature conditions. Candidates:
- Zabat et al. 2018 (sauerkraut at about 21 °C, *Foods* 7:77)
- Kimchi time series (Jung et al. 2011; Kim et al. 2022, 500 days at 4 °C)
- USDA-ARS Raleigh cucumber papers (Breidt, McFeeters, Johanningsmeier)
- Classic sauerkraut temperature and salt studies (Pederson & Albury, NY Agric. Exp.
  Station; locate and verify)

- **Output:** `data/processed/veg/` tidy CSV `(source_id, vegetable, salt_pct, temp_C,
  time_h, variable, value, unit)`, with a provenance README per source.
- **Done when:** at least 3 temperatures and 2 salt levels are covered, or the gap is
  documented.

### P4.2 Buffer model for vegetable juice 🔎
- **Output:** pH as a function of lactic and acetic acid, following the Breidt group's
  buffer-capacity models for cucumber juice. Parameterize cabbage separately if data exist.
- **Done when:** predicted pH vs. measured pH, given measured acids, is within about
  0.15 units on the harvested data.

### P4.3 Succession model
- **Organisms:**
  - heterofermentative LAB `H` (*Leuconostoc*/*Weissella*, early, less acid-tolerant)
  - homofermentative LAB `O` (*L. plantarum*, late, acid-tolerant)
- **States:** sugar `S`, lactic `La`, acetic `Ac`, ethanol, CO₂.
- **Behaviour targets from the literature:**
  - Low salt (about 2 %) and cool (about 18 °C) favour heterofermenters
  - High salt (about 3.5 %) and warm (about 32 °C) favour homofermenters
  - Finished kraut reaches about 1.4–1.6 % acid (as lactic)
- **Done when:** the qualitative targets are met and the model is fitted to P4.1 data
  with P1.5/P1.6.

### P4.4 Safety layer
- **Output:** an *E. coli* inactivation submodel (log-linear or Weibull) whose rate
  depends on pH, undissociated acid and T. Calibrate it to Breidt et al. data on
  O157:H7 survival in fermented cucumber brine; in their fermented-juice medium at
  pH 3.9, a 5-log reduction took 1–16 days depending on conditions.
- **Done when:** the submodel outputs `days_to_5log(salt, T)` with an uncertainty band,
  and the docs clearly label it educational.

### P4.5 Kitchen calculator v1
- **Output:** a Pluto notebook with inputs (salt %, room T, optional fridge date) and
  outputs: predicted pH and acidity curves, day pH ≤ 4.6, day reaching target
  acidity, day of 5-log reduction, all with bands.
- **Done when:** it runs from a clean environment, and the warnings are visible on the page.

---

## Phase 5 — Later

- **Kitchen data loop.** Log jar height (rubber band + photo) and pH-strip readings.
  Fit a hierarchical model: literature posterior becomes the prior, which kitchen data
  update. A Pluto or web form for logging.
- **Calculator outside Julia.** Port the fitted sourdough and vegetable models to a
  small static web page. The models are small enough for plain JS, or Rust + diffsol
  compiled to WASM (see TOOLING.md).
- **Community dynamics.** Fit generalized Lotka–Volterra or consumer-resource models to
  sequencing time series (kimchi, sauerkraut, sourdough). This needs
  compositional-data methods, since amplicon data are relative abundances.
