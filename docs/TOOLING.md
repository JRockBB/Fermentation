# Tooling assessment

Survey as of October 2026. Re-check versions before pinning.

## What this project needs

1. **Small ODE systems:** 3–10 states, 5–30 parameters, sometimes mildly stiff
   (populations near zero, sharp acid inhibition).
2. **Events and inputs:** feeding/dilution events, time-varying temperature.
3. **Point fits:** multistart maximum likelihood with bounds.
4. **Bayesian fits:** NUTS with ODE gradients. The cocoa paper is Bayesian, and kitchen
   predictions need uncertainty bands. Hierarchical models come later.
5. **Interactive exploration:** plots, notebooks with sliders for the kitchen calculators.
6. **Coding-agent friendliness:** stable APIs, plenty of examples, clear errors.
7. **Optional later:** a small calculator that runs outside the research stack (web page / phone).

## Julia (chosen)

| Library | Role | Notes |
|---|---|---|
| **OrdinaryDiffEq.jl** | ODE solving | `Tsit5` (non-stiff), `Rodas5P` / `FBDF` (stiff). Callbacks handle feedings (`PresetTimeCallback`) and positivity. Load this rather than the full DifferentialEquations.jl. |
| **Optimization.jl** (+ OptimizationOptimJL) | MLE | LBFGS with ForwardDiff gradients is plenty at this model size. |
| **Turing.jl** | Bayesian | SciML documents the Turing + ODE workflow directly. DiffEqBayes.jl offers a simplified wrapper that can also call Stan or DynamicHMC. |
| **PEtab.jl** | Multi-condition fitting (optional) | Links ODE models to measurement data across experimental conditions, with multistart via Optim, Ipopt or Fides, plus priors/posteriors. A strong fit for the 23-dataset cocoa benchmark. Requires a ModelingToolkit or Catalyst model, so adopt it only if that pays off. |
| **ModelingToolkit.jl / Catalyst.jl** | Symbolic modeling (deferred) | Catalyst gives a reaction-network DSL on MTK and benchmarks well (PLOS Comp Bio 2023). But MTK v10 renamed core APIs (e.g. `structural_simplify` → `mtkcompile`). MTK v11 (Dec 2025) split into an MIT **ModelingToolkitBase** and an AGPL-dependent **ModelingToolkit**. API churn trips up agents, and our systems are small, so we write RHS functions by hand. If needed later, prefer ModelingToolkitBase. |
| **Pigeons.jl** | Multimodal posteriors (optional) | Parallel tempering. Use it if NUTS struggles with the multi-population cocoa model. |
| **Pluto.jl** | Notebooks / calculators | Reactive sliders, good for the kitchen calculators. |
| **CairoMakie** or **Plots.jl** | Figures | Either works. |
| **CSV.jl, DataFrames.jl, TOML (stdlib)** | Data and parameters | |

**Gaps.**
- No mature Julia package covers food predictive microbiology (cardinal models, Baranyi
  primary models). We implement those ourselves (P1.1); the code is small.
- Julia start-up and compile time is the main friction. Mitigate with OrdinaryDiffEq
  instead of DifferentialEquations, and a persistent REPL or Pluto session.

## Python (reference and cross-checks)

| Library | Fit for us | Notes |
|---|---|---|
| **SciPy `solve_ivp`** | Prototyping only | Simple, but gives no gradients and is slow inside MCMC loops. |
| **Diffrax (JAX) + NumPyro** | Strongest Python alternative | The solver is part of the JAX graph, so NUTS gets exact gradients. JIT compile makes repeated solves fast. |
| **Stan via cmdstanpy** | Best cross-check | Mature built-in ODE solvers (RK45, BDF, adjoint). Stan's sampler is a trusted reference. Use it to confirm one Julia posterior (e.g. cocoa baseline). |
| **PyMC (+ sunode / nutpie)** | Viable | nutpie is a fast NUTS sampler for PyMC and Stan models, built on the Rust crate nuts-rs. |
| **pyPESTO + AMICI** | Viable, heavier | Systems-biology estimation toolbox: multistart, profiles, sampling, AMICI's C++ simulation and sensitivities. Shares the PEtab problem format with PEtab.jl, so the same problem can be run in both. |
| **Catalax** | Interesting, young | JAX framework for biochemical ODEs with neural ODEs and HMC inference. |
| COBRApy / MICOM | Out of scope | Genome-scale metabolic models; only if we go down the yogurt/community-FBA path. |

## R (reference only)

- **biogrowth** (Garre et al., Wageningen / Cartagena): primary and secondary predictive
  microbiology models, fitting and prediction under static and dynamic conditions, MCMC.
  It is the best existing reference for our P1.1 secondary models. **GPL-3:** compare
  numbers against it, do not copy code.

## Rust (deferred)

| Library | Notes |
|---|---|
| **diffsol** | Pure-Rust ODE/DAE solver: BDF, SDIRK/ESDIRK, explicit RK. Has forward sensitivities, adjoints, and a DSL (DiffSL). Designed to be called from Python and R. |
| **nuts-rs** | Rust NUTS implementation behind nutpie. |
| **argmin** | Optimization framework. |

**Assessment.** The pieces exist, but nothing ties model, likelihood, gradients and
sampler together the way Turing or Stan do. We would hand-wire every likelihood.
That's not worth it for the research core. Rust's realistic niche here is Phase 5: a
fitted model compiled to WASM for an offline calculator page. Even then, these models
are small enough that plain JavaScript may be simpler.

## Scorecard

Scores are 1–5, higher is better, against this project's needs.

| Criterion | Julia (SciML + Turing) | Python (Diffrax + NumPyro) | Python (Stan) | Rust (diffsol + nuts-rs) |
|---|---|---|---|---|
| ODE solver quality, stiffness, events | 5 | 4 | 4 | 4 |
| Gradients through solver | 5 | 5 | 5 | 3 (manual wiring) |
| Bayesian workflow | 4 | 5 | 5 | 2 |
| Multistart MLE / profiles | 4 | 3 | 3 | 2 |
| Interactive notebooks | 5 (Pluto) | 4 (Jupyter) | 3 | 1 |
| API stability for agents | 3 (avoid MTK churn) | 4 | 5 | 3 |
| Your preference | 5 | 3 | 3 | 2 |
| Deployable calculator | 2 | 2 | 1 | 5 (WASM) |

## Decision

- **Core:** Julia with OrdinaryDiffEq.jl, Optimization.jl and Turing.jl, using hand-written RHS functions.
- **Cross-check:** reproduce the cocoa baseline posterior once in Stan (cmdstanpy) to
  validate the Julia pipeline.
- **Reference:** R biogrowth for secondary-model numerics.
- **Revisit** PEtab.jl (with ModelingToolkitBase) at P2.5 if multi-condition fitting
  gets painful. Revisit Rust/WASM at Phase 5.

## Sources

- SciML Turing + ODE tutorial: https://docs.sciml.ai/SciMLSensitivity/stable/examples/bayesian/turing_bayesian
- SciML parameter estimation overview: https://docs.sciml.ai/Overview/dev/highlevels/inverse_problems/
- Catalyst paper (PLOS Comp Bio 2023): https://pmc.ncbi.nlm.nih.gov/articles/PMC10584191
- ModelingToolkit v10 release notes: https://newreleases.io/project/github/SciML/ModelingToolkit.jl/release/v10.0.0
- ModelingToolkit v11 licensing post: https://sciml.ai/news/2026/07/27/ModelingToolkitv11/index.html
- PEtab.jl optimizers: https://sebapersson.github.io/PEtab.jl/stable/tutorials/parameter_estimation/optimizers
- PEtab.jl / SBMLImporter.jl preprint: https://www.biorxiv.org/content/10.1101/2025.04.30.651378.full.pdf
- pyPESTO (Bioinformatics 2023): https://pmc.ncbi.nlm.nih.gov/articles/PMC10689677
- nutpie: https://github.com/pymc-devs/nutpie
- Diffrax + NumPyro in practice: https://arxiv.org/pdf/2606.24966
- Catalax: https://pypi.org/project/catalax/
- biogrowth (CRAN): https://www.rdocumentation.org/packages/biogrowth/versions/1.0.8
- diffsol: https://docs.rs/crate/diffsol/latest/source/src/lib.rs
