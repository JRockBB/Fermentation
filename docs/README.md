# FermentModels

This document records the original project outline. See the
[repository README](../README.md) for current setup and implementation status.

Mechanistic models of food fermentations, built for two audiences at once:

1. **The kitchen** — practical predictions for home sourdough and vegetable ferments
   ("I fed 1:5:5 at 21 °C, when does it peak?", "2 % salt at 20 °C, when is my kraut
   safe and when is it tart?").
2. **The notebook** — an academic reproduction and extension of a published cocoa bean
   fermentation model, which doubles as the test bed for the shared modeling engine.

All three tracks share one engine: microbial populations consume sugars, produce acids
(and ethanol / CO₂), and are slowed by temperature, salt, pH and accumulated acid. Only
the species, substrates and parameter values change between foods.

## Ground rules

- **No home data to start.** Every parameter is taken from, or fitted to, published data.
  Each value carries its source (see `AGENTS.md` → provenance). Kitchen data can be added
  later to recalibrate.
- **Mechanistic first.** Small, interpretable ODE systems (3–10 states). Machine learning
  only as a later add-on, if ever.
- **Predictions are ballpark.** Real flours, cabbages and wild microbes differ from any
  published system. Outputs should always carry uncertainty bands.
- **Not food-safety certification.** The vegetable safety estimates are educational. The
  pH ≤ 4.6 threshold and published inactivation data inform them; they do not replace
  tested recipes or a pH meter.

## Tracks

| Track | Purpose | Status |
|---|---|---|
| **Core engine** | Secondary models (temperature, pH, salt, acid inhibition), acid/pH chemistry, ODE assembly, events, fitting | planned |
| **Cocoa** (academic) | Reproduce Moreno-Zambrano et al. 2018; benchmark on 23 literature datasets; extend | planned |
| **Sourdough** (kitchen) | Starter peak time and sourness vs. feed ratio and temperature | planned |
| **Vegetables** (kitchen) | Time to safe pH and to target acidity vs. salt % and temperature | planned |

Work order: core → cocoa (validates the engine against real data with released code) →
sourdough → vegetables. Details in [`PLAN.md`](PLAN.md).

## Language and libraries

Julia, built on OrdinaryDiffEq.jl (solving), Optimization.jl (point fits) and Turing.jl
(Bayesian fits). Python, Rust and R libraries are used only as cross-checks or
references. Reasoning in [`TOOLING.md`](TOOLING.md).

## Repository layout (target)

```
FermentModels/
├── README.md
├── AGENTS.md                 # instructions for coding agents (symlink as CLAUDE.md if needed)
├── Project.toml / Manifest.toml
├── docs/
│   ├── PLAN.md               # phased, agent-sized tasks
│   ├── TOOLING.md            # library assessment
│   └── SOURCES.md            # papers & datasets, licenses, trust levels
├── src/
│   ├── FermentModels.jl
│   ├── core/                 # secondary models, chemistry, assembly, events, fitting
│   ├── cocoa/
│   ├── sourdough/
│   └── vegetables/
├── params/                   # TOML parameter sets, every value with provenance
├── data/
│   ├── raw/                  # downloaded or digitized, never edited
│   └── processed/            # cleaned, tidy CSV, produced by scripts
├── scripts/                  # data fetching, digitizing notes, figure generation
├── notebooks/                # Pluto notebooks (exploration, interactive calculators)
└── test/
```

## Quickstart (once scaffolded)

```julia
julia --project -e 'using Pkg; Pkg.instantiate()'
julia --project -e 'using Pkg; Pkg.test()'
```
