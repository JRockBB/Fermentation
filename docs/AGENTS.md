# Instructions for coding agents

Read `README.md` and `docs/PLAN.md` first. Work on one task ID from the plan at a time
and state which one at the top of your PR or summary.

## Hard rules

1. **Never invent parameter values.** Every number in `params/` needs a source: DOI,
   table/figure/equation, and units. If a value is not available, leave it as
   `missing` with a `TODO` note and say so in your summary. Fitted values must cite the
   fit script and dataset that produced them.
2. **Never edit `data/raw/`.** Raw files are as downloaded or as digitized. All
   cleaning happens in scripts that write to `data/processed/`.
3. **Record provenance for digitized data.** Digitized points (e.g. WebPlotDigitizer)
   go in `data/raw/digitized/<source-id>/` with a `README.md` giving the paper, figure
   number, axis calibration and the person or agent who digitized them.
4. **Respect licenses.** Cocoa data is CC0. R `biogrowth` is GPL-3: use it as a
   numerical reference only and do not copy its code into this repo.
5. **Every model change ships with tests** (see Testing).

## Julia conventions

- Use **OrdinaryDiffEq.jl**, not the full DifferentialEquations.jl, to keep load times down.
- **Do not use ModelingToolkit.jl or Catalyst.jl** unless the task explicitly says so.
  Write right-hand-side functions by hand (`f!(du, u, p, t)`). MTK had breaking
  releases in v10 and v11, and code written against older versions is often wrong.
- Parameters are passed as `NamedTuple`s or ComponentArrays, loaded from TOML in `params/`.
- Pure functions for secondary models (`γ_T(T, p)`, `γ_pH(pH, p)` …) so they can be
  unit-tested and plotted alone.
- Make types AD-friendly (no hard-coded `Float64` in RHS), so ForwardDiff works through
  the solver.
- Pin package versions in `Project.toml` `[compat]`.

## Units (internal)

| Quantity | Unit |
|---|---|
| time | hours (h); convert to days only in user-facing output |
| temperature | °C |
| concentrations | mmol/L (acids, sugars) unless a source model uses g/L, in which case convert at load time |
| populations | log10 CFU/g or CFU/mL, state the basis in the parameter file |
| salt | % w/w NaCl on total vegetable + salt mass |

Convert at the boundary (loading params, printing results), never inside RHS functions.

## Parameter file format (`params/*.toml`)

```toml
[lab_mu_max]
value = 0.62
units = "1/h"
source = "doi:10.1128/aem.64.7.2616-2623.1998"
location = "Table 2"
notes = "L. sanfranciscensis LTH2581, 30 °C, pH 5.6"
confidence = "high"   # high | medium | low | placeholder
```

(The value above is illustrative only. Extract the real number from the paper.)

## Testing

- **Unit tests:** each secondary model returns 0 outside its cardinal limits, reaches 1
  at its optimum, and stays within [0, 1].
- **Invariants:** no negative concentrations; carbon or mass balance holds within
  tolerance where the model defines stoichiometry.
- **Synthetic recovery:** simulate with known parameters, add noise, fit, and check that
  the true values lie inside the 95 % interval.
- **Regression:** published reproductions (e.g. cocoa fits) are pinned as reference
  trajectories with tolerances.

Run with `julia --project -e 'using Pkg; Pkg.test()'`.

## When you are unsure

Stop and write the question into your summary rather than guessing. This applies
especially to biology (which organism does what) and to which paper's equations to follow.
