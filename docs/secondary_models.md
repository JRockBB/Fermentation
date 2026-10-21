# Secondary growth models (P1.1)

The implementation in `src/core/secondary.jl` uses equations from the papers,
without copying code from GPL references. Temperature uses °C, pH is dimensionless,
and each response scales an independently supplied optimum growth rate.

| API | Equation source |
|---|---|
| `γ_T_ctmi(T; Tmin, Topt, Tmax)` | Rosso et al. (1995), equation 6, reproducing the 1993 CTMI |
| `γ_pH_cardinal(pH; pHmin, pHopt, pHmax)` | Rosso et al. (1995), equation 7 |
| `γ_T_ganzle(T; a, b, c, Tref)` | Gänzle et al. (1998), Table 1 model 1a and Appendix |
| `γ_acid(HA; MIC, α=1)` | Gänzle et al. (1998), Table 1 model 4 |
| `γ_salt(NaCl; max, α=1)` | Same power-law inhibition family applied to salt mass fraction |

Sources: [Rosso et al., doi:10.1128/aem.61.2.610-616.1995](https://pbil.univ-lyon1.fr/members/lobry/articles/AEM_1995_61_610.pdf)
and [Gänzle et al., doi:10.1128/aem.64.7.2616-2623.1998](https://journals.asm.org/doi/10.1128/AEM.64.7.2616-2623.1998).

The Gänzle formula in the initial plan had the wrong exponential sign. The paper
uses `a*x^b*exp(-c*x)` with positive coefficients and `x = Tref - T`.
Normalizing at `x = b/c` cancels `a`. This function has an upper cardinal limit,
but no finite lower limit; no artificial lower cutoff is introduced.

For CTMI, ordered cardinal temperatures alone are insufficient. Algebraically,
when the optimum lies below the midpoint the denominator crosses zero inside the
growth interval. The implementation rejects those parameter sets. Inhibition
models require nonnegative concentrations and positive limits and exponents.
All formulas return zero at and above their upper growth/inhibition limit.
Clipped boundaries are nonsmooth; derivative tests use interior points.

The salt function assumes maximal growth at zero salt and does not represent
salt-stimulated growth. Published ionic-strength parameters cannot be substituted
for salt mass percentages without a documented conversion and calibration.
No organism-specific numerical parameters are supplied here.

To regenerate the synthetic shape-review plots:

```sh
julia --project=scripts -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=scripts scripts/secondary_figures.jl
```

Outputs: `docs/figures/secondary_models.png` and `.svg`. Plotting dependencies
live in the separate `scripts/` environment.
