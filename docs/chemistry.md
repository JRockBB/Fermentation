# Acid and pH chemistry (P1.2)

`undissociated_fraction(pH, pKa)` implements monoprotic acid speciation, and
`undissociated_acid(total, pH; pKa)` returns HA in the input concentration units.
`henderson_hasselbalch(acid, base; pKa)` is the simple fallback when both HA and
its conjugate base are known.

`pH(acids, buffer_params)` solves an ideal monoprotic buffer plus added neutral
acids. Its charge-balance structure follows the acid-buffer case of equation 3
in [Breidt and Skinner (2022), doi:10.4315/JFP-22-068](https://www.ars.usda.gov/ARSUserFiles/60702500/Fermented%20Veg%20Manuscripts/p450.pdf).
Their full food-specific model uses multiple effective buffers, measured
titration data, and ionic-strength corrections. Those calibrations are not
included in this generic single-buffer implementation.

All concentrations are mmol/L. Let `h` denote pH and define:

```text
H(h)  = 10^(3-h)
OH(h) = 10^(3+h-pKw)
f(h, pKa) = 1/(1 + 10^(pKa-h))
Q = B*f(pH0, pKa_buffer) + OH(pH0) - H(pH0)
0 = H(h) + Q - OH(h) - B*f(h, pKa_buffer) - sum(C_i*f(h, pKa_i))
```

`Q` represents the fixed counterion charge needed to reproduce the supplied
initial pH before the listed fermentation acids are added. The residual strictly
decreases with pH, giving a unique root when bracketed. The caller supplies pKas
and pKw for the relevant conditions; there are no hidden food-specific constants.
Acid salts, polyprotic acids, ionic activities, and changing buffer composition
need extensions before use with the corresponding food data.

The buffer parameters are a NamedTuple with `total`, `pKa`, `pH0`, `pKw`, and
`acid_pKas`. Alternatively, `buffer_from_capacity` converts a measured **local,
buffer-only** capacity at `pH0` (mmol/L per pH unit) into a total monoprotic buffer
concentration. It is not the same as an integrated titration to a chosen endpoint.

The root is bracketed by bisection and refined by Newton steps. Refinement also
propagates the implicit derivative through the solution; tests compare ForwardDiff
with finite differences and a hand-derived sensitivity. Titration tests invert
the charge balance at known target pHs and check that adding acid lowers pH.
