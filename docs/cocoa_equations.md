# Cocoa baseline transcription (P2.2)

This is a reviewable transcription of the archived implementation, with paper
equation references and explicit distinctions between physical and scaled
coordinates. The Julia implementation and its outstanding reproduction check
are documented in [P2.3](cocoa_baseline.md).

Primary sources:

- Moreno-Zambrano et al. (2018), doi:10.1098/rsos.180964,
  [paper PDF](https://upload.wikimedia.org/wikipedia/commons/5/50/A_mathematical_model_of_cocoa_bean_fermentation.pdf),
  Tables 1–2, equations (2.1)–(2.8), and sections 2.2 and 2.4.
- [Archived Stan model](../data/raw/cocoa/CocoaModel-master/simulation/model.stan),
  whose bytes are pinned in [the inventory](cocoa_archive_files.csv).
- The three archived R drivers and `results/table.R`, described in
  [the original-materials audit](cocoa_original.md).
- [Supplement, version 2](../data/raw/cocoa/rsos180964_si_001.pdf),
  doi:10.6084/m9.figshare.7152614.v2, section S3, pages 7–10.

The equations below were checked against the Stan source, the paper's text
extraction, and rendered supplement S3 pages 7–10. S3 confirms the conversion
factors below and identifies a rescaling error in the archived table script.
The main paper's printed unit labels were available only through text extraction;
their disagreement with the equations is retained explicitly.

The repository's cocoa "baseline" means the final eight-state model, called
**M4** in supplement S2 (Table S2, page 4). It is distinct from the supplement's
earlier seven-state model named M1.

## State coordinates and initial conditions

Time is in hours. Physical state concentrations are mg of the named material per
g of pulp; microbial states are dry-biomass equivalents, not CFU or log CFU.

| Index | Symbol here | Source abbreviation | Material |
|---:|---|---|---|
| 1 | `G` | Glc | Glucose |
| 2 | `F` | Fru | Fructose |
| 3 | `E` | EtOH | Ethanol |
| 4 | `L` | LA | Lactic acid |
| 5 | `A` | Ac | Acetic acid |
| 6 | `Xy` | Y | Yeast biomass |
| 7 | `Xl` | LAB | Lactic-acid-bacterial biomass |
| 8 | `Xa` | AAB | Acetic-acid-bacterial biomass |

The initial state is fixed at `t0=0`. Each number below is copied from line 22 of
the named R driver, in the preceding state order and source units:

```text
camu.R:
[51.963, 57.741, 0, 0, 0, 0.029180401, 0.007868827, 3.36634e-6]

p_box1.R:
[55.482, 49.669, 0, 0, 1.379, 0.175424909, 0.003047264, 2.86522e-6]

p_box2.R:
[42.936, 67.249, 0.318, 0, 0.942, 0.004754351, 0.000221264, 3.86508e-6]
```

The fitting arrays exclude this initial row. Observation times and sample counts
are listed in [P2.1](cocoa_original.md#datasets-actually-fitted). The verification
script reads these numerical literals without evaluating R.

## Growth and mortality rates

Use the source parameter names below. In this section they denote **physical**
parameters; the conversion from sampled normalized parameters appears later.
Each rate has units mg of its microbial group per g pulp per hour.

```text
v1 = mu1 * G/(ks1 + G) * Xy
v2 = mu2 * F/(ks2 + F) * Xy
v3 = mu3 * G/(ks3 + G) * Xl
v4 = mu4 * E/(ks4 + E) * Xa
v5 = mu5 * L/(ks5*Xa + L) * Xa

v6 = k1 * Xy * E
v7 = k2 * Xl * L
v8 = k3 * Xa * A^2
```

These correspond to paper Table 1 and the terms in `model.stan:26–33`.
`v1`–`v4` are Monod routes; `v5` is Contois. The two yeast routes and the two
AAB routes are additive. They must not become a product of substrate limitation
factors. At `Xa=L=0`, the literal Contois expression is `0/0`; a Julia extension
should define the nonnegative continuous extension as zero and test it. This is
a numerical boundary convention, not an additional biological mechanism.

## Eight state equations

| Paper equation | State derivative | Stan RHS line |
|---|---|---:|
| (2.1) | `dG/dt = -yc1*v1 - yc2*v3` | 26 |
| (2.2) | `dF/dt = -yc3*v2` | 27 |
| (2.3) | `dE/dt = yc4*v1 + yc5*v2 - yc6*v4` | 28 |
| (2.4) | `dL/dt = yc7*v3 - yc8*v5` | 29 |
| (2.5) | `dA/dt = yc9*v3 + yc10*v4 + yc11*v5` | 30 |
| (2.6) | `dXy/dt = v1 + v2 - v6` | 31 |
| (2.7) | `dXl/dt = v3 - v7` | 32 |
| (2.8) | `dXa/dt = v4 + v5 - v8` | 33 |

The archived implementation treats these yield coefficients as independent
positive parameters. Its state set omits products such as CO2 and dead biomass;
the sum of these eight concentrations is not a specified conservation law.
No temperature, pH, lag, feeding, or oxygen state appears in these equations.

## Parameter dictionary and units

This mapping follows `model.stan:88–112` and paper Tables 1–2. Units are written
dimensionally from the equations and source state basis. The final column is
the index in the Stan solver vector `theta`.

| Source name | Meaning | Physical units | `theta` |
|---|---|---|---:|
| `mu1` | Yeast growth on glucose | h^-1 | 12 |
| `mu2` | Yeast growth on fructose | h^-1 | 13 |
| `mu3` | LAB growth on glucose | h^-1 | 14 |
| `mu4` | AAB growth on ethanol | h^-1 | 15 |
| `mu5` | AAB growth on lactate | h^-1 | 16 |
| `ks1` | Yeast/glucose Monod constant | mg glucose / g pulp | 17 |
| `ks2` | Yeast/fructose Monod constant | mg fructose / g pulp | 18 |
| `ks3` | LAB/glucose Monod constant | mg glucose / g pulp | 19 |
| `ks4` | AAB/ethanol Monod constant | mg ethanol / g pulp | 20 |
| `ks5` | AAB/lactate Contois constant | mg lactate / mg AAB | 21 |
| `k1` | Ethanol-dependent yeast mortality | g pulp / (mg ethanol h) | 22 |
| `k2` | Lactate-dependent LAB mortality | g pulp / (mg lactate h) | 23 |
| `k3` | Squared-acetate-dependent AAB mortality | (g pulp)^2 / ((mg acetate)^2 h) | 24 |
| `yc1` | Glucose consumed per yeast growth via `v1` | mg glucose / mg yeast | 1 |
| `yc2` | Glucose consumed per LAB growth via `v3` | mg glucose / mg LAB | 2 |
| `yc3` | Fructose consumed per yeast growth via `v2` | mg fructose / mg yeast | 3 |
| `yc4` | Ethanol produced per yeast growth via `v1` | mg ethanol / mg yeast | 4 |
| `yc5` | Ethanol produced per yeast growth via `v2` | mg ethanol / mg yeast | 5 |
| `yc6` | Ethanol consumed per AAB growth via `v4` | mg ethanol / mg AAB | 6 |
| `yc7` | Lactate produced per LAB growth via `v3` | mg lactate / mg LAB | 7 |
| `yc8` | Lactate consumed per AAB growth via `v5` | mg lactate / mg AAB | 8 |
| `yc9` | Acetate produced per LAB growth via `v3` | mg acetate / mg LAB | 9 |
| `yc10` | Acetate produced per AAB growth via `v4` | mg acetate / mg AAB | 10 |
| `yc11` | Acetate produced per AAB growth via `v5` | mg acetate / mg AAB | 11 |

**Unit discrepancies:** the text extracted from printed Table 2
labels the Contois constant as mg lactate/g pulp, although multiplication by AAB
concentration in the denominator requires a lactate-to-AAB mass ratio. Its
mortality labels omit pulp-mass factors required by the source concentration
basis. The table above records dimensional deductions rather than silently
copying those labels. Supplement S3, page 8, independently confirms the
Contois conversion `ks5_z*s_L/s_Xa` and the mortality conversions used below.

## Normalization and inverse parameter conversion

For each dataset separately let `s_i = max(x0_i, x_1i, ..., x_Ti)` and
`z_i = x_i/s_i`. This is exactly the R `scl` calculation and Stan's transformed
data block. Denote a sampled coefficient in normalized coordinates by a suffix
`_z`. Substitution into `dx_i/dt = s_i * dz_i/dt` gives:

| Coefficient | Conversion to physical coordinates |
|---|---|
| `mu1`–`mu5` | unchanged |
| `ks1`, `ks3` | `ks_z * s_G` |
| `ks2` | `ks2_z * s_F` |
| `ks4` | `ks4_z * s_E` |
| `ks5` | `ks5_z * s_L/s_Xa` |
| `k1` | `k1_z / s_E` |
| `k2` | `k2_z / s_L` |
| `k3` | `k3_z / s_A^2` |
| `yc1` | `yc1_z * s_G/s_Xy` |
| `yc2` | `yc2_z * s_G/s_Xl` |
| `yc3` | `yc3_z * s_F/s_Xy` |
| `yc4`, `yc5` | `yc_z * s_E/s_Xy` |
| `yc6` | `yc6_z * s_E/s_Xa` |
| `yc7` | `yc7_z * s_L/s_Xl` |
| `yc8` | `yc8_z * s_L/s_Xa` |
| `yc9` | `yc9_z * s_A/s_Xl` |
| `yc10`, `yc11` | `yc_z * s_A/s_Xa` |

These algebraic deductions from the archived differential equations agree with
the authors' conversion formulas in supplement S3, pages 8–9.
`table.R` applies the same conversions except that its `yc3` lines use `s_G`
where the fructose equation requires `s_F`. Preserve the literal script for
reproducibility, and keep this discrepancy visible when comparing parameter
tables. **S3 page 9 explicitly uses the fructose maximum for this coefficient**,
confirming that the table script disagrees with the published derivation. The
Julia conversion should follow S3; comparisons to literal `table.R` results
must label the difference.

## Priors and observation model

The sampler works in normalized state coordinates. All 24 kinetic parameters
have independent positive-truncated `Normal(0.5,0.3)` priors. One common
`sigma_z > 0` has a half-Cauchy prior with scale 1. These values come directly
from `model.stan:117–144`, also described by paper equation (2.13). Support is
positive and unbounded above. The source's absolute-Normal random initialization
is distinct from these priors.

At each fitted time and state:

```text
data_z[t,i] ~ Normal(solution_z(t)[i], sigma_z)
```

The mean is the **integrated state**, as implemented by `x_hat` in Stan, not the
RHS derivative. This distinction matters because the prose around paper
equations (2.10)–(2.12) reuses an `f` symbol for model prediction. `sigma_z` is
estimated jointly with the 24 kinetic coefficients. The normalized negative log
likelihood includes `n*log(sigma_z)` with `n=T*8`; dropping this term would change
the inferred noise scale. In physical coordinates each state's observation SD
is `s_i*sigma_z`, rather than one shared SD in mg/g pulp.

Fit each of the three datasets separately. The initial row sets the ODE initial
condition and contributes no likelihood term. Original solver and sampler
settings are listed in [P2.1](cocoa_original.md#what-the-original-programs-do).

## Implementation and remaining source checks

[P2.3](cocoa_baseline.md) implements all five additive growth routes, three
mortality routes, source state ordering, and physical/normalized parameter
conversion. Tests cover the eight RHS equations, analytical limits, positivity,
coordinate equivalence, and solver derivatives. Supplement Tables S4–S6
(pages 11–13) have been transcribed with provenance into three parameter files.
[P2.4](cocoa_fitting.md) adds the inferred noise scale, exact source priors, and
fit reporting without altering the source data or equations.

These marginal summaries do not recover joint parameter correlations. A
trajectory evaluated at the reported parameter means is a plug-in prediction;
for this nonlinear model it need not equal the posterior mean trajectory.
P2.3 comparisons must identify which prediction is being reproduced, and
marginal quantiles must not be treated as joint posterior draws when building
uncertainty bands.
