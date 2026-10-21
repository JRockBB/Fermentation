# Original cocoa materials (P2.1)

Inspected 2026-10-19. The archive supports Moreno-Zambrano, Grimbs, Ullrich and
Hütt (2018), *A mathematical model of cocoa bean fermentation*,
[doi:10.1098/rsos.180964](https://doi.org/10.1098/rsos.180964).
The [data record](https://zenodo.org/records/4935403) is
doi:10.5061/dryad.321d33v. The archive README declares CC0-1.0.

## Files and integrity

`data/raw/cocoa/CocoaModel.zip` contains 12,401 bytes. Its SHA-256 is
`9cc6c82ae32c7360408d0bef072e80462df882a004d7aa8975c52f301568c40a`;
its MD5 agrees with the data record. Download provenance remains in
`CocoaModel.zip.provenance.toml`. The files below were extracted byte-for-byte
under `data/raw/cocoa/CocoaModel-master/`; no source code or data were edited.

The ZIP has seven files and two explicit directory entries (`results/` and
`simulation/`). Paths in this table are relative to `CocoaModel-master/`.

| File | Bytes | Contents |
|---|---:|---|
| `INSTRUCTIONS.md` | 3,550 | Dataset identities, execution instructions, references, and scaling notes. |
| `README.txt` | 771 | Purpose, author contact, and CC0 declaration. |
| `simulation/camu.R` | 3,902 | Embedded Camu data, initialization, scaling, and one RStan fit. |
| `simulation/p_box1.R` | 3,718 | Embedded Papalexandratou box 1 data and the same fitting workflow. |
| `simulation/p_box2.R` | 3,713 | Embedded Papalexandratou box 2 data and the same fitting workflow. |
| `simulation/model.stan` | 4,005 | Shared eight-state ODE, priors, unknown observation scale, likelihood, and pointwise log likelihood. |
| `results/table.R` | 9,551 | Load fitted draws, rescale parameters, summarize means/SDs, compare the two boxes, and print a LaTeX table. |

There are no separate CSV datasets, saved posterior draws, compiled binaries,
figures, package lockfiles, or supplementary PDFs. The 23-dataset benchmark from
the 2022 paper is a separate acquisition task (P2.5).

P2.2 additionally acquired `rsos180964_si_001.pdf` separately from the publisher:
doi:10.6084/m9.figshare.7152614.v2, file 13237406, 6,149,683 bytes, CC BY 4.0.
Its provenance sidecar records the download and hash. It is not a ZIP member.
`scripts/fetch/cocoa_supplement.jl` verifies/retrieves this exact version.

[Per-file hashes](cocoa_archive_files.csv) pin the extracted bytes. Run:

```sh
julia --project scripts/verify_cocoa.jl
```

This verifies the archive, provenance, complete member list, byte lengths, and
hashes. It reads only numeric literals from the three R data blocks, checks their
dimensions, and reports scaling maxima. It does not execute R or write data.
Git attributes disable line-ending conversion for raw files.

If only the ZIP is present in a fresh copy, extract it into `data/raw/cocoa/`
without overwriting existing files, then run the verifier. For example, in
PowerShell, `Expand-Archive data/raw/cocoa/CocoaModel.zip data/raw/cocoa` creates
the original directory tree; omit `-Force` to preserve existing files.

## Datasets actually fitted

Every script stores `x0`, `t0`, `ts`, and `x` directly in its data block (lines
21–34). `x0` is fixed as the initial state, and `x` contains subsequent
observations. R's `structure(..., .Dim=c(T,8))` is column-major.

| Script | Original data source | Times including t0 | Likelihood observations |
|---|---|---:|---:|
| `camu.R` | Camu et al. 2007, doi:10.1128/AEM.02189-06 | 17 | 16 × 8 = 128 |
| `p_box1.R` | Papalexandratou et al. 2011, doi:10.1016/j.fm.2011.06.003; box 1 | 14 | 13 × 8 = 104 |
| `p_box2.R` | Same publication; box 2 | 14 | 13 × 8 = 104 |

All arrays use the state order **glucose, fructose, ethanol, lactic acid, acetic
acid, yeast, LAB, AAB**. There are 336 likelihood observations and 24 fixed
initial values. There are no missing-value markers in these embedded arrays.

The Camu grid in hours is
`0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 84, 96, 120, 144`.
Both box grids are
`0, 6, 12, 24, 30, 36, 48, 54, 60, 72, 84, 96, 120, 144`.

The paper identifies the Camu trial as a Ghana heap fermentation and the two
Brazil trials as wooden-box fermentations (§2.1). Time is in hours and source
concentrations are **mg per g of pulp**, including the microbial biomass
converted from CFU (§2.2, §2.4.2, figures 3–5). The conversion assumptions are
15, 1.25, and 0.28 pg/CFU for yeast, LAB, and AAB respectively. These are source
assumptions, not independently measured conversion factors for these datasets.
See the [paper PDF](https://upload.wikimedia.org/wikipedia/commons/5/50/A_mathematical_model_of_cocoa_bean_fermentation.pdf).

No conversion to mmol/L has been attempted: mg/g pulp cannot become a volumetric
concentration without further information. Raw embedded data also do not supply
measurement standard deviations, temperature schedules, pH observations, or
the original log-CFU arrays.

## What the original programs do

The R scripts use `rstan`, clear the R workspace, and require replacing a literal
`setwd("path")`. `table.R` also uses `xtable` and `dplyr`; its helper credits
Margossian and Gillespie. No executable environment is pinned in the archive.
The paper cites RStan 2.14.1 (reference 38), which identifies a historical
reference rather than a reproducible package lockfile. We inspected the source;
we have not run its original RStan sampler.

1. Each dataset is fitted independently. The R driver computes `scl` as the
   maximum of each column, including its initial value (`x1 = rbind(x0,x)`).
2. Stan divides both initial conditions and observations by those maxima
   (`model.stan:46–57`). The likelihood and kinetic parameters therefore operate
   on these dataset-specific dimensionless state coordinates.
3. The RHS (`model.stan:18–37`) has five substrate-dependent growth routes:
   yeast on glucose and fructose, LAB on glucose, and AAB on ethanol and lactate.
   Lactate uptake by AAB has a Contois denominator; the others use Monod
   denominators. Ethanol, lactate, and squared acetate enter the yeast, LAB, and
   AAB mortality terms respectively. There is no explicit temperature or pH
   forcing in this baseline.
4. Stan estimates 24 kinetic coefficients: five `mu`, five `ks`, three `k`, and
   eleven `yc`, plus one shared positive `sigma` per fitted dataset. Initial
   states are data, not estimated parameters.
5. All 24 kinetic priors are `Normal(0.5,0.3)` restricted to positive support;
   `sigma` has a positive `Cauchy(0,1)` prior (`model.stan:58–84,117–144`). These
   are not log-normal priors and are not bounded above at one. Observations have
   independent Gaussian likelihoods on the scaled coordinates.
6. The solver is Stan's `integrate_ode_rk45`, with relative and absolute
   tolerances `1e-6` and maximum step count `1e8` (`model.stan:114`). Each R driver
   requests four chains, 3,000 iterations including 1,000 warmup, seed 100785,
   `adapt_delta=0.985`, initial step size 0.01, and maximum tree depth 12. Its
   initialization draws absolute Normal(0.5,0.2) values; this is an initialization
   rule, distinct from the model priors.
7. Generated quantities provide a `T × 8` pointwise `log_lik` array. The archive
   does not include scripts calculating LOO/WAIC or producing the paper's plots.
8. Each driver saves `output_camu.Rsave`, `output_pbox1.Rsave`, or
   `output_pbox2.Rsave`. These files are absent from the archive. `table.R`
   expects them, rescales posterior draws, computes parameter means and SDs,
   and prints a table including a standardized difference between boxes.

## Reproduction issues and implementation status

- **Unknown observation noise:** P2.4 now supports the original inferred shared
  scale, its half-Cauchy prior, and the full Gaussian normalization. Analytical
  MLE/profile tests and an independent numerical posterior check cover this
  extension. See [the fitting workflow](cocoa_fitting.md).
- **Source-unit basis:** retain mg/g pulp and the documented CFU-to-biomass
  assumptions when transcribing P2.2. Any boundary conversion needs explicit
  dimensional support and provenance.
- **Confirmed rescaling inconsistency:** `table.R:46,77,108` multiplies `yc3`
  by `glu/yst`, while `yc3` multiplies the fructose uptake term in
  `model.stan:27` (via `theta[3]`). Supplement S3 page 9 explicitly specifies
  `fru/yst`, confirming the discrepancy. The [P2.2 transcription](cocoa_equations.md)
  follows S3 and preserves the original script untouched.
- **Missing reference draws:** supplementary Tables S4–S6 have been transcribed
  into sourced normalized parameter files. Joint posterior samples and plotting
  code remain absent. Digitized central curves expose a discrepancy with
  simulations at printed mean parameters; see [P2.3](cocoa_baseline.md).
- **Engine mapping:** P2.3 implements additive uptake routes and product-dependent
  mortality explicitly, reusing the common parameter and inference tools.

P2.1 establishes the archive contents and fitting workflow. Equation transcription
is in P2.2, baseline implementation and parameter extraction in P2.3, and fit
reproduction checks in P2.4. The remaining discrepancies are kept explicit.
