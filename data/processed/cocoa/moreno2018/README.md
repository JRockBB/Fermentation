# Cocoa observations and source-curve comparison

`scripts/prepare_cocoa.jl` generates the three tidy observation tables from the
untouched R numeric literals in `data/raw/cocoa/CocoaModel-master/simulation/`.
`provenance.toml` records input checksums, the script, units, initial states, and
normalizing maxima. No values are rounded, cleaned, interpolated, or converted
to volumetric units. Initial rows are marked `initial` and must not enter the
likelihood. Subsequent rows are marked `observation` (336 values total).

All values use mg/g pulp, including dry-biomass equivalents for microbes. The
source CFU conversions and their limitations are in `docs/cocoa_original.md`.
No temperature, pH, or measurement-standard-deviation series is supplied.

`scripts/compare_cocoa_figures.jl` generates `figure_comparison.csv` and its
TOML summary from the published parameter means and the digitized central
curves under `data/raw/digitized/moreno2018-figures/`. It reports the existing
reproduction gap against the preselected 2% axis-range tolerance; it does not
change parameter values to improve agreement. The raw digitization README
records attribution, calibration, source PDF checksum, and extraction method.

The `fits/` directory contains saved fit configurations and diagnostic reports.
Its README distinguishes completed, interrupted, and failed runs. Only fits
whose reports pass convergence and reproduction checks may be called accepted
reproductions. Predictive bands use retained joint posterior draws, never
independently sampled marginal table summaries.
