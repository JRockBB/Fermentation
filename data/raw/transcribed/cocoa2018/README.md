# Published posterior summaries

Manually transcribed by Codex on 2026-10-19 from rendered Tables S4, S5, and S6
(pages 11, 12, and 13) of the publisher's supplement, version 2:
doi:10.6084/m9.figshare.7152614.v2, supporting Moreno-Zambrano et al. (2018),
doi:10.1098/rsos.180964. License: CC BY 4.0; attribution to the four article authors.

The preserved PDF is `data/raw/cocoa/rsos180964_si_001.pdf`, SHA-256
`10a337b696e50ca417e915ac90e22064e04a65533b89fb3300fc5b7542084a85`.
The CSV records the printed mean, posterior standard deviation (not SE mean),
and 2.5%, 50%, and 97.5% marginal quantiles. All values retain the printed
three-decimal precision. A printed zero quantile reflects rounding; it does not
change the strictly positive prior support. No axis calibration applies to
these numerical tables. No points were read from a trajectory figure.

Parameter names map top-to-bottom to the archived Stan names: mu1-mu5,
ks1-ks5, k1-k3, yc1-yc11, sigma. `docs/cocoa_equations.md` gives the biological
symbols and mapping. Camu, box 1, and box 2 correspond to Tables S4-S6.

These are **normalized-coordinate** posterior summaries, using each trial's
own state maxima. mu1-mu5 and k1-k3 have units h^-1 in these coordinates;
ks1-ks5, yc1-yc11, and sigma are dimensionless. They are neither physical
mg/g-pulp parameter values nor joint posterior draws. Correlations cannot be
reconstructed from this table.

`scripts/prepare_cocoa.jl` creates the sourced parameter TOMLs from the mean
column and the processed observation CSVs from the untouched R data blocks.
