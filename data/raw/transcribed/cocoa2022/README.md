# P2.5 source transcription

Source: Moreno-Zambrano, Ullrich and Hütt (2022),
doi:10.1098/rsos.210274. Supplement version 1:
doi:10.6084/m9.figshare.19149492.v1 (CC BY 4.0).
The exact PDF is preserved under `data/raw/cocoa_benchmark/`.

Codex transcribed the MI(0) rows of Tables S3 and S4 (page 14), S5 and S6
(page 15), and S7 (page 16) on 2026-10-20. These are printed numbers,
not digitized plot coordinates; no axis calibration applies. Blank cells
represent printed dashes, not zeros. Success means a numeric entry is
reported for that baseline fit, not a newly assessed convergence result.
BMA weights are rounded to three decimals and depend on the other fitted
models. A printed zero is a rounded weight, not proof of zero probability.

The first transcription is immutable. A second visual check at original
render resolution found three transcription mistakes. They are recorded in
`baseline_diagnostics_corrections.csv`; the preparation script checks the
old values and applies these corrections only to processed output. The
source PDF itself has not been corrected or changed.

The paper's main text describes bulk/tail ESS >100, while supplement
Section 4 describes 400. The printed successful baseline entries for
`brwb4` (tail ESS 288) and `ecpb1` (bulk ESS 253) are preserved; the
loader does not silently reclassify the paper's results.

No PSIS-LOO numbers were transcribed: Figure S1 contains plotted results,
and the archived plotting script computes them from absent `.Rsave`
outputs. MLE negative log likelihood and training RMSE must not be
presented as reproductions of PSIS-LOO.
