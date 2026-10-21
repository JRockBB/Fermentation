# Published cocoa prediction curves

Digitized by Codex on 2026-10-19 and visually checked against all 24 panels.
Source: Moreno-Zambrano, Grimbs, Ullrich and Huett (2018),
*A mathematical model of cocoa bean fermentation*, doi:10.1098/rsos.180964,
Figures 3-5. The article and figures are CC BY 4.0; attribution to these authors.

The accessible figure copies are in Mauricio Moreno-Zambrano's thesis,
*Mathematical Modeling of Cocoa Bean Fermentation*, chapter 2, Figures 2.3-2.5,
printed pages 27-29 (PDF pages 45-47): <https://d-nb.info/1240902832/34>.
The thesis PDF was downloaded on 2026-10-19 and has SHA-256
`d5e63c6eca4c2a82fea9e13f3bb79d456681dfea8a567637632e5613370a6a99`.
The full thesis is not redistributed here; the three PNGs contain only the
article's figure panels and axes. Their crop origins in the full-page renders
are camu `(1560,7306)`, box1 `(1084,6047)`, and box2 `(1520,4706)` pixels.

## Calibration and extraction

The built-in Windows PDF renderer produced 9000 x 12729 pixel pages. Coordinates
in the CSVs refer to these full pages, with origin at the top left, not the PNG
crops. Both axes are linear. `calibration.csv` records the pixel centers of
the 0/144-hour ticks and the zero/highest-labeled concentration ticks for every
panel. The highest labeled tick need not equal the top of the plot, and zero
need not equal its bottom. Units are hours and mg/g pulp, including biomass.

`scripts/digitize_cocoa_figures.ps1` identifies the **solid central red curve**
at 12, 24, 36, 48, 72, 96, and 120 hours. It uses connected red components in a
61-pixel horizontal window to distinguish the solid curve from dashed credible
limits. A component must span at least 55 columns and intersect the center;
ambiguous or occluded points are omitted. All 168 requested points were found.
Blue-circle overlays were visually inspected against the source in all panels.
Experimental black points and dashed interval limits are not digitized.

At this resolution a line is several pixels wide. Values have floating-point
precision for reproducible coordinate conversion, not that measurement accuracy.
Allow approximately 0.5% of the calibrated vertical range for digitization.
The comparison tolerance of 2% of that range was chosen **before** computing
model differences. It must not be enlarged to hide discrepancies.

The source central curve is piecewise linear between displayed simulation
times. These seven comparison times coincide with the archived observation
times for each trial. The archive does not identify the exact posterior summary
used for that central curve or provide its posterior draws.

Reproduction requires the pinned thesis PDF and Windows. Render PDF pages
44-46 (zero-based) with `scripts/render_cocoa_source.ps1`, then run
`scripts/digitize_cocoa_figures.ps1 -ImageDirectory <renders> -OutputDirectory <new-directory>`.
The digitizer refuses to overwrite existing CSVs. Never overwrite these raw
files; perform further transformations under `data/processed/`.
