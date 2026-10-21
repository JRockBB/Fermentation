# Processed data

Store reproducible outputs under `<food>/<source-id>/`. Each output must identify
its raw input and checksum, the processing script, transformations, and units.
Never repair an input in `data/raw/`; change the processing script and regenerate
the output instead.

Time series should use tidy CSV columns `source_id, experiment, time_h, variable,
value, unit`, with environmental metadata (temperature, salt, population basis)
in additional columns or a source README. Preserve missing observations as
missing, not zero. Convert time to hours and concentration to mmol/L at this
boundary when a documented molecular weight and measurement basis permit it.
