# Original data

Store downloaded bytes under `<food>/<source-or-file>` and keep the accompanying
`.provenance.toml` recording the original URL, UTC acquisition time, SHA-256, size,
source, and license. Never edit or overwrite acquired files. Use
`julia --project scripts/fetch/cocoa.jl` for the first cocoa archive. Rerunning it
verifies the existing bytes and metadata without changing the acquisition time.

Archive extraction must preserve original contents in a new directory; do not
execute downloaded code automatically. Cleaning and unit conversion belong in
scripts that write to `data/processed/`.

Digitized observations belong in `digitized/<source-id>/` with a README identifying
the paper, figure, axis calibration, digitization method, and person or agent who
created the points. Keep the original point exports unchanged.
