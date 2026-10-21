"""Run from the repository root. Intentionally rejects pilot1's runtime logs."""
import csv
import hashlib
from pathlib import Path
import tomllib

root = Path("data/processed/cocoa_benchmark/fits/bayes_brwb3_budget_pilot1")
experiment = tomllib.loads((root / "experiment.toml").read_text())
issues = []

def rows(relative):
    path = root / relative
    assert hashlib.sha256(path.read_bytes()).hexdigest() == experiment["artifact_sha256"][relative]
    with path.open(newline="") as file:
        return list(csv.DictReader(file))

for scope in ("solve", "output_interval"):
    runtime = rows(f"{scope}/brwb3/runtime.csv")
    for chain in range(1, 5):
        checkpoints = rows(f"{scope}/brwb3/chain_{chain:03}.csv")
        events = [r for r in runtime if int(r["chain"]) == chain and int(r["retained_draw"]) > 0]
        assert len(events) == len(checkpoints) == 25
        for event, kept in zip(events, checkpoints):
            if (event["divergent"], event["tree_depth"]) != (kept["divergent"], kept["tree_depth"]):
                issues.append((scope, chain, int(kept["draw"]), kept["tree_depth"], event["tree_depth"]))
print("Runtime discrepancies (scope, chain, draw, checkpoint depth, runtime depth):", issues)
assert not issues, "Runtime/checkpoint diagnostics differ; this pilot is preserved with a failed runtime check."
