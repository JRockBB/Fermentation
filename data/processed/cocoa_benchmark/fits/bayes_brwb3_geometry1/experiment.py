"""Descriptive log-coordinate correlations of an audited, unaccepted pilot."""
import csv
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
import statistics
import sys
import tomllib

ROOT = Path.cwd()
HISTORY = ROOT / "data/processed/cocoa_benchmark/fits"
PARENT = HISTORY / "bayes_brwb3_ghhp4_v2"
OUTPUT = HISTORY / "bayes_brwb3_geometry1"
AUDIT = HISTORY / "bayes_brwb3_ghhp4_v2.audit.toml"
assert not OUTPUT.exists()

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def toml(path):
    return tomllib.loads(path.read_text(encoding="utf-8"))

def save_csv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

audit = toml(AUDIT)
assert audit["integrity_passed"]
files = ["settings.toml", "brwb3/input.toml", "brwb3/report.toml", "brwb3/priors.toml",
         "brwb3/chain_001.csv", "brwb3/chain_002.csv"]
hashes = {name: digest(PARENT / name) for name in files}
assert all(audit["parent_sha256"][name] == value for name, value in hashes.items())
settings, report = toml(PARENT / "settings.toml"), toml(PARENT / "brwb3/report.toml")
assert report["status"] == "posterior_diagnostics_failed" and not report["diagnostics_passed"]
data = toml(PARENT / "brwb3/input.toml")
assert digest(ROOT / data["source_file"]) == data["source_sha256"]
for entry in toml(PARENT / "brwb3/priors.toml").values():
    assert digest(ROOT / entry["file"]) == entry["sha256"]
chains = {}
names = None
for chain in range(1, settings["chains"] + 1):
    with (PARENT / f"brwb3/chain_{chain:03}.csv").open(newline="") as file:
        reader = csv.DictReader(file)
        columns = reader.fieldnames[2:-2]
        names = columns if names is None else names
        assert columns == names and len(names) == 25
        values = list(reader)
    assert len(values) == settings["samples_per_chain"]
    assert [(int(row["chain"]), int(row["draw"])) for row in values] == [(chain, i + 1) for i in range(len(values))]
    chains[f"chain_{chain}"] = {}
    for name in names:
        x = [float(row[name]) for row in values]
        assert all(math.isfinite(v) and v > 0 for v in x)
        chains[f"chain_{chain}"][name] = list(map(math.log, x))
groups = {**chains, "pooled": {name: list(itertools.chain.from_iterable(chain[name] for chain in chains.values())) for name in names}}
assert len(groups["pooled"][names[0]]) == report["retained_draws"] == 2000
correlations, scales, lookup = [], [], {}
for group, columns in groups.items():
    for name in names:
        x = columns[name]
        assert statistics.stdev(x) > 0
        scales.append({"group": group, "parameter": name, "samples": len(x),
                       "mean_log": statistics.mean(x), "sd_log": statistics.stdev(x),
                       "minimum_log": min(x), "maximum_log": max(x)})
    for a, b in itertools.combinations(names, 2):
        correlation = statistics.correlation(columns[a], columns[b])
        assert -1 <= correlation <= 1
        correlations.append({"group": group, "parameter_a": a, "parameter_b": b, "pearson_log": correlation})
        lookup[group, a, b] = correlation
pairs = []
for a, b in itertools.combinations(names, 2):
    per_chain = [lookup[group, a, b] for group in chains]
    pairs.append({"parameter_a": a, "parameter_b": b, "pooled_pearson_log": lookup["pooled", a, b],
                  "chain_1_pearson_log": per_chain[0], "chain_2_pearson_log": per_chain[1],
                  "minimum_absolute_within_chain": min(map(abs, per_chain)),
                  "same_sign_in_both_chains": per_chain[0] * per_chain[1] > 0})
pairs.sort(key=lambda row: row["minimum_absolute_within_chain"], reverse=True)
assert all(digest(PARENT / name) == value for name, value in hashes.items())
OUTPUT.mkdir()
shutil.copyfile(__file__, OUTPUT / "experiment.py")
save_csv(OUTPUT / "correlations.csv", correlations[0].keys(), correlations)
save_csv(OUTPUT / "parameter_scales.csv", scales[0].keys(), scales)
save_csv(OUTPUT / "ranked_pairs.csv", pairs[0].keys(), pairs)
metadata = {
    "completed_utc": datetime.now(timezone.utc).isoformat(), "trial": "brwb3",
    "parent_directory": str(PARENT.relative_to(ROOT)).replace("\\", "/"),
    "parent_sha256": hashes, "prior_audit_sha256": digest(AUDIT), "source_sha256": data["source_sha256"],
    "script_sha256": digest(OUTPUT / "experiment.py"), "python_version": sys.version,
    "parent_accepted": False, "diagnostics_accept_fit": False,
    "chains": len(chains), "draws": 2000, "parameters": names, "correlations": len(correlations),
    "method": "Pearson correlations and sample SDs of log(parameter) in the fixed archived normalized coordinates, within each chain and pooled. All 300 pairs ranked by the smaller absolute within-chain correlation; signs retained.",
    "units": "Correlations are dimensionless; log coordinates use the archived source-normalized numerical coefficients.",
    "limitations": "Descriptive statistics of an unaccepted pilot, not population posterior estimates, identifiability proof, or a divergence diagnosis. No ODE/model replay, new sampling, parameter-file change, or change to acceptance. Correlation does not establish a cause of slow warmup. Current implementation hash guards for model replay remain unchanged.",
    "output_sha256": {name: digest(OUTPUT / name) for name in ("correlations.csv", "parameter_scales.csv", "ranked_pairs.csv")},
}
(OUTPUT / "settings.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
print(json.dumps(pairs[:8], indent=2))
