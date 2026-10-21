"""Independent CSV/TOML integrity check; no posterior refitting or relabelling."""
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import tomllib

ROOT = Path.cwd()
PARENT = ROOT / "data/processed/cocoa_benchmark/fits/bayes_brwb3_budget_pilot2"
OUTPUT = PARENT.with_suffix(".runtime-audit.toml")
SCRIPT = PARENT.with_suffix(".runtime-audit.py")
assert not OUTPUT.exists() and not SCRIPT.exists()


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def toml(path):
    return tomllib.loads(path.read_text(encoding="utf-8"))


def rows(path):
    with path.open(newline="", encoding="utf-8") as file:
        return list(csv.DictReader(file))


experiment = toml(PARENT / "experiment.toml")
assert experiment["all_common_draws_identical"] and not experiment["fit_accepted"]
assert digest(PARENT / "experiment.jl") == experiment["script_sha256"]
parents = {name: digest(PARENT / name) for name in experiment["artifact_sha256"]}
assert parents == experiment["artifact_sha256"]
for name, expected in experiment["environment_sha256"].items():
    assert digest(ROOT / name) == expected

runtime_count = 0
draw_count = 0
for scope in ("solve", "output_interval"):
    case = PARENT / scope
    trial = case / "brwb3"
    settings, report = toml(case / "settings.toml"), toml(trial / "report.toml")
    assert settings["solver_maxiters_scope"] == scope
    assert settings["solver_max_step_attempts"] == 500
    assert report["status"] == "posterior_diagnostics_failed" and not report["diagnostics_passed"]
    assert report["minimum_ess"] == 400 and report["rhat_threshold"] == 1.01
    for name, expected in settings["implementation_sha256"].items():
        assert digest(ROOT / name) == expected
    audit = toml(PARENT / f"{scope}.audit.toml")
    assert audit["integrity_passed"] and audit["solver_budget_scope"] == scope
    assert audit["trials"]["brwb3"]["likelihood_replay_complete"]
    for name, expected in audit["parent_sha256"].items():
        assert digest(case / name) == expected
    runtime = rows(trial / "runtime.csv")
    assert len(runtime) == 4 * (50 + 25)
    elapsed = [float(row["elapsed_seconds"]) for row in runtime]
    assert all(math.isfinite(t) and t >= 0 for t in elapsed)
    assert elapsed == sorted(elapsed)
    divergences = depth_hits = 0
    for chain in range(1, 5):
        checkpoint = rows(trial / f"chain_{chain:03}.csv")
        events = [row for row in runtime if int(row["chain"]) == chain]
        assert len(checkpoint) == 25 and len(events) == 75
        assert [int(row["draw"]) for row in checkpoint] == list(range(1, 26))
        previous = 0.0
        for iteration, event in enumerate(events, 1):
            assert int(event["iteration"]) == iteration
            phase = "initialization" if iteration == 1 else "warmup" if iteration <= 50 else "sampling"
            assert event["phase"] == phase
            assert int(event["retained_draw"]) == max(0, iteration - 50)
            elapsed = float(event["elapsed_seconds"])
            difference = float(event["since_previous_completion_seconds"])
            assert elapsed >= previous and math.isclose(elapsed - previous, difference, abs_tol=1e-9)
            previous = elapsed
            if iteration == 1:
                assert event["divergent"] == event["tree_depth"] == ""
            else:
                assert event["divergent"] in ("true", "false")
                assert 0 <= int(event["tree_depth"]) <= 8
            if iteration > 50:
                kept = checkpoint[iteration - 51]
                assert int(kept["chain"]) == chain
                assert kept["divergent"] == event["divergent"]
                assert kept["tree_depth"] == event["tree_depth"]
                divergences += kept["divergent"] == "true"
                depth_hits += int(kept["tree_depth"]) >= 8
                parameters = [float(value) for name, value in kept.items()
                              if name not in ("chain", "draw", "divergent", "tree_depth")]
                assert len(parameters) == 25 and all(math.isfinite(p) and p > 0 for p in parameters)
        draw_count += len(checkpoint)
        runtime_count += len(events)
    assert divergences == report["divergences"] == 0
    assert depth_hits == report["depth_hits"] == 80

for name in ("parameters.csv", "predictive.csv", "predictive_selection.csv",
             "likelihood_observations.csv", "log_likelihood.csv",
             *(f"chain_{chain:03}.csv" for chain in range(1, 5))):
    assert digest(PARENT / "solve/brwb3" / name) == digest(PARENT / "output_interval/brwb3" / name)
assert all(digest(PARENT / name) == expected for name, expected in parents.items())
shutil.copyfile(__file__, SCRIPT)
summary = {
    "integrity_passed": True,
    "runtime_events_checked": runtime_count,
    "retained_draws_checked": draw_count,
    "paired_draws_identical": draw_count // 2,
    "warmup_excluded_from_draws": True,
    "sampling_diagnostics_passed": False,
    "parent_experiment_sha256": digest(PARENT / "experiment.toml"),
    "audit_script_sha256": digest(SCRIPT),
    "limitation": "Integrity and identical output do not establish MCMC convergence or explain the previous long-run delay.",
}
OUTPUT.write_text("\n".join(f"{key} = {json.dumps(value)}" for key, value in summary.items()) + "\n", encoding="utf-8")
print(json.dumps(summary, indent=2))
