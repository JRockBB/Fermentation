"""Audit the declared solver scan without modifying its parent or result files."""
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import tomllib

root = Path.cwd()
fits = root / "data/processed/cocoa_benchmark/fits"
directory = fits / "bayes_brwb3_runtime_solvers1"
hash_checks = 0


def toml(path):
    return tomllib.loads(path.read_text(encoding="utf-8"))


def csv_rows(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_hash(path, expected):
    global hash_checks
    assert digest(path) == expected, path
    hash_checks += 1


def near(actual, expected):
    assert math.isclose(actual, expected, abs_tol=1e-13, rel_tol=1e-11), (actual, expected)


metadata = toml(directory / "experiment.toml")
assert metadata["status"] == "complete"
assert not metadata["diagnostics_accept_fit"]
parent = root / metadata["parent_directory"]
check_hash(parent / "experiment.toml", metadata["parent_experiment_sha256"])
for file, expected in metadata["parent_artifact_sha256"].items():
    check_hash(parent / file, expected)
for file, expected in metadata["artifact_sha256"].items():
    check_hash(directory / file, expected)
for group in ("implementation_sha256", "environment_sha256"):
    for file, expected in metadata[group].items():
        check_hash(root / file, expected)
check_hash(directory / "experiment.jl", metadata["script_sha256"])

selection = csv_rows(directory / "selection.csv")
assert selection == csv_rows(parent / "runtime_diagnostics/brwb3/selection.csv")
assert len(selection) == 12
cases = [case["name"] for case in metadata["cases"]]
comparisons = csv_rows(directory / "comparison.csv")
timings = csv_rows(directory / "gradient_timings.csv")
assert len(comparisons) == metadata["profile_cases"] == len(selection) * len(cases)
assert len(timings) == metadata["timed_gradients"]
gates = metadata["comparison_gates"]
accuracy = {case: True for case in cases}
ratios = {case: [] for case in cases}
expensive_fast = {case: True for case in cases}
expensive = []
completed = 0

for point in selection:
    chain, draw = int(point["chain"]), int(point["draw"])
    checkpoint = csv_rows(parent / f"fit/brwb3/chain_{chain:03}.csv")[draw - 1]
    reports = {case: toml(directory / f"chain_{chain}_draw_{draw}_{case}.toml") for case in cases}
    fine = [reports[case] for case in ("tsit5_1e10", "rodas5p_1e10")]
    times = {
        case: sorted((row for row in timings if int(row["chain"]) == chain
                      and int(row["draw"]) == draw and row["case"] == case),
                     key=lambda row: int(row["round"])) for case in cases
    }
    for case, report in reports.items():
        assert (report["chain"], report["draw"], report["reasons"]) == (chain, draw, point["reasons"])
        profile, numerics = report["profile"], report["numerics"]
        assert not profile["diagnostics_accept_fit"] and not numerics["diagnostics_accept_fit"]
        assert all(float(checkpoint[key]) == value for key, value in profile["parameters"].items())
        if report["status"] != "complete":
            accuracy[case] = False
            expensive_fast[case] = False
            continue
        completed += 1
        assert profile["parameter_order"] == numerics["parameter_order"]
        assert len(profile["gradient_log"]) == 25
        assert len(profile["prediction_normalized"]) == 96
        for key in ("nll", "gradient_log", "prediction_normalized"):
            assert profile[key] == numerics[key]
        if case == "tsit5_1e8":
            old = toml(parent / f"runtime_diagnostics/brwb3/chain_{chain}_draw_{draw}.toml")
            for key in ("nll", "gradient_log", "prediction_normalized", "accepted_ode_steps",
                        "rejected_ode_steps", "rhs_evaluations", "parameters"):
                assert profile[key] == old[key]
        assert [int(row["round"]) for row in times[case]] == [1, 2, 3, 4, 5]
        for row in times[case]:
            assert all(math.isfinite(float(row[key])) and float(row[key]) >= 0
                       for key in ("seconds", "allocated_bytes", "gc_seconds"))
        for fd in numerics["finite_differences"]:
            near(fd["max_scaled_error"], max(abs(a - b) / max(1., abs(b))
                 for a, b in zip(fd["gradient"], profile["gradient_log"])))
        if any(ref["status"] != "complete" for ref in fine):
            accuracy[case] = False
            expensive_fast[case] = False
            continue
        prediction_error = max(abs(a - b) for ref in fine
            for a, b in zip(profile["prediction_normalized"], ref["profile"]["prediction_normalized"]))
        nll_error = max(abs(profile["nll"] - ref["profile"]["nll"]) for ref in fine)
        gradient_error = max(abs(a - b) / max(1., abs(b)) for ref in fine
            for a, b in zip(profile["gradient_log"], ref["profile"]["gradient_log"]))
        fd_error = max(fd["max_scaled_error"] for fd in numerics["finite_differences"])
        numerical_ok = (prediction_error <= gates["max_normalized_prediction_difference"]
            and nll_error <= gates["max_absolute_nll_difference"]
            and gradient_error <= gates["max_scaled_gradient_difference"]
            and fd_error <= gates["max_fd_scaled_error"])
        accuracy[case] &= numerical_ok
        measured = [float(row["seconds"]) for row in times[case]]
        baseline = [float(row["seconds"]) for row in times["tsit5_1e8"]]
        ratio = statistics.median(measured) / statistics.median(baseline)
        ratios[case].append(ratio)
        faster = sum(a < b for a, b in zip(measured, baseline))
        rows = [row for row in comparisons if int(row["chain"]) == chain
                and int(row["draw"]) == draw and row["case"] == case]
        assert len(rows) == 1
        row = rows[0]
        assert row["numerical_screen_passed"] == str(numerical_ok).lower()
        for key, expected in (
            ("max_prediction_difference_to_fine", prediction_error),
            ("max_absolute_nll_difference_to_fine", nll_error),
            ("max_scaled_gradient_difference_to_fine", gradient_error),
            ("max_fd_scaled_error", fd_error), ("time_ratio_to_tsit5_1e8", ratio),
            ("median_gradient_seconds", statistics.median(measured)), ("faster_rounds", faster)):
            near(float(row[key]), expected)
        for key in ("accepted_ode_steps", "rejected_ode_steps", "rhs_evaluations"):
            assert int(row[key]) == profile[key]
        if chain == 4 and draw in (24, 25):
            expensive_fast[case] &= (ratio <= gates["max_expensive_state_time_ratio"]
                and faster >= gates["minimum_faster_rounds_at_each_expensive_state"])
            expensive.append(row)

reference_ok = accuracy["tsit5_1e10"] and accuracy["rodas5p_1e10"]
candidates = {case: reference_ok and accuracy[case] and expensive_fast[case]
              and len(ratios[case]) == 12 and statistics.median(ratios[case]) <= gates["max_median_paired_time_ratio"]
              for case in ("rodas5p_1e8", "rodas5p_1e10")}
assert metadata["completed_cases"] == completed
assert metadata["numerical_screen_passed"] == accuracy
assert metadata["fine_reference_screen_passed"] == reference_ok
assert metadata["pilot_candidate_screen_passed"] == candidates
for case, values in ratios.items():
    if values:
        near(metadata["median_paired_gradient_time_ratio"][case], statistics.median(values))
summary = {
    "integrity_passed": True, "diagnostics_accept_fit": False, "hash_checks": hash_checks,
    "profile_cases": len(comparisons), "completed_cases": completed, "timed_gradients": len(timings),
    "fine_reference_screen_passed": reference_ok, "numerical_screen_passed": accuracy,
    "pilot_candidate_screen_passed": candidates,
    "median_paired_gradient_time_ratio": metadata["median_paired_gradient_time_ratio"],
    "worst_numerics": {case: {key: max(float(row[key]) for row in comparisons
        if row["case"] == case and row["status"] == "complete") for key in (
        "max_prediction_difference_to_fine", "max_absolute_nll_difference_to_fine",
        "max_scaled_gradient_difference_to_fine", "max_fd_scaled_error")}
        for case in cases if any(row["case"] == case and row["status"] == "complete" for row in comparisons)},
    "expensive_states": expensive,
    "experiment_sha256": digest(directory / "experiment.toml"),
    "audit_script_sha256": digest(Path(__file__)),
}
destination = fits / "bayes_brwb3_runtime_solvers1.audit.json"
with destination.open("x", encoding="utf-8") as stream:
    json.dump(summary, stream, indent=2)
print(json.dumps(summary, indent=2))
