"""Audit complete or softly stopped adaptation runs without changing their artifacts."""
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import tomllib

root = Path.cwd()
history = root / "data/processed/cocoa_benchmark/fits"
directory = history / "bayes_brwb3_adaptation_pilot1"
checked_hashes = 0


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def toml(path):
    return tomllib.loads(path.read_text(encoding="utf-8"))


def rows(path):
    with path.open(encoding="utf-8", newline="") as stream:
        result = list(csv.DictReader(stream))
    assert all(None not in row and all(value is not None for value in row.values()) for row in result), path
    return result


def check_hash(path, expected):
    global checked_hashes
    assert digest(path) == expected, path
    checked_hashes += 1


def boolean(value):
    assert value in ("true", "false"), value
    return value == "true"


def near(a, b):
    if not (math.isfinite(a) and math.isfinite(b)):
        assert a == b or (math.isnan(a) and math.isnan(b)), (a, b)
        return
    assert math.isclose(a, b, abs_tol=1e-9, rel_tol=1e-10), (a, b)


def quantile(values, probability):
    ordered = sorted(values)
    index = (len(ordered) - 1) * probability
    lo, hi = math.floor(index), math.ceil(index)
    return ordered[lo] + (index - lo) * (ordered[hi] - ordered[lo])


experiment = toml(directory / "experiment.toml")
assert experiment["status"] in ("complete", "sampling_stopped"), experiment["status"]
complete = experiment["status"] == "complete"
receipt_path = history / "bayes_brwb3_adaptation_pilot1.process.json"
receipt = json.loads(receipt_path.read_text(encoding="utf-8-sig"))
assert receipt["status"] == "exited" and receipt["exit_code"] == 0, receipt
assert receipt["watchdog_seconds"] == experiment["process_watchdog_seconds"] == 2700
for stream in ("stdout", "stderr"):
    check_hash(history / receipt[stream], receipt[stream + "_sha256"])
check_hash(directory / "experiment.jl", experiment["script_sha256"])
check_hash(directory / "watchdog.ps1", experiment["watchdog_script_sha256"])
for file, expected in experiment["artifact_sha256"].items():
    check_hash(directory / file, expected)
for group in ("implementation_sha256", "environment_sha256"):
    for file, expected in experiment[group].items():
        check_hash(root / file, expected)
for prefix in ("reference", "numerical_comparison"):
    parent = root / experiment[prefix + "_directory"]
    check_hash(parent / "experiment.toml", experiment[prefix + "_experiment_sha256"])
    for file, expected in experiment[prefix + "_artifact_sha256"].items():
        check_hash(parent / file, expected)

fit = directory / "fit"
settings = toml(fit / "settings.toml")
trial_path = fit / "brwb3"
report = toml(trial_path / "report.toml")
reference_fit = root / experiment["reference_directory"] / "fit"
reference_input = toml(reference_fit / "brwb3/input.toml")
source_input = toml(trial_path / "input.toml")
assert source_input == reference_input
assert toml(trial_path / "priors.toml") == toml(reference_fit / "brwb3/priors.toml")
check_hash(root / source_input["source_file"], source_input["source_sha256"])
for key in ("samples_per_chain", "warmup_per_chain", "chains", "max_depth", "target_accept",
            "abstol", "reltol", "solver_maxiters_scope", "solver_max_step_attempts",
            "max_sampling_seconds", "stop_on_divergence", "chain_execution", "julia_threads"):
    assert settings[key] == experiment[key], key
assert settings["solver"] == "tsit5" and settings["runtime_schema_version"] == 2
assert settings["seed"] == experiment["base_seed"] == 141085
assert settings["predictive_draws"] == 50
assert settings["trial_seeds"]["brwb3"] == experiment["trial_seed"] == 141100
assert settings["trial_solver_maxiters"]["brwb3"] == experiment["solver_maxiters"] == 10000
warmup, samples, chains = settings["warmup_per_chain"], settings["samples_per_chain"], settings["chains"]
assert (warmup, samples, chains) == (250, 100, 4)
checkpoints = {}
retained = []
for chain in range(1, chains + 1):
    assert (trial_path / f"chain_{chain:03}.csv").read_text(encoding="utf-8").splitlines()[0] == (reference_fit / f"brwb3/chain_{chain:03}.csv").read_text(encoding="utf-8").splitlines()[0]
    saved = rows(trial_path / f"chain_{chain:03}.csv")
    retained.append(len(saved))
    assert len(saved) <= samples
    for draw, row in enumerate(saved, 1):
        assert (int(row["chain"]), int(row["draw"])) == (chain, draw)
        assert all(math.isfinite(float(value)) and float(value) > 0 for key, value in row.items()
                   if key not in ("chain", "draw", "divergent", "tree_depth"))
        boolean(row["divergent"])
        assert 0 <= int(row["tree_depth"]) <= settings["max_depth"]
        checkpoints[chain, draw] = row

runtime = rows(trial_path / "runtime.csv")
assert (trial_path / "runtime.csv").read_text(encoding="utf-8").splitlines()[0] == (reference_fit / "brwb3/runtime.csv").read_text(encoding="utf-8").splitlines()[0]
stats = ("divergent", "tree_depth", "n_steps", "step_size", "acceptance_rate", "hamiltonian_energy",
         "hamiltonian_energy_error", "max_hamiltonian_energy_error")
calls, last_times = [0] * chains, [0.] * chains
global_elapsed = 0.
kept_events = {}
for row in runtime:
    chain, iteration, draw = (int(row[key]) for key in ("chain", "iteration", "retained_draw"))
    assert 1 <= chain <= chains and iteration == calls[chain - 1] + 1
    assert iteration <= warmup + samples
    phase = "initialization" if iteration == 1 else "warmup" if iteration <= warmup else "sampling"
    assert row["phase"] == phase and draw == max(0, iteration - warmup)
    elapsed, gap = float(row["elapsed_seconds"]), float(row["since_previous_completion_seconds"])
    assert math.isfinite(elapsed) and math.isfinite(gap) and gap >= 0 and elapsed >= global_elapsed
    near(elapsed - last_times[chain - 1], gap)
    if iteration == 1:
        assert all(row[key] == "" for key in stats)
    else:
        divergent = boolean(row["divergent"])
        assert 0 <= int(row["tree_depth"]) <= settings["max_depth"] and int(row["n_steps"]) > 0
        step_size, acceptance = float(row["step_size"]), float(row["acceptance_rate"])
        assert math.isfinite(step_size) and step_size > 0
        assert 0 <= acceptance <= 1 or (divergent and not math.isfinite(acceptance))
        assert all(math.isfinite(float(row[key])) for key in stats[5:]) or divergent
    if draw:
        saved = checkpoints[chain, draw]
        assert row["divergent"] == saved["divergent"] and row["tree_depth"] == saved["tree_depth"]
        kept_events[chain, draw] = row
    calls[chain - 1], last_times[chain - 1], global_elapsed = iteration, elapsed, elapsed
assert retained == [max(0, count - warmup) for count in calls]
assert len(kept_events) == len(checkpoints)

if complete:
    assert calls == [warmup + samples] * chains and retained == [samples] * chains
    assert report["status"] in ("posterior_diagnostics_passed", "posterior_diagnostics_failed")
    assert experiment["fit_accepted"] == report["diagnostics_passed"]
    assert experiment["joint_diagnostics_available"]
    raw_statistics = rows(trial_path / "sampler_statistics.csv")
    assert len(raw_statistics) == len(checkpoints)
    for index, row in enumerate(raw_statistics):
        chain, draw = index // samples + 1, index % samples + 1
        assert (int(row["chain"]), int(row["draw"])) == (chain, draw)
        event = kept_events[chain, draw]
        for key in stats:
            if key == "divergent":
                assert boolean(row[key]) == boolean(event[key])
            else:
                left, right = float(row[key]), float(event[key])
                assert left == right or (math.isnan(left) and math.isnan(right))
    for artifact in (directory / "fit.audit.toml", directory / "runtime_diagnostics/settings.toml"):
        audited = toml(artifact)
        for file, expected in audited["parent_sha256"].items():
            check_hash(fit / file, expected)
    mixing = rows(directory / "mixing_comparison.csv")
    assert [row["case"] for row in mixing] == ["reference_25", "longer_full_100", "longer_first_25", "longer_last_25"]
    for row, parent_report in ((mixing[0], toml(reference_fit / "brwb3/report.toml")), (mixing[1], report)):
        for key in ("max_rhat", "min_ess_bulk", "min_ess_tail", "divergences", "depth_hits"):
            near(float(row[key]), float(parent_report[key]))
    for row in mixing:
        total = int(row["total_draws"])
        assert total == 4 * int(row["draws_per_chain"])
        near(float(row["depth_hit_fraction"]), int(row["depth_hits"]) / total)
        near(float(row["min_bulk_ess_per_draw"]), float(row["min_ess_bulk"]) / total)
        near(float(row["min_tail_ess_per_draw"]), float(row["min_ess_tail"]) / total)
else:
    assert report["status"] == "sampling_stopped" and not report["diagnostics_passed"]
    assert not experiment["fit_accepted"] and not experiment["joint_diagnostics_available"]
    assert calls == report["completed_sampler_calls"] == experiment["completed_sampler_calls"]
    assert retained == report["partial_retained_draws"] == experiment["partial_retained_draws"]
    for actual, expected in zip(last_times, report["last_completion_seconds"]):
        near(actual, expected)
    assert report["stop_reason"] == experiment["stop_reason"]
    if "divergence" in report["stop_reason"]:
        assert any(boolean(row["divergent"]) for row in checkpoints.values())
    else:
        assert "time budget" in report["stop_reason"] and global_elapsed >= settings["max_sampling_seconds"]
    for file in ("sampler_statistics.csv", "parameters.csv", "predictive.csv", "log_likelihood.csv"):
        assert not (trial_path / file).exists()
    assert not any(key in report for key in ("max_rhat", "min_ess_bulk", "min_ess_tail"))
    mixing = []

phase_summaries = []
for chain in range(1, chains + 1):
    for phase in ("initialization", "warmup", "sampling"):
        selected = [row for row in runtime if int(row["chain"]) == chain and row["phase"] == phase]
        if not selected:
            continue
        duration = [float(row["since_previous_completion_seconds"]) for row in selected]
        transitions = [row for row in selected if row["n_steps"]]
        steps = sum(int(row["n_steps"]) for row in transitions)
        phase_summaries.append({
            "chain": chain, "phase": phase, "calls": len(selected), "median_seconds": statistics.median(duration),
            "p95_seconds": quantile(duration, .95), "max_seconds": max(duration), "leapfrog_steps": steps,
            "seconds_per_leapfrog": sum(float(row["since_previous_completion_seconds"]) for row in transitions) / steps if steps else None,
            "median_leapfrog_steps": statistics.median(int(row["n_steps"]) for row in transitions) if transitions else None,
            "min_step_size": min(float(row["step_size"]) for row in transitions) if transitions else None,
            "max_step_size": max(float(row["step_size"]) for row in transitions) if transitions else None,
            "divergences": sum(boolean(row["divergent"]) for row in transitions),
            "depth_hits": sum(int(row["tree_depth"]) >= settings["max_depth"] for row in transitions),
        })

result = {
    "integrity_passed": True, "diagnostics_accept_fit": False, "outcome": experiment["status"],
    "hash_checks": checked_hashes, "runtime_events_checked": len(runtime), "retained_draws_checked": len(checkpoints),
    "completed_sampler_calls": calls, "retained_per_chain": retained,
    "retained_divergences": sum(boolean(row["divergent"]) for row in checkpoints.values()),
    "retained_depth_hits": sum(int(row["tree_depth"]) >= settings["max_depth"] for row in checkpoints.values()),
    "joint_diagnostics_available": complete, "raw_sampler_correspondence_checked": complete,
    "fit_sampling_accepted": experiment["fit_accepted"], "stop_reason": experiment.get("stop_reason"),
    "phase_summaries": phase_summaries, "mixing_comparison": mixing,
    "experiment_sha256": digest(directory / "experiment.toml"), "process_receipt_sha256": digest(receipt_path),
    "audit_script_sha256": digest(Path(__file__)),
    "limitations": "Partial records establish completion order and internal checkpoint/runtime correspondence only; raw transitions are unavailable after a stop. No joint R-hat/ESS or posterior claim for partial chains. Completion gaps include scheduling and GC; retained states do not reconstruct Hamiltonian paths.",
}
destination = history / "bayes_brwb3_adaptation_pilot1.audit.json"
with destination.open("x", encoding="utf-8") as stream:
    json.dump(result, stream, indent=2)
print(json.dumps(result, indent=2))
