#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shlex
import subprocess
import sys
import time
import signal
from collections import namedtuple
from datetime import date
from pathlib import Path

from trsm_inputs import json_safe
from trsm_parallel import (THREAD_ENVIRONMENT, available_cpus, cancellation_signals,
                           effective_jobs, group_alive, job_limit, signal_group,
                           stop_processes, worker_environment)
from trsm_scan_campaign import (CampaignLock, CampaignStateError, atomic_write_json,
                                checkpoint_path, configuration_fingerprint,
                                load_json, utc_now, _file_digest)

CAMPAIGN_SCHEMA = "trsm_seed_campaign_v1"
OPERATIONAL_OPTIONS = {"jobs", "heartbeat_seconds", "checkpoint_every", "shutdown_grace_seconds",
                       "resume", "aggregate_only", "campaign_dir"}
EXECUTION_OPTIONS = OPERATIONAL_OPTIONS - {"resume", "aggregate_only", "campaign_dir"}
PATH_OPTIONS = {"generator_script", "run_cwd", "python_executable", "ewpt_plot_output",
                "ewpt_executable", "ewpt_minima_executable"}


RUN_TAG_PREFIX = "13.6"
RUN_MG5 = False
STRENGTH_PRIORITY = ("nucl", "perc", "compl", "crit")
STRENGTH_RE = re.compile(
    r"fopt_strength_(?P<kind>[a-z_]+)_(?P<index>\d+):.*?"
    r"ew_jump/T=(?P<value>[-+0-9.eE]+)"
)

Strength = namedtuple(
    "Strength",
    ["value", "kind", "transition_index", "summary_path"],
    defaults=[math.nan, "", "", ""],
)
SeedResult = namedtuple(
    "SeedResult",
    [
        "seed",
        "returncode",
        "viable_count",
        "ewpt_runs",
        "best_ew_jump_over_T",
        "best_strength_kind",
        "point_output",
        "log_path",
        "ewpt_workdir",
        "elapsed_seconds",
        "command",
        "best_point",
        "status", "draw_count", "evo_thc_count", "error",
    ],
    defaults=[0, 0, 0, math.nan, "", "", "", "", 0.0, [], None, "complete", 0, 0, ""],
)
CampaignResult = namedtuple(
    "CampaignResult",
    ["seed_results", "combined_points", "best_points"],
)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Run generate_trsm_points.py for a range of seeds and aggregate viable/EWPT results."
    )
    parser.add_argument("--seed-start", type=int)
    parser.add_argument("--nseeds", type=int)
    parser.add_argument("--nrandom", type=int, default=100)
    parser.add_argument("--nrandom-count-evo-thc", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--checkpoint-every", type=int, default=1)
    parser.add_argument("--jobs", type=job_limit, default=4, help="Concurrent scans, or auto; capped by available CPUs.")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--resume", action="store_true", help="Continue unfinished scans in an existing campaign.")
    mode.add_argument("--aggregate-only", action="store_true", help="Rebuild products from committed rows without running physics.")
    parser.add_argument("--shutdown-grace-seconds", type=float, default=30.0)
    parser.add_argument(
        "--heartbeat-seconds",
        type=float,
        default=30.0,
        help="Print a running-status heartbeat this often while waiting for seed jobs; use 0 to disable.",
    )
    parser.add_argument("--campaign-dir", type=Path, required=True)
    parser.add_argument(
        "--generator-script",
        type=Path,
        default=Path(__file__).resolve().with_name("generate_trsm_points.py"),
        help="Path to generate_trsm_points.py; mainly useful for tests.",
    )
    parser.add_argument(
        "--run-cwd",
        type=Path,
        default=None,
        help="Optional parent for isolated seed directories; defaults to <campaign-dir>/seeds.",
    )
    parser.add_argument(
        "--python-executable",
        type=Path,
        default=default_python_executable(),
        help="Python interpreter used for seed subprocesses; defaults to the active virtualenv Python if available.",
    )
    parser.add_argument(
        "--write-dm-failed",
        action="store_true",
        help="Forward --write-dm-failed to generate_trsm_points.py for each seed.",
    )
    parser.add_argument("--run-ewpt", action="store_true")
    parser.add_argument("--run-ewpt-on-dm-failed", action="store_true")
    parser.add_argument("--ewpt-require-eq418", action="store_true")
    parser.add_argument("--ewpt-thigh", type=float, default=300.0)
    parser.add_argument("--ewpt-plot-phases", action="store_true")
    parser.add_argument("--ewpt-plot-output", type=Path)
    parser.add_argument("--ewpt-plot-format", choices=["png", "pdf", "both"], default="both")
    parser.add_argument("--ewpt-executable", type=Path)
    parser.add_argument("--ewpt-minima-executable", type=Path)
    parser.add_argument("--ewpt-sym-threshold", type=float, default=1.0)
    parser.add_argument("--ewpt-w1-threshold", type=float, default=5.0)
    parser.add_argument("--ewpt-wx-threshold", type=float, default=1.0)
    parser.add_argument("--ewpt-ws-threshold", type=float, default=1.0)
    parser.add_argument("--generator-extra-arg", action="append", default=[], help=argparse.SUPPRESS)
    argv = list(sys.argv[1:] if argv is None else argv)
    args = parser.parse_args(argv)
    args._provided = {action.dest for action in parser._actions
                      if any(token.split("=", 1)[0] in action.option_strings for token in argv)}

    if not (args.resume or args.aggregate_only) and (args.seed_start is None or args.nseeds is None):
        parser.error("fresh campaigns require --seed-start and --nseeds")
    if args.nseeds is not None and args.nseeds <= 0:
        parser.error("--nseeds must be positive")
    if args.nrandom < 0:
        parser.error("--nrandom must be non-negative")
    if args.checkpoint_every < 1:
        parser.error("--checkpoint-every must be positive")
    if not math.isfinite(args.shutdown_grace_seconds) or args.shutdown_grace_seconds < 0:
        parser.error("--shutdown-grace-seconds must be finite and non-negative")
    if not math.isfinite(args.heartbeat_seconds) or args.heartbeat_seconds < 0:
        parser.error("--heartbeat-seconds must be non-negative")
    # Accept old configurations which forwarded these through the escape hatch.
    extra = list(args.generator_extra_arg)
    cleaned = []
    index = 0
    while index < len(extra):
        option, separator, value = extra[index].partition("=")
        if option == "--checkpoint-every":
            if not separator:
                index += 1
                if index == len(extra): parser.error("--checkpoint-every needs a value")
                value = extra[index]
            if "checkpoint_every" not in args._provided:
                args.checkpoint_every = int(value)
        elif option == "--nrandom-count-evo-thc":
            if "nrandom_count_evo_thc" not in args._provided:
                args.nrandom_count_evo_thc = True
        elif option in {"--nrandom", "--resume-from", "--output-manifest", "--preflight",
                        "--ewpt-workdir", "--ewpt-multithreading", "--no-ewpt-multithreading",
                        "--run-mg5", "--ewpt-plot-output"}:
            parser.error(f"{option} must not be supplied through --generator-extra-arg in a campaign")
        else:
            cleaned.append(extra[index])
        index += 1
    args.generator_extra_arg = cleaned
    if args.checkpoint_every < 1:
        parser.error("--checkpoint-every must be positive")
    if args.ewpt_plot_output is not None and (args.ewpt_plot_output.is_absolute() or ".." in args.ewpt_plot_output.parts):
        parser.error("--ewpt-plot-output must be relative to each point directory")
    return args


def default_python_executable(environ=None):
    if environ is None:
        environ = os.environ
    virtual_env = environ.get("VIRTUAL_ENV")
    if virtual_env:
        candidate = Path(virtual_env) / "bin" / "python"
        if candidate.exists():
            return candidate
    return Path(sys.executable)


def expand_seeds(seed_start, nseeds):
    return list(range(seed_start, seed_start + nseeds))


def run_tag(seed, run_date=None):
    if run_date is None:
        run_date = date.today()
    return f"{RUN_TAG_PREFIX}-{run_date.strftime('%Y%m%d')}-{seed}-{RUN_MG5}_vxzero"


def seed_output_path(run_cwd, seed, run_date=None):
    return Path(run_cwd) / "output" / f"trsm_points_{run_tag(seed, run_date)}.dat"


def seed_dm_failed_output_path(run_cwd, seed, run_date=None):
    return Path(run_cwd) / "output" / f"trsm_points_{run_tag(seed, run_date)}_dm_failed.dat"


def build_generator_command(args, seed, generator_script, ewpt_workdir):
    command = [
        str(args.python_executable),
        "-u",
        str(Path(generator_script)),
        str(seed),
        "--nrandom",
        str(args.nrandom),
    ]
    if args.write_dm_failed:
        command.append("--write-dm-failed")
    if args.run_ewpt_on_dm_failed:
        command.append("--run-ewpt-on-dm-failed")
    if args.run_ewpt:
        command.append("--run-ewpt")
    if args.run_ewpt or args.run_ewpt_on_dm_failed:
        command.extend(["--ewpt-workdir", str(ewpt_workdir)])
        if args.ewpt_require_eq418:
            command.append("--ewpt-require-eq418")
        command.extend(["--ewpt-thigh", str(args.ewpt_thigh)])
        if args.ewpt_plot_phases:
            command.append("--ewpt-plot-phases")
        if args.ewpt_plot_output is not None:
            command.extend(["--ewpt-plot-output", str(args.ewpt_plot_output)])
        command.extend(["--ewpt-plot-format", args.ewpt_plot_format])
        if args.ewpt_executable is not None:
            command.extend(["--ewpt-executable", str(args.ewpt_executable)])
        if args.ewpt_minima_executable is not None:
            command.extend(["--ewpt-minima-executable", str(args.ewpt_minima_executable)])
        command.extend(
            [
                "--ewpt-sym-threshold",
                str(args.ewpt_sym_threshold),
                "--ewpt-w1-threshold",
                str(args.ewpt_w1_threshold),
                "--ewpt-wx-threshold",
                str(args.ewpt_wx_threshold),
                "--ewpt-ws-threshold",
                str(args.ewpt_ws_threshold),
            ]
        )
    command.extend(args.generator_extra_arg)
    if args.nrandom_count_evo_thc:
        command.append("--nrandom-count-evo-thc")
    command.extend(["--checkpoint-every", str(args.checkpoint_every), "--no-ewpt-multithreading"])
    return command


def read_viable_rows(path, seed):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="ascii", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        return [{"seed": str(seed), **row} for row in reader]


def best_strength_from_summary(path):
    matches = []
    text = Path(path).read_text(encoding="ascii", errors="replace")
    for match in STRENGTH_RE.finditer(text):
        kind = match.group("kind")
        if kind not in STRENGTH_PRIORITY:
            continue
        matches.append(
            Strength(
                value=float(match.group("value")),
                kind=kind,
                transition_index=match.group("index"),
                summary_path=str(path),
            )
        )
    return select_best_strength(matches)


def strength_sort_key(strength):
    if strength is None or math.isnan(strength.value):
        return (len(STRENGTH_PRIORITY), 0.0)
    return (STRENGTH_PRIORITY.index(strength.kind), -strength.value)


def select_best_strength(strengths):
    strengths = [strength for strength in strengths if strength is not None and not math.isnan(strength.value)]
    if not strengths:
        return Strength()
    return sorted(strengths, key=strength_sort_key)[0]


def read_trsm_input(path):
    path = Path(path)
    if not path.exists():
        return {}
    with path.open(encoding="ascii", newline="") as stream:
        rows = list(csv.reader(stream, delimiter="\t"))
    if len(rows) < 2:
        return {}
    header = ["point_index"] + rows[0][1:]
    return dict(zip(header, rows[1]))


def collect_best_point_records(seed, seed_ewpt_dir, max_draw=None):
    for summary_path in Path(seed_ewpt_dir).glob("point_*/ewpt_summary.txt"):
        index = summary_path.parent.name.removeprefix("point_")
        if not index.isdigit() or (max_draw is not None and int(index) > max_draw):
            continue
        strength = best_strength_from_summary(summary_path)
        if math.isnan(strength.value):
            continue
        point_info = read_trsm_input(summary_path.parent / "TRSM_Input.tsv")
        yield (
            {
                "seed": str(seed),
                "point_index": point_info.get("point_index", summary_path.parent.name),
                "strength_kind": strength.kind,
                "transition_index": strength.transition_index,
                "best_ew_jump_over_T": format_float(strength.value),
                "ewpt_summary": str(summary_path),
                "m1": point_info.get("m1", ""),
                "m2": point_info.get("m2", ""),
                "m3": point_info.get("m3", ""),
                "vs": point_info.get("vs", ""),
                "a12": point_info.get("a12", ""),
                "lx": point_info.get("lx", ""),
                "lphix": point_info.get("lphix", ""),
                "lsx": point_info.get("lsx", ""),
            }
        )


def best_record(records):
    return min(
        records,
        key=lambda record: (
            STRENGTH_PRIORITY.index(record["strength_kind"]),
            -float(record["best_ew_jump_over_T"]),
            str(record["point_index"]),
        ),
        default=None,
    )


def format_float(value):
    if value is None or math.isnan(value):
        return "nan"
    return f"{value:.12g}"


def configuration(args):
    return {key: str(value) if isinstance(value, Path) else value
            for key, value in vars(args).items()
            if not key.startswith("_") and key not in OPERATIONAL_OPTIONS}


def canonicalize_paths(args):
    args.campaign_dir = args.campaign_dir.expanduser().resolve()
    for key in PATH_OPTIONS - {"ewpt_plot_output"}:
        value = getattr(args, key)
        if value is not None:
            # Keep the venv interpreter path: resolving its symlink loses the venv.
            setattr(args, key, Path(os.path.abspath(value.expanduser())))
    # Resolve forwarded file paths before moving into private working directories.
    path_flags = {"--micromegas-main", "--dm-limit-table", "--ewpt-executable", "--ewpt-minima-executable"}
    extra = iter(args.generator_extra_arg)
    normalized = []
    for token in extra:
        option, separator, value = token.partition("=")
        if option in path_flags:
            if not separator:
                try:
                    value = next(extra)
                except StopIteration as error:
                    raise CampaignStateError(f"{option} needs a path") from error
            if option.removeprefix("--").replace("-", "_") in args._provided:
                continue
            normalized.append(option + "=" + str(Path(value).expanduser().resolve()))
        else:
            normalized.append(token)
    args.generator_extra_arg = normalized


def hydrate_configuration(args, state):
    for key, value in state.get("execution", {}).items():
        if key in EXECUTION_OPTIONS and key not in args._provided:
            setattr(args, key, value)
    saved = state["configuration"]
    for key, value in saved.items():
        current = getattr(args, key)
        if isinstance(current, Path):
            current = str(current)
        if key in args._provided and current != value:
            raise CampaignStateError(f"Cannot change saved campaign option {key!r} during resume/aggregation")
        setattr(args, key, Path(value) if key in PATH_OPTIONS and value is not None else value)


def run_preflight(args):
    command = build_generator_command(args, args.seed_start, args.generator_script,
                                      args.campaign_dir / "preflight-ewpt") + ["--preflight"]
    print("Checking scan runtime before launching workers...", flush=True)
    env = worker_environment(args.campaign_dir / "preflight-cache")
    with (args.campaign_dir / "preflight.log").open("a", encoding="utf-8") as log:
        log.write("command: " + shlex.join(command) + "\n")
        with cancellation_signals() as stop:
            process = subprocess.Popen(command, cwd=args.campaign_dir, env=env,
                                       stdout=subprocess.PIPE, stderr=log, text=True,
                                       start_new_session=True)
            deadline = time.monotonic() + 180
            try:
                while True:
                    try:
                        stdout, _ = process.communicate(timeout=.1)
                        break
                    except subprocess.TimeoutExpired:
                        if stop["signal"] is not None:
                            raise CampaignStateError("Runtime preflight interrupted; campaign remains resumable")
                        if time.monotonic() >= deadline:
                            raise CampaignStateError("Runtime preflight timed out; see preflight.log")
            except BaseException:
                stop_processes([process], args.shutdown_grace_seconds)
                raise
        log.write(stdout)
    if process.returncode:
        raise CampaignStateError(f"Runtime preflight failed; see {args.campaign_dir / 'preflight.log'}")
    receipts = [line.removeprefix("TRSM_PREFLIGHT ") for line in stdout.splitlines()
                if line.startswith("TRSM_PREFLIGHT ")]
    if len(receipts) != 1:
        raise CampaignStateError("Generator did not return a unique TRSM_PREFLIGHT receipt")
    receipt = json.loads(receipts[0])
    receipt["generator_sha256"] = hashlib.sha256(args.generator_script.read_bytes()).hexdigest()
    return receipt


def new_seed_record(seed, args):
    directory = args.campaign_dir
    return {"seed": seed, "status": "queued", "attempts": 0,
            "cwd": str((args.run_cwd or directory / "seeds") / f"seed_{seed}"),
            "manifest": str(directory / "manifests" / f"seed_{seed}.json"),
            "ewpt_workdir": str(directory / "ewpt" / f"seed_{seed}"),
            "log_path": str(directory / "logs" / f"seed_{seed}.log"),
            "point_output": "", "returncode": None, "elapsed_seconds": 0.0,
            "draw_count": 0, "evo_thc_count": 0, "error": ""}


def seed_checkpoint(record):
    """Find authoritative state; never guess a run tag or adopt loose TSVs."""
    manifest = Path(record["manifest"])
    output = record.get("point_output")
    if manifest.is_file():
        payload = load_json(manifest, "seed output manifest")
        try:
            output = payload["outputs"]["main"]["path"]
        except (KeyError, TypeError) as error:
            raise CampaignStateError(f"Malformed output manifest: {manifest}") from error
    if not output:
        # Recover a crash between the initial checkpoint and manifest writes.
        matches = list((Path(record["cwd"]) / "output").glob("*.checkpoint.json"))
        if len(matches) > 1:
            raise CampaignStateError(f"Ambiguous checkpoints for seed {record['seed']}")
        if matches:
            output = load_json(matches[0], "seed checkpoint").get("scan_path")
    if not output:
        if any((Path(record["cwd"]) / "output").glob("*.dat")):
            raise CampaignStateError(f"Seed {record['seed']} has output without a recoverable checkpoint")
        return None
    output = Path(output).resolve()
    if not output.is_relative_to(Path(record["cwd"]).resolve()):
        raise CampaignStateError(f"Seed output is outside its isolated directory: {output}")
    checkpoint = load_json(checkpoint_path(output), "seed checkpoint")
    if checkpoint.get("seed") != record["seed"]:
        raise CampaignStateError(f"Checkpoint seed does not match seed {record['seed']}")
    for key in ("draw_count", "evo_thc_count", "viable_count"):
        value = checkpoint.get(key)
        if type(value) is not int or value < 0 or (key != "draw_count" and value > checkpoint["draw_count"]):
            raise CampaignStateError(f"Invalid checkpoint counter {key} for seed {record['seed']}")
    if checkpoint.get("outputs", {}).get("main", {}).get("path") != str(output):
        raise CampaignStateError(f"Checkpoint output path mismatch for seed {record['seed']}")
    record.update(point_output=str(output), draw_count=checkpoint["draw_count"],
                  evo_thc_count=checkpoint["evo_thc_count"])
    return checkpoint


def target_complete(checkpoint, args):
    if checkpoint is None:
        return False
    if checkpoint.get("target") != args.nrandom or checkpoint.get("count_evo_thc") != args.nrandom_count_evo_thc:
        raise CampaignStateError("Seed checkpoint target/counting mode differs from the campaign")
    key = "evo_thc_count" if args.nrandom_count_evo_thc else "draw_count"
    return checkpoint["status"] == "complete" and checkpoint[key] == args.nrandom


def validate_prefix(path, expected):
    """Validate saved boundary hashes without reading or modifying the tail."""
    size = expected.get("size")
    if type(size) is not int or size < 0:
        raise CampaignStateError(f"Invalid committed output size: {path}")
    if not Path(path).is_file():
        if size == 0 and not expected.get("exists"):
            return
        raise CampaignStateError(f"Committed output is missing: {path}")
    if Path(path).stat().st_size < size:
        raise CampaignStateError(f"Output is shorter than its checkpoint: {path}")
    if size:
        sample = min(size, 4096)
        for key, offset in (("head_sha256", 0), ("tail_sha256", size - sample)):
            if _file_digest(path, offset, sample) != expected.get(key):
                raise CampaignStateError(f"Committed output {key} mismatch: {path}")
        with Path(path).open("rb") as stream:
            stream.seek(size - 1)
            if stream.read(1) != b"\n":
                raise CampaignStateError(f"Committed output ends in a partial line: {path}")


def committed_lines(path, expected):
    remaining = expected["size"]
    if not remaining:
        return
    with Path(path).open("rb") as stream:
        while remaining:
            line = stream.readline(remaining)
            if not line:
                raise CampaignStateError(f"Output shrank during aggregation: {path}")
            remaining -= len(line)
            yield line.decode("ascii")


def committed_rows(path, checkpoint):
    expected = checkpoint["outputs"]["main"]
    validate_prefix(path, expected)
    reader = csv.DictReader(committed_lines(path, expected), delimiter="\t")
    count = 0
    for row in reader:
        if None in row or None in row.values():
            raise CampaignStateError(f"Malformed committed TSV row: {path}")
        count += 1
        yield row
    if count != expected["data_rows"]:
        raise CampaignStateError(f"Committed row count differs from checkpoint: {path}")


def assert_no_live_workers(state):
    for record in state["seeds"]:
        pid = record.get("pid")
        if pid and group_alive(pid):
            raise CampaignStateError(f"Seed {record['seed']} still owns live process group {pid}; let it finish before resuming")
        output = record.get("point_output")
        if output:
            # Reuse the generator's conservative same-host lock rules.
            with CampaignLock(output, ["campaign", "check-lock"]):
                pass


def launch_seed(record, args):
    checkpoint = seed_checkpoint(record)
    if checkpoint:
        target_complete(checkpoint, args)  # check counting-mode consistency
        command = [str(args.python_executable), "-u", str(args.generator_script),
                   "--resume-from", record["point_output"], "--nrandom", str(args.nrandom),
                   "--checkpoint-every", str(args.checkpoint_every)]
    else:
        command = build_generator_command(args, record["seed"], args.generator_script,
                                          Path(record["ewpt_workdir"]))
        command += ["--output-manifest", record["manifest"]]
    cwd = Path(record["cwd"])
    cwd.mkdir(parents=True, exist_ok=True)
    Path(record["manifest"]).parent.mkdir(parents=True, exist_ok=True)
    log_path = Path(record["log_path"])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log = log_path.open("a", encoding="utf-8", buffering=1)
    log.write(f"\nattempt: {record['attempts'] + 1} at {utc_now()}\ncommand: {shlex.join(command)}\n\n")
    log.flush()
    try:
        process = subprocess.Popen(command, cwd=cwd, env=worker_environment(cwd / "cache"),
                                   stdout=log, stderr=subprocess.STDOUT, text=True,
                                   start_new_session=True)
    except BaseException:
        log.close()
        raise
    record.update(status="running", pid=process.pid, attempts=record["attempts"] + 1,
                  command=command, error="", started_utc=utc_now(), returncode=None)
    return {"process": process, "log": log, "start": time.monotonic(), "record": record}


def finish_seed(job, args, interrupted=False):
    if job.get("finished"):
        return
    record, process = job["record"], job["process"]
    code = process.wait()
    elapsed = time.monotonic() - job["start"]
    job["log"].write(f"\nreturncode: {code}\nelapsed_seconds: {elapsed:.6f}\n")
    job["log"].close()
    job["finished"] = True
    record.update(returncode=code, elapsed_seconds=record["elapsed_seconds"] + elapsed,
                  finished_utc=utc_now(), status="interrupted" if interrupted else "failed")
    try:
        checkpoint = seed_checkpoint(record)
        if checkpoint:
            validate_prefix(record["point_output"], checkpoint["outputs"]["main"])
        if code == 0 and target_complete(checkpoint, args):
            record["status"] = "complete"
            record.pop("pid", None)
        elif not interrupted:
            record["error"] = f"Worker exited {code} without a complete target checkpoint"
    except (OSError, ValueError, KeyError, CampaignStateError) as error:
        record["error"] = str(error)
        record["status"] = "failed"
    if not group_alive(process.pid):
        record.pop("pid", None)


def persist_state(args, state):
    state["updated_utc"] = utc_now()
    atomic_write_json(args.campaign_dir / "campaign_state.json", state)


def heartbeat(state):
    counts = {name: sum(row["status"] == name for row in state["seeds"])
              for name in ("queued", "running", "complete", "failed", "interrupted")}
    for record in state["seeds"]:
        if record["status"] == "running":
            try:
                seed_checkpoint(record)
            except (OSError, ValueError, KeyError, CampaignStateError):
                pass  # initialization/recovery can be between its atomic writes
    print("[running] " + " ".join(f"{key}={value}" for key, value in counts.items()) +
          f" raw_draws={sum(row['draw_count'] for row in state['seeds'])}"
          f" evo/thc={sum(row['evo_thc_count'] for row in state['seeds'])}", flush=True)


def supervise(args, state):
    pending = iter(row for row in state["seeds"] if row["status"] == "queued")
    active = {}
    exhausted = False
    shutdown = False
    started = time.monotonic()
    next_heartbeat = started + args.heartbeat_seconds
    with cancellation_signals() as stop:
        try:
            while active or not exhausted:
                if stop["signal"] is not None and not shutdown:
                    shutdown = True
                    exhausted = True
                    print("Stopping campaign; allowing workers to checkpoint...", flush=True)
                    deadline = time.monotonic() + args.shutdown_grace_seconds
                    for job in active.values():
                        signal_group(job["process"].pid, signal.SIGTERM)
                        job["deadline"] = deadline
                while not exhausted and not shutdown and stop["signal"] is None and len(active) < state["effective_jobs"]:
                    try:
                        record = next(pending)
                    except StopIteration:
                        exhausted = True
                        break
                    try:
                        print(f"[launch] seed={record['seed']}", flush=True)
                        active[record["seed"]] = launch_seed(record, args)
                    except (OSError, ValueError, KeyError, CampaignStateError) as error:
                        record.update(status="failed", error=str(error), returncode=1)
                    persist_state(args, state)
                for seed, job in list(active.items()):
                    process = job["process"]
                    code = process.poll()
                    # A finished leader can leave native descendants behind.
                    if code is not None and group_alive(process.pid) and "deadline" not in job:
                        signal_group(process.pid, signal.SIGTERM)
                        job["deadline"] = time.monotonic() + args.shutdown_grace_seconds
                    if "deadline" in job and time.monotonic() >= job["deadline"]:
                        signal_group(process.pid, signal.SIGKILL)
                    if code is not None and (not group_alive(process.pid) or time.monotonic() >= job.get("deadline", math.inf)):
                        finish_seed(job, args, interrupted=shutdown)
                        del active[seed]
                        persist_state(args, state)
                        done = sum(row["status"] in ("complete", "failed", "interrupted") for row in state["seeds"])
                        print(f"[{done}/{len(state['seeds'])} done] seed={seed} status={job['record']['status']} "
                              f"draws={job['record']['draw_count']} evo/thc={job['record']['evo_thc_count']} "
                              f"elapsed={time.monotonic() - started:.1f}s", flush=True)
                if args.heartbeat_seconds and time.monotonic() >= next_heartbeat:
                    heartbeat(state)
                    persist_state(args, state)
                    next_heartbeat = time.monotonic() + args.heartbeat_seconds
                if active:
                    time.sleep(0.05)
        finally:
            if active:
                try:
                    stop_processes([job["process"] for job in active.values()], args.shutdown_grace_seconds)
                finally:
                    for job in active.values():
                        if job["process"].poll() is not None:
                            finish_seed(job, args, interrupted=True)
                        else:
                            job["record"].update(status="interrupted", error="Worker could not be stopped; its writer lock remains authoritative")
                            job["log"].close()
                    persist_state(args, state)
    return shutdown


def seed_result(record, checkpoint):
    best = None
    ewpt_runs = 0
    if checkpoint is not None:
        # Ignore uncommitted and quarantined point directories.
        for summary in Path(record["ewpt_workdir"]).glob("point_*/ewpt_summary.txt"):
            index = summary.parent.name.removeprefix("point_")
            if index.isdigit() and int(index) <= checkpoint["draw_count"]:
                ewpt_runs += 1
        best = best_record(collect_best_point_records(record["seed"], record["ewpt_workdir"], checkpoint["draw_count"]))
    return SeedResult(seed=record["seed"], returncode=record.get("returncode"),
                      viable_count=checkpoint["viable_count"] if checkpoint else 0,
                      ewpt_runs=ewpt_runs,
                      best_ew_jump_over_T=float(best["best_ew_jump_over_T"]) if best else math.nan,
                      best_strength_kind=best["strength_kind"] if best else "",
                      point_output=record["point_output"], log_path=record["log_path"],
                      ewpt_workdir=record["ewpt_workdir"], elapsed_seconds=record["elapsed_seconds"],
                      command=record.get("command", []), best_point=best, status=record["status"],
                      draw_count=checkpoint["draw_count"] if checkpoint else 0,
                      evo_thc_count=checkpoint["evo_thc_count"] if checkpoint else 0,
                      error=record.get("error", ""))


def seed_result_row(result):
    row = dict(result._asdict())
    row.pop("command")
    row.pop("best_point")
    row["best_ew_jump_over_T"] = format_float(result.best_ew_jump_over_T)
    return row


def write_tsv(path, rows, fieldnames):
    path = Path(path)
    temp = path.with_name("." + path.name + ".tmp")
    try:
        with temp.open("w", encoding="ascii", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({field: row.get(field, "") for field in fieldnames})
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


class TSVRows:
    """Re-iterable aggregate product without retaining the campaign in RAM."""
    def __init__(self, path, count):
        self.path, self.count = path, count

    def __iter__(self):
        with self.path.open(encoding="ascii", newline="") as stream:
            yield from csv.DictReader(stream, delimiter="\t")

    def __len__(self):
        return self.count


def aggregate(args, state):
    directory = args.campaign_dir
    snapshots = []
    results = []
    for record in sorted(state["seeds"], key=lambda row: row["seed"]):
        try:
            checkpoint = seed_checkpoint(record)
            if checkpoint:
                # Fully validate this seed before publishing any of its rows.
                for _ in committed_rows(record["point_output"], checkpoint):
                    pass
                if target_complete(checkpoint, args) and record["status"] != "running":
                    record.update(status="complete", returncode=0, error="")
                elif record["status"] == "complete":
                    raise CampaignStateError("Completed seed no longer has a complete checkpoint")
                snapshots.append((record, checkpoint))
            elif record["status"] == "complete":
                raise CampaignStateError("Completed seed has no checkpoint")
            if record["status"] == "running":
                record["status"] = "interrupted"
        except (OSError, ValueError, KeyError, CampaignStateError) as error:
            checkpoint = None
            record.update(status="failed", returncode=1, error=str(error))
        results.append(seed_result(record, checkpoint))
    names = ("ewpt_baryo_candidate", "ewpt_gw_candidate", "dm_subset", "vacuum_tree_global",
             "rg_bfb", "rg_unitarity", "experimental_subset", "flavour")
    counts = {key: {"True": 0, "False": 0, "unassessed": 0} for key in names}
    combined_path = directory / "combined_points.tsv"
    candidate_path = directory / "candidate_points.tsv"
    tmp_combined = combined_path.with_name(".combined_points.tsv.tmp")
    tmp_candidates = candidate_path.with_name(".candidate_points.tsv.tmp")
    row_count = 0
    candidate_count = 0
    fields = None
    try:
        with tmp_combined.open("w", encoding="ascii", newline="") as combined, tmp_candidates.open("w", encoding="ascii", newline="") as candidates:
            for record, checkpoint in snapshots:
                reader = csv.DictReader(committed_lines(record["point_output"], checkpoint["outputs"]["main"]), delimiter="\t")
                if not reader.fieldnames:
                    continue
                current_fields = ["seed"] + reader.fieldnames
                if fields is None:
                    fields = current_fields
                    writers = [csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n") for stream in (combined, candidates)]
                    for writer in writers:
                        writer.writeheader()
                elif current_fields != fields:
                    raise CampaignStateError("Cannot combine scan outputs with different column layouts")
                for row in reader:
                    row = {"seed": str(record["seed"]), **row}
                    writers[0].writerow(row)
                    row_count += 1
                    for name in names:
                        value = row.get(name)
                        counts[name][value if value in ("True", "False") else "unassessed"] += 1
                    if any(row.get(name) == "True" for name in names[:2]):
                        writers[1].writerow(row)
                        candidate_count += 1
            if fields is None:
                combined.write("seed\n")
                candidates.write("seed\n")
        os.replace(tmp_combined, combined_path)
        os.replace(tmp_candidates, candidate_path)
    finally:
        tmp_combined.unlink(missing_ok=True)
        tmp_candidates.unlink(missing_ok=True)
    best_rows = sorted((result.best_point for result in results if result.best_point),
                       key=lambda row: (STRENGTH_PRIORITY.index(row["strength_kind"]),
                                        -float(row["best_ew_jump_over_T"]), int(row["seed"])))
    best_fields = ["seed", "point_index", "strength_kind", "transition_index", "best_ew_jump_over_T",
                   "ewpt_summary", "m1", "m2", "m3", "vs", "a12", "lx", "lphix", "lsx"]
    write_tsv(directory / "best_points.tsv", best_rows, best_fields)
    summaries = [seed_result_row(result) for result in results]
    write_tsv(directory / "campaign_summary.tsv", summaries, list(summaries[0]))
    atomic_write_json(directory / "candidate_counts.json", counts)
    complete = all(result.status == "complete" for result in results)
    state["status"] = "complete" if complete else "incomplete"
    persist_state(args, state)
    payload = {"metadata": {**state["configuration"], "constraint_version": "trsm_constraints_v2",
                            "jobs": args.jobs, "effective_jobs": state.get("effective_jobs", 0),
                            "complete": complete, "committed_rows_only": True,
                            "best_points_semantics": "Legacy EW-jump ranking; use candidate flags for v2 selection"},
               "totals": {"draw_count": sum(result.draw_count for result in results),
                          "evo_thc_count": sum(result.evo_thc_count for result in results),
                          "viable_count": sum(result.viable_count for result in results),
                          "ewpt_runs": sum(result.ewpt_runs for result in results),
                          "failed_seeds": sum(result.status == "failed" for result in results),
                          "completed_seeds": sum(result.status == "complete" for result in results),
                          "combined_rows": row_count, "candidate_rows": candidate_count},
               "incomplete_seeds": [result.seed for result in results if result.status != "complete"],
               "seeds": [{**seed_result_row(result), "command": result.command, "best_point": result.best_point} for result in results],
               "best_points": best_rows}
    atomic_write_json(directory / "campaign_summary.json", json_safe(payload))
    return CampaignResult(results, TSVRows(combined_path, row_count), best_rows)


def run_campaign(args):
    canonicalize_paths(args)
    directory = args.campaign_dir
    directory.mkdir(parents=True, exist_ok=True)
    with CampaignLock(directory / "campaign_state.json", sys.argv):
        state_path = directory / "campaign_state.json"
        if args.resume or args.aggregate_only:
            state = load_json(state_path, "campaign state")
            if state.get("schema") != CAMPAIGN_SCHEMA:
                raise CampaignStateError("Unsupported campaign state schema")
            if state.get("hostname") != __import__("socket").gethostname():
                raise CampaignStateError("This campaign belongs to another host; cross-host resume is not supported")
            if state.get("configuration_fingerprint") != configuration_fingerprint(state["configuration"]):
                raise CampaignStateError("Saved campaign configuration fingerprint does not match")
            hydrate_configuration(args, state)
            if [record["seed"] for record in state["seeds"]] != expand_seeds(args.seed_start, args.nseeds):
                raise CampaignStateError("Saved seed list differs from the campaign configuration")
            for record in state["seeds"]:
                try:
                    seed_checkpoint(record)
                except (OSError, ValueError, KeyError, CampaignStateError):
                    pass
            assert_no_live_workers(state)
        else:
            if state_path.exists() or (directory / "campaign_summary.json").exists():
                raise CampaignStateError("Campaign already exists; use --resume or --aggregate-only")
            state = {"schema": CAMPAIGN_SCHEMA, "hostname": __import__("socket").gethostname(),
                     "created_utc": utc_now(), "configuration": configuration(args),
                     "configuration_fingerprint": configuration_fingerprint(configuration(args)),
                     "execution": {key: getattr(args, key) for key in EXECUTION_OPTIONS},
                     "status": "prepared", "seeds": [new_seed_record(seed, args)
                         for seed in expand_seeds(args.seed_start, args.nseeds)]}
            for record in state["seeds"]:
                cwd = Path(record["cwd"])
                if cwd.exists() and any(cwd.iterdir()):
                    raise CampaignStateError(f"Seed directory is not empty: {cwd}")
            persist_state(args, state)
        if args.aggregate_only:
            return aggregate(args, state)
        receipt = run_preflight(args)
        fingerprint = configuration_fingerprint(receipt)
        if state.get("runtime_fingerprint") not in (None, fingerprint):
            raise CampaignStateError("Runtime/source fingerprint changed; use the original installation or start a new campaign")
        state["runtime_fingerprint"] = fingerprint
        state["runtime_receipt"] = receipt
        for record in state["seeds"]:
            try:
                checkpoint = seed_checkpoint(record)
                if checkpoint:
                    validate_prefix(record["point_output"], checkpoint["outputs"]["main"])
                if target_complete(checkpoint, args):
                    record.update(status="complete", returncode=0, error="")
                elif record["status"] == "complete":
                    raise CampaignStateError("Completed seed lost its checkpoint")
                else:
                    record["status"] = "queued"
            except (OSError, ValueError, KeyError, CampaignStateError) as error:
                record.update(status="failed", returncode=1, error=str(error))
        remaining = sum(record["status"] == "queued" for record in state["seeds"])
        state["effective_jobs"] = effective_jobs(args.jobs, remaining)
        state["execution"] = {key: getattr(args, key) for key in EXECUTION_OPTIONS}
        state["thread_environment"] = THREAD_ENVIRONMENT
        state["status"] = "running"
        persist_state(args, state)
        print(f"Starting TRSM seed campaign: nseeds={args.nseeds} requested_jobs={args.jobs} "
              f"available_cpus={available_cpus()} effective_jobs={state['effective_jobs']} "
              f"target={args.nrandom} count={'evo/thc' if args.nrandom_count_evo_thc else 'draws'}", flush=True)
        if remaining:
            supervise(args, state)
        return aggregate(args, state)


def run_seed(seed, args):
    """Run one isolated seed; primarily a small integration-test entry point."""
    canonicalize_paths(args)
    record = new_seed_record(seed, args)
    job = launch_seed(record, args)
    try:
        job["process"].wait()
    finally:
        stop_processes([job["process"]], args.shutdown_grace_seconds)
    finish_seed(job, args)
    return seed_result(record, seed_checkpoint(record))


def main(argv=None):
    try:
        result = run_campaign(parse_args(argv))
    except (CampaignStateError, OSError, ValueError) as error:
        print(f"Campaign error: {error}", file=sys.stderr)
        return 1
    return 0 if all(seed.status == "complete" for seed in result.seed_results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
