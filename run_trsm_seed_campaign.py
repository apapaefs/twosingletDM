#!/usr/bin/env python3

import argparse
import csv
import json
import math
import os
import re
import shlex
import subprocess
import sys
import time
from collections import namedtuple
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from datetime import date
from pathlib import Path


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
    ],
    defaults=[0, 0, 0, math.nan, "", "", "", "", 0.0, [], None],
)
CampaignResult = namedtuple(
    "CampaignResult",
    ["seed_results", "combined_points", "best_points"],
)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Run generate_trsm_points.py for a range of seeds and aggregate viable/EWPT results."
    )
    parser.add_argument("--seed-start", type=int, required=True)
    parser.add_argument("--nseeds", type=int, required=True)
    parser.add_argument("--nrandom", type=int, default=100)
    parser.add_argument("--jobs", type=int, default=4)
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
        default=Path.cwd(),
        help="Working directory for seed subprocesses; generate_trsm_points.py writes output/ here.",
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
    args = parser.parse_args(argv)

    if args.nseeds <= 0:
        parser.error("--nseeds must be positive")
    if args.nrandom < 0:
        parser.error("--nrandom must be non-negative")
    if args.jobs <= 0:
        parser.error("--jobs must be positive")
    if args.heartbeat_seconds < 0:
        parser.error("--heartbeat-seconds must be non-negative")
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
    if args.run_ewpt:
        command.extend(["--run-ewpt", "--ewpt-workdir", str(ewpt_workdir)])
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
    return command


def read_viable_rows(path, seed):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
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


def collect_best_point_records(seed, seed_ewpt_dir):
    records = []
    for summary_path in sorted(Path(seed_ewpt_dir).glob("point_*/ewpt_summary.txt")):
        strength = best_strength_from_summary(summary_path)
        if math.isnan(strength.value):
            continue
        point_info = read_trsm_input(summary_path.parent / "TRSM_Input.tsv")
        records.append(
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
    return records


def best_record(records):
    if not records:
        return None
    return sorted(
        records,
        key=lambda record: (
            STRENGTH_PRIORITY.index(record["strength_kind"]),
            -float(record["best_ew_jump_over_T"]),
        ),
    )[0]


def format_float(value):
    if value is None or math.isnan(value):
        return "nan"
    return f"{value:.12g}"


def write_seed_log(log_path, command, completed, elapsed):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(
        "\n".join(
            [
                "command: " + shlex.join(command),
                f"returncode: {completed.returncode}",
                f"elapsed_seconds: {elapsed:.6f}",
                "",
                "stdout:",
                completed.stdout or "",
                "",
                "stderr:",
                completed.stderr or "",
                "",
            ]
        ),
        encoding="utf-8",
    )


def run_subprocess_to_log(command, cwd, log_path):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    with log_path.open("w", encoding="utf-8", buffering=1) as log:
        log.write("command: " + shlex.join(command) + "\n\n")
        log.write("subprocess output:\n")
        log.flush()
        process = subprocess.Popen(
            command,
            cwd=cwd,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        returncode = process.wait()
        elapsed = time.monotonic() - start
        log.write("\n")
        log.write(f"returncode: {returncode}\n")
        log.write(f"elapsed_seconds: {elapsed:.6f}\n")
    return returncode, elapsed


def run_seed(seed, args):
    campaign_dir = Path(args.campaign_dir).expanduser().resolve()
    run_cwd = Path(args.run_cwd).expanduser().resolve()
    run_date = date.today()
    log_path = campaign_dir / "logs" / f"seed_{seed}.log"
    ewpt_workdir = campaign_dir / "ewpt" / f"seed_{seed}"
    point_output = seed_output_path(run_cwd, seed, run_date)
    if point_output.exists():
        point_output.unlink()
    dm_failed_output = seed_dm_failed_output_path(run_cwd, seed, run_date)
    if args.write_dm_failed and dm_failed_output.exists():
        dm_failed_output.unlink()
    generator_script = Path(args.generator_script).expanduser().resolve()
    command = build_generator_command(args, seed, generator_script, ewpt_workdir)

    returncode, elapsed = run_subprocess_to_log(command, run_cwd, log_path)

    viable_rows = read_viable_rows(point_output, seed)
    summaries = sorted(ewpt_workdir.glob("point_*/ewpt_summary.txt"))
    best_points = collect_best_point_records(seed, ewpt_workdir)
    best = best_record(best_points)
    best_value = float(best["best_ew_jump_over_T"]) if best is not None else math.nan
    best_kind = best["strength_kind"] if best is not None else ""

    return SeedResult(
        seed=seed,
        returncode=returncode,
        viable_count=len(viable_rows),
        ewpt_runs=len(summaries),
        best_ew_jump_over_T=best_value,
        best_strength_kind=best_kind,
        point_output=str(point_output),
        log_path=str(log_path),
        ewpt_workdir=str(ewpt_workdir),
        elapsed_seconds=elapsed,
        command=command,
        best_point=best,
    )


def run_seed_jobs(args, seeds, run_seed_func=run_seed, executor_cls=ThreadPoolExecutor):
    with executor_cls(max_workers=args.jobs) as executor:
        future_to_seed = {}
        for seed in seeds:
            print(f"[launch] seed={seed}", flush=True)
            future_to_seed[executor.submit(run_seed_func, seed, args)] = seed

        pending = set(future_to_seed)
        while pending:
            timeout = args.heartbeat_seconds if args.heartbeat_seconds > 0 else None
            done, pending = wait(pending, timeout=timeout, return_when=FIRST_COMPLETED)
            if not done:
                active = ",".join(str(future_to_seed[future]) for future in sorted(pending, key=lambda item: future_to_seed[item]))
                print(f"[running] active={len(pending)} seeds={active}", flush=True)
                continue
            for future in done:
                yield future.result()


def write_tsv(path, rows, fieldnames):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def combine_viable_points(seed_results):
    rows = []
    for result in sorted(seed_results, key=lambda item: item.seed):
        rows.extend(read_viable_rows(result.point_output, result.seed))
    return rows


def seed_result_row(result):
    return {
        "seed": result.seed,
        "returncode": result.returncode,
        "viable_count": result.viable_count,
        "ewpt_runs": result.ewpt_runs,
        "best_strength_kind": result.best_strength_kind,
        "best_ew_jump_over_T": format_float(result.best_ew_jump_over_T),
        "point_output": result.point_output,
        "log_path": result.log_path,
        "ewpt_workdir": result.ewpt_workdir,
        "elapsed_seconds": f"{result.elapsed_seconds:.6f}",
    }


def progress_line(done, total, result, elapsed):
    return (
        f"[{done}/{total} done] seed={result.seed} viable={result.viable_count} "
        f"ewpt_runs={result.ewpt_runs} "
        f"best_ew_jump/T={format_float(result.best_ew_jump_over_T)} "
        f"elapsed={elapsed:.1f}s"
    )


def write_campaign_outputs(args, seed_results):
    campaign_dir = Path(args.campaign_dir).expanduser().resolve()
    campaign_dir.mkdir(parents=True, exist_ok=True)

    summary_rows = [seed_result_row(result) for result in sorted(seed_results, key=lambda item: item.seed)]
    summary_fields = [
        "seed",
        "returncode",
        "viable_count",
        "ewpt_runs",
        "best_strength_kind",
        "best_ew_jump_over_T",
        "point_output",
        "log_path",
        "ewpt_workdir",
        "elapsed_seconds",
    ]
    write_tsv(campaign_dir / "campaign_summary.tsv", summary_rows, summary_fields)

    combined_rows = combine_viable_points(seed_results)
    combined_fields = list(combined_rows[0].keys()) if combined_rows else ["seed"]
    write_tsv(campaign_dir / "combined_points.tsv", combined_rows, combined_fields)

    best_rows = [
        result.best_point for result in sorted(seed_results, key=lambda item: (
            STRENGTH_PRIORITY.index(item.best_strength_kind)
            if item.best_strength_kind in STRENGTH_PRIORITY
            else len(STRENGTH_PRIORITY),
            -item.best_ew_jump_over_T if not math.isnan(item.best_ew_jump_over_T) else 0.0,
            item.seed,
        ))
        if result.best_point is not None
    ]
    best_fields = [
        "seed",
        "point_index",
        "strength_kind",
        "transition_index",
        "best_ew_jump_over_T",
        "ewpt_summary",
        "m1",
        "m2",
        "m3",
        "vs",
        "a12",
        "lx",
        "lphix",
        "lsx",
    ]
    write_tsv(campaign_dir / "best_points.tsv", best_rows, best_fields)

    payload = {
        "metadata": {
            "seed_start": args.seed_start,
            "nseeds": args.nseeds,
            "nrandom": args.nrandom,
            "jobs": args.jobs,
            "campaign_dir": str(campaign_dir),
            "generator_script": str(Path(args.generator_script).expanduser().resolve()),
            "run_cwd": str(Path(args.run_cwd).expanduser().resolve()),
            "write_dm_failed": args.write_dm_failed,
            "run_ewpt": args.run_ewpt,
            "ewpt_require_eq418": args.ewpt_require_eq418,
            "ewpt_thigh": args.ewpt_thigh,
        },
        "totals": {
            "viable_count": sum(result.viable_count for result in seed_results),
            "ewpt_runs": sum(result.ewpt_runs for result in seed_results),
            "failed_seeds": sum(1 for result in seed_results if result.returncode != 0),
        },
        "seeds": [
            {
                **seed_result_row(result),
                "command": result.command,
                "best_point": result.best_point,
            }
            for result in sorted(seed_results, key=lambda item: item.seed)
        ],
        "best_points": best_rows,
    }
    (campaign_dir / "campaign_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=True) + "\n",
        encoding="ascii",
    )
    return combined_rows, best_rows


def run_campaign(args):
    campaign_dir = Path(args.campaign_dir).expanduser().resolve()
    campaign_dir.mkdir(parents=True, exist_ok=True)
    seeds = expand_seeds(args.seed_start, args.nseeds)
    start = time.monotonic()
    seed_results = []
    print(
        f"Starting TRSM seed campaign: seeds={seeds[0]}..{seeds[-1]} "
        f"nseeds={len(seeds)} jobs={args.jobs} logs={campaign_dir / 'logs'}",
        flush=True,
    )
    for result in run_seed_jobs(args, seeds):
        seed_results.append(result)
        print(progress_line(len(seed_results), len(seeds), result, time.monotonic() - start), flush=True)
    combined_rows, best_rows = write_campaign_outputs(args, seed_results)
    return CampaignResult(
        seed_results=sorted(seed_results, key=lambda item: item.seed),
        combined_points=combined_rows,
        best_points=best_rows,
    )


def main(argv=None):
    args = parse_args(argv)
    run_campaign(args)


if __name__ == "__main__":
    main()
