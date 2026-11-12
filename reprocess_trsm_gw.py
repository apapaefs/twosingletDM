#!/usr/bin/env python3
"""Run CalcGW on saved EWPT-eligible points, without generating a new scan.

Streams the large campaign TSV; commits each point to SQLite. --max-points is
an optional pilot budget: a partial run only writes its checkpoint, never a
misleading finished TSV. Resume continues the same point identities.
"""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import csv
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import sqlite3

from trsm_gw import DEFAULT_GW_EXECUTABLE, GW_COLUMNS, run_calcgw


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def records(path):
    csv.field_size_limit(64 * 1024 * 1024)
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        header = reader.fieldnames
        required = {"M2", "M3", "vs", "vx", "a12", "lX", "lPhiX", "lSX", "ewpt_eligible"}
        if not header or len(header) != len(set(header)) or not required.issubset(header):
            raise ValueError("input needs unique columns, vx=0 TRSM parameters and stored ewpt_eligible")
        if set(header) & set(GW_COLUMNS):
            raise ValueError("input already contains GW columns; use the original source with --resume")
        yield header
        for index, row in enumerate(reader, 1):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"malformed input row {index}")
            if row["ewpt_eligible"] not in ("True", "False", "nan", "None", "null", ""):
                raise ValueError(f"invalid ewpt_eligible on row {index}")
            yield index, row


def reprocess(source, output, *, executable=DEFAULT_GW_EXECUTABLE, thigh=1000.,
              trans_temp="perc", vwall=.95, epsturb=.1, timeout=3600.,
              jobs=1, max_points=None, resume=False, dry_run=False, runner=run_calcgw):
    source, output, executable = (Path(p).expanduser().resolve() for p in (source, output, executable))
    if jobs < 1 or (max_points is not None and max_points < 1):
        raise ValueError("jobs and max-points must be positive")
    if not all(math.isfinite(v) for v in (thigh, vwall, epsturb, timeout)) or min(thigh, timeout) <= 0:
        raise ValueError("temperatures/timeouts must be positive and options finite")
    if vwall not in (-1, -2) and not 0 < vwall <= 1:
        raise ValueError("vwall must lie in (0,1] or be -1/-2 for BSMPT estimates")
    if epsturb != -1 and not 0 <= epsturb <= 1:
        raise ValueError("epsturb must lie in [0,1] or be -1 for the BSMPT upper bound")
    if trans_temp not in ("nucl_approx", "nucl", "perc", "compl"):
        raise ValueError("invalid transition temperature")
    sidecar = output.with_suffix(".metadata.json")
    checkpoint = Path(str(output) + ".partial")
    lock_path = Path(str(output) + ".lock")
    work = Path(str(output) + ".work")
    temporary = Path(str(output) + ".tmp")
    metadata_temp = Path(str(sidecar) + ".tmp")
    source_sidecar = source.with_suffix(".metadata.json")
    if {source, source_sidecar} & {output, sidecar, checkpoint, lock_path, work, temporary, metadata_temp} \
            or source.is_relative_to(work) or source_sidecar.is_relative_to(work):
        raise ValueError("input and output/checkpoint paths must not overlap")
    rows = records(source)
    header = next(rows)
    counts = Counter(rows=0, eligible=0)
    for _, row in rows:
        counts["rows"] += 1
        counts["eligible"] += row["ewpt_eligible"] == "True"
    print(f"Saved scan: {counts['rows']:,} rows; {counts['eligible']:,} EWPT-eligible for CalcGW", flush=True)
    if dry_run:
        return dict(counts)
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ValueError(f"build the existing BSMPT CalcGW target first: {executable}")
    root = Path(__file__).resolve().parent
    options = dict(executable=str(executable), thigh=thigh, trans_temp=trans_temp,
                   vwall=vwall, epsturb=epsturb, timeout=timeout)
    config = {"schema": "trsm_calcgw_v1", "input": str(source), "input_sha256": sha256(source),
              "input_metadata_sha256": sha256(source_sidecar) if source_sidecar.exists() else None,
              "output": str(output), "selection": "stored ewpt_eligible == True", "options": options,
              "executable_sha256": sha256(executable),
              "sources_sha256": {p: sha256(root / p) for p in
                  ("reprocess_trsm_gw.py", "trsm_gw.py", "test_trsm_ewpt.py", "trsm_inputs.py", "config/sm-inputs-v2.json")}}
    metadata = json.loads(source_sidecar.read_text()) if source_sidecar.exists() else {}
    if not isinstance(metadata, dict):
        raise ValueError("source metadata must be an object")
    output.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError("another CalcGW enrichment is writing this output") from error
        if output.exists():
            saved = json.loads(sidecar.read_text()).get("gw_enrichment", {}) if sidecar.exists() else {}
            if resume and saved.get("configuration") == config and saved.get("output_sha256") == sha256(output):
                return saved["counts"]
            raise FileExistsError(f"refusing to overwrite output: {output}")
        if sidecar.exists():
            saved = json.loads(sidecar.read_text()).get("gw_enrichment", {})
            if not resume or saved.get("configuration") != config:
                raise FileExistsError(f"refusing to overwrite metadata: {sidecar}")
        if resume != checkpoint.exists():
            raise ValueError("checkpoint exists: add --resume" if checkpoint.exists() else "no checkpoint to resume")
        if work.exists() and not resume:
            raise FileExistsError(f"refusing to reuse existing work directory: {work}")
        connection = sqlite3.connect(checkpoint)
        try:
            if resume:
                saved = json.loads(connection.execute("SELECT value FROM config").fetchone()[0])
                if saved != config:
                    raise ValueError("checkpoint does not match input, code, executable, or physics options")
            else:
                with connection:
                    connection.execute("CREATE TABLE config (value TEXT)")
                    connection.execute("INSERT INTO config VALUES (?)", (json.dumps(config, sort_keys=True),))
                    connection.execute("CREATE TABLE results (row_index INTEGER PRIMARY KEY, value TEXT)")
            done = {index for index, in connection.execute("SELECT row_index FROM results")}
            new_count = 0

            def save_completed(pending):
                nonlocal new_count
                finished, _ = wait(pending, return_when=FIRST_COMPLETED)
                for future in finished:
                    index = pending.pop(future)
                    result = future.result()
                    payload = json.dumps(result, sort_keys=True, allow_nan=False)
                    with connection:
                        connection.execute("INSERT INTO results VALUES (?, ?)", (index, payload))
                    done.add(index)
                    new_count += 1
                    print(f"CalcGW {len(done):,}/{counts['eligible']:,}; row {index}: {result['gw_status']}", flush=True)

            with ThreadPoolExecutor(max_workers=jobs) as pool:
                pending = {}
                submitted = 0
                rows = records(source)
                next(rows)
                for index, row in rows:
                    if row["ewpt_eligible"] != "True" or index in done:
                        continue
                    if max_points is not None and submitted >= max_points:
                        break
                    pending[pool.submit(runner, row, index, work / f"row-{index:07d}", options)] = index
                    submitted += 1
                    if len(pending) >= jobs:
                        save_completed(pending)
                while pending:
                    save_completed(pending)
            counts.update({"processed": len(done), "newly_processed": new_count})
            for value, in connection.execute("SELECT value FROM results"):
                result = json.loads(value)
                counts["status_" + result["gw_status"]] += 1
                counts["successful_transitions"] += result["gw_successful_transition_count"]
            if len(done) != counts["eligible"]:
                print("Pilot finished; checkpoint saved. Repeat without --max-points and with --resume to finish.", flush=True)
                return dict(counts, complete=False)
            if sha256(source) != config["input_sha256"] or (
                sha256(source_sidecar) if source_sidecar.exists() else None
            ) != config["input_metadata_sha256"]:
                raise ValueError("source changed during CalcGW; refusing to publish a mixed output")
            rows = records(source)
            next(rows)
            with temporary.open("w", encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=header + list(GW_COLUMNS), delimiter="\t", lineterminator="\n")
                writer.writeheader()
                for index, row in rows:
                    saved = connection.execute("SELECT value FROM results WHERE row_index=?", (index,)).fetchone()
                    result = json.loads(saved[0]) if saved else {"gw_status": "not_selected"}
                    row.update({key: json.dumps(value, separators=(",", ":"), allow_nan=False)
                                if isinstance(value, (dict, list)) else value for key, value in result.items()})
                    writer.writerow(row)
                stream.flush()
                os.fsync(stream.fileno())
            counts["complete"] = True
            metadata["gw_enrichment"] = {"configuration": config, "counts": dict(counts),
                                          "output_sha256": sha256(temporary)}
            metadata_temp.write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
            os.replace(metadata_temp, sidecar)
            os.replace(temporary, output)
            return dict(counts)
        finally:
            connection.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--executable", type=Path, default=DEFAULT_GW_EXECUTABLE)
    parser.add_argument("--thigh", type=float, default=1000., help="Match the saved campaign (default 1000 GeV)")
    parser.add_argument("--trans-temp", choices=("nucl_approx", "nucl", "perc", "compl"), default="perc")
    parser.add_argument("--vwall", type=float, default=.95)
    parser.add_argument("--epsturb", type=float, default=.1)
    parser.add_argument("--timeout", type=float, default=3600., help="Seconds per point")
    parser.add_argument("--jobs", type=int, default=1, help="Concurrent single-threaded CalcGW processes")
    parser.add_argument("--max-points", type=int, help="Optional new-point pilot budget; resume afterward")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = vars(parser.parse_args())
    source = args.pop("input")
    print(json.dumps(reprocess(source, **args), indent=2))


if __name__ == "__main__":
    main()
