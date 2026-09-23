#!/usr/bin/env python3

"""Run BSMPT EWPT checks on points stored in an existing TRSM scan.

The source TSV is never modified.  Rows pass the EWPT selection when every
stored non-DM constraint passes; the stored aggregate DM result is preserved
and counted but deliberately does not gate BSMPT.  Progress is transactionally
checkpointed in SQLite and the requested output TSV appears only after every
input row has been accounted for.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sqlite3
import sys
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Mapping, Sequence

from dm_thermal_relic_diagnostic import (
    RESONANCE_COLUMNS,
    THERMAL_VEV_COLUMNS, THERMAL_COUPLING_COLUMNS, thermal_input_updates,
    resonance_proximity_updates,
    thermal_vev_updates,
)
from ewpt_entry_criterion import EW_ENTRY_COLUMNS, ew_entry_updates
from ewpt_x_history import X_HISTORY_COLUMNS, x_history_updates


CHECKPOINT_SCHEMA_VERSION = 1
REPROCESSOR_SCHEMA = "trsm_ewpt_reprocessing_v2"
from trsm_inputs import M1 as M1_GEV, PHYSICS_VERSION
from ewpt_assessment import STATUS_COLUMNS, status_updates
NON_DM_CONSTRAINT_COLUMNS = ("thc", "hb", "hs", "ewpo", "wmass")
POINT_COLUMNS = ("M2", "M3", "vs", "vx", "a12", "lX", "lPhiX", "lSX")
REQUIRED_INPUT_COLUMNS = POINT_COLUMNS + NON_DM_CONSTRAINT_COLUMNS + ("evo","dm",)
EWPT_COLUMNS = (
    "ewpt_ew_true_over_T",
    "ewpt_ew_jump_over_T",
    "ewpt_global_phase_path",
    "ewpt_has_x_broken",
    "ewpt_ew_step_index",
    "ewpt_status",
    "ewpt_error",
    *EW_ENTRY_COLUMNS,
    *X_HISTORY_COLUMNS,
    *THERMAL_VEV_COLUMNS,
    *RESONANCE_COLUMNS,
    *THERMAL_COUPLING_COLUMNS, *STATUS_COLUMNS, "ewpt_constraint_version",
)
EWPT_STRENGTH_PRIORITY = ("nucl", "perc", "compl", "crit")
ELIGIBLE_OUTCOMES = {
    "success",
    "failed",
    "existing",
    "skipped_eq418",
}


class EWPTReprocessingError(RuntimeError):
    """Raised when an existing scan cannot be reprocessed safely."""


@dataclass(frozen=True)
class RowEvaluation:
    updates: Mapping[str, object]
    outcome: str
    dm_passed: bool


@dataclass(frozen=True)
class ReprocessingResult:
    output: Path
    workdir: Path
    provenance: Path
    metadata: Path | None
    counts: Mapping[str, int]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def json_safe(value):
    if isinstance(value,float) and not math.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def one_line_error(error: BaseException, max_length: int = 1000) -> str:
    message = " ".join(str(error).split())
    if len(message) > max_length:
        return message[:max_length] + "... [truncated]"
    return message


def optional_text(value: object) -> str:
    text = "" if value is None else str(value).strip()
    if text.lower() in {"", "nan", "none", "null"}:
        return ""
    return text


def finite_number(value: object) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def require_float(row: Mapping[str, str], column: str, row_number: int) -> float:
    value = finite_number(row.get(column))
    if value is None:
        raise EWPTReprocessingError(
            f"input row {row_number} has an invalid {column!r} value"
        )
    return value


def require_bool(row: Mapping[str, str], column: str, row_number: int) -> bool:
    value = row.get(column)
    if value == "True":
        return True
    if value == "False":
        return False
    raise EWPTReprocessingError(
        f"input row {row_number} has {column}={value!r}; expected True or False"
    )


def existing_ewpt_attempt(row: Mapping[str, str]) -> bool:
    if optional_text(row.get("ewpt_status")):
        return True
    if optional_text(row.get("ewpt_error")):
        return True
    if finite_number(row.get("ewpt_ew_true_over_T")) is not None:
        return True
    if finite_number(row.get("ewpt_ew_jump_over_T")) is not None:
        return True
    if finite_number(row.get("ewpt_ew_entry_true_over_T")) is not None:
        return True
    if optional_text(row.get("ewpt_global_phase_path")):
        return True
    if finite_number(row.get("ewpt_ew_step_index")) is not None:
        return True
    return str(row.get("ewpt_has_x_broken", "")).strip() in {"True", "False"}


def empty_ewpt_updates() -> dict[str, object]:
    return {
        "ewpt_ew_true_over_T": None,
        "ewpt_ew_jump_over_T": None,
        **dict.fromkeys(EW_ENTRY_COLUMNS),
        "ewpt_global_phase_path": None,
        "ewpt_has_x_broken": None,
        "ewpt_ew_step_index": None,
        "ewpt_status": None,
        "ewpt_error": None,
        **dict.fromkeys(X_HISTORY_COLUMNS),
        **dict.fromkeys(THERMAL_VEV_COLUMNS),
        **dict.fromkeys(RESONANCE_COLUMNS),
    }


def select_primary_strength(payload: Mapping[str, object]) -> Mapping[str, object] | None:
    strengths = payload.get("transition_strengths") or []
    for temperature_kind in EWPT_STRENGTH_PRIORITY:
        candidates = []
        for strength in strengths:
            if not isinstance(strength, Mapping):
                continue
            if strength.get("temperature_kind") != temperature_kind:
                continue
            value = finite_number(strength.get("ew_true_over_T"))
            if value is not None:
                candidates.append((value, strength))
        if candidates:
            return max(candidates, key=lambda item: item[0])[1]
    return None


def updates_from_payload(
    payload: Mapping[str, object], *, w1_threshold: float = 5.0,
    freezeout_temperature: float | None = None,
    m2: float | None = None, m3: float | None = None, point=None,
) -> dict[str, object]:
    updates = empty_ewpt_updates()
    updates.update(status_updates(payload))
    updates["ewpt_constraint_version"]=PHYSICS_VERSION
    updates["ewpt_error"] = ""

    strength = select_primary_strength(payload)
    if strength is not None:
        updates["ewpt_ew_true_over_T"] = finite_number(
            strength.get("ew_true_over_T")
        )
        updates["ewpt_ew_jump_over_T"] = finite_number(
            strength.get("ew_jump_over_T")
        )
    updates.update(ew_entry_updates(payload, w1_threshold=w1_threshold))
    updates.update(x_history_updates(payload, freezeout_temperature))
    updates.update(thermal_vev_updates(payload, freezeout_temperature))
    updates.update(resonance_proximity_updates(m2, m3, freezeout_temperature))
    if point is not None:
        updates.update(thermal_input_updates(payload, point))

    minimatracer = payload.get("minimatracer") or {}
    if isinstance(minimatracer, Mapping):
        raw_path = minimatracer.get("global_phase_path") or []
        if isinstance(raw_path, str):
            labels = [part.strip() for part in raw_path.split("->") if part.strip()]
        else:
            labels = [str(part).strip() for part in raw_path if str(part).strip()]
        if labels:
            updates["ewpt_global_phase_path"] = " -> ".join(labels)
            updates["ewpt_has_x_broken"] = any(
                "X_BROKEN" in label for label in labels
            )
        step = minimatracer.get("ew_step_index")
        if step is not None:
            updates["ewpt_ew_step_index"] = step
    return updates


class ProductionEvaluator:
    """Adapter around :mod:`test_trsm_ewpt` for one stored scan row."""

    def __init__(
        self,
        ewpt_module,
        config,
        workdir: Path,
        *,
        require_eq418: bool = False,
        rerun_existing: bool = False,
    ):
        self.ewpt = ewpt_module
        self.config = config
        self.workdir = Path(workdir)
        self.require_eq418 = require_eq418
        self.rerun_existing = rerun_existing

    def __call__(self, row: Mapping[str, str], point_index: int) -> RowEvaluation:
        row_number = point_index + 1
        def verdict(name):
            value=row.get(name)
            if value not in (True,False,"True","False",None,"nan","None","null",""):
                raise EWPTReprocessingError(f"input row {row_number} has invalid {name}={value!r}")
            return True if value in (True,"True") else False if value in (False,"False") else None
        dm_passed=verdict("dm")
        known_exclusion=any(verdict(k) is False for k in ("hb","hs","ewpo"))
        mass=finite_number(row.get("M2"))
        known_exclusion |= mass is not None and 133<=mass<=999 and verdict("wmass") is False
        if verdict("thc") is not True or known_exclusion:
            return RowEvaluation({}, "ineligible", dm_passed)
        if not self.rerun_existing and row.get("ewpt_constraint_version")==PHYSICS_VERSION and existing_ewpt_attempt(row):
            return RowEvaluation({}, "existing", dm_passed)

        m2 = require_float(row, "M2", row_number)
        m3 = require_float(row, "M3", row_number)
        vs = require_float(row, "vs", row_number)
        vx = require_float(row, "vx", row_number)
        a12 = require_float(row, "a12", row_number)
        lx = require_float(row, "lX", row_number)
        lphix = require_float(row, "lPhiX", row_number)
        lsx = require_float(row, "lSX", row_number)
        if not math.isclose(vx, 0.0, rel_tol=0.0, abs_tol=1.0e-12):
            raise EWPTReprocessingError(
                f"input row {row_number} has vx={vx}; only vx=0 scans are supported"
            )

        point = self.ewpt.TRSMEWPTPoint(
            index=point_index,
            m1=M1_GEV,
            m2=m2,
            m3=m3,
            vs=vs,
            a12=a12,
            lx=lx,
            lphix=lphix,
            lsx=lsx,
        )
        eq418_row = {
            "m1": M1_GEV,
            "m2": m2,
            "m3": m3,
            "vs": vs,
            "a12": a12,
            "lx": lx,
            "lphix": lphix,
            "lsx": lsx,
        }
        if self.require_eq418:
            check = self.ewpt.check_eq_4_18(eq418_row)
            if not check.satisfied:
                print(f"Row {point_index}: Eq. 4.18 failed; skipping BSMPT")
                return RowEvaluation({}, "skipped_eq418", dm_passed)

        point_workdir = self.workdir / f"point_{point_index:06d}"
        print(
            f"Row {point_index}: non-DM constraints pass, dm={dm_passed}; "
            f"running BSMPT in {point_workdir}"
        )
        try:
            result = self.ewpt.run_trsm_ewpt(
                point,
                config=self.config,
                workdir=point_workdir,
                keep_files=True,
            )
        except Exception as error:
            point_workdir.mkdir(parents=True, exist_ok=True)
            message = one_line_error(error)
            (point_workdir / "ewpt_error.txt").write_text(
                message + "\n", encoding="utf-8"
            )
            updates = empty_ewpt_updates()
            updates["ewpt_status"] = "failed"
            updates["ewpt_error"] = message
            print(f"Row {point_index}: BSMPT failed; continuing: {message}")
            return RowEvaluation(updates, "failed", dm_passed)

        summary = self.ewpt.summarize_result(result)
        payload = self.ewpt.result_to_json(result)
        point_workdir.mkdir(parents=True, exist_ok=True)
        (point_workdir / "ewpt_summary.txt").write_text(
            summary + "\n", encoding="utf-8"
        )
        (point_workdir / "ewpt_result.json").write_text(
            json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
            encoding="utf-8",
        )
        print(summary)
        return RowEvaluation(
            updates_from_payload(
                payload,
                w1_threshold=self.config.w1_threshold,
                freezeout_temperature=finite_number(row.get("dm_freezeout_temperature_GeV")),
                m2=m2,
                m3=m3, point=row,
            ),
            "success",
            dm_passed,
        )


def read_input(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise EWPTReprocessingError(f"input is empty: {path}") from error
        if not header or any(not column for column in header):
            raise EWPTReprocessingError("input header contains an empty column name")
        duplicates = sorted({name for name in header if header.count(name) > 1})
        if duplicates:
            raise EWPTReprocessingError(
                f"input has duplicate columns: {', '.join(duplicates)}"
            )
        missing = [name for name in REQUIRED_INPUT_COLUMNS if name not in header]
        if missing:
            raise EWPTReprocessingError(
                f"input is missing columns: {', '.join(missing)}"
            )
        rows = []
        for row_number, values in enumerate(reader, start=2):
            if len(values) != len(header):
                raise EWPTReprocessingError(
                    f"input row {row_number} has {len(values)} fields; "
                    f"expected {len(header)}"
                )
            rows.append(dict(zip(header, values)))
    if not rows:
        raise EWPTReprocessingError(f"input has a header but no data rows: {path}")
    return header, rows


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while True:
            chunk = stream.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def input_identity(
    path: Path, header: Sequence[str], row_count: int
) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "sha256": sha256_file(path),
        "size": path.stat().st_size,
        "header": list(header),
        "row_count": int(row_count),
    }


def output_header(input_header: Sequence[str]) -> list[str]:
    result = list(input_header)
    result.extend(column for column in EWPT_COLUMNS if column not in result)
    return result


def format_value(value: object) -> str:
    if value is None:
        return "nan"
    if type(value) is bool:
        return "True" if value else "False"
    return str(value)


def merged_payload(
    input_row: Mapping[str, str],
    updates: Mapping[str, object],
    header: Sequence[str],
) -> str:
    fields = []
    for column in header:
        if column in updates:
            fields.append(format_value(updates[column]))
        else:
            fields.append(input_row.get(column, "nan"))
    return "\t".join(fields)


def checkpoint_path(output: Path) -> Path:
    return output.with_name(output.name + ".partial")


def provenance_path(output: Path) -> Path:
    return output.with_name(output.stem + ".ewpt-reprocess.json")


def scan_metadata_path(path: Path) -> Path:
    return path.with_suffix(".metadata.json")


def create_checkpoint(path: Path, metadata: Mapping[str, object]) -> sqlite3.Connection:
    connection = sqlite3.connect(path)
    try:
        with connection:
            connection.execute(
                "CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)"
            )
            connection.execute(
                "CREATE TABLE rows ("
                "row_index INTEGER PRIMARY KEY, payload TEXT NOT NULL, "
                "outcome TEXT NOT NULL, dm_passed INTEGER)"
            )
            connection.executemany(
                "INSERT INTO metadata(key, value) VALUES (?, ?)",
                [
                    (
                        key,
                        json.dumps(value, sort_keys=True, separators=(",", ":")),
                    )
                    for key, value in metadata.items()
                ],
            )
    except Exception:
        connection.close()
        raise
    return connection


def open_checkpoint(
    path: Path, expected_metadata: Mapping[str, object]
) -> tuple[sqlite3.Connection, int]:
    try:
        connection = sqlite3.connect(path)
        stored = {
            key: json.loads(value)
            for key, value in connection.execute("SELECT key, value FROM metadata")
        }
        if stored != dict(expected_metadata):
            raise EWPTReprocessingError(
                "checkpoint does not match the current input, output, workdir, "
                "or EWPT configuration"
            )
        count, minimum, maximum = connection.execute(
            "SELECT COUNT(*), MIN(row_index), MAX(row_index) FROM rows"
        ).fetchone()
        if count and (minimum != 0 or maximum != count - 1):
            raise EWPTReprocessingError("checkpoint rows are not contiguous")
        return connection, int(count)
    except (sqlite3.DatabaseError, json.JSONDecodeError, EWPTReprocessingError) as error:
        try:
            connection.close()
        except UnboundLocalError:
            pass
        if isinstance(error, EWPTReprocessingError):
            raise
        raise EWPTReprocessingError(f"invalid checkpoint {path}: {error}") from error


def checkpoint_counts(connection: sqlite3.Connection) -> Counter:
    return Counter(
        {
            outcome: int(count)
            for outcome, count in connection.execute(
                "SELECT outcome, COUNT(*) FROM rows GROUP BY outcome"
            )
        }
    )


def export_completed_output(
    connection: sqlite3.Connection,
    output: Path,
    header: Sequence[str],
    expected_rows: int,
) -> None:
    temporary = output.with_name("." + output.name + ".complete.tmp")
    count = connection.execute("SELECT COUNT(*) FROM rows").fetchone()[0]
    if count != expected_rows:
        raise EWPTReprocessingError(
            f"checkpoint contains {count} rows; expected {expected_rows}"
        )
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        stream.write("\t".join(header) + "\n")
        for expected_index, (row_index, payload) in enumerate(
            connection.execute(
                "SELECT row_index, payload FROM rows ORDER BY row_index"
            )
        ):
            if row_index != expected_index:
                raise EWPTReprocessingError("checkpoint rows are not contiguous")
            stream.write(payload + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, output)


def write_json_atomic(path: Path, payload: Mapping[str, object]) -> None:
    temporary = path.with_name("." + path.name + ".tmp")
    temporary.write_text(
        json.dumps(json_safe(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def load_scan_metadata(input_path: Path, output: Path) -> dict | None:
    source = scan_metadata_path(input_path)
    if not source.is_file():
        return None
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise EWPTReprocessingError(
            f"cannot read source scan metadata {source}: {error}"
        ) from error
    if not isinstance(payload, dict):
        raise EWPTReprocessingError(f"source scan metadata is not an object: {source}")
    destination = scan_metadata_path(output)
    if destination.resolve() == source.resolve():
        raise EWPTReprocessingError(
            "input and output names resolve to the same metadata sidecar"
        )
    history = payload.get("postprocessing_history", [])
    if not isinstance(history, list):
        raise EWPTReprocessingError(
            f"source scan metadata has invalid postprocessing_history: {source}"
        )
    return payload


def write_scan_metadata(
    input_path: Path,
    output: Path,
    identity: Mapping[str, object],
    history_entry: Mapping[str, object],
    payload: dict | None,
) -> Path | None:
    if payload is None:
        return None
    destination = scan_metadata_path(output)
    payload["scan_file"] = output.name
    payload["source_scan"] = {
        "path": str(input_path),
        "sha256": identity["sha256"],
    }
    history = payload.setdefault("postprocessing_history", [])
    history.append(dict(history_entry))
    write_json_atomic(destination, payload)
    return destination


def ensure_fresh_workdir(workdir: Path, resume: bool) -> None:
    if resume:
        workdir.mkdir(parents=True, exist_ok=True)
        return
    if workdir.exists() and any(workdir.iterdir()):
        raise FileExistsError(
            f"EWPT workdir is not empty; choose a new directory or use --resume: "
            f"{workdir}"
        )
    workdir.mkdir(parents=True, exist_ok=True)


def reprocess(
    input_path: Path,
    output: Path,
    evaluator: Callable[[Mapping[str, str], int], RowEvaluation],
    *,
    workdir: Path,
    settings: Mapping[str, object],
    command_line: Sequence[str] = (),
    resume: bool = False,
    checkpoint_every: int = 1,
) -> ReprocessingResult:
    input_path = input_path.expanduser().resolve()
    output = output.expanduser().resolve()
    workdir = workdir.expanduser().resolve()
    if input_path == output:
        raise EWPTReprocessingError("input and output paths must differ")
    if not input_path.is_file():
        raise FileNotFoundError(f"input scan not found: {input_path}")
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    if checkpoint_every <= 0:
        raise ValueError("checkpoint_every must be positive")

    header, rows = read_input(input_path)
    identity = input_identity(input_path, header, len(rows))
    final_header = output_header(header)
    source_metadata = load_scan_metadata(input_path, output)
    output.parent.mkdir(parents=True, exist_ok=True)
    ensure_fresh_workdir(workdir, resume)
    partial = checkpoint_path(output)
    provenance = provenance_path(output)
    destination_metadata = scan_metadata_path(output)
    if not resume:
        for path in (provenance, destination_metadata):
            if path.exists():
                raise FileExistsError(f"refusing to overwrite existing sidecar: {path}")

    checkpoint_metadata = {
        "schema_version": CHECKPOINT_SCHEMA_VERSION,
        "reprocessor_schema": REPROCESSOR_SCHEMA,
        "input_identity": identity,
        "output_header": final_header,
        "output_path": str(output),
        "workdir": str(workdir),
        "settings": json_safe(dict(settings)),
    }
    if resume:
        if not partial.is_file():
            raise FileNotFoundError(f"checkpoint not found: {partial}")
        connection, completed = open_checkpoint(partial, checkpoint_metadata)
    else:
        if partial.exists():
            raise FileExistsError(
                f"checkpoint already exists; use --resume or move it aside: {partial}"
            )
        connection = create_checkpoint(partial, checkpoint_metadata)
        completed = 0
    if completed > len(rows):
        connection.close()
        raise EWPTReprocessingError("checkpoint has more rows than the input")
    if completed:
        print(f"Resuming at input row {completed + 1:,}/{len(rows):,}")

    started = utc_now()
    try:
        for batch_start in range(completed, len(rows), checkpoint_every):
            batch_stop = min(batch_start + checkpoint_every, len(rows))
            connection.execute("BEGIN")
            try:
                for index in range(batch_start, batch_stop):
                    result = evaluator(rows[index], index + 1)
                    if result.outcome not in ELIGIBLE_OUTCOMES | {"ineligible"}:
                        raise EWPTReprocessingError(
                            f"row {index + 2} produced unknown outcome {result.outcome!r}"
                        )
                    updates = dict(result.updates)
                    if finite_number(rows[index].get("vx")) == 0.0:
                        updates.update(resonance_proximity_updates(
                            rows[index].get("M2"), rows[index].get("M3"),
                            rows[index].get("dm_freezeout_temperature_GeV"),
                            widths=(rows[index].get("dm_h1_width_GeV"),rows[index].get("dm_h2_width_GeV")),
                        ))
                    payload = merged_payload(rows[index], updates, final_header)
                    connection.execute(
                        "INSERT INTO rows(row_index, payload, outcome, dm_passed) "
                        "VALUES (?, ?, ?, ?)",
                        (index, payload, result.outcome, int(result.dm_passed) if result.dm_passed is not None else None),
                    )
                connection.commit()
            except BaseException:
                connection.rollback()
                raise
            counts = checkpoint_counts(connection)
            print(
                f"Checkpointed {batch_stop:,}/{len(rows):,} rows "
                f"(success={counts['success']}, failed={counts['failed']}, "
                f"existing={counts['existing']}, ineligible={counts['ineligible']})"
            )

        if input_identity(input_path, header, len(rows)) != identity:
            raise EWPTReprocessingError("input changed while EWPT reprocessing was running")
        export_completed_output(connection, output, final_header, len(rows))
        counts = checkpoint_counts(connection)
        dm_pass, dm_fail = connection.execute(
            "SELECT "
            "SUM(CASE WHEN outcome != 'ineligible' AND dm_passed = 1 THEN 1 ELSE 0 END), "
            "SUM(CASE WHEN outcome != 'ineligible' AND dm_passed = 0 THEN 1 ELSE 0 END) "
            "FROM rows"
        ).fetchone()
    finally:
        connection.close()

    completed_at = utc_now()
    for outcome in sorted(ELIGIBLE_OUTCOMES | {"ineligible"}):
        counts.setdefault(outcome, 0)
    counts["eligible_dm_pass"] = int(dm_pass or 0)
    counts["eligible_dm_fail"] = int(dm_fail or 0)
    counts["input_rows"] = len(rows)
    counts["eligible"] = sum(counts[name] for name in ELIGIBLE_OUTCOMES)
    history_entry = {
        "schema": REPROCESSOR_SCHEMA,
        "started_utc": started,
        "completed_utc": completed_at,
        "command_line": list(command_line),
        "settings": json_safe(dict(settings)),
        "counts": dict(counts),
        "workdir": str(workdir),
    }
    write_json_atomic(
        provenance,
        {
            **history_entry,
            "input": identity,
            "output": str(output),
        },
    )
    metadata = write_scan_metadata(
        input_path,
        output,
        identity,
        history_entry,
        source_metadata,
    )
    partial.unlink()
    print(f"Wrote {len(rows):,} rows to {output}")
    print(
        "EWPT selection: "
        f"eligible={counts['eligible']:,} "
        f"(dm pass={counts['eligible_dm_pass']:,}, "
        f"dm fail={counts['eligible_dm_fail']:,}); "
        f"success={counts['success']:,}, failed={counts['failed']:,}, "
        f"existing={counts['existing']:,}, "
        f"Eq.4.18 skipped={counts['skipped_eq418']:,}"
    )
    return ReprocessingResult(
        output=output,
        workdir=workdir,
        provenance=provenance,
        metadata=metadata,
        counts=dict(counts),
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run BSMPT on stored vx=0 scan rows that pass all non-DM "
            "constraints, without regenerating points or requiring dm=True."
        )
    )
    parser.add_argument("input", type=Path, help="Existing TRSM scan TSV")
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="New TSV receiving the preserved rows and updated EWPT columns",
    )
    parser.add_argument(
        "--ewpt-workdir",
        type=Path,
        help=(
            "Directory for point_NNNNNN BSMPT outputs; defaults to "
            "<output-stem>_ewpt beside the output TSV"
        ),
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume the transactionally checkpointed <output>.partial file",
    )
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=1,
        metavar="N",
        help="Commit a restartable checkpoint after N input rows (default: 1)",
    )
    parser.add_argument(
        "--rerun-existing-ewpt",
        action="store_true",
        help="Rerun rows that already contain an EWPT attempt instead of preserving it",
    )
    parser.add_argument(
        "--ewpt-executable",
        type=Path,
        help="CalcTemps executable passed to the EWPT runner",
    )
    parser.add_argument(
        "--ewpt-minima-executable",
        type=Path,
        help="MinimaTracer executable passed to the EWPT runner",
    )
    parser.add_argument("--ewpt-thigh", type=float, default=300.0)
    parser.add_argument("--ewpt-multistepmode", default="default")
    parser.add_argument("--ewpt-plot-phases", action="store_true")
    parser.add_argument("--ewpt-plot-output", type=Path)
    parser.add_argument(
        "--ewpt-plot-format",
        choices=["png", "pdf", "both"],
        default="both",
    )
    parser.add_argument("--ewpt-require-eq418", action="store_true")
    parser.add_argument("--ewpt-sym-threshold", type=float, default=1.0)
    parser.add_argument("--ewpt-w1-threshold", type=float, default=5.0)
    parser.add_argument("--ewpt-wx-threshold", type=float, default=1.0)
    parser.add_argument("--ewpt-ws-threshold", type=float, default=1.0)
    args = parser.parse_args(argv)
    if args.checkpoint_every < 1:
        parser.error("--checkpoint-every must be at least 1")
    if not math.isfinite(args.ewpt_thigh) or args.ewpt_thigh <= 0.0:
        parser.error("--ewpt-thigh must be positive and finite")
    return args


def run(
    argv: Sequence[str] | None = None,
    *,
    ewpt_module=None,
) -> ReprocessingResult:
    args = parse_args(argv)
    if ewpt_module is None:
        import test_trsm_ewpt as ewpt_module

    output = args.output.expanduser().resolve()
    workdir = args.ewpt_workdir
    if workdir is None:
        workdir = output.parent / f"{output.stem}_ewpt"
    workdir = Path(workdir).expanduser().resolve()
    config_kwargs = {
        "multistepmode": args.ewpt_multistepmode,
        "thigh": args.ewpt_thigh,
        "plot_phases": args.ewpt_plot_phases,
        "plot_output": args.ewpt_plot_output,
        "plot_format": args.ewpt_plot_format,
        "sym_threshold": args.ewpt_sym_threshold,
        "w1_threshold": args.ewpt_w1_threshold,
        "wx_threshold": args.ewpt_wx_threshold,
        "ws_threshold": args.ewpt_ws_threshold,
    }
    if args.ewpt_executable is not None:
        config_kwargs["executable"] = args.ewpt_executable
    if args.ewpt_minima_executable is not None:
        config_kwargs["minima_executable"] = args.ewpt_minima_executable
    config = ewpt_module.EWPTConfig(**config_kwargs)
    from trsm_constraint_profile import physics_manifest
    from trsm_micromegas import default_micromegas_main
    settings = {
        "physics_version":PHYSICS_VERSION,
        "physics_manifest":physics_manifest(str(default_micromegas_main()),str(getattr(config,"executable",getattr(ewpt_module,"DEFAULT_EXECUTABLE",Path("/unavailable/CalcTemps")))),str(getattr(config,"minima_executable",None) or Path(getattr(config,"executable",getattr(ewpt_module,"DEFAULT_EXECUTABLE",Path("/unavailable/CalcTemps")))).with_name("MinimaTracer"))),
        "selection": "evo & thc & hb & hs & ewpo & wmass; dm recorded but not required",
        "rerun_existing_ewpt": args.rerun_existing_ewpt,
        "ewpt_require_eq418": args.ewpt_require_eq418,
        **config_kwargs,
    }
    evaluator = ProductionEvaluator(
        ewpt_module,
        config,
        workdir,
        require_eq418=args.ewpt_require_eq418,
        rerun_existing=args.rerun_existing_ewpt,
    )
    command_line = [Path(sys.argv[0]).name, *(sys.argv[1:] if argv is None else argv)]
    return reprocess(
        args.input,
        output,
        evaluator,
        workdir=workdir,
        settings=settings,
        command_line=command_line,
        resume=args.resume,
        checkpoint_every=args.checkpoint_every,
    )


def main() -> None:
    run()


if __name__ == "__main__":
    main()
