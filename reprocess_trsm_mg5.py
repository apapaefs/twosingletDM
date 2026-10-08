#!/usr/bin/env python3
"""Add missing associated-production rates to a separate saved scan copy.

Only stored full-viability points enter MadGraph. Existing finite rates and
all other assessments are retained. A SQLite checkpoint commits each native
process result independently, so --resume never repeats a completed result.
"""

import argparse
import csv
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import sqlite3
import sys

from generate_mg5_trsm_xsecs import get_mg5_xsec, mg5_runtime_receipt
from trsm_mg5_rates import (
    ASSOCIATED_PROCESSES, DERIVED_RATE_COLUMNS, derive_mg5_rates,
    finite_number, stored_full_viability,
)

LAMBDA_COLUMNS = ("K111", "K112", "K113", "K123", "K122", "K1111", "K1112", "K1113", "K133", "K233")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_rows(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        header = reader.fieldnames
        if not header or any(not name for name in header) or len(set(header)) != len(header):
            raise ValueError("input needs a nonempty, unique TSV header")
        yield header
        for index, row in enumerate(reader, 1):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"malformed input data row {index}")
            yield index, row


def input_identity(source):
    sidecar = source.with_suffix(".metadata.json")
    return {"path": str(source), "sha256": sha256(source),
            "metadata_sha256": sha256(sidecar) if sidecar.exists() else None}


def validate_source_energy(metadata, requested):
    """Prevent summing retained rates with a newly computed different energy."""
    evidence = []

    def add(value, unit, origin):
        number = finite_number(value)
        scale = {"tev": 1.0, "gev": .001}.get(str(unit).strip().lower())
        if number is None or number <= 0 or scale is None:
            raise ValueError(f"invalid source collider energy at {origin}")
        value_tev = number * scale
        if not math.isclose(value_tev, requested, rel_tol=1e-12, abs_tol=0.0):
            raise ValueError(f"source collider energy at {origin} is {value_tev:g} TeV; "
                             f"requested {requested:g} TeV would mix existing and new rates")
        evidence.append({"source": origin, "energy_TeV": value_tev})

    for entry in metadata.get("fixed_parameters", []):
        if isinstance(entry, dict) and entry.get("variable") == "sqrt(s)":
            add(entry.get("value"), entry.get("unit"), "fixed_parameters.sqrt(s)")
    for name, container in (("", metadata), ("configuration.", metadata.get("configuration", {})),
                            ("config.", metadata.get("config", {}))):
        fixed = container.get("fixed", {}) if isinstance(container, dict) else {}
        if isinstance(fixed, dict) and "Energy" in fixed:
            add(fixed["Energy"], "TeV", name + "fixed.Energy")
    prior = metadata.get("mg5_enrichment", {})
    prior = prior.get("configuration", {}) if isinstance(prior, dict) else {}
    if isinstance(prior, dict) and "energy_TeV" in prior:
        add(prior["energy_TeV"], "TeV", "mg5_enrichment.configuration.energy_TeV")
    return {"status": "verified_source_metadata" if evidence else "assumed_from_requested_energy",
            "energy_TeV": requested, "evidence": evidence}


def native_rate(process, run_name, row, energy, get_xsec):
    def number(key):
        value = finite_number(row.get(key))
        if value is None:
            raise ValueError(f"missing or nonfinite stored MG5 input {key}")
        return value

    lambdas = [number(key) for key in LAMBDA_COLUMNS]
    masses = [number("M2"), number("M3")]
    widths = [number("w1"), number("w2"), number("w3")]
    if any(value <= 0 for value in masses) or any(value < 0 for value in widths):
        raise ValueError("MG5 requires positive masses and nonnegative stored widths")
    value = get_xsec(process, run_name, lambdas,
                     number("k1"), number("k2"), number("k3"),
                     masses[0], widths[1], masses[1], widths[2],
                     ecm=energy, w1=widths[0], k233=lambdas[-1])
    rate = finite_number(value)
    if rate is None or rate < 0:
        raise ValueError(f"MG5 returned invalid cross section {value!r}")
    return rate


def configuration(source, output, processes, energy, runtime_receipt):
    root = Path(__file__).resolve().parent
    files = ("reprocess_trsm_mg5.py", "trsm_mg5_rates.py", "generate_mg5_trsm_xsecs.py", "tools/madloop_build.py")
    return {"schema": "trsm_mg5_enrichment_v1", "input": input_identity(source),
            "output": str(output), "processes": list(processes), "energy_TeV": energy,
            "selection": "stored_full_viability_v1", "runtime": runtime_receipt(processes),
            "sources_sha256": {name: sha256(root / name) for name in files}}


def open_checkpoint(path, config, resume):
    if path.exists() and not resume:
        raise FileExistsError(f"checkpoint exists; repeat with --resume: {path}")
    if resume and not path.exists():
        raise FileNotFoundError(f"checkpoint not found: {path}")
    connection = sqlite3.connect(path)
    try:
        if resume:
            saved = json.loads(connection.execute("SELECT value FROM metadata WHERE key='configuration'").fetchone()[0])
            if saved != config:
                raise ValueError("checkpoint does not match input, output, code, runtime, or MG5 options")
        else:
            with connection:
                connection.execute("CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
                connection.execute("CREATE TABLE rates (row_index INTEGER, process TEXT, value REAL NOT NULL, PRIMARY KEY(row_index, process))")
                connection.execute("INSERT INTO metadata VALUES ('configuration', ?)", (json.dumps(config, sort_keys=True),))
        return connection
    except Exception:
        connection.close()
        raise


def reprocess(source, output, *, processes=ASSOCIATED_PROCESSES, energy=13.6,
              resume=False, dry_run=False, get_xsec=get_mg5_xsec,
              runtime_receipt=mg5_runtime_receipt):
    """Serialize writers to the output/checkpoint while preserving the input."""
    source, output = Path(source).expanduser().resolve(), Path(output).expanduser().resolve()
    if source == output:
        raise ValueError("input and output paths must differ")
    source_metadata = source.with_suffix(".metadata.json")
    sidecar = output.with_suffix(".metadata.json")
    auxiliary_paths = {output, sidecar, output.with_name("." + output.name + ".complete.tmp"),
                       sidecar.with_name("." + sidecar.name + ".complete.tmp")}
    auxiliary_paths.update(output.with_name(output.name + suffix) for suffix in (".partial", ".lock"))
    if {source, source_metadata.resolve()} & {path.resolve() for path in auxiliary_paths}:
        raise ValueError("input or its metadata collides with an output path")
    options = dict(processes=processes, energy=energy, resume=resume, dry_run=dry_run,
                   get_xsec=get_xsec, runtime_receipt=runtime_receipt)
    if dry_run:
        return _reprocess(source, output, **options)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.with_name(output.name + ".lock").open("a+") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError("another MG5 enrichment is writing this output") from error
        return _reprocess(source, output, **options)


def _reprocess(source, output, *, processes, energy, resume, dry_run, get_xsec,
               runtime_receipt):
    processes = tuple(dict.fromkeys(processes))
    if not processes or any(process not in ASSOCIATED_PROCESSES for process in processes):
        raise ValueError("select supported associated-production processes")
    if not math.isfinite(energy) or energy <= 0:
        raise ValueError("energy must be finite and positive")
    checkpoint = output.with_name(output.name + ".partial")
    sidecar = output.with_suffix(".metadata.json")
    source_metadata = source.with_suffix(".metadata.json")
    metadata = json.loads(source_metadata.read_text()) if source_metadata.exists() else {}
    if not isinstance(metadata, dict):
        raise ValueError("source metadata must be a JSON object")
    source_energy = validate_source_energy(metadata, energy)
    config = configuration(source, output, processes, energy, runtime_receipt)
    config["source_energy"] = source_energy
    if output.exists():
        if resume and sidecar.exists():
            saved = json.loads(sidecar.read_text())
            enrichment = saved.get("mg5_enrichment", {})
            if enrichment.get("configuration") == config and enrichment.get("output_sha256") == sha256(output):
                return enrichment["counts"]
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    if sidecar.exists():
        saved = json.loads(sidecar.read_text())
        if not resume or saved.get("mg5_enrichment", {}).get("configuration") != config:
            raise FileExistsError(f"refusing to overwrite existing metadata: {sidecar}")
    initial_rows = read_rows(source)
    header = next(initial_rows)
    counts = {"rows": 0, "eligible": 0, "requested_missing_rates": 0}
    for _, row in initial_rows:
        counts["rows"] += 1
        if stored_full_viability(row):
            counts["eligible"] += 1
            counts["requested_missing_rates"] += sum(
                finite_number(row.get(f"mg5_xsec_{process}_pb")) is None for process in processes)
    if not counts["rows"]:
        raise ValueError("input has no data rows")
    print(f"Saved scan: {counts['rows']:,} rows; {counts['eligible']:,} full viable; "
          f"{counts['requested_missing_rates']:,} missing requested rates", flush=True)
    if dry_run:
        return counts
    output.parent.mkdir(parents=True, exist_ok=True)
    connection = open_checkpoint(checkpoint, config, resume)
    namespace = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()[:20]
    temporary = output.with_name("." + output.name + ".complete.tmp")
    metadata_temp = sidecar.with_name("." + sidecar.name + ".complete.tmp")
    try:
        rows = read_rows(source)
        next(rows)
        for index, row in rows:
            if not stored_full_viability(row):
                continue
            for process in processes:
                if finite_number(row.get(f"mg5_xsec_{process}_pb")) is not None:
                    continue
                if connection.execute("SELECT 1 FROM rates WHERE row_index=? AND process=?", (index, process)).fetchone():
                    continue
                try:
                    rate = native_rate(process, f"enrich-{namespace}-row{index}", row, energy, get_xsec)
                    with connection:
                        connection.execute("INSERT INTO rates VALUES (?, ?, ?)", (index, process, rate))
                except Exception as error:
                    raise RuntimeError(f"MG5 enrichment failed at data row {index}, process {process}: {error}") from error
                print(f"Checkpointed data row {index}: {process} = {rate:.8g} pb", flush=True)
        fields = header + [key for key in (
            *(f"mg5_xsec_{process}_pb" for process in processes), *DERIVED_RATE_COLUMNS) if key not in header]
        with temporary.open("w", encoding="utf-8", newline="") as target:
            writer = csv.DictWriter(target, fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            rows = read_rows(source)
            next(rows)
            updates = iter(connection.execute("SELECT row_index, process, value FROM rates ORDER BY row_index, process"))
            current = next(updates, None)
            written = 0
            for index, row in rows:
                while current is not None and current[0] == index:
                    _, process, value = current
                    row[f"mg5_xsec_{process}_pb"] = str(value)
                    current = next(updates, None)
                for key, value in derive_mg5_rates(row).items():
                    if finite_number(row.get(key)) is None:
                        row[key] = str(value)
                writer.writerow({key: row.get(key, "nan") for key in fields})
                written += 1
            target.flush()
            os.fsync(target.fileno())
        if current is not None or written != counts["rows"] or input_identity(source) != config["input"]:
            raise ValueError("input changed during enrichment; original output was not published")
        metadata.update(scan_file=output.name, output_role="mg5_enrichment")
        metadata["mg5_enrichment"] = {"configuration": config, "counts": counts,
                                       "output_sha256": sha256(temporary)}
        with metadata_temp.open("w", encoding="utf-8") as target:
            json.dump(metadata, target, indent=2, sort_keys=True)
            target.write("\n")
            target.flush()
            os.fsync(target.fileno())
        os.replace(metadata_temp, sidecar)
        os.replace(temporary, output)
    finally:
        connection.close()
        temporary.unlink(missing_ok=True)
        metadata_temp.unlink(missing_ok=True)
    checkpoint.unlink()
    return counts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mg5-process", action="append", choices=ASSOCIATED_PROCESSES)
    parser.add_argument("--energy", type=float, default=13.6, help="Collider energy in TeV (default: 13.6)")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true", help="Validate runtime and report missing rates without running MG5")
    args = parser.parse_args(argv)
    # Saved-point enrichment is serial. Avoid inheriting an installation's
    # unrestricted core count; an explicit caller environment still wins.
    os.environ.setdefault("TRSM_MG5_CORES", "1")
    try:
        reprocess(args.input, args.output, processes=args.mg5_process or ASSOCIATED_PROCESSES,
                  energy=args.energy, resume=args.resume, dry_run=args.dry_run)
    except (OSError, ValueError, RuntimeError, sqlite3.Error) as error:
        print(f"MG5 enrichment stopped: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
