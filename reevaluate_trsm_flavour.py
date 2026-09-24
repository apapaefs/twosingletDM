#!/usr/bin/env python3
"""Append current flavour assessments to a separate copy of a saved vx=0 scan."""

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import tempfile

from scan_output import format_output_value
from trsm_flavour import FLAVOUR_COLUMNS, flavour_configuration, flavour_updates_from_row
from trsm_scan_campaign import atomic_write_json


def file_hash(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def reevaluate(source, output):
    source, output = Path(source).resolve(), Path(output).resolve()
    sidecar = output.with_suffix(".metadata.json")
    if source == output:
        raise ValueError("input and output must be different files")
    if output.exists() or sidecar.exists():
        raise FileExistsError(f"refusing to overwrite output or metadata: {output}")
    configuration = flavour_configuration()
    identity = file_hash(source)
    metadata = {}
    source_metadata = source.with_suffix(".metadata.json")
    if source_metadata.is_file():
        metadata = json.loads(source_metadata.read_text())
        if not isinstance(metadata, dict):
            raise ValueError("input metadata must be a JSON object")
    output.parent.mkdir(parents=True, exist_ok=True)
    counts = Counter()
    temporary = None
    try:
        with source.open(encoding="ascii", newline="") as stream, tempfile.NamedTemporaryFile(
                mode="w", encoding="ascii", newline="", dir=output.parent,
                prefix=f".{output.name}.", suffix=".tmp", delete=False) as target:
            temporary = Path(target.name)
            reader = csv.DictReader(stream, delimiter="\t")
            header = reader.fieldnames
            if not header or any(not key for key in header) or len(header) != len(set(header)):
                raise ValueError("input needs a nonempty, unique TSV header")
            fields = header + [key for key in FLAVOUR_COLUMNS if key not in header]
            writer = csv.DictWriter(target, fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row_number, row in enumerate(reader, start=2):
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(f"malformed TSV row {row_number}")
                updates = flavour_updates_from_row(row)
                counts[updates["flavour_status"]] += 1
                row.update({key: format_output_value(value) for key, value in updates.items()})
                writer.writerow(row)
            target.flush()
            os.fsync(target.fileno())
        if file_hash(source) != identity:
            raise ValueError("input changed during reevaluation; use a completed scan or a snapshot")
        if metadata.get("schema") not in ("trsm_scan_metadata_v1", "trsm_scan_metadata_v2"):
            metadata = {"source_metadata": metadata}
        metadata.setdefault("schema", "trsm_scan_metadata_v1")
        metadata.setdefault("variable_ranges", [])
        metadata.update(scan_file=output.name, output_role="flavour_reevaluation")
        metadata["flavour"] = configuration
        metadata["flavour_reevaluation"] = {
            "source": str(source), "source_sha256": identity,
            "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "rows": sum(counts.values()), "status_counts": dict(counts),
            "decay_inputs": "reconstructed with current Python decay and canonical portal helpers",
            "other_assessments": "retained from source; not reevaluated",
        }
        atomic_write_json(sidecar, metadata)
        os.replace(temporary, output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return dict(counts)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    counts = reevaluate(args.input, args.output)
    print(f"Saved {args.output}: {sum(counts.values())} rows; {counts}")


if __name__ == "__main__":
    main()
