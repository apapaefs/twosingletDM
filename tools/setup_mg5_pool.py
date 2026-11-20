#!/usr/bin/env python3
"""Create independent, compiled MadEvent copies while the source is idle."""

import argparse
from contextlib import ExitStack
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from generate_mg5_trsm_xsecs import MGLocation, ProcLocation, mg5_runtime_receipt
from trsm_mg5_pool import SCHEMA, MANIFEST


def clone_process(source, destination):
    """Copy writable files; preserve internal links and external library links."""
    source, destination = Path(source).resolve(), Path(destination).resolve()

    def ignore(directory, names):
        top = Path(directory) == source
        excluded = set()
        for name in names:
            if (name == "__pycache__" or name.endswith((".pyc", ".dcmd"))
                    or name.startswith(".trsm") or name == "ME5_debug"
                    or (top and name in {"Events", "HTML", "madevent.tar.gz"})
                    or re.fullmatch(r"G\d+(?:_\d+)?", name)):
                excluded.add(name)
        return excluded

    shutil.copytree(source, destination, symlinks=True, ignore=ignore)
    for folder in ("Events", "HTML"):
        (destination / folder).mkdir(exist_ok=True)
    for directory, dirs, files in os.walk(destination, followlinks=False):
        for name in dirs + files:
            link = Path(directory) / name
            if not link.is_symlink():
                continue
            original = source / link.relative_to(destination)
            target = Path(os.path.abspath(original.parent / os.readlink(original)))
            if not target.is_relative_to(source) and target.resolve().is_relative_to(source):
                target = target.resolve()
            if target.is_relative_to(source):
                mapped = destination / target.relative_to(source)
                new_target = os.path.relpath(mapped, link.parent)
            else:
                if not original.relative_to(source).parts[0] == "lib":
                    raise ValueError(f"External non-library link cannot be pooled: {original}")
                if not target.exists():
                    raise ValueError(f"Missing external library: {target}")
                new_target = str(target)
            link.unlink()
            link.symlink_to(new_target)


def setup_pool(source_root, destination, processes, slots):
    source_root = Path(source_root).expanduser().resolve()
    destination = Path(destination).expanduser().resolve()
    if slots < 1:
        raise ValueError("The pool needs at least one slot per process")
    if destination.exists():
        raise FileExistsError(f"Pool destination already exists: {destination}")
    if destination.is_relative_to(source_root):
        raise ValueError("Put the pool outside the source MG5 installation")
    processes = list(dict.fromkeys(processes))
    if not processes or any(p not in ProcLocation for p in processes):
        raise ValueError("Select supported MG5 processes")
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".mg5-pool-building-", dir=destination.parent))
    try:
        template = staging / "template"
        template.mkdir()
        with ExitStack() as stack:
            for process in processes:
                source = source_root / ProcLocation[process]
                stream = stack.enter_context((source / ".trsm_madevent.lock").open("a+"))
                try:
                    fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError as error:
                    raise RuntimeError(f"Source MG5 process is active: {source}") from error
            shutil.copy2(source_root / "VERSION", template / "VERSION")
            for process in processes:
                clone_process(source_root / ProcLocation[process], template / ProcLocation[process])
        receipt = mg5_runtime_receipt(processes, template, include_pool=False)
        # The template is an exact copy of the source's compiled processes.
        receipt["location"] = str(source_root)
        for process in processes:
            receipt["processes"][process]["directory"] = str(source_root / ProcLocation[process].rstrip("/"))
        locations = [f"slots/{i:03d}" for i in range(slots)]
        for slot in locations:
            root = staging / slot
            root.mkdir(parents=True)
            shutil.copy2(template / "VERSION", root / "VERSION")
            for process in processes:
                clone_process(template / ProcLocation[process], root / ProcLocation[process])
            copied = mg5_runtime_receipt(processes, root, include_pool=False)
            if copied["sources_sha256"] != receipt["sources_sha256"]:
                raise ValueError(f"Source/model hash mismatch in {slot}")
            print(f"Prepared {slot} ({len(processes)} processes)", flush=True)
        data = {"schema": SCHEMA, "created_utc": datetime.now(timezone.utc).isoformat(),
                "source_root": str(source_root), "source_receipt": receipt,
                "processes": processes, "slots": locations,
                "copy_policy": "private regular files, internal relative links, external library links; no historical events or integration directories"}
        (staging / MANIFEST).write_text(json.dumps(data, indent=2) + "\n")
        staging.rename(destination)
        return data
    except BaseException:
        shutil.rmtree(staging)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(MGLocation))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--slots", type=int, default=16)
    parser.add_argument("--process", action="append", choices=tuple(ProcLocation))
    args = parser.parse_args()
    data = setup_pool(args.source, args.output,
                      args.process or ["gg_heta0", "pp_eta0Z", "gg_eta0Z"], args.slots)
    print(json.dumps({"pool": str(args.output.resolve()), "slots_per_process": len(data["slots"]),
                      "processes": data["processes"]}, indent=2))


if __name__ == "__main__":
    main()
