"""Independent MadEvent workspaces, leased separately for each process.

The pool changes scheduling only. Input rounding, cards, integration settings,
and the legacy per-directory MadEvent lock are owned by the existing driver.
"""

from contextlib import contextmanager
from functools import lru_cache
import fcntl
import hashlib
import json
import os
from pathlib import Path
import time


SCHEMA = "trsm_mg5_pool_v1"
MANIFEST = "pool.json"


@lru_cache(maxsize=8)
def _manifest(root, mtime_ns):
    path = Path(root) / MANIFEST
    data = json.loads(path.read_text())
    if data.get("schema") != SCHEMA:
        raise ValueError(f"Unsupported MG5 pool manifest: {path}")
    slots = data.get("slots")
    if not isinstance(slots, list) or not slots or len(slots) != len(set(slots)):
        raise ValueError(f"Invalid MG5 pool slots: {path}")
    for slot in slots:
        if not isinstance(slot, str) or Path(slot).is_absolute():
            raise ValueError(f"Invalid MG5 pool slot path: {slot!r}")
        location = (Path(root) / slot).resolve()
        if not location.is_relative_to(Path(root)) or location == Path(root):
            raise ValueError(f"MG5 pool slot leaves the pool: {slot!r}")
    return data


def configured_pool(source_root):
    value = os.environ.get("TRSM_MG5_POOL_ROOT")
    if not value:
        return None
    root = Path(value).expanduser().resolve()
    path = root / MANIFEST
    data = _manifest(str(root), path.stat().st_mtime_ns)
    if Path(data["source_root"]).resolve() != Path(source_root).expanduser().resolve():
        raise ValueError("MG5 pool was prepared from a different source installation")
    return root, data


def validate_pool(source_receipt, process_locations, receipt_function):
    """Verify every copy against the unchanged reference process/model sources."""
    configured = configured_pool(source_receipt["location"])
    if configured is None:
        return None
    root, data = configured
    processes = list(process_locations)
    if not set(processes).issubset(data["processes"]):
        raise ValueError("MG5 pool does not contain every requested process")
    expected = source_receipt["sources_sha256"]
    recorded = data["source_receipt"]["sources_sha256"]
    if any(recorded.get(name) != digest for name, digest in expected.items()):
        raise ValueError("MG5 source installation differs from the pool's recorded source")
    for slot in data["slots"]:
        receipt = receipt_function(processes, root / slot, include_pool=False)
        if receipt["sources_sha256"] != expected:
            raise ValueError(f"MG5 pool source/model mismatch in {slot}")
    return {"root": str(root), "slots_per_process": len(data["slots"]),
            "processes": processes, "manifest_sha256": hashlib.sha256(
                (root / MANIFEST).read_bytes()).hexdigest(),
            "locking": "exclusive_per_process_slot"}


@contextmanager
def process_workspace(process, source_root, process_locations):
    configured = configured_pool(source_root)
    if configured is None:
        yield source_root
        return
    root, data = configured
    if process not in data["processes"]:
        raise ValueError(f"Process {process!r} is absent from the MG5 pool")
    slots = data["slots"]
    # Do not consume the scan's random-number stream when scheduling jobs.
    offset = os.getpid() % len(slots)
    slots = slots[offset:] + slots[:offset]
    started = time.monotonic()
    while True:
        for slot in slots:
            location = root / slot
            directory = location / process_locations[process]
            if not (directory / "bin/madevent").is_file():
                raise FileNotFoundError(f"Incomplete MG5 pool process: {directory}")
            with (directory / ".trsm_pool_slot.lock").open("a+") as stream:
                try:
                    fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError:
                    continue
                acquired = time.monotonic()
                print(f"MG5_POOL acquired process={process} slot={slot} "
                      f"wait_seconds={acquired-started:.3f}", flush=True)
                try:
                    yield location
                finally:
                    fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
                    print(f"MG5_POOL released process={process} slot={slot} "
                          f"elapsed_seconds={time.monotonic()-acquired:.3f}", flush=True)
                return
        time.sleep(0.1 + (os.getpid() % 17) * 0.005)
