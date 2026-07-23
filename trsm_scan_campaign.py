"""Crash-safe state management for resumable TRSM scan campaigns."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import socket
import uuid
from datetime import datetime, timezone
from pathlib import Path


CHECKPOINT_SCHEMA = "trsm_scan_checkpoint_v1"
SAMPLING_ALGORITHM_VERSION = "trsm_random_vxzero_v1"


class CampaignStateError(RuntimeError):
    """Raised when saved campaign state is unsafe or inconsistent."""


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def timestamp_token():
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def checkpoint_path(scan_path):
    return Path(scan_path).with_suffix(".checkpoint.json")


def lock_path(scan_path):
    return Path(scan_path).with_suffix(".lock")


def atomic_write_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("w", encoding="utf-8") as stream:
            json.dump(payload, stream, indent=2, sort_keys=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def load_json(path, description):
    path = Path(path)
    try:
        with path.open(encoding="utf-8") as stream:
            payload = json.load(stream)
    except FileNotFoundError as error:
        raise CampaignStateError(f"Missing {description}: {path}") from error
    except (OSError, json.JSONDecodeError) as error:
        raise CampaignStateError(f"Cannot read {description} {path}: {error}") from error
    if not isinstance(payload, dict):
        raise CampaignStateError(f"Invalid {description} {path}: expected a JSON object")
    return payload


def encode_rng_state(state):
    if isinstance(state, tuple):
        return [encode_rng_state(item) for item in state]
    if isinstance(state, list):
        return [encode_rng_state(item) for item in state]
    if state is None or isinstance(state, (int, float, str, bool)):
        return state
    raise CampaignStateError(f"Unsupported RNG-state value: {type(state).__name__}")


def decode_rng_state(state):
    if isinstance(state, list):
        return tuple(decode_rng_state(item) for item in state)
    if state is None or isinstance(state, (int, float, str, bool)):
        return state
    raise CampaignStateError(f"Invalid JSON RNG-state value: {type(state).__name__}")


def configuration_fingerprint(configuration):
    canonical = json.dumps(
        configuration,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("ascii")
    return hashlib.sha256(canonical).hexdigest()


def _file_digest(path, start, length):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        stream.seek(start)
        remaining = length
        while remaining:
            chunk = stream.read(min(remaining, 1024 * 1024))
            if not chunk:
                break
            digest.update(chunk)
            remaining -= len(chunk)
    return digest.hexdigest()


def inspect_tsv(path):
    """Return a compact integrity record for an append-only TSV output."""
    path = Path(path)
    if not path.exists():
        return {
            "path": str(path.resolve()),
            "exists": False,
            "size": 0,
            "data_rows": 0,
            "ends_with_newline": True,
            "head_sha256": None,
            "tail_sha256": None,
        }

    size = path.stat().st_size
    if size == 0:
        return {
            "path": str(path.resolve()),
            "exists": True,
            "size": 0,
            "data_rows": 0,
            "ends_with_newline": True,
            "head_sha256": None,
            "tail_sha256": None,
        }

    newline_count = 0
    final_byte = b""
    with path.open("rb") as stream:
        while True:
            chunk = stream.read(1024 * 1024)
            if not chunk:
                break
            newline_count += chunk.count(b"\n")
            final_byte = chunk[-1:]
    sample = min(size, 4096)
    return {
        "path": str(path.resolve()),
        "exists": True,
        "size": size,
        "data_rows": max(0, newline_count - 1),
        "ends_with_newline": final_byte == b"\n",
        "head_sha256": _file_digest(path, 0, sample),
        "tail_sha256": _file_digest(path, size - sample, sample),
    }


def output_records(output_paths):
    return {
        role: inspect_tsv(path)
        for role, path in sorted(output_paths.items())
    }


def advance_output_records(output_paths, previous_records):
    """Update integrity records by reading only bytes appended since the last draw."""
    updated = {}
    for role, path in sorted(output_paths.items()):
        path = Path(path)
        previous = previous_records.get(role)
        if previous is None:
            raise CampaignStateError(f"No previous output record for role {role!r}")
        old_size = int(previous["size"])
        if not path.exists():
            if old_size:
                raise CampaignStateError(f"Output disappeared during scan: {path}")
            updated[role] = inspect_tsv(path)
            continue
        new_size = path.stat().st_size
        if new_size < old_size:
            raise CampaignStateError(f"Output shrank during scan: {path}")
        appended_newlines = 0
        final_byte = b""
        if new_size > old_size:
            with path.open("rb") as stream:
                stream.seek(old_size)
                while True:
                    chunk = stream.read(1024 * 1024)
                    if not chunk:
                        break
                    appended_newlines += chunk.count(b"\n")
                    final_byte = chunk[-1:]
            if final_byte != b"\n":
                raise CampaignStateError(
                    f"Output append did not finish with a newline: {path}"
                )
        sample = min(new_size, 4096)
        data_rows = int(previous["data_rows"]) + appended_newlines
        if old_size == 0 and new_size > 0:
            data_rows -= 1
        updated[role] = {
            "path": str(path.resolve()),
            "exists": True,
            "size": new_size,
            "data_rows": max(0, data_rows),
            "ends_with_newline": (
                previous["ends_with_newline"] if new_size == old_size else True
            ),
            "head_sha256": (
                _file_digest(path, 0, sample) if sample else None
            ),
            "tail_sha256": (
                _file_digest(path, new_size - sample, sample) if sample else None
            ),
        }
    return updated


def backup_file(path, reason):
    path = Path(path)
    if not path.exists():
        return None
    base = path.with_name(f"{path.name}.{reason}-{timestamp_token()}")
    backup = base
    suffix = 1
    while backup.exists():
        backup = base.with_name(f"{base.name}-{suffix}")
        suffix += 1
    shutil.copy2(path, backup)
    return backup


def repair_incomplete_tsv_tail(path):
    """Back up and remove a non-newline-terminated legacy TSV tail."""
    path = Path(path)
    record = inspect_tsv(path)
    if record["size"] == 0 or record["ends_with_newline"]:
        return None
    backup = backup_file(path, "incomplete-tail-backup")
    with path.open("r+b") as stream:
        stream.seek(0, os.SEEK_END)
        position = stream.tell() - 1
        while position >= 0:
            stream.seek(position)
            if stream.read(1) == b"\n":
                stream.truncate(position + 1)
                stream.flush()
                os.fsync(stream.fileno())
                return backup
            position -= 1
        stream.truncate(0)
        stream.flush()
        os.fsync(stream.fileno())
    return backup


def reconcile_outputs(output_paths, expected_records):
    """Roll append-only files back to the last checkpoint after making backups."""
    backups = []
    for role, path in sorted(output_paths.items()):
        path = Path(path)
        if role not in expected_records:
            raise CampaignStateError(
                f"Checkpoint does not describe the {role!r} output"
            )
        expected = expected_records[role]
        if expected.get("path") != str(path.resolve()):
            raise CampaignStateError(
                f"Checkpoint {role} output path does not match {path.resolve()}"
            )
        current = inspect_tsv(path)
        expected_size = int(expected.get("size", -1))
        if expected_size < 0:
            raise CampaignStateError(
                f"Checkpoint has an invalid size for the {role!r} output"
            )
        if current["size"] < expected_size:
            raise CampaignStateError(
                f"{role} output is shorter than its checkpoint "
                f"({current['size']} < {expected_size} bytes): {path}"
            )
        if current["size"] > expected_size:
            backup = backup_file(path, "uncheckpointed-tail-backup")
            if backup is not None:
                backups.append(backup)
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("r+b") as stream:
                stream.truncate(expected_size)
                stream.flush()
                os.fsync(stream.fileno())

        restored = inspect_tsv(path)
        if restored["size"] != expected_size:
            raise CampaignStateError(f"Could not restore {role} output to checkpoint")
        for key in (
            "exists",
            "data_rows",
            "ends_with_newline",
            "head_sha256",
            "tail_sha256",
        ):
            if restored.get(key) != expected.get(key):
                raise CampaignStateError(
                    f"{role} output does not match checkpoint field {key}: {path}"
                )
    return backups


def quarantine_uncheckpointed_point_dirs(workdir, last_completed_draw):
    workdir = Path(workdir)
    if not workdir.exists():
        return []
    quarantined = []
    stamp = timestamp_token()
    for child in sorted(workdir.iterdir()):
        if not child.is_dir() or not child.name.startswith("point_"):
            continue
        index_text = child.name.removeprefix("point_")
        if not index_text.isdigit() or int(index_text) <= int(last_completed_draw):
            continue
        destination = child.with_name(f"{child.name}.interrupted-{stamp}")
        suffix = 1
        while destination.exists():
            destination = child.with_name(
                f"{child.name}.interrupted-{stamp}-{suffix}"
            )
            suffix += 1
        child.rename(destination)
        quarantined.append(destination)
    return quarantined


def _pid_is_alive(pid):
    try:
        os.kill(int(pid), 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


class CampaignLock:
    """Atomic same-host writer lock with conservative stale-lock handling."""

    def __init__(self, scan_path, command_line):
        self.path = lock_path(scan_path)
        self.command_line = list(command_line)
        self.hostname = socket.gethostname()
        self.pid = os.getpid()
        self.token = uuid.uuid4().hex
        self.acquired = False
        self.archived_stale_lock = None

    def _payload(self):
        return {
            "schema": "trsm_scan_lock_v1",
            "hostname": self.hostname,
            "pid": self.pid,
            "token": self.token,
            "started_utc": utc_now(),
            "command_line": self.command_line,
        }

    def _archive_stale_same_host_lock(self):
        existing = load_json(self.path, "campaign lock")
        host = existing.get("hostname")
        pid = existing.get("pid")
        if host != self.hostname:
            raise CampaignStateError(
                f"Campaign is locked by host {host!r}; refusing to assume it is stale: "
                f"{self.path}"
            )
        if not isinstance(pid, int):
            raise CampaignStateError(
                f"Campaign lock has no valid PID; refusing to remove it: {self.path}"
            )
        if _pid_is_alive(pid):
            raise CampaignStateError(
                f"Campaign already has a live writer (PID {pid} on {host}): "
                f"{self.path}"
            )
        archived = self.path.with_name(
            f"{self.path.name}.stale-{timestamp_token()}"
        )
        suffix = 1
        while archived.exists():
            archived = self.path.with_name(
                f"{self.path.name}.stale-{timestamp_token()}-{suffix}"
            )
            suffix += 1
        os.replace(self.path, archived)
        self.archived_stale_lock = archived

    def acquire(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        for attempt in range(2):
            try:
                descriptor = os.open(
                    self.path,
                    os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                    0o644,
                )
            except FileExistsError:
                if attempt:
                    raise CampaignStateError(
                        f"Could not acquire campaign lock: {self.path}"
                    )
                self._archive_stale_same_host_lock()
                continue
            with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
                json.dump(self._payload(), stream, indent=2, sort_keys=False)
                stream.write("\n")
                stream.flush()
                os.fsync(stream.fileno())
            self.acquired = True
            return self
        raise CampaignStateError(f"Could not acquire campaign lock: {self.path}")

    def release(self):
        if not self.acquired:
            return
        try:
            existing = load_json(self.path, "campaign lock")
            if existing.get("token") == self.token:
                self.path.unlink()
        except CampaignStateError:
            pass
        finally:
            self.acquired = False

    def __enter__(self):
        return self.acquire()

    def __exit__(self, exc_type, exc_value, traceback):
        self.release()
        return False
