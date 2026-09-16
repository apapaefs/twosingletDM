#!/usr/bin/env python3
"""Stop a live evo/thc-counted TRSM scan at a completed-row target."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import signal
import socket
import subprocess
import sys
import time
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


RECEIPT_SCHEMA = "trsm_scan_stop_receipt_v1"
LOCK_SCHEMA = "trsm_scan_stop_lock_v1"


class StopMonitorError(RuntimeError):
    """Raised when stopping a scan cannot be done safely."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def nonnegative_float(value: str) -> float:
    parsed = float(value)
    if parsed < 0.0:
        raise argparse.ArgumentTypeError("must be nonnegative")
    return parsed


def metadata_path(scan_path: Path) -> Path:
    return scan_path.with_suffix(".metadata.json")


def default_receipt_path(scan_path: Path, target: int) -> Path:
    return scan_path.with_name(f"{scan_path.stem}.stop-at-{target}.json")


def default_lock_path(scan_path: Path, target: int) -> Path:
    return scan_path.with_name(f"{scan_path.stem}.stop-at-{target}.lock")


def atomic_write_json(path: Path, payload: dict) -> None:
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


def load_campaign_metadata(scan_path: Path) -> tuple[Path, dict, int]:
    path = metadata_path(scan_path)
    try:
        with path.open(encoding="utf-8") as stream:
            payload = json.load(stream)
    except FileNotFoundError as error:
        raise StopMonitorError(f"Missing scan metadata: {path}") from error
    except (OSError, json.JSONDecodeError) as error:
        raise StopMonitorError(f"Cannot read scan metadata {path}: {error}") from error

    if not isinstance(payload, dict):
        raise StopMonitorError(f"Invalid scan metadata {path}: expected an object")
    options = payload.get("options")
    if not isinstance(options, dict):
        raise StopMonitorError(f"Invalid scan metadata {path}: missing options")
    if options.get("nrandom_count_evo_thc") is not True:
        raise StopMonitorError(
            "This monitor requires a scan run with --nrandom-count-evo-thc"
        )
    if options.get("write_evo_thc_points") is not True:
        raise StopMonitorError(
            "This monitor requires a scan run with --write-evo-thc-points"
        )
    seed = payload.get("seed", options.get("seed"))
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise StopMonitorError(f"Invalid or missing scan seed in {path}")
    recorded_scan = payload.get("scan_file")
    if recorded_scan and Path(recorded_scan).name != scan_path.name:
        raise StopMonitorError(
            f"Metadata describes {recorded_scan!r}, not {scan_path.name!r}"
        )
    return path, payload, seed


class EvoThcFileCounter:
    """Incrementally count complete TSV rows that explicitly pass evo and thc."""

    def __init__(self, path: Path):
        self.path = Path(path).expanduser().resolve()
        self.count = 0
        self._offset = 0
        self._pending = b""
        self._device: int | None = None
        self._inode: int | None = None
        self._column_count = 0
        self._evo_index = -1
        self._thc_index = -1
        self._open_and_load()

    @property
    def pending_bytes(self) -> int:
        return len(self._pending)

    def _open_and_load(self) -> None:
        try:
            stream = self.path.open("rb")
        except OSError as error:
            raise StopMonitorError(f"Cannot open scan output {self.path}: {error}") from error
        with stream:
            stat = os.fstat(stream.fileno())
            self._device = stat.st_dev
            self._inode = stat.st_ino
            header = stream.readline()
            if not header.endswith(b"\n"):
                raise StopMonitorError(
                    f"Scan output has no complete TSV header: {self.path}"
                )
            try:
                columns = header.rstrip(b"\r\n").decode("utf-8").split("\t")
            except UnicodeDecodeError as error:
                raise StopMonitorError(
                    f"Scan output header is not UTF-8: {self.path}"
                ) from error
            if len(columns) != len(set(columns)):
                raise StopMonitorError(f"Scan output has duplicate TSV columns: {self.path}")
            missing = [name for name in ("evo", "thc") if name not in columns]
            if missing:
                raise StopMonitorError(
                    f"Scan output lacks required columns {missing}: {self.path}"
                )
            self._column_count = len(columns)
            self._evo_index = columns.index("evo")
            self._thc_index = columns.index("thc")
            self._offset = stream.tell()
            self._consume_stream(stream)

    def _consume_line(self, line: bytes) -> None:
        line = line.rstrip(b"\r")
        if not line:
            raise StopMonitorError(
                f"Blank TSV row at data row {self.count + 1} in {self.path}"
            )
        fields = line.split(b"\t")
        if len(fields) != self._column_count:
            raise StopMonitorError(
                f"Malformed TSV data row {self.count + 1} in {self.path}: "
                f"expected {self._column_count} fields, found {len(fields)}"
            )
        evo = fields[self._evo_index]
        thc = fields[self._thc_index]
        if evo != b"True" or thc != b"True":
            raise StopMonitorError(
                f"Data row {self.count + 1} is not an evo/thc-passing point "
                f"(evo={evo!r}, thc={thc!r})"
            )
        self.count += 1

    def _consume_stream(self, stream) -> None:
        while True:
            chunk = stream.read(1024 * 1024)
            if not chunk:
                break
            data = self._pending + chunk
            lines = data.split(b"\n")
            self._pending = lines.pop()
            for line in lines:
                self._consume_line(line)
        self._offset = stream.tell()

    def update(self) -> int:
        try:
            stream = self.path.open("rb")
        except OSError as error:
            raise StopMonitorError(f"Cannot reopen scan output {self.path}: {error}") from error
        with stream:
            stat = os.fstat(stream.fileno())
            if (stat.st_dev, stat.st_ino) != (self._device, self._inode):
                raise StopMonitorError(
                    f"Scan output was replaced while being monitored: {self.path}"
                )
            if stat.st_size < self._offset:
                raise StopMonitorError(
                    f"Scan output shrank while being monitored: {self.path}"
                )
            stream.seek(self._offset)
            self._consume_stream(stream)
        return self.count


@dataclass(frozen=True)
class ProcessIdentity:
    pid: int
    pgid: int
    start_time: str
    command: str


def parse_ps_line(line: str) -> ProcessIdentity:
    fields = line.strip().split(None, 2)
    if len(fields) != 3:
        raise StopMonitorError(f"Cannot parse process record: {line!r}")
    try:
        pid = int(fields[0])
        pgid = int(fields[1])
    except ValueError as error:
        raise StopMonitorError(f"Cannot parse process IDs: {line!r}") from error
    remainder = fields[2]
    if len(remainder) < 25:
        raise StopMonitorError(f"Cannot parse process start time: {line!r}")
    start_time = remainder[:24]
    command = remainder[24:].strip()
    if not command:
        raise StopMonitorError(f"Process record has no command: {line!r}")
    return ProcessIdentity(pid, pgid, start_time, command)


def ps_environment() -> dict[str, str]:
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    return environment


def get_process_identity(pid: int) -> ProcessIdentity | None:
    result = subprocess.run(
        ["ps", "-p", str(pid), "-o", "pid=,pgid=,lstart=,command="],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=ps_environment(),
    )
    if result.returncode != 0 or not result.stdout.strip():
        return None
    return parse_ps_line(result.stdout)


def list_process_identities() -> list[ProcessIdentity]:
    result = subprocess.run(
        ["ps", "-axo", "pid=,pgid=,lstart=,command="],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=ps_environment(),
    )
    return [
        parse_ps_line(line)
        for line in result.stdout.splitlines()
        if line.strip()
    ]


def scan_process_mismatch(process: ProcessIdentity, seed: int) -> str | None:
    try:
        tokens = shlex.split(process.command)
    except ValueError as error:
        return f"cannot parse command line: {error}"
    generator_indices = [
        index
        for index, token in enumerate(tokens)
        if Path(token).name == "generate_trsm_points.py"
    ]
    if len(generator_indices) != 1:
        return "command is not a single generate_trsm_points.py process"
    index = generator_indices[0]
    if index + 1 >= len(tokens) or tokens[index + 1] != str(seed):
        return f"generator seed is not {seed}"
    required = ("--nrandom-count-evo-thc", "--write-evo-thc-points")
    missing = [option for option in required if option not in tokens]
    if missing:
        return f"generator command lacks {', '.join(missing)}"
    return None


def select_scan_process(seed: int, requested_pid: int | None) -> ProcessIdentity:
    if requested_pid is not None:
        process = get_process_identity(requested_pid)
        if process is None:
            raise StopMonitorError(f"PID {requested_pid} is not running")
        mismatch = scan_process_mismatch(process, seed)
        if mismatch:
            raise StopMonitorError(f"PID {requested_pid} is unsafe to signal: {mismatch}")
        return process

    matches = [
        process
        for process in list_process_identities()
        if scan_process_mismatch(process, seed) is None
    ]
    if not matches:
        raise StopMonitorError(
            f"No live evo/thc-counted generate_trsm_points.py process "
            f"was found for seed {seed}"
        )
    if len(matches) > 1:
        details = ", ".join(str(process.pid) for process in matches)
        raise StopMonitorError(
            f"Multiple matching processes were found ({details}); pass --pid explicitly"
        )
    return matches[0]


def same_process(original: ProcessIdentity, current: ProcessIdentity | None) -> bool:
    return current is not None and (
        current.pid,
        current.pgid,
        current.start_time,
        current.command,
    ) == (
        original.pid,
        original.pgid,
        original.start_time,
        original.command,
    )


def pid_is_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


class MonitorLock:
    def __init__(self, path: Path):
        self.path = path
        self.hostname = socket.gethostname()
        self.pid = os.getpid()
        self.token = uuid.uuid4().hex
        self.acquired = False

    def _payload(self) -> dict:
        return {
            "schema": LOCK_SCHEMA,
            "hostname": self.hostname,
            "pid": self.pid,
            "token": self.token,
            "started_utc": utc_now(),
        }

    def acquire(self) -> "MonitorLock":
        self.path.parent.mkdir(parents=True, exist_ok=True)
        for attempt in range(2):
            try:
                descriptor = os.open(
                    self.path,
                    os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                    0o644,
                )
            except FileExistsError:
                try:
                    with self.path.open(encoding="utf-8") as stream:
                        existing = json.load(stream)
                except (OSError, json.JSONDecodeError) as error:
                    raise StopMonitorError(
                        f"Cannot validate existing monitor lock {self.path}: {error}"
                    ) from error
                host = existing.get("hostname")
                pid = existing.get("pid")
                if host != self.hostname:
                    raise StopMonitorError(
                        f"Monitor lock belongs to another host ({host!r}): {self.path}"
                    )
                if isinstance(pid, int) and pid_is_alive(pid):
                    raise StopMonitorError(
                        f"A live stop monitor already exists as PID {pid}: {self.path}"
                    )
                if attempt:
                    raise StopMonitorError(f"Cannot acquire monitor lock: {self.path}")
                archived = self.path.with_name(
                    f"{self.path.name}.stale-{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}"
                )
                os.replace(self.path, archived)
                continue
            with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
                json.dump(self._payload(), stream, indent=2)
                stream.write("\n")
                stream.flush()
                os.fsync(stream.fileno())
            self.acquired = True
            return self
        raise StopMonitorError(f"Cannot acquire monitor lock: {self.path}")

    def release(self) -> None:
        if not self.acquired:
            return
        try:
            with self.path.open(encoding="utf-8") as stream:
                existing = json.load(stream)
            if existing.get("token") == self.token:
                self.path.unlink()
        except (FileNotFoundError, OSError, json.JSONDecodeError):
            pass
        self.acquired = False

    def __enter__(self) -> "MonitorLock":
        return self.acquire()

    def __exit__(self, exc_type, exc_value, traceback) -> bool:
        self.release()
        return False


def signal_scope(process: ProcessIdentity) -> str:
    if process.pgid == process.pid and process.pgid != os.getpgrp():
        return "process-group"
    return "process"


def send_stop_signal(
    process: ProcessIdentity,
    signal_number: int,
    scope: str,
) -> None:
    if scope == "process-group":
        os.killpg(process.pgid, signal_number)
    else:
        os.kill(process.pid, signal_number)


def stop_scope_is_alive(process: ProcessIdentity, scope: str) -> bool:
    try:
        if scope == "process-group":
            os.killpg(process.pgid, 0)
        else:
            os.kill(process.pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Monitor a --nrandom-count-evo-thc/--write-evo-thc-points scan "
            "and stop its verified generator process at a completed-row target."
        )
    )
    parser.add_argument("scan_file", type=Path, help="Main TRSM TSV output")
    parser.add_argument(
        "--target",
        type=positive_int,
        default=20000,
        help="Completed evo/thc point target (default: 20000)",
    )
    parser.add_argument(
        "--pid",
        type=positive_int,
        help="Expected generator PID; otherwise discover it from the saved seed",
    )
    parser.add_argument(
        "--poll-interval",
        type=nonnegative_float,
        default=0.1,
        help="Seconds between append checks (default: 0.1)",
    )
    parser.add_argument(
        "--grace-seconds",
        type=nonnegative_float,
        default=60.0,
        help="How long to wait for the process group to exit (default: 60)",
    )
    parser.add_argument(
        "--signal",
        choices=("TERM", "INT"),
        default="TERM",
        help="Signal sent at the target (default: TERM)",
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Validate the output and process, report the current count, then exit",
    )
    parser.add_argument(
        "--receipt",
        type=Path,
        help="Status JSON path (default: adjacent target-specific file)",
    )
    return parser


def run(args: argparse.Namespace) -> int:
    scan_path = args.scan_file.expanduser().resolve()
    if not scan_path.is_file():
        raise StopMonitorError(f"Scan output does not exist: {scan_path}")
    metadata_file, metadata, seed = load_campaign_metadata(scan_path)
    counter = EvoThcFileCounter(scan_path)
    process = select_scan_process(seed, args.pid)
    receipt_path = (
        args.receipt.expanduser().resolve()
        if args.receipt is not None
        else default_receipt_path(scan_path, args.target)
    )

    print(
        f"Validated seed-{seed} scan PID {process.pid}: "
        f"{counter.count}/{args.target} complete evo/thc points",
        flush=True,
    )
    if args.check_only:
        return 0

    stop_signal = signal.SIGTERM if args.signal == "TERM" else signal.SIGINT
    scope = signal_scope(process)
    receipt = {
        "schema": RECEIPT_SCHEMA,
        "status": "armed",
        "scan_file": str(scan_path),
        "metadata_file": str(metadata_file),
        "seed": seed,
        "target_evo_thc_count": args.target,
        "initial_evo_thc_count": counter.count,
        "last_observed_evo_thc_count": counter.count,
        "pending_partial_bytes": counter.pending_bytes,
        "monitor_pid": os.getpid(),
        "monitor_host": socket.gethostname(),
        "generator_pid": process.pid,
        "generator_pgid": process.pgid,
        "generator_start_time": process.start_time,
        "generator_command": process.command,
        "signal": args.signal,
        "signal_scope": scope,
        "started_utc": utc_now(),
        "metadata_requested_points": metadata.get("requested_points"),
    }
    lock = MonitorLock(default_lock_path(scan_path, args.target))
    with lock:
        atomic_write_json(receipt_path, receipt)
        last_reported = counter.count
        try:
            while counter.count < args.target:
                if not pid_is_alive(process.pid):
                    receipt.update(
                        {
                            "status": "generator_exited_before_target",
                            "last_observed_evo_thc_count": counter.update(),
                            "finished_utc": utc_now(),
                        }
                    )
                    atomic_write_json(receipt_path, receipt)
                    raise StopMonitorError(
                        f"Generator PID {process.pid} exited at "
                        f"{counter.count}/{args.target} points"
                    )
                time.sleep(args.poll_interval)
                counter.update()
                if counter.count != last_reported:
                    print(
                        f"evo/thc points: {counter.count}/{args.target}",
                        flush=True,
                    )
                    last_reported = counter.count

            current = get_process_identity(process.pid)
            if not same_process(process, current):
                receipt.update(
                    {
                        "status": "process_identity_changed",
                        "last_observed_evo_thc_count": counter.count,
                        "finished_utc": utc_now(),
                    }
                )
                atomic_write_json(receipt_path, receipt)
                raise StopMonitorError(
                    f"PID {process.pid} no longer has its validated identity; "
                    "refusing to signal it"
                )

            receipt.update(
                {
                    "status": "signal_sending",
                    "last_observed_evo_thc_count": counter.count,
                    "pending_partial_bytes": counter.pending_bytes,
                    "signal_sent_utc": utc_now(),
                }
            )
            atomic_write_json(receipt_path, receipt)
            send_stop_signal(process, stop_signal, scope)
            print(
                f"Sent SIG{args.signal} to {scope} "
                f"{process.pgid if scope == 'process-group' else process.pid} "
                f"at {counter.count} complete points",
                flush=True,
            )

            deadline = time.monotonic() + args.grace_seconds
            while stop_scope_is_alive(process, scope) and time.monotonic() < deadline:
                time.sleep(min(0.2, max(0.01, args.grace_seconds)))
            counter.update()
            alive = stop_scope_is_alive(process, scope)
            receipt.update(
                {
                    "status": "signal_sent_process_alive" if alive else "stopped",
                    "last_observed_evo_thc_count": counter.count,
                    "pending_partial_bytes": counter.pending_bytes,
                    "finished_utc": utc_now(),
                }
            )
            atomic_write_json(receipt_path, receipt)
            if alive:
                raise StopMonitorError(
                    f"The {scope} is still alive {args.grace_seconds:g} seconds "
                    "after the stop signal; no stronger signal was sent"
                )
            print(
                f"Generator stopped; final complete evo/thc count is {counter.count}. "
                f"Receipt: {receipt_path}",
                flush=True,
            )
            return 0
        except KeyboardInterrupt:
            receipt.update(
                {
                    "status": "monitor_interrupted",
                    "last_observed_evo_thc_count": counter.update(),
                    "pending_partial_bytes": counter.pending_bytes,
                    "finished_utc": utc_now(),
                }
            )
            atomic_write_json(receipt_path, receipt)
            raise


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(args)
    except StopMonitorError as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    sys.exit(main())
