#!/usr/bin/env python3
"""Check or repair the compilation artifacts needed by MadLoop initialization.

Usage: python tools/madloop_build.py [--repair] /path/to/MG5/process [...]
No cards, physics sources, events, or campaign checkpoints are changed.
"""

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys


class MadLoopBuildError(ValueError):
    pass


def madloop_subprocesses(process):
    process = Path(process).expanduser().resolve()
    root = process / "SubProcesses"
    if not root.is_dir():
        raise MadLoopBuildError(f"No generated SubProcesses directory: {root}")
    paths = []
    for path in sorted(root.glob("PV*")):
        if not path.is_dir():
            continue
        if path.resolve().parent != root.resolve():
            raise MadLoopBuildError(f"MadLoop directory leaves the selected process: {path}")
        for name in ("makefile", "loop_matrix.f", "check_sa.f"):
            if not (path / name).is_file():
                raise MadLoopBuildError(f"Incomplete MadLoop source directory: {path / name}")
        paths.append(path)
    characteristics = root / "proc_characteristics"
    if (not paths and characteristics.is_file()
            and re.search(r"(?im)^\s*loop_induced\s*=\s*true\s*$", characteristics.read_text())):
        raise MadLoopBuildError(f"Loop-induced process has no MadLoop subprocesses: {process}")
    return paths


def polynomial_modules(directory):
    """Get module names from the generated source, without assuming an ML5 prefix."""
    source = Path(directory) / "polynomial.f"
    if not source.is_file():
        # Non-optimized MadLoop exports do not have polynomial routines.
        return []
    names = re.findall(r"(?im)^\s*module\s+(?!procedure\b)([a-z][a-z0-9_]*)\s*$",
                       source.read_text())
    if not names:
        raise MadLoopBuildError(f"No Fortran module declaration in {source}")
    return [source.with_name(name.lower() + ".mod") for name in dict.fromkeys(names)]


def _nonempty(path):
    return path.is_file() and path.stat().st_size > 0


def validate_madloop(process):
    """Read-only preflight; a parent-directory module alone is insufficient."""
    artifacts, missing = [], []
    for directory in madloop_subprocesses(process):
        check = directory / "check"
        modules = polynomial_modules(directory)
        artifacts.extend([check, *modules])
        if not _nonempty(check) or not os.access(check, os.X_OK):
            missing.append(check)
        missing.extend(path for path in modules if not _nonempty(path))
    if missing:
        command = shlex.join(["python", "tools/madloop_build.py", "--repair",
                              str(Path(process).expanduser().resolve())])
        raise MadLoopBuildError(
            "Missing MadLoop initialization artifacts: " + ", ".join(map(str, missing))
            + f". With scans using this process stopped, run: {command}")
    return artifacts


def build_madloop_checks(process, run):
    """Compile each check target after the integration binary/OLP library.

    The caller must own a fresh installation or hold the MadEvent process lock.
    `run` executes a command list and raises on a nonzero return code.
    """
    for directory in madloop_subprocesses(process):
        command = ["make", "-j1", "-C", str(directory)]
        if any(not _nonempty(path) for path in polynomial_modules(directory)):
            # OLP_static builds from SubProcesses/, leaving .mod there but .o
            # in PV*/. Make does not track .mod as an output: force the local
            # compilation even when polynomial.o is already up to date.
            command.extend(["-W", "polynomial.f"])
        command.append("check")
        run(command)
    return validate_madloop(process)


@contextmanager
def repair_lock(process):
    """Use the same lock as generate_mg5_trsm_xsecs._madevent_lock."""
    with (process / ".trsm_madevent.lock").open("a+") as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise MadLoopBuildError(f"MG5 process is in use; repair was not started: {process}") from error
        try:
            yield
        finally:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def repair_madloop(process):
    process = Path(process).expanduser().resolve()
    # Validate the selection before creating a lock or log.
    if not madloop_subprocesses(process):
        return []
    with repair_lock(process):
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        log_path = process / f"madloop-build-{stamp}.log"
        env = dict(os.environ)
        for key in ("MAKEFLAGS", "GNUMAKEFLAGS", "MFLAGS", "MAKELEVEL", "MAKEOVERRIDES"):
            env.pop(key, None)
        with log_path.open("x") as log:
            print(f"Building MadLoop checks for {process}; log: {log_path}", flush=True)

            def run(command):
                log.write("\n$ " + shlex.join(command) + "\n")
                log.flush()
                result = subprocess.run(command, cwd=process, env=env,
                                        stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT)
                if result.returncode:
                    raise MadLoopBuildError(
                        f"MadLoop compilation failed (exit {result.returncode}); inspect {log_path}")

            return build_madloop_checks(process, run)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("process", type=Path, nargs="+", help="Generated MG5 process directory")
    parser.add_argument("--repair", action="store_true", help="Build missing artifacts under the MG5 writer lock")
    args = parser.parse_args(argv)
    failed = False
    for process in args.process:
        try:
            artifacts = repair_madloop(process) if args.repair else validate_madloop(process)
            print(f"OK: {process} ({len(artifacts)} MadLoop initialization artifacts)")
        except (OSError, MadLoopBuildError) as error:
            print(f"ERROR: {error}", file=sys.stderr)
            failed = True
    return int(failed)


if __name__ == "__main__":
    sys.exit(main())
