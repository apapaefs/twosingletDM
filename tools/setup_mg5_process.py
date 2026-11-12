#!/usr/bin/env python3
"""Add one compiled process to an existing MG5 runtime without replacing it."""

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from generate_mg5_trsm_xsecs import MGLocation, ProcLocation, mg5_runtime_receipt
from tools.bootstrap_runtime import (
    SetupError, compile_mg5_process, process_card, sha256, write_json,
)


def configuration(root):
    result = {}
    for line in (root / "input/mg5_configuration.txt").read_text().splitlines():
        key, separator, value = line.partition("#")[0].partition("=")
        if separator:
            result[key.strip()] = value.strip()
    return result


def compiler(value, fallback):
    value = value if value and value != "None" else fallback
    resolved = shutil.which(value)
    if resolved is None:
        raise SetupError(f"Compiler not found: {value}")
    return resolved


def generation_card(spec, **settings):
    # Unlike a fresh bootstrap, this runtime may also serve saved campaigns.
    # MG5 otherwise saves configuration-changing 'set' commands automatically.
    return "\n".join(
        line + " --no_save" if line.startswith("set ") else line
        for line in process_card(spec, **settings).splitlines()
    ) + "\n"


def validate_definition(target, spec):
    card = target / "Cards/proc_card_mg5.dat"
    if not card.is_file():
        raise SetupError(f"Existing process has no generation card: {target}")
    commands = [" ".join(line.partition("#")[0].split())
                for line in card.read_text().splitlines()]
    generated = [line for line in commands if line.startswith(("generate ", "add process "))]
    imports = [line for line in commands if line.startswith("import model ")]
    if (generated != ["generate " + spec["generate"]]
            or not imports or imports[-1] != "import model loop_sm_twoscalar_generic"):
        raise SetupError(f"Refusing to replace a different process/model at {target}")


def model_hashes(model):
    files = sorted(model.rglob("*.py")) + sorted(model.glob("restrict*.dat"))
    if not files or not (model / "parameters.py").is_file():
        raise SetupError(f"Missing loop_sm_twoscalar_generic UFO: {model}")
    return {str(path.relative_to(model)): sha256(path) for path in files}


@contextmanager
def setup_lock(root):
    with (root / ".trsm-process-setup.lock").open("a+") as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise SetupError("Another MG5 process setup is running") from error
        try:
            yield
        finally:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def setup_process(root, process, *, jobs=1):
    root = Path(root).expanduser().resolve()
    launcher = root / "bin/mg5_aMC"
    if not launcher.is_file():
        raise SetupError(f"Missing MG5 launcher: {launcher}")
    spec = json.loads((ROOT / "config/runtime-sources-v2.json").read_text())["mg5_processes"][process]
    target = root / spec["directory"]
    if spec["directory"] != ProcLocation[process].rstrip("/"):
        raise SetupError("Process configuration disagrees with the scan interface")
    model = root / "models/loop_sm_twoscalar_generic"
    identity = {"process": process, "definition": spec, "model_sha256": model_hashes(model),
                "mg5_version_sha256": sha256(root / "VERSION")}
    receipt_path = target / "trsm-process-setup.json"
    with setup_lock(root):
        if target.exists():
            validate_definition(target, spec)
            if receipt_path.is_file():
                saved = json.loads(receipt_path.read_text())
                if saved["identity"] != identity:
                    raise SetupError(f"Existing process source identity differs: {target}")
                receipt = mg5_runtime_receipt([process], mgloc=root)
                if receipt["sources_sha256"] != saved["runtime"]["sources_sha256"]:
                    raise SetupError(f"Existing generated process sources changed: {target}")
            else:
                exported = target / "bin/internal/ufomodel"
                # A matching card alone cannot establish which UFO was used.
                for name, digest in identity["model_sha256"].items():
                    if name.endswith(".py") and (
                        not (exported / name).is_file() or sha256(exported / name) != digest
                    ):
                        raise SetupError(f"Existing process UFO differs from the runtime: {target}")
                receipt = mg5_runtime_receipt([process], mgloc=root)
            print(f"Validated existing process: {target}", flush=True)
            return receipt

        settings = configuration(root)
        fc = compiler(os.environ.get("FC") or settings.get("fortran_compiler"), "gfortran")
        cxx = compiler(os.environ.get("CXX") or settings.get("cpp_compiler"), "c++")
        collier = settings.get("collier")
        if collier and collier not in ("None", "auto"):
            collier = Path(collier)
            collier = collier if collier.is_absolute() else root / collier
        else:
            collier = next((path for path in (root.parent / "collier", root / "HEPTools/collier",
                                             root / "HEPTools/lib")
                            if (path / "libcollier.a").is_file()), None)
        if collier is None or not (collier / "libcollier.a").is_file():
            raise SetupError("No compiled COLLIER library found in this runtime")

        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        logs = root / "trsm-process-setup-logs"
        logs.mkdir(exist_ok=True)
        card = logs / f"{process}-{stamp}.mg5"
        card.write_text(generation_card(spec, jobs=jobs, fc=fc, cxx=cxx, collier=collier))
        log_path = card.with_suffix(".log")
        print(f"Generating and compiling {process}; log: {log_path}", flush=True)
        env = dict(os.environ)
        for key in ("MAKEFLAGS", "GNUMAKEFLAGS", "MFLAGS", "MAKELEVEL", "MAKEOVERRIDES"):
            env.pop(key, None)
        with log_path.open("x") as log:
            def run(command, cwd=None):
                command = [str(part) for part in command]
                log.write("\n$ " + shlex.join(command) + "\n")
                log.flush()
                result = subprocess.run(command, cwd=cwd or root, env=env,
                                        stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT)
                if result.returncode:
                    raise SetupError(f"Command failed ({result.returncode}); inspect {log_path}")

            run([sys.executable, launcher, card], cwd=root)
            if not (target / "bin/madevent").is_file():
                raise SetupError(f"MG5 did not generate {process}; inspect {log_path}")
            validate_definition(target, spec)
            compile_mg5_process(target, sys.executable, run)
        receipt = mg5_runtime_receipt([process], mgloc=root)
        write_json(receipt_path, {"identity": identity, "runtime": receipt,
                                 "created_utc": datetime.now(timezone.utc).isoformat(),
                                 "compilers": {"FC": fc, "CXX": cxx}, "collier": str(collier),
                                 "generation_card": str(card), "build_log": str(log_path)})
        print(f"Process ready: {target}", flush=True)
        return receipt


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mg5-process", required=True, choices=tuple(ProcLocation))
    parser.add_argument("--mg5-location", type=Path, default=Path(MGLocation),
                        help="Defaults to TRSM_MG5_LOCATION or the usual sibling MG5 installation")
    parser.add_argument("--jobs", type=int, default=1)
    args = parser.parse_args(argv)
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    try:
        setup_process(args.mg5_location, args.mg5_process, jobs=args.jobs)
    except (SetupError, OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"Process setup stopped: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
