#!/usr/bin/env python3
"""Repair MG5 survey iterations and loop-induced helicity sampling settings.

This maintenance command is for campaigns retaining their original scan checkout.
New scans request three iterations and exact loop helicity sums directly.
Use --repair to apply;
without it, the command only checks the selected generated process directories.
"""

import argparse
import ast
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from tools.madloop_build import repair_lock


GUARD = """        # TRSM: keep the survey maximum consistent with its minimum.
        if options['maxiter'] < options['miniter']:
            logger.warning('Survey iterations %s are below minimum %s; using the minimum.',
                           options['maxiter'], options['miniter'])
            options['maxiter'] = options['miniter']

"""


def guarded_source(source):
    """Insert after the split-grid branch, whose one-iteration chunks are valid."""
    if GUARD in source:
        return source
    if "# TRSM: keep the survey maximum" in source:
        raise ValueError("A different TRSM survey guard is already installed")
    tree = ast.parse(source)
    functions = [node for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)
                 and node.name == "write_parameter"
                 and [arg.arg for arg in node.args.args] == ["self", "parralelization", "Pdirs"]]
    if len(functions) != 1:
        raise ValueError("Unsupported MG5 survey source: expected one write_parameter method")
    function = functions[0]
    branches = [node for node in function.body if isinstance(node, ast.If)
                and isinstance(node.test, ast.UnaryOp) and isinstance(node.test.op, ast.Not)
                and isinstance(node.test.operand, ast.Name) and node.test.operand.id == "Pdirs"]
    lines = source.splitlines(keepends=True)
    body = "".join(lines[function.lineno - 1:function.end_lineno])
    if (len(branches) != 1 or "options['miniter'] = 1" not in body
            or "'miniter': self.min_iterations" not in body):
        raise ValueError("Unsupported MG5 survey source: iteration handling differs")
    insertion = branches[0].lineno - 1
    result = "".join(lines[:insertion]) + GUARD + "".join(lines[insertion:])
    compile(result, "gen_ximprove.py", "exec")
    return result


def exact_helicity_card(source):
    """Keep card formatting and all physics inputs; change only helicity sampling."""
    pattern = re.compile(r"(?m)^([ \t]*)(\S+)([ \t]*=[ \t]*"
                         r"(nhel(?:_survey|_refine)?)(?=[ \t!#]|$))")
    matches = list(pattern.finditer(source))
    names = [match[4] for match in matches]
    if (names.count("nhel") != 1 or len(names) != len(set(names))
            or any(match[2] not in ("0", "1") for match in matches)):
        raise ValueError("Unsupported MG5 run card: expected one nhel = 0 or 1 setting")
    return pattern.sub(lambda match: match[1] + "0" + match[3], source)


def pending_updates(process):
    """Validate every input before any replacement; return changed files only."""
    characteristics = (process / "SubProcesses/proc_characteristics").read_text()
    matches = re.findall(r"(?im)^\s*loop_induced\s*=\s*(true|false)\s*$", characteristics)
    if len(matches) != 1:
        raise ValueError("Unsupported MG5 process: missing or ambiguous loop_induced flag")
    changes = [("bin/internal/gen_ximprove.py", guarded_source)]
    if matches[0].lower() == "true":
        changes.append(("Cards/run_card.dat", exact_helicity_card))
        if (process / "Cards/run_card_default.dat").exists():
            changes.append(("Cards/run_card_default.dat", exact_helicity_card))
    updates = {}
    for relative, transform in changes:
        source = process / relative
        if not source.resolve().is_relative_to(process):
            raise ValueError(f"Refusing to change a source outside the selected process: {source}")
        original = source.read_bytes()
        updated = transform(original.decode("utf-8")).encode("utf-8")
        if updated != original:
            updates[relative] = (original, updated)
    return updates


def _replace_bytes(source, content):
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=source.parent, prefix=".survey-", delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        temporary.chmod(source.stat().st_mode & 0o777)
        temporary.replace(source)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def repair_survey(process):
    process = Path(process).expanduser().resolve()
    with repair_lock(process):
        updates = pending_updates(process)
        if not updates:
            return None
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        backup = process / ".trsm-maintenance" / ("survey-minimum-" + stamp)
        backup.mkdir(parents=True)
        files = []
        for relative, (original, updated) in updates.items():
            saved = backup / relative
            saved.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(process / relative, saved)
            files.append({"source": str(process / relative), "backup": str(saved),
                          "old_sha256": hashlib.sha256(original).hexdigest(),
                          "new_sha256": hashlib.sha256(updated).hexdigest()})
        receipt = {"schema": 1, "repair": "survey_minimum_and_exact_loop_helicities", "utc": stamp,
                   "files": files,
                   "tool_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        (backup / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
        for relative, (_, updated) in updates.items():
            _replace_bytes(process / relative, updated)
        return backup / "receipt.json"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("process", nargs="+", type=Path, help="Generated MG5 process directory")
    parser.add_argument("--repair", action="store_true", help="Guard survey iterations and set exact loop helicity sums")
    args = parser.parse_args(argv)
    failed = False
    for process in args.process:
        try:
            if args.repair:
                receipt = repair_survey(process)
                print(f"OK: {process}; " + (f"repair receipt: {receipt}" if receipt else "repair already installed"))
            else:
                selected = process.expanduser().resolve()
                with repair_lock(selected):
                    changes = pending_updates(selected)
                if changes:
                    raise ValueError(f"Survey repair needed: {process} ({', '.join(changes)}); "
                                     "rerun with --repair while scans are stopped")
                print(f"OK: {process}")
        except (OSError, ValueError) as error:
            print(f"ERROR: {error}", file=sys.stderr)
            failed = True
    return int(failed)


if __name__ == "__main__":
    sys.exit(main())
