"""Locate shared DM modules independently of the shell's working directory."""
import os
from pathlib import Path
import sys


def repository_root():
    if sys.version_info < (3, 10):
        raise SystemExit('The shared TRSM modules require Python 3.10 or newer.')
    explicit = os.environ.get('TRSM_REPO_ROOT')
    candidates = [Path(explicit).expanduser().resolve()] if explicit else Path(__file__).resolve().parents
    for root in candidates:
        if (root / 'test_trsm_DM.py').is_file() and (root / 'trsm_inputs.py').is_file():
            sys.path.insert(0, str(root))
            return root
    raise SystemExit(
        'Cannot locate the shared TRSM Python modules (test_trsm_DM.py, trsm_inputs.py).\n'
        f'Example script: {Path(__file__).resolve()}\n'
        f'Python: {sys.executable}\n'
        'Run these scripts from a complete twosingletDM checkout, or set '
        'TRSM_REPO_ROOT to that checkout. No pip package named test_trsm_DM is needed.\n'
        'See DM/steer_example/README.md, Troubleshooting.'
    )
