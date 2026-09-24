#!/usr/bin/env python3
"""Check the standalone DM example from any working directory."""
import argparse
import importlib.util
import subprocess
import sys

if sys.version_info < (3, 10):
    raise SystemExit('The shared TRSM modules require Python 3.10 or newer.')

from run_single_point import REPO_ROOT, add_physics_arguments, configure
import test_trsm_DM


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--no-native', action='store_true', help='Check Python imports and data without micrOMEGAs')
    add_physics_arguments(parser)
    args = parser.parse_args()
    print('Python:', sys.executable, sys.version.split()[0])
    print('Repository:', REPO_ROOT)
    print('DM module:', test_trsm_DM.__file__)
    for name in ('numpy', 'matplotlib'):
        print(f'Plot dependency {name}:', 'available' if importlib.util.find_spec(name) else 'not installed')
    try:
        configure(args, replay=args.no_native)
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f'Runtime check failed: {error}\nSee DM/README.md for installation commands.', file=sys.stderr)
        return 1
    print('Direct-detection limit:', args.si_table.model_id if args.si_table else args.limit_model)
    print('micrOMEGAs executable:', args.micromegas_main)
    print('Native capability check:', 'not requested' if args.no_native else 'passed')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
