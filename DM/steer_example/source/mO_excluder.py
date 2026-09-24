#!/usr/bin/env python3
"""Compatibility entry point: assess complete raw logs using the v2 provider.

Example: mO_excluder.py --card MO_inp1.dat --micromegas-output OUT_mO_1
The original positional interface retains its historical relic/DD/gamma policy.
Use complete logs for v2 assessment including solver status and optional CMB.
"""
from pathlib import Path
import sys
import os
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_single_point import main

if __name__ == '__main__':
    if len(sys.argv) == 2 and sys.argv[1] in ('--plot-dirdet-limits', '--plot-indirect-limits'):
        from plot.plot_limits import main as plot_limits
        output = os.environ.get('MO_EXCLUDER_OUTPUT_DIR', str(Path(__file__).resolve().parents[1] / 'output'))
        kind = 'direct' if sys.argv[1] == '--plot-dirdet-limits' else 'indirect'
        sys.exit(plot_limits([kind, '--output-dir', output]))
    if len(sys.argv) > 1 and not sys.argv[1].startswith('-'):
        from legacy_excluder import main as legacy_main
        try:
            sys.exit(legacy_main(sys.argv[1:]))
        except (OSError, ValueError) as error:
            sys.exit(str(error))
    sys.exit(main())
