#!/usr/bin/env python3
"""Reproduce four native CMB benchmarks and check the exposed comparison arithmetic."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

from run_single_point import sha256, write_json, write_tsv
from trsm_micromegas import normalize_micromegas_version

BASE = Path(__file__).resolve().parent
QUANTITIES = ('dm_omega', 'dm_xf', 'dm_freezeout_temperature_GeV', 'dm_dir_det',
              'dm_cmb_ratio_raw', 'dm_cmb_abundance_fraction', 'dm_cmb_ratio')


def validate(output, version, *, rtol=2e-6):
    rows = json.loads((output / 'results.json').read_text())
    reference = json.loads((BASE / 'benchmarks/cmb-reference-v2.json').read_text())
    checks = []
    def check(label, condition):
        checks.append({'check': label, 'passed': bool(condition)})
    check('benchmark input checksum', sha256(BASE / 'benchmarks/cmb-points.dat') == reference['points_sha256'])
    expected_rows = reference['backends'][version]
    check('four benchmark indices', [r['index'] for r in rows] == [r['index'] for r in expected_rows])
    for row, expected in zip(rows, expected_rows):
        label = f"point {row['index']}"
        raw = (output / f"OUT_mO_{row['index']}").read_text()
        lines = [line.removeprefix('TRSM_PlanckCMB_v1 ') for line in raw.splitlines()
                 if line.startswith('TRSM_PlanckCMB_v1 ')]
        check(label + ' exactly one native CMB result', len(lines) == 1)
        payload = json.loads(lines[0]) if len(lines) == 1 else {}
        omega = row.get('dm_omega')
        raw_ratio = payload.get('ratio_raw')
        valid = all(isinstance(v, (int, float)) and math.isfinite(v) and v >= 0
                    for v in (omega, raw_ratio)) and payload.get('status') == 'ok'
        check(label + ' finite successful native CMB result', valid)
        if valid:
            # Deliberately do not call assess_cmb_limit: independent arithmetic on raw evidence.
            fraction = min(1., omega / .12)
            ratio = fraction * fraction * raw_ratio
            check(label + ' xi=min(1,Omega/0.12)', row.get('dm_cmb_abundance_fraction') == fraction)
            check(label + ' raw Planck ratio preserved', row.get('dm_cmb_ratio_raw') == raw_ratio)
            check(label + ' ratio=raw*xi^2', math.isclose(row.get('dm_cmb_ratio') if row.get('dm_cmb_ratio') is not None else math.nan, ratio, rel_tol=2e-14, abs_tol=1e-30))
            check(label + ' strict exclusion at ratio>1', row.get('dm_cmb_excluded') is (ratio > 1))
        check(label + ' benchmark DM mass', row.get('MX') == expected['MX'])
        check(label + ' darkOmega success', row.get('dm_solver_error') == 0)
        for stage in ('relic', 'indirect'):
            check(label + f' loop hook in {stage}', (row.get(f'dm_loop_hook_{stage}_calls') or 0) > 0)
        inputs = json.loads(row.get('dm_actual_inputs') or '{}')
        check(label + ' h1 pole 125.09 GeV', inputs.get('Mh') == 125.09)
        check(label + ' virtual W/Z enabled', inputs.get('VWdecay') == inputs.get('VZdecay') == 1)
        for key in QUANTITIES:
            value = row.get(key)
            check(label + ' reference ' + key,
                  value is not None and math.isclose(value, expected[key], rel_tol=rtol, abs_tol=1e-30))
        for key in ('dm_cmb_excluded', 'dm_passed'):
            check(label + ' reference ' + key, row.get(key) is expected[key])
    report = {'schema': 'trsm_example_cmb_validation_v1', 'backend_version': version,
              'reference_relative_tolerance': rtol, 'reference_absolute_tolerance': 1e-30,
              'reference_sha256': sha256(BASE / 'benchmarks/cmb-reference-v2.json'),
              'passed': all(c['passed'] for c in checks), 'checks': checks,
              'scope': 'Native benchmark regression plus independent abundance/rescaling arithmetic; '
                       'not an independent relic solver or Planck likelihood calculation.'}
    write_json(output / 'validation.json', report)
    write_tsv(output / 'cmb-comparison.tsv', [{key: r.get(key) for key in
        ('index', 'MX', *QUANTITIES, 'dm_cmb_excluded', 'dm_passed')} for r in rows])
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True, help='New directory, or existing results with --check-only')
    parser.add_argument('--micromegas-version', choices=('6', '7', '6.1.15', '7.1.4'), default='7.1.4')
    parser.add_argument('--micromegas-main', type=Path)
    parser.add_argument('--check-only', action='store_true', help='Validate already produced benchmark logs/results')
    parser.add_argument('--plot', action='store_true', help='Also write cmb_ratios.png (requires matplotlib)')
    args = parser.parse_args(argv)
    version = normalize_micromegas_version(args.micromegas_version)
    if not args.check_only:
        command = [sys.executable, str(BASE / 'run_scan.py'), '--input', str(BASE / 'benchmarks/cmb-points.dat'),
                   '--output-dir', str(args.output_dir), '--micromegas-version', version, '--planck-cmb']
        if args.micromegas_main:
            command += ['--micromegas-main', str(args.micromegas_main)]
        # A failed calculation should leave its raw evidence and an explicit failing report.
        completed = subprocess.run(command)
        if completed.returncode and not (args.output_dir / 'results.json').is_file():
            return completed.returncode
    report = validate(args.output_dir, version)
    failed = [c['check'] for c in report['checks'] if not c['passed']]
    print(f"CMB validation: {len(report['checks']) - len(failed)}/{len(report['checks'])} checks pass")
    for failure in failed:
        print('FAIL:', failure)
    print('Comparison:', args.output_dir.resolve() / 'cmb-comparison.tsv')
    if args.plot:
        subprocess.run([sys.executable, str(BASE / 'plot/plot_cmb.py'), str(args.output_dir)], check=True)
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
