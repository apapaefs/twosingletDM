#!/usr/bin/env python3
"""A standalone, inspectable front end to the production v2 DM assessment."""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from dataclasses import asdict, dataclass

from _bootstrap import repository_root
REPO_ROOT = repository_root()

from test_trsm_DM import DMPoint, test_dm, write_micromegas_card
from trsm_cmb import add_cmb_arguments, cmb_configuration, cmb_diagnostics, require_cmb_capability
from trsm_direct_detection import DEFAULT_LIMIT_MODEL, DEFAULT_LIMIT_TABLE, load_si_limit_table
from trsm_inputs import (ABUNDANCE_REFERENCE, PHYSICS_VERSION, RELIC_UPPER_LIMIT,
                         json_safe, micromegas_sm_inputs, sm_inputs)
from trsm_micromegas import (DEFAULT_MICROMEGAS_VERSION, default_micromegas_main,
                            normalize_micromegas_version, require_v2_capability)

CARD_FIELDS = (('LX', 'lx'), ('LHX', 'lhx'), ('LSX', 'lsx'), ('MX', 'mx'),
               ('vevs', 'vevs'), ('SinT', 'sint'), ('Mh2', 'mh2'))
LEGACY_COLUMNS = ('index', 'LX', 'LHX', 'LSX', 'MX', 'vevs', 'SinT', 'Mh2',
                  'MDM', 'Omega', 'DirDet', 'DirDetLimit', 'DirDetBaseLimit',
                  'IndirAvailable', 'IndirEnergy', 'IndirFlux', 'IndirLimit', 'IndirRatio')
SCHEMA = 'trsm_steer_example_v2'


@dataclass(frozen=True)
class PointInput:
    index: int
    lx: float
    lhx: float
    lsx: float
    mx: float
    vevs: float
    sint: float
    mh2: float

    def __post_init__(self):
        if not isinstance(self.index, int) or self.index < 0:
            raise ValueError('Point index must be a nonnegative integer')
        if not all(math.isfinite(getattr(self, key)) for _, key in CARD_FIELDS):
            raise ValueError('All input parameters must be finite')
        if self.mx <= 0 or self.mh2 <= 0 or self.vevs == 0 or abs(self.sint) >= 1:
            raise ValueError('Require MX, Mh2 > 0, vevs != 0 and |SinT| < 1')

    def dm_point(self):
        return DMPoint(self.lx, self.lhx, self.lsx, self.mx, self.vevs,
                       math.asin(self.sint), self.mh2)

    def columns(self):
        return {'index': self.index, **{label: getattr(self, key) for label, key in CARD_FIELDS}}


def read_card(card_path, index=None):
    card_path = Path(card_path)
    if index is None:
        match = re.fullmatch(r'MO_inp(\d+)\.dat', card_path.name)
        if not match:
            raise ValueError('Use --index with a card not named MO_inp<index>.dat')
        index = int(match[1])
    values = {}
    for line in card_path.read_text().splitlines():
        fields = line.split('#', 1)[0].split()
        if not fields:
            continue
        if len(fields) != 2 or fields[0] in values:
            raise ValueError(f'Malformed or duplicate card entry: {line}')
        values[fields[0]] = float(fields[1])
    allowed = {name for name, _ in CARD_FIELDS} | set(micromegas_sm_inputs())
    if set(values) - allowed:
        raise ValueError(f'Unsupported card overrides: {sorted(set(values) - allowed)}')
    for name, value in micromegas_sm_inputs().items():
        if name in values and not math.isclose(values[name], value, rel_tol=2e-14):
            raise ValueError(f'{name} conflicts with the shared v2 SM inputs ({value:.17g})')
    missing = [name for name, _ in CARD_FIELDS if name not in values]
    if missing:
        raise ValueError(f'Missing card inputs: {missing}')
    return PointInput(index, **{key: values[name] for name, key in CARD_FIELDS})


def read_points(path):
    points, seen = [], set()
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        fields = line.split('#', 1)[0].split()
        if not fields:
            continue
        if len(fields) != 8:
            raise ValueError(f'{path}:{number}: expected index and seven parameters')
        point = PointInput(int(fields[0]), *map(float, fields[1:]))
        if point.index in seen:
            raise ValueError(f'Duplicate point index {point.index}')
        seen.add(point.index)
        points.append(point)
    if not points:
        raise ValueError(f'No points in {path}')
    return points


def write_card(point, card_path):
    write_micromegas_card(point.dm_point(), Path(card_path))


def add_physics_arguments(parser):
    parser.add_argument('--micromegas-version', choices=('6', '7', '6.1.15', '7.1.4'),
                        default=DEFAULT_MICROMEGAS_VERSION)
    parser.add_argument('--micromegas-main', type=Path, help='Path to the rebuilt v2 TRSM/main')
    add_cmb_arguments(parser)
    parser.add_argument('--no-rescale', action='store_true', help='Use unit DM abundance for all signals')
    parser.add_argument('--relic-upper-limit', type=float, default=RELIC_UPPER_LIMIT,
                        help='Relic-density upper cut (default: %(default)s); abundance reference remains 0.12')
    parser.add_argument('--limit-model', default=DEFAULT_LIMIT_MODEL,
                        choices=(DEFAULT_LIMIT_MODEL, 'legacy-output', 'lz2025-source'))
    parser.add_argument('--limit-table', type=Path, help='Explicit normalized SI table (default: LZ WS2024)')
    parser.add_argument('--timeout', type=float, default=600, help='Seconds allowed per native point')


def configure(args, *, replay=False):
    args.micromegas_version = normalize_micromegas_version(args.micromegas_version)
    if args.planck_cmb is None:
        args.planck_cmb = args.micromegas_version == '7.1.4'
    if args.timeout <= 0 or not math.isfinite(args.timeout):
        raise ValueError('--timeout must be finite and positive')
    if args.relic_upper_limit <= 0 or not math.isfinite(args.relic_upper_limit):
        raise ValueError('--relic-upper-limit must be finite and positive')
    if args.limit_table and args.limit_model != DEFAULT_LIMIT_MODEL:
        raise ValueError('Use either --limit-table or a historical --limit-model')
    args.si_table = (load_si_limit_table(args.limit_table or DEFAULT_LIMIT_TABLE)
                     if args.limit_model == DEFAULT_LIMIT_MODEL else None)
    args.micromegas_main = Path(args.micromegas_main or os.environ.get('MICROMEGAS_MAIN') or
                               default_micromegas_main(args.micromegas_version)).expanduser().resolve()
    args.capabilities = None
    if not replay:
        args.capabilities = require_v2_capability(args.micromegas_main)
        if args.planck_cmb:
            args._cmb_driver = require_cmb_capability(args.micromegas_main)
    return args


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def metadata(args, *, replay=False):
    revision = subprocess.run(['git', '-C', str(REPO_ROOT), 'rev-parse', 'HEAD'],
                              capture_output=True, text=True)
    sources = ['test_trsm_DM.py', 'trsm_cmb.py', 'trsm_direct_detection.py',
               'trsm_inputs.py', 'trsm_micromegas.py', 'config/sm-inputs-v2.json',
               'DM/steer_example/run_single_point.py', 'DM/steer_example/run_scan.py',
               'DM/steer_example/cutflow.py']
    return {
        'schema': SCHEMA, 'physics_version': PHYSICS_VERSION,
        'source_commit': revision.stdout.strip() if revision.returncode == 0 else None,
        'source_sha256': {name: sha256(REPO_ROOT / name) for name in sources},
        'mode': 'raw_output_replay' if replay else 'native',
        'backend_version_requested': args.micromegas_version,
        'executable': None if replay else str(args.micromegas_main),
        'executable_sha256': None if replay else sha256(args.micromegas_main),
        'driver_capabilities': args.capabilities, 'sm_inputs': sm_inputs(),
        'relic_upper_limit': args.relic_upper_limit, 'abundance_reference': ABUNDANCE_REFERENCE,
        'rescale': not args.no_rescale, 'planck_cmb': {
            **cmb_configuration(args), 'rescale': not args.no_rescale,
            'abundance_rescaling': '1' if args.no_rescale else 'min(1, Omega_h2 / 0.12)^2'},
        'direct_detection': args.si_table.metadata() if args.si_table else {'model': args.limit_model},
        'legacy_columns': list(LEGACY_COLUMNS), 'timeout_seconds': args.timeout,
    }


def write_json(path, value):
    Path(path).write_text(json.dumps(json_safe(value), indent=2, allow_nan=False) + '\n')


def tsv_value(value):
    if value is None:
        return 'nan'
    if isinstance(value, bool):
        return '1' if value else '0'
    if isinstance(value, float):
        return format(value, '.17g')
    return str(value)


def write_tsv(path, rows):
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, columns, delimiter='\t')
        writer.writeheader()
        writer.writerows({key: tsv_value(row.get(key)) for key in columns} for row in rows)


def legacy_values(row):
    keys = ('dm_mdm', 'dm_omega', 'dm_dir_det', 'dm_dir_det_limit', 'dm_lux_base_limit',
            'dm_indirect_available', 'dm_indirect_energy', 'dm_indirect_flux',
            'dm_indirect_limit', 'dm_indirect_ratio')
    return [row[key] for key in LEGACY_COLUMNS[:8]] + [row.get(key) for key in keys]


def legacy_line(row, *, scan=False):
    values = legacy_values(row)
    return '\t'.join(map(tsv_value, values[:11] if scan else values)) + '\n'


class NativeRunner:
    """One isolated writable CalcHEP cache per process/runner, reused across points."""
    def __init__(self, args):
        self.args = args
        self.directory = tempfile.TemporaryDirectory(prefix='trsm_example_')

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.directory.cleanup()

    def run(self, card):
        command = [str(self.args.micromegas_main), str(Path(card).resolve())]
        if self.args.planck_cmb:
            command.append('--planck-cmb')
        invocation = {'command': command, 'returncode': None, 'error': None,
                      'runtime_directory': self.directory.name}
        try:
            done = subprocess.run(command, cwd=self.directory.name,
                                  env=dict(os.environ, TRSM_RUNTIME_DIR=self.directory.name),
                                  capture_output=True, text=True, timeout=self.args.timeout)
            invocation['returncode'] = done.returncode
            if done.returncode:
                invocation['error'] = f'micrOMEGAs exited with status {done.returncode}'
            return done.stdout, done.stderr, invocation
        except subprocess.TimeoutExpired as error:
            invocation['error'] = f'micrOMEGAs exceeded {self.args.timeout:g} seconds'
            def decode(text):
                return text.decode(errors='replace') if isinstance(text, bytes) else text or ''
            return decode(error.stdout), decode(error.stderr), invocation
        except OSError as error:
            invocation['error'] = str(error)
            return '', '', invocation


def evaluate_point(point, args, output_dir, *, runner=None, raw_output=None, raw_source=None):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = str(point.index)
    card = output_dir / f'MO_inp{stem}.dat'
    result_path = output_dir / f'result_{stem}.json'
    if result_path.exists() or (output_dir / f'OUT_mO_{stem}').exists() or card.exists():
        raise ValueError(f'Output for point {stem} already exists in {output_dir}; use a new directory')
    write_card(point, card)
    if raw_output is None:
        raw_output, stderr, invocation = runner.run(card)
    else:
        stderr = ''
        invocation = {'mode': 'raw_output_replay', 'source': str(raw_source) if raw_source else None,
                      'returncode': None, 'error': None}
    invocation['card_sha256'] = sha256(card)
    (output_dir / f'OUT_mO_{stem}').write_text(raw_output)
    (output_dir / f'ERR_mO_{stem}').write_text(stderr)
    passed, _info, diagnostics = test_dm(**asdict(point.dm_point()), raw_output=raw_output,
                                      limit_model=args.limit_model, limit_table=args.si_table,
                                      relic_upper_limit=args.relic_upper_limit,
                                      rescale=not args.no_rescale, planck_cmb=args.planck_cmb)
    if not args.planck_cmb:
        diagnostics.update(cmb_diagnostics(enabled=False))
    if invocation['error']:
        passed = None
        diagnostics.update(dm_calculation_status='execution_error',
                           dm_assessment_reason=invocation['error'], dm_xf=None,
                           dm_freezeout_temperature_GeV=None, dm_direct_detection_available=False,
                           dm_indirect_available=False)
        for name in ('relic', 'direct_detection', 'indirect_detection'):
            diagnostics[f'dm_{name}_excluded'] = None
        diagnostics.update(cmb_diagnostics(enabled=args.planck_cmb, reason=invocation['error']))
    # Keep a failed/missing CMB diagnosis visible even when another constraint excludes.
    row = json_safe({**point.columns(), 'dm_passed': passed, **diagnostics})
    write_json(result_path, {'schema': SCHEMA, 'point': row, 'invocation': invocation,
                             'raw_output_sha256': hashlib.sha256(raw_output.encode()).hexdigest()})
    write_tsv(output_dir / f'result_{stem}.tsv', [row])
    (output_dir / f'DM_data_{stem}').write_text(legacy_line(row))
    # Keep the directory layout consumed by the original shell/plot workflow.
    for subdir, filename, content in (
        ('OUT_mO', f'OUT_mO_{stem}', raw_output),
        ('DM_data', f'DM_data_{stem}', legacy_line(row)),
    ):
        (output_dir / subdir).mkdir(exist_ok=True)
        (output_dir / subdir / filename).write_text(content)
    (output_dir / f'scan_result_{stem}.dat').write_text(legacy_line(row, scan=True))
    return row


def print_summary(row):
    def verdict(value):
        return 'unassessed' if value is None else 'excluded' if value else 'pass'
    print(f"Point {row['index']}: Omega={row.get('dm_omega')}  Xf={row.get('dm_xf')}  Tf={row.get('dm_freezeout_temperature_GeV')} GeV")
    for label, key in (('Relic density', 'relic'), ('Direct detection', 'direct_detection'),
                       ('Gamma lines', 'indirect_detection')):
        print(f"  {label}: {verdict(row.get('dm_' + key + '_excluded'))}")
    print(f"  Planck CMB: raw={row.get('dm_cmb_ratio_raw')}  xi={row.get('dm_cmb_abundance_fraction')}  rescaled={row.get('dm_cmb_ratio')}  {verdict(row.get('dm_cmb_excluded')) if row.get('dm_cmb_enabled') else 'disabled'} ({row['dm_cmb_status']})")
    print('  DM:', 'unassessed' if row['dm_passed'] is None else 'pass' if row['dm_passed'] else 'excluded',
          '|', row.get('dm_assessment_reason', ''))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--card', type=Path)
    parser.add_argument('--index', type=int)
    for _, key in CARD_FIELDS:
        parser.add_argument('--' + key, type=float)
    parser.add_argument('--micromegas-output', type=Path, help='Reassess a saved raw log without running a backend')
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parent / 'single_point_output')
    add_physics_arguments(parser)
    args = parser.parse_args(argv)
    if args.card and any(getattr(args, key) is not None for _, key in CARD_FIELDS):
        parser.error('Use --card or explicit point parameters, not both')
    if not args.card and any(getattr(args, key) is None for _, key in CARD_FIELDS):
        parser.error('Provide --card or all seven point parameters')
    return args


def main(argv=None):
    try:
        args = parse_args(argv)
        replay = args.micromegas_output is not None
        configure(args, replay=replay)
        point = (read_card(args.card, args.index) if args.card else
                 PointInput(1 if args.index is None else args.index,
                            **{key: getattr(args, key) for _, key in CARD_FIELDS}))
        if replay:
            row = evaluate_point(point, args, args.output_dir,
                                 raw_output=args.micromegas_output.read_text(), raw_source=args.micromegas_output)
        else:
            with NativeRunner(args) as runner:
                row = evaluate_point(point, args, args.output_dir, runner=runner)
        write_json(args.output_dir / f'metadata_{point.index}.json', metadata(args, replay=replay))
        print_summary(row)
        print('Saved:', args.output_dir.resolve() / f'result_{point.index}.json')
        return 1 if row.get('dm_calculation_status') != 'success' or (args.planck_cmb and not row['dm_cmb_available']) else 0
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f'Error: {error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
