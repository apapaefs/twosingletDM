#!/usr/bin/env python3
"""Evaluate an oks.dat list or cards with relic, DD, gamma-line and optional CMB cuts."""
import argparse
from contextlib import nullcontext, redirect_stdout, redirect_stderr
from pathlib import Path
import sys
import subprocess

from run_single_point import (NativeRunner, add_physics_arguments, configure, evaluate_point,
                              legacy_line, metadata, print_summary, read_card, read_points,
                              tsv_value, write_json, write_tsv)
from cutflow import stages


class Tee:
    def __init__(self, terminal, log):
        self.terminal, self.log = terminal, log

    def write(self, text):
        self.log.write(text)
        self.log.flush()
        return self.terminal.write(text)

    def flush(self):
        self.terminal.flush()
        self.log.flush()


def raw_output_path(directory, index):
    for path in (directory / f'OUT_mO_{index}', directory / 'OUT_mO' / f'OUT_mO_{index}'):
        if path.is_file():
            return path
    raise ValueError(f'Missing raw output for point {index} in {directory}')


def aggregate(output, rows):
    """Rebuild tables from evaluated records; never append stale/duplicate rows."""
    output = Path(output)
    write_tsv(output / 'results.tsv', rows)
    write_json(output / 'results.json', rows)
    (output / 'oks.dat').write_text(''.join('\t'.join(tsv_value(row[key]) for key in
        ('index', 'LX', 'LHX', 'LSX', 'MX', 'vevs', 'SinT', 'Mh2')) + '\n' for row in rows))
    # Legacy plots have no status column; only assessed solver results enter them.
    valid = [r for r in rows if r.get('dm_calculation_status') == 'success']
    (output / 'scan_results.dat').write_text(''.join(legacy_line(r, scan=True) for r in valid))
    def passes(r, name):
        return r.get('dm_' + name + '_excluded') is False
    def fails(r, name):
        return r.get('dm_' + name + '_excluded') is True
    subsets = {
        'allall': lambda r: r['dm_passed'] is True,
        'dmexcl': lambda r: r['dm_passed'] is False,
        'dmunassessed': lambda r: r['dm_passed'] is None,
        'relic_pass': lambda r: passes(r, 'relic'),
        'relic_strict': lambda r: passes(r, 'relic') and .119 <= r['dm_omega'] <= .121,
        'omexcl': lambda r: fails(r, 'relic'),
        'luxpass': lambda r: passes(r, 'direct_detection'),
        'luxexcl': lambda r: fails(r, 'direct_detection'),
        'all_dirpass': lambda r: passes(r, 'relic') and passes(r, 'direct_detection'),
        'all_indirpass': lambda r: all(passes(r, name) for name in ('relic', 'direct_detection', 'indirect_detection')),
        'omgpass_dirfail': lambda r: passes(r, 'relic') and fails(r, 'direct_detection'),
        'indirpass': lambda r: passes(r, 'indirect_detection') and r.get('dm_indirect_available'),
        'indirexcl': lambda r: fails(r, 'indirect_detection'),
        'indir_caughtit': lambda r: passes(r, 'relic') and passes(r, 'direct_detection') and fails(r, 'indirect_detection'),
        'dir_caughtit': lambda r: passes(r, 'relic') and fails(r, 'direct_detection') and passes(r, 'indirect_detection'),
        'dir_indir_caughtit': lambda r: passes(r, 'relic') and fails(r, 'direct_detection') and fails(r, 'indirect_detection'),
        'cmbpass': lambda r: passes(r, 'cmb'),
        'cmbexcl': lambda r: fails(r, 'cmb'),
        'cmbunassessed': lambda r: r.get('dm_cmb_enabled') and r.get('dm_cmb_excluded') is None,
        'cmb_caughtit': lambda r: fails(r, 'cmb') and all(passes(r, name) for name in ('relic', 'direct_detection', 'indirect_detection')),
    }
    counts = {}
    for name, selected in subsets.items():
        selected_rows = [r for r in rows if selected(r)]
        (output / (name + '.dat')).write_text(''.join(legacy_line(r) for r in selected_rows))
        counts[name] = len(selected_rows)
    write_json(output / 'counts.json', {'points': len(rows), **counts})
    write_json(output / 'cutflow.json', list(stages(rows)))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--input', type=Path, help='Eight-column oks.dat')
    group.add_argument('--cards', type=Path, help='Directory of MO_inp<index>.dat cards')
    parser.add_argument('--output-dir', type=Path, required=True, help='New, empty output directory')
    parser.add_argument('--raw-output-dir', type=Path, help='Replay flat or historical OUT_mO/ logs without a backend')
    parser.add_argument('--log-file', type=Path, help='Terminal log (default: <output-dir>/MOrun.log)')
    add_physics_arguments(parser)
    args = parser.parse_args(argv)
    try:
        replay = args.raw_output_dir is not None
        configure(args, replay=replay)
        points = (read_points(args.input) if args.input else
                  sorted((read_card(p) for p in args.cards.glob('MO_inp*.dat')), key=lambda p: p.index))
        if not points or len({p.index for p in points}) != len(points):
            raise ValueError('Require at least one point with unique indices')
        if args.output_dir.exists() and any(args.output_dir.iterdir()):
            raise ValueError('Output directory is not empty; use a new directory')
        if replay:
            for point in points:
                raw_output_path(args.raw_output_dir, point.index)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        write_json(args.output_dir / 'metadata.json', metadata(args, replay=replay))
        rows = []
        log_path = args.log_file or args.output_dir / 'MOrun.log'
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open('a') as log, redirect_stdout(Tee(sys.stdout, log)), redirect_stderr(Tee(sys.stderr, log)):
            with nullcontext() if replay else NativeRunner(args) as runner:
                for point in points:
                    raw_path = raw_output_path(args.raw_output_dir, point.index) if replay else None
                    row = evaluate_point(point, args, args.output_dir, runner=runner,
                                         raw_output=raw_path.read_text() if replay else None, raw_source=raw_path)
                    rows.append(row)
                    aggregate(args.output_dir, rows)
                    print_summary(row)
                    sys.stdout.flush()
        print('Saved:', args.output_dir.resolve() / 'results.tsv')
        return int(any(r.get('dm_calculation_status') != 'success' or
                       (args.planck_cmb and not r['dm_cmb_available']) for r in rows))
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f'Error: {error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
