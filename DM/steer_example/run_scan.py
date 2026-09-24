#!/usr/bin/env python3
"""Run an oks.dat list or a directory of cards, retaining CMB and failure diagnostics."""
import argparse
from contextlib import nullcontext
from pathlib import Path
import sys
import subprocess

from run_single_point import (NativeRunner, add_physics_arguments, configure, evaluate_point,
                              legacy_line, metadata, print_summary, read_card, read_points,
                              tsv_value, write_json, write_tsv)


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


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--input', type=Path, help='Eight-column oks.dat')
    group.add_argument('--cards', type=Path, help='Directory of MO_inp<index>.dat cards')
    parser.add_argument('--output-dir', type=Path, required=True, help='New, empty output directory')
    parser.add_argument('--raw-output-dir', type=Path, help='Replay flat OUT_mO_<index> logs without a backend')
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
                if not (args.raw_output_dir / f'OUT_mO_{point.index}').is_file():
                    raise ValueError(f'Missing raw output for point {point.index}')
        args.output_dir.mkdir(parents=True, exist_ok=True)
        write_json(args.output_dir / 'metadata.json', metadata(args, replay=replay))
        rows = []
        with nullcontext() if replay else NativeRunner(args) as runner:
            for point in points:
                raw_path = args.raw_output_dir / f'OUT_mO_{point.index}' if replay else None
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
