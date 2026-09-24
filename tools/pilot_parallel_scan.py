#!/usr/bin/env python3
"""Compare a small serial/parallel native EWPT campaign and record sampled resources.

Default ranges are a known EWPT-eligible smoke region, not a production profile.
Use --config to measure the intended production profile instead.
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from trsm_scan_campaign import atomic_write_json


def process_snapshot(pgids):
    """RSS is summed over a worker's Python/native process group (may share pages)."""
    if not pgids:
        return {}, {}
    fields = 'pid=,pgid=,rss=,comm='
    output = subprocess.check_output(['ps', '-axo', fields], text=True)
    groups, names = {}, {}
    for line in output.splitlines():
        pid, group, rss, command = line.strip().split(None, 3)
        if int(group) in pgids:
            groups[int(group)] = groups.get(int(group), 0) + int(rss) * 1024
            names[int(pid)] = Path(command).name
    native = {pid: name for pid, name in names.items() if name in ('CalcTemps', 'MinimaTracer', 'PhaseProbe')}
    threads = {}
    if native:
        if sys.platform == 'darwin':
            listing = subprocess.run(['ps', '-M', '-p', ','.join(map(str, native)), '-o', 'pid='],
                                     capture_output=True, text=True)
            counts = {}
            for line in listing.stdout.splitlines():
                values = line.split()
                # Darwin -M prepends its USER/PID thread table even with -o.
                token = values[0] if values and values[0].isdigit() else (values[1] if len(values) > 1 else '')
                if token.isdigit():
                    pid = int(token); counts[pid] = counts.get(pid, 0) + 1
        else:
            listing = subprocess.run(['ps', '-p', ','.join(map(str, native)), '-o', 'pid=,nlwp='],
                                     capture_output=True, text=True)
            counts = {int(parts[0]): int(parts[1]) for line in listing.stdout.splitlines()
                      if len(parts := line.split()) == 2}
        for pid, count in counts.items():
            if pid in native:
                threads[native[pid]] = max(threads.get(native[pid], 0), count)
    return groups, threads


def measure(command, directory, log_path, sample_seconds, timeout):
    start = time.monotonic()
    peak_group_rss = peak_total_rss = 0
    peak_threads = {}
    samples = 0
    measurement_error = None
    with log_path.open('w') as log:
        process = subprocess.Popen(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        try:
            while process.poll() is None:
                if time.monotonic() - start > timeout:
                    raise TimeoutError(f'Pilot exceeded {timeout:g} seconds; see {log_path}')
                state_path = directory / 'campaign_state.json'
                if state_path.is_file() and measurement_error is None:
                    state = json.loads(state_path.read_text())
                    pgids = {row['pid'] for row in state['seeds'] if row.get('pid') and row['status'] == 'running'}
                    try:
                        rss, threads = process_snapshot(pgids)
                        if rss:
                            samples += 1
                            peak_group_rss = max(peak_group_rss, max(rss.values()))
                            peak_total_rss = max(peak_total_rss, sum(rss.values()))
                        for name, count in threads.items():
                            peak_threads[name] = max(peak_threads.get(name, 0), count)
                    except (OSError, subprocess.SubprocessError, ValueError) as error:
                        measurement_error = str(error)
                time.sleep(sample_seconds)
        finally:
            if process.poll() is None:
                process.terminate()
                process.wait(timeout=60)
    elapsed = time.monotonic() - start
    if process.returncode:
        raise RuntimeError(f'Pilot failed ({process.returncode}); see {log_path}')
    summary = json.loads((directory / 'campaign_summary.json').read_text())
    return {'elapsed_seconds': elapsed, 'totals': summary['totals'],
            'draws_per_second': summary['totals']['draw_count'] / elapsed,
            'evo_thc_per_second': summary['totals']['evo_thc_count'] / elapsed,
            'peak_worker_group_rss_bytes': peak_group_rss,
            'peak_simultaneous_worker_rss_bytes': peak_total_rss,
            'native_max_threads': peak_threads, 'resource_samples': samples,
            'resource_measurement_error': measurement_error,
            'campaign_directory_bytes': sum(p.stat().st_size for p in directory.rglob('*') if p.is_file())}


def compare(serial, parallel):
    from itertools import zip_longest
    exact = ('seed', 'point_index', 'M2', 'M3', 'vs', 'a12', 'lX', 'lPhiX', 'lSX',
             'evo', 'thc', 'experimental_subset', 'dm_subset', 'vacuum_tree_global',
             'rg_bfb', 'rg_unitarity', 'ewpt_baryo_candidate', 'ewpt_gw_candidate',
             'ewpt_execution_status', 'ewpt_equilibrium_status')
    numeric = {'dm_omega': (2e-6, 1e-30), 'dm_xf': (2e-6, 1e-30),
               'dm_dir_det': (2e-6, 1e-30), 'higgstools_hs_chi2': (1e-8, 1e-10),
               'ewpt_ew_entry_temperature_GeV': (0, .003),
               'ewpt_ew_entry_jump_over_T': (1e-3, 1e-10),
               'ewpt_gw_max_field_jump_over_T': (1e-3, 1e-10)}
    count = ewpt = 0
    with (serial / 'combined_points.tsv').open() as first, (parallel / 'combined_points.tsv').open() as second:
        for a, b in zip_longest(csv.DictReader(first, delimiter='\t'), csv.DictReader(second, delimiter='\t')):
            if a is None or b is None: raise AssertionError('Serial/parallel row counts differ')
            for key in exact:
                if a.get(key) != b.get(key): raise AssertionError(f'{key} differs: {a.get(key)} / {b.get(key)}')
            for key, (relative, absolute) in numeric.items():
                x, y = a.get(key), b.get(key)
                if x == y: continue
                try: x, y = float(x), float(y)
                except (TypeError, ValueError): raise AssertionError(f'{key} availability differs')
                if not (math.isnan(x) and math.isnan(y)) and not math.isclose(x, y, rel_tol=relative, abs_tol=absolute):
                    raise AssertionError(f'{key} differs: {x} / {y}')
            if a.get('ewpt_eligible') == 'True':
                for row in (a, b):
                    if row.get('ewpt_calctemps_returncode') != '0' or row.get('ewpt_minimatracer_returncode') != '0':
                        raise AssertionError('Eligible pilot point did not complete both EWPT executables')
                ewpt += 1
            count += 1
    if not ewpt: raise AssertionError('Pilot did not exercise EWPT; choose an eligible smoke region')
    return {'matched_rows': count, 'matched_ewpt_rows': ewpt}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--config', type=Path)
    parser.add_argument('--seed-start', type=int, default=917319)
    parser.add_argument('--nseeds', type=int, default=2)
    parser.add_argument('--nrandom', type=int, default=1)
    parser.add_argument('--jobs', type=int, default=2)
    parser.add_argument('--sample-seconds', type=float, default=.5)
    parser.add_argument('--timeout-seconds', type=float, default=3600)
    args = parser.parse_args()
    if min(args.nseeds, args.nrandom, args.jobs, args.sample_seconds, args.timeout_seconds) <= 0:
        parser.error('Counts, intervals and timeout must be positive')
    destination = args.output_dir.expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=False)
    if args.config:
        config = json.loads(args.config.read_text())
    else:
        config = {'constraint_version': 'trsm_constraints_v2', 'seed_start': args.seed_start,
                  'nseeds': args.nseeds, 'points_per_seed': args.nrandom, 'jobs': args.jobs,
                  'nrandom_count_evo_thc': True, 'ewpt_thigh_GeV': 300,
                  'generator_arguments': ['--independent-m3', '--m2-min', '180', '--m2-max', '220',
                    '--m3-min', '50', '--m3-max', '55', '--vs-min', '290', '--vs-max', '310',
                    '--k1-min', '.99', '--k1-max', '1', '--lx-min', '.05', '--lx-max', '.2', '--no-print-info']}
    config_path = destination / 'configuration.json'
    atomic_write_json(config_path, config)
    report = {'profile': 'supplied' if args.config else 'EWPT-eligible smoke region',
              'memory_metric': 'Sampled summed RSS including Python and native descendants; shared pages can be counted more than once',
              'runs': {}}
    for label, jobs in (('serial', 1), ('parallel', args.jobs)):
        directory = destination / label
        command = [sys.executable, str(ROOT / 'tools/run_next_scan.py'), '--config', str(config_path),
                   '--campaign-dir', str(directory), '--seed-start', str(args.seed_start), '--nseeds', str(args.nseeds),
                   '--nrandom', str(args.nrandom), '--nrandom-count-evo-thc', '--jobs', str(jobs), '--run']
        print(f'Running {label} pilot: {args.nseeds} scans, target {args.nrandom}, jobs {jobs}', flush=True)
        report['runs'][label] = measure(command, directory, destination / (label + '.log'),
                                        args.sample_seconds, args.timeout_seconds)
        atomic_write_json(destination / 'report.json', report)
    report['comparison'] = compare(destination / 'serial', destination / 'parallel')
    for run in report['runs'].values():
        if run['resource_measurement_error'] or not run['resource_samples']:
            raise RuntimeError('Process metrics unavailable; rerun where ps can inspect owned workers')
        if not run['native_max_threads'] or any(count > 1 for count in run['native_max_threads'].values()):
            raise AssertionError('Native EWPT thread counts were not confirmed to be one')
    report['status'] = 'passed'
    atomic_write_json(destination / 'report.json', report)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
