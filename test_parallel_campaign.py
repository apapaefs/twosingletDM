"""Integration tests for bounded processes, recovery and committed aggregation."""
import contextlib
import csv
import io
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import textwrap
import time
import tracemalloc
import unittest
from unittest import mock

import run_trsm_seed_campaign as campaign
from trsm_parallel import THREAD_ENVIRONMENT, effective_jobs
from trsm_scan_campaign import CampaignStateError, atomic_write_json, inspect_tsv
from tools import run_next_scan
from tools import pilot_parallel_scan

ROOT = Path(__file__).resolve().parent

FAKE = r'''
import argparse, csv, json, os, random, signal, subprocess, sys, time
from pathlib import Path
from trsm_scan_campaign import (CampaignLock, atomic_write_json, inspect_tsv,
    reconcile_outputs, encode_rng_state, decode_rng_state)
p = argparse.ArgumentParser()
p.add_argument('seed', nargs='?', type=int)
p.add_argument('--nrandom', type=int, default=2)
p.add_argument('--checkpoint-every', type=int, default=1)
p.add_argument('--nrandom-count-evo-thc', action='store_true')
p.add_argument('--no-ewpt-multithreading', action='store_true')
p.add_argument('--output-manifest', type=Path)
p.add_argument('--resume-from', type=Path)
p.add_argument('--preflight', action='store_true')
p.add_argument('--sleep', type=float, default=.005)
p.add_argument('--child', action='store_true')
p.add_argument('--fail-seed', type=int)
p.add_argument('--run-ewpt', action='store_true')
p.add_argument('--ewpt-workdir')
a, unused = p.parse_known_args()
if a.preflight:
    print('TRSM_PREFLIGHT {"fixture": "deterministic"}')
    sys.exit(0)
start = time.time()
if a.resume_from:
    output = a.resume_from
    saved = json.loads(output.with_suffix('.metadata.json').read_text())
    a.seed = saved['seed']
    a.nrandom_count_evo_thc = saved['count_evo_thc']
    a.output_manifest = Path(saved['manifest'])
    a.sleep = saved['sleep']
    a.child = saved['child']
    a.fail_seed = saved['fail_seed']
    a.no_ewpt_multithreading = saved['single_thread']
else:
    output = Path.cwd() / 'output' / 'fixture.dat'
output.parent.mkdir(exist_ok=True)
checkpoint_path = output.with_suffix('.checkpoint.json')
def interrupt(signum, frame):
    raise KeyboardInterrupt()
signal.signal(signal.SIGTERM, interrupt)
with CampaignLock(output, sys.argv):
    rng = random.Random(a.seed)
    if a.resume_from:
        cp = json.loads(checkpoint_path.read_text())
        reconcile_outputs({'main': output}, cp['outputs'])
        rng.setstate(decode_rng_state(cp['rng_state']))
    else:
        output.write_text('point_index\tcoordinate\tevo\tthc\ttheory_strict_subset\tdm_subset\n')
        cp = dict(seed=a.seed, scan_path=str(output), draw_count=0, evo_thc_count=0,
                  viable_count=0, count_evo_thc=a.nrandom_count_evo_thc, target=a.nrandom)
        atomic_write_json(output.with_suffix('.metadata.json'), dict(seed=a.seed,
            count_evo_thc=a.nrandom_count_evo_thc, manifest=str(a.output_manifest),
            sleep=a.sleep, child=a.child, fail_seed=a.fail_seed,
            single_thread=a.no_ewpt_multithreading))
    def save(status):
        cp.update(status=status, outputs={'main': inspect_tsv(output)}, rng_state=encode_rng_state(rng.getstate()))
        atomic_write_json(checkpoint_path, cp)
        atomic_write_json(a.output_manifest, {'outputs': {'main': {'path': str(output)}}})
    save('running')
    Path('environment.json').write_text(json.dumps({key: os.environ[key] for key in
        ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
         'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS', 'TRSM_HIGHS_THREADS')}))
    assert a.no_ewpt_multithreading
    if a.child:
        child = subprocess.Popen([sys.executable, '-c',
            'import signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(60)'])
        Path('child.pid').write_text(str(child.pid))
    if a.seed == a.fail_seed:
        save('failed')
        sys.exit(7)
    try:
        while cp['evo_thc_count' if a.nrandom_count_evo_thc else 'draw_count'] < a.nrandom:
            index = cp['draw_count'] + 1
            evo, thc = index % 2 == 0, index % 3 != 0
            with output.open('a', newline='') as stream:
                csv.writer(stream, delimiter='\t', lineterminator='\n').writerow(
                    [index, rng.random(), evo, thc, False, False])
            cp['draw_count'] = index
            cp['evo_thc_count'] += evo and thc
            save('running')
            time.sleep(a.sleep)
        save('complete')
    except KeyboardInterrupt:
        save('interrupted')
        sys.exit(130)
Path('timing.json').write_text(json.dumps({'start': start, 'end': time.time()}))
'''


class ParallelCampaignTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.generator = self.root / 'fake.py'
        self.generator.write_text('import sys\nsys.path.insert(0, ' + repr(str(ROOT)) + ')\n' + textwrap.dedent(FAKE))

    def tearDown(self):
        self.temporary.cleanup()

    def args(self, name='campaign', *extra):
        return campaign.parse_args(['--campaign-dir', str(self.root / name), '--seed-start', '71',
            '--nseeds', '3', '--nrandom', '3', '--nrandom-count-evo-thc', '--jobs', '2',
            '--generator-script', str(self.generator), '--python-executable', sys.executable,
            '--heartbeat-seconds', '0', '--shutdown-grace-seconds', '.25', *extra])

    def run_campaign(self, args):
        with contextlib.redirect_stdout(io.StringIO()):
            return campaign.run_campaign(args)

    def resume(self, name='campaign', *extra):
        return self.run_campaign(campaign.parse_args(['--campaign-dir', str(self.root / name), '--resume', *extra]))

    def test_capacity_and_input_validation(self):
        with mock.patch('trsm_parallel.available_cpus', return_value=192):
            self.assertEqual(effective_jobs('auto', 200), 192)
            self.assertEqual(effective_jobs(300, 200), 192)
            self.assertEqual(effective_jobs(7, 200), 7)
            self.assertEqual(effective_jobs('auto', 2), 2)
            self.assertEqual(effective_jobs(1, 0), 0)
        with contextlib.redirect_stderr(io.StringIO()):
            for value in ('0', '-2', 'invalid'):
                with self.assertRaises(SystemExit): self.args('bad', '--jobs', value)

    def test_parallel_matches_serial_and_counts_only_evo_thc(self):
        serial = self.run_campaign(self.args('serial', '--jobs', '1'))
        parallel = self.run_campaign(self.args('parallel', '--generator-extra-arg=--sleep', '--generator-extra-arg=.04'))
        self.assertEqual(list(serial.combined_points), list(parallel.combined_points))
        for result in parallel.seed_results:
            self.assertEqual(result.status, 'complete')
            self.assertEqual(result.evo_thc_count, 3)
            self.assertEqual(result.draw_count, 8)
        self.assertEqual(len(parallel.combined_points), 24)
        rows = list(parallel.combined_points)
        self.assertTrue(all(row['theory_strict_subset'] == 'False' for row in rows))
        events = []
        for path in (self.root / 'parallel' / 'seeds').glob('seed_*/timing.json'):
            timing = json.loads(path.read_text())
            events += [(timing['start'], 1), (timing['end'], -1)]
            env = json.loads(path.with_name('environment.json').read_text())
            self.assertTrue(all(env[key] == '1' for key in THREAD_ENVIRONMENT))
        active = peak = 0
        for _, delta in sorted(events):
            active += delta
            peak = max(peak, active)
        self.assertLessEqual(peak, 2)
        if campaign.available_cpus() >= 2: self.assertEqual(peak, 2)

    def test_raw_draw_mode_keeps_accurate_evo_count(self):
        result = self.run_campaign(self.args('raw', '--no-nrandom-count-evo-thc', '--nrandom', '3'))
        self.assertTrue(all(row.draw_count == 3 and row.evo_thc_count == 1 for row in result.seed_results))

    def test_completed_resume_is_noop_and_aggregation_needs_no_runtime(self):
        first = self.run_campaign(self.args())
        directory = self.root / 'campaign'
        logs = {p: p.read_bytes() for p in (directory / 'logs').glob('*')}
        second = self.resume()
        self.assertEqual(list(first.combined_points), list(second.combined_points))
        self.assertTrue(all(path.read_bytes() == content for path, content in logs.items()))
        self.generator.unlink()
        result = self.run_campaign(campaign.parse_args(['--campaign-dir', str(directory), '--aggregate-only']))
        self.assertEqual(len(result.combined_points), 24)

    def test_missing_manifest_recovers_from_checkpoint(self):
        self.run_campaign(self.args())
        directory = self.root / 'campaign'
        state = json.loads((directory / 'campaign_state.json').read_text())
        for record in state['seeds']:
            Path(record['manifest']).unlink()
            record['point_output'] = ''
        atomic_write_json(directory / 'campaign_state.json', state)
        result = self.resume()
        self.assertTrue(all(row.status == 'complete' for row in result.seed_results))
        self.assertEqual(len(result.combined_points), 24)

    def test_failure_is_isolated_and_incomplete_exit_is_nonzero(self):
        result = self.run_campaign(self.args('fail', '--generator-extra-arg=--fail-seed', '--generator-extra-arg=72'))
        self.assertEqual([row.status for row in result.seed_results], ['complete', 'failed', 'complete'])
        self.assertEqual(result.seed_results[1].returncode, 7)
        self.assertEqual(len(result.combined_points), 16)
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(campaign.main(['--campaign-dir', str(self.root / 'fail'), '--aggregate-only']), 1)

    def test_missing_or_corrupt_checkpoint_excludes_only_affected_seed(self):
        self.run_campaign(self.args())
        directory = self.root / 'campaign'
        state = json.loads((directory / 'campaign_state.json').read_text())
        Path(state['seeds'][0]['point_output']).with_suffix('.checkpoint.json').unlink()
        output = Path(state['seeds'][1]['point_output'])
        data = output.read_bytes()
        output.write_bytes(data.replace(b'coordinate', b'coordinatX'))
        result = self.run_campaign(campaign.parse_args(['--campaign-dir', str(directory), '--aggregate-only']))
        self.assertEqual([row.status for row in result.seed_results], ['failed', 'failed', 'complete'])
        self.assertEqual(len(result.combined_points), 8)

    def test_partial_aggregate_reads_committed_prefix_without_changing_tail(self):
        self.run_campaign(self.args())
        directory = self.root / 'campaign'
        state = json.loads((directory / 'campaign_state.json').read_text())
        output = Path(state['seeds'][0]['point_output'])
        with output.open('ab') as stream: stream.write(b'999\tuncommitted\tTrue')
        before = output.read_bytes()
        checkpoint = json.loads(output.with_suffix('.checkpoint.json').read_text())
        checkpoint['status'] = 'interrupted'
        atomic_write_json(output.with_suffix('.checkpoint.json'), checkpoint)
        state['seeds'][0]['status'] = 'interrupted'
        atomic_write_json(directory / 'campaign_state.json', state)
        args = campaign.parse_args(['--campaign-dir', str(directory), '--aggregate-only'])
        a = self.run_campaign(args)
        product = (directory / 'combined_points.tsv').read_bytes()
        b = self.run_campaign(args)
        self.assertEqual(len(a.combined_points), 24)
        self.assertEqual(len(b.combined_points), 24)
        self.assertEqual(product, (directory / 'combined_points.tsv').read_bytes())
        self.assertEqual(output.read_bytes(), before)
        self.assertEqual(a.seed_results[0].status, 'interrupted')

    def test_rejects_configuration_and_source_changes_and_fresh_overwrite(self):
        self.run_campaign(self.args())
        with self.assertRaises(CampaignStateError): self.resume('campaign', '--nrandom', '8')
        with self.assertRaises(CampaignStateError): self.run_campaign(self.args())
        self.generator.write_text(self.generator.read_text() + '\n# changed code\n')
        with self.assertRaises(CampaignStateError): self.resume()

    def test_relative_provider_paths_are_resolved_before_isolation(self):
        args = self.args('paths', '--generator-extra-arg=--micromegas-main',
                         '--generator-extra-arg=relative/TRSM/main')
        campaign.canonicalize_paths(args)
        self.assertIn('--micromegas-main=' + str((Path.cwd() / 'relative/TRSM/main').resolve()),
                      args.generator_extra_arg)

    def test_preflight_failure_launches_no_seeds(self):
        args = self.args('preflight', '--python-executable', str(self.root / 'missing-python'))
        with self.assertRaises(OSError): self.run_campaign(args)
        state = json.loads((self.root / 'preflight/campaign_state.json').read_text())
        self.assertTrue(all(row['status'] == 'queued' and row['attempts'] == 0 for row in state['seeds']))
        self.assertFalse((self.root / 'preflight/seeds').exists())

    def test_configuration_corruption_is_rejected(self):
        self.run_campaign(self.args())
        path = self.root / 'campaign/campaign_state.json'
        state = json.loads(path.read_text())
        state['configuration']['nrandom'] = 9
        atomic_write_json(path, state)
        with self.assertRaisesRegex(CampaignStateError, 'configuration fingerprint'):
            self.resume()

    def test_resume_retains_execution_settings_unless_overridden(self):
        self.run_campaign(self.args('settings', '--checkpoint-every', '5'))
        self.resume('settings', '--jobs', '1')
        state = json.loads((self.root / 'settings/campaign_state.json').read_text())
        self.assertEqual(state['execution']['checkpoint_every'], 5)
        self.assertEqual(state['execution']['shutdown_grace_seconds'], .25)
        self.assertEqual(state['execution']['jobs'], 1)

    def test_stop_native_descendant_and_resume_matches_reference(self):
        directory = self.root / 'stopped'
        arguments = ['--campaign-dir', str(directory), '--seed-start', '71', '--nseeds', '3',
            '--nrandom', '3', '--nrandom-count-evo-thc', '--jobs', '1', '--generator-script', str(self.generator),
            '--python-executable', sys.executable, '--heartbeat-seconds', '0', '--shutdown-grace-seconds', '.3',
            '--generator-extra-arg=--child', '--generator-extra-arg=--sleep', '--generator-extra-arg=.06']
        with (self.root / 'supervisor.log').open('w') as log:
            process = subprocess.Popen([sys.executable, str(ROOT / 'run_trsm_seed_campaign.py'), *arguments], stdout=log, stderr=log)
            try:
                deadline = time.monotonic() + 15
                cp_path = directory / 'seeds/seed_71/output/fixture.checkpoint.json'
                while time.monotonic() < deadline:
                    if cp_path.exists() and json.loads(cp_path.read_text())['draw_count'] >= 1: break
                    if process.poll() is not None: self.fail((self.root / 'supervisor.log').read_text())
                    time.sleep(.02)
                else: self.fail('worker did not start')
                # A second supervisor must reject the active campaign.
                with self.assertRaises(CampaignStateError): self.resume('stopped')
                child = int((directory / 'seeds/seed_71/child.pid').read_text())
                process.send_signal(signal.SIGTERM)
                self.assertEqual(process.wait(timeout=10), 1)
                with self.assertRaises(ProcessLookupError): os.kill(child, 0)
            finally:
                if process.poll() is None:
                    process.send_signal(signal.SIGTERM)
                    process.wait(timeout=10)
        state = json.loads((directory / 'campaign_state.json').read_text())
        self.assertEqual([row['status'] for row in state['seeds']], ['interrupted', 'queued', 'queued'])
        reference = self.run_campaign(self.args('reference', '--jobs', '1'))
        resumed = self.resume('stopped', '--jobs', '2')
        self.assertEqual(list(reference.combined_points), list(resumed.combined_points))
        log_text = (directory / 'logs/seed_71.log').read_text()
        self.assertIn('attempt: 1', log_text)
        self.assertIn('attempt: 2', log_text)

    def test_aggregation_memory_is_bounded(self):
        args = self.args('large', '--nseeds', '1', '--nrandom', '12000', '--no-nrandom-count-evo-thc')
        campaign.canonicalize_paths(args)
        args.campaign_dir.mkdir()
        record = campaign.new_seed_record(71, args)
        output = Path(record['cwd']) / 'output' / 'fixture.dat'
        output.parent.mkdir(parents=True)
        with output.open('w') as stream:
            stream.write('point_index\tfiller\n')
            filler = 'x' * 4096
            for i in range(12000): stream.write(f'{i + 1}\t{filler}\n')
        cp = dict(seed=71, scan_path=str(output), status='complete', target=12000,
                  count_evo_thc=False, draw_count=12000, evo_thc_count=0, viable_count=0,
                  outputs={'main': inspect_tsv(output)})
        atomic_write_json(output.with_suffix('.checkpoint.json'), cp)
        record.update(point_output=str(output), status='complete', returncode=0)
        state = {'seeds': [record], 'configuration': campaign.configuration(args)}
        tracemalloc.start()
        try:
            result = campaign.aggregate(args, state)
            _, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        self.assertEqual(len(result.combined_points), 12000)
        self.assertLess(peak, 8 * 1024 * 1024)


class WrapperTests(unittest.TestCase):
    def test_cli_overrides_config_and_evo_is_explicit(self):
        args = run_next_scan.parse_args(['--campaign-dir', '/tmp/example', '--seed-start', '19',
            '--nseeds', '7', '--nrandom', '13', '--nrandom-count-evo-thc', '--jobs', 'auto', '--checkpoint-every', '5'])
        parsed = campaign.parse_args(run_next_scan.build_command(args)[2:])
        self.assertEqual((parsed.seed_start, parsed.nseeds, parsed.nrandom, parsed.jobs, parsed.checkpoint_every), (19, 7, 13, 'auto', 5))
        self.assertTrue(parsed.nrandom_count_evo_thc)
        self.assertTrue(parsed.run_ewpt)

    def test_resume_does_not_read_config_or_add_physics_overrides(self):
        args = run_next_scan.parse_args(['--campaign-dir', '/tmp/example', '--config', '/does/not/exist', '--resume', '--jobs', '3'])
        command = run_next_scan.build_command(args)
        self.assertIn('--resume', command)
        self.assertNotIn('--nrandom', command)
        self.assertNotIn('--generator-extra-arg', command)

    def test_old_config_defaults_to_raw_draws(self):
        args = run_next_scan.parse_args(['--campaign-dir', '/tmp/example'])
        parsed = campaign.parse_args(run_next_scan.build_command(args)[2:])
        self.assertFalse(parsed.nrandom_count_evo_thc)
        self.assertEqual(parsed.nrandom, 2500)

    def test_legacy_extra_count_flag_and_cli_override(self):
        with tempfile.TemporaryDirectory() as temp:
            config = Path(temp) / 'profile.json'
            data = json.loads((ROOT / 'config/next-scan-v2.json').read_text())
            data['generator_arguments'].append('--nrandom-count-evo-thc')
            config.write_text(json.dumps(data))
            for extra, expected in (([], True), (['--no-nrandom-count-evo-thc'], False)):
                args = run_next_scan.parse_args(['--config', str(config), '--campaign-dir', '/tmp/example', *extra])
                parsed = campaign.parse_args(run_next_scan.build_command(args)[2:])
                self.assertEqual(parsed.nrandom_count_evo_thc, expected)


class ResourceMeasurementTests(unittest.TestCase):
    def test_darwin_thread_table_and_group_rss(self):
        rows = '11 11 1024 /venv/bin/python\n12 11 2048 /build/bin/CalcTemps\n99 99 9999 unrelated\n'
        threads = 'USER PID TT CPU\nuser 12 ?? 50.0 R 12\nuser 12 ?? 50.0 R 12\n'
        with mock.patch.object(pilot_parallel_scan.sys, 'platform', 'darwin'), \
             mock.patch.object(pilot_parallel_scan.subprocess, 'check_output', return_value=rows), \
             mock.patch.object(pilot_parallel_scan.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0, threads, '')):
            rss, measured = pilot_parallel_scan.process_snapshot({11})
        self.assertEqual(rss, {11: 3072 * 1024})
        self.assertEqual(measured, {'CalcTemps': 2})

    def test_linux_thread_table(self):
        with mock.patch.object(pilot_parallel_scan.sys, 'platform', 'linux'), \
             mock.patch.object(pilot_parallel_scan.subprocess, 'check_output', return_value='11 11 1024 /bin/MinimaTracer\n'), \
             mock.patch.object(pilot_parallel_scan.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0, '11 1\n', '')):
            self.assertEqual(pilot_parallel_scan.process_snapshot({11})[1], {'MinimaTracer': 1})


if __name__ == '__main__':
    unittest.main()
