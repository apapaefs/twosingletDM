"""Automatic CalcGW gating, checkpoint recovery and scan compatibility."""
import contextlib
import io
import json
import os
from pathlib import Path
import socket
import sys
import tempfile
import unittest
from unittest import mock

import run_trsm_seed_campaign as campaign
from reprocess_trsm_gw import reprocess
from test_trsm_gw import native_row, write_scan
from tools import run_campaign_calcgw as followup
from tools import run_next_scan
from trsm_scan_campaign import CampaignLock, CampaignStateError

ROOT = Path(__file__).resolve().parent


class CampaignCalcGWTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.state = dict(schema='trsm_seed_campaign_v1', hostname=socket.gethostname(),
                          status='complete', seeds=[dict(seed=1, status='complete')],
                          configuration=dict(python_executable=sys.executable, ewpt_thigh=1000.),
                          execution=dict(jobs=2))
        self.save_state()
        write_scan(self.directory / 'combined_points.tsv')

    def save_state(self):
        (self.directory / 'campaign_state.json').write_text(json.dumps(self.state))

    def test_profiles_enable_followup_and_allow_explicit_disable(self):
        for name in ('odysseus-m2-m3-4-1000.json', 'odysseus-m2-4-200-m3-4-65.json'):
            options = ['--config', str(ROOT / 'config' / name), '--campaign-dir', str(self.directory)]
            for extra, enabled in (([], True), (['--no-run-calcgw'], False)):
                args = run_next_scan.parse_args(options + extra + ['--calcgw-jobs', '7'])
                parsed = campaign.parse_args(run_next_scan.build_command(args)[2:])
                self.assertEqual(parsed.run_calcgw, enabled)
                self.assertEqual(parsed.calcgw_jobs, 7)

    def test_followup_is_operational_and_survives_resume_without_config(self):
        options = ['--campaign-dir', str(self.directory), '--seed-start', '1', '--nseeds', '1']
        original = campaign.parse_args(options)
        enabled = campaign.parse_args(options + ['--run-calcgw', '--calcgw-jobs', '7'])
        self.assertEqual(campaign.configuration(original), campaign.configuration(enabled))
        saved = dict(configuration=campaign.configuration(original),
                     execution=dict(run_calcgw=True, calcgw_jobs=7))
        args = run_next_scan.parse_args(['--campaign-dir', str(self.directory), '--resume',
                                        '--config', '/missing/profile.json'])
        resumed = campaign.parse_args(run_next_scan.build_command(args)[2:])
        campaign.hydrate_configuration(resumed, saved)
        self.assertTrue(resumed.run_calcgw)
        self.assertEqual(resumed.calcgw_jobs, 7)
        disabled = campaign.parse_args(['--campaign-dir', str(self.directory), '--resume', '--no-run-calcgw'])
        campaign.hydrate_configuration(disabled, saved)
        self.assertFalse(disabled.run_calcgw)

    def test_scan_failure_and_aggregate_only_never_launch_calcgw(self):
        for status, extra, expected in (('failed', [], 1), ('complete', ['--aggregate-only'], 0)):
            args = ['--campaign-dir', str(self.directory), '--run-calcgw', *extra]
            if not extra:
                args += ['--seed-start', '1', '--nseeds', '1']
            result = mock.Mock(seed_results=[mock.Mock(status=status)])
            with mock.patch.object(campaign, 'run_campaign', return_value=result), \
                 mock.patch.object(followup, 'run_followup') as run:
                self.assertEqual(campaign.main(args), expected)
                run.assert_not_called()

    def test_success_launches_followup_and_propagates_its_exit_code(self):
        result = mock.Mock(seed_results=[mock.Mock(status='complete')])
        with mock.patch.object(campaign, 'run_campaign', return_value=result), \
             mock.patch.object(followup, 'run_followup', return_value=7) as run:
            self.assertEqual(campaign.main(['--campaign-dir', str(self.directory), '--seed-start', '1',
                                           '--nseeds', '1', '--run-calcgw', '--calcgw-jobs', '3']), 7)
            run.assert_called_once_with(self.directory, jobs=3)

    def test_wait_does_not_use_output_until_supervisor_releases_lock(self):
        lock = self.directory / 'campaign_state.lock'
        lock.write_text(json.dumps(dict(hostname=socket.gethostname(), pid=os.getpid())))
        with mock.patch.object(followup.time, 'sleep', side_effect=lambda _: lock.unlink()) as sleep:
            self.assertEqual(followup.completed_scan(self.directory, wait=True), self.state)
        sleep.assert_called_once()

    def test_incomplete_or_missing_output_never_starts_calcgw(self):
        for status in ('running', 'incomplete'):
            self.state['status'] = status
            self.save_state()
            with mock.patch.object(followup, 'run_command') as run, \
                 self.assertRaisesRegex(ValueError, 'not complete'):
                followup.run_followup(self.directory, wait=True)
            run.assert_not_called()
        self.state['status'] = 'complete'
        self.save_state()
        (self.directory / 'combined_points.tsv').unlink()
        with self.assertRaisesRegex(ValueError, 'no combined_points'):
            followup.completed_scan(self.directory)

    def test_completion_during_status_check_is_not_reported_as_interruption(self):
        self.state['status'] = 'running'
        self.save_state()
        lock = self.directory / 'campaign_state.lock'
        lock.write_text(json.dumps(dict(hostname=socket.gethostname(), pid=os.getpid())))

        def finish_scan(pid):
            self.state['status'] = 'complete'
            self.save_state()
            lock.unlink()
            return False

        with mock.patch.object(followup, '_pid_is_alive', side_effect=finish_scan):
            self.assertEqual(followup.completed_scan(self.directory, wait=True)['status'], 'complete')

    def test_duplicate_followup_cannot_start_a_second_writer(self):
        with CampaignLock(self.directory / 'calcgw_followup.json', ['test']):
            with mock.patch.object(followup, 'run_command') as run, self.assertRaises(CampaignStateError):
                followup.run_followup(self.directory)
            run.assert_not_called()

    def test_actual_subprocess_resumes_pilot_without_recomputing_points(self):
        executable = self.directory / 'CalcGW'
        calls = self.directory / 'native_calls.txt'
        executable.write_text(
            '#!' + sys.executable + '\nimport csv, json, sys\nfrom pathlib import Path\n'
            'args = dict(arg[2:].split("=", 1) for arg in sys.argv[1:])\n'
            'with open(args["input"]) as f:\n'
            '    rows = list(csv.reader(f, delimiter="\\t"))\n'
            'index = int(rows[1][0])\n'
            'row = json.loads(' + repr(json.dumps(native_row())) + ')\n'
            'row["index"] = index\n'
            'with open(args["output"], "w") as f:\n'
            '    w = csv.DictWriter(f, fieldnames=list(row), delimiter="\\t")\n'
            '    w.writeheader(); w.writerow(row)\n'
            'with open(' + repr(str(calls)) + ', "a") as f: f.write(str(index) + "\\n")\n')
        executable.chmod(0o755)
        self.state['configuration']['ewpt_executable'] = str(executable.with_name('CalcTemps'))
        self.save_state()
        source = self.directory / 'combined_points.tsv'
        output = self.directory / 'combined_points_gw.tsv'
        original = source.read_bytes()
        with contextlib.redirect_stdout(io.StringIO()):
            pilot = reprocess(source, output, executable=executable, thigh=1000., max_points=1)
            self.assertFalse(pilot['complete'])
            self.assertEqual(followup.run_followup(self.directory, jobs=2), 0)
            self.assertEqual(followup.run_followup(self.directory, jobs=2), 0)
        self.assertEqual(sorted(calls.read_text().splitlines()), ['1', '3', '4'])
        self.assertEqual(source.read_bytes(), original)
        progress = json.loads((self.directory / 'calcgw_followup.json').read_text())
        self.assertEqual(progress['status'], 'complete')
        self.assertEqual(progress['counts']['processed'], 3)
        self.assertIn('--resume', progress['command'])


if __name__ == '__main__':
    unittest.main()
