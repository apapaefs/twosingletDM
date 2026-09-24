"""Regression tests for the collaborator-facing standalone DM workflow."""
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

EXAMPLE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXAMPLE))
import run_single_point as single
import run_scan
from test_trsm_DM import test_dm as production_dm
from trsm_inputs import ABUNDANCE_REFERENCE, RELIC_UPPER_LIMIT, micromegas_sm_inputs

RAW = '''Masses of odd sector Particles:
~X : MX = 50 ||
==== Calculation of relic density =====
Xf=25 Omega=0.06 darkOmega_error=0
~X[~X]-nucleon cross sections[pb]:
 proton SI 1e-15 [1e-15] SD 0 [0]
 neutron SI 1e-15 [1e-15] SD 0 [0]

TRSM_inputs_v2 {"MX":50,"Mh":125.09,"Mh2":380,"width_h1":0.004,"width_h2":0.1}
TRSM_loop_hook_v2 relic_calls=30 indirect_calls=2
FermiLAT_line_channel A A: E_gamma=50[GeV], sigmaV=1e-30[cm^3 s^-1], N_gamma*sigmaV=2e-30[cm^3 s^-1], Phi_R16=1e-20[cm^-2 s^-1]
TRSM_PlanckCMB_v1 {"status":"ok","ratio_raw":8,"sigma_v_cm3_s":1e-26}
'''


class ExampleTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.point = single.PointInput(1, .1, .05, .15, 50., 200., -.14943813247359922, 380.)
        self.card = self.root / 'MO_inp1.dat'
        single.write_card(self.point, self.card)

    def args(self, *extra):
        args = single.parse_args(['--card', str(self.card), *extra])
        return single.configure(args, replay=True)

    def evaluate(self, raw=RAW, *extra, name='out'):
        return single.evaluate_point(self.point, self.args(*extra), self.root / name,
                                     raw_output=raw)

    def test_defaults_and_precision(self):
        args = self.args()
        self.assertEqual(args.micromegas_version, '7.1.4')
        self.assertTrue(args.planck_cmb)
        self.assertEqual(RELIC_UPPER_LIMIT, .121)
        values = dict(line.split() for line in self.card.read_text().splitlines())
        for key, value in micromegas_sm_inputs().items():
            self.assertEqual(float(values[key]), value)
        restored = single.read_card(self.card)
        self.assertAlmostEqual(restored.sint, self.point.sint, places=16)
        self.assertFalse(self.args('--micromegas-version', '6').planck_cmb)
        self.assertTrue(self.args('--micromegas-version', '6', '--planck-cmb').planck_cmb)

    def test_explicit_sm_conflict_is_not_silently_ignored(self):
        self.card.write_text(self.card.read_text().replace('125.09', '125.1'))
        with self.assertRaisesRegex(ValueError, 'shared v2'):
            single.read_card(self.card)

    def test_cmb_rescaling_independent_arithmetic_and_pipeline_agreement(self):
        row = self.evaluate()
        self.assertEqual(row['dm_cmb_ratio_raw'], 8)
        self.assertEqual(row['dm_cmb_abundance_fraction'], .5)
        self.assertEqual(row['dm_cmb_ratio'], 8 * (.06 / .12) ** 2)
        self.assertTrue(row['dm_cmb_excluded'])
        self.assertFalse(row['dm_passed'])
        self.assertEqual(row['dm_freezeout_temperature_GeV'], 2.)
        passed, _, diagnostics = production_dm(**single.asdict(self.point.dm_point()),
                                              raw_output=RAW, planck_cmb=True, limit_table=self.args().si_table)
        for key, value in single.json_safe(diagnostics).items():
            self.assertEqual(row[key], value, key)
        self.assertEqual(row['dm_passed'], passed)

    def test_strict_cmb_boundary_and_abundance_cap(self):
        for i, (omega, ratio, expected, excluded) in enumerate([
            (.06, 4., 1., False), (.06, 4.00001, 1.0000025, True),
            (.24, .5, .5, False), (0., 8., 0., False)]):
            raw = RAW.replace('Omega=0.06', f'Omega={omega}').replace('"ratio_raw":8', f'"ratio_raw":{ratio}')
            row = self.evaluate(raw, name=f'case{i}')
            self.assertAlmostEqual(row['dm_cmb_ratio'], expected)
            self.assertEqual(row['dm_cmb_excluded'], excluded)

    def test_no_rescaling_applies_to_all_signals(self):
        scaled = self.evaluate(name='scaled')
        full = self.evaluate(RAW, '--no-rescale', name='full')
        self.assertEqual(full['dm_cmb_ratio'], 8)
        self.assertEqual(full['dm_cmb_abundance_fraction'], 1)
        self.assertAlmostEqual(scaled['dm_dir_det_limit'], 2 * full['dm_dir_det_limit'])
        self.assertAlmostEqual(scaled['dm_indirect_ratio'], .25 * full['dm_indirect_ratio'])

    def test_cmb_disabled_missing_and_error_are_distinct(self):
        no_cmb = '\n'.join(line for line in RAW.splitlines() if not line.startswith('TRSM_PlanckCMB'))
        missing = self.evaluate(no_cmb, name='missing')
        self.assertEqual(missing['dm_cmb_status'], 'missing_output')
        self.assertIsNone(missing['dm_cmb_excluded'])
        self.assertIsNone(missing['dm_passed'])
        disabled = self.evaluate(no_cmb, '--no-planck-cmb', name='disabled')
        self.assertEqual(disabled['dm_cmb_status'], 'disabled')
        self.assertTrue(disabled['dm_passed'])
        error = self.evaluate(no_cmb + '\nTRSM_PlanckCMB_v1 {"status":"error","reason":"calcSpectrum_error"}\n', name='error')
        self.assertEqual(error['dm_cmb_reason'], 'calcSpectrum_error')
        self.assertIsNone(error['dm_cmb_excluded'])

    def test_solver_failure_and_legacy_logs_remain_unassessed(self):
        for name, raw in [('failed', RAW.replace('darkOmega_error=0', 'darkOmega_error=7')),
                          ('legacy', RAW.replace('darkOmega_error=0', '')),
                          ('invalid', RAW.replace('Omega=0.06', 'Omega=nan'))]:
            row = self.evaluate(raw, name=name)
            self.assertIsNone(row['dm_passed'])
            self.assertIsNone(row['dm_freezeout_temperature_GeV'])
            self.assertIsNone(row['dm_cmb_excluded'])

    def test_native_failure_cannot_pass_using_partial_output(self):
        class FailedRunner:
            def run(self, card):
                return RAW, 'native error', {'returncode': 9, 'error': 'failed'}
        row = single.evaluate_point(self.point, self.args(), self.root / 'failed', runner=FailedRunner())
        self.assertIsNone(row['dm_passed'])
        self.assertIsNone(row['dm_cmb_excluded'])
        self.assertEqual((self.root / 'failed/ERR_mO_1').read_text(), 'native error')

    def test_native_command_enables_cmb_and_uses_isolated_directory(self):
        args = self.args()
        fake = subprocess.CompletedProcess([], 0, RAW, 'warning')
        with single.NativeRunner(args) as runner, patch.object(single.subprocess, 'run', return_value=fake) as run:
            raw, stderr, invocation = runner.run(self.card)
            call = run.call_args
            self.assertEqual(call.args[0][-1], '--planck-cmb')
            self.assertEqual(call.kwargs['cwd'], call.kwargs['env']['TRSM_RUNTIME_DIR'])
            self.assertNotEqual(call.kwargs['cwd'], str(self.root))
            self.assertEqual(stderr, 'warning')
            self.assertEqual(invocation['returncode'], 0)

    def test_timeout_retains_partial_stdout_and_stderr(self):
        args = self.args()
        timeout = subprocess.TimeoutExpired('main', 600, output=b'partial result', stderr=b'failure')
        with single.NativeRunner(args) as runner, patch.object(single.subprocess, 'run', side_effect=timeout):
            raw, stderr, invocation = runner.run(self.card)
        self.assertEqual(raw, 'partial result')
        self.assertEqual(stderr, 'failure')
        self.assertIn('exceeded', invocation['error'])

    def test_dd_coverage_and_relic_boundary(self):
        # A valid low-mass CMB result must not invent DD table coverage.
        low = RAW.replace('MX = 50', 'MX = 5').replace('"MX":50', '"MX":5').replace('"ratio_raw":8', '"ratio_raw":0.1')
        row = self.evaluate(low, name='low')
        self.assertIsNone(row['dm_direct_detection_excluded'])
        self.assertFalse(row['dm_cmb_excluded'])
        self.assertIsNone(row['dm_passed'])
        for i, (omega, excluded) in enumerate([(.121, False), (.121000000001, True)]):
            row = self.evaluate(RAW.replace('Omega=0.06', f'Omega={omega}'), name=f'omega{i}')
            self.assertEqual(row['dm_relic_excluded'], excluded)

    def test_serialization_column_alignment_and_cmb_subsets(self):
        row = self.evaluate()
        run_scan.aggregate(self.root / 'out', [row])
        saved = json.loads((self.root / 'out/result_1.json').read_text())
        self.assertEqual(saved['point'], row)
        with (self.root / 'out/results.tsv').open() as stream:
            table = list(csv.DictReader(stream, delimiter='\t'))
        self.assertEqual(table[0]['dm_cmb_excluded'], '1')
        self.assertEqual(len((self.root / 'out/scan_results.dat').read_text().split()), 11)
        self.assertEqual(len((self.root / 'out/cmbexcl.dat').read_text().split()), 18)
        counts = json.loads((self.root / 'out/counts.json').read_text())
        self.assertEqual(counts['cmb_caughtit'], 1)
        self.assertEqual(counts['allall'] + counts['dmexcl'] + counts['dmunassessed'], 1)
        with self.assertRaisesRegex(ValueError, 'already exists'):
            self.evaluate()

    def test_batch_cli_replays_complete_logs_and_preserves_inputs(self):
        rawdir = self.root / 'raw'
        rawdir.mkdir()
        (rawdir / 'OUT_mO_1').write_text(RAW)
        oks = self.root / 'oks.dat'
        oks.write_text('# header\n1 .1 .05 .15 50 200 -.14943813247359922 380\n')
        out = self.root / 'batch'
        cmd = [sys.executable, str(EXAMPLE / 'run_scan.py'), '--input', str(oks),
               '--raw-output-dir', str(rawdir), '--output-dir', str(out)]
        result = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(json.loads((out / 'results.json').read_text())[0]['dm_cmb_excluded'])
        self.assertEqual(single.read_points(out / 'oks.dat'), [self.point])
        self.assertNotEqual(subprocess.run(cmd, capture_output=True).returncode, 0)

    def test_duplicate_indices_rejected(self):
        oks = self.root / 'duplicate.dat'
        oks.write_text('1 .1 .05 .15 50 200 .1 380\n' * 2)
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            single.read_points(oks)


if __name__ == '__main__':
    unittest.main()
