"""Exercise the installation and historical collaborator entry points."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

EXAMPLE = Path(__file__).resolve().parents[1]
ROOT = EXAMPLE.parents[1]
sys.path.insert(0, str(EXAMPLE))
from cutflow import stages
from test_run_single_point import RAW
import run_single_point as single


class CompatibilityTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.env = dict(os.environ)
        self.env.pop('PYTHONPATH', None)
        self.env.pop('TRSM_REPO_ROOT', None)
        self.env['MPLCONFIGDIR'] = str(self.directory / 'cache')

    def python(self, script, *args, cwd=None, env=None):
        return subprocess.run([sys.executable, str(EXAMPLE / script), *map(str, args)],
                              cwd=cwd or self.directory, env=env or self.env,
                              capture_output=True, text=True, timeout=30)

    def test_scripts_import_without_pythonpath_from_each_historical_directory(self):
        for script in ('run_single_point.py', 'run_scan.py', 'generate_oks.py',
                       'source/write_mo.py', 'source/mO_excluder.py', 'validate_cmb.py',
                       'check_installation.py'):
            with self.subTest(script=script):
                result = self.python(script, '--help', cwd=(EXAMPLE / script).parent)
                self.assertEqual(result.returncode, 0, result.stderr)
        result = self.python('check_installation.py', '--no-native')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(str(ROOT / 'test_trsm_DM.py'), result.stdout)

    def test_missing_checkout_has_actionable_diagnostic_and_explicit_path_works(self):
        detached = self.directory / 'example'
        detached.mkdir()
        for name in ('run_single_point.py', '_bootstrap.py'):
            shutil.copy2(EXAMPLE / name, detached)
        command = [sys.executable, str(detached / 'run_single_point.py'), '--help']
        result = subprocess.run(command, cwd=self.directory, env=self.env, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('TRSM_REPO_ROOT', result.stderr)
        self.assertNotIn('Traceback', result.stderr)
        result = subprocess.run(command, cwd=self.directory, env=dict(self.env, TRSM_REPO_ROOT=str(ROOT)), capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_sourcing_setup_errors_does_not_exit_or_change_caller_state(self):
        for shell in ('bash', 'zsh'):
            if not shutil.which(shell):
                continue
            for arguments in ('', ' /not/a/micromegas/directory'):
                with self.subTest(shell=shell, args=arguments):
                    script = 'old_pwd=$PWD; old_options=$-; source "$1"' + arguments + '\n' + (
                        'rc=$?; test "$rc" = 2 || exit 9; '
                        'test "$PWD" = "$old_pwd" || exit 10; '
                        'test "$-" = "$old_options" || exit 11; echo shell-survived')
                    result = subprocess.run([shell, '-c', script, shell, str(ROOT / 'DM/setup_micromegas.sh')],
                                            cwd=self.directory, capture_output=True, text=True)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertIn('shell-survived', result.stdout)

    def test_sourced_setup_resolves_its_own_path_and_rejects_existing_install(self):
        install = self.directory / 'micromegas_7.1.4'
        (install / 'TRSM').mkdir(parents=True)
        sentinel = install / 'TRSM/keep'
        sentinel.write_text('existing installation')
        for shell in ('sh', 'bash', 'zsh'):
            if not shutil.which(shell):
                continue
            command = [shell, str(ROOT / 'DM/setup_micromegas.sh'), str(install)]
            if shell != 'sh':
                command = [shell, '-c', 'source "$1" "$2"', shell, command[1], command[2]]
            result = subprocess.run(command, cwd=self.directory, capture_output=True, text=True)
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertIn('Refusing to overwrite', result.stderr)
            self.assertEqual(sentinel.read_text(), 'existing installation')

    def test_positional_excluder_preserves_historical_arithmetic_and_files(self):
        points = self.directory / 'oks.dat'
        points.write_text('7 .1 .02 .03 50 300 .1 200\n')
        output = self.directory / 'legacy'
        env = dict(self.env, MO_EXCLUDER_OKS_FILE=str(points), MO_EXCLUDER_OUTPUT_DIR=str(output))
        result = self.python('source/mO_excluder.py', '7', '50', '.06', '1e-15', '50', '1e-20', env=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        values = list(map(float, (output / 'DM_data/DM_data_7').read_text().split()))
        self.assertEqual(len(values), 18)
        self.assertEqual((output / 'DM_data/DM_data').read_bytes(), (output / 'DM_data/DM_data_7').read_bytes())
        self.assertAlmostEqual(values[11] / values[12], .121 / .06)
        self.assertTrue((output / 'allall.dat').read_text())
        self.assertEqual((output / 'dmexcl.dat').read_text(), '')
        self.assertIn('solver status and CMB were not supplied', result.stdout)

    def test_original_card_shell_workflow_and_nested_raw_replay(self):
        points = self.directory / 'oks.dat'
        points.write_text('1 .1 .05 .15 50 200 -.14943813247359922 380\n')
        cards = self.directory / 'cards'
        for _ in range(2):
            result = self.python('source/write_mo.py', '--input', points, '--output-dir', cards)
            self.assertEqual(result.returncode, 0, result.stderr)
        raw = self.directory / 'raw/OUT_mO'
        raw.mkdir(parents=True)
        (raw / 'OUT_mO_1').write_text(RAW)
        output = self.directory / 'out'
        result = subprocess.run(['bash', str(EXAMPLE / 'run/MOrun.sh'), '--raw-output-dir', str(raw.parent), '--no-planck-cmb'],
                                cwd=EXAMPLE / 'run', env=dict(self.env, PYTHON=sys.executable,
                                    CARD_DIR=str(cards), OUTPUT_DIR=str(output)), capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('Relic density: pass', (output / 'MOrun.log').read_text())
        self.assertEqual((output / 'OUT_mO/OUT_mO_1').read_text(), RAW)
        self.assertEqual(len((output / 'DM_data/DM_data_1').read_text().split()), 18)
        self.assertEqual(len(json.loads((output / 'cutflow.json').read_text())), 4)
        self.assertEqual(json.loads((output / 'results.json').read_text())[0]['dm_cmb_status'], 'disabled')

    def test_card_regeneration_preserves_previous_cards_and_excludes_stale_indices(self):
        points, cards = self.directory / 'oks.dat', self.directory / 'cards'
        points.write_text('1 .1 .02 .03 50 300 .1 200\n2 .1 .02 .03 60 300 .1 200\n')
        result = self.python('source/write_mo.py', '--input', points, '--output-dir', cards)
        self.assertEqual(result.returncode, 0, result.stderr)
        old = (cards / 'MO_inp1.dat').read_bytes()
        points.write_text('1 .1 .02 .03 55 300 .1 200\n')
        result = self.python('source/write_mo.py', '--input', points, '--output-dir', cards)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([p.name for p in cards.glob('MO_inp*.dat')], ['MO_inp1.dat'])
        backup, = cards.glob('.previous-cards-*')
        self.assertEqual((backup / 'MO_inp1.dat').read_bytes(), old)
        self.assertTrue((backup / 'MO_inp2.dat').is_file())
        self.assertEqual(single.read_card(cards / 'MO_inp1.dat').mx, 55)

    def test_cutflow_keeps_unknown_and_disabled_constraints_distinct(self):
        def row(index, **overrides):
            return dict(index=index, dm_calculation_status='success', dm_relic_excluded=False,
                        dm_direct_detection_excluded=False, dm_indirect_detection_excluded=False,
                        dm_cmb_enabled=True, dm_cmb_excluded=False) | overrides
        rows = [row(1), row(2, dm_relic_excluded=True), row(3, dm_direct_detection_excluded=True),
                row(4, dm_indirect_detection_excluded=True), row(5, dm_cmb_excluded=True),
                row(6, dm_cmb_excluded=None), row(7, dm_cmb_enabled=False, dm_cmb_excluded=None)]
        steps = list(stages(rows))
        self.assertEqual([step['passed'] for step in steps], [7, 6, 5, 4, 2])
        self.assertEqual(steps[-1]['unassessed'], 1)
        self.assertEqual(steps[-1]['excluded'], 4)
        self.assertEqual(steps[-1]['indices'], [1, 7])

    def test_plot_driver_handles_empty_selections(self):
        import importlib.util
        if not importlib.util.find_spec('matplotlib') or not importlib.util.find_spec('numpy'):
            self.skipTest('Optional plotting dependencies are not installed')
        points = self.directory / 'oks.dat'
        points.write_text('1 .1 .05 .15 50 200 -.14943813247359922 380\n')
        raw = self.directory / 'raw'
        raw.mkdir()
        (raw / 'OUT_mO_1').write_text(RAW.replace('Omega=0.06', 'Omega=0.2'))
        output = self.directory / 'out'
        result = self.python('run_scan.py', '--input', points, '--raw-output-dir', raw,
                             '--output-dir', output, '--no-planck-cmb')
        self.assertEqual(result.returncode, 0, result.stderr)
        # Keep this test's historical outplots away from user plots in the repo.
        shutil.copytree(EXAMPLE / 'plot', self.directory / 'plot', ignore=shutil.ignore_patterns('outplots', '__pycache__'))
        shutil.copy2(EXAMPLE / 'cutflow.py', self.directory)
        result = subprocess.run([sys.executable, str(self.directory / 'plot/steer_plots.py'),
                                 '--outdir', str(output), '--xvar', 'MX', '--yvar', 'LHX'],
                                cwd=self.directory, env=self.env, capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue((output / 'constraint_cutflow.png').is_file())
        self.assertFalse((output / 'cmb_ratios.png').exists())
        self.assertIn('no points in this subset', result.stdout)

    def test_relic_upper_limit_override_does_not_change_abundance_reference(self):
        point = single.PointInput(1, .1, .05, .15, 50., 200., -.14943813247359922, 380.)
        card = self.directory / 'MO_inp1.dat'
        single.write_card(point, card)
        args = single.configure(single.parse_args(['--card', str(card), '--relic-upper-limit', '.123']), replay=True)
        row = single.evaluate_point(point, args, self.directory / 'custom',
                                    raw_output=RAW.replace('Omega=0.06', 'Omega=0.122'))
        self.assertFalse(row['dm_relic_excluded'])
        self.assertEqual(row['dm_relic_upper_limit'], .123)
        self.assertEqual(row['dm_cmb_abundance_fraction'], 1)


if __name__ == '__main__':
    unittest.main()
