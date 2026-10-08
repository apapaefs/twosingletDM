import tempfile
import json
import os
import subprocess
import sys
import textwrap
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from generate_mg5_trsm_xsecs import ProcLocation, _survey_iterations, drive_mg, mg5_runtime_receipt
from mg5_process_runner import run_mg5_processes, selected_mg5_processes


class TestMG5ProcessRunner(unittest.TestCase):
    def test_requested_associated_production_processes_are_configured(self):
        self.assertEqual(ProcLocation["gg_heta0"], "gg_heta0/")
        self.assertEqual(ProcLocation["pp_eta0Z"], "pp_eta0Z/")

    def test_selected_processes_skips_unavailable_with_warning(self):
        buffer = StringIO()

        with redirect_stdout(buffer):
            selected = selected_mg5_processes(
                ["hh", "missing", "hhh"],
                proc_location={"hh": "gg_hh_twoscalar/", "hhh": "gg_hhh_twoscalar/"},
            )

        self.assertEqual(selected, ["hh", "hhh"])
        self.assertIn("Warning: MG5 process 'missing' is not available", buffer.getvalue())

    def test_run_mg5_processes_returns_cross_section_dictionary(self):
        calls = []

        def fake_get_xsec(process, run_tag, lambdas, k1, k2, k3, m2, w2, m3, w3, ecm):
            calls.append((process, run_tag, ecm))
            return {"hh": 1.25, "hhh": 0.031}[process]

        buffer = StringIO()
        with redirect_stdout(buffer):
            xsecs = run_mg5_processes(
                ["hh", "missing", "hhh"],
                "SCAN13.6",
                [1, 2, 3],
                0.9,
                0.1,
                0.0,
                300.0,
                1.0,
                500.0,
                2.0,
                13.6,
                get_xsec=fake_get_xsec,
                proc_location={"hh": "gg_hh_twoscalar/", "hhh": "gg_hhh_twoscalar/"},
            )

        self.assertEqual(xsecs, {"hh": 1.25, "hhh": 0.031})
        self.assertEqual(calls, [("hh", "SCAN13.6", 13.6), ("hhh", "SCAN13.6", 13.6)])
        self.assertIn("MG5 hh xsec [pb] = 1.25", buffer.getvalue())

    def test_run_mg5_processes_warns_when_selection_is_empty(self):
        buffer = StringIO()

        with redirect_stdout(buffer):
            xsecs = run_mg5_processes(
                [],
                "SCAN13.6",
                [1, 2, 3],
                0.9,
                0.1,
                0.0,
                300.0,
                1.0,
                500.0,
                2.0,
                13.6,
                get_xsec=lambda *args, **kwargs: self.fail("MG5 should not run"),
            )

        self.assertEqual(xsecs, {})
        self.assertIn("MG5ProcessesToRun is empty", buffer.getvalue())

    def test_optional_h1_width_and_k233_reach_madgraph_interface(self):
        calls = []

        def fake_get_xsec(*args, **kwargs):
            calls.append(kwargs)
            return 0.5

        xsecs = run_mg5_processes(
            ["gg_heta0"],
            "SCAN13.6-point4",
            list(range(10)),
            0.99,
            0.1,
            0.0,
            300.0,
            1.0,
            50.0,
            0.0,
            13.6,
            w1=0.0042,
            k233=25.0,
            get_xsec=fake_get_xsec,
        )

        self.assertEqual(xsecs, {"gg_heta0": 0.5})
        self.assertEqual(
            calls,
            [{"ecm": 13.6, "w1": 0.0042, "k233": 25.0}],
        )

    def test_madevent_command_sets_associated_production_parameters(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            mg_root = Path(tmpdir)
            process_dir = mg_root / "gg_heta0"
            madevent = process_dir / "bin" / "madevent"
            madevent.parent.mkdir(parents=True)
            madevent.write_text("", encoding="ascii")
            characteristics = process_dir / "SubProcesses/proc_characteristics"
            characteristics.parent.mkdir()
            characteristics.write_text("loop_induced = True\n", encoding="ascii")
            captured = {}

            def fake_run(args, **kwargs):
                command_file = Path(args[1])
                captured["text"] = command_file.read_text(encoding="ascii")
                run_name = next(line.split()[1] for line in captured["text"].splitlines()
                                if line.startswith("generate_events "))
                lhefile = (
                    process_dir
                    / "Events"
                    / run_name
                    / "unweighted_events.lhe.gz"
                )
                lhefile.parent.mkdir(parents=True)
                lhefile.write_bytes(b"placeholder")
                return SimpleNamespace(returncode=0, stdout="")

            lambdas = [[str(value) for value in range(1, 11)]]
            with patch.dict(os.environ, {"TRSM_MG5_CORES": "1"}), patch(
                "generate_mg5_trsm_xsecs.subprocess.run",
                side_effect=fake_run,
            ):
                count = drive_mg(
                    "gg_heta0",
                    "unit",
                    mg_root,
                    0.99,
                    0.1,
                    0.0,
                    lambdas,
                    300.0,
                    1.2,
                    50.0,
                    0.0,
                    1,
                    1,
                    ecm=13.6,
                    w1=0.0042,
                    k233=25.0,
                )

        self.assertEqual(count, 1)
        self.assertTrue(captured["text"].startswith("set nb_core 1\nset run_mode 0\n"))
        self.assertIn("--iterations=3", captured["text"])
        self.assertIn("set nhel 0\n", captured["text"])
        self.assertIn("set ebeam1 6800.0", captured["text"])
        self.assertIn("set Meta 300.0", captured["text"])
        self.assertIn("set Weta 1.2", captured["text"])
        self.assertIn("set Miota 50.0", captured["text"])
        self.assertIn("set Wiota 0.0", captured["text"])
        self.assertIn("set WH 0.0042", captured["text"])
        self.assertIn("set kap133 9", captured["text"])
        self.assertIn("set kap233 25.0", captured["text"])

    def test_tree_level_run_preserves_helicity_settings(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            process = root / "pp_eta0Z"
            madevent = process / "bin/madevent"
            madevent.parent.mkdir(parents=True)
            madevent.touch()
            characteristics = process / "SubProcesses/proc_characteristics"
            characteristics.parent.mkdir()
            characteristics.write_text("loop_induced = False\n")
            with patch("generate_mg5_trsm_xsecs.subprocess.run", return_value=SimpleNamespace(
                    returncode=1, stdout="stop after capturing commands")):
                with self.assertRaises(RuntimeError):
                    drive_mg("pp_eta0Z", "tree", root, .99, .1, 0.,
                             [[str(i) for i in range(10)]], 300., 1., 50., 0., 1, 1)
            command = (process / "mg5_pp_eta0Z_lambdavar_runtree.dcmd").read_text()
            self.assertIn("--iterations=3", command)
            self.assertNotIn("set nhel", command)

    def test_survey_meets_madevent_minimum_for_all_process_types(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            process_dir = Path(tmpdir)
            characteristics = process_dir / "SubProcesses" / "proc_characteristics"
            characteristics.parent.mkdir(parents=True)
            for content in (None, "loop_induced = True\n", "loop_induced = False\n", "unrelated = True\n"):
                with self.subTest(characteristics=content):
                    if content is not None:
                        characteristics.write_text(content, encoding="ascii")
                    self.assertGreaterEqual(_survey_iterations(process_dir), 3)

    def test_zero_survey_without_event_file_still_fails_instead_of_recording_zero(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            madevent = root / "gg_heta0/bin/madevent"
            madevent.parent.mkdir(parents=True)
            madevent.touch()
            with patch("generate_mg5_trsm_xsecs.subprocess.run", return_value=SimpleNamespace(
                    returncode=0, stdout="Survey return zero cross section.")):
                with self.assertRaisesRegex(RuntimeError, "Survey return zero cross section"):
                    drive_mg("gg_heta0", "zero", root, .99, .1, 0.,
                             [[str(i) for i in range(10)]], 300., 1., 50., 0., 1, 1)

    def test_preflight_checks_compiled_process_and_ignores_mutable_cards(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            files = {'VERSION': '3.5.15', 'pp_eta0Z/bin/madevent': '#!/bin/sh\n',
                     'pp_eta0Z/Cards/proc_card_mg5.dat': 'generate p p > eta0 z\n',
                     'pp_eta0Z/SubProcesses/subproc.mg': 'P0_test\n',
                     'pp_eta0Z/SubProcesses/P0_test/madevent': '#!/bin/sh\n',
                     'pp_eta0Z/bin/internal/ufomodel/parameters.py': 'parameters = []\n'}
            for name, content in files.items():
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(content)
                if path.name == 'madevent': path.chmod(0o755)
            receipt = mg5_runtime_receipt(['pp_eta0Z'], root)
            self.assertEqual(list(receipt['processes']), ['pp_eta0Z'])
            (root / 'pp_eta0Z/Cards/param_card.dat').write_text('parameters changed during run')
            self.assertEqual(receipt, mg5_runtime_receipt(['pp_eta0Z'], root))
            model = root / 'pp_eta0Z/bin/internal/ufomodel/parameters.py'
            model.write_text('different model')
            self.assertNotEqual(receipt, mg5_runtime_receipt(['pp_eta0Z'], root))
            (root / 'pp_eta0Z/SubProcesses/P0_test/madevent').unlink()
            with self.assertRaisesRegex(ValueError, 'compiled MG5 subprocess'):
                mg5_runtime_receipt(['pp_eta0Z'], root)
            with self.assertRaises(FileNotFoundError):
                mg5_runtime_receipt(['gg_heta0'], root)

    def test_preflight_rejects_missing_madloop_module_and_ignores_rebuilds_in_receipt(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            process = root / "gg_heta0"
            pv = process / "SubProcesses/PV5_test"
            files = {
                "VERSION": "3.5.15",
                "gg_heta0/bin/madevent": "#!/bin/sh\n",
                "gg_heta0/Cards/proc_card_mg5.dat": "generate g g > h eta0 [noborn=QCD]\n",
                "gg_heta0/SubProcesses/subproc.mg": "P0_test\n",
                "gg_heta0/SubProcesses/P0_test/madevent": "binary",
                "gg_heta0/SubProcesses/PV5_test/makefile": "makefile",
                "gg_heta0/SubProcesses/PV5_test/loop_matrix.f": "loop source",
                "gg_heta0/SubProcesses/PV5_test/check_sa.f": "check source",
                "gg_heta0/SubProcesses/PV5_test/polynomial.f": "      MODULE ML5_5_CONSTANTS\n",
                "gg_heta0/SubProcesses/PV5_test/check": "check binary",
                "gg_heta0/SubProcesses/ml5_5_constants.mod": "parent module",
            }
            for name, content in files.items():
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(content)
                path.chmod(0o755)
            with self.assertRaisesRegex(ValueError, "Missing MadLoop.*ml5_5_constants.mod"):
                mg5_runtime_receipt(["gg_heta0"], root)
            module = pv / "ml5_5_constants.mod"
            module.write_text("local module")
            receipt = mg5_runtime_receipt(["gg_heta0"], root)
            module.write_text("rebuilt module")
            (pv / "check").write_text("rebuilt check binary")
            self.assertEqual(receipt, mg5_runtime_receipt(["gg_heta0"], root))
            self.assertNotIn("madloop_artifacts", receipt["processes"]["gg_heta0"])

    def test_concurrent_workers_serialize_cards_and_keep_their_cross_sections(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            process = root / 'pp_eta0Z'
            madevent = process / 'bin/madevent'
            madevent.parent.mkdir(parents=True)
            madevent.write_text('#!' + sys.executable + '\n' + textwrap.dedent('''\
                import gzip, json, sys, time
                from pathlib import Path
                lines = Path(sys.argv[1]).read_text().splitlines()
                assert lines[:2] == ['set nb_core 1', 'set run_mode 0']
                run = next(line.split()[1] for line in lines if line.startswith('generate_events '))
                mass = next(line.split()[2] for line in lines if line.startswith('set Meta '))
                guard = Path('active')
                with guard.open('x') as stream:
                    stream.write(run)
                with Path('timeline.jsonl').open('a') as stream:
                    stream.write(json.dumps(['start', run]) + '\\n')
                card = Path('shared_param_card')
                card.write_text(mass)
                time.sleep(.25)
                assert card.read_text() == mass
                lhe = Path('Events') / run / 'unweighted_events.lhe.gz'
                lhe.parent.mkdir(parents=True)
                with gzip.open(lhe, 'wt') as stream:
                    stream.write('Integrated weight : ' + mass + '\\n')
                with Path('timeline.jsonl').open('a') as stream:
                    stream.write(json.dumps(['end', run]) + '\\n')
                guard.unlink()
                '''))
            madevent.chmod(0o755)
            code = ('from generate_mg5_trsm_xsecs import get_mg5_xsec; import sys; '
                    'value=get_mg5_xsec("pp_eta0Z", sys.argv[1], list(range(1,11)), '
                    '.99,.1,0,float(sys.argv[2]),1.,50.,0.,ecm=13.6,w1=.004,k233=2.); '
                    'print("RESULT", value)')
            env = dict(os.environ, TRSM_MG5_LOCATION=str(root), TRSM_MG5_CORES='1',
                       PYTHONPATH=str(Path(__file__).resolve().parent))
            children = [subprocess.Popen([sys.executable, '-B', '-c', code, tag, str(mass)],
                        env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                        for tag, mass in [('first', 300), ('second', 400)]]
            try:
                for child, mass in zip(children, (300, 400)):
                    output, _ = child.communicate(timeout=15)
                    self.assertEqual(child.returncode, 0, output)
                    self.assertIn(f'RESULT {mass}.0', output)
            finally:
                for child in children:
                    if child.poll() is None:
                        child.kill()
                    child.wait()
            timeline = [json.loads(line) for line in (process / 'timeline.jsonl').read_text().splitlines()]
            self.assertEqual([row[0] for row in timeline], ['start', 'end', 'start', 'end'])
            self.assertEqual(timeline[0][1], timeline[1][1])
            self.assertEqual(timeline[2][1], timeline[3][1])


if __name__ == "__main__":
    unittest.main()
