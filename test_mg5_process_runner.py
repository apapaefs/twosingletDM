import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from generate_mg5_trsm_xsecs import ProcLocation, _survey_iterations, drive_mg
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
            captured = {}

            def fake_run(args, **kwargs):
                command_file = Path(args[1])
                captured["text"] = command_file.read_text(encoding="ascii")
                run_name = captured["text"].splitlines()[0].split()[1]
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
            with patch(
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
        self.assertIn("--iterations=3", captured["text"])
        self.assertIn("set ebeam1 6800.0", captured["text"])
        self.assertIn("set Meta 300.0", captured["text"])
        self.assertIn("set Weta 1.2", captured["text"])
        self.assertIn("set Miota 50.0", captured["text"])
        self.assertIn("set Wiota 0.0", captured["text"])
        self.assertIn("set WH 0.0042", captured["text"])
        self.assertIn("set kap133 9", captured["text"])
        self.assertIn("set kap233 25.0", captured["text"])

    def test_loop_induced_process_keeps_single_fast_survey_iteration(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            process_dir = Path(tmpdir)
            characteristics = process_dir / "SubProcesses" / "proc_characteristics"
            characteristics.parent.mkdir(parents=True)
            characteristics.write_text("loop_induced = True\n", encoding="ascii")

            self.assertEqual(_survey_iterations(process_dir), 1)


if __name__ == "__main__":
    unittest.main()
