import csv
import gc
import importlib.util
import math
import os
import subprocess
import shutil
import tempfile
from pathlib import Path
import unittest


SCRIPT_PATH = Path(__file__).resolve().parent / "test_trsm_ewpt.py"
BSMPT_FIXTURE = Path(__file__).resolve().parents[1] / "BSMPT" / "test.output_new.csv"


def load_module():
    spec = importlib.util.spec_from_file_location("test_trsm_ewpt", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestTRSMEWPT(unittest.TestCase):
    def test_writes_trsm_input_with_expected_header_and_values(self):
        ewpt = load_module()
        point = ewpt.TRSMEWPTPoint(
            index=7,
            m1=125.09,
            m2=300.0,
            m3=400.0,
            vs=200.0,
            a12=0.2,
            lx=0.1,
            lphix=0.05,
            lsx=0.05,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            infile = Path(tmpdir) / "TRSM_Input.tsv"
            ewpt.write_trsm_input(infile, point)
            with infile.open(encoding="ascii") as stream:
                rows = list(csv.reader(stream, delimiter="\t"))

        self.assertEqual(rows[0], ["", "m1", "m2", "m3", "vs", "a12", "lx", "lphix", "lsx"])
        self.assertEqual(rows[1], ["7", "125.09", "300.0", "400.0", "200.0", "0.2", "0.1", "0.05", "0.05"])

    def test_builds_calctemps_command_from_config(self):
        ewpt = load_module()
        config = ewpt.EWPTConfig(
            executable=Path("/tmp/CalcTemps"),
            multistepmode="2",
            use_multithreading=True,
        )

        command = ewpt.build_calctemps_command(
            config,
            input_path=Path("/tmp/input.tsv"),
            output_path=Path("/tmp/output.csv"),
        )

        self.assertEqual(command[:6], ["/tmp/CalcTemps", "trsm", "/tmp/input.tsv", "/tmp/output.csv", "2", "2"])
        self.assertIn("--multistepmode=2", command)
        self.assertIn("--usemultithreading=true", command)

    def test_parses_calctemps_output_with_numeric_coercion_and_history(self):
        ewpt = load_module()

        rows = ewpt.parse_calctemps_output(BSMPT_FIXTURE)

        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["index"], 1.0)
        self.assertTrue(math.isclose(row["m2"], 300.0))
        self.assertTrue(math.isclose(row["T_crit_0"], 162.05976502096183))
        self.assertEqual(row["status_crit_0"], "success")
        self.assertEqual(row["transition_history"], "0-(0)->1")

    def test_eq_4_18_check_passes_for_bsmpt_fixture(self):
        ewpt = load_module()
        row = ewpt.parse_calctemps_output(BSMPT_FIXTURE)[0]

        check = ewpt.check_eq_4_18(row)

        self.assertTrue(check.satisfied)
        self.assertTrue(check.conditions["c_phiphi > 0"])
        self.assertTrue(check.conditions["C_phix > 0"])
        self.assertTrue(check.conditions["D > 0"])
        self.assertTrue(math.isclose(check.couplings["c_phiphi"], 0.0766285972235586))

    def test_eq_4_18_check_fails_when_coupling_matrix_is_not_positive(self):
        ewpt = load_module()
        row = {
            "m1": 125.09,
            "m2": 300.0,
            "m3": 400.0,
            "vs": 200.0,
            "a12": 0.2,
            "lx": -0.1,
            "lphix": 0.05,
            "lsx": 0.05,
        }

        check = ewpt.check_eq_4_18(row)

        self.assertFalse(check.satisfied)
        self.assertFalse(check.conditions["c_xx > 0"])
        self.assertFalse(check.conditions["C_phix > 0"])

    def test_summary_prints_eq_4_18_status(self):
        ewpt = load_module()
        row = ewpt.parse_calctemps_output(BSMPT_FIXTURE)[0]
        result = ewpt.EWPTResult(
            rows=[row],
            input_path=Path("/tmp/TRSM_Input.tsv"),
            output_path=Path("/tmp/test.output_new.csv"),
            command=[],
            returncode=0,
        )

        summary = ewpt.summarize_result(result)

        self.assertIn("eq_4_18_satisfied: True", summary)
        self.assertIn("D > 0=True", summary)

    def test_run_trsm_ewpt_uses_fake_runner_and_parses_output(self):
        ewpt = load_module()
        point = ewpt.TRSMEWPTPoint(m2=300.0, m3=400.0, vs=200.0, a12=0.2, lx=0.1, lphix=0.05, lsx=0.05)
        calls = []

        def fake_runner(command, check, capture_output, text, cwd):
            calls.append((command, check, capture_output, text, cwd))
            output_path = Path(command[3])
            output_path.write_text(
                "\t".join(["", "m1", "m2", "status_crit_0", "T_crit_0", "transition_history"])
                + "\n"
                + "\t".join(["1", "125.09", "300", "success", "162.5", "0-(0)->1"])
                + "\n",
                encoding="ascii",
            )
            return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

        with tempfile.TemporaryDirectory() as tmpdir:
            result = ewpt.run_trsm_ewpt(point, workdir=Path(tmpdir), runner=fake_runner)

            self.assertEqual(len(calls), 1)
            self.assertEqual(result.returncode, 0)
            self.assertTrue(result.input_path.exists())
            self.assertTrue(result.output_path.exists())
            self.assertEqual(result.rows[0]["status_crit_0"], "success")
            self.assertTrue(math.isclose(result.rows[0]["T_crit_0"], 162.5))

    def test_relative_workdir_uses_absolute_paths_for_calctemps(self):
        ewpt = load_module()
        point = ewpt.TRSMEWPTPoint(m2=300.0, m3=400.0, vs=200.0, a12=0.2, lx=0.1, lphix=0.05, lsx=0.05)

        def fake_runner(command, check, capture_output, text, cwd):
            self.assertTrue(Path(command[2]).is_absolute())
            self.assertTrue(Path(command[3]).is_absolute())
            self.assertTrue(Path(cwd).is_absolute())
            output_path = Path(command[3])
            output_path.write_text(
                "\t".join(["", "m1", "m2", "status_crit_0"])
                + "\n"
                + "\t".join(["1", "125.09", "300", "success"])
                + "\n",
                encoding="ascii",
            )
            return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

        with tempfile.TemporaryDirectory() as tmpdir:
            original_cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                result = ewpt.run_trsm_ewpt(
                    point,
                    workdir=Path("ewpt_results"),
                    runner=fake_runner,
                )
            finally:
                os.chdir(original_cwd)

        self.assertTrue(result.output_path.is_absolute())
        self.assertEqual(result.rows[0]["status_crit_0"], "success")

    def test_keep_files_leaves_generated_files_after_result_is_deleted(self):
        ewpt = load_module()
        point = ewpt.TRSMEWPTPoint(m2=300.0, m3=400.0, vs=200.0, a12=0.2, lx=0.1, lphix=0.05, lsx=0.05)

        def fake_runner(command, check, capture_output, text, cwd):
            output_path = Path(command[3])
            output_path.write_text(
                "\t".join(["", "m1", "m2", "status_crit_0"])
                + "\n"
                + "\t".join(["1", "125.09", "300", "success"])
                + "\n",
                encoding="ascii",
            )
            return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

        result = ewpt.run_trsm_ewpt(point, keep_files=True, runner=fake_runner)
        run_dir = result.output_path.parent
        output_path = result.output_path

        try:
            self.assertTrue(output_path.exists())
            del result
            gc.collect()
            self.assertTrue(output_path.exists())
        finally:
            shutil.rmtree(run_dir, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
