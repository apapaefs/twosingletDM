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
from unittest import mock


SCRIPT_PATH = Path(__file__).resolve().parent / "test_trsm_ewpt.py"
BSMPT_FIXTURE = Path(__file__).resolve().parents[1] / "BSMPT" / "test.output_new.csv"


def load_module():
    spec = importlib.util.spec_from_file_location("test_trsm_ewpt", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_minimatracer_fixture(path):
    header = [
        "",
        "m1",
        "status_nlo_stability",
        "status_ewsr",
        "status_tracing",
        "Temp_0",
        "w1(Temp_0)",
        "wx(Temp_0)",
        "ws(Temp_0)",
        "Veff(Temp_0)",
        "Temp_1",
        "w1(Temp_1)",
        "wx(Temp_1)",
        "ws(Temp_1)",
        "Veff(Temp_1)",
        "runtime",
    ]
    rows = [
        ["1", "125.09", "success", "ew_sym_res", "success", "10", "0", "0", "5", "-10", "nan", "nan", "nan", "nan", "nan", "0.1"],
        ["1", "125.09", "success", "ew_sym_res", "success", "0", "20", "0", "6", "-1", "5", "8", "0", "5", "-7", "0.1"],
        ["1", "125.09", "success", "ew_sym_res", "success", "10", "0", "0", "4", "-9", "6", "9", "0", "5", "-8", "0.1"],
    ]
    with Path(path).open("w", encoding="ascii", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def command_option(command, name):
    prefix = f"--{name}="
    for argument in command:
        if argument.startswith(prefix):
            return argument[len(prefix) :]
    raise AssertionError(f"missing command option {prefix}")


def synthetic_phase_traces(ewpt):
    return [
        ewpt.PhaseTrace(
            index=0,
            samples=[
                ewpt.PhaseSample(temp=80.0, w1=0.0, wx=0.0, ws=0.0, veff=-100.0),
                ewpt.PhaseSample(temp=100.0, w1=0.0, wx=0.0, ws=0.0, veff=-120.0),
            ],
        ),
        ewpt.PhaseTrace(
            index=1,
            samples=[
                ewpt.PhaseSample(temp=40.0, w1=0.0, wx=0.0, ws=50.0, veff=-400.0),
                ewpt.PhaseSample(temp=80.0, w1=0.0, wx=0.0, ws=40.0, veff=-200.0),
            ],
        ),
        ewpt.PhaseTrace(
            index=2,
            samples=[
                ewpt.PhaseSample(temp=0.0, w1=246.0, wx=0.0, ws=100.0, veff=-800.0),
                ewpt.PhaseSample(temp=50.0, w1=180.0, wx=0.0, ws=80.0, veff=-500.0),
            ],
        ),
    ]


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
            use_multithreading=True,
        )

        command = ewpt.build_calctemps_command(
            config,
            input_path=Path("/tmp/input.tsv"),
            output_path=Path("/tmp/output.csv"),
        )

        self.assertEqual(command[0], "/tmp/CalcTemps")
        self.assertIn("--model=trsm", command)
        self.assertIn("--input=/tmp/input.tsv", command)
        self.assertIn("--output=/tmp/output.csv", command)
        self.assertIn("--firstline=2", command)
        self.assertIn("--lastline=2", command)
        self.assertIn("--multistepmode=default", command)
        self.assertIn("--usemultithreading=true", command)
        self.assertIn("--thigh=300.0", command)

    def test_builds_minimatracer_command_from_config(self):
        ewpt = load_module()
        config = ewpt.EWPTConfig(
            executable=Path("/tmp/bin/CalcTemps"),
            minima_executable=Path("/tmp/bin/MinimaTracer"),
            multistepmode="2",
            use_multithreading=True,
            thigh=500.0,
        )

        command = ewpt.build_minimatracer_command(
            config,
            input_path=Path("/tmp/input.tsv"),
            output_prefix=Path("/tmp/minima"),
        )

        self.assertEqual(command[0], "/tmp/bin/MinimaTracer")
        self.assertIn("--model=trsm", command)
        self.assertIn("--input=/tmp/input.tsv", command)
        self.assertIn("--output=/tmp/minima", command)
        self.assertIn("--firstline=2", command)
        self.assertIn("--lastline=2", command)
        self.assertIn("--multistepmode=2", command)
        self.assertIn("--usemultithreading=true", command)
        self.assertIn("--thigh=500.0", command)

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

    def test_calculates_fopt_strengths_from_calctemps_vev_jumps(self):
        ewpt = load_module()
        row = ewpt.parse_calctemps_output(BSMPT_FIXTURE)[0]

        strengths = ewpt.calculate_fopt_strengths(row)
        by_key = {
            (strength.transition_index, strength.temperature_kind): strength
            for strength in strengths
        }
        nucl = by_key[(0, "nucl")]
        expected_ew_jump = abs(row["w1_nucl_true_0"] - row["w1_nucl_false_0"]) / row["T_nucl_0"]
        expected_field_jump = math.sqrt(
            (row["w1_nucl_true_0"] - row["w1_nucl_false_0"]) ** 2
            + (row["wx_nucl_true_0"] - row["wx_nucl_false_0"]) ** 2
            + (row["ws_nucl_true_0"] - row["ws_nucl_false_0"]) ** 2
        ) / row["T_nucl_0"]

        self.assertTrue(math.isclose(nucl.ew_jump_over_T, expected_ew_jump))
        self.assertTrue(math.isclose(nucl.field_jump_over_T, expected_field_jump))
        self.assertTrue(math.isclose(nucl.ew_true_over_T, abs(row["w1_nucl_true_0"]) / row["T_nucl_0"]))
        self.assertEqual(nucl.false_vev["ws"], row["ws_nucl_false_0"])
        self.assertEqual(nucl.true_vev["w1"], row["w1_nucl_true_0"])

    def test_parses_minimatracer_output_with_nan_and_duplicate_temperature_handling(self):
        ewpt = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "minima_1.tsv"
            write_minimatracer_fixture(path)

            traces = ewpt.parse_minimatracer_output(path)

        self.assertEqual([trace.index for trace in traces], [0, 1])
        self.assertEqual([sample.temp for sample in traces[0].samples], [0.0, 10.0])
        self.assertTrue(math.isclose(traces[0].samples[1].veff, -10.0))
        self.assertEqual([sample.temp for sample in traces[1].samples], [5.0, 6.0])

    def test_reconstructs_and_classifies_global_branch_on_observed_union_grid(self):
        ewpt = load_module()

        analysis = ewpt.analyze_minimatracer_phases(
            synthetic_phase_traces(ewpt),
            ewpt.PhaseClassificationThresholds(),
        )

        self.assertEqual(analysis.temperature_grid, [0.0, 40.0, 50.0, 80.0, 100.0])
        self.assertEqual(analysis.global_phase_path, ["SYM", "SINGLET_S", "EW"])
        self.assertEqual(analysis.ew_step_index, 1)
        cooling_labels = [point.label for point in sorted(analysis.global_branch, key=lambda point: point.temp, reverse=True)]
        self.assertEqual(cooling_labels[0], "SYM")
        self.assertEqual(cooling_labels[-1], "EW")

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
        self.assertIn("fopt_strength_nucl_0:", summary)
        self.assertIn("ew_jump/T=", summary)
        self.assertIn("field_jump/T=", summary)

    def test_summary_prints_all_calctemps_transition_indices(self):
        ewpt = load_module()
        row = {
            "status_tracing": "success",
            "status_coex_pairs": "success",
            "status_crit_0": "success",
            "T_crit_0": 523.7,
            "status_nucl_0": "success",
            "T_nucl_0": 523.5,
            "status_crit_1": "success",
            "T_crit_1": 159.5,
            "status_nucl_1": "not_met",
            "T_nucl_1": math.nan,
            "status_perc_1": "success",
            "T_perc_1": 159.4,
            "transition_history": "0-(0)->1-(1)->2",
        }
        result = ewpt.EWPTResult(
            rows=[row],
            input_path=Path("/tmp/TRSM_Input.tsv"),
            output_path=Path("/tmp/test.output_new.csv"),
            command=[],
            returncode=0,
        )

        summary = ewpt.summarize_result(result)

        self.assertIn("status_crit_1: success", summary)
        self.assertIn("T_crit_1: 159.5", summary)
        self.assertIn("status_nucl_1: not_met", summary)
        self.assertIn("T_nucl_1: nan", summary)
        self.assertIn("status_perc_1: success", summary)

    def test_json_includes_transition_strengths(self):
        ewpt = load_module()
        row = ewpt.parse_calctemps_output(BSMPT_FIXTURE)[0]
        result = ewpt.EWPTResult(
            rows=[row],
            input_path=Path("/tmp/TRSM_Input.tsv"),
            output_path=Path("/tmp/test.output_new.csv"),
            command=[],
            returncode=0,
        )

        payload = ewpt.result_to_json(result)

        self.assertIn("transition_strengths", payload)
        self.assertGreater(len(payload["transition_strengths"]), 0)
        self.assertEqual(payload["transition_strengths"][0]["transition_index"], 0)

    def test_summary_prints_minimatracer_diagnosis(self):
        ewpt = load_module()
        row = ewpt.parse_calctemps_output(BSMPT_FIXTURE)[0]
        analysis = ewpt.analyze_minimatracer_phases(
            synthetic_phase_traces(ewpt),
            ewpt.PhaseClassificationThresholds(),
        )
        result = ewpt.EWPTResult(
            rows=[row],
            input_path=Path("/tmp/TRSM_Input.tsv"),
            output_path=Path("/tmp/test.output_new.csv"),
            command=[],
            returncode=0,
            minima_output_path=Path("/tmp/minima_1.tsv"),
            minima_analysis=analysis,
            plot_paths=(Path("/tmp/minima_phases.png"), Path("/tmp/minima_phases.pdf")),
        )

        summary = ewpt.summarize_result(result)

        self.assertIn("minima_output: /tmp/minima_1.tsv", summary)
        self.assertIn("phase_plot: /tmp/minima_phases.png", summary)
        self.assertIn("phase_plot: /tmp/minima_phases.pdf", summary)
        self.assertIn("global_phase_path: SYM -> SINGLET_S -> EW", summary)
        self.assertIn("ew_step_index: 1", summary)

    def test_run_trsm_ewpt_uses_fake_runner_and_parses_output(self):
        ewpt = load_module()
        point = ewpt.TRSMEWPTPoint(m2=300.0, m3=400.0, vs=200.0, a12=0.2, lx=0.1, lphix=0.05, lsx=0.05)
        calls = []

        def fake_runner(command, check, capture_output, text, cwd):
            calls.append((command, check, capture_output, text, cwd))
            if Path(command[0]).name == "MinimaTracer":
                write_minimatracer_fixture(str(Path(command_option(command, "output"))) + "_1.tsv")
            else:
                output_path = Path(command_option(command, "output"))
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

            self.assertEqual([Path(call[0][0]).name for call in calls], ["MinimaTracer", "CalcTemps"])
            self.assertEqual(result.returncode, 0)
            self.assertTrue(result.input_path.exists())
            self.assertTrue(result.output_path.exists())
            self.assertTrue(result.minima_output_path.exists())
            self.assertEqual(result.rows[0]["status_crit_0"], "success")
            self.assertTrue(math.isclose(result.rows[0]["T_crit_0"], 162.5))
            self.assertEqual(result.minima_analysis.global_phase_path, ["SINGLET_S", "EW"])

    def test_relative_workdir_uses_absolute_paths_for_calctemps(self):
        ewpt = load_module()
        point = ewpt.TRSMEWPTPoint(m2=300.0, m3=400.0, vs=200.0, a12=0.2, lx=0.1, lphix=0.05, lsx=0.05)

        def fake_runner(command, check, capture_output, text, cwd):
            self.assertTrue(Path(command_option(command, "input")).is_absolute())
            self.assertTrue(Path(command_option(command, "output")).is_absolute())
            self.assertTrue(Path(cwd).is_absolute())
            if Path(command[0]).name == "MinimaTracer":
                write_minimatracer_fixture(str(Path(command_option(command, "output"))) + "_1.tsv")
            else:
                output_path = Path(command_option(command, "output"))
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
            if Path(command[0]).name == "MinimaTracer":
                write_minimatracer_fixture(str(Path(command_option(command, "output"))) + "_1.tsv")
            else:
                output_path = Path(command_option(command, "output"))
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

    @unittest.skipIf(importlib.util.find_spec("matplotlib") is None, "matplotlib is not installed")
    def test_plot_minimatracer_phases_writes_png_and_pdf(self):
        ewpt = load_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            output_base = Path(tmpdir) / "phase_plot"
            paths = ewpt.plot_minimatracer_phases(
                synthetic_phase_traces(ewpt),
                output_base,
                plot_format="both",
            )

            self.assertEqual([path.suffix for path in paths], [".png", ".pdf"])
            for path in paths:
                self.assertTrue(path.exists())
                self.assertGreater(path.stat().st_size, 0)

    @unittest.skipIf(importlib.util.find_spec("matplotlib") is None, "matplotlib is not installed")
    def test_plot_minimatracer_phases_uses_requested_thigh_as_x_axis_limit(self):
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt

        ewpt = load_module()
        recorded_limits = []
        original_subplots = plt.subplots

        def subplots_spy(*args, **kwargs):
            fig, ax = original_subplots(*args, **kwargs)
            original_set_xlim = ax.set_xlim

            def set_xlim_spy(*limits, **set_xlim_kwargs):
                recorded_limits.append(limits)
                return original_set_xlim(*limits, **set_xlim_kwargs)

            ax.set_xlim = set_xlim_spy
            return fig, ax

        with tempfile.TemporaryDirectory() as tmpdir:
            with mock.patch("matplotlib.pyplot.subplots", side_effect=subplots_spy):
                ewpt.plot_minimatracer_phases(
                    synthetic_phase_traces(ewpt),
                    Path(tmpdir) / "phase_plot",
                    plot_format="png",
                    thigh=500.0,
                )

        self.assertIn(((0.0, 500.0),), recorded_limits)


if __name__ == "__main__":
    unittest.main()
