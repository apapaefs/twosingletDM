import importlib.util
import contextlib
import io
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from types import SimpleNamespace


SCRIPT_PATH = Path(__file__).resolve().parent / "generate_trsm_points.py"


def install_stub_modules():
    stubs = {}

    generate_trsm_info = types.ModuleType("generate_trsm_info")
    generate_trsm_info.__all__ = []
    stubs["generate_trsm_info"] = generate_trsm_info

    for name in [
        "test_trsm_evolution",
        "test_trsm_higgstools",
        "trsm_kstoalphas",
        "generate_mg5_trsm_xsecs",
        "test_trsm_theory_constraints",
        "singlet_EWPO",
    ]:
        module = types.ModuleType(name)
        module.__all__ = []
        stubs[name] = module

    mg5_process_runner = types.ModuleType("mg5_process_runner")
    mg5_process_runner.run_mg5_processes = lambda *args, **kwargs: {}
    stubs["mg5_process_runner"] = mg5_process_runner

    scan_output = types.ModuleType("scan_output")
    scan_output.write_valid_point = lambda *args, **kwargs: None
    stubs["scan_output"] = scan_output

    test_trsm_dm = types.ModuleType("test_trsm_DM")
    test_trsm_dm.print_dm_info = lambda *args, **kwargs: None
    test_trsm_dm.test_dm = lambda *args, **kwargs: (False, {}, {})
    stubs["test_trsm_DM"] = test_trsm_dm

    prettytable = types.ModuleType("prettytable")

    class PrettyTable:
        def __init__(self, *args, **kwargs):
            self.rows = []

        def add_row(self, row):
            self.rows.append(row)

        def __str__(self):
            return "\n".join(str(row) for row in self.rows)

    prettytable.PrettyTable = PrettyTable
    stubs["prettytable"] = prettytable

    originals = {}
    for name, module in stubs.items():
        originals[name] = sys.modules.get(name)
        sys.modules[name] = module
    return originals


def restore_stub_modules(originals):
    for name, module in originals.items():
        if module is None:
            sys.modules.pop(name, None)
        else:
            sys.modules[name] = module


def load_generator_module(argv=None):
    originals = install_stub_modules()
    old_argv = sys.argv[:]
    sys.argv = [str(SCRIPT_PATH)] + list(argv or ["--nrandom", "0"])
    module_name = f"generate_trsm_points_under_test_{id(argv)}"
    try:
        spec = importlib.util.spec_from_file_location(module_name, SCRIPT_PATH)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.argv = old_argv
        restore_stub_modules(originals)


def ewpt_args(**overrides):
    values = {
        "run_ewpt": True,
        "ewpt_executable": Path("/tmp/CalcTemps"),
        "ewpt_minima_executable": Path("/tmp/MinimaTracer"),
        "ewpt_thigh": 800.0,
        "ewpt_multistepmode": "default",
        "ewpt_workdir": None,
        "ewpt_plot_phases": False,
        "ewpt_plot_output": None,
        "ewpt_plot_format": "both",
        "ewpt_require_eq418": False,
        "ewpt_sym_threshold": 1.0,
        "ewpt_w1_threshold": 5.0,
        "ewpt_wx_threshold": 1.0,
        "ewpt_ws_threshold": 1.0,
        "run_ewpt_on_dm_failed": False,
        "write_dm_failed": False,
        "write_all_points": False,
        "higgstools_details": False,
        "higgstools_top": 5,
        "save_higgstools_details": True,
        "resonantDM1": False,
        "resonantDM2": False,
        "m1": 125.09,
        "m3": None,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def viable_point_info():
    return {
        "M2": 380.0,
        "M3": 500.0,
        "vs": 200.0,
        "a12": -0.15,
        "lX": 0.10,
        "lPhiX": 0.044,
        "lSX": 0.15,
    }


class FakeEWPTModule:
    def __init__(self, eq418_satisfied=True):
        self.calls = []
        self.eq418_satisfied = eq418_satisfied
        self.eq418_checks = []

    class TRSMEWPTPoint:
        def __init__(self, **kwargs):
            self.__dict__.update(kwargs)

    class EWPTConfig:
        def __init__(self, **kwargs):
            self.__dict__.update(kwargs)

    def run_trsm_ewpt(self, point, config, workdir, keep_files):
        self.calls.append(
            {
                "point": point,
                "config": config,
                "workdir": Path(workdir),
                "keep_files": keep_files,
            }
        )
        return SimpleNamespace(tag="fake-result")

    def summarize_result(self, result):
        return "output: fake.csv\nglobal_phase_path: SINGLET_S -> EW"

    def result_to_json(self, result):
        return {
            "tag": result.tag,
            "minimatracer": {
                "global_phase_path": ["SINGLET_S", "X_BROKEN", "EW_X_BROKEN", "EW"],
                "ew_step_index": 2,
            },
            "transition_strengths": [
                {
                    "temperature_kind": "crit",
                    "transition_index": 0,
                    "ew_true_over_T": 0.4,
                },
                {
                    "temperature_kind": "nucl",
                    "transition_index": 0,
                    "ew_true_over_T": 0.9,
                },
            ],
        }

    def check_eq_4_18(self, row):
        self.eq418_checks.append(row)
        return SimpleNamespace(
            satisfied=self.eq418_satisfied,
            conditions={"D > 0": self.eq418_satisfied},
        )


class FailingEWPTModule(FakeEWPTModule):
    def run_trsm_ewpt(self, point, config, workdir, keep_files):
        self.calls.append(
            {
                "point": point,
                "config": config,
                "workdir": Path(workdir),
                "keep_files": keep_files,
            }
        )
        raise RuntimeError("BSMPT exploded")


class TestGenerateTRSMPointsEWPT(unittest.TestCase):
    def test_scan_output_includes_ewpt_columns(self):
        import scan_output

        self.assertIn("higgstools_hb_selected_limits", scan_output.output_columns({}))
        self.assertIn("higgstools_hb_top_obs", scan_output.output_columns({}))
        self.assertIn("higgstools_hs_chi2", scan_output.output_columns({}))
        self.assertIn("higgstools_hs_delta_chi2", scan_output.output_columns({}))
        self.assertIn("higgstools_hs_top_chi2", scan_output.output_columns({}))
        self.assertIn("ewpt_status", scan_output.output_columns({}))
        self.assertIn("ewpt_error", scan_output.output_columns({}))
        self.assertIn("dm_indirect_available", scan_output.output_columns({}))
        self.assertIn("dm_indirect_ratio", scan_output.output_columns({}))
        self.assertIn("dm_indirect_detection_excluded", scan_output.output_columns({}))
        self.assertIn("ewpt_ew_true_over_T", scan_output.output_columns({}))
        self.assertIn("ewpt_global_phase_path", scan_output.output_columns({}))
        self.assertIn("ewpt_has_x_broken", scan_output.output_columns({}))
        self.assertIn("ewpt_ew_step_index", scan_output.output_columns({}))

    def test_parse_args_defaults_to_not_running_ewpt(self):
        generator = load_generator_module()

        args = generator.parse_args([])

        self.assertFalse(args.run_ewpt)
        self.assertFalse(args.write_all_points)
        self.assertTrue(args.save_higgstools_details)

    def test_parse_args_accepts_ewpt_options(self):
        generator = load_generator_module()

        args = generator.parse_args(
            [
                "--run-ewpt",
                "--ewpt-thigh",
                "1000",
                "--ewpt-workdir",
                "ewpt-out",
                "--ewpt-plot-phases",
                "--ewpt-plot-format",
                "png",
                "--ewpt-require-eq418",
            ]
        )

        self.assertTrue(args.run_ewpt)
        self.assertEqual(args.ewpt_thigh, 1000.0)
        self.assertEqual(args.ewpt_workdir, Path("ewpt-out"))
        self.assertTrue(args.ewpt_plot_phases)
        self.assertEqual(args.ewpt_plot_format, "png")
        self.assertTrue(args.ewpt_require_eq418)

    def test_parse_args_accepts_higgstools_detail_options(self):
        generator = load_generator_module()

        args = generator.parse_args(
            ["--higgstools-details", "--higgstools-top", "7", "--no-save-higgstools-details"]
        )

        self.assertTrue(args.higgstools_details)
        self.assertEqual(args.higgstools_top, 7)
        self.assertFalse(args.save_higgstools_details)

    def test_parse_args_accepts_resonant_dm_options(self):
        generator = load_generator_module()

        args = generator.parse_args(["--resonantDM1"])

        self.assertTrue(args.resonantDM1)
        self.assertFalse(args.resonantDM2)

    def test_parse_args_rejects_both_resonant_dm_options(self):
        generator = load_generator_module()

        with self.assertRaises(SystemExit):
            generator.parse_args(["--resonantDM1", "--resonantDM2"])

    def test_explicit_point_can_omit_m3_when_resonant_dm_is_enabled(self):
        generator = load_generator_module()

        args = generator.parse_args(
            [
                "--m2",
                "7.526",
                "--vs",
                "768.5",
                "--a12",
                "0.3078",
                "--lx",
                "0.4426",
                "--lphix",
                "8.696e-05",
                "--lsx",
                "0.0002258",
                "--resonantDM1",
            ]
        )

        self.assertTrue(generator.has_explicit_point(args))
        self.assertAlmostEqual(generator.effective_m3(args, 7.526), 250.18)

    def test_effective_m3_uses_resonant_dm_modes(self):
        generator = load_generator_module()

        base = ewpt_args(resonantDM1=False, resonantDM2=False, m1=125.09, m3=600.0)
        self.assertEqual(generator.effective_m3(base, 200.0), 600.0)

        dm1 = ewpt_args(resonantDM1=True, resonantDM2=False, m1=125.09, m3=None)
        self.assertAlmostEqual(generator.effective_m3(dm1, 200.0), 250.18)

        dm2 = ewpt_args(resonantDM1=False, resonantDM2=True, m1=125.09, m3=None)
        self.assertEqual(generator.effective_m3(dm2, 200.0), 400.0)

    def test_parse_args_accepts_write_all_points(self):
        generator = load_generator_module()

        args = generator.parse_args(["--write-all-points"])

        self.assertTrue(args.write_all_points)

    def test_parse_args_accepts_dm_failed_output_option(self):
        generator = load_generator_module()

        args = generator.parse_args(["--write-dm-failed"])

        self.assertTrue(args.write_dm_failed)

    def test_parse_args_accepts_ewpt_on_dm_failed_option(self):
        generator = load_generator_module()

        args = generator.parse_args(["--run-ewpt-on-dm-failed"])

        self.assertTrue(args.run_ewpt_on_dm_failed)

    def test_dm_failed_option_writes_sidecar_for_otherwise_allowed_point(self):
        generator = load_generator_module()
        records = []

        def fake_write(path, point_info, mg5xsecs=None):
            records.append((Path(path), dict(point_info)))

        generator.write_valid_point_file = fake_write
        generator.RunTag = "unit"
        generator.OutputDir = "output/"
        generator.cli_args = SimpleNamespace(
            run_ewpt=False,
            write_dm_failed=True,
            run_ewpt_on_dm_failed=False,
        )
        generator.np = SimpleNamespace(sin=lambda value: value)
        for name in [
            "Mz",
            "Mw",
            "Delta_S_central_wU",
            "Delta_T_central_wU",
            "Delta_U_central_wU",
            "errS_wU",
            "errT_wU",
            "errU_wU",
            "covST_wU",
            "covSU_wU",
            "covTU_wU",
            "pred",
            "H1",
            "H2",
            "H3",
        ]:
            setattr(generator, name, 0)
        generator.generate_lams = lambda *args, **kwargs: (
            200.0,
            0.0,
            380.0,
            500.0,
            -0.15,
            0.0,
            0.0,
            1.0,
            2.0,
            3.0,
            11.0,
            12.0,
            13.0,
            123.0,
            122.0,
            1111.0,
            1112.0,
            1113.0,
            133.0,
            0.99,
            -0.1,
            0.0,
            {},
            {},
            {},
            0.0,
            0.0,
            0.0,
        )
        generator.check_EWPO_wU = lambda *args, **kwargs: True
        generator.check_wmass_tania = lambda *args, **kwargs: True
        generator.analyze_parampoint = lambda *args, **kwargs: (True, True)
        generator.theory_constraints_vxzero = lambda *args, **kwargs: True
        generator.test_evo_vxzero = lambda *args, **kwargs: True
        generator.test_dm = lambda *args, **kwargs: (
            False,
            {},
            {"dm_mdm": 750.0, "dm_relic_excluded": True},
        )

        result = generator.evaluate_trsm_point_vxzero(
            123,
            380.0,
            500.0,
            200.0,
            -0.15,
            0.10,
            0.05,
            0.15,
            point_index=9,
        )

        self.assertEqual(result, 0)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0][0], Path("output/trsm_points_unit_dm_failed.dat"))
        self.assertEqual(records[0][1]["dm"], False)
        self.assertEqual(records[0][1]["dm_mdm"], 750.0)

    def test_ewpt_hook_skips_when_flag_disabled_or_point_not_viable(self):
        generator = load_generator_module()
        fake_ewpt = FakeEWPTModule()

        disabled = ewpt_args(run_ewpt=False)
        self.assertIsNone(
            generator.run_ewpt_if_requested(
                point_info=viable_point_info(),
                passed=True,
                args=disabled,
                point_index=1,
                ewpt_module=fake_ewpt,
            )
        )

        enabled = ewpt_args(run_ewpt=True)
        self.assertIsNone(
            generator.run_ewpt_if_requested(
                point_info=viable_point_info(),
                passed=False,
                args=enabled,
                point_index=1,
                ewpt_module=fake_ewpt,
            )
        )
        self.assertEqual(fake_ewpt.calls, [])

    def test_ewpt_hook_runs_for_dm_failed_point_when_explicitly_enabled(self):
        generator = load_generator_module()
        fake_ewpt = FakeEWPTModule()
        point_info = viable_point_info()

        with tempfile.TemporaryDirectory() as tmpdir:
            result = generator.run_ewpt_if_requested(
                point_info=point_info,
                passed=False,
                args=ewpt_args(
                    run_ewpt=False,
                    run_ewpt_on_dm_failed=True,
                    ewpt_workdir=Path(tmpdir),
                ),
                point_index=7,
                ewpt_module=fake_ewpt,
                allow_dm_failed=True,
            )

        self.assertEqual(result.tag, "fake-result")
        self.assertEqual(len(fake_ewpt.calls), 1)
        self.assertEqual(fake_ewpt.calls[0]["point"].index, 7)
        self.assertEqual(point_info["ewpt_ew_true_over_T"], 0.9)
        self.assertEqual(
            point_info["ewpt_global_phase_path"],
            "SINGLET_S -> X_BROKEN -> EW_X_BROKEN -> EW",
        )
        self.assertTrue(point_info["ewpt_has_x_broken"])
        self.assertEqual(point_info["ewpt_ew_step_index"], 2)

    def test_ewpt_hook_skips_when_eq418_prefilter_fails(self):
        generator = load_generator_module()
        fake_ewpt = FakeEWPTModule(eq418_satisfied=False)

        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            result = generator.run_ewpt_if_requested(
                point_info=viable_point_info(),
                passed=True,
                args=ewpt_args(ewpt_require_eq418=True),
                point_index=4,
                ewpt_module=fake_ewpt,
            )

        self.assertIsNone(result)
        self.assertEqual(fake_ewpt.calls, [])
        self.assertEqual(fake_ewpt.eq418_checks[0]["m2"], 380.0)
        self.assertEqual(fake_ewpt.eq418_checks[0]["lphix"], 0.044)
        output = stdout.getvalue()
        self.assertIn("Viable EWPT candidate:", output)
        self.assertIn("m2=380.0", output)
        self.assertIn("m3=500.0", output)
        self.assertIn("lphix=0.044", output)
        self.assertIn("failed Eq. 4.18", output)

    def test_ewpt_hook_runs_for_viable_point_and_writes_summary_files(self):
        generator = load_generator_module()
        fake_ewpt = FakeEWPTModule()

        with tempfile.TemporaryDirectory() as tmpdir:
            result = generator.run_ewpt_if_requested(
                point_info=viable_point_info(),
                passed=True,
                args=ewpt_args(ewpt_workdir=Path(tmpdir)),
                point_index=3,
                ewpt_module=fake_ewpt,
            )

            workdir = Path(tmpdir) / "point_000003"
            summary_path = workdir / "ewpt_summary.txt"
            json_path = workdir / "ewpt_result.json"

            self.assertEqual(result.tag, "fake-result")
            self.assertEqual(len(fake_ewpt.calls), 1)
            self.assertEqual(fake_ewpt.eq418_checks, [])
            self.assertEqual(fake_ewpt.calls[0]["point"].m2, 380.0)
            self.assertEqual(fake_ewpt.calls[0]["point"].m3, 500.0)
            self.assertEqual(fake_ewpt.calls[0]["point"].lx, 0.10)
            self.assertEqual(fake_ewpt.calls[0]["config"].thigh, 800.0)
            self.assertTrue(fake_ewpt.calls[0]["keep_files"])
            self.assertTrue(summary_path.exists())
            self.assertIn("global_phase_path", summary_path.read_text(encoding="utf-8"))
            self.assertEqual(
                json.loads(json_path.read_text(encoding="utf-8"))["tag"],
                "fake-result",
            )

    def test_ewpt_hook_records_primary_ew_true_over_T_on_point_info(self):
        generator = load_generator_module()
        fake_ewpt = FakeEWPTModule()
        point_info = viable_point_info()

        with tempfile.TemporaryDirectory() as tmpdir:
            generator.run_ewpt_if_requested(
                point_info=point_info,
                passed=True,
                args=ewpt_args(ewpt_workdir=Path(tmpdir)),
                point_index=3,
                ewpt_module=fake_ewpt,
            )

        self.assertEqual(point_info["ewpt_ew_true_over_T"], 0.9)
        self.assertEqual(point_info["ewpt_status"], "success")
        self.assertEqual(
            point_info["ewpt_global_phase_path"],
            "SINGLET_S -> X_BROKEN -> EW_X_BROKEN -> EW",
        )
        self.assertTrue(point_info["ewpt_has_x_broken"])
        self.assertEqual(point_info["ewpt_ew_step_index"], 2)

    def test_ewpt_hook_records_failure_instead_of_aborting_scan(self):
        generator = load_generator_module()
        fake_ewpt = FailingEWPTModule()
        point_info = viable_point_info()

        with tempfile.TemporaryDirectory() as tmpdir:
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                result = generator.run_ewpt_if_requested(
                    point_info=point_info,
                    passed=True,
                    args=ewpt_args(ewpt_workdir=Path(tmpdir)),
                    point_index=5,
                    ewpt_module=fake_ewpt,
                )

            workdir = Path(tmpdir) / "point_000005"
            error_path = workdir / "ewpt_error.txt"

            self.assertIsNone(result)
            self.assertEqual(len(fake_ewpt.calls), 1)
            self.assertEqual(point_info["ewpt_status"], "failed")
            self.assertIn("BSMPT exploded", point_info["ewpt_error"])
            self.assertTrue(error_path.exists())
            self.assertIn("BSMPT exploded", error_path.read_text(encoding="utf-8"))
            self.assertIn("EWPT analysis failed", stdout.getvalue())

    def test_viable_point_writes_record_after_ewpt_strength_hook(self):
        generator = load_generator_module()
        rows = []

        def fake_write(path, point_info, mg5xsecs=None):
            rows.append(dict(point_info))

        generator.write_valid_point_file = fake_write
        generator.RunTag = "unit"
        generator.OutputDir = "output/"
        generator.cli_args = ewpt_args()
        generator.cli_args.write_dm_failed = False
        generator.np = SimpleNamespace(sin=lambda value: value)
        for name in [
            "Mz",
            "Mw",
            "Delta_S_central_wU",
            "Delta_T_central_wU",
            "Delta_U_central_wU",
            "errS_wU",
            "errT_wU",
            "errU_wU",
            "covST_wU",
            "covSU_wU",
            "covTU_wU",
            "pred",
            "H1",
            "H2",
            "H3",
        ]:
            setattr(generator, name, 0)
        generator.generate_lams = lambda *args, **kwargs: (
            200.0,
            0.0,
            380.0,
            500.0,
            -0.15,
            0.0,
            0.0,
            1.0,
            2.0,
            3.0,
            11.0,
            12.0,
            13.0,
            123.0,
            122.0,
            1111.0,
            1112.0,
            1113.0,
            133.0,
            0.99,
            -0.1,
            0.0,
            {},
            {},
            {},
            0.0,
            0.0,
            0.0,
        )
        generator.check_EWPO_wU = lambda *args, **kwargs: True
        generator.check_wmass_tania = lambda *args, **kwargs: True
        generator.analyze_parampoint = lambda *args, **kwargs: (True, True)
        generator.theory_constraints_vxzero = lambda *args, **kwargs: True
        generator.test_evo_vxzero = lambda *args, **kwargs: True
        generator.test_dm = lambda *args, **kwargs: (
            True,
            {},
            {"dm_mdm": 750.0, "dm_relic_excluded": False},
        )

        def fake_ewpt(point_info, passed, args, point_index):
            point_info["ewpt_ew_true_over_T"] = 1.23
            point_info["ewpt_global_phase_path"] = "SINGLET_S -> EW"
            point_info["ewpt_has_x_broken"] = False
            point_info["ewpt_ew_step_index"] = 1

        generator.run_ewpt_if_requested = fake_ewpt

        result = generator.evaluate_trsm_point_vxzero(
            123,
            380.0,
            500.0,
            200.0,
            -0.15,
            0.10,
            0.05,
            0.15,
            point_index=9,
        )

        self.assertEqual(result, 1)
        self.assertEqual(rows[0]["ewpt_ew_true_over_T"], 1.23)
        self.assertEqual(rows[0]["ewpt_global_phase_path"], "SINGLET_S -> EW")
        self.assertFalse(rows[0]["ewpt_has_x_broken"])
        self.assertEqual(rows[0]["ewpt_ew_step_index"], 1)

    def test_higgstools_details_are_forwarded_to_analyze_parampoint(self):
        generator = load_generator_module()
        calls = []

        generator.write_valid_point_file = lambda *args, **kwargs: None
        generator.RunTag = "unit"
        generator.OutputDir = "output/"
        generator.cli_args = ewpt_args(
            run_ewpt=False,
            higgstools_details=True,
            higgstools_top=7,
        )
        generator.np = SimpleNamespace(sin=lambda value: value)
        for name in [
            "Mz",
            "Mw",
            "Delta_S_central_wU",
            "Delta_T_central_wU",
            "Delta_U_central_wU",
            "errS_wU",
            "errT_wU",
            "errU_wU",
            "covST_wU",
            "covSU_wU",
            "covTU_wU",
            "pred",
            "H1",
            "H2",
            "H3",
        ]:
            setattr(generator, name, 0)
        generator.generate_lams = lambda *args, **kwargs: (
            200.0,
            0.0,
            380.0,
            500.0,
            -0.15,
            0.0,
            0.0,
            1.0,
            2.0,
            3.0,
            11.0,
            12.0,
            13.0,
            123.0,
            122.0,
            1111.0,
            1112.0,
            1113.0,
            133.0,
            0.99,
            -0.1,
            0.0,
            {},
            {},
            {},
            0.0,
            0.0,
            0.0,
        )
        generator.check_EWPO_wU = lambda *args, **kwargs: True
        generator.check_wmass_tania = lambda *args, **kwargs: True
        generator.theory_constraints_vxzero = lambda *args, **kwargs: True
        generator.test_evo_vxzero = lambda *args, **kwargs: True
        generator.test_dm = lambda *args, **kwargs: (
            True,
            {},
            {"dm_mdm": 750.0, "dm_relic_excluded": False},
        )

        def fake_analyze(*args, **kwargs):
            calls.append(kwargs)
            return True, True

        generator.analyze_parampoint = fake_analyze

        result = generator.evaluate_trsm_point_vxzero(
            123,
            380.0,
            500.0,
            200.0,
            -0.15,
            0.10,
            0.05,
            0.15,
            point_index=9,
        )

        self.assertEqual(result, 1)
        self.assertEqual(
            calls,
            [{"print_details": True, "details_top": 7, "return_details": True}],
        )

    def test_higgstools_details_are_saved_to_point_record_by_default(self):
        generator = load_generator_module()
        rows = []

        def fake_write(path, point_info, mg5xsecs=None):
            rows.append(dict(point_info))

        generator.write_valid_point_file = fake_write
        generator.RunTag = "unit"
        generator.OutputDir = "output/"
        generator.cli_args = ewpt_args(run_ewpt=False)
        generator.np = SimpleNamespace(sin=lambda value: value)
        for name in [
            "Mz",
            "Mw",
            "Delta_S_central_wU",
            "Delta_T_central_wU",
            "Delta_U_central_wU",
            "errS_wU",
            "errT_wU",
            "errU_wU",
            "covST_wU",
            "covSU_wU",
            "covTU_wU",
            "pred",
            "H1",
            "H2",
            "H3",
        ]:
            setattr(generator, name, 0)
        generator.generate_lams = lambda *args, **kwargs: (
            200.0,
            0.0,
            380.0,
            500.0,
            -0.15,
            0.0,
            0.0,
            1.0,
            2.0,
            3.0,
            11.0,
            12.0,
            13.0,
            123.0,
            122.0,
            1111.0,
            1112.0,
            1113.0,
            133.0,
            0.99,
            -0.1,
            0.0,
            {},
            {},
            {},
            0.0,
            0.0,
            0.0,
        )
        generator.check_EWPO_wU = lambda *args, **kwargs: True
        generator.check_wmass_tania = lambda *args, **kwargs: True
        generator.theory_constraints_vxzero = lambda *args, **kwargs: True
        generator.test_evo_vxzero = lambda *args, **kwargs: True
        generator.test_dm = lambda *args, **kwargs: (
            True,
            {},
            {"dm_mdm": 750.0, "dm_relic_excluded": False},
        )

        details = {
            "higgsbounds": {
                "selected_limits": {"H2": {"obsRatio": 1.1, "reference": "OPAL"}},
                "top_observed_ratios": [{"obsRatio": 1.1, "reference": "OPAL"}],
            },
            "higgssignals": {
                "chi2": 158.0,
                "delta_chi2": 5.5,
                "top_chi2_contributors": [{"chisq": 10.0, "reference": "CMS"}],
            },
        }

        def fake_analyze(*args, **kwargs):
            self.assertTrue(kwargs["return_details"])
            return True, True, details

        generator.analyze_parampoint = fake_analyze

        result = generator.evaluate_trsm_point_vxzero(
            123,
            380.0,
            500.0,
            200.0,
            -0.15,
            0.10,
            0.05,
            0.15,
            point_index=9,
        )

        self.assertEqual(result, 1)
        self.assertEqual(rows[0]["higgstools_hs_chi2"], 158.0)
        self.assertEqual(rows[0]["higgstools_hs_delta_chi2"], 5.5)
        self.assertIn('"H2"', rows[0]["higgstools_hb_selected_limits"])
        self.assertIn('"OPAL"', rows[0]["higgstools_hb_top_obs"])
        self.assertIn('"CMS"', rows[0]["higgstools_hs_top_chi2"])

    def test_dm_failed_point_runs_ewpt_and_writes_strength_sidecar_when_enabled(self):
        generator = load_generator_module()
        rows = []
        ewpt_calls = []

        def fake_write(path, point_info, mg5xsecs=None):
            rows.append((Path(path), dict(point_info)))

        generator.write_valid_point_file = fake_write
        generator.RunTag = "unit"
        generator.OutputDir = "output/"
        generator.cli_args = ewpt_args(run_ewpt=False, run_ewpt_on_dm_failed=True)
        generator.np = SimpleNamespace(sin=lambda value: value)
        for name in [
            "Mz",
            "Mw",
            "Delta_S_central_wU",
            "Delta_T_central_wU",
            "Delta_U_central_wU",
            "errS_wU",
            "errT_wU",
            "errU_wU",
            "covST_wU",
            "covSU_wU",
            "covTU_wU",
            "pred",
            "H1",
            "H2",
            "H3",
        ]:
            setattr(generator, name, 0)
        generator.generate_lams = lambda *args, **kwargs: (
            200.0,
            0.0,
            380.0,
            500.0,
            -0.15,
            0.0,
            0.0,
            1.0,
            2.0,
            3.0,
            11.0,
            12.0,
            13.0,
            123.0,
            122.0,
            1111.0,
            1112.0,
            1113.0,
            133.0,
            0.99,
            -0.1,
            0.0,
            {},
            {},
            {},
            0.0,
            0.0,
            0.0,
        )
        generator.check_EWPO_wU = lambda *args, **kwargs: True
        generator.check_wmass_tania = lambda *args, **kwargs: True
        generator.analyze_parampoint = lambda *args, **kwargs: (True, True)
        generator.theory_constraints_vxzero = lambda *args, **kwargs: True
        generator.test_evo_vxzero = lambda *args, **kwargs: True
        generator.test_dm = lambda *args, **kwargs: (
            False,
            {},
            {"dm_mdm": 750.0, "dm_relic_excluded": True},
        )

        def fake_ewpt(point_info, passed, args, point_index, **kwargs):
            ewpt_calls.append((passed, kwargs.get("allow_dm_failed")))
            point_info["ewpt_ew_true_over_T"] = 0.77
            point_info["ewpt_global_phase_path"] = "SINGLET_S -> X_BROKEN -> EW"
            point_info["ewpt_has_x_broken"] = True
            point_info["ewpt_ew_step_index"] = 2

        generator.run_ewpt_if_requested = fake_ewpt

        result = generator.evaluate_trsm_point_vxzero(
            123,
            380.0,
            500.0,
            200.0,
            -0.15,
            0.10,
            0.05,
            0.15,
            point_index=9,
        )

        self.assertEqual(result, 0)
        self.assertEqual(ewpt_calls, [(False, True)])
        self.assertEqual(rows[0][0], Path("output/trsm_points_unit_dm_failed.dat"))
        self.assertEqual(rows[0][1]["ewpt_ew_true_over_T"], 0.77)
        self.assertEqual(rows[0][1]["ewpt_global_phase_path"], "SINGLET_S -> X_BROKEN -> EW")
        self.assertTrue(rows[0][1]["ewpt_has_x_broken"])
        self.assertEqual(rows[0][1]["ewpt_ew_step_index"], 2)

    def test_write_all_points_records_failed_constraint_point_and_runs_ewpt(self):
        generator = load_generator_module()
        rows = []
        ewpt_calls = []

        def fake_write(path, point_info, mg5xsecs=None):
            rows.append((Path(path), dict(point_info)))

        generator.write_valid_point_file = fake_write
        generator.RunTag = "unit"
        generator.OutputDir = "output/"
        generator.cli_args = ewpt_args(run_ewpt=False, write_all_points=True)
        generator.np = SimpleNamespace(sin=lambda value: value)
        for name in [
            "Mz",
            "Mw",
            "Delta_S_central_wU",
            "Delta_T_central_wU",
            "Delta_U_central_wU",
            "errS_wU",
            "errT_wU",
            "errU_wU",
            "covST_wU",
            "covSU_wU",
            "covTU_wU",
            "pred",
            "H1",
            "H2",
            "H3",
        ]:
            setattr(generator, name, 0)
        generator.generate_lams = lambda *args, **kwargs: (
            200.0,
            0.0,
            380.0,
            500.0,
            -0.15,
            0.0,
            0.0,
            1.0,
            2.0,
            3.0,
            11.0,
            12.0,
            13.0,
            123.0,
            122.0,
            1111.0,
            1112.0,
            1113.0,
            133.0,
            0.99,
            -0.1,
            0.0,
            {},
            {},
            {},
            0.0,
            0.0,
            0.0,
        )
        generator.check_EWPO_wU = lambda *args, **kwargs: True
        generator.check_wmass_tania = lambda *args, **kwargs: True
        generator.analyze_parampoint = lambda *args, **kwargs: (False, True)
        generator.theory_constraints_vxzero = lambda *args, **kwargs: True
        generator.test_evo_vxzero = lambda *args, **kwargs: True
        generator.test_dm = lambda *args, **kwargs: (
            True,
            {},
            {"dm_mdm": 750.0, "dm_relic_excluded": False},
        )

        def fake_ewpt(point_info, passed, args, point_index, **kwargs):
            ewpt_calls.append((passed, kwargs.get("allow_any_failed")))
            point_info["ewpt_ew_true_over_T"] = 0.66
            point_info["ewpt_global_phase_path"] = "SYM -> EW"
            point_info["ewpt_has_x_broken"] = False
            point_info["ewpt_ew_step_index"] = 0

        generator.run_ewpt_if_requested = fake_ewpt

        result = generator.evaluate_trsm_point_vxzero(
            123,
            380.0,
            500.0,
            200.0,
            -0.15,
            0.10,
            0.05,
            0.15,
            point_index=9,
        )

        self.assertEqual(result, 0)
        self.assertEqual(ewpt_calls, [(False, True)])
        self.assertEqual(rows[0][0], Path("output/trsm_points_unit.dat"))
        self.assertFalse(rows[0][1]["hb"])
        self.assertEqual(rows[0][1]["ewpt_ew_true_over_T"], 0.66)
        self.assertEqual(rows[0][1]["ewpt_global_phase_path"], "SYM -> EW")


if __name__ == "__main__":
    unittest.main()
