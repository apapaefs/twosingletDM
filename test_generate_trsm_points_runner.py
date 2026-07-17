import importlib.util
import contextlib
import io
import json
import math
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
    generate_trsm_info.PORTAL_CONVENTION_ID = "trsm_vxzero_canonical_v1"

    def scalar_to_identical_scalar_width(daughter_mass, parent_mass, coupling):
        if parent_mass <= 2.0 * daughter_mass:
            return 0.0
        beta = math.sqrt(1.0 - 4.0 * daughter_mass**2 / parent_mass**2)
        return coupling**2 * beta / (8.0 * math.pi * parent_mass)

    def vxzero_portal_couplings(lphix, lsx, vs, a12, v=246.0):
        c12 = math.cos(a12)
        s12 = math.sin(a12)
        return (
            0.5 * (lphix * v * c12 - lsx * vs * s12),
            0.5 * (lphix * v * s12 + lsx * vs * c12),
        )

    def vxzero_invisible_decay_info(m1, m2, m3, k133, k233, base_w1, base_w2):
        gamma1 = scalar_to_identical_scalar_width(m3, m1, k133)
        gamma2 = scalar_to_identical_scalar_width(m3, m2, k233)
        w1 = base_w1 + gamma1
        w2 = base_w2 + gamma2
        return {
            "h1_h3h3_width": gamma1,
            "h1_h3h3_br": gamma1 / w1 if w1 else 0.0,
            "h2_h3h3_width": gamma2,
            "h2_h3h3_br": gamma2 / w2 if w2 else 0.0,
            "w1": w1,
            "w2": w2,
        }

    generate_trsm_info.scalar_to_identical_scalar_width = scalar_to_identical_scalar_width
    generate_trsm_info.vxzero_portal_couplings = vxzero_portal_couplings
    generate_trsm_info.vxzero_invisible_decay_info = vxzero_invisible_decay_info
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
        "write_dm_failed_to_main": False,
        "write_all_points": False,
        "write_evo_thc_points": False,
        "print_info": False,
        "higgstools_details": False,
        "higgstools_top": 5,
        "save_higgstools_details": True,
        "resonantDM1": False,
        "resonantDM2": False,
        "approximate_resonantDM": False,
        "independent_m3": False,
        "delta_res": None,
        "nrandom_count_evo_thc": False,
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
        self.assertIn("dm_indirect_channels_seen", scan_output.output_columns({}))
        self.assertIn("dm_indirect_channels_used", scan_output.output_columns({}))
        self.assertIn("dm_indirect_ratio", scan_output.output_columns({}))
        self.assertIn("dm_indirect_detection_excluded", scan_output.output_columns({}))
        self.assertIn("dm_limit_model", scan_output.output_columns({}))
        self.assertIn("dm_rescale", scan_output.output_columns({}))
        self.assertIn("ewpt_ew_true_over_T", scan_output.output_columns({}))
        self.assertIn("ewpt_global_phase_path", scan_output.output_columns({}))
        self.assertIn("ewpt_has_x_broken", scan_output.output_columns({}))
        self.assertIn("ewpt_ew_step_index", scan_output.output_columns({}))

    def test_vxzero_low_m3_is_recorded_as_stable_with_finite_outputs(self):
        import generate_trsm_info

        result = generate_trsm_info.get_point_info(
            246.0,
            200.0,
            0.0,
            125.09,
            380.0,
            8.5,
            -0.15,
            0.0,
            0.0,
            False,
            lX=0.10,
            lPhiX=0.05,
            lSX=0.15,
        )
        h3_brs = result[16]
        paramsubs = result[17]
        xs13_h3 = result[13]
        xs136_h3 = result[20]

        self.assertTrue(all(math.isfinite(float(value)) for value in h3_brs))
        self.assertTrue(all(float(value) == 0.0 for value in h3_brs))
        self.assertEqual(float(paramsubs["width3"]), 0.0)
        self.assertEqual(float(xs13_h3), 0.0)
        self.assertEqual(float(xs136_h3), 0.0)

    def test_parse_args_defaults_to_not_running_ewpt(self):
        generator = load_generator_module()

        args = generator.parse_args([])

        self.assertFalse(args.run_ewpt)
        self.assertFalse(args.write_all_points)
        self.assertFalse(args.write_evo_thc_points)
        self.assertFalse(args.print_info)
        self.assertTrue(args.save_higgstools_details)
        self.assertFalse(args.independent_m3)

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

    def test_parse_args_accepts_approximate_resonant_dm_option(self):
        generator = load_generator_module()

        args = generator.parse_args(["--approximate-resonantDM", "--delta-res", "15"])

        self.assertTrue(args.approximate_resonantDM)
        self.assertEqual(args.delta_res, 15.0)

    def test_parse_args_accepts_independent_m3_option(self):
        generator = load_generator_module()

        args = generator.parse_args(["--independent-m3"])

        self.assertTrue(args.independent_m3)

    def test_parse_args_rejects_independent_m3_with_other_mass_modes(self):
        generator = load_generator_module()

        conflicting_modes = [
            ["--resonantDM1"],
            ["--resonantDM2"],
            ["--approximate-resonantDM", "--delta-res", "10"],
        ]
        for mode in conflicting_modes:
            with self.subTest(mode=mode):
                with self.assertRaises(SystemExit):
                    generator.parse_args(["--independent-m3", *mode])

    def test_parse_args_rejects_independent_m3_in_explicit_point_mode(self):
        generator = load_generator_module()

        with self.assertRaises(SystemExit):
            generator.parse_args(
                [
                    "--independent-m3",
                    "--m2", "380",
                    "--m3", "100",
                    "--vs", "200",
                    "--a12", "-0.15",
                    "--lx", "0.10",
                    "--lphix", "0.05",
                    "--lsx", "0.15",
                ]
            )

    def test_parse_args_rejects_approximate_resonant_dm_without_delta(self):
        generator = load_generator_module()

        with self.assertRaises(SystemExit):
            generator.parse_args(["--approximate-resonantDM"])

    def test_parse_args_rejects_approximate_resonant_dm_with_exact_mode(self):
        generator = load_generator_module()

        with self.assertRaises(SystemExit):
            generator.parse_args(
                ["--approximate-resonantDM", "--delta-res", "15", "--resonantDM2"]
            )

    def test_parse_args_rejects_negative_delta_res(self):
        generator = load_generator_module()

        with self.assertRaises(SystemExit):
            generator.parse_args(["--approximate-resonantDM", "--delta-res", "-1"])

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
        self.assertAlmostEqual(generator.effective_m3(args, 7.526), 62.545)

    def test_effective_m3_uses_resonant_dm_modes(self):
        generator = load_generator_module()

        base = ewpt_args(resonantDM1=False, resonantDM2=False, m1=125.09, m3=600.0)
        self.assertEqual(generator.effective_m3(base, 200.0), 600.0)

        dm1 = ewpt_args(resonantDM1=True, resonantDM2=False, m1=125.09, m3=None)
        self.assertAlmostEqual(generator.effective_m3(dm1, 200.0), 62.545)

        dm2 = ewpt_args(resonantDM1=False, resonantDM2=True, m1=125.09, m3=None)
        self.assertEqual(generator.effective_m3(dm2, 200.0), 100.0)

    def test_resonant_dm2_sampler_respects_the_configured_m3_range(self):
        generator = load_generator_module()
        calls = []

        def fake_uniform(low, high):
            calls.append((low, high))
            return 200.0

        generator.random.uniform = fake_uniform
        args = ewpt_args(resonantDM2=True)

        m2, m3 = generator.sample_random_masses(args)

        self.assertEqual((m2, m3), (200.0, 100.0))
        self.assertEqual(
            calls,
            [
                (
                    max(generator.m2_min, 2.0 * generator.m3_min),
                    min(generator.m2_max, 2.0 * generator.m3_max),
                )
            ],
        )
        self.assertGreaterEqual(m3, generator.m3_min)
        self.assertLessEqual(m3, generator.m3_max)

    def test_approximate_resonant_mass_sampler_samples_both_branches(self):
        generator = load_generator_module()

        uniform_calls = []
        samples = iter([300.0, 604.0])
        generator.random.random = lambda: 0.25

        def fake_uniform(low, high):
            uniform_calls.append((low, high))
            return next(samples)

        generator.random.uniform = fake_uniform

        m2, m3 = generator.sample_approximate_resonant_masses(10.0)

        self.assertEqual((m2, m3), (604.0, 300.0))
        self.assertLessEqual(abs(m2 - 2.0 * m3), 10.0)
        self.assertEqual(uniform_calls, [(8, 505.0), (590.0, 610.0)])

        uniform_calls.clear()
        samples = iter([300.0, 62.0])
        generator.random.random = lambda: 0.75

        m2, m3 = generator.sample_approximate_resonant_masses(10.0)

        self.assertEqual((m2, m3), (300.0, 62.0))
        self.assertLessEqual(abs(2.0 * m3 - 125.09), 10.0)
        self.assertEqual(
            uniform_calls,
            [
                (generator.m2_min, generator.m2_max),
                ((125.09 - 10.0) / 2.0, (125.09 + 10.0) / 2.0),
            ],
        )

    def test_independent_m3_sampler_uses_full_configured_range(self):
        generator = load_generator_module()
        calls = []
        samples = iter([900.0, 100.0])

        def fake_uniform(low, high):
            calls.append((low, high))
            return next(samples)

        generator.random.uniform = fake_uniform
        args = ewpt_args(independent_m3=True)

        m2, m3 = generator.sample_random_masses(args)

        self.assertEqual((m2, m3), (900.0, 100.0))
        self.assertEqual(
            calls,
            [
                (generator.m2_min, generator.m2_max),
                (generator.m3_min, generator.m3_max),
            ],
        )
        self.assertLess(m3, m2 + generator.mhiggs)

    def test_default_mass_sampler_stays_inside_both_configured_ranges(self):
        generator = load_generator_module()
        calls = []
        samples = iter([800.0, 950.0])

        def fake_uniform(low, high):
            calls.append((low, high))
            return next(samples)

        generator.random.uniform = fake_uniform

        m2, m3 = generator.sample_random_masses(ewpt_args())

        self.assertEqual((m2, m3), (800.0, 950.0))
        self.assertEqual(
            calls,
            [
                (
                    generator.m2_min,
                    min(generator.m2_max, generator.m3_max - generator.mhiggs),
                ),
                (m2 + generator.mhiggs, generator.m3_max),
            ],
        )
        self.assertLessEqual(m2, generator.m2_max)
        self.assertLessEqual(m3, generator.m3_max)
        self.assertGreaterEqual(m3, m2 + generator.mhiggs)

    def test_independent_m3_random_scan_passes_sampled_pair_to_evaluator(self):
        generator = load_generator_module(["--independent-m3", "--nrandom", "1"])
        captured = []

        generator.nrandom = 1
        generator.cli_args = ewpt_args(independent_m3=True)
        generator.random.seed = lambda *args, **kwargs: None
        generator.random.uniform = lambda low, high: (
            generator.math.cos(0.2)
            if low == generator.k1_min and high == generator.k1_max
            else (low + high) / 2.0
        )
        generator.sample_random_masses = lambda args: (900.0, 100.0)
        generator.randsign = lambda: 1
        generator.np = SimpleNamespace(
            arccos=generator.math.acos,
            sqrt=generator.math.sqrt,
        )
        generator.round_sig = lambda value, _digits: value
        generator.tqdm = lambda iterable: iterable

        def fake_evaluate(myseed, m2, m3, *args, **kwargs):
            captured.append((m2, m3))
            return 1

        generator.evaluate_trsm_point_vxzero = fake_evaluate

        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            result = generator.run_random_vxzero_scan()

        self.assertEqual(result, 1)
        self.assertEqual(captured, [(900.0, 100.0)])
        self.assertIn("Sampling M3 independently", stdout.getvalue())
        self.assertIn(
            "Including h1/h2 -> h3 h3 invisible widths in HiggsTools inputs",
            stdout.getvalue(),
        )

    def test_k133_k233_inverse_map_round_trips_lambdas(self):
        generator = load_generator_module()
        vs = 500.0
        a12 = 0.2
        k133 = 3.0
        k233 = -4.0

        lphix, lsx = generator.k133_k233_to_lambdas(k133, k233, vs, a12)
        actual_k133, actual_k233 = generator.lambdas_to_k133_k233(
            lphix,
            lsx,
            vs,
            a12,
        )

        self.assertAlmostEqual(actual_k133, k133)
        self.assertAlmostEqual(actual_k233, k233)

    def test_valid_point_info_records_derived_k233(self):
        generator = load_generator_module()
        vs = 500.0
        a12 = 0.2
        lphix = 0.1
        lsx = -0.2

        _k133, expected_k233 = generator.lambdas_to_k133_k233(
            lphix,
            lsx,
            vs,
            a12,
        )
        point_info = generator.valid_point_info(
            300.0,
            600.0,
            vs,
            0.0,
            a12,
            0.0,
            0.0,
            0.3,
            lphix,
            lsx,
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
            0.1,
            0.0,
            True,
            True,
            True,
            True,
            True,
            True,
            True,
        )

        self.assertAlmostEqual(point_info["K233"], expected_k233)
        self.assertTrue(point_info["higgs_invisible_widths_included"])
        self.assertEqual(point_info["portal_convention"], "trsm_vxzero_canonical_v1")
        self.assertEqual(
            point_info["micromegas_model_convention"],
            "trsm_vxzero_canonical_v1",
        )
        self.assertEqual(point_info["h1_h3h3_width"], 0.0)
        self.assertEqual(point_info["h2_h3h3_width"], 0.0)

    def test_parse_args_accepts_write_all_points(self):
        generator = load_generator_module()

        args = generator.parse_args(["--write-all-points"])

        self.assertTrue(args.write_all_points)

    def test_parse_args_accepts_write_evo_thc_points(self):
        generator = load_generator_module()

        args = generator.parse_args(["--write-evo-thc-points"])

        self.assertTrue(args.write_evo_thc_points)

    def test_parse_args_rejects_write_all_and_write_evo_thc_together(self):
        generator = load_generator_module()

        with self.assertRaises(SystemExit):
            generator.parse_args(["--write-all-points", "--write-evo-thc-points"])

    def test_parse_args_accepts_count_evo_thc_option(self):
        generator = load_generator_module()

        args = generator.parse_args(["--nrandom-count-evo-thc"])

        self.assertTrue(args.nrandom_count_evo_thc)

    def test_parse_args_accepts_dm_failed_output_option(self):
        generator = load_generator_module()

        args = generator.parse_args(["--write-dm-failed"])

        self.assertTrue(args.write_dm_failed)

    def test_parse_args_accepts_dm_failed_main_output_option(self):
        generator = load_generator_module()

        args = generator.parse_args(["--write-dm-failed-to-main"])

        self.assertTrue(args.write_dm_failed_to_main)

    def test_parse_args_accepts_print_info_option(self):
        generator = load_generator_module()

        args = generator.parse_args(["--print-info"])

        self.assertTrue(args.print_info)

    def test_parse_args_accepts_k133_k233_scan_mode(self):
        generator = load_generator_module()

        args = generator.parse_args(["--scan-k133-k233"])

        self.assertTrue(args.scan_k133_k233)

    def test_parse_args_accepts_k133_k233_log_scan_mode(self):
        generator = load_generator_module()

        args = generator.parse_args(["--scan-k133-k233-log"])

        self.assertTrue(args.scan_k133_k233_log)

    def test_k133_k233_linear_and_log_ranges_are_separate(self):
        generator = load_generator_module()

        self.assertEqual(generator.K133_min, 1e-4)
        self.assertEqual(generator.K133_max, 8.0)
        self.assertEqual(generator.K233_min, 1e-4)
        self.assertEqual(generator.K233_max, 8.0)
        self.assertEqual(generator.K133_pow_min, -3)
        self.assertEqual(generator.K133_pow_max, 3)
        self.assertEqual(generator.K233_pow_min, -3)
        self.assertEqual(generator.K233_pow_max, 3)

    def test_parse_args_accepts_ewpt_on_dm_failed_option(self):
        generator = load_generator_module()

        args = generator.parse_args(["--run-ewpt-on-dm-failed"])

        self.assertTrue(args.run_ewpt_on_dm_failed)

    def test_k133_k233_scan_mode_derives_lambdas_for_random_points(self):
        generator = load_generator_module(["--scan-k133-k233", "--nrandom", "1"])
        captured = []
        samples = iter(
            [
                0.9800665778412416,  # k1 = cos(0.2)
                10.0,  # m2
                600.0,  # m3
                500.0,  # vs
                0.3,  # lX
                3.0,  # K133
                1.0,  # K133 sign factor
                4.0,  # K233
                -1.0,  # K233 sign factor
            ]
        )

        generator.nrandom = 1
        generator.cli_args = ewpt_args(
            scan_k133_k233=True,
            resonantDM1=False,
            resonantDM2=False,
        )
        generator.random.seed = lambda *args, **kwargs: None
        generator.random.uniform = lambda *args, **kwargs: next(samples)
        generator.randsign = lambda: 1
        generator.np = SimpleNamespace(
            arccos=generator.math.acos,
            sqrt=generator.math.sqrt,
        )
        generator.round_sig = lambda value, _digits: value
        generator.tqdm = lambda iterable: iterable

        def fake_evaluate(myseed, m2, m3, vs, a12, lx, lphix, lsx, **kwargs):
            captured.append(
                {
                    "m2": m2,
                    "m3": m3,
                    "vs": vs,
                    "a12": a12,
                    "lx": lx,
                    "lphix": lphix,
                    "lsx": lsx,
                    **kwargs,
                }
            )
            return 1

        generator.evaluate_trsm_point_vxzero = fake_evaluate

        result = generator.run_random_vxzero_scan()

        expected_lphix, expected_lsx = generator.k133_k233_to_lambdas(
            3.0,
            -4.0,
            500.0,
            0.2,
        )
        self.assertEqual(result, 1)
        self.assertEqual(len(captured), 1)
        self.assertAlmostEqual(captured[0]["a12"], 0.2)
        self.assertEqual(captured[0]["lx"], 0.3)
        self.assertAlmostEqual(captured[0]["lphix"], expected_lphix)
        self.assertAlmostEqual(captured[0]["lsx"], expected_lsx)

    def test_k133_k233_log_scan_mode_derives_lambdas_for_random_points(self):
        generator = load_generator_module(["--scan-k133-k233-log", "--nrandom", "1"])
        captured = []
        samples = iter(
            [
                0.9800665778412416,  # k1 = cos(0.2)
                10.0,  # m2
                600.0,  # m3
                500.0,  # vs
                0.3,  # lX
                1.0,  # K133 power -> 10 GeV
                -2.0,  # K233 power -> 0.01 GeV
            ]
        )
        signs = iter([1, -1, 1, 1])

        generator.nrandom = 1
        generator.cli_args = ewpt_args(
            scan_k133_k233_log=True,
            resonantDM1=False,
            resonantDM2=False,
        )
        generator.random.seed = lambda *args, **kwargs: None
        generator.random.uniform = lambda *args, **kwargs: next(samples)
        generator.randsign = lambda: next(signs)
        generator.np = SimpleNamespace(
            arccos=generator.math.acos,
            sqrt=generator.math.sqrt,
        )
        generator.round_sig = lambda value, _digits: value
        generator.tqdm = lambda iterable: iterable

        def fake_evaluate(myseed, m2, m3, vs, a12, lx, lphix, lsx, **kwargs):
            captured.append(
                {
                    "m2": m2,
                    "m3": m3,
                    "vs": vs,
                    "a12": a12,
                    "lx": lx,
                    "lphix": lphix,
                    "lsx": lsx,
                    **kwargs,
                }
            )
            return 1

        generator.evaluate_trsm_point_vxzero = fake_evaluate

        result = generator.run_random_vxzero_scan()

        expected_lphix, expected_lsx = generator.k133_k233_to_lambdas(
            10.0,
            -0.01,
            500.0,
            0.2,
        )
        self.assertEqual(result, 1)
        self.assertEqual(len(captured), 1)
        self.assertAlmostEqual(captured[0]["a12"], 0.2)
        self.assertEqual(captured[0]["lx"], 0.3)
        self.assertAlmostEqual(captured[0]["lphix"], expected_lphix)
        self.assertAlmostEqual(captured[0]["lsx"], expected_lsx)

    def test_approximate_resonant_scan_uses_sampled_mass_pair(self):
        generator = load_generator_module(
            ["--approximate-resonantDM", "--delta-res", "10", "--nrandom", "1"]
        )
        captured = []

        generator.nrandom = 1
        generator.cli_args = ewpt_args(
            approximate_resonantDM=True,
            delta_res=10.0,
        )
        generator.random.seed = lambda *args, **kwargs: None
        generator.random.uniform = lambda low, high: (
            0.9800665778412416 if low == generator.k1_min and high == generator.k1_max else low
        )
        generator.sample_approximate_resonant_masses = (
            lambda delta, m1=125.09: (604.0, 300.0)
        )
        generator.randsign = lambda: 1
        generator.np = SimpleNamespace(
            arccos=generator.math.acos,
            sqrt=generator.math.sqrt,
        )
        generator.round_sig = lambda value, _digits: value
        generator.tqdm = lambda iterable: iterable

        def fake_evaluate(myseed, m2, m3, vs, a12, lx, lphix, lsx, **kwargs):
            captured.append(
                {
                    "m2": m2,
                    "m3": m3,
                    "vs": vs,
                    "a12": a12,
                    "lx": lx,
                    "lphix": lphix,
                    "lsx": lsx,
                    **kwargs,
                }
            )
            return 1

        generator.evaluate_trsm_point_vxzero = fake_evaluate

        result = generator.run_random_vxzero_scan()

        self.assertEqual(result, 1)
        self.assertEqual(len(captured), 1)
        self.assertEqual(captured[0]["m2"], 604.0)
        self.assertEqual(captured[0]["m3"], 300.0)

    def test_random_scan_count_evo_thc_runs_until_target_count(self):
        generator = load_generator_module(["--nrandom-count-evo-thc", "--nrandom", "2"])
        captured = []
        statuses = iter(
            [
                {"viable": 0, "evo_thc": False},
                {"viable": 0, "evo_thc": True},
                {"viable": 1, "evo_thc": True},
            ]
        )

        generator.nrandom = 2
        generator.cli_args = ewpt_args(nrandom_count_evo_thc=True)
        generator.random.seed = lambda *args, **kwargs: None
        generator.random.uniform = lambda low, high: (
            0.9800665778412416 if low == generator.k1_min and high == generator.k1_max else low
        )
        generator.randsign = lambda: 1
        generator.np = SimpleNamespace(
            arccos=generator.math.acos,
            sqrt=generator.math.sqrt,
        )
        generator.round_sig = lambda value, _digits: value
        generator.tqdm = lambda iterable, **kwargs: iterable

        def fake_evaluate(myseed, m2, m3, vs, a12, lx, lphix, lsx, **kwargs):
            captured.append(dict(kwargs))
            return next(statuses)

        generator.evaluate_trsm_point_vxzero = fake_evaluate

        result = generator.run_random_vxzero_scan()

        self.assertEqual(result, 1)
        self.assertEqual(len(captured), 3)
        self.assertEqual([call["point_index"] for call in captured], [1, 2, 3])
        self.assertTrue(all(call["return_status"] for call in captured))

    def test_print_info_count_evo_thc_scan_reports_progress_counters(self):
        generator = load_generator_module(["--nrandom-count-evo-thc", "--print-info", "--nrandom", "2"])
        statuses = iter(
            [
                {"viable": 0, "evo_thc": False},
                {"viable": 0, "evo_thc": True},
                {"viable": 1, "evo_thc": True},
            ]
        )

        generator.nrandom = 2
        generator.cli_args = ewpt_args(nrandom_count_evo_thc=True, print_info=True)
        generator.random.seed = lambda *args, **kwargs: None
        generator.random.uniform = lambda low, high: (
            0.9800665778412416 if low == generator.k1_min and high == generator.k1_max else low
        )
        generator.randsign = lambda: 1
        generator.np = SimpleNamespace(
            arccos=generator.math.acos,
            sqrt=generator.math.sqrt,
        )
        generator.round_sig = lambda value, _digits: value
        generator.tqdm = lambda iterable, **kwargs: iterable
        generator.evaluate_trsm_point_vxzero = lambda *args, **kwargs: next(statuses)

        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            result = generator.run_random_vxzero_scan()

        output = stdout.getvalue()
        self.assertEqual(result, 1)
        self.assertIn("Progress: draws=1 evo/thc=0/2 viable=0", output)
        self.assertIn("Progress: draws=2 evo/thc=1/2 viable=0", output)
        self.assertIn("Progress: draws=3 evo/thc=2/2 viable=1", output)

    def test_evaluate_vxzero_status_reports_evo_thc_before_later_failures(self):
        generator = load_generator_module()

        generator.cli_args = ewpt_args(nrandom_count_evo_thc=True, run_ewpt=False)
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

        def fail_if_dm_runs(*args, **kwargs):
            raise AssertionError("DM should not run after a later non-DM failure")

        generator.test_dm = fail_if_dm_runs

        status = generator.evaluate_trsm_point_vxzero(
            123,
            380.0,
            500.0,
            200.0,
            -0.15,
            0.10,
            0.05,
            0.15,
            return_status=True,
        )

        self.assertEqual(status, {"viable": 0, "evo_thc": True})

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

    def test_dm_failed_main_option_writes_otherwise_allowed_point_to_main_output(self):
        generator = load_generator_module()
        records = []

        def fake_write(path, point_info, mg5xsecs=None):
            records.append((Path(path), dict(point_info)))

        generator.write_valid_point_file = fake_write
        generator.RunTag = "unit"
        generator.OutputDir = "output/"
        generator.cli_args = ewpt_args(
            run_ewpt=False,
            write_dm_failed_to_main=True,
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
            {
                "dm_mdm": 750.0,
                "dm_relic_excluded": True,
                "dm_limit_model": "lz2025-source",
                "dm_rescale": True,
            },
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
        self.assertEqual(records[0][0], Path("output/trsm_points_unit.dat"))
        self.assertEqual(records[0][1]["dm"], False)
        self.assertEqual(records[0][1]["dm_mdm"], 750.0)
        self.assertEqual(records[0][1]["dm_limit_model"], "lz2025-source")
        self.assertTrue(records[0][1]["dm_rescale"])

    def test_print_info_flag_prints_dm_failed_point_details(self):
        generator = load_generator_module()

        generator.write_valid_point_file = lambda *args, **kwargs: None
        generator.RunTag = "unit"
        generator.OutputDir = "output/"
        generator.cli_args = SimpleNamespace(
            run_ewpt=False,
            write_dm_failed=True,
            run_ewpt_on_dm_failed=False,
            print_info=True,
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
            "DM check: Fail\n  Reason: test failure",
            {"dm_mdm": 750.0, "dm_relic_excluded": True},
        )
        generator.print_dm_info = lambda info: print(info)

        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
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

        output = stdout.getvalue()
        expected_k233 = generator.lambdas_to_k133_k233(0.05, 0.15, 200.0, -0.15)[1]
        self.assertEqual(result, 0)
        self.assertIn("['M2', 380.0]", output)
        self.assertIn(f"['K233', {expected_k233}]", output)
        self.assertIn("['dm', 'Fail']", output)
        self.assertIn("DM check: Fail", output)
        self.assertIn("Reason: test failure", output)

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
            [
                {
                    "print_details": True,
                    "details_top": 7,
                    "return_details": True,
                    "h1_direct_invisible_width": 0.0,
                    "h2_direct_invisible_width": 0.0,
                }
            ],
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

    def test_write_evo_thc_points_records_evo_thc_point_with_later_failure(self):
        generator = load_generator_module()
        rows = []

        def fake_write(path, point_info, mg5xsecs=None):
            rows.append((Path(path), dict(point_info)))

        generator.write_valid_point_file = fake_write
        generator.RunTag = "unit"
        generator.OutputDir = "output/"
        generator.cli_args = ewpt_args(run_ewpt=False, write_evo_thc_points=True)
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
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0][0], Path("output/trsm_points_unit.dat"))
        self.assertTrue(rows[0][1]["evo"])
        self.assertTrue(rows[0][1]["thc"])
        self.assertFalse(rows[0][1]["hb"])

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
