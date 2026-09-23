import importlib.util
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


SCRIPT_PATH = Path(__file__).resolve().parent / "plot_trsm_constraint_suite.py"


def load_module():
    spec = importlib.util.spec_from_file_location("plot_trsm_constraint_suite", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


HEADER = [
    "M2",
    "M3",
    "a12",
    "K133",
    "K233",
    "lSX",
    "evo",
    "thc",
    "hb",
    "hs",
    "ewpo",
    "wmass",
    "dm",
    "dm_omega",
    "dm_relic_upper_limit",
    "dm_dir_det",
    "dm_dir_det_limit",
    "dm_relic_excluded",
    "dm_direct_detection_excluded",
    "dm_indirect_available",
    "dm_indirect_ratio",
    "dm_indirect_detection_excluded",
    "higgstools_hs_delta_chi2",
    "ewpt_ew_true_over_T",
    "w2",
]

ROWS = [
    [200, 325, -0.1, -100, 100, -2, True, True, False, False, True, False, False, 0.2, 0.1, 2e-9, 1e-9, True, True, False, 0, False, 25, "nan", 1e-8],
    [250, 500, 0.2, -10, 10, -1, True, True, True, True, True, False, True, 0.05, 0.1, 0.5e-9, 1e-9, False, False, False, 0, False, 16, "nan", 1e-6],
    [300, 600, -0.3, -1, 1, -0.1, True, True, True, True, True, True, False, 0.3, 0.1, 0.25e-9, 1e-9, True, False, True, 0.25, False, 9, "nan", 1e-4],
    [350, 700, 0.4, -0.1, 0.1, 0, True, True, True, True, True, True, True, 0.02, 0.1, 0.1e-9, 1e-9, False, False, True, 0.75, False, 4, "nan", 1e-2],
    [400, 800, -0.5, 0.1, -0.1, 0.1, True, True, False, True, True, True, False, 0.08, 0.1, 4e-9, 1e-9, False, True, True, 1.25, True, 1, "nan", 1e-1],
    [450, 900, 0.6, 1, -1, 1, True, True, True, False, True, True, True, 0.01, 0.1, 0.2e-9, 1e-9, False, False, False, 0, False, 0.25, "nan", 1],
    [500, 1000, -0.7, 10, -10, 2, True, True, True, True, True, True, False, 0.4, 0.1, 3e-9, 1e-9, True, True, True, 0.5, False, 0, "nan", 10],
    [550, 1100, 0.8, 100, -100, 10, True, True, True, True, True, True, True, 0.1, 0.1, 1e-9, 1e-9, False, False, True, 1.0, False, "nan", "nan", 100],
]

BSMPT_EXTRA_HEADER = [
    "ewpt_ew_jump_over_T",
    "ewpt_status",
    "ewpt_global_phase_path",
    "ewpt_has_x_broken",
    "ewpt_ew_step_index",
]
BSMPT_HEADER = HEADER + BSMPT_EXTRA_HEADER

RATE_EXTRA_HEADER = [
    "xs136_lo_h1_pb",
    "xs136_lo_h2_pb",
    "h1_h2h2_br",
    "h2_h1h1_br",
    "h1_h3h3_br",
    "h2_h3h3_br",
    "mg5_xsec_gg_heta0_pb",
    "mg5_xsec_pp_eta0Z_pb",
]
RATE_HEADER = HEADER + RATE_EXTRA_HEADER


def rate_fixture_rows():
    rows = []
    for index, original in enumerate(ROWS, start=1):
        rows.append(
            list(original)
            + [
                20.0 / index,
                1.0 / index,
                0.02 * index,
                0.03 * index,
                0.04,
                0.05,
                0.2 / index,
                0.1 / index,
            ]
        )
    return rows


SIGNAL_EXTRA_HEADER = ["k2", "h2_h3h3_br"]
SIGNAL_HEADER = HEADER + SIGNAL_EXTRA_HEADER


def bsmpt_fixture_rows():
    payloads = [
        ("nan", "nan", "nan", "nan", "nan", "nan"),
        ("nan", "nan", "failed", "nan", "nan", "nan"),
        ("nan", "nan", "success", "SYM -> EW", False, 0),
        (0.75, 0.65, "success", "SYM -> SINGLET_S -> EW", False, 1),
        (1.25, 1.10, "success", "SYM -> X_BROKEN -> EW_X_BROKEN", True, 1),
        (2.0, 1.80, "success", "SYM -> EW", False, 0),
        (0.4, 0.20, "success", "SYM -> MIXED -> EW", False, 1),
        ("nan", "nan", "nan", "nan", "nan", "nan"),
    ]
    rows = []
    for original, payload in zip(ROWS, payloads):
        strength, jump, status, phase_path, has_x_broken, ew_step = payload
        row = list(original)
        row[HEADER.index("ewpt_ew_true_over_T")] = strength
        row.extend([jump, status, phase_path, has_x_broken, ew_step])
        rows.append(row)
    return rows


def signal_fixture_rows():
    specifications = [
        (200.0, 50.0, 0.20, 0.10, 1.0),
        (500.0, 100.0, -0.10, 0.50, 20.0),
        (900.0, 300.0, 0.05, 0.80, 135.0),
    ]
    rows = []
    for m2, m3, k2, branching_ratio, width in specifications:
        row = list(ROWS[-1])
        row[HEADER.index("M2")] = m2
        row[HEADER.index("M3")] = m3
        row[HEADER.index("w2")] = width
        for column in ("evo", "thc", "hb", "hs", "ewpo", "wmass", "dm"):
            row[HEADER.index(column)] = True
        row[HEADER.index("dm_relic_excluded")] = False
        row[HEADER.index("dm_direct_detection_excluded")] = False
        row[HEADER.index("dm_indirect_available")] = True
        row[HEADER.index("dm_indirect_detection_excluded")] = False
        row.extend([k2, branching_ratio])
        rows.append(row)
    return rows


def write_fixture(path, header=None, rows=None):
    if header is None:
        header = HEADER
    if rows is None:
        rows = ROWS
    lines = ["\t".join(str(value) for value in header)]
    lines.extend("\t".join(str(value) for value in row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


class TestPlotTRSMConstraintSuite(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plotter = load_module()

    def load_fixture(self, directory):
        path = Path(directory) / "points.tsv"
        write_fixture(path)
        return self.plotter.load_scan(path)

    def load_bsmpt_fixture(self, directory):
        path = Path(directory) / "points_bsmpt.tsv"
        write_fixture(path, BSMPT_HEADER, bsmpt_fixture_rows())
        return self.plotter.load_scan(path)

    def load_rate_fixture(self, directory):
        path = Path(directory) / "points_rates.tsv"
        write_fixture(path, RATE_HEADER, rate_fixture_rows())
        return self.plotter.load_scan(path)

    def load_signal_fixture(self, directory):
        path = Path(directory) / "points_signal.tsv"
        write_fixture(path, SIGNAL_HEADER, signal_fixture_rows())
        return self.plotter.load_scan(path)

    def test_strict_bool_accepts_only_canonical_tokens(self):
        self.assertIs(self.plotter.strict_bool("True"), True)
        self.assertIs(self.plotter.strict_bool("False"), False)
        for value in ["true", "FALSE", "1", "0", "yes", "", " True ", "nan"]:
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "expected 'True' or 'False'"):
                    self.plotter.strict_bool(value, "dm", 9)

    def test_nullable_boolean_accepts_nan_and_legacy_blank(self):
        self.assertIs(self.plotter.strict_nullable_bool("True"), True)
        self.assertIs(self.plotter.strict_nullable_bool("False"), False)
        self.assertIsNone(self.plotter.strict_nullable_bool("nan"))
        self.assertIsNone(self.plotter.strict_nullable_bool(""))
        for value in ["NaN", "none", " nan "]:
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "expected 'True' or 'False'"):
                    self.plotter.strict_nullable_bool(
                        value, "dm_relic_excluded", 384
                    )

    def test_load_scan_rejects_missing_column_and_invalid_boolean(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            missing_path = Path(tmpdir) / "missing.tsv"
            missing_header = [name for name in HEADER if name != "M3"]
            missing_rows = [
                [value for index, value in enumerate(row) if HEADER[index] != "M3"]
                for row in ROWS
            ]
            write_fixture(missing_path, missing_header, missing_rows)
            with self.assertRaisesRegex(ValueError, "Missing required input columns: M3"):
                self.plotter.load_scan(missing_path)

            invalid_path = Path(tmpdir) / "invalid.tsv"
            invalid_rows = [list(row) for row in ROWS]
            invalid_rows[0][HEADER.index("dm")] = "TRUE"
            write_fixture(invalid_path, HEADER, invalid_rows)
            with self.assertRaisesRegex(ValueError, "dm on row 2"):
                self.plotter.load_scan(invalid_path)

    def test_historical_blank_results_remain_unavailable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "historical.tsv"
            rows = [list(row) for row in ROWS]
            rows[0][HEADER.index("dm_direct_detection_excluded")] = ""
            rows[0][HEADER.index("ewpt_ew_true_over_T")] = ""
            write_fixture(path, HEADER, rows)
            data = self.plotter.load_scan(path)
            self.assertTrue(np.isnan(data.f("ewpt_ew_true_over_T")[0]))
            self.assertFalse(data.b("dm_direct_detection_excluded")[0])
            rows[0][HEADER.index("M2")] = ""
            write_fixture(path, HEADER, rows)
            with self.assertRaisesRegex(ValueError, "M2 on row 2"):
                self.plotter.load_scan(path)

    def test_rate_columns_are_derived_and_rendered_on_log_scale(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_rate_fixture(tmpdir)

        np.testing.assert_allclose(
            data.f("mono_higgs_xsec_pb"),
            data.f("mg5_xsec_gg_heta0_pb") * data.f("h2_h3h3_br"),
        )
        expected_cascade = (
            data.f("xs136_lo_h2_pb")
            * data.f("h2_h1h1_br")
            * 2.0
            * data.f("h1_h3h3_br")
            * (1.0 - data.f("h1_h3h3_br"))
        )
        np.testing.assert_allclose(
            data.f("xsec_h2_h1h1_one_h1_invisible_pb"), expected_cascade
        )
        spec = self.plotter.PLOT_BY_STEM["52_mono_higgs_xsec_vs_m2"]
        self.assertIsNone(self.plotter.spec_unavailable_reason(data, spec))
        fig, ax = self.plotter.plt.subplots()
        try:
            self.plotter.render_rate_xy(ax, data, spec)
            self.assertEqual(ax.get_yscale(), "log")
            self.assertGreater(len(ax.collections), 0)
        finally:
            self.plotter.plt.close(fig)

        no_dm_spec = self.plotter.PLOT_BY_STEM[
            "60_mono_higgs_xsec_no_dm_vs_m2"
        ]
        self.assertIsNone(
            self.plotter.spec_unavailable_reason(data, no_dm_spec)
        )
        fig, ax = self.plotter.plt.subplots()
        try:
            self.plotter.render_rate_xy(ax, data, no_dm_spec)
            self.assertEqual(ax.get_yscale(), "log")
            self.assertIn("All non-DM constraints", ax.get_title())
        finally:
            self.plotter.plt.close(fig)

        paths = self.plotter.expected_figure_paths(Path("plots"), "both", data=data)
        self.assertEqual(len(paths), 132)
        self.assertTrue(
            any(path.stem == "dashboard_scalar_cascade_rates" for path in paths)
        )
        self.assertTrue(any(path.stem == "dashboard_mg5_mono_rates" for path in paths))
        self.assertTrue(
            any(
                path.stem == "dashboard_scalar_cascade_rates_no_dm"
                for path in paths
            )
        )
        self.assertTrue(
            any(path.stem == "dashboard_mg5_mono_rates_no_dm" for path in paths)
        )

    def test_scan_metadata_is_loaded_and_rendered_with_configured_ranges(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_path = tmpdir / "points.tsv"
            write_fixture(input_path)
            metadata_path = input_path.with_suffix(".metadata.json")
            metadata = {
                "schema": "trsm_scan_metadata_v1",
                "created_utc": "2026-07-18T12:00:00+00:00",
                "seed": 888,
                "requested_points": 10000,
                "stopping_rule": "Stop after evo/thc target",
                "output_selection": "Points passing evo and thc",
                "mass_sampling": {
                    "mode": "independent_m3",
                    "description": "M2 and M3 are sampled independently.",
                },
                "portal_sampling": {
                    "mode": "k133_k233_log",
                    "description": "Signed logarithmic K scan.",
                },
                "portal_convention": "trsm_vxzero_canonical_v1",
                "variable_ranges": [
                    {
                        "variable": "M2",
                        "column": "M2",
                        "configured_min": 4,
                        "configured_max": 1000,
                        "effective_min": 4,
                        "effective_max": 1000,
                        "unit": "GeV",
                        "sampling": "uniform",
                        "note": "",
                    },
                    {
                        "variable": "K133",
                        "column": "K133",
                        "configured_min": -1000,
                        "configured_max": 1000,
                        "effective_min": -1000,
                        "effective_max": 1000,
                        "unit": "GeV",
                        "sampling": "log-uniform magnitude with random sign",
                        "note": "|K133| >= 0.001 GeV",
                    },
                ],
                "fixed_parameters": [
                    {
                        "variable": "vx",
                        "value": 0,
                        "unit": "GeV",
                        "note": "Dark-matter branch",
                    }
                ],
                "command_line": [
                    "generate_trsm_points.py",
                    "888",
                    "--independent-m3",
                ],
                "options": {
                    "independent_m3": True,
                    "nrandom": 10000,
                    "run_ewpt": True,
                    "ewpt_thigh": 1000,
                },
            }
            metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

            data = self.plotter.load_scan(input_path)
            rendered = self.plotter.scan_information_html(data)

        self.assertEqual(data.metadata["seed"], 888)
        self.assertEqual(data.metadata_source, metadata_path)
        self.assertIn("Configured scan metadata", rendered)
        self.assertIn("independent_m3", rendered)
        self.assertIn("4 – 1000", rendered)
        self.assertIn("200 – 550", rendered)
        self.assertIn("generate_trsm_points.py 888 --independent-m3", rendered)
        self.assertIn("Dark-matter branch", rendered)
        self.assertIn("BSMPT requested", rendered)
        self.assertIn("BSMPT high temperature", rendered)

    def test_explicit_invalid_scan_metadata_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_path = tmpdir / "points.tsv"
            metadata_path = tmpdir / "broken.json"
            write_fixture(input_path)
            metadata_path.write_text("not-json", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "Could not read scan metadata"):
                self.plotter.load_scan(input_path, metadata_path)

    def test_composite_masks_and_four_way_categories(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_fixture(tmpdir)

        np.testing.assert_array_equal(data.b("theory"), np.ones(8, dtype=bool))
        np.testing.assert_array_equal(
            data.b("experimental"),
            [False, False, True, True, False, False, True, True],
        )
        np.testing.assert_array_equal(
            data.b("dm"), [False, True, False, True, False, True, False, True]
        )
        np.testing.assert_array_equal(
            data.b("non_dm_viability"),
            [False, False, True, True, False, False, True, True],
        )
        np.testing.assert_array_equal(
            data.b("full_viability"),
            [False, False, False, True, False, False, False, True],
        )
        self.assertEqual(
            data.derived["fourway"].tolist(),
            [
                "neither",
                "DM only",
                "experimental only",
                "both",
                "neither",
                "DM only",
                "experimental only",
                "both",
            ],
        )

    def test_bsmpt_results_are_detected_and_classified(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_bsmpt_fixture(tmpdir)

        self.assertTrue(self.plotter.has_bsmpt_results(data))
        self.assertEqual(data.s("ewpt_status")[0], "")
        self.assertEqual(
            data.derived["bsmpt_status"].tolist(),
            [
                "not run",
                "failed",
                "success / no selected FOPT",
                "legacy v/T ≤ 1 diagnostic",
                "legacy v/T > 1 diagnostic",
                "legacy v/T > 1 diagnostic",
                "legacy v/T ≤ 1 diagnostic",
                "not run",
            ],
        )
        self.assertEqual(
            data.derived["bsmpt_phase"].tolist(),
            [
                "not run",
                "failed",
                "direct EW",
                "singlet-assisted",
                "X-broken",
                "direct EW",
                "other / multistep",
                "not run",
            ],
        )
        self.assertEqual(
            data.derived["bsmpt_step"].tolist(),
            [
                "not run",
                "failed",
                "step 0",
                "step 1",
                "step 1",
                "step 0",
                "step 1",
                "not run",
            ],
        )
        self.assertEqual(np.count_nonzero(data.b("bsmpt_attempted")), 6)
        self.assertEqual(np.count_nonzero(data.b("bsmpt_success")), 5)
        self.assertEqual(np.count_nonzero(data.b("bsmpt_failed")), 1)
        self.assertEqual(np.count_nonzero(data.b("bsmpt_selected_fopt")), 4)
        self.assertEqual(np.count_nonzero(data.b("bsmpt_weak_fopt")), 2)
        self.assertEqual(np.count_nonzero(data.b("bsmpt_strong_fopt")), 2)
        self.assertEqual(np.count_nonzero(data.b("bsmpt_x_broken_path")), 1)
        self.assertEqual(np.count_nonzero(data.b("bsmpt_direct_ew_entry")), 2)
        self.assertEqual(np.count_nonzero(data.b("bsmpt_multistep_ew_entry")), 3)

    def test_legacy_scan_has_no_bsmpt_attempts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_fixture(tmpdir)

        self.assertFalse(self.plotter.has_bsmpt_results(data))
        self.assertTrue(np.all(data.s("ewpt_status") == ""))
        self.assertTrue(np.all(data.derived["bsmpt_status"] == "not run"))
        self.assertTrue(np.all(data.derived["bsmpt_phase"] == "not run"))
        self.assertTrue(np.all(data.derived["bsmpt_step"] == "not run"))

    def test_bsmpt_nullable_boolean_is_strict_when_present(self):
        rows = bsmpt_fixture_rows()
        rows[2][BSMPT_HEADER.index("ewpt_has_x_broken")] = "yes"
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "bad_bsmpt_bool.tsv"
            write_fixture(path, BSMPT_HEADER, rows)
            with self.assertRaisesRegex(
                ValueError, "ewpt_has_x_broken on row 4"
            ):
                self.plotter.load_scan(path)

    def test_each_composite_input_flag_is_required(self):
        base_row = list(ROWS[-1])
        for column in ["evo", "thc", "hb", "hs", "ewpo", "wmass"]:
            with self.subTest(column=column), tempfile.TemporaryDirectory() as tmpdir:
                row = list(base_row)
                row[HEADER.index(column)] = False
                path = Path(tmpdir) / "one.tsv"
                write_fixture(path, rows=[row])
                data = self.plotter.load_scan(path)
                if column in {"evo", "thc"}:
                    self.assertFalse(data.b("theory")[0])
                    self.assertTrue(data.b("experimental")[0])
                else:
                    self.assertTrue(data.b("theory")[0])
                    self.assertFalse(data.b("experimental")[0])
                self.assertFalse(data.b("full_viability")[0])

    def test_cumulative_constraint_masks_match_collaborator_sequence(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_fixture(tmpdir)

        masks = self.plotter.cumulative_constraint_masks(data)
        self.assertEqual(
            list(masks),
            ["all", "hb", "hb_hs", "hb_hs_wmass", "hb_hs_wmass_dm"],
        )
        self.assertEqual(
            [int(np.count_nonzero(mask)) for mask in masks.values()],
            [8, 6, 5, 4, 2],
        )
        previous = masks["all"]
        for mask in list(masks.values())[1:]:
            self.assertTrue(np.all(mask <= previous))
            previous = mask

        no_ewpo_row = list(ROWS[-1])
        no_ewpo_row[HEADER.index("ewpo")] = False
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "no-ewpo.tsv"
            write_fixture(path, rows=[no_ewpo_row])
            no_ewpo_data = self.plotter.load_scan(path)
        self.assertFalse(no_ewpo_data.b("experimental")[0])
        self.assertTrue(
            no_ewpo_data.b("dm")[0]
            and self.plotter.cumulative_constraint_masks(no_ewpo_data)[
                "hb_hs_wmass_dm"
            ][0]
        )

    def test_cmb_failures_and_legacy_unassessed_summary(self):
        from trsm_cmb import CMB_COLUMNS, CMBSignal, assess_cmb_limit, cmb_diagnostics
        records = []
        for signal in (CMBSignal(True, 8, "ok", ""), CMBSignal(), CMBSignal(True, 0, "ok", "")):
            row = list(ROWS[1])
            result = assess_cmb_limit(signal, .05)
            row[HEADER.index("dm")] = result.passed
            records.append(row + ["nan" if value is None else value for value in cmb_diagnostics(result).values()])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cmb.tsv"
            write_fixture(path, header=HEADER + list(CMB_COLUMNS), rows=records)
            data = self.plotter.load_scan(path)
            legacy = self.load_fixture(tmp)
        self.assertEqual(data.derived["dm_failure"].tolist(), ["CMB only", "CMB unavailable", "pass"])
        summary = {row.metric: row for row in self.plotter.build_summary(data, [])}
        self.assertEqual(summary["dm_component_mismatch"].count, 0)
        self.assertEqual(summary["cmb_excluded"].count, 1)
        self.assertEqual(summary["cmb_unavailable"].count, 1)
        self.assertEqual(summary["cmb_unassessed"].count, 0)
        legacy_summary = {row.metric: row for row in self.plotter.build_summary(legacy, [])}
        self.assertEqual(legacy_summary["cmb_unassessed"].count, len(legacy))

    def test_dm_failure_and_indirect_categories(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_fixture(tmpdir)

        self.assertEqual(
            data.derived["dm_failure"].tolist(),
            [
                "relic + direct",
                "pass",
                "relic only",
                "pass",
                "direct only",
                "pass",
                "relic + direct",
                "pass",
            ],
        )
        self.assertEqual(
            data.derived["indirect"].tolist(),
            [
                "unavailable",
                "unavailable",
                "allowed",
                "allowed",
                "excluded",
                "unavailable",
                "allowed",
                "allowed",
            ],
        )
        indirect = data.f("indirect_ratio_available")
        self.assertTrue(math.isnan(indirect[0]))
        self.assertTrue(math.isnan(indirect[1]))
        self.assertEqual(indirect[2], 0.25)
        self.assertNotEqual(data.derived["indirect"][0], "allowed")

    def test_unavailable_dm_component_results_are_never_counted_as_passing(self):
        unavailable = list(ROWS[0])
        unavailable[HEADER.index("dm")] = False
        for column in (
            "dm_relic_excluded",
            "dm_direct_detection_excluded",
            "dm_indirect_available",
            "dm_indirect_detection_excluded",
        ):
            unavailable[HEADER.index(column)] = "nan"
        for column in (
            "dm_omega",
            "dm_relic_upper_limit",
            "dm_dir_det",
            "dm_dir_det_limit",
            "dm_indirect_ratio",
        ):
            unavailable[HEADER.index(column)] = "nan"

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "unavailable.tsv"
            write_fixture(path, rows=[unavailable])
            data = self.plotter.load_scan(path)
            summary = {
                row.metric: row for row in self.plotter.build_summary(data)
            }

        self.assertFalse(data.b("dm_result_available")[0])
        self.assertFalse(data.b("relic_available")[0])
        self.assertFalse(data.b("direct_available")[0])
        self.assertFalse(data.b("relic_pass")[0])
        self.assertFalse(data.b("direct_pass")[0])
        self.assertEqual(data.derived["dm_failure"][0], "DM unavailable")
        self.assertEqual(data.derived["indirect"][0], "DM unavailable")
        relic_categories, relic_styles = self.plotter.category_styles(
            data, "relic_pass"
        )
        self.assertEqual(relic_categories[0], "unavailable")
        self.assertIn("unavailable", relic_styles)
        self.assertEqual(summary["dm_result_available"].count, 0)
        self.assertEqual(summary["dm_result_unavailable"].count, 1)
        self.assertEqual(summary["relic_pass"].count, 0)
        self.assertEqual(summary["relic_pass"].denominator, 0)
        self.assertEqual(summary["direct_pass"].count, 0)
        self.assertEqual(summary["dm_failure_DM_unavailable"].count, 1)
        self.assertEqual(summary["indirect_dm_unavailable"].count, 1)
        self.assertEqual(summary["dm_component_mismatch"].denominator, 0)

    def test_safe_ratio_and_positive_log10(self):
        numerator = np.array([2.0, 5.0, 0.0, 1.0, np.nan, 1.0])
        denominator = np.array([1.0, 10.0, 1.0, 0.0, 1.0, np.inf])
        ratios = self.plotter.safe_ratio(numerator, denominator)
        np.testing.assert_allclose(ratios[:3], [2.0, 0.5, 0.0])
        self.assertTrue(np.all(np.isnan(ratios[3:])))

        logs = self.plotter.positive_log10(np.array([0.1, 1.0, 10.0, 0.0, -1.0, np.nan]))
        np.testing.assert_allclose(logs[:3], [-1.0, 0.0, 1.0])
        self.assertTrue(np.all(np.isnan(logs[3:])))

    def test_yr4_signal_grid_and_interpolation(self):
        grid = self.plotter.load_yr4_cross_section_grid()
        self.assertEqual(len(grid.mass_gev), 114)
        self.assertEqual((grid.mass_gev[0], grid.mass_gev[-1]), (10.0, 3000.0))
        index_125 = int(np.flatnonzero(grid.mass_gev == 125.0)[0])
        index_500 = int(np.flatnonzero(grid.mass_gev == 500.0)[0])
        self.assertAlmostEqual(grid.ggf_pb[index_125], 45.142)
        self.assertAlmostEqual(grid.vbf_pb[index_125], 4.237)
        self.assertAlmostEqual(grid.ggf_pb[index_500], 5.0558)
        self.assertAlmostEqual(grid.vbf_pb[index_500], 0.54126)

        interpolated = self.plotter.interpolate_yr4_cross_sections(
            np.array([5.0, 125.0, 500.0, 3001.0]),
            grid,
        )
        self.assertTrue(np.isnan(interpolated["ggf_pb"][0]))
        self.assertAlmostEqual(interpolated["ggf_pb"][1], 45.142)
        self.assertAlmostEqual(interpolated["vbf_pb"][2], 0.54126)
        self.assertTrue(np.isnan(interpolated["vbf_pb"][3]))

        midpoint = self.plotter.interpolate_no_extrapolation(
            np.array([15.0]),
            np.array([10.0, 20.0]),
            np.array([100.0, 25.0]),
            log_y=True,
        )
        self.assertAlmostEqual(midpoint[0], 50.0)

    def test_signal_observables_use_k2_squared_and_strict_thresholds(self):
        grid = self.plotter.load_yr4_cross_section_grid()
        result = self.plotter.derive_signal_observables(
            np.array([500.0, 500.0, 200.0, 5.0]),
            np.array([100.0, 100.0, 100.0, 1.0]),
            np.array([0.2, -0.2, 0.2, 0.2]),
            np.array([5.0, 5.0, 1.0, 1.0]),
            np.array([0.5, 0.5, 0.5, 0.5]),
            np.array([True, True, True, True]),
            grid,
        )
        expected_ggf = 1000.0 * 5.0558 * 0.2**2 * 0.5
        expected_vbf = 1000.0 * 0.54126 * 0.2**2 * 0.5
        self.assertAlmostEqual(result["signal_ggf_rate_fb"][0], expected_ggf)
        self.assertAlmostEqual(result["signal_vbf_rate_fb"][0], expected_vbf)
        self.assertAlmostEqual(
            result["signal_dominant_rate_fb"][0],
            expected_ggf + expected_vbf,
        )
        self.assertAlmostEqual(
            result["signal_dominant_rate_fb"][0],
            result["signal_dominant_rate_fb"][1],
        )
        self.assertEqual(
            result["signal_width_category"].tolist(),
            ["intermediate", "intermediate", "narrow", "broad"],
        )
        self.assertTrue(result["signal_viable_open"][0])
        self.assertTrue(result["signal_viable_open"][1])
        self.assertFalse(result["h2_h3h3_kinematically_open"][2])
        self.assertFalse(result["signal_viable_open"][2])
        self.assertFalse(result["signal_yr4_grid_available"][3])
        self.assertFalse(result["signal_viable_open"][3])
        self.assertLess(
            result["signal_ggf_rate_low_fb"][0],
            result["signal_ggf_rate_fb"][0],
        )
        self.assertGreater(
            result["signal_ggf_rate_high_fb"][0],
            result["signal_ggf_rate_fb"][0],
        )

    def test_signal_fixture_derives_full_viable_rates_and_width_categories(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_signal_fixture(tmpdir)

        self.assertTrue(self.plotter.has_signal_results(data))
        self.assertIsNone(self.plotter.signal_availability_reason(data))
        self.assertEqual(np.count_nonzero(data.b("signal_viable_open")), 3)
        self.assertEqual(
            data.derived["signal_width_category"].tolist(),
            ["narrow", "intermediate", "broad"],
        )
        self.assertTrue(np.all(data.f("signal_dominant_rate_fb") > 0.0))
        np.testing.assert_allclose(
            data.f("signal_raw_hllhc_events"),
            data.f("signal_dominant_rate_fb")
            * self.plotter.HL_LHC_LUMINOSITY_FB,
        )

    def test_legacy_fixture_explains_signal_plot_omission(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_fixture(tmpdir)

        self.assertFalse(self.plotter.has_signal_results(data))
        reason = self.plotter.signal_availability_reason(data)
        self.assertIn("k2", reason)
        self.assertIn("h2_h3h3_br", reason)

    def test_fixture_ratios_match_expected_values(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_fixture(tmpdir)

        np.testing.assert_allclose(data.f("relic_ratio"), [2, 0.5, 3, 0.2, 0.8, 0.1, 4, 1])
        np.testing.assert_allclose(data.f("direct_ratio"), [2, 0.5, 0.25, 0.1, 4, 0.2, 3, 1])
        expected_indirect = [np.nan, np.nan, 0.25, 0.75, 1.25, np.nan, 0.5, 1.0]
        np.testing.assert_allclose(data.f("indirect_ratio_available"), expected_indirect, equal_nan=True)

    def test_finite_mask_filters_each_required_array(self):
        first = np.array([1.0, np.nan, np.inf, -2.0])
        second = np.array([3.0, 4.0, 5.0, np.nan])
        np.testing.assert_array_equal(
            self.plotter.finite_mask(first, second),
            [True, False, False, False],
        )

    def test_robust_norms_handle_one_sided_constant_and_signed_values(self):
        threshold = self.plotter.robust_threshold_norm(
            np.array([-2.0, -0.5, 0.0, 0.5, 3.0, np.nan, np.inf]),
            center=0.0,
        )
        self.assertLess(threshold.vmin, 0.0)
        self.assertEqual(threshold.vcenter, 0.0)
        self.assertGreater(threshold.vmax, 0.0)

        for values in [np.array([0.2, 0.4, np.nan]), np.array([0.0, 0.0])]:
            with self.subTest(values=values):
                norm = self.plotter.robust_threshold_norm(values, center=0.0)
                self.assertLess(norm.vmin, 0.0)
                self.assertGreater(norm.vmax, 0.0)

        signed = self.plotter.robust_symlog_norm(
            np.array([-100.0, -1.0, 0.0, 1.0, 100.0, np.nan, np.inf])
        )
        self.assertLess(signed.vmin, 0.0)
        self.assertGreater(signed.vmax, 0.0)
        self.assertGreater(signed.linthresh, 0.0)
        mapped = signed(np.array([-1.0, 0.0, 1.0]))
        self.assertLess(mapped[0], mapped[1])
        self.assertLess(mapped[1], mapped[2])

        logarithmic = self.plotter.robust_log_norm_to_one(
            np.array([1.0e-6, 1.0e-3, 0.2, np.nan, np.inf])
        )
        self.assertGreater(logarithmic.vmin, 0.0)
        self.assertGreaterEqual(logarithmic.vmax, 1.0)
        self.assertLess(logarithmic(1.0e-3), logarithmic(1.0))

        ticks = self.plotter.sparse_symlog_ticks(signed)
        self.assertIn(0.0, ticks)
        self.assertLessEqual(len([tick for tick in ticks if tick > 0.0]), 3)
        self.assertEqual(ticks, sorted(ticks))

        bsmpt_spec = self.plotter.PLOT_BY_STEM[
            "31_bsmpt_ew_true_over_t_m2_m3"
        ]
        bsmpt_norm = self.plotter.norm_for(
            bsmpt_spec, np.array([0.4, 0.9, 1.0, 1.5, 2.0])
        )
        self.assertEqual(
            bsmpt_norm.vcenter,
            self.plotter.BSMPT_STRONG_EWPT_THRESHOLD,
        )
        jump_spec = self.plotter.PLOT_BY_STEM[
            "31b_bsmpt_ew_jump_over_t_m2_m3"
        ]
        jump_norm = self.plotter.norm_for(
            jump_spec, np.array([0.2, 0.8, 1.0, 1.4])
        )
        self.assertEqual(
            jump_norm.vcenter,
            self.plotter.BSMPT_STRONG_EWPT_THRESHOLD,
        )

    def test_continuous_marker_groups_share_one_norm(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_fixture(tmpdir)
        spec = self.plotter.PLOT_BY_STEM["12_log10_relic_ratio_m2_m3"]
        fig, ax = self.plotter.plt.subplots()
        try:
            self.plotter.render_continuous_mass(fig, ax, data, spec)
            scatter_norms = [collection.norm for collection in ax.collections]
            self.assertEqual(len(scatter_norms), 4)
            self.assertTrue(all(norm is scatter_norms[0] for norm in scatter_norms))
        finally:
            self.plotter.plt.close(fig)

    def test_h2_width_map_uses_dm_markers_and_logarithmic_color(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_fixture(tmpdir)
        spec = self.plotter.PLOT_BY_STEM["23_h2_width_dm_status_m2_m3"]
        fig, ax = self.plotter.plt.subplots()
        try:
            self.plotter.render_continuous_mass(fig, ax, data, spec)
            scatter_norms = [collection.norm for collection in ax.collections]
            self.assertEqual(len(scatter_norms), 2)
            self.assertTrue(all(norm is scatter_norms[0] for norm in scatter_norms))
            self.assertIsInstance(scatter_norms[0], self.plotter.LogNorm)
            legend_labels = [text.get_text() for text in ax.get_legend().get_texts()]
            self.assertTrue(any(label.startswith("Fail:") for label in legend_labels))
            self.assertTrue(any(label.startswith("Pass:") for label in legend_labels))
            self.assertFalse(
                any(text.get_text().startswith("finite N") for text in ax.texts)
            )
        finally:
            self.plotter.plt.close(fig)

    def test_h2_width_map_is_skipped_for_legacy_input_without_w2(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "legacy.tsv"
            write_fixture(
                path,
                header=HEADER[:-1],
                rows=[row[:-1] for row in ROWS],
            )
            data = self.plotter.load_scan(path)
        spec = self.plotter.PLOT_BY_STEM["23_h2_width_dm_status_m2_m3"]
        fig, ax = self.plotter.plt.subplots()
        try:
            with self.assertRaisesRegex(
                self.plotter.PlotUnavailable, "w2 column is unavailable"
            ):
                self.plotter.render_continuous_mass(fig, ax, data, spec)
        finally:
            self.plotter.plt.close(fig)

    def test_resonance_plots_apply_requested_selections_and_center_zero(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_fixture(tmpdir)

        cases = (
            (
                "38_k133_vs_m3_experimental_resonance",
                "experimental",
                4,
            ),
            (
                "39_k133_vs_m3_relic_pass_resonance",
                "relic_pass",
                5,
            ),
            (
                "40_k233_vs_m3_all_resonance",
                "all",
                8,
            ),
        )
        norm_bounds = []
        axis_bounds = []
        for stem, selection, expected_count in cases:
            with self.subTest(stem=stem):
                spec = self.plotter.PLOT_BY_STEM[stem]
                fig, ax = self.plotter.plt.subplots()
                try:
                    self.plotter.render_resonance_xy(fig, ax, data, spec)
                    self.assertEqual(len(ax.collections), 1)
                    points = ax.collections[0]
                    self.assertEqual(len(points.get_offsets()), expected_count)
                    self.assertIsInstance(points.norm, self.plotter.SymLogNorm)
                    self.assertLess(points.norm.vmin, 0.0)
                    self.assertGreater(points.norm.vmax, 0.0)
                    self.assertLess(
                        sum(points.cmap(points.norm(0.0))[:3]), 0.5
                    )
                    self.assertEqual(ax.get_yscale(), "symlog")
                    norm_bounds.append((points.norm.vmin, points.norm.vmax))
                    axis_bounds.append((ax.get_xlim(), ax.get_ylim()))
                    mask = self.plotter.selection_mask(data, selection)
                    np.testing.assert_allclose(
                        np.sort(np.asarray(points.get_array(), dtype=float)),
                        np.sort(data.f("m2_minus_2m3")[mask]),
                    )
                    self.assertIn(f"N={expected_count}", ax.get_title())
                finally:
                    self.plotter.plt.close(fig)
        self.assertTrue(all(bounds == norm_bounds[0] for bounds in norm_bounds))
        self.assertTrue(
            all(
                np.allclose(xlimits, axis_bounds[0][0])
                and np.allclose(ylimits, axis_bounds[0][1])
                for xlimits, ylimits in axis_bounds
            )
        )

    def test_portal_mass_maps_color_only_selected_points(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_fixture(tmpdir)

        cases = (
            ("43_k133_experimental_m2_m3", "experimental", "K133", 4),
            ("44_k133_dm_m2_m3", "dm", "K133", 4),
            ("45_k233_experimental_m2_m3", "experimental", "K233", 4),
            ("46_k233_dm_m2_m3", "dm", "K233", 4),
        )
        norm_bounds = {}
        for stem, selection, value_name, expected_count in cases:
            with self.subTest(stem=stem):
                spec = self.plotter.PLOT_BY_STEM[stem]
                fig, ax = self.plotter.plt.subplots()
                try:
                    self.plotter.render_selected_continuous_mass(
                        fig, ax, data, spec
                    )
                    self.assertEqual(len(ax.collections), 2)
                    background, colored = ax.collections
                    self.assertEqual(len(background.get_offsets()), len(data))
                    self.assertEqual(len(colored.get_offsets()), expected_count)
                    self.assertIsInstance(colored.norm, self.plotter.SymLogNorm)
                    norm_bounds.setdefault(value_name, []).append(
                        (colored.norm.vmin, colored.norm.vmax)
                    )
                    mask = self.plotter.selection_mask(data, selection)
                    np.testing.assert_allclose(
                        np.sort(np.asarray(colored.get_array(), dtype=float)),
                        np.sort(data.f(value_name)[mask]),
                    )
                    self.assertIn("gray: all stored points", ax.get_title())
                finally:
                    self.plotter.plt.close(fig)
        for bounds in norm_bounds.values():
            self.assertEqual(bounds[0], bounds[1])

    def test_bsmpt_status_strength_and_phase_renderers(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_bsmpt_fixture(tmpdir)

        status_spec = self.plotter.PLOT_BY_STEM["30_bsmpt_status_m2_m3"]
        fig, ax = self.plotter.plt.subplots()
        try:
            self.plotter.render_categorical_mass(ax, data, status_spec)
            labels = [text.get_text() for text in ax.get_legend().get_texts()]
            self.assertEqual(len(ax.collections), 5)
            self.assertTrue(any(label.startswith("Not run: 2") for label in labels))
            self.assertTrue(
                any(label.startswith("BSMPT failed: 1") for label in labels)
            )
            self.assertIn(r"$v_{\rm EW,true}(T_*)/T_*\geq1$", ax.get_title())
        finally:
            self.plotter.plt.close(fig)

        strength_spec = self.plotter.PLOT_BY_STEM[
            "31_bsmpt_ew_true_over_t_m2_m3"
        ]
        fig, ax = self.plotter.plt.subplots()
        try:
            self.plotter.render_continuous_mass(fig, ax, data, strength_spec)
            scatter_norms = [collection.norm for collection in ax.collections]
            self.assertEqual(len(scatter_norms), 4)
            self.assertTrue(all(norm is scatter_norms[0] for norm in scatter_norms))
            self.assertEqual(
                scatter_norms[0].vcenter,
                self.plotter.BSMPT_STRONG_EWPT_THRESHOLD,
            )
            labels = [text.get_text() for text in ax.get_legend().get_texts()]
            self.assertEqual(len(labels), 4)
            self.assertFalse(any(label.startswith("Not run:") for label in labels))
        finally:
            self.plotter.plt.close(fig)

        jump_spec = self.plotter.PLOT_BY_STEM[
            "31b_bsmpt_ew_jump_over_t_m2_m3"
        ]
        fig, ax = self.plotter.plt.subplots()
        try:
            self.plotter.render_continuous_mass(fig, ax, data, jump_spec)
            self.assertEqual(len(ax.collections), 4)
            self.assertIn(
                r"$\Delta v_{\rm EW}(T_*)/T_*$",
                [axis.get_ylabel() for axis in fig.axes],
            )
        finally:
            self.plotter.plt.close(fig)

    def test_bsmpt_strength_xy_and_count_renderers(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_bsmpt_fixture(tmpdir)

        strength_spec = self.plotter.PLOT_BY_STEM["34_bsmpt_strength_vs_m2"]
        fig, ax = self.plotter.plt.subplots()
        try:
            self.plotter.render_bsmpt_strength_xy(ax, data, strength_spec)
            self.assertEqual(len(ax.collections), 4)
            self.assertEqual(ax.get_yscale(), "log")
            self.assertTrue(
                any(
                    np.allclose(line.get_ydata(), [1.0, 1.0])
                    for line in ax.lines
                )
            )
        finally:
            self.plotter.plt.close(fig)

        count_spec = self.plotter.PLOT_BY_STEM["36_bsmpt_counts"]
        fig, ax = self.plotter.plt.subplots()
        try:
            self.plotter.render_bsmpt_bars(ax, data, count_spec)
            self.assertEqual(len(ax.patches), 8)
            widths = [patch.get_width() for patch in ax.patches]
            self.assertEqual(widths[:3], [6, 5, 1])
        finally:
            self.plotter.plt.close(fig)

    def test_bsmpt_plots_are_unavailable_without_attempts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_fixture(tmpdir)
        spec = self.plotter.PLOT_BY_STEM["30_bsmpt_status_m2_m3"]
        fig, ax = self.plotter.plt.subplots()
        try:
            with self.assertRaisesRegex(
                self.plotter.PlotUnavailable,
                "BSMPT was not run",
            ):
                self.plotter.render_spec(fig, ax, data, spec)
        finally:
            self.plotter.plt.close(fig)

    def test_ew_jump_plot_is_omitted_for_legacy_bsmpt_scan(self):
        legacy_header = [
            column
            for column in BSMPT_HEADER
            if column != "ewpt_ew_jump_over_T"
        ]
        jump_index = BSMPT_HEADER.index("ewpt_ew_jump_over_T")
        legacy_rows = [
            row[:jump_index] + row[jump_index + 1 :]
            for row in bsmpt_fixture_rows()
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "legacy_bsmpt.tsv"
            write_fixture(path, legacy_header, legacy_rows)
            data = self.plotter.load_scan(path)

        stems = self.plotter.figure_stems_for_data(data)
        self.assertIn("31_bsmpt_ew_true_over_t_m2_m3", stems)
        self.assertNotIn("31b_bsmpt_ew_jump_over_t_m2_m3", stems)
        self.assertNotIn("31d_bsmpt_ew_entry_status_m2_m3", stems)
        self.assertNotIn("34b_bsmpt_ew_entry_strength_vs_m2", stems)
        self.assertNotIn("35b_bsmpt_ew_entry_strength_vs_m3", stems)
        self.assertNotIn("35c_bsmpt_selected_vs_ew_entry_jump", stems)
        self.assertNotIn("31f_bsmpt_gw_status_m2_m3", stems)

    def test_ew_entry_status_plot_shows_runs_without_critical_entry(self):
        extra = [
            "ewpt_ew_entry_jump_over_T",
            "ewpt_baryo_candidate",
        ]
        rows = [row + ["nan"] * len(extra) for row in bsmpt_fixture_rows()]
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "no_completed_entry.tsv"
            write_fixture(path, BSMPT_HEADER + extra, rows)
            data = self.plotter.load_scan(path)
        stems = self.plotter.figure_stems_for_data(data)
        self.assertIn("31d_bsmpt_ew_entry_status_m2_m3", stems)
        self.assertNotIn("31c_bsmpt_ew_entry_jump_over_t_m2_m3", stems)
        self.assertNotIn("34b_bsmpt_ew_entry_strength_vs_m2", stems)
        self.assertNotIn("35b_bsmpt_ew_entry_strength_vs_m3", stems)
        self.assertNotIn("35c_bsmpt_selected_vs_ew_entry_jump", stems)

    def test_ew_entry_plot_and_summary_use_the_new_criterion(self):
        extra = [
            "ewpt_ew_entry_jump_over_T",
            "ewpt_ew_entry_percolated",
            "ewpt_ew_entry_completed",
            "ewpt_baryo_candidate",
            "ewpt_ew_entry_nucl_jump_over_T",
            "ewpt_ew_entry_perc_jump_over_T",
            "ewpt_gw_crit_field_jump_over_T",
            "ewpt_gw_nucl_field_jump_over_T",
            "ewpt_gw_perc_field_jump_over_T",
            "ewpt_gw_max_field_jump_over_T",
            "ewpt_gw_candidate",
        ]
        rows = [row + ["nan"] * len(extra) for row in bsmpt_fixture_rows()]
        def set_values(index, **values):
            for column, value in values.items():
                rows[index][len(BSMPT_HEADER) + extra.index(column)] = value

        set_values(2, ewpt_gw_crit_field_jump_over_T=1.2,
                   ewpt_gw_max_field_jump_over_T=1.2,
                   ewpt_gw_candidate=True)
        set_values(3, ewpt_ew_entry_jump_over_T=0.2,
                   ewpt_ew_entry_percolated=True, ewpt_ew_entry_completed=True,
                   ewpt_baryo_candidate=False,
                   ewpt_ew_entry_nucl_jump_over_T=0.3,
                   ewpt_ew_entry_perc_jump_over_T=0.4,
                   ewpt_gw_crit_field_jump_over_T=0.2,
                   ewpt_gw_nucl_field_jump_over_T=0.3,
                   ewpt_gw_perc_field_jump_over_T=0.4,
                   ewpt_gw_max_field_jump_over_T=0.4, ewpt_gw_candidate=False)
        set_values(4, ewpt_ew_entry_jump_over_T=1.2,
                   ewpt_ew_entry_percolated=False, ewpt_ew_entry_completed=False,
                   ewpt_baryo_candidate=True,
                   ewpt_gw_crit_field_jump_over_T=1.2,
                   ewpt_gw_max_field_jump_over_T=1.2, ewpt_gw_candidate=True)
        set_values(5, ewpt_ew_entry_jump_over_T=0.5,
                   ewpt_ew_entry_percolated=True, ewpt_ew_entry_completed=False,
                   ewpt_baryo_candidate=False,
                   ewpt_ew_entry_nucl_jump_over_T=1.5,
                   ewpt_ew_entry_perc_jump_over_T=1.8,
                   ewpt_gw_crit_field_jump_over_T=0.5,
                   ewpt_gw_nucl_field_jump_over_T=1.5,
                   ewpt_gw_perc_field_jump_over_T=1.8,
                   ewpt_gw_max_field_jump_over_T=1.8, ewpt_gw_candidate=True)
        set_values(6, ewpt_ew_entry_jump_over_T=0.3,
                   ewpt_ew_entry_percolated=False, ewpt_ew_entry_completed=False,
                   ewpt_baryo_candidate=False,
                   ewpt_gw_crit_field_jump_over_T=0.3,
                   ewpt_gw_max_field_jump_over_T=0.3, ewpt_gw_candidate=False)
        rows[3][HEADER.index("ewpt_ew_true_over_T")] = 2.5
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "ew_entry.tsv"
            write_fixture(path, BSMPT_HEADER + extra, rows)
            data = self.plotter.load_scan(path)

        stems = self.plotter.figure_stems_for_data(data)
        self.assertIn("31c_bsmpt_ew_entry_jump_over_t_m2_m3", stems)
        self.assertIn("31d_bsmpt_ew_entry_status_m2_m3", stems)
        self.assertIn("31e_bsmpt_gw_max_jump_over_t_m2_m3", stems)
        self.assertIn("31f_bsmpt_gw_status_m2_m3", stems)
        self.assertIn("34b_bsmpt_ew_entry_strength_vs_m2", stems)
        self.assertIn("35b_bsmpt_ew_entry_strength_vs_m3", stems)
        self.assertIn("35c_bsmpt_selected_vs_ew_entry_jump", stems)
        self.assertIn("35d_bsmpt_gw_max_jump_vs_m3", stems)
        self.assertIn("35e_bsmpt_ew_entry_temperature_jumps", stems)
        self.assertIn("35f_bsmpt_gw_temperature_jumps", stems)
        self.assertEqual(
            data.derived["bsmpt_ew_entry"].tolist(),
            [
                "not run",
                "failed",
                "no recorded EW entry",
                "EW entry weak",
                "baryogenesis candidate",
                "EW entry weak",
                "EW entry weak",
                "not run",
            ],
        )
        self.assertEqual(
            data.derived["bsmpt_gw"].tolist(),
            ["not run", "failed", "GW candidate", "GW weak", "GW candidate",
             "GW candidate", "GW weak", "not run"],
        )
        summary = {row.metric: row for row in self.plotter.build_summary(data)}
        self.assertEqual(summary["bsmpt_ew_entry_fopt_identified"].count, 4)
        self.assertEqual(summary["bsmpt_ew_entry_percolated"].count, 2)
        self.assertEqual(summary["bsmpt_ew_entry_completed"].count, 1)
        self.assertEqual(summary["bsmpt_baryogenesis_candidates"].count, 1)
        self.assertEqual(summary["bsmpt_baryogenesis_candidates"].denominator, 4)
        self.assertEqual(summary["bsmpt_gw_candidates"].count, 3)
        self.assertEqual(summary["bsmpt_gw_candidates"].denominator, 5)
        self.assertEqual(len(self.plotter.bsmpt_bar_metrics(data)), 10)
        spec = self.plotter.PLOT_BY_STEM[
            "31c_bsmpt_ew_entry_jump_over_t_m2_m3"
        ]
        fig, ax = self.plotter.plt.subplots()
        try:
            self.plotter.render_continuous_mass(fig, ax, data, spec)
            self.assertTrue(ax.collections)
        finally:
            self.plotter.plt.close(fig)
        for stem in (
            "34b_bsmpt_ew_entry_strength_vs_m2",
            "35b_bsmpt_ew_entry_strength_vs_m3",
        ):
            mass_spec = self.plotter.PLOT_BY_STEM[stem]
            fig, ax = self.plotter.plt.subplots()
            try:
                self.plotter.render_bsmpt_strength_xy(ax, data, mass_spec)
                self.assertEqual(len(ax.collections), 2)
                self.assertEqual(ax.get_yscale(), "log")
                self.assertIn("EW entry at", ax.get_title())
            finally:
                self.plotter.plt.close(fig)
        comparison = self.plotter.PLOT_BY_STEM[
            "35c_bsmpt_selected_vs_ew_entry_jump"
        ]
        fig, ax = self.plotter.plt.subplots()
        try:
            self.plotter.render_spec(fig, ax, data, comparison)
            self.assertEqual(ax.get_xscale(), "log")
            self.assertEqual(ax.get_yscale(), "log")
            self.assertTrue(any("EW jump" in note.get_text() for note in ax.texts))
        finally:
            self.plotter.plt.close(fig)
        status_spec = self.plotter.PLOT_BY_STEM[
            "31d_bsmpt_ew_entry_status_m2_m3"
        ]
        fig, ax = self.plotter.plt.subplots()
        try:
            self.plotter.render_categorical_mass(ax, data, status_spec)
            self.assertEqual(len(ax.collections), 5)
            labels = [text.get_text() for text in ax.get_legend().get_texts()]
            self.assertTrue(any("EW entry" in label for label in labels))
        finally:
            self.plotter.plt.close(fig)

        for stem in ("35e_bsmpt_ew_entry_temperature_jumps",
                     "35f_bsmpt_gw_temperature_jumps"):
            fig, ax = self.plotter.plt.subplots()
            try:
                self.plotter.render_spec(fig, ax, data, self.plotter.PLOT_BY_STEM[stem])
                self.assertEqual(ax.get_xscale(), "log")
                self.assertEqual(ax.get_yscale(), "log")
                self.assertGreaterEqual(len(ax.collections), 2)
            finally:
                self.plotter.plt.close(fig)

    def test_freezeout_x_window_plot_uses_sampled_global_window(self):
        extra = [
            "dm_freezeout_temperature_GeV",
            "ewpt_x_broken_min_T_GeV",
            "ewpt_x_broken_max_T_GeV",
            "ewpt_x_broken_at_or_after_freezeout",
        ]
        rows = [row + ["nan"] * len(extra) for row in bsmpt_fixture_rows()]
        rows[3][-4:] = [10, 30, 70, False]
        rows[4][-4:] = [40, 30, 70, True]
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "freezeout.tsv"
            write_fixture(path, BSMPT_HEADER + extra, rows)
            data = self.plotter.load_scan(path)

        spec = self.plotter.PLOT_BY_STEM["35g_freezeout_vs_x_broken_window"]
        self.assertIsNone(self.plotter.spec_unavailable_reason(data, spec))
        summary = {row.metric: row for row in self.plotter.build_summary(data)}
        self.assertEqual(summary["dm_freezeout_temperature_available"].count, 2)
        self.assertEqual(summary["bsmpt_x_freezeout_comparison_assessed"].count, 2)
        self.assertEqual(summary["bsmpt_x_broken_at_or_after_freezeout"].count, 1)
        fig, ax = self.plotter.plt.subplots()
        try:
            self.plotter.render_spec(fig, ax, data, spec)
            self.assertEqual(ax.get_xscale(), "symlog")
            self.assertEqual(ax.get_yscale(), "symlog")
            self.assertEqual(len(ax.collections), 4)
            self.assertEqual(sum(len(item.get_offsets()) for item in ax.collections), 4)
        finally:
            self.plotter.plt.close(fig)

    def test_cumulative_renderer_uses_nested_styles_and_skips_missing_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_fixture(tmpdir)

        spec = self.plotter.PLOT_BY_STEM[
            "26_cumulative_constraints_m2_a12"
        ]
        fig, ax = self.plotter.plt.subplots()
        try:
            self.plotter.render_cumulative_xy(ax, data, spec)
            self.assertEqual(len(ax.collections), 5)
            labels = [text.get_text() for text in ax.get_legend().get_texts()]
            self.assertEqual(len(labels), 5)
            self.assertTrue(labels[0].startswith("All stored: 8 (100.00%)"))
            self.assertTrue(labels[-1].endswith("2 (25.00%)"))
            self.assertIn("EWPO not applied", ax.get_title())
        finally:
            self.plotter.plt.close(fig)

        missing_spec = self.plotter.PLOT_BY_STEM[
            "25_cumulative_constraints_m2_vs"
        ]
        fig, ax = self.plotter.plt.subplots()
        try:
            with self.assertRaisesRegex(
                self.plotter.PlotUnavailable, "vs column is unavailable"
            ):
                self.plotter.render_cumulative_xy(ax, data, missing_spec)
        finally:
            self.plotter.plt.close(fig)

    def test_mass_guides_include_both_requested_relations(self):
        guide_rows = [list(ROWS[index]) for index in range(3)]
        for row, (m2, m3) in zip(
            guide_rows,
            ((100.0, 100.0), (1000.0, 100.0), (100.0, 1000.0)),
        ):
            row[HEADER.index("M2")] = m2
            row[HEADER.index("M3")] = m3
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "guides.tsv"
            write_fixture(path, rows=guide_rows)
            data = self.plotter.load_scan(path)
        fig, ax = self.plotter.plt.subplots()
        try:
            self.plotter.draw_mass_guides(ax, data, annotate=True)
            lines_by_color = {line.get_color(): line for line in ax.lines}

            m2_eq_2m3 = lines_by_color[
                self.plotter.M2_EQ_2M3_GUIDE_COLOR
            ]
            np.testing.assert_allclose(
                m2_eq_2m3.get_ydata(),
                0.5 * m2_eq_2m3.get_xdata(),
            )
            self.assertEqual(m2_eq_2m3.get_linestyle(), "-.")

            m3_eq_2m2 = lines_by_color[
                self.plotter.M3_EQ_2M2_GUIDE_COLOR
            ]
            np.testing.assert_allclose(
                m3_eq_2m2.get_ydata(),
                2.0 * m3_eq_2m2.get_xdata(),
            )
            self.assertEqual(m3_eq_2m2.get_linestyle(), "--")

            annotation = "\n".join(text.get_text() for text in ax.texts)
            self.assertIn(r"$M_2=2M_3$", annotation)
            self.assertIn(r"$M_3=2M_2$", annotation)
        finally:
            self.plotter.plt.close(fig)

    def test_registry_and_expected_paths_are_unique(self):
        self.assertEqual(len(self.plotter.PLOT_SPECS), 81)
        self.assertEqual(len(self.plotter.DASHBOARDS), 12)
        stems = self.plotter.all_figure_stems()
        self.assertEqual(len(stems), 93)
        self.assertEqual(len(set(stems)), 93)
        self.assertTrue(any("bsmpt" in stem for stem in stems))
        self.assertTrue(any("signal" in stem for stem in stems))
        self.assertEqual(
            [spec.stem for spec in self.plotter.PLOT_SPECS],
            [
                "01_dm_experimental_fourway_m2_m3",
                "02_dm_status_m2_m3",
                "03_experimental_status_m2_m3",
                "04_full_viability_status_m2_m3",
                "05_higgsbounds_status_m2_m3",
                "06_higgssignals_status_m2_m3",
                "07_wmass_status_m2_m3",
                "08_relic_density_status_m2_m3",
                "09_direct_detection_status_m2_m3",
                "10_indirect_detection_status_m2_m3",
                "11_dm_failure_modes_m2_m3",
                "12_log10_relic_ratio_m2_m3",
                "13_log10_direct_detection_ratio_m2_m3",
                "14_indirect_ratio_m2_m3",
                "15_abs_a12_m2_m3",
                "16_k233_m2_m3",
                "17_lsx_m2_m3",
                "18_higgssignals_delta_chi2_m2_m3",
                "19_m2_vs_abs_a12",
                "20_k133_vs_k233",
                "21_relic_vs_direct_ratio",
                "22_constraint_counts",
                "23_h2_width_dm_status_m2_m3",
                "24_cumulative_constraints_m2_m3",
                "25_cumulative_constraints_m2_vs",
                "26_cumulative_constraints_m2_a12",
                "27_cumulative_constraints_m3_lx",
                "28_cumulative_constraints_m3_lphix",
                "29_cumulative_constraints_m3_lsx",
                "30_bsmpt_status_m2_m3",
                "31_bsmpt_ew_true_over_t_m2_m3",
                "31b_bsmpt_ew_jump_over_t_m2_m3",
                "31c_bsmpt_ew_entry_jump_over_t_m2_m3",
                "31d_bsmpt_ew_entry_status_m2_m3",
                "31e_bsmpt_gw_max_jump_over_t_m2_m3",
                "31f_bsmpt_gw_status_m2_m3",
                "32_bsmpt_phase_history_m2_m3",
                "33_bsmpt_ew_entry_step_m2_m3",
                "34_bsmpt_strength_vs_m2",
                "34b_bsmpt_ew_entry_strength_vs_m2",
                "35_bsmpt_strength_vs_m3",
                "35b_bsmpt_ew_entry_strength_vs_m3",
                "35c_bsmpt_selected_vs_ew_entry_jump",
                "35d_bsmpt_gw_max_jump_vs_m3",
                "35e_bsmpt_ew_entry_temperature_jumps",
                "35f_bsmpt_gw_temperature_jumps",
                "35g_freezeout_vs_x_broken_window",
                "36_bsmpt_counts",
                "37_k133_vs_m3_all_resonance",
                "38_k133_vs_m3_experimental_resonance",
                "39_k133_vs_m3_relic_pass_resonance",
                "40_k233_vs_m3_all_resonance",
                "41_k233_vs_m3_experimental_resonance",
                "42_k233_vs_m3_relic_pass_resonance",
                "43_k133_experimental_m2_m3",
                "44_k133_dm_m2_m3",
                "45_k233_experimental_m2_m3",
                "46_k233_dm_m2_m3",
                "47_m2_vs_a12",
                "48_h2_h1h1_one_h1_invisible_xsec_vs_m2",
                "49_h2_h1h1_one_h1_invisible_xsec_vs_m3",
                "50_h1_h2h2_one_h2_invisible_xsec_vs_m2",
                "51_h1_h2h2_one_h2_invisible_xsec_vs_m3",
                "52_mono_higgs_xsec_vs_m2",
                "53_mono_higgs_xsec_vs_m3",
                "54_mono_z_xsec_vs_m2",
                "55_mono_z_xsec_vs_m3",
                "56_h2_h1h1_one_h1_invisible_xsec_no_dm_vs_m2",
                "57_h2_h1h1_one_h1_invisible_xsec_no_dm_vs_m3",
                "58_h1_h2h2_one_h2_invisible_xsec_no_dm_vs_m2",
                "59_h1_h2h2_one_h2_invisible_xsec_no_dm_vs_m3",
                "60_mono_higgs_xsec_no_dm_vs_m2",
                "61_mono_higgs_xsec_no_dm_vs_m3",
                "62_mono_z_xsec_no_dm_vs_m2",
                "63_mono_z_xsec_no_dm_vs_m3",
                "64_signal_rate_m2_m3",
                "65_signal_rates_vs_m2",
                "66_signal_rate_vs_m3",
                "67_signal_k2sq_vs_br",
                "68_signal_width_fraction_vs_rate",
                "69_signal_rate_vs_direct_detection",
            ],
        )
        self.assertEqual(
            self.plotter.DASHBOARDS,
            {
                "dashboard_status_summary": (
                    "01_dm_experimental_fourway_m2_m3",
                    "04_full_viability_status_m2_m3",
                    "03_experimental_status_m2_m3",
                    "05_higgsbounds_status_m2_m3",
                    "06_higgssignals_status_m2_m3",
                    "07_wmass_status_m2_m3",
                ),
                "dashboard_dm_summary": (
                    "02_dm_status_m2_m3",
                    "08_relic_density_status_m2_m3",
                    "09_direct_detection_status_m2_m3",
                    "10_indirect_detection_status_m2_m3",
                ),
                "dashboard_diagnostic_summary": (
                    "12_log10_relic_ratio_m2_m3",
                    "13_log10_direct_detection_ratio_m2_m3",
                    "15_abs_a12_m2_m3",
                    "16_k233_m2_m3",
                    "17_lsx_m2_m3",
                    "18_higgssignals_delta_chi2_m2_m3",
                ),
                "dashboard_cumulative_constraint_summary": (
                    "24_cumulative_constraints_m2_m3",
                    "25_cumulative_constraints_m2_vs",
                    "26_cumulative_constraints_m2_a12",
                    "27_cumulative_constraints_m3_lx",
                    "28_cumulative_constraints_m3_lphix",
                    "29_cumulative_constraints_m3_lsx",
                ),
                "dashboard_bsmpt_summary": (
                    "30_bsmpt_status_m2_m3",
                    "31c_bsmpt_ew_entry_jump_over_t_m2_m3",
                    "31d_bsmpt_ew_entry_status_m2_m3",
                    "31e_bsmpt_gw_max_jump_over_t_m2_m3",
                    "31f_bsmpt_gw_status_m2_m3",
                    "32_bsmpt_phase_history_m2_m3",
                ),
                "dashboard_portal_resonance_summary": (
                    "37_k133_vs_m3_all_resonance",
                    "38_k133_vs_m3_experimental_resonance",
                    "39_k133_vs_m3_relic_pass_resonance",
                    "40_k233_vs_m3_all_resonance",
                    "41_k233_vs_m3_experimental_resonance",
                    "42_k233_vs_m3_relic_pass_resonance",
                ),
                "dashboard_portal_mass_plane_summary": (
                    "43_k133_experimental_m2_m3",
                    "44_k133_dm_m2_m3",
                    "45_k233_experimental_m2_m3",
                    "46_k233_dm_m2_m3",
                ),
                "dashboard_scalar_cascade_rates": (
                    "48_h2_h1h1_one_h1_invisible_xsec_vs_m2",
                    "49_h2_h1h1_one_h1_invisible_xsec_vs_m3",
                    "50_h1_h2h2_one_h2_invisible_xsec_vs_m2",
                    "51_h1_h2h2_one_h2_invisible_xsec_vs_m3",
                ),
                "dashboard_mg5_mono_rates": (
                    "52_mono_higgs_xsec_vs_m2",
                    "53_mono_higgs_xsec_vs_m3",
                    "54_mono_z_xsec_vs_m2",
                    "55_mono_z_xsec_vs_m3",
                ),
                "dashboard_scalar_cascade_rates_no_dm": (
                    "56_h2_h1h1_one_h1_invisible_xsec_no_dm_vs_m2",
                    "57_h2_h1h1_one_h1_invisible_xsec_no_dm_vs_m3",
                    "58_h1_h2h2_one_h2_invisible_xsec_no_dm_vs_m2",
                    "59_h1_h2h2_one_h2_invisible_xsec_no_dm_vs_m3",
                ),
                "dashboard_mg5_mono_rates_no_dm": (
                    "60_mono_higgs_xsec_no_dm_vs_m2",
                    "61_mono_higgs_xsec_no_dm_vs_m3",
                    "62_mono_z_xsec_no_dm_vs_m2",
                    "63_mono_z_xsec_no_dm_vs_m3",
                ),
                "dashboard_signal_summary": (
                    "64_signal_rate_m2_m3",
                    "65_signal_rates_vs_m2",
                    "66_signal_rate_vs_m3",
                    "67_signal_k2sq_vs_br",
                    "68_signal_width_fraction_vs_rate",
                    "69_signal_rate_vs_direct_detection",
                ),
            },
        )
        indirect = self.plotter.PLOT_BY_STEM["14_indirect_ratio_m2_m3"]
        self.assertEqual(indirect.value, "indirect_ratio_available")
        self.assertEqual(indirect.norm_kind, "log_to_one")
        k233 = self.plotter.PLOT_BY_STEM["16_k233_m2_m3"]
        self.assertEqual(k233.value, "abs_K233")
        self.assertEqual(k233.norm_kind, "log_to_one")
        self.assertEqual(k233.cmap, "viridis")
        h2_width = self.plotter.PLOT_BY_STEM["23_h2_width_dm_status_m2_m3"]
        self.assertEqual(h2_width.scheme, "dm")
        self.assertEqual(h2_width.value, "w2")
        self.assertEqual(h2_width.norm_kind, "log")
        bsmpt_strength = self.plotter.PLOT_BY_STEM[
            "31_bsmpt_ew_true_over_t_m2_m3"
        ]
        self.assertTrue(bsmpt_strength.requires_bsmpt)
        self.assertEqual(bsmpt_strength.norm_kind, "threshold1")
        resonance = self.plotter.PLOT_BY_STEM[
            "38_k133_vs_m3_experimental_resonance"
        ]
        self.assertEqual(resonance.kind, "resonance_xy")
        self.assertEqual(resonance.selection, "experimental")
        self.assertEqual(resonance.value, "m2_minus_2m3")
        self.assertEqual(resonance.cmap, "trsm_resonance")
        portal_mass = self.plotter.PLOT_BY_STEM[
            "46_k233_dm_m2_m3"
        ]
        self.assertEqual(portal_mass.kind, "selected_continuous_mass")
        self.assertEqual(portal_mass.selection, "dm")
        self.assertEqual(portal_mass.value, "K233")
        signed_mixing = self.plotter.PLOT_BY_STEM["47_m2_vs_a12"]
        self.assertEqual(signed_mixing.x, "M2")
        self.assertEqual(signed_mixing.y, "a12")
        scalar_cascade = self.plotter.PLOT_BY_STEM[
            "48_h2_h1h1_one_h1_invisible_xsec_vs_m2"
        ]
        self.assertEqual(scalar_cascade.kind, "rate_xy")
        self.assertEqual(scalar_cascade.selection, "full_viability")
        mono_higgs = self.plotter.PLOT_BY_STEM["52_mono_higgs_xsec_vs_m2"]
        self.assertEqual(mono_higgs.kind, "rate_xy")
        self.assertEqual(mono_higgs.selection, "full_viability")
        self.assertEqual(mono_higgs.y, "mono_higgs_xsec_pb")
        mono_higgs_no_dm = self.plotter.PLOT_BY_STEM[
            "60_mono_higgs_xsec_no_dm_vs_m2"
        ]
        self.assertEqual(mono_higgs_no_dm.selection, "non_dm_viability")
        signal_rate = self.plotter.PLOT_BY_STEM["64_signal_rate_m2_m3"]
        self.assertTrue(signal_rate.requires_signal)
        self.assertEqual(signal_rate.value, "signal_dominant_rate_fb")

        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_fixture(tmpdir)
        np.testing.assert_allclose(data.f("abs_K233"), np.abs(data.f("K233")))
        np.testing.assert_allclose(
            data.f("m2_minus_2m3"), data.f("M2") - 2.0 * data.f("M3")
        )
        _categories, binary_styles = self.plotter.category_styles(data, "dm")
        self.assertNotEqual(binary_styles["fail"].marker, binary_styles["pass"].marker)

        paths = self.plotter.expected_figure_paths(Path("plots"), "both")
        self.assertEqual(len(paths), 186)
        self.assertEqual(len(set(paths)), 186)
        self.assertEqual(sum(path.suffix == ".png" for path in paths), 93)
        self.assertEqual(sum(path.suffix == ".pdf" for path in paths), 93)

        legacy_paths = self.plotter.expected_figure_paths(
            Path("plots"), "both", data=data
        )
        self.assertEqual(len(legacy_paths), 92)
        self.assertFalse(any("bsmpt" in path.stem for path in legacy_paths))

        with tempfile.TemporaryDirectory() as tmpdir:
            bsmpt_data = self.load_bsmpt_fixture(tmpdir)
        bsmpt_paths = self.plotter.expected_figure_paths(
            Path("plots"), "both", data=bsmpt_data
        )
        self.assertEqual(len(bsmpt_paths), 110)
        self.assertTrue(any(path.stem == "dashboard_bsmpt_summary" for path in bsmpt_paths))

        with tempfile.TemporaryDirectory() as tmpdir:
            signal_data = self.load_signal_fixture(tmpdir)
        signal_paths = self.plotter.expected_figure_paths(
            Path("plots"), "both", data=signal_data
        )
        self.assertEqual(len(signal_paths), 106)
        self.assertTrue(
            any(path.stem == "dashboard_signal_summary" for path in signal_paths)
        )

    def test_summary_contains_categories_omissions_and_high_mass_tail(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_fixture(tmpdir)
            rows = self.plotter.build_summary(data)
            by_metric = {row.metric: row for row in rows}

            self.assertEqual(by_metric["total_rows"].count, 8)
            self.assertEqual(by_metric["experimental"].count, 4)
            self.assertEqual(by_metric["non_dm_viability"].count, 4)
            self.assertEqual(by_metric["dm"].count, 4)
            self.assertEqual(by_metric["full_viability"].count, 2)
            self.assertEqual(by_metric["cumulative_selection_all"].count, 8)
            self.assertEqual(by_metric["cumulative_selection_hb"].count, 6)
            self.assertEqual(by_metric["cumulative_selection_hb_hs"].count, 5)
            self.assertEqual(
                by_metric["cumulative_selection_hb_hs_wmass"].count, 4
            )
            self.assertEqual(
                by_metric["cumulative_selection_hb_hs_wmass_dm"].count, 2
            )
            self.assertIn(
                "EWPO intentionally not applied",
                by_metric["cumulative_selection_hb_hs_wmass_dm"].note,
            )
            self.assertEqual(by_metric["fourway_neither"].count, 2)
            self.assertEqual(by_metric["fourway_DM_only"].count, 2)
            self.assertEqual(by_metric["fourway_experimental_only"].count, 2)
            self.assertEqual(by_metric["fourway_both"].count, 2)
            self.assertEqual(by_metric["dm_failure_pass"].count, 4)
            self.assertEqual(by_metric["dm_failure_relic_only"].count, 1)
            self.assertEqual(by_metric["dm_failure_direct_only"].count, 1)
            self.assertEqual(by_metric["dm_failure_relic_and_direct"].count, 2)
            self.assertEqual(by_metric["indirect_unavailable"].count, 3)
            self.assertEqual(by_metric["indirect_allowed"].count, 4)
            self.assertEqual(by_metric["indirect_excluded"].count, 1)
            self.assertEqual(by_metric["m3_above_nominal_max"].count, 1)
            self.assertIn("m2_in_reversed_sampling_domain", by_metric)
            expected_invisible_open = np.count_nonzero(
                (2.0 * data.f("M3") < self.plotter.SM_LIKE_HIGGS_MASS_GEV)
                | (2.0 * data.f("M3") < data.f("M2"))
            )
            self.assertEqual(
                by_metric["higgs_invisible_decay_open"].count,
                expected_invisible_open,
            )
            self.assertEqual(
                by_metric["higgs_invisible_decay_open_but_unmodelled"].count,
                expected_invisible_open,
            )
            self.assertEqual(by_metric["bsmpt_attempted"].count, 0)
            self.assertIn("plots omitted", by_metric["bsmpt_attempted"].note)
            self.assertEqual(by_metric["ewpt_finite"].count, 0)
            self.assertIn("omitted_plot_evo", by_metric)
            self.assertIn("omitted_plot_thc", by_metric)
            self.assertIn("omitted_plot_ewpo", by_metric)

            summary_path = Path(tmpdir) / "constraint_summary.tsv"
            self.plotter.write_summary(summary_path, rows)
            text = summary_path.read_text(encoding="ascii")
            self.assertTrue(text.startswith("metric\tcount\tdenominator\tpercent\tnote\n"))
            self.assertIn("full_viability\t2\t8\t25", text)

    def test_bsmpt_summary_counts_and_denominators(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_bsmpt_fixture(tmpdir)
            by_metric = {
                row.metric: row for row in self.plotter.build_summary(data)
            }

        self.assertEqual(
            (by_metric["bsmpt_attempted"].count, by_metric["bsmpt_attempted"].denominator),
            (6, 8),
        )
        self.assertEqual(
            (by_metric["bsmpt_success"].count, by_metric["bsmpt_success"].denominator),
            (5, 6),
        )
        self.assertEqual(
            (by_metric["bsmpt_failed"].count, by_metric["bsmpt_failed"].denominator),
            (1, 6),
        )
        self.assertEqual(by_metric["bsmpt_success_no_selected_fopt"].count, 1)
        self.assertEqual(
            (by_metric["bsmpt_selected_fopt"].count, by_metric["bsmpt_selected_fopt"].denominator),
            (4, 5),
        )
        self.assertEqual(
            (
                by_metric["bsmpt_selected_strong_fopt"].count,
                by_metric["bsmpt_selected_strong_fopt"].denominator,
            ),
            (2, 4),
        )
        self.assertEqual(by_metric["bsmpt_phase_history_available"].count, 5)
        self.assertEqual(by_metric["bsmpt_direct_ew_entry"].count, 2)
        self.assertEqual(by_metric["bsmpt_multistep_ew_entry"].count, 3)
        self.assertEqual(by_metric["bsmpt_x_broken_path"].count, 1)
        self.assertIn(
            "not an additional scan constraint",
            by_metric["bsmpt_selected_strong_fopt"].note,
        )

    def test_signal_summary_records_rates_widths_and_yr4_provenance(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_signal_fixture(tmpdir)
            by_metric = {
                row.metric: row for row in self.plotter.build_summary(data)
            }

        self.assertEqual(
            (
                by_metric["signal_full_viable_h2_to_h3h3_open"].count,
                by_metric["signal_full_viable_h2_to_h3h3_open"].denominator,
            ),
            (3, 3),
        )
        self.assertEqual(
            by_metric["signal_full_viable_rate_available"].count,
            3,
        )
        self.assertEqual(by_metric["signal_width_narrow"].count, 1)
        self.assertEqual(by_metric["signal_width_intermediate"].count, 1)
        self.assertEqual(by_metric["signal_width_broad"].count, 1)
        self.assertEqual(by_metric["yr4_signal_table_rows"].count, 114)
        self.assertIn(
            self.plotter.YR4_SIGNAL_REPOSITORY_COMMIT,
            by_metric["yr4_signal_table_rows"].note,
        )
        self.assertIn(
            "range",
            by_metric["signal_full_viable_rate_available"].note,
        )

    def test_mass_guide_summary_uses_strict_regions(self):
        h2_open = list(ROWS[0])
        h2_open[HEADER.index("M2")] = 500.0
        h2_open[HEADER.index("M3")] = 100.0

        m3_above_2m2 = list(ROWS[1])
        m3_above_2m2[HEADER.index("M2")] = 100.0
        m3_above_2m2[HEADER.index("M3")] = 500.0

        at_threshold = list(ROWS[2])
        at_threshold[HEADER.index("M2")] = 200.0
        at_threshold[HEADER.index("M3")] = 400.0

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "thresholds.tsv"
            write_fixture(path, rows=[h2_open, m3_above_2m2, at_threshold])
            data = self.plotter.load_scan(path)
            by_metric = {
                row.metric: row for row in self.plotter.build_summary(data)
            }

        h2_metric = by_metric["h2_to_h3h3_kinematically_open"]
        h3_metric = by_metric["m3_above_2m2_reference"]
        self.assertEqual((h2_metric.count, h2_metric.denominator), (1, 3))
        self.assertEqual((h3_metric.count, h3_metric.denominator), (1, 3))
        self.assertIn("M2 = 2*M3", h2_metric.note)
        self.assertIn("M3 = 2*M2", h3_metric.note)
        self.assertIn("not an h3 decay threshold", h3_metric.note)

    def test_invisible_width_summary_uses_optional_provenance(self):
        header = HEADER + ["higgs_invisible_widths_included"]
        closed = list(ROWS[0]) + [True]
        open_modelled = list(ROWS[1]) + [True]
        open_modelled[HEADER.index("M3")] = 40.0
        open_unmodelled = list(ROWS[2]) + [False]
        open_unmodelled[HEADER.index("M3")] = 50.0

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "provenance.tsv"
            write_fixture(
                path,
                header=header,
                rows=[closed, open_modelled, open_unmodelled],
            )
            data = self.plotter.load_scan(path)
            rows = self.plotter.build_summary(data)

        by_metric = {row.metric: row for row in rows}
        self.assertEqual(by_metric["higgs_invisible_decay_open"].count, 2)
        self.assertEqual(
            by_metric["higgs_invisible_decay_open_but_unmodelled"].count,
            1,
        )

    def test_cli_defaults_and_overrides(self):
        args = self.plotter.parse_args(["points.tsv"])
        self.assertEqual(args.format, "both")
        self.assertEqual(args.dpi, 200)
        self.assertIsNone(args.scan_metadata)
        self.assertTrue(str(args.output_dir).endswith("plots/points_constraints"))

        args = self.plotter.parse_args(
            [
                "points.tsv",
                "--format",
                "png",
                "--dpi",
                "72",
                "--output-dir",
                "out",
                "--scan-metadata",
                "scan.json",
            ]
        )
        self.assertEqual(args.format, "png")
        self.assertEqual(args.dpi, 72)
        self.assertEqual(args.output_dir, Path("out"))
        self.assertEqual(args.scan_metadata, Path("scan.json"))

    def test_run_orchestrates_all_figures_and_summary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_path = tmpdir / "points.tsv"
            output_dir = tmpdir / "out"
            write_fixture(input_path)

            def fake_standalone(data, spec, target, plot_format, dpi):
                paths = []
                for extension in self.plotter.extensions_for_format(plot_format):
                    path = target / f"{spec.stem}.{extension}"
                    path.write_bytes(b"figure")
                    paths.append(path)
                return paths

            def fake_dashboard(data, stem, plot_stems, target, plot_format, dpi):
                paths = []
                for extension in self.plotter.extensions_for_format(plot_format):
                    path = target / f"{stem}.{extension}"
                    path.write_bytes(b"figure")
                    paths.append(path)
                return paths

            with mock.patch.object(
                self.plotter, "render_standalone", side_effect=fake_standalone
            ), mock.patch.object(
                self.plotter, "render_dashboard", side_effect=fake_dashboard
            ):
                paths = self.plotter.run(
                    [
                        str(input_path),
                        "--output-dir",
                        str(output_dir),
                        "--format",
                        "both",
                        "--dpi",
                        "72",
                    ]
                )

            self.assertEqual(len(paths), 92)
            self.assertEqual(len(set(paths)), 92)
            self.assertTrue((output_dir / "constraint_summary.tsv").exists())
            index_path = output_dir / "index.html"
            self.assertTrue(index_path.exists())
            index_text = index_path.read_text(encoding="utf-8")
            self.assertIn("TRSM constraint plot suite", index_text)
            self.assertIn("Observed-range fallback", index_text)
            self.assertIn("Observed stored rows only", index_text)
            self.assertIn("200 – 550", index_text)
            self.assertIn('href="constraint_summary.tsv"', index_text)
            self.assertIn(
                'src="01_dm_experimental_fourway_m2_m3.png"', index_text
            )
            self.assertIn(
                'href="dashboard_status_summary.pdf"', index_text
            )
            self.assertIn("No stored BSMPT evaluations.", index_text)
            for path in paths:
                self.assertIn(f'href="{path.name}"', index_text)
            loaded = self.plotter.load_scan(input_path)
            self.assertEqual(
                {path.name for path in paths},
                {
                    path.name
                    for path in self.plotter.expected_figure_paths(
                        output_dir, "both", data=loaded
                    )
                },
            )

    def test_run_adds_bsmpt_figures_only_when_results_exist(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_path = tmpdir / "points_bsmpt.tsv"
            output_dir = tmpdir / "out"
            write_fixture(input_path, BSMPT_HEADER, bsmpt_fixture_rows())

            def fake_standalone(data, spec, target, plot_format, dpi):
                paths = []
                for extension in self.plotter.extensions_for_format(plot_format):
                    path = target / f"{spec.stem}.{extension}"
                    path.write_bytes(b"figure")
                    paths.append(path)
                return paths

            def fake_dashboard(data, stem, plot_stems, target, plot_format, dpi):
                paths = []
                for extension in self.plotter.extensions_for_format(plot_format):
                    path = target / f"{stem}.{extension}"
                    path.write_bytes(b"figure")
                    paths.append(path)
                return paths

            with mock.patch.object(
                self.plotter, "render_standalone", side_effect=fake_standalone
            ), mock.patch.object(
                self.plotter, "render_dashboard", side_effect=fake_dashboard
            ):
                paths = self.plotter.run(
                    [
                        str(input_path),
                        "--output-dir",
                        str(output_dir),
                        "--format",
                        "both",
                        "--dpi",
                        "72",
                    ]
                )

            self.assertEqual(len(paths), 110)
            self.assertTrue(
                (output_dir / "dashboard_bsmpt_summary.png").exists()
            )
            self.assertTrue((output_dir / "30_bsmpt_status_m2_m3.pdf").exists())
            index_text = (output_dir / "index.html").read_text(encoding="utf-8")
            self.assertIn("BSMPT phase-transition candidates", index_text)
            self.assertIn("6 attempted BSMPT evaluations", index_text)
            self.assertIn('src="dashboard_bsmpt_summary.png"', index_text)
            data = self.plotter.load_scan(input_path)
            self.assertEqual(
                {path.name for path in paths},
                {
                    path.name
                    for path in self.plotter.expected_figure_paths(
                        output_dir, "both", data=data
                    )
                },
            )

    def test_plot_index_handles_pdf_only_and_skipped_figures(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            data = self.load_fixture(tmpdir)
            pdf_path = tmpdir / "dashboard_status_summary.pdf"
            pdf_path.write_bytes(b"pdf")
            skipped = [
                ("14_indirect_ratio_m2_m3", "No finite values < threshold")
            ]
            summary_rows = self.plotter.build_summary(data, skipped)
            index_path = tmpdir / "index.html"
            self.plotter.write_plot_index(
                index_path, data, [pdf_path], summary_rows, skipped
            )

            text = index_path.read_text(encoding="utf-8")
            self.assertIn('href="dashboard_status_summary.pdf"', text)
            self.assertNotIn('src="dashboard_status_summary.png"', text)
            self.assertIn("Preview unavailable for PDF-only output.", text)
            self.assertIn("No finite values &lt; threshold", text)
            self.assertIn("skipped_figure_14_indirect_ratio_m2_m3", text)

    def test_plot_index_includes_signal_definition_and_provenance(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            data = self.load_signal_fixture(tmpdir)
            preview = tmpdir / "dashboard_signal_summary.png"
            preview.write_bytes(b"png")
            summary = self.plotter.build_summary(data)
            index_path = tmpdir / "index.html"
            self.plotter.write_plot_index(
                index_path,
                data,
                [preview],
                summary,
            )

            text = index_path.read_text(encoding="utf-8")
            self.assertIn("Full-viability collider signal plots", text)
            self.assertIn("k2^2 * sigma_P^YR4", text)
            self.assertIn("10--3000 GeV", text)
            self.assertIn(self.plotter.YR4_SIGNAL_SOURCE_URL, text)
            self.assertIn('src="dashboard_signal_summary.png"', text)
            self.assertIn("NWA", text)
            self.assertIn("before acceptance", text)

    def test_headline_png_smoke(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            data = self.load_fixture(tmpdir)
            spec = self.plotter.PLOT_BY_STEM["01_dm_experimental_fourway_m2_m3"]
            paths = self.plotter.render_standalone(
                data,
                spec,
                tmpdir,
                plot_format="png",
                dpi=72,
            )

            self.assertEqual(paths, [tmpdir / f"{spec.stem}.png"])
            self.assertTrue(paths[0].exists())
            self.assertGreater(paths[0].stat().st_size, 1024)
            self.assertEqual(paths[0].read_bytes()[:8], b"\x89PNG\r\n\x1a\n")

    def test_cumulative_plot_png_smoke(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            data = self.load_fixture(tmpdir)
            spec = self.plotter.PLOT_BY_STEM[
                "26_cumulative_constraints_m2_a12"
            ]
            paths = self.plotter.render_standalone(
                data,
                spec,
                tmpdir,
                plot_format="png",
                dpi=72,
            )

            self.assertEqual(paths, [tmpdir / f"{spec.stem}.png"])
            self.assertTrue(paths[0].exists())
            self.assertGreater(paths[0].stat().st_size, 1024)
            self.assertEqual(paths[0].read_bytes()[:8], b"\x89PNG\r\n\x1a\n")

    def test_signal_dashboard_png_smoke(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            data = self.load_signal_fixture(tmpdir)
            stem = "dashboard_signal_summary"
            paths = self.plotter.render_dashboard(
                data,
                stem,
                self.plotter.DASHBOARDS[stem],
                tmpdir,
                plot_format="png",
                dpi=72,
            )

            self.assertEqual(paths, [tmpdir / f"{stem}.png"])
            self.assertTrue(paths[0].exists())
            self.assertGreater(paths[0].stat().st_size, 4096)
            self.assertEqual(paths[0].read_bytes()[:8], b"\x89PNG\r\n\x1a\n")

    def test_bsmpt_dashboard_png_smoke(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            data = self.load_bsmpt_fixture(tmpdir)
            stem = "dashboard_bsmpt_summary"
            paths = self.plotter.render_dashboard(
                data,
                stem,
                self.plotter.DASHBOARDS[stem],
                tmpdir,
                plot_format="png",
                dpi=72,
            )

            self.assertEqual(paths, [tmpdir / f"{stem}.png"])
            self.assertTrue(paths[0].exists())
            self.assertGreater(paths[0].stat().st_size, 4096)
            self.assertEqual(paths[0].read_bytes()[:8], b"\x89PNG\r\n\x1a\n")


if __name__ == "__main__":
    unittest.main()
