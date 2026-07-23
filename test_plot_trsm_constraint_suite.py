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
    "ewpt_status",
    "ewpt_global_phase_path",
    "ewpt_has_x_broken",
    "ewpt_ew_step_index",
]
BSMPT_HEADER = HEADER + BSMPT_EXTRA_HEADER


def bsmpt_fixture_rows():
    payloads = [
        ("nan", "nan", "nan", "nan", "nan"),
        ("nan", "failed", "nan", "nan", "nan"),
        ("nan", "success", "SYM -> EW", False, 0),
        (0.75, "success", "SYM -> SINGLET_S -> EW", False, 1),
        (1.25, "success", "SYM -> X_BROKEN -> EW_X_BROKEN", True, 1),
        (2.0, "success", "SYM -> EW", False, 0),
        (0.4, "success", "SYM -> MIXED -> EW", False, 1),
        ("nan", "nan", "nan", "nan", "nan"),
    ]
    rows = []
    for original, payload in zip(ROWS, payloads):
        strength, status, phase_path, has_x_broken, ew_step = payload
        row = list(original)
        row[HEADER.index("ewpt_ew_true_over_T")] = strength
        row.extend([status, phase_path, has_x_broken, ew_step])
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

    def test_strict_bool_accepts_only_canonical_tokens(self):
        self.assertIs(self.plotter.strict_bool("True"), True)
        self.assertIs(self.plotter.strict_bool("False"), False)
        for value in ["true", "FALSE", "1", "0", "yes", "", " True ", "nan"]:
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "expected 'True' or 'False'"):
                    self.plotter.strict_bool(value, "dm", 9)

    def test_nullable_dm_boolean_accepts_only_canonical_nan(self):
        self.assertIs(self.plotter.strict_nullable_bool("True"), True)
        self.assertIs(self.plotter.strict_nullable_bool("False"), False)
        self.assertIsNone(self.plotter.strict_nullable_bool("nan"))
        for value in ["NaN", "none", "", " nan "]:
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
                "selected weak FOPT",
                "selected strong FOPT",
                "selected strong FOPT",
                "selected weak FOPT",
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
        self.assertEqual(len(self.plotter.PLOT_SPECS), 36)
        self.assertEqual(len(self.plotter.DASHBOARDS), 5)
        stems = self.plotter.all_figure_stems()
        self.assertEqual(len(stems), 41)
        self.assertEqual(len(set(stems)), 41)
        self.assertTrue(any("bsmpt" in stem for stem in stems))
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
                "32_bsmpt_phase_history_m2_m3",
                "33_bsmpt_ew_entry_step_m2_m3",
                "34_bsmpt_strength_vs_m2",
                "35_bsmpt_strength_vs_m3",
                "36_bsmpt_counts",
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
                    "31_bsmpt_ew_true_over_t_m2_m3",
                    "32_bsmpt_phase_history_m2_m3",
                    "33_bsmpt_ew_entry_step_m2_m3",
                    "34_bsmpt_strength_vs_m2",
                    "36_bsmpt_counts",
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

        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_fixture(tmpdir)
        np.testing.assert_allclose(data.f("abs_K233"), np.abs(data.f("K233")))
        _categories, binary_styles = self.plotter.category_styles(data, "dm")
        self.assertNotEqual(binary_styles["fail"].marker, binary_styles["pass"].marker)

        paths = self.plotter.expected_figure_paths(Path("plots"), "both")
        self.assertEqual(len(paths), 82)
        self.assertEqual(len(set(paths)), 82)
        self.assertEqual(sum(path.suffix == ".png" for path in paths), 41)
        self.assertEqual(sum(path.suffix == ".pdf" for path in paths), 41)

        legacy_paths = self.plotter.expected_figure_paths(
            Path("plots"), "both", data=data
        )
        self.assertEqual(len(legacy_paths), 66)
        self.assertFalse(any("bsmpt" in path.stem for path in legacy_paths))

        with tempfile.TemporaryDirectory() as tmpdir:
            bsmpt_data = self.load_bsmpt_fixture(tmpdir)
        bsmpt_paths = self.plotter.expected_figure_paths(
            Path("plots"), "both", data=bsmpt_data
        )
        self.assertEqual(len(bsmpt_paths), 82)
        self.assertTrue(any(path.stem == "dashboard_bsmpt_summary" for path in bsmpt_paths))

    def test_summary_contains_categories_omissions_and_high_mass_tail(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_fixture(tmpdir)
            rows = self.plotter.build_summary(data)
            by_metric = {row.metric: row for row in rows}

            self.assertEqual(by_metric["total_rows"].count, 8)
            self.assertEqual(by_metric["experimental"].count, 4)
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

            self.assertEqual(len(paths), 66)
            self.assertEqual(len(set(paths)), 66)
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

            self.assertEqual(len(paths), 82)
            self.assertTrue(
                (output_dir / "dashboard_bsmpt_summary.png").exists()
            )
            self.assertTrue((output_dir / "30_bsmpt_status_m2_m3.pdf").exists())
            index_text = (output_dir / "index.html").read_text(encoding="utf-8")
            self.assertIn("BSMPT electroweak phase-transition plots", index_text)
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
