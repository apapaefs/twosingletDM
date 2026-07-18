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
]

ROWS = [
    [200, 325, -0.1, -100, 100, -2, True, True, False, False, True, False, False, 0.2, 0.1, 2e-9, 1e-9, True, True, False, 0, False, 25, "nan"],
    [250, 500, 0.2, -10, 10, -1, True, True, True, True, True, False, True, 0.05, 0.1, 0.5e-9, 1e-9, False, False, False, 0, False, 16, "nan"],
    [300, 600, -0.3, -1, 1, -0.1, True, True, True, True, True, True, False, 0.3, 0.1, 0.25e-9, 1e-9, True, False, True, 0.25, False, 9, "nan"],
    [350, 700, 0.4, -0.1, 0.1, 0, True, True, True, True, True, True, True, 0.02, 0.1, 0.1e-9, 1e-9, False, False, True, 0.75, False, 4, "nan"],
    [400, 800, -0.5, 0.1, -0.1, 0.1, True, True, False, True, True, True, False, 0.08, 0.1, 4e-9, 1e-9, False, True, True, 1.25, True, 1, "nan"],
    [450, 900, 0.6, 1, -1, 1, True, True, True, False, True, True, True, 0.01, 0.1, 0.2e-9, 1e-9, False, False, False, 0, False, 0.25, "nan"],
    [500, 1000, -0.7, 10, -10, 2, True, True, True, True, True, True, False, 0.4, 0.1, 3e-9, 1e-9, True, True, True, 0.5, False, 0, "nan"],
    [550, 1100, 0.8, 100, -100, 10, True, True, True, True, True, True, True, 0.1, 0.1, 1e-9, 1e-9, False, False, True, 1.0, False, "nan", "nan"],
]


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

    def test_strict_bool_accepts_only_canonical_tokens(self):
        self.assertIs(self.plotter.strict_bool("True"), True)
        self.assertIs(self.plotter.strict_bool("False"), False)
        for value in ["true", "FALSE", "1", "0", "yes", "", " True ", "nan"]:
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "expected 'True' or 'False'"):
                    self.plotter.strict_bool(value, "dm", 9)

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
                "options": {"independent_m3": True, "nrandom": 10000},
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

    def test_registry_and_expected_paths_are_unique(self):
        self.assertEqual(len(self.plotter.PLOT_SPECS), 22)
        self.assertEqual(len(self.plotter.DASHBOARDS), 3)
        stems = self.plotter.all_figure_stems()
        self.assertEqual(len(stems), 25)
        self.assertEqual(len(set(stems)), 25)
        self.assertFalse(any("ewpt" in stem for stem in stems))
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
            },
        )
        indirect = self.plotter.PLOT_BY_STEM["14_indirect_ratio_m2_m3"]
        self.assertEqual(indirect.value, "indirect_ratio_available")
        self.assertEqual(indirect.norm_kind, "log_to_one")

        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_fixture(tmpdir)
        _categories, binary_styles = self.plotter.category_styles(data, "dm")
        self.assertNotEqual(binary_styles["fail"].marker, binary_styles["pass"].marker)

        paths = self.plotter.expected_figure_paths(Path("plots"), "both")
        self.assertEqual(len(paths), 50)
        self.assertEqual(len(set(paths)), 50)
        self.assertEqual(sum(path.suffix == ".png" for path in paths), 25)
        self.assertEqual(sum(path.suffix == ".pdf" for path in paths), 25)

    def test_summary_contains_categories_omissions_and_high_mass_tail(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = self.load_fixture(tmpdir)
            rows = self.plotter.build_summary(data)
            by_metric = {row.metric: row for row in rows}

            self.assertEqual(by_metric["total_rows"].count, 8)
            self.assertEqual(by_metric["experimental"].count, 4)
            self.assertEqual(by_metric["dm"].count, 4)
            self.assertEqual(by_metric["full_viability"].count, 2)
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
            self.assertEqual(by_metric["ewpt_finite"].count, 0)
            self.assertIn("omitted_plot_evo", by_metric)
            self.assertIn("omitted_plot_thc", by_metric)
            self.assertIn("omitted_plot_ewpo", by_metric)

            summary_path = Path(tmpdir) / "constraint_summary.tsv"
            self.plotter.write_summary(summary_path, rows)
            text = summary_path.read_text(encoding="ascii")
            self.assertTrue(text.startswith("metric\tcount\tdenominator\tpercent\tnote\n"))
            self.assertIn("full_viability\t2\t8\t25", text)

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

            self.assertEqual(len(paths), 50)
            self.assertEqual(len(set(paths)), 50)
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
            for path in paths:
                self.assertIn(f'href="{path.name}"', index_text)
            self.assertEqual(
                {path.name for path in paths},
                {path.name for path in self.plotter.expected_figure_paths(output_dir, "both")},
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


if __name__ == "__main__":
    unittest.main()
