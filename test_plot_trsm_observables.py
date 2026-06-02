import importlib.util
import math
import tempfile
import unittest
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().parent / "plot_trsm_observables.py"


def load_module():
    spec = importlib.util.spec_from_file_location("plot_trsm_observables", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_sample_points(path):
    path.write_text(
        "\n".join(
            [
                "M2\tM3\tvs\ta12\tlX\tlPhiX\tlSX\tewpt_ew_true_over_T\tdm_omega\tdm_dir_det\tdm_relic_excluded",
                "100\t400\t50\t-0.10\t0.8\t0.3\t0.2\t0.4\t0.02\t1.0e-11\tFalse",
                "200\t500\t60\t0.20\t0.9\t0.2\t0.1\t0.8\t0.10\t2.0e-11\tTrue",
                "300\t600\t70\t-0.30\t1.0\t0.1\t0.05\tnan\t0.30\t3.0e-11\tTrue",
            ]
        )
        + "\n",
        encoding="ascii",
    )


class TestPlotTRSMObservables(unittest.TestCase):
    def test_load_rows_and_filter_numeric_plot_points(self):
        plotter = load_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            points_path = Path(tmpdir) / "points.tsv"
            write_sample_points(points_path)

            rows = plotter.load_rows(points_path)
            filtered = plotter.rows_with_finite_columns(
                rows,
                ["M2", "ewpt_ew_true_over_T", "dm_omega"],
            )

        self.assertEqual(len(rows), 3)
        self.assertEqual(len(filtered), 2)
        self.assertEqual(filtered[0]["M2"], "100")

    def test_derived_observables_are_added_as_columns(self):
        plotter = load_module()

        rows = [
            {
                "M2": "100",
                "M3": "400",
                "vs": "50",
                "a12": "-0.25",
                "lX": "0.8",
                "lPhiX": "0.3",
                "lSX": "0.2",
                "dm_omega": "0.01",
            },
            {
                "M2": "100",
                "M3": "0",
                "vs": "50",
                "a12": "nan",
                "lX": "0.8",
                "lPhiX": "0.3",
                "lSX": "0.2",
                "dm_omega": "-1",
            },
        ]
        plotter.add_derived_observables(rows)

        self.assertEqual(rows[0]["M2_over_M3"], "0.25")
        self.assertEqual(rows[0]["abs_a12"], "0.25")
        self.assertEqual(rows[0]["log10_dm_omega"], "-2")
        expected_lphi = (
            plotter.DEFAULT_M1**2 * math.cos(-0.25) ** 2
            + 100.0**2 * math.sin(-0.25) ** 2
        ) / (2.0 * plotter.SM_VEV**2)
        self.assertAlmostEqual(float(rows[0]["lPhi"]), expected_lphi)
        self.assertEqual(rows[1]["M2_over_M3"], "nan")
        self.assertEqual(rows[1]["log10_dm_omega"], "nan")
        self.assertEqual(rows[1]["lPhi"], "nan")

    def test_eq418_observables_match_ewpt_formula(self):
        plotter = load_module()

        rows = [
            {
                "M2": "160",
                "M3": "500",
                "vs": "50",
                "a12": "-0.18",
                "lX": "0.8",
                "lPhiX": "0.3",
                "lSX": "0.28",
            }
        ]
        plotter.add_derived_observables(rows)
        row = rows[0]

        m1 = plotter.DEFAULT_M1
        m2 = 160.0
        vs = 50.0
        a12 = -0.18
        lphi = (
            m1**2 * math.cos(a12) ** 2
            + m2**2 * math.sin(a12) ** 2
        ) / (2.0 * plotter.SM_VEV**2)
        ls = (
            m2**2 * math.cos(a12) ** 2
            + m1**2 * math.sin(a12) ** 2
        ) / (2.0 * vs**2)
        lphis = ((-m1**2 + m2**2) * math.cos(a12) * math.sin(a12)) / (
            plotter.SM_VEV * vs
        )
        c_phi = 0.5 * lphi
        c_x = 0.5 * 0.8
        c_s = 0.5 * ls
        c_phix = 0.5 * 0.3
        c_phis = 0.5 * lphis
        c_xs = 0.5 * 0.28
        expected = {
            "eq418_c_phiphi": c_phi,
            "eq418_c_xx": c_x,
            "eq418_c_ss": c_s,
            "eq418_C_phix": c_phi * c_x - c_phix**2,
            "eq418_C_phis": c_phi * c_s - c_phis**2,
            "eq418_C_xs": c_x * c_s - c_xs**2,
            "eq418_D": (
                c_phi * c_x * c_s
                + 2.0 * c_phix * c_phis * c_xs
                - c_phis**2 * c_x
                - c_phi * c_xs**2
                - c_phix**2 * c_s
            ),
        }

        for column, value in expected.items():
            self.assertAlmostEqual(float(row[column]), value)

    def test_eq418_plot_group_resolves_to_separate_plots(self):
        plotter = load_module()

        specs = plotter.resolve_plot_specs(
            plotter.parse_args(["points.tsv", "--plot", "eq418_all"])
        )

        self.assertEqual(
            [spec.x for spec in specs],
            [
                "eq418_c_phiphi",
                "eq418_c_xx",
                "eq418_c_ss",
                "eq418_C_phix",
                "eq418_C_phis",
                "eq418_C_xs",
                "eq418_D",
            ],
        )
        self.assertTrue(all(spec.y == "ewpt_ew_true_over_T" for spec in specs))

    def test_eq418_log_plot_group_resolves_to_symlog_x_plots(self):
        plotter = load_module()

        specs = plotter.resolve_plot_specs(
            plotter.parse_args(["points.tsv", "--plot", "eq418_all_log"])
        )

        self.assertEqual(len(specs), 7)
        self.assertTrue(all(spec.x_symlog for spec in specs))
        self.assertTrue(all(not spec.x_log for spec in specs))
        self.assertTrue(
            all(spec.output_stem.endswith("_symlogx") for spec in specs)
        )
        self.assertEqual(specs[-1].x, "eq418_D")

    def test_log_plot_rows_drop_nonpositive_log_axis_values(self):
        plotter = load_module()

        rows = [
            {"x": "1.0", "y": "0.1"},
            {"x": "0.0", "y": "0.2"},
            {"x": "-2.0", "y": "0.3"},
            {"x": "nan", "y": "0.4"},
        ]
        spec = plotter.PlotSpec(name="test", x="x", y="y", x_log=True)

        filtered = plotter.rows_for_plot(rows, spec)

        self.assertEqual(filtered, [rows[0]])

    def test_symlog_plot_rows_keep_negative_axis_values(self):
        plotter = load_module()

        rows = [
            {"x": "1.0", "y": "0.1"},
            {"x": "0.0", "y": "0.2"},
            {"x": "-2.0", "y": "0.3"},
            {"x": "nan", "y": "0.4"},
        ]
        spec = plotter.PlotSpec(name="test", x="x", y="y", x_symlog=True)

        filtered = plotter.rows_for_plot(rows, spec)

        self.assertEqual(filtered, rows[:3])

    def test_cli_override_can_request_symlog_axis(self):
        plotter = load_module()

        spec = plotter.resolve_plot_spec(
            plotter.parse_args(["points.tsv", "--x-symlog"])
        )

        self.assertTrue(spec.x_symlog)
        self.assertFalse(spec.x_log)

    def test_load_rows_adds_derived_observables_by_default(self):
        plotter = load_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            points_path = Path(tmpdir) / "points.tsv"
            write_sample_points(points_path)

            rows = plotter.load_rows(points_path)

        self.assertIn("M2_over_M3", rows[0])
        self.assertEqual(rows[0]["M2_over_M3"], "0.25")

    def test_default_preset_is_ewpt_strength_vs_m2(self):
        plotter = load_module()

        spec = plotter.resolve_plot_spec(
            plotter.parse_args(["points.tsv", "--plot", "ewpt_vs_M2"])
        )

        self.assertEqual(spec.x, "M2")
        self.assertEqual(spec.y, "ewpt_ew_true_over_T")
        self.assertEqual(spec.output_stem, "ewpt_ew_true_over_T_vs_M2")

    def test_cli_overrides_make_new_scatter_plot_without_editing_code(self):
        plotter = load_module()

        spec = plotter.resolve_plot_spec(
            plotter.parse_args(
                [
                    "points.tsv",
                    "--x",
                    "M2_over_M3",
                    "--y",
                    "dm_omega",
                    "--color-by",
                    "log10_dm_omega",
                    "--size-by",
                    "M2",
                    "--marker-by",
                    "dm_relic_excluded",
                    "--marker",
                    "s",
                    "--size",
                    "42",
                    "--output-stem",
                    "custom",
                ]
            )
        )

        self.assertEqual(spec.x, "M2_over_M3")
        self.assertEqual(spec.y, "dm_omega")
        self.assertEqual(spec.color_by, "log10_dm_omega")
        self.assertEqual(spec.size_by, "M2")
        self.assertEqual(spec.marker_by, "dm_relic_excluded")
        self.assertEqual(spec.marker, "s")
        self.assertEqual(spec.size, 42.0)
        self.assertEqual(spec.output_stem, "custom")

    def test_plot_writes_png_when_matplotlib_is_available(self):
        if importlib.util.find_spec("matplotlib") is None:
            self.skipTest("matplotlib is not installed")
        plotter = load_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            points_path = tmpdir / "points.tsv"
            write_sample_points(points_path)

            paths = plotter.run(
                [
                    str(points_path),
                    "--output-dir",
                    str(tmpdir),
                    "--format",
                    "png",
                ]
            )

            self.assertEqual(paths, [tmpdir / "ewpt_ew_true_over_T_vs_M2.png"])
            self.assertTrue(paths[0].exists())


if __name__ == "__main__":
    unittest.main()
