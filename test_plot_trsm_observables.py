import importlib.util
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
                "M2\tM3\ta12\tewpt_ew_true_over_T\tdm_omega\tdm_dir_det\tdm_relic_excluded",
                "100\t400\t-0.10\t0.4\t0.02\t1.0e-11\tFalse",
                "200\t500\t0.20\t0.8\t0.10\t2.0e-11\tTrue",
                "300\t600\t-0.30\tnan\t0.30\t3.0e-11\tTrue",
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
            {"M2": "100", "M3": "400", "a12": "-0.25", "dm_omega": "0.01"},
            {"M2": "100", "M3": "0", "a12": "nan", "dm_omega": "-1"},
        ]
        plotter.add_derived_observables(rows)

        self.assertEqual(rows[0]["M2_over_M3"], "0.25")
        self.assertEqual(rows[0]["abs_a12"], "0.25")
        self.assertEqual(rows[0]["log10_dm_omega"], "-2")
        self.assertEqual(rows[1]["M2_over_M3"], "nan")
        self.assertEqual(rows[1]["log10_dm_omega"], "nan")

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
