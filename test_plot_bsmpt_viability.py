"""Regression checks for BSMPT star meaning and full-viability selections."""

from pathlib import Path
import tempfile
import unittest

import numpy as np
from matplotlib.markers import MarkerStyle

import plot_trsm_constraint_suite as plots
from test_plot_trsm_constraint_suite import BSMPT_HEADER, bsmpt_fixture_rows, write_fixture


def fixture(directory, *, entry_columns=True):
    extras = [
        "ewpt_ew_entry_jump_over_T", "ewpt_baryo_candidate",
        "ewpt_ew_entry_nucl_jump_over_T", "ewpt_ew_entry_perc_jump_over_T",
        "ewpt_gw_max_field_jump_over_T", "ewpt_gw_candidate",
        "ewpt_gw_crit_field_jump_over_T", "ewpt_gw_nucl_field_jump_over_T",
        "ewpt_gw_perc_field_jump_over_T", "dm_freezeout_temperature_GeV",
        "ewpt_x_broken_min_T_GeV", "ewpt_x_broken_max_T_GeV",
        "ewpt_x_broken_at_or_after_freezeout",
    ]
    header = BSMPT_HEADER + extras
    rows = []
    for index, entry in enumerate((0.5, 1.0, 1.2, "nan", 1.5, 2.0)):
        row = list(bsmpt_fixture_rows()[3])
        row[header.index("M2")] = 200 + index * 100
        row[header.index("M3")] = 20 + index * 10
        row[header.index("ewpt_ew_true_over_T")] = 2.0 + index
        row[header.index("ewpt_ew_jump_over_T")] = 0.1 + index
        row[header.index("ewpt_ew_step_index")] = 2
        row[header.index("dm")] = index != 4
        row[header.index("ewpt_status")] = "failed" if index == 5 else "success"
        row.extend([entry, index in (2, 4), 0.6 + index, 0.7 + index,
                    2.5 + index, True, 2.1 + index, 2.2 + index, 2.3 + index,
                    10 + index, 20 + index, 30 + index, False])
        rows.append(row)
    if not entry_columns:
        keep = [i for i, name in enumerate(header)
                if not name.startswith("ewpt_ew_entry_") and name != "ewpt_baryo_candidate"]
        header = [header[i] for i in keep]
        rows = [[row[i] for i in keep] for row in rows]
    path = Path(directory) / "bsmpt.tsv"
    write_fixture(path, header, rows)
    return plots.load_scan(path)


def star_count(ax):
    style = MarkerStyle("*")
    vertices = style.get_path().transformed(style.get_transform()).vertices
    return sum(len(collection.get_offsets()) for collection in ax.collections
               if collection.get_paths()
               and collection.get_paths()[0].vertices.shape == vertices.shape
               and np.allclose(collection.get_paths()[0].vertices, vertices))


class TestBSMPTFullViability(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.data = fixture(self.directory.name)

    def test_strict_entry_mask_and_stored_flags(self):
        np.testing.assert_array_equal(self.data.b("bsmpt_strong_ew_entry"),
                                      [False, False, True, False, True, False])
        np.testing.assert_array_equal(self.data.b("ewpt_baryo_candidate"),
                                      [False, False, True, False, True, False])
        self.assertTrue(np.all(self.data.b("ewpt_gw_candidate")))

    def test_all_bsmpt_specs_have_full_viability_counterparts(self):
        originals = [spec for spec in plots.PLOT_SPECS
                     if spec.requires_bsmpt and spec.selection is None]
        self.assertEqual(len(originals), 19)
        for spec in originals:
            clone = plots.PLOT_BY_STEM[spec.stem + "_full_viability"]
            self.assertEqual(clone.selection, "full_viability")
            self.assertEqual(clone.kind, spec.kind)
        self.assertEqual(plots.DASHBOARDS["dashboard_bsmpt_summary_full_viability"],
                         tuple(stem + "_full_viability" for stem in
                               plots.DASHBOARDS["dashboard_bsmpt_summary"]))

    def test_every_scatter_family_has_only_strong_entry_stars(self):
        for spec in plots.BSMPT_FULL_VIABILITY_SPECS:
            if spec.kind == "bsmpt_bars":
                continue
            with self.subTest(stem=spec.stem):
                self.assertIsNone(plots.spec_unavailable_reason(self.data, spec))
                fig, ax = plots.plt.subplots()
                try:
                    plots.render_spec(fig, ax, self.data, spec)
                    expected = 2 if spec.kind in {
                        "bsmpt_temperature_comparison_xy", "x_window_freezeout_xy"
                    } else 1
                    self.assertEqual(star_count(ax), expected)
                    # The only excluded point is itself strong. A missing selection
                    # would double the star count in every scatter family.
                    if spec.kind in {"categorical_mass", "continuous_mass"}:
                        for collection in ax.collections:
                            self.assertFalse(np.any(collection.get_offsets()[:, 0] == 600))
                    fig.canvas.draw()
                finally:
                    plots.plt.close(fig)

    def test_legacy_strength_and_gw_diagnostics_never_infer_stars(self):
        data = fixture(self.directory.name, entry_columns=False)
        self.assertFalse(np.any(data.b("bsmpt_strong_ew_entry")))
        for stem in ("30_bsmpt_status_m2_m3", "31f_bsmpt_gw_status_m2_m3",
                     "33_bsmpt_ew_entry_step_m2_m3", "34_bsmpt_strength_vs_m2"):
            fig, ax = plots.plt.subplots()
            try:
                plots.render_spec(fig, ax, data, plots.PLOT_BY_STEM[stem])
                self.assertEqual(star_count(ax), 0)
                self.assertTrue(any("classification unavailable" in text.get_text()
                                    for text in ax.texts))
            finally:
                plots.plt.close(fig)

    def test_counts_and_denominators_use_selected_population(self):
        metrics = plots.bsmpt_bar_metrics(self.data, "full_viability")
        self.assertEqual(metrics[0][1], 5)
        self.assertEqual(metrics[1][1], 4)
        self.assertEqual(metrics[2][1], 1)
        self.assertEqual(next(count for label, count, *_ in metrics
                              if label.startswith("EW entry:")), 1)
        fig, ax = plots.plt.subplots()
        try:
            plots.render_bsmpt_bars(ax, self.data,
                plots.PLOT_BY_STEM["36_bsmpt_counts_full_viability"])
            self.assertIn("selected N = 5", ax.get_xlabel())
            self.assertIn("5 (100.00%)", [text.get_text() for text in ax.texts])
        finally:
            plots.plt.close(fig)

    def test_empty_selection_is_unavailable_for_every_counterpart(self):
        self.data.derived["full_viability"][:] = False
        for spec in plots.BSMPT_FULL_VIABILITY_SPECS:
            self.assertIsNotNone(plots.spec_unavailable_reason(self.data, spec))
        self.assertFalse(plots.dashboard_available(self.data,
            plots.DASHBOARDS["dashboard_bsmpt_summary_full_viability"]))

    def test_values_outside_selection_do_not_make_plot_available(self):
        selected = self.data.b("full_viability")
        for key in ("ewpt_ew_entry_jump_over_T", "ewpt_ew_entry_nucl_jump_over_T",
                    "ewpt_ew_entry_perc_jump_over_T", "ewpt_gw_max_field_jump_over_T",
                    "dm_freezeout_temperature_GeV"):
            self.data.floats[key][selected] = np.nan
        for stem in ("31c_bsmpt_ew_entry_jump_over_t_m2_m3",
                     "34b_bsmpt_ew_entry_strength_vs_m2",
                     "35c_bsmpt_selected_vs_ew_entry_jump",
                     "35e_bsmpt_ew_entry_temperature_jumps",
                     "35g_freezeout_vs_x_broken_window"):
            self.assertIsNotNone(plots.spec_unavailable_reason(
                self.data, plots.PLOT_BY_STEM[stem + "_full_viability"]))


if __name__ == "__main__":
    unittest.main()
