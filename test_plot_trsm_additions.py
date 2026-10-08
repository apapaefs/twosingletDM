"""Regression coverage for cumulative DM and associated-production diagnostics."""

import tempfile
import unittest
from pathlib import Path

import numpy as np

import plot_trsm_constraint_suite as plots
from test_plot_trsm_constraint_suite import HEADER, ROWS, write_fixture
from trsm_mg5_rates import derive_mg5_rates


class PlotAdditionsTests(unittest.TestCase):
    def load_rows(self, updates):
        records = []
        for index, changes in enumerate(updates):
            record = dict(zip(HEADER, ROWS[-1]))
            record.update(M2=200 + index, M3=20 + index)
            record.update(changes)
            records.append(record)
        columns = list(dict.fromkeys(name for record in records for name in record))
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "points.tsv"
            write_fixture(source, columns, [[r.get(c, "nan") for c in columns] for r in records])
            return plots.load_scan(source)

    def test_cumulative_stages_require_known_verdicts_and_preserve_coverage(self):
        common = dict(dm_cmb_enabled=False, dm_cmb_available=False, dm_cmb_excluded="nan")
        cases = [
            {}, {"dm_relic_excluded": True}, {"dm_direct_detection_excluded": True},
            {"dm_indirect_detection_excluded": True},
            dict(dm_cmb_enabled=True, dm_cmb_available=True, dm_cmb_excluded=True),
            dict(dm_cmb_enabled=True, dm_cmb_available=False),
            {"dm_indirect_available": False},
            {"dm_indirect_detection_excluded": "nan"},
            {"dm_relic_excluded": "nan"}, {"flavour": False},
            dict(evo=False, constraint_version="trsm_constraints_v2", experimental_subset=True),
            {"evo": False},
            dict(dm_cmb_enabled=True, dm_cmb_available=True, dm_cmb_excluded="nan"),
        ]
        data = self.load_rows([common | case for case in cases])
        masks = plots.cumulative_dm_masks(data)
        self.assertEqual([np.count_nonzero(mask) for mask in masks.values()], [11, 9, 8, 6, 3])
        np.testing.assert_array_equal(np.flatnonzero(masks["cmb"]), [0, 6, 10])
        for before, after in zip(list(masks.values()), list(masks.values())[1:]):
            self.assertFalse(np.any(after & ~before))
        summary = {row.metric: row for row in plots.build_summary(data)}
        self.assertEqual(summary["cumulative_dm_cmb"].denominator, 11)
        self.assertEqual(summary["cumulative_dm_vs_stored_full_mismatch"].count, 8)
        self.assertEqual(summary["cumulative_dm_cmb_unavailable"].count, 2)

    def test_legacy_cmb_stage_is_identity_and_missing_baseline_is_unavailable(self):
        data = self.load_rows([{}, {"dm_indirect_available": False}])
        masks = plots.cumulative_dm_masks(data)
        np.testing.assert_array_equal(masks["indirect"], masks["cmb"])
        self.assertIn("not enabled", plots.dm_stage_label(data, "cmb"))
        unavailable = self.load_rows([{"flavour": "nan"}])
        spec = plots.PLOT_BY_STEM["75_cumulative_dm_m2_m3"]
        self.assertIn("no non-DM", plots.spec_unavailable_reason(unavailable, spec))

    def test_missing_gluon_rate_is_not_a_zero_total_and_valid_zeros_survive(self):
        data = self.load_rows([
            dict(mg5_xsec_pp_eta0Z_pb=tree, mg5_xsec_gg_eta0Z_pb=gluon, h2_h3h3_br=br)
            for tree, gluon, br in ((2, 3, .5), (2, "nan", .5), (0, 0, .5),
                                     ("nan", 3, .5), (-1, 3, .5), (2, 3, 1.2), (2, 3, 0))
        ])
        np.testing.assert_allclose(data.f("mg5_xsec_pp_eta0Z_total_pb"),
                                   [5, np.nan, 0, np.nan, np.nan, 5, 5], equal_nan=True)
        np.testing.assert_allclose(data.f("mono_z_total_xsec_pb"),
                                   [2.5, np.nan, 0, np.nan, np.nan, np.nan, 0], equal_nan=True)
        np.testing.assert_allclose(data.f("mono_z_xsec_pb"),
                                   [1, 1, 0, np.nan, np.nan, np.nan, 0], equal_nan=True)

    def test_stored_finite_rates_preserved_and_nan_rates_backfilled(self):
        data = self.load_rows([
            dict(mg5_xsec_pp_eta0Z_pb=2, mg5_xsec_gg_eta0Z_pb=3, h2_h3h3_br=.5,
                 mono_z_xsec_pb=old, mono_z_total_xsec_pb=total)
            for old, total in ((9, 12), ("nan", "nan"))
        ])
        np.testing.assert_allclose(data.f("mono_z_xsec_pb"), [9, 1])
        np.testing.assert_allclose(data.f("mono_z_total_xsec_pb"), [12, 2.5])

    def test_vectorized_rates_match_production_helper(self):
        cases = ((2, 3, .5), (2, np.nan, .5), (0, 0, 1), (-1, 2, .5),
                 (2, 3, 1.1), (1e308, 1e308, .5))
        raw = {
            "mg5_xsec_pp_eta0Z_pb": np.asarray([c[0] for c in cases]),
            "mg5_xsec_gg_eta0Z_pb": np.asarray([c[1] for c in cases]),
            "h2_h3h3_br": np.asarray([c[2] for c in cases]),
        }
        vectorized = plots.derive_associated_rates(raw)
        for index in range(len(cases)):
            expected = derive_mg5_rates({column: values[index] for column, values in raw.items()})
            for column, values in vectorized.items():
                np.testing.assert_allclose(values[index], expected.get(column, np.nan), equal_nan=True)

    def test_new_renderers_select_expected_points(self):
        data = self.load_rows([
            dict(mg5_xsec_pp_eta0Z_pb=2, mg5_xsec_gg_eta0Z_pb=3, h2_h3h3_br=.5),
            dict(mg5_xsec_pp_eta0Z_pb=2, mg5_xsec_gg_eta0Z_pb="nan", h2_h3h3_br=.5),
            dict(dm=False, dm_relic_excluded=True,
                 mg5_xsec_pp_eta0Z_pb=2, mg5_xsec_gg_eta0Z_pb=3, h2_h3h3_br=.5),
        ])
        self.assertIsNotNone(plots.signal_availability_reason(data))
        for stem in ("80_zh2_production_vs_m2", "81_zh2_production_vs_m3"):
            spec = plots.PLOT_BY_STEM[stem]
            self.assertIsNone(plots.spec_unavailable_reason(data, spec))
            fig, ax = plots.plt.subplots()
            try:
                plots.render_zh2_rates_xy(ax, data, spec)
                self.assertEqual([len(c.get_offsets()) for c in ax.collections], [2, 1, 1])
                self.assertEqual(ax.get_yscale(), "log")
                fig.canvas.draw()
            finally:
                plots.plt.close(fig)
        spec = plots.PLOT_BY_STEM["76_cumulative_dm_m3_k133"]
        fig, ax = plots.plt.subplots()
        try:
            plots.render_cumulative_dm_xy(ax, data, spec)
            self.assertEqual([len(c.get_offsets()) for c in ax.collections], [3, 2, 2, 2, 2])
            self.assertEqual(ax.get_yscale(), "symlog")
            fig.canvas.draw()
        finally:
            plots.plt.close(fig)

    def test_html_groups_full_viability_mg5_independently_of_yr4(self):
        data = self.load_rows([dict(mg5_xsec_gg_heta0_pb=4, mg5_xsec_pp_eta0Z_pb=2,
                                    mg5_xsec_gg_eta0Z_pb=3, h2_h3h3_br=.5,
                                    ewpt_status="success", ewpt_ew_true_over_T=2)])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "index.html"
            plots.write_plot_index(path, data, [], plots.build_summary(data))
            html = path.read_text()
        signal = html.split('<section id="signal">', 1)[1].split('<section id="bsmpt">', 1)[0]
        standalone = html.split('<section id="standalone">', 1)[1]
        for stem in ("52_mono_higgs_xsec_vs_m2", "55_mono_z_xsec_vs_m3",
                     "84_mono_z_total_xsec_vs_m2", "dashboard_mg5_mono_rates"):
            self.assertEqual(html.count(f'id="{stem}"'), 1)
            self.assertIn(f'id="{stem}"', signal)
            self.assertNotIn(f'id="{stem}"', standalone)
        self.assertIn('id="60_mono_higgs_xsec_no_dm_vs_m2"', standalone)
        self.assertIn("YR4 signal plots unavailable", signal)
        self.assertIn('id="bsmpt-full-viability"', html)
        self.assertEqual(html.count('id="30_bsmpt_status_m2_m3_full_viability"'), 1)


if __name__ == "__main__":
    unittest.main()
