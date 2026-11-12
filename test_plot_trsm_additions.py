"""Regression coverage for cumulative DM and associated-production diagnostics."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

import plot_trsm_constraint_suite as plots
from test_plot_trsm_constraint_suite import HEADER, ROWS, write_fixture
from trsm_mg5_rates import derive_mg5_rates


class PlotAdditionsTests(unittest.TestCase):
    def load_rows(self, updates, metadata=None, **kwargs):
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
            if metadata is not None:
                source.with_suffix(".metadata.json").write_text(json.dumps(metadata))
            return plots.load_scan(source, **kwargs)

    def test_portal_intersection_uses_versioned_experiment_and_aggregate_dm(self):
        data = self.load_rows([
            {}, {"dm": False}, {"hb": False},
            dict(constraint_version="trsm_constraints_v2", experimental_subset=True, hb=False),
            dict(constraint_version="trsm_constraints_v2", experimental_subset=False),
            dict(thc=False, flavour=False),
        ])
        np.testing.assert_array_equal(np.flatnonzero(data.b("experimental_dm")), [0, 3, 5])
        np.testing.assert_array_equal(np.flatnonzero(data.b("full_viability")), [0, 3])
        spec = plots.PLOT_BY_STEM["86_k133_experimental_dm_m2_m3"]
        fig, ax = plots.plt.subplots()
        try:
            plots.render_spec(fig, ax, data, spec)
            self.assertEqual([len(c.get_offsets()) for c in ax.collections], [6, 3])
            fig.canvas.draw()
        finally:
            plots.plt.close(fig)

    def test_portal_dashboard_couplings_are_not_swapped(self):
        data = self.load_rows([dict(K133=-.03, K233=21), dict(K133=.07, K233=-8)])
        stems = plots.DASHBOARDS["dashboard_portal_mass_plane_summary_experimental_dm"]
        for stem, column in zip(stems, ("K133", "K233")):
            spec = plots.PLOT_BY_STEM[stem]
            self.assertEqual(spec.value, column)
            fig, ax = plots.plt.subplots()
            try:
                plots.render_spec(fig, ax, data, spec)
                colored = ax.collections[1]
                expected = dict(zip(data.f("M2"), data.f(column)))
                for xy, value in zip(colored.get_offsets(), colored.get_array()):
                    self.assertEqual(value, expected[xy[0]])
            finally:
                plots.plt.close(fig)

    def test_determinant_matches_eq418_and_displays_both_signs(self):
        from test_trsm_ewpt import check_eq_4_18
        common = dict(vs=300, a12=.1, lX=.3, lPhiX=.01, lSX=.02, ewpt_status="assessed",
                      ewpt_ew_entry_jump_over_T=1.4, ewpt_baryo_candidate=True,
                      ewpt_gw_max_field_jump_over_T=2, ewpt_gw_candidate=True)
        data = self.load_rows([common, common | dict(lPhiX=2), common | dict(lX=0, lPhiX=0, lSX=0)])
        expected = [check_eq_4_18(dict(m1=plots.INPUT_M1, m2=data.f("M2")[i], vs=300, a12=.1,
                                     lx=data.f("lX")[i], lphix=data.f("lPhiX")[i], lsx=data.f("lSX")[i])).couplings["D"]
                    for i in range(3)]
        np.testing.assert_allclose(data.f("ewpt_determinant_D"), expected, rtol=1e-14)
        self.assertGreater(expected[0], 0)
        self.assertLess(expected[1], 0)
        self.assertEqual(expected[2], 0)
        spec = plots.DETERMINANT_SPECS[0]
        fig, ax = plots.plt.subplots()
        try:
            plots.render_spec(fig, ax, data, spec)
            self.assertEqual(ax.get_xscale(), "symlog")
            self.assertEqual(sum(len(c.get_offsets()) for c in ax.collections), 3)
            fig.canvas.draw()
        finally:
            plots.plt.close(fig)

    def test_determinants_match_mass_basis_portal_identity(self):
        from plot_trsm_observables import eq418_couplings
        from test_trsm_ewpt import check_eq_4_18

        common = dict(M2=600, vs=200, lX=.1, lPhiX=0, lSX=0)
        data = self.load_rows([
            common | dict(a12=0),
            common | dict(a12=.2),
            common | dict(a12=-.2, lPhiX=.05, lSX=.1),
            common | dict(a12=.2, lPhiX=2, lSX=.3),
            common | dict(a12=.2, lX=0),
        ])
        m1, v = plots.INPUT_M1, plots.INPUT_VEV
        m2, vs = data.f("M2"), data.f("vs")
        c, s = np.cos(data.f("a12")), np.sin(data.f("a12"))
        k133 = (data.f("lPhiX") * v * c - data.f("lSX") * vs * s) / 2
        k233 = (data.f("lPhiX") * v * s + data.f("lSX") * vs * c) / 2
        # Schur complement expressed in physical masses and h_i X^2 coefficients.
        expected = (m1**2 * m2**2 / (4 * v**2 * vs**2)
                    * (data.f("lX") - 2 * k133**2 / m1**2 - 2 * k233**2 / m2**2))
        np.testing.assert_allclose(data.f("ewpt_determinant_D"), expected, rtol=1e-12, atol=1e-14)
        for i, value in enumerate(expected):
            row = {key: data.f(key)[i] for key in ("M2", "vs", "a12", "lX", "lPhiX", "lSX")}
            with self.subTest(point=i):
                np.testing.assert_allclose(eq418_couplings(row)["eq418_D"], value, rtol=1e-12, atol=1e-14)
                check = check_eq_4_18(dict(m1=m1, m2=row["M2"], vs=row["vs"], a12=row["a12"],
                                           lx=row["lX"], lphix=row["lPhiX"], lsx=row["lSX"]))
                np.testing.assert_allclose(check.couplings["D"], value, rtol=1e-12, atol=1e-14)

    def test_calcgw_uses_paired_transitions_and_parent_selection(self):
        from test_trsm_gw import native_row
        from trsm_gw import parse_gw_row
        native = native_row()
        native["h2Omega_2_sw_3"] = 0
        result = parse_gw_row(native)
        entry = dict(gw_status="success", gw_transitions=json.dumps(result["gw_transitions"]))
        data = self.load_rows([entry, entry | dict(dm=False)])
        self.assertEqual(len(data.gw["alpha"]), 4)
        np.testing.assert_array_equal(data.gw["row_index"], [0, 0, 1, 1])
        spec = plots.PLOT_BY_STEM["93_calcgw_alpha_beta_full_viability"]
        self.assertIsNone(plots.spec_unavailable_reason(data, spec))
        fig, ax = plots.plt.subplots()
        try:
            plots.render_spec(fig, ax, data, spec)
            np.testing.assert_allclose(ax.collections[0].get_offsets(), [[.04, 200], [.7, 30]])
            self.assertEqual(ax.get_xscale(), "log")
            self.assertEqual(ax.get_yscale(), "log")
            fig.canvas.draw()
        finally:
            plots.plt.close(fig)
        spec = plots.PLOT_BY_STEM["96_calcgw_sound_wave_amplitude_full_viability"]
        x, y, selected, valid = plots.gw_plot_values(data, spec)
        self.assertEqual(np.count_nonzero(selected), 2)
        self.assertEqual(np.count_nonzero(valid), 1)
        unavailable = self.load_rows([{}])
        self.assertIn("no successful CalcGW", plots.spec_unavailable_reason(unavailable, spec))
        with self.assertRaisesRegex(ValueError, "invalid gw_transitions"):
            self.load_rows([dict(gw_transitions="{bad}")])

    def test_calcgw_stars_follow_34b_parent_points_in_every_panel(self):
        from test_trsm_gw import native_row
        from trsm_gw import parse_gw_row
        transitions = parse_gw_row(native_row())["gw_transitions"]
        entry = dict(gw_status="success", gw_transitions=json.dumps(transitions),
                     ewpt_status="success", ewpt_ew_entry_jump_over_T=2)
        data = self.load_rows([
            entry, entry | dict(dm=False),
            entry | dict(ewpt_ew_entry_jump_over_T=1),
            entry | dict(ewpt_ew_entry_jump_over_T="nan"),
            entry | dict(ewpt_status="failed"),
        ])
        np.testing.assert_array_equal(data.b("bsmpt_strong_ew_entry"), [True, True, False, False, False])
        for base in plots.CALCGW_SPECS + plots.CALCGW_SNR_ZOOM_SPECS:
            for suffix, star_masses in (("", [200, 200, 201, 201]),
                                        ("_full_viability", [200, 200])):
                spec = plots.PLOT_BY_STEM[base.stem + suffix]
                with self.subTest(stem=spec.stem):
                    fig, ax = plots.plt.subplots()
                    try:
                        plots.render_spec(fig, ax, data, spec)
                        self.assertEqual(len(ax.collections), 2)
                        stars = ax.collections[1]
                        self.assertEqual(len(stars.get_offsets()), len(star_masses))
                        self.assertEqual(len(stars.get_paths()[0].vertices), 11)
                        if spec.x == "M2":
                            np.testing.assert_allclose(stars.get_offsets()[:, 0], star_masses)
                        self.assertIs(ax.collections[0].norm, stars.norm)
                        if spec.y == "snr_lisa":
                            np.testing.assert_allclose(ax.lines[0].get_ydata(), [5, 5])
                        if spec.y_limits is not None:
                            self.assertEqual(ax.get_ylim(), (1, 100))
                        fig.canvas.draw()
                    finally:
                        plots.plt.close(fig)

    def test_calcgw_zoom_preserves_values_and_reports_points_in_view(self):
        from test_trsm_gw import native_row
        from trsm_gw import parse_gw_row
        transitions = parse_gw_row(native_row())["gw_transitions"]
        transitions[0]["snr_lisa"] = .1
        transitions[1]["snr_lisa"] = 200
        data = self.load_rows([dict(gw_status="success", gw_transitions=json.dumps(transitions),
                                   ewpt_status="success", ewpt_ew_entry_jump_over_T=2)])
        spec = plots.CALCGW_SNR_ZOOM_SPECS[0]
        fig, ax = plots.plt.subplots()
        try:
            plots.render_spec(fig, ax, data, spec)
            np.testing.assert_allclose(ax.collections[0].get_offsets()[:, 1], [.1, 200])
            self.assertEqual(ax.get_ylim(), (1, 100))
            self.assertIn("In view: transitions = 0; points = 0", ax.texts[-1].get_text())
            self.assertIn("(points: 0)", ax.get_legend().get_texts()[0].get_text())
        finally:
            plots.plt.close(fig)

    def test_hh_normalization_missing_values_and_full_viability_selection(self):
        data = self.load_rows([
            {"mg5_xsec_hh_pb": value} for value in (.01617, .1617, 0, "nan", -1, "inf")
        ] + [dict(mg5_xsec_hh_pb=1.617, dm=False)])
        np.testing.assert_allclose(data.f("hh_enhancement"),
                                   [1, 10, 0, np.nan, np.nan, np.nan, 100], equal_nan=True)
        summary = {row.metric: row for row in plots.build_summary(data)}
        self.assertEqual(summary["hh_full_viability_rate_available"].count, 3)
        self.assertEqual(summary["hh_full_viability_rate_available"].denominator, 6)
        self.assertEqual(summary["hh_enhanced_above_sm"].count, 1)
        for stem in plots.DASHBOARDS["dashboard_hh_enhancement"]:
            spec = plots.PLOT_BY_STEM[stem]
            self.assertIsNone(plots.spec_unavailable_reason(data, spec))
            fig, ax = plots.plt.subplots()
            try:
                plots.render_spec(fig, ax, data, spec)
                self.assertEqual(len(ax.collections[0].get_offsets()), 2)
                self.assertIn("16.17 fb", ax.get_title())
                self.assertIn("unavailable = 3", ax.texts[-1].get_text())
                if spec.kind == "hh_enhancement_xy":
                    np.testing.assert_allclose(ax.collections[0].get_offsets()[:, 1], [1, 10])
                    self.assertEqual(ax.get_yscale(), "log")
                    np.testing.assert_allclose(ax.lines[0].get_ydata(), [1, 1])
                fig.canvas.draw()
            finally:
                plots.plt.close(fig)

    def test_hh_never_substitutes_factorized_resonance_and_handles_empty_rates(self):
        spec = plots.PLOT_BY_STEM["88_hh_enhancement_m2_m3"]
        data = self.load_rows([dict(xs136_lo_h2_pb=10, h2_h1h1_br=.5)])
        self.assertIn("--mg5-process hh", plots.spec_unavailable_reason(data, spec))
        self.assertNotIn("hh_enhancement", data.derived)
        for changes in ({"mg5_xsec_hh_pb": "nan"}, {"mg5_xsec_hh_pb": 0},
                        dict(mg5_xsec_hh_pb=.1, dm=False)):
            with self.subTest(changes=changes):
                data = self.load_rows([changes])
                self.assertIn("no finite positive", plots.spec_unavailable_reason(data, spec))

    def test_hh_energy_receipt_and_custom_reference(self):
        # Legacy enrichment receipts have no plotting metadata schema.
        for metadata in (
            {"mg5_enrichment": {"configuration": {"energy_TeV": 13}}},
            {"fixed_parameters": [{"variable": "sqrt(s)", "value": 13000, "unit": "GeV"}]},
            {"configuration": {"fixed": {"Energy": 13}}},
        ):
            data = self.load_rows([dict(mg5_xsec_hh_pb=.1)], metadata=metadata)
            self.assertIn("require 13.6 TeV", data.hh_error)
            self.assertNotIn("hh_enhancement", data.derived)
        data = self.load_rows([dict(mg5_xsec_hh_pb=.1)], hh_sm_xsec_pb=.02,
                              metadata={"mg5_enrichment": {"configuration": {"energy_TeV": 13.6}}})
        self.assertIsNone(data.hh_error)
        np.testing.assert_allclose(data.f("hh_enhancement"), [5])
        for value in (0, -1, float("nan"), float("inf")):
            with self.assertRaisesRegex(ValueError, "finite and positive"):
                self.load_rows([dict(mg5_xsec_hh_pb=.1)], hh_sm_xsec_pb=value)

    def test_targeted_dashboard_run_includes_panels_without_other_suites(self):
        selected = ("dashboard_portal_mass_plane_summary_experimental_dm", "dashboard_hh_enhancement")
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "points.tsv"
            write_fixture(source, HEADER + ["mg5_xsec_hh_pb"], [list(row) + [.01617] for row in ROWS])
            def standalone(data, spec, *args):
                return [Path(directory) / (spec.stem + ".png")]
            def dashboard(data, stem, *args):
                return [Path(directory) / (stem + ".png")]
            with mock.patch.object(plots, "render_standalone", side_effect=standalone), \
                 mock.patch.object(plots, "render_dashboard", side_effect=dashboard), \
                 mock.patch("plot_trsm_v2.render_suite") as v2:
                paths = plots.run([str(source), "--output-dir", directory, "--format", "png",
                                   "--plots", *selected])
            v2.assert_not_called()
            self.assertEqual({p.stem for p in paths},
                             set(selected) | {stem for name in selected for stem in plots.DASHBOARDS[name]})
            html = (Path(directory) / "index.html").read_text()
            self.assertIn("16.17 fb", html)
            self.assertIn("continuum or its interference", html)

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
