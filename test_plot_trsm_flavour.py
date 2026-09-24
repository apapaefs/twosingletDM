import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

import plot_trsm_constraint_suite as plots
from test_plot_trsm_constraint_suite import HEADER, ROWS, write_fixture
from trsm_flavour import FLAVOUR_COLUMNS, generated_flavour_updates
from test_trsm_flavour import base_brs


def mixed_fixture(path):
    header = HEADER + [c for c in FLAVOUR_COLUMNS if c not in HEADER]
    rows = []
    cases = [(5, .01, .1, .8), (5.2, 1., 1., 0.), (5.4, 1., 0., 1.),
             (5.6, 1., .1, .8), (100, .1, .1, .8), (6, .1, None, None), (6.5, 0., 0., 0.)]
    for mass, mixing, mu, tau in cases:
        dm_mass = 20 if mass == 100 else 1
        point = dict(zip(HEADER, ROWS[3]))
        point.update(M2=mass, M3=dm_mass)
        brs = None if mu is None else base_brs(mu, tau)
        point.update(generated_flavour_updates(mass, dm_mass, (1, mixing, 0),
                                              (None, brs, None), (None, 1, 0)))
        rows.append(["nan" if point.get(key) is None else point[key] for key in header])
    write_fixture(path, header=header, rows=rows)


class TestFlavourPlots(unittest.TestCase):
    def test_legacy_data_is_unassessed_and_cannot_pass_full_viability(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "old.tsv"
            index = HEADER.index("flavour")
            write_fixture(path, header=HEADER[:index] + HEADER[index+1:],
                          rows=[r[:index] + r[index+1:] for r in ROWS])
            data = plots.load_scan(path)
        self.assertFalse(np.any(data.b("flavour_available")))
        self.assertFalse(np.any(data.b("full_viability")))
        self.assertFalse(np.any(data.b("non_dm_viability")))
        self.assertEqual(np.count_nonzero(data.b("pre_flavour_full_viability")), 2)
        summary = {s.metric: s for s in plots.build_summary(data)}
        self.assertEqual(summary["flavour_unassessed"].count, len(ROWS))
        for stem in ("70_flavour_status_m2_m3", "71_flavour_status_lowmass"):
            self.assertIsNone(plots.spec_unavailable_reason(data, plots.PLOT_BY_STEM[stem]))
            fig, ax = plots.plt.subplots()
            try:
                plots.render_spec(fig, ax, data, plots.PLOT_BY_STEM[stem])
            finally:
                plots.plt.close(fig)
        self.assertIsNotNone(plots.spec_unavailable_reason(data, plots.PLOT_BY_STEM["73_flavour_lepton_products"]))

    def test_categories_and_full_selection_use_stored_nullable_verdicts(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "mixed.tsv"
            mixed_fixture(path)
            data = plots.load_scan(path)
        self.assertEqual(data.s("flavour_category").tolist(),
                         ["passed", "mumu", "tautau", "both", "outside_coverage", "unassessed", "zero_signal"])
        np.testing.assert_array_equal(data.b("full_viability"), [True, False, False, False, True, False, True])
        summary = {s.metric: s for s in plots.build_summary(data)}
        self.assertEqual(summary["flavour_removed_otherwise_viable"].count, 3)
        self.assertEqual(summary["flavour_unassessed_otherwise_viable"].count, 1)
        self.assertEqual(data.f("flavour_h2_k2sq")[1], 1.)

    def test_each_new_renderer_and_zero_annotations(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "mixed.tsv"
            mixed_fixture(path)
            data = plots.load_scan(path)
            for spec in plots.PLOT_SPECS:
                if not spec.kind.startswith("flavour_"):
                    continue
                with self.subTest(spec=spec.stem):
                    self.assertIsNone(plots.spec_unavailable_reason(data, spec))
                    fig, ax = plots.plt.subplots()
                    try:
                        plots.render_spec(fig, ax, data, spec)
                        fig.canvas.draw()
                        if spec.kind == "flavour_products":
                            self.assertEqual(len(fig.axes), 2)
                            for axis in fig.axes:
                                self.assertEqual(axis.get_yscale(), "log")
                                self.assertIn("Zero predictions omitted", axis.get_title())
                        if spec.kind == "flavour_zoom":
                            self.assertEqual(ax.get_xlim(), (4., 9.2))
                    finally:
                        plots.plt.close(fig)

    def test_uniform_passing_status_is_not_omitted(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "pass.tsv"
            write_fixture(path)
            data = plots.load_scan(path)
        spec = plots.PLOT_BY_STEM["70_flavour_status_m2_m3"]
        self.assertIsNone(plots.spec_unavailable_reason(data, spec))
        fig, ax = plots.plt.subplots()
        try:
            plots.render_spec(fig, ax, data, spec)
            self.assertIn("Passing: 8", [t.get_text() for t in ax.get_legend().get_texts()])
        finally:
            plots.plt.close(fig)


if __name__ == "__main__":
    unittest.main()
