import unittest

import Higgs.predictions as HP
import numpy as np

from generate_trsm_info import generate_lams, vxzero_invisible_decay_info, vxzero_portal_couplings
from test_trsm_higgstools import H1, H2, H3, analyze_parampoint, ensure_sum_unit, pred


def base_br_arrays():
    h1_brs = np.zeros(12, dtype=float)
    h1_brs[0] = 1.0
    h1_brs[-1] = 0.004

    h2_brs = np.zeros(13, dtype=float)
    h2_brs[0] = 0.5
    h2_brs[11] = 0.5
    h2_brs[-1] = 2.0
    h3_brs = np.zeros(15, dtype=float)
    return h1_brs, h2_brs, h3_brs


class TestHiggsToolsInvisibleWidths(unittest.TestCase):
    def analyze(self, h1_width, h2_width):
        h1_brs, h2_brs, h3_brs = base_br_arrays()
        return analyze_parampoint(
            pred,
            H1,
            H2,
            H3,
            125.09,
            300.0,
            50.0,
            1.0,
            0.1,
            0.0,
            h1_brs,
            h2_brs,
            h3_brs,
            h1_direct_invisible_width=h1_width,
            h2_direct_invisible_width=h2_width,
        )

    def test_direct_invisible_width_rescales_visible_and_cascade_brs(self):
        self.analyze(0.001, 0.5)

        self.assertAlmostEqual(H1.totalWidth(), 0.005)
        self.assertAlmostEqual(H1.br(HP.Decay.directInv), 0.2)
        self.assertAlmostEqual(H1.br(HP.Decay.bb), 0.8)
        self.assertAlmostEqual(H2.totalWidth(), 2.5)
        self.assertAlmostEqual(H2.br(HP.Decay.directInv), 0.2)
        self.assertAlmostEqual(H2.br("H1", "H1"), 0.4)

    def test_explicit_zero_clears_previous_direct_invisible_widths(self):
        self.analyze(0.001, 0.5)
        self.analyze(0.0, 0.0)

        self.assertAlmostEqual(H1.totalWidth(), 0.004)
        self.assertAlmostEqual(H2.totalWidth(), 2.0)
        self.assertEqual(H1.br(HP.Decay.directInv), 0.0)
        self.assertEqual(H2.br(HP.Decay.directInv), 0.0)

    def test_invalid_direct_invisible_widths_are_rejected(self):
        h1_brs, h2_brs, h3_brs = base_br_arrays()
        common = (
            pred,
            H1,
            H2,
            H3,
            125.09,
            300.0,
            50.0,
            1.0,
            0.1,
            0.0,
            h1_brs,
            h2_brs,
            h3_brs,
        )
        with self.assertRaises(ValueError):
            analyze_parampoint(*common, h1_direct_invisible_width=-0.1)
        with self.assertRaises(ValueError):
            analyze_parampoint(*common, h2_direct_invisible_width=float("nan"))

    def test_nonrounded_base_total_is_preserved_exactly(self):
        h1_brs, h2_brs, h3_brs = base_br_arrays()
        base_width = 0.004123456789123
        invisible_width = 0.000987654321987
        h1_brs[-1] = base_width

        normalized = ensure_sum_unit(h1_brs)
        self.assertEqual(normalized[-1], base_width)
        analyze_parampoint(
            pred,
            H1,
            H2,
            H3,
            125.09,
            300.0,
            50.0,
            1.0,
            0.1,
            0.0,
            h1_brs,
            h2_brs,
            h3_brs,
            h1_direct_invisible_width=invisible_width,
        )

        self.assertEqual(H1.totalWidth(), base_width + invisible_width)
        self.assertEqual(
            H1.br(HP.Decay.directInv),
            invisible_width / (base_width + invisible_width),
        )

    def test_pure_invisible_h2_clears_prior_nonzero_state(self):
        # Seed H2 with visible, cascade, and direct-invisible decays first.
        self.analyze(0.001, 0.5)
        self.assertGreater(H2.br(HP.Decay.bb), 0.0)
        self.assertGreater(H2.br("H1", "H1"), 0.0)

        generated = generate_lams(
            1,
            300.0,
            50.0,
            500.0,
            0.0,
            0.0,
            0.0,
            0.0,
            False,
            lX=0.1,
            lPhiX=0.0,
            lSX=0.1,
        )
        w1, w2 = generated[7], generated[8]
        k133 = generated[18]
        k1, k2, k3 = generated[19:22]
        h1_brs, h2_brs, h3_brs = generated[22:25]
        _canonical_k133, k233 = vxzero_portal_couplings(0.0, 0.1, 500.0, 0.0)
        invisible = vxzero_invisible_decay_info(
            125.09,
            300.0,
            50.0,
            k133,
            k233,
            h1_brs[-1],
            h2_brs[-1],
        )

        self.assertTrue(np.all(np.isfinite(h2_brs)))
        self.assertTrue(np.all(h2_brs == 0.0))
        self.assertEqual(w1, invisible["w1"])
        self.assertEqual(w2, invisible["w2"])
        analyze_parampoint(
            pred,
            H1,
            H2,
            H3,
            125.09,
            300.0,
            50.0,
            k1,
            k2,
            k3,
            h1_brs,
            h2_brs,
            h3_brs,
            h1_direct_invisible_width=invisible["h1_h3h3_width"],
            h2_direct_invisible_width=invisible["h2_h3h3_width"],
        )

        self.assertEqual(H2.totalWidth(), w2)
        self.assertEqual(H2.br(HP.Decay.directInv), 1.0)
        self.assertEqual(H2.br(HP.Decay.bb), 0.0)
        self.assertEqual(H2.br("H1", "H1"), 0.0)


if __name__ == "__main__":
    unittest.main()
