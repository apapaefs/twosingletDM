import math
import unittest

import numpy as np

from generate_trsm_info import (
    BR_interpolators_SM,
    HIGGSTOOLS_YR4_LOWMASS_SWITCH_GEV,
    PORTAL_CONVENTION_ID,
    Gam_h2_to_h1h1,
    calc_h2_BRs,
    generate_lams,
    higgstools_yr4_lowmass_br_widths,
    scalar_to_identical_scalar_width,
    vxzero_invisible_decay_info,
    vxzero_portal_couplings,
)


class TestCanonicalInvisibleWidths(unittest.TestCase):
    def test_lowmass_reference_matches_higgstools_at_campaign_row(self):
        values = higgstools_yr4_lowmass_br_widths(8.114)
        expected = np.array(
            [
                0.0,
                0.35115206577,
                0.0017078770609,
                0.43290941348,
                0.0034772359181,
                0.0,
                0.21070729756,
                4.6108495634e-05,
                0.0,
                1.7977728154e-09,
                6.0329023875e-10,
                3.4784141706e-05,
            ]
        )

        self.assertTrue(np.all(np.isfinite(values)))
        np.testing.assert_allclose(values, expected, rtol=1.0e-11, atol=1.0e-18)

    def test_campaign_low_m2_row_has_finite_base_widths(self):
        result = generate_lams(
            64,
            8.114,
            338.8,
            221.9,
            0.0,
            0.1415,
            0.0,
            0.0,
            False,
            lX=0.4781,
            lPhiX=0.0035056476325440854,
            lSX=-0.0006017021941833114,
        )
        h2_brs = result[23]

        self.assertTrue(np.all(np.isfinite(h2_brs)))
        self.assertGreater(h2_brs[-1], 0.0)
        self.assertAlmostEqual(result[8], h2_brs[-1])

    def test_lowmass_switch_preserves_legacy_high_mass_interpolators(self):
        boundary = HIGGSTOOLS_YR4_LOWMASS_SWITCH_GEV
        below = higgstools_yr4_lowmass_br_widths(np.nextafter(boundary, 0.0))
        legacy = np.asarray(
            [float(interpolator(boundary)) for interpolator in BR_interpolators_SM]
        )
        at_boundary = calc_h2_BRs(
            BR_interpolators_SM,
            125.09,
            boundary,
            0.2,
            0.0,
        )

        # The tracked pre-existing table and HiggsTools YR4 source differ
        # slightly at 20 GeV.  Keep the checkpoint-compatible legacy branch at
        # the boundary, while bounding and documenting the one-sided step.
        np.testing.assert_array_equal(at_boundary[:-2], legacy[:-1])
        self.assertAlmostEqual(at_boundary[-1], 0.2**2 * legacy[-1])
        self.assertLess(np.max(np.abs(below[:-1] - legacy[:-1])), 3.6e-4)
        self.assertLess(abs(below[-1] - legacy[-1]), 6.0e-7)

    def test_portal_couplings_pure_lsx_benchmark(self):
        k133, k233 = vxzero_portal_couplings(0.0, 0.1, 500.0, 0.0)

        self.assertEqual(PORTAL_CONVENTION_ID, "trsm_vxzero_canonical_v1")
        self.assertAlmostEqual(k133, 0.0)
        self.assertAlmostEqual(k233, 25.0)

    def test_identical_scalar_width_benchmark_and_sign_invariance(self):
        expected = 25.0**2 * math.sqrt(1.0 - 4.0 * 50.0**2 / 300.0**2) / (
            8.0 * math.pi * 300.0
        )

        self.assertAlmostEqual(
            scalar_to_identical_scalar_width(50.0, 300.0, 25.0),
            expected,
        )
        self.assertAlmostEqual(
            scalar_to_identical_scalar_width(50.0, 300.0, -25.0),
            expected,
        )
        self.assertAlmostEqual(Gam_h2_to_h1h1(50.0, 300.0, 25.0), expected)

    def test_identical_scalar_width_is_zero_at_and_below_threshold(self):
        self.assertEqual(scalar_to_identical_scalar_width(50.0, 100.0, 25.0), 0.0)
        self.assertEqual(scalar_to_identical_scalar_width(51.0, 100.0, 25.0), 0.0)

    def test_invisible_decay_info_adds_to_base_widths(self):
        info = vxzero_invisible_decay_info(
            125.09,
            300.0,
            50.0,
            0.0,
            25.0,
            0.004,
            1.2,
        )

        expected_gamma2 = scalar_to_identical_scalar_width(50.0, 300.0, 25.0)
        self.assertEqual(info["h1_h3h3_width"], 0.0)
        self.assertAlmostEqual(info["h2_h3h3_width"], expected_gamma2)
        self.assertAlmostEqual(info["w1"], 0.004)
        self.assertAlmostEqual(info["w2"], 1.2 + expected_gamma2)
        self.assertAlmostEqual(
            info["h2_h3h3_br"],
            expected_gamma2 / (1.2 + expected_gamma2),
        )

    def test_generate_lams_returns_physical_totals_but_base_br_arrays(self):
        result = generate_lams(
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
        w1, w2 = result[7], result[8]
        k133 = result[18]
        h1_brs, h2_brs = result[22], result[23]
        _expected_k133, expected_k233 = vxzero_portal_couplings(
            0.0, 0.1, 500.0, 0.0
        )
        info = vxzero_invisible_decay_info(
            125.09,
            300.0,
            50.0,
            k133,
            expected_k233,
            h1_brs[-1],
            h2_brs[-1],
        )

        self.assertAlmostEqual(w1, info["w1"])
        self.assertAlmostEqual(w2, info["w2"])
        self.assertGreater(w2, h2_brs[-1])
        self.assertEqual(h2_brs[-1], 0.0)
        self.assertTrue(np.all(np.isfinite(h2_brs)))
        self.assertTrue(np.all(h2_brs[:-1] == 0.0))

    def test_non_vxzero_widths_remain_the_base_widths(self):
        result = generate_lams(
            1,
            300.0,
            600.0,
            500.0,
            100.0,
            0.1,
            0.2,
            0.3,
            False,
        )

        self.assertAlmostEqual(result[7], result[22][-1])
        self.assertAlmostEqual(result[8], result[23][-1])


if __name__ == "__main__":
    unittest.main()
