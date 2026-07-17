import math
import unittest

from singlet_EWPO import (
    WMASS_MH2_MAX,
    WMASS_MH2_MIN,
    check_wmass_tania,
    sinthetalim_wmass,
)


class TestWmassConstraint(unittest.TestCase):
    def test_uses_tabulated_interpolated_limit(self):
        self.assertAlmostEqual(float(sinthetalim_wmass(380.0)), 0.214)
        self.assertTrue(check_wmass_tania(380.0, 0.214))
        self.assertFalse(check_wmass_tania(380.0, 0.215))

    def test_is_symmetric_under_mixing_sign(self):
        self.assertTrue(check_wmass_tania(380.0, 0.20))
        self.assertTrue(check_wmass_tania(380.0, -0.20))
        self.assertFalse(check_wmass_tania(380.0, 0.22))
        self.assertFalse(check_wmass_tania(380.0, -0.22))

    def test_passes_with_warning_outside_table_mass_range(self):
        self.assertTrue(check_wmass_tania(WMASS_MH2_MIN, 0.0))
        self.assertTrue(check_wmass_tania(WMASS_MH2_MAX, 0.0))
        with self.assertWarnsRegex(RuntimeWarning, "W-mass constraint not applied"):
            self.assertTrue(check_wmass_tania(WMASS_MH2_MIN - 1.0, 1.0))
        with self.assertWarnsRegex(RuntimeWarning, "W-mass constraint not applied"):
            self.assertTrue(check_wmass_tania(WMASS_MH2_MAX + 1.0, 1.0))

    def test_rejects_non_finite_inputs(self):
        self.assertFalse(check_wmass_tania(math.nan, 0.0))
        self.assertFalse(check_wmass_tania(380.0, math.nan))
        self.assertFalse(check_wmass_tania(math.inf, 0.0))


if __name__ == "__main__":
    unittest.main()
