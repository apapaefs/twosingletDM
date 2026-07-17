import math
import unittest

from test_trsm_theory_constraints import (
    theory_constraints,
    theory_constraints_vxzero,
)


class TestTrsmTheoryConstraints(unittest.TestCase):
    def test_general_model_accepts_a_bounded_perturbative_point(self):
        self.assertTrue(
            theory_constraints(
                200.0,
                300.0,
                300.0,
                400.0,
                0.1,
                0.05,
                0.02,
            )
        )

    def test_vxzero_accepts_a_bounded_perturbative_point(self):
        self.assertTrue(
            theory_constraints_vxzero(
                200.0,
                300.0,
                100.0,
                0.1,
                0.1,
                0.001,
                0.001,
            )
        )

    def test_vxzero_rejects_large_negative_portal(self):
        # The former one-sided bounds accepted this point because -10 < 8*pi
        # and boundedness from below was not checked.
        self.assertFalse(
            theory_constraints_vxzero(
                200.0,
                300.0,
                100.0,
                0.1,
                0.1,
                -10.0,
                0.001,
            )
        )

    def test_vxzero_rejects_negative_diagonal_quartic(self):
        self.assertFalse(
            theory_constraints_vxzero(
                200.0,
                300.0,
                100.0,
                0.1,
                -0.1,
                0.001,
                0.001,
            )
        )

    def test_constraints_reject_nonfinite_or_zero_vevs(self):
        self.assertFalse(
            theory_constraints_vxzero(
                0.0,
                300.0,
                100.0,
                0.1,
                0.1,
                0.001,
                0.001,
            )
        )
        self.assertFalse(
            theory_constraints_vxzero(
                200.0,
                300.0,
                100.0,
                math.nan,
                0.1,
                0.001,
                0.001,
            )
        )
        self.assertFalse(theory_constraints(200.0, 0.0, 300.0, 400.0, 0.1, 0.1, 0.1))


if __name__ == "__main__":
    unittest.main()
