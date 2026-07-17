import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

import test_trsm_evolution as evolution


class TestTrsmEvolutionSolver(unittest.TestCase):
    @patch("test_trsm_evolution.solve_ivp")
    def test_solver_requests_only_the_test_scale(self, solve_ivp):
        initial_conditions = [0.1] * 13
        solve_ivp.return_value = SimpleNamespace(
            success=True,
            y=np.ones((13, 1)),
        )

        result = evolution.solv_eq(
            91.0,
            1000.0,
            0.0001,
            initial_conditions,
            900.0,
        )

        self.assertEqual(result, 1)
        solve_ivp.assert_called_once_with(
            evolution.rhs,
            (91.0, 1000.0),
            initial_conditions,
            t_eval=[900.0],
        )

    @patch("test_trsm_evolution.solve_ivp")
    def test_solver_failure_is_rejected_without_reading_a_solution(self, solve_ivp):
        solve_ivp.return_value = SimpleNamespace(success=False)

        self.assertEqual(
            evolution.solv_eq(91.0, 1000.0, 0.0001, [0.1] * 13, 900.0),
            0,
        )

    @patch("test_trsm_evolution.solve_ivp")
    def test_nonfinite_solution_is_rejected(self, solve_ivp):
        values = np.ones((13, 1))
        values[5, 0] = np.nan
        solve_ivp.return_value = SimpleNamespace(success=True, y=values)

        self.assertEqual(
            evolution.solv_eq(91.0, 1000.0, 0.0001, [0.1] * 13, 900.0),
            0,
        )

    @patch("test_trsm_evolution.solve_ivp")
    def test_test_scale_outside_integration_interval_is_rejected(self, solve_ivp):
        self.assertEqual(
            evolution.solv_eq(91.0, 1000.0, 0.0001, [0.1] * 13, 1200.0),
            0,
        )
        solve_ivp.assert_not_called()


if __name__ == "__main__":
    unittest.main()
