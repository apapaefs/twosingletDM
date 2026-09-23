import json
import math
import unittest
from pathlib import Path

from dm_thermal_relic_diagnostic import (
    resonance_proximity_updates,
    thermal_vev_updates,
)


def history(*rows):
    return {"minimatracer": {"global_branch": [
        {"temp": temp, "w1": ew, "ws": singlet, "phase_index": phase}
        for temp, ew, singlet, phase in rows
    ]}}


class TestThermalRelicDiagnostic(unittest.TestCase):
    def test_smooth_window_interpolates_within_one_phase(self):
        payload = history(
            (0, 246, 300, 0),
            (5, 245, 299, 0),
            (10, 240, 290, 0),
            (20, 220, 270, 0),
            (30, 200, 250, 0),
        )
        result = thermal_vev_updates(payload, 7.5)
        self.assertEqual(result["dm_relic_thermal_vev_window_low_T_GeV"], 3.75)
        self.assertEqual(result["dm_relic_thermal_vev_window_high_T_GeV"], 15)
        self.assertAlmostEqual(result["dm_relic_thermal_ew_vev_Tf_over_T0"], 242.5 / 246)
        self.assertAlmostEqual(result["dm_relic_thermal_s_vev_Tf_over_T0"], 294.5 / 300)
        self.assertAlmostEqual(result["dm_relic_thermal_ew_vev_max_fractional_shift"], 16 / 246)
        self.assertAlmostEqual(result["dm_relic_thermal_s_vev_max_fractional_shift"], 20 / 300)
        self.assertIs(result["dm_relic_thermal_phase_boundary_bracket_overlaps_window"], False)
        self.assertIs(result["dm_relic_thermal_vev_shift_ge_10pct"], False)

    def test_phase_boundary_is_not_interpolated_or_called_low_sensitivity(self):
        payload = history(
            (0, 246, 300, 0),
            (8, 245, 299, 0),
            (12, 244, 298, 1),
            (20, 243, 297, 1),
        )
        result = thermal_vev_updates(payload, 10)
        self.assertIsNone(result["dm_relic_thermal_ew_vev_Tf_over_T0"])
        self.assertIs(result["dm_relic_thermal_phase_boundary_bracket_overlaps_window"], True)
        self.assertIsNone(result["dm_relic_thermal_vev_shift_ge_10pct"])

    def test_large_shift_and_missing_coverage(self):
        payload = history((0, 246, 300, 0), (5, 240, 290, 0), (10, 210, 250, 0))
        result = thermal_vev_updates(payload, 5)
        self.assertIs(result["dm_relic_thermal_vev_shift_ge_10pct"], True)
        mild = history((0, 246, 300, 0), (5, 245, 299, 0), (10, 244, 298, 0))
        partial = thermal_vev_updates(mild, 8)
        self.assertIsNone(partial["dm_relic_thermal_phase_boundary_bracket_overlaps_window"])
        self.assertIsNone(partial["dm_relic_thermal_vev_shift_ge_10pct"])
        self.assertIsNone(thermal_vev_updates({}, 5)["dm_relic_thermal_vev_shift_ge_10pct"])
        unknown_phase = history((0, 246, 300, 0), (5, 245, 299, None), (10, 244, 298, 0))
        self.assertIsNone(thermal_vev_updates(unknown_phase, 5)["dm_relic_thermal_vev_shift_ge_10pct"])

    def test_resonance_gaps_keep_sign_and_use_freezeout_scale(self):
        result = resonance_proximity_updates(120, 60, 5)
        self.assertAlmostEqual(result["dm_resonance_h1_mass_gap_GeV"], 5.09)
        self.assertEqual(result["dm_resonance_h2_mass_gap_GeV"], 0)
        self.assertAlmostEqual(result["dm_resonance_h1_abs_gap_over_Tf"], abs(5.09) / 5)
        self.assertEqual(result["dm_resonance_h2_abs_gap_over_Tf"], 0)
        self.assertEqual(result["dm_resonance_nearest_mediator"], "h2")
        missing_tf = resonance_proximity_updates(120, 60)
        self.assertEqual(missing_tf["dm_resonance_nearest_mediator"], "h2")
        self.assertIsNone(missing_tf["dm_resonance_h2_abs_gap_over_Tf"])

    def test_archived_point_223615_is_stable_on_equilibrium_branch(self):
        archive = Path(__file__).resolve().parent / "benchmarks/v2/point_223615/legacy-thermal-history.json"
        payload = json.loads(archive.read_text(encoding="utf-8"))
        result = thermal_vev_updates(payload, 23.51 / 31.2)
        self.assertTrue(math.isclose(result["dm_relic_thermal_ew_vev_Tf_over_T0"], 1, rel_tol=1e-5))
        self.assertTrue(math.isclose(result["dm_relic_thermal_s_vev_Tf_over_T0"], 1, rel_tol=1e-5))
        self.assertIs(result["dm_relic_thermal_vev_shift_ge_10pct"], False)


if __name__ == "__main__":
    unittest.main()
