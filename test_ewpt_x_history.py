import json
import unittest

from ewpt_x_history import x_history_updates


def payload(*points):
    return {"minimatracer": {"global_branch": [
        {"temp": temperature, "label": label}
        for temperature, label in points
    ]}}


class TestXHistory(unittest.TestCase):
    def setUp(self):
        self.history = payload(
            (0, "EW"), (20, "EW"), (30, "EW_X_BROKEN"),
            (70, "EW_X_BROKEN"), (80, "EW"), (100, "SYM"),
        )

    def test_global_x_window_and_freezeout_overlap(self):
        early = x_history_updates(self.history, 40)
        self.assertEqual(early["ewpt_x_broken_min_T_GeV"], 30)
        self.assertEqual(early["ewpt_x_broken_max_T_GeV"], 70)
        self.assertEqual(json.loads(early["ewpt_x_broken_intervals_GeV"]), [[30, 70]])
        self.assertEqual(early["ewpt_x_final_restoration_low_T_GeV"], 20)
        self.assertEqual(early["ewpt_x_final_restoration_high_T_GeV"], 30)
        self.assertEqual(early["ewpt_x_phase_at_freezeout"], "broken")
        self.assertIs(early["ewpt_x_broken_at_or_after_freezeout"], True)
        self.assertIs(early["dm_relic_z2_freezeout_compatible"], False)

        late = x_history_updates(self.history, 10)
        self.assertEqual(late["ewpt_x_phase_at_freezeout"], "unbroken")
        self.assertIs(late["ewpt_x_broken_at_or_after_freezeout"], False)
        self.assertIs(late["dm_relic_z2_freezeout_compatible"], True)

    def test_boundary_and_later_x_breaking_are_distinguished(self):
        boundary = x_history_updates(self.history, 25)
        self.assertEqual(boundary["ewpt_x_phase_at_freezeout"], "boundary_unresolved")
        self.assertIsNone(boundary["ewpt_x_broken_at_or_after_freezeout"])
        self.assertIsNone(boundary["dm_relic_z2_freezeout_compatible"])

        before_breaking = x_history_updates(self.history, 90)
        self.assertEqual(before_breaking["ewpt_x_phase_at_freezeout"], "unbroken")
        self.assertIs(before_breaking["ewpt_x_broken_at_or_after_freezeout"], True)
        self.assertIs(before_breaking["dm_relic_z2_freezeout_compatible"], False)

    def test_local_x_minimum_does_not_count_without_global_x_branch(self):
        history = payload((0, "EW"), (10, "EW"), (100, "SYM"))
        history["minimatracer"]["phase_traces"] = [
            {"samples": [{"temp": 0, "wx": 60}]},
        ]
        updates = x_history_updates(history, 5)
        self.assertEqual(json.loads(updates["ewpt_x_broken_intervals_GeV"]), [])
        self.assertIs(updates["ewpt_x_broken_at_or_after_freezeout"], False)
        self.assertIs(updates["dm_relic_z2_freezeout_compatible"], True)
        self.assertIsNone(updates["ewpt_x_final_restoration_low_T_GeV"])

    def test_multiple_windows_and_missing_inputs(self):
        history = payload(
            (0, "EW"), (10, "X_BROKEN"), (20, "EW"),
            (40, "X_BROKEN"), (50, "SYM"),
        )
        updates = x_history_updates(history, None)
        self.assertEqual(json.loads(updates["ewpt_x_broken_intervals_GeV"]), [[10, 10], [40, 40]])
        self.assertEqual(updates["ewpt_x_final_restoration_low_T_GeV"], 0)
        self.assertEqual(updates["ewpt_x_final_restoration_high_T_GeV"], 10)
        self.assertIsNone(updates["ewpt_x_phase_at_freezeout"])
        self.assertIsNone(updates["dm_relic_z2_freezeout_compatible"])
        self.assertIsNone(x_history_updates({})["ewpt_x_broken_intervals_GeV"])

    def test_x_broken_at_zero_has_no_final_restoration(self):
        history = payload((0, "X_BROKEN"), (10, "X_BROKEN"), (20, "EW"))
        updates = x_history_updates(history, 5)
        self.assertIsNone(updates["ewpt_x_final_restoration_low_T_GeV"])
        self.assertIsNone(updates["ewpt_x_final_restoration_high_T_GeV"])
        self.assertIs(updates["dm_relic_z2_freezeout_compatible"], False)


if __name__ == "__main__":
    unittest.main()
