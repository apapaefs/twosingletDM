import unittest

from ewpt_entry_criterion import ew_entry_updates
from reprocess_trsm_ewpt import updates_from_payload


def transition(index, kind, temperature, false_vev, true_vev):
    false = dict(zip(("w1", "wx", "ws"), false_vev))
    true = dict(zip(("w1", "wx", "ws"), true_vev))
    return {
        "transition_index": index,
        "temperature_kind": kind,
        "status": "success",
        "temperature": temperature,
        "false_vev": false,
        "true_vev": true,
        "ew_true_over_T": abs(true["w1"]) / temperature,
        "ew_jump_over_T": abs(true["w1"] - false["w1"]) / temperature,
    }


class TestEWPTCandidateCriteria(unittest.TestCase):
    def test_weak_ew_entry_not_replaced_by_later_large_true_vev(self):
        payload = {"transition_strengths": [
            transition(1, "crit", 150, (0, 0, 0), (25, 0, 0)),
            transition(1, "perc", 140, (0, 0, 0), (28, 0, 0)),
            transition(1, "compl", 139, (0, 0, 0), (29, 0, 0)),
            transition(2, "nucl", 80, (230, 0, 0), (232, 0, 0)),
        ]}
        updates = ew_entry_updates(payload)
        self.assertEqual(updates["ewpt_ew_entry_temperature_kind"], "crit")
        self.assertAlmostEqual(updates["ewpt_ew_entry_jump_over_T"], 25 / 150)
        self.assertFalse(updates["ewpt_baryo_candidate"])
        self.assertFalse(updates["ewpt_gw_candidate"])
        self.assertTrue(updates["ewpt_ew_entry_completed"])
        reprocessed = updates_from_payload(payload)
        self.assertGreater(reprocessed["ewpt_ew_true_over_T"], 2)
        self.assertFalse(reprocessed["ewpt_baryo_candidate"])

    def test_critical_only_strong_ew_entry_is_selected(self):
        payload = {"transition_strengths": [
            transition(0, "crit", 100, (0, 0, 0), (125, 0, 0))
        ]}
        updates = ew_entry_updates(payload)
        self.assertEqual(updates["ewpt_ew_entry_jump_over_T"], 1.25)
        self.assertTrue(updates["ewpt_baryo_candidate"])
        self.assertTrue(updates["ewpt_gw_candidate"])
        self.assertFalse(updates["ewpt_ew_entry_percolated"])
        self.assertFalse(updates["ewpt_ew_entry_completed"])
        self.assertIsNone(updates["ewpt_ew_entry_nucl_jump_over_T"])

    def test_nucleation_and_percolation_raise_jump_without_baryo_reclassification(self):
        payload = {"transition_strengths": [
            transition(0, "crit", 100, (0, 0, 0), (60, 0, 0)),
            transition(0, "nucl", 60, (0, 0, 0), (90, 0, 0)),
            transition(0, "perc", 50, (0, 0, 0), (100, 0, 0)),
        ]}
        updates = ew_entry_updates(payload)
        self.assertFalse(updates["ewpt_baryo_candidate"])
        self.assertTrue(updates["ewpt_gw_candidate"])
        self.assertEqual(updates["ewpt_ew_entry_nucl_jump_over_T"], 1.5)
        self.assertEqual(updates["ewpt_ew_entry_perc_jump_over_T"], 2)
        self.assertEqual(updates["ewpt_gw_max_temperature_kind"], "perc")
        self.assertEqual(updates["ewpt_gw_max_field_jump_over_T"], 2)

    def test_singlet_only_jump_from_zero_is_gw_candidate(self):
        payload = {"transition_strengths": [
            transition(1, "crit", 100, (0, 0, 0), (0, 0, 120))
        ]}
        updates = ew_entry_updates(payload)
        self.assertFalse(updates["ewpt_baryo_candidate"])
        self.assertIsNone(updates["ewpt_ew_entry_jump_over_T"])
        self.assertTrue(updates["ewpt_gw_candidate"])
        self.assertAlmostEqual(updates["ewpt_gw_crit_field_jump_over_T"], 1.2)

    def test_no_critical_entry_but_nucleation_gw_candidate(self):
        payload = {"transition_strengths": [
            transition(1, "nucl", 50, (0, 0, 0), (0, 55, 0))
        ]}
        updates = ew_entry_updates(payload)
        self.assertFalse(updates["ewpt_baryo_candidate"])
        self.assertTrue(updates["ewpt_gw_candidate"])
        self.assertIsNone(updates["ewpt_gw_crit_field_jump_over_T"])
        self.assertEqual(updates["ewpt_gw_max_temperature_kind"], "nucl")

    def test_percolation_alone_can_trigger_gw_flag(self):
        payload = {"transition_strengths": [
            transition(1, "crit", 100, (0, 0, 0), (0, 40, 0)),
            transition(1, "perc", 50, (0, 0, 0), (0, 60, 0)),
        ]}
        updates = ew_entry_updates(payload)
        self.assertFalse(updates["ewpt_baryo_candidate"])
        self.assertTrue(updates["ewpt_gw_candidate"])
        self.assertEqual(updates["ewpt_gw_crit_field_jump_over_T"], 0.4)
        self.assertEqual(updates["ewpt_gw_perc_field_jump_over_T"], 1.2)
        self.assertEqual(updates["ewpt_gw_max_temperature_kind"], "perc")

    def test_strongest_critical_ew_entry_is_selected_even_without_completion(self):
        payload = {"transition_strengths": [
            transition(0, "crit", 200, (0, 0, 0), (50, 0, 0)),
            transition(1, "crit", 100, (0, 0, 0), (130, 0, 0)),
            transition(0, "compl", 80, (0, 0, 0), (100, 0, 0)),
        ]}
        updates = ew_entry_updates(payload)
        self.assertEqual(updates["ewpt_ew_entry_transition_index"], 1)
        self.assertTrue(updates["ewpt_baryo_candidate"])
        self.assertFalse(updates["ewpt_ew_entry_completed"])
        self.assertIsNone(
            ew_entry_updates(payload, w1_threshold=140)["ewpt_ew_entry_jump_over_T"]
        )

    def test_threshold_is_strict_and_failed_rows_are_ignored(self):
        payload = {"transition_strengths": [
            transition(0, "crit", 100, (0, 0, 0), (100, 0, 0)),
            {**transition(1, "crit", 100, (0, 0, 0), (200, 0, 0)), "status": "failed"},
        ]}
        updates = ew_entry_updates(payload)
        self.assertFalse(updates["ewpt_baryo_candidate"])
        self.assertFalse(updates["ewpt_gw_candidate"])
        self.assertEqual(updates["ewpt_ew_entry_jump_over_T"], 1)

    def test_missing_temperature_status_uses_finite_transition(self):
        critical = transition(0, "crit", 100, (0, 0, 0), (120, 0, 0))
        critical["status"] = None
        updates = ew_entry_updates({"transition_strengths": [critical]})
        self.assertTrue(updates["ewpt_baryo_candidate"])
        self.assertTrue(updates["ewpt_gw_candidate"])

    def test_no_transition_is_negative_on_a_successful_run(self):
        updates = ew_entry_updates({"transition_strengths": []})
        self.assertFalse(updates["ewpt_baryo_candidate"])
        self.assertFalse(updates["ewpt_gw_candidate"])
        self.assertIsNone(updates["ewpt_ew_entry_jump_over_T"])


if __name__ == "__main__":
    unittest.main()
