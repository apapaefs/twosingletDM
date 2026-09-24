import csv
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from flavour import flavour_observables as ups
from generate_trsm_info import generate_lams
from reevaluate_trsm_flavour import reevaluate
from trsm_flavour import (FLAVOUR_COLUMNS, assess_scalar, combine_flavour,
                          flavour_configuration, flavour_updates_from_row,
                          generated_flavour_updates)


def base_brs(mu=.1, tau=.8, width=1.):
    return [0., tau, mu] + [0.] * 9 + [width]


def row(m2=5, m3=20, angle=.4, **updates):
    result = dict(M2=str(m2), M3=str(m3), vx="0", vs="300", a12=str(angle),
                  lX=".1", lPhiX=".01", lSX=".001", k2="", w2="historic-width",
                  dm="True", ewpt_baryo_candidate="True")
    result.update(updates)
    return result


class TestFlavourPhysics(unittest.TestCase):
    def test_frozen_collaborator_curve_values_and_normalization(self):
        self.assertAlmostEqual(ups.babar_mumu_limit(5), 6.4958934e-7, delta=1.e-14)
        self.assertAlmostEqual(min(ups.belle_tautau_limit(5), ups.babar_tautau_limit(5)),
                               4.7304496e-6, delta=1.e-13)
        expected = (.0239 * 1.1663787e-5 * 4.18**2 /
                    (math.sqrt(2) * math.pi / 137.035999084) * .2**2 *
                    (1 - (5/9.46040)**2) * (2/3) * (1 - (5/9.46040)**6))
        self.assertAlmostEqual(ups.BRUpsgammahi(5, .2), expected, delta=1.e-18)

    def test_checker_and_adapter_agree_across_boundaries(self):
        for mass in (.211, .212, 3.549999, 3.55, 3.8, 3.8343, 4, 5, 9.2, 9.200001, 125.09):
            for coupling in (0, -.4, .8, .2 + .3j):
                with self.subTest(mass=mass, coupling=coupling):
                    expected = ups.CheckUpsLeptonBounds(mass, coupling, .1, .8)
                    result = assess_scalar("h2", mass, coupling, base_brs(), 1.)
                    self.assertEqual(result["passed"], bool(expected))
                    self.assertIs(type(expected), int)

    def test_exact_limit_passes_and_both_channels_are_recorded(self):
        with patch.object(ups, "BRUpsgammahi", return_value=1.), \
                patch.object(ups, "belle_mumu_limit", return_value=.5), \
                patch.object(ups, "babar_mumu_limit", return_value=.8), \
                patch.object(ups, "belle_tautau_limit", return_value=.5), \
                patch.object(ups, "babar_tautau_limit", return_value=.8):
            self.assertEqual(ups.CheckUpsLeptonBounds(5, 1., .5, .5), 1)
            for mu, tau, failing in ((np.nextafter(.5, 1), .5, {"mumu"}),
                                     (.5, .6, {"tautau"}), (.6, .6, {"mumu", "tautau"})):
                result = ups.EvaluateUpsLeptonBounds(5, 1., mu, tau)
                self.assertFalse(result["passed"])
                self.assertEqual(set(result["channels"]), {"mumu", "tautau"})
                self.assertEqual({k for k, c in result["channels"].items() if c["excluded"]}, failing)
                self.assertEqual(result["channels"]["mumu"]["experiment"], "Belle")

    def test_sign_and_complex_phase_invariance(self):
        for coupling in (.3, -.3, .3j):
            self.assertEqual(ups.EvaluateUpsLeptonBounds(5, coupling, .1, .8),
                             ups.EvaluateUpsLeptonBounds(5, .3, .1, .8))

    def test_coverage_uses_csv_extent_without_extrapolation(self):
        result = ups.EvaluateUpsLeptonBounds(3.8, .1, .1, .8)
        self.assertTrue(result["channels"]["mumu"]["covered"])
        self.assertFalse(result["channels"]["tautau"]["covered"])
        self.assertEqual(result["channels"]["tautau"]["status"], "outside_curve")
        self.assertIsNone(result["channels"]["tautau"]["ratio"])
        self.assertFalse(ups.EvaluateUpsLeptonBounds(9.20001, 1., 1., 1.)["channels"]["mumu"]["covered"])

    def test_physical_width_suppresses_visible_products_once(self):
        plain = assess_scalar("h2", 5, .8, base_brs(), 1.)
        invisible = assess_scalar("h2", 5, .8, base_brs(), 100.)
        self.assertFalse(plain["passed"])
        self.assertTrue(invisible["passed"])
        for channel in ("mumu", "tautau"):
            self.assertAlmostEqual(invisible["channels"][channel]["prediction"] * 100,
                                   plain["channels"][channel]["prediction"])
        tiny = assess_scalar("h2", 5, .1, base_brs(width=1.e-15), 1.e-15)
        self.assertEqual(tiny["channels"]["mumu"]["br"], .1)

    def test_zero_mixing_and_stable_h3_need_no_decay_table(self):
        zero = assess_scalar("h3", .3, 0, stable=True)
        self.assertTrue(zero["passed"])
        self.assertEqual(zero["status"], "zero_signal")
        self.assertEqual(zero["channels"]["mumu"]["prediction"], 0)
        self.assertIsNone(assess_scalar("h3", 5, .1, stable=True)["passed"])

    def test_invalid_inputs_are_unassessed_and_known_exclusion_wins(self):
        for mass, coupling, brs, width in (
            (-1, .1, base_brs(), 1), (math.nan, .1, base_brs(), 1),
            (5, math.inf, base_brs(), 1), (5, .1, None, 1),
            (5, .1, base_brs(mu=math.nan), 1), (5, .1, base_brs(mu=1.1), 1),
            (5, .1, base_brs(width=0), 1), (5, .1, base_brs(), .5),
        ):
            with self.subTest(mass=mass, coupling=coupling, width=width):
                self.assertIsNone(assess_scalar("h2", mass, coupling, brs, width)["passed"])
        failed = assess_scalar("h2", 5, 1, base_brs(), 1)
        unknown = assess_scalar("h3", 5, .1)
        result = combine_flavour([failed, unknown])
        self.assertIs(result["flavour"], False)
        self.assertIn("h3", result["flavour_reason"])
        self.assertNotIn("Infinity", result["flavour_details"])
        self.assertNotIn("NaN", result["flavour_details"])

    def test_saved_parameters_match_generation_with_open_and_closed_invisible_decay(self):
        for m2, m3 in ((4, 20), (5, 20), (5, 1), (9.2, 4), (10, 4)):
            saved = row(m2=m2, m3=m3)
            result = generate_lams(1, m2, m3, 300, 0, .4, 0, 0, False,
                                   lX=.1, lPhiX=.01, lSX=.001)
            generated = generated_flavour_updates(m2, m3, result[19:22], result[22:25], result[7:10])
            reevaluated = flavour_updates_from_row(saved)
            self.assertEqual(generated, reevaluated)

    def test_unsupported_saved_rows_do_not_pass(self):
        for saved in (row(m2=3.8), row(vx="1"), row(a12="nan"), row(lPhiX="")):
            self.assertIsNone(flavour_updates_from_row(saved)["flavour"])
        self.assertTrue(flavour_updates_from_row(dict(M2="200", M3="20", vx="0"))["flavour"])
        self.assertTrue(flavour_updates_from_row(row(m2=.3, angle=0))["flavour"])

    def test_preflight_checks_data_and_records_all_fingerprints(self):
        config = flavour_configuration()
        self.assertEqual(len(config["curve_coverage_GeV"]), 4)
        for name in ups.LIMIT_FILES:
            self.assertEqual(len(config["sources_sha256"]["flavour/" + name]), 64)
        self.assertIn("trsm_flavour.py", config["sources_sha256"])
        with tempfile.TemporaryDirectory() as tmp, patch.object(ups, "LIMIT_DATA_DIR", Path(tmp)):
            (Path(tmp) / ups.LIMIT_FILES[0]).write_text("mass,limit\n4,0\n5,1e-6\n")
            with self.assertRaisesRegex(ValueError, "Invalid flavour"):
                flavour_configuration()
        ups._load_limit_curve.cache_clear()


class TestFlavourReevaluation(unittest.TestCase):
    def test_changed_or_legacy_flavour_prescription_rejects_resume_before_output(self):
        from test_generate_trsm_points_runner import load_generator_module
        from trsm_scan_campaign import CampaignStateError
        generator = load_generator_module()
        generator.cli_args.resume_from = Path("unused-saved-scan.tsv")
        configuration = flavour_configuration()
        for saved in ({}, {"flavour": {**configuration, "method": "other"}}):
            generator.cli_args.resume_metadata = saved
            with patch.object(generator, "write_valid_point_file") as writer:
                with self.assertRaisesRegex(CampaignStateError, "Flavour prescription differs"):
                    generator.main()
                writer.assert_not_called()

    def test_separate_output_preserves_other_columns_and_records_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, target = Path(tmp) / "old.tsv", Path(tmp) / "new.tsv"
            original_rows = [row(), row(m2=3.8), row(m2=200)]
            with source.open("w") as stream:
                writer = csv.DictWriter(stream, original_rows[0], delimiter="\t")
                writer.writeheader()
                writer.writerows(original_rows)
            original_bytes = source.read_bytes()
            counts = reevaluate(source, target)
            self.assertEqual(sum(counts.values()), 3)
            self.assertEqual(source.read_bytes(), original_bytes)
            with target.open() as stream:
                result = list(csv.DictReader(stream, delimiter="\t"))
            for before, after in zip(original_rows, result):
                self.assertEqual(before, {key: after[key] for key in before})
                self.assertTrue(set(FLAVOUR_COLUMNS) <= set(after))
            self.assertEqual(result[1]["flavour"], "nan")
            metadata = json.loads(target.with_suffix(".metadata.json").read_text())
            self.assertEqual(metadata["flavour_reevaluation"]["rows"], 3)
            self.assertEqual(metadata["schema"], "trsm_scan_metadata_v1")
            with self.assertRaises(FileExistsError):
                reevaluate(source, target)
            with self.assertRaises(ValueError):
                reevaluate(source, source)

    def test_malformed_tsv_does_not_publish_partial_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, target = Path(tmp) / "old.tsv", Path(tmp) / "new.tsv"
            source.write_text("M2\tM3\tvx\n200\t20\t0\n300\t0\n")
            with self.assertRaisesRegex(ValueError, "malformed"):
                reevaluate(source, target)
            self.assertFalse(target.exists())
            self.assertFalse(target.with_suffix(".metadata.json").exists())


if __name__ == "__main__":
    unittest.main()
