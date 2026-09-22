import contextlib
import io
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock

from trsm_cmb import (
    CMB_COLUMNS, CMB_METHOD, CMB_OUTPUT_PREFIX, CMBSignal, assess_cmb_limit,
    cmb_diagnostics, parse_cmb_signal, require_cmb_capability,
)
from test_trsm_DM import parse_micromegas_output, summarize_dm_result, test_dm
from test_trsm_DM_unit import MICROMEGAS_OUTPUT_WITH_RESCALED_INDIRECT
from test_generate_trsm_points_runner import load_generator_module
from trsm_scan_campaign import configuration_fingerprint
from scan_output import output_columns


def signal_output(ratio=1.0, rate=1e-26):
    return CMB_OUTPUT_PREFIX + json.dumps({"status": "ok", "ratio_raw": ratio, "sigma_v_cm3_s": rate}) + "\n"


def capable_driver(path, **overrides):
    capability = dict(method=CMB_METHOD, pann_limit_cm3_s_GeV=3.2e-28,
                      spectrum_key=7, spectra_flag=0, vrot_km_s=220, vz_decay=0, vw_decay=0)
    capability.update(overrides)
    payload = json.dumps({"schema": "trsm_driver_capabilities_v1", "planck_cmb": capability})
    path.write_text(f"#!{sys.executable}\nprint({payload!r})\n")
    path.chmod(0o755)
    return path


class TestCMBPhysics(unittest.TestCase):
    def test_boundary_and_abundance_rescaling(self):
        for raw, omega, fraction, ratio, excluded in (
            (1, .12, 1, 1, False), (1.000001, .12, 1, 1.000001, True),
            (4, .06, .5, 1, False), (4.1, .06, .5, 1.025, True),
            (2, .121, 1, 2, True), (1e300, 0, 0, 0, False),
            (0, .12, 1, 0, False),
        ):
            with self.subTest(raw=raw, omega=omega):
                result = assess_cmb_limit(parse_cmb_signal(signal_output(raw)), omega)
                self.assertEqual(result.fraction, fraction)
                self.assertAlmostEqual(result.ratio, ratio)
                self.assertIs(result.excluded, excluded)
                self.assertIs(result.passed, not excluded)

    def test_missing_invalid_and_failed_results_are_unavailable_even_at_zero_abundance(self):
        for output in (
            "", signal_output(-1), signal_output(math.nan), signal_output(math.inf),
            signal_output(None), signal_output(True), signal_output(0, -1),
            signal_output(0, math.nan), signal_output() * 2,
            CMB_OUTPUT_PREFIX + "broken json", CMB_OUTPUT_PREFIX + "[]",
            CMB_OUTPUT_PREFIX + '{"status":"error","reason":"calcSpectrum_error"}',
        ):
            with self.subTest(output=output):
                signal = parse_cmb_signal(output)
                self.assertFalse(signal.available)
                self.assertTrue(signal.reason)
                result = assess_cmb_limit(signal, 0)
                self.assertFalse(result.passed)
                self.assertIsNone(result.excluded)
                self.assertIsNone(cmb_diagnostics(result)["dm_cmb_ratio"])

    def test_enabled_cmb_gates_aggregate_without_losing_other_diagnostics(self):
        base = MICROMEGAS_OUTPUT_WITH_RESCALED_INDIRECT
        for extra, passed, available in ((signal_output(100), False, True),
                                         (signal_output(0), True, True), ("", False, False)):
            result = parse_micromegas_output(base + extra)
            self.assertTrue(summarize_dm_result(None, result).passed)
            summary = summarize_dm_result(None, result, planck_cmb=True)
            self.assertIs(summary.passed, passed)
            verdict, info, fields = test_dm(.2, .1, .5, 1000, 500, .3, 500,
                                           raw_output=base + extra, planck_cmb=True)
            self.assertIs(verdict, passed)
            self.assertEqual(fields["dm_omega"], .049)
            self.assertIs(fields["dm_cmb_available"], available)
            self.assertFalse(fields["dm_direct_detection_excluded"])
            self.assertIn("CMB", info)

    def test_disabled_does_not_add_columns_or_change_results(self):
        base = MICROMEGAS_OUTPUT_WITH_RESCALED_INDIRECT
        before = summarize_dm_result(None, parse_micromegas_output(base))
        after = summarize_dm_result(None, parse_micromegas_output(base + signal_output(1e6)))
        self.assertEqual(before.passed, after.passed)
        self.assertIsNone(after.cmb_limit)
        self.assertFalse(set(CMB_COLUMNS) & set(output_columns({})))
        self.assertEqual(output_columns({}, planck_cmb=True), output_columns({}) + list(CMB_COLUMNS))

    def test_core_failure_marks_cmb_unavailable(self):
        passed, _, fields = test_dm(.2, .1, .5, 1000, 500, .3, 500,
                                    raw_output="no valid DM output", planck_cmb=True)
        self.assertFalse(passed)
        self.assertTrue(fields["dm_cmb_enabled"])
        self.assertFalse(fields["dm_cmb_available"])
        self.assertEqual(fields["dm_cmb_status"], "dm_unavailable")


class TestCMBScanCompatibility(unittest.TestCase):
    def test_version_defaults_and_explicit_overrides(self):
        generator = load_generator_module()
        for options, expected in (([], False), (["--micromegas-version", "7"], True),
                                  (["--micromegas-version", "7", "--no-planck-cmb"], False),
                                  (["--micromegas-version", "6", "--planck-cmb"], True)):
            self.assertIs(generator.parse_args(options).planck_cmb, expected)

    def test_capability_query_rejects_old_or_incompatible_driver(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "main"
            self.assertEqual(require_cmb_capability(capable_driver(path))["spectrum_key"], 7)
            for overrides in ({"spectrum_key": 0}, {"method": "other"}, {"vrot_km_s": math.nan}):
                capable_driver(path, **overrides)
                with self.assertRaisesRegex(ValueError, "rebuilt TRSM driver"):
                    require_cmb_capability(path)
            path.write_text(f"#!{sys.executable}\nprint('old TRSM driver')\n")
            generator = load_generator_module(["--micromegas-version", "7", "--micromegas-main", str(path)])
            generator.reset_output = Mock()
            with self.assertRaisesRegex(generator.CampaignStateError, "rebuilt TRSM driver"):
                generator.main()
            generator.reset_output.assert_not_called()

    def test_legacy_version_seven_resume_preserves_fingerprint_and_columns(self):
        generator = load_generator_module()
        with tempfile.TemporaryDirectory() as tmp:
            scan = Path(tmp) / "legacy.dat"
            args = generator.parse_args(["123", "--micromegas-version", "7", "--no-planck-cmb"])
            del args.planck_cmb
            before = generator.immutable_scan_configuration(args)
            metadata = generator.build_scan_metadata(args, "legacy", scan)
            del metadata["planck_cmb"]
            generator.write_scan_metadata_file(scan.with_suffix(".metadata.json"), metadata)
            resumed = generator.parse_args(["--resume-from", str(scan)])
            self.assertFalse(resumed.planck_cmb)
            self.assertEqual(configuration_fingerprint(before), configuration_fingerprint(generator.immutable_scan_configuration(resumed)))
            self.assertEqual(output_columns({}, planck_cmb=resumed.planck_cmb), output_columns({}))
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                generator.parse_args(["--resume-from", str(scan), "--planck-cmb"])

    def test_cmb_method_is_part_of_new_fingerprint(self):
        generator = load_generator_module()
        args = generator.parse_args(["--micromegas-version", "7"])
        args._cmb_driver = {"spectra_flag": 0}
        original = configuration_fingerprint(generator.immutable_scan_configuration(args))
        args._cmb_driver = {"spectra_flag": 1}
        self.assertNotEqual(original, configuration_fingerprint(generator.immutable_scan_configuration(args)))
        args.planck_cmb = False
        self.assertNotEqual(original, configuration_fingerprint(generator.immutable_scan_configuration(args)))

    def test_changed_driver_settings_reject_scan_resume_before_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            driver = capable_driver(Path(tmp) / "main")
            generator = load_generator_module(["--micromegas-version", "7", "--micromegas-main", str(driver)])
            generator.cli_args._cmb_driver = require_cmb_capability(driver)
            metadata = generator.build_scan_metadata(generator.cli_args, "scan", Path(tmp) / "scan.dat")
            metadata["planck_cmb"]["driver"]["vrot_km_s"] = 230
            generator.cli_args.resume_from = Path(tmp) / "scan.dat"
            generator.cli_args.resume_metadata = metadata
            generator.reset_output = Mock()
            with self.assertRaisesRegex(generator.CampaignStateError, "CMB settings differ"):
                generator.main()
            generator.reset_output.assert_not_called()


if __name__ == "__main__":
    unittest.main()
