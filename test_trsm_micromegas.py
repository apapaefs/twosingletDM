import contextlib
import io
import math
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

from test_generate_trsm_points_runner import load_generator_module
from trsm_micromegas import default_micromegas_main, micromegas_configuration
from trsm_scan_campaign import configuration_fingerprint


class TestMicromegasSelection(unittest.TestCase):
    @patch.dict(os.environ)
    def test_default_and_version_aliases(self):
        os.environ.pop("TRSM_RUNTIME_ROOT", None)
        generator = load_generator_module()
        for value, version in (("6", "6.1.15"), ("7", "7.1.4"), ("7.1.4", "7.1.4")):
            args = generator.parse_args(["--micromegas-version", value])
            self.assertEqual(args.micromegas_version, version)
            self.assertEqual(
                Path(micromegas_configuration(args)["executable"]),
                Path(__file__).resolve().parents[1] / "runtime-v2" / f"micromegas_{version}" / "TRSM/main",
            )
        self.assertEqual(generator.cli_args.micromegas_version, "7.1.4")
        self.assertEqual(default_micromegas_main(), default_micromegas_main("7"))
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            generator.parse_args(["--micromegas-version", "8"])

    def test_version_separates_run_tags_and_metadata(self):
        old = load_generator_module(["123", "--micromegas-version", "6"])
        new = load_generator_module(["123", "--micromegas-version", "7"])
        self.assertNotEqual(new.RunTag,old.RunTag)
        self.assertIn("-mo6.1.15",old.RunTag)
        self.assertIn("cmb-planck2018",new.RunTag)
        metadata = new.build_scan_metadata(new.cli_args, new.RunTag, Path("scan.dat"))
        self.assertEqual(metadata["micromegas"]["version"], "7.1.4")
        self.assertEqual(metadata["micromegas"]["executable"], str(default_micromegas_main("7")))

    def test_selected_executable_reaches_dm_provider(self):
        generator = load_generator_module(["123", "--micromegas-version", "7"])
        generator.generate_lams = lambda *args, **kwargs: (
            200.0, 0.0, 380.0, 500.0, -0.15, 0.0, 0.0,
            1.0, 2.0, 0.0, 11.0, 12.0, 13.0, 123.0, 122.0,
            1111.0, 1112.0, 1113.0, 133.0, 0.99, -0.1, 0.0,
            {}, {}, {}, 0.5, 0.25, 0.0,
        )
        generator.np = math
        for name in (
            "pred", "H1", "H2", "H3", "Mz", "Mw",
            "Delta_S_central_wU", "Delta_T_central_wU", "Delta_U_central_wU",
            "errS_wU", "errT_wU", "errU_wU", "covST_wU", "covSU_wU", "covTU_wU",
        ):
            setattr(generator, name, 0)
        for name in ("check_EWPO_wU", "check_wmass_tania", "theory_constraints_vxzero", "test_evo_vxzero"):
            setattr(generator, name, lambda *args: True)
        generator.analyze_parampoint = lambda *args, **kwargs: (True, True)
        generator.test_dm = Mock(return_value=(False, {}, {}))
        generator.evaluate_trsm_point_vxzero(123, 380.0, 500.0, 200.0, -0.15, 0.1, 0.05, 0.15)
        self.assertEqual(
            generator.test_dm.call_args.kwargs["micromegas_main"],
            default_micromegas_main("7"),
        )

    def test_resume_restores_version_and_custom_executable(self):
        generator = load_generator_module()
        with tempfile.TemporaryDirectory() as tmp:
            scan = Path(tmp) / "scan.dat"
            executable = Path(tmp) / "custom/main"
            args = generator.parse_args([
                "123", "--micromegas-version", "7", "--micromegas-main", str(executable),
            ])
            before = generator.immutable_scan_configuration(args)
            generator.write_scan_metadata_file(
                scan.with_suffix(".metadata.json"),
                generator.build_scan_metadata(args, "scan", scan),
            )
            resumed = generator.parse_args(["--resume-from", str(scan)])
            self.assertEqual(resumed.micromegas_version, "7.1.4")
            self.assertEqual(resumed.micromegas_main, executable.resolve())
            self.assertEqual(before, generator.immutable_scan_configuration(resumed))
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                generator.parse_args(["--resume-from", str(scan), "--micromegas-version", "6"])

    def test_legacy_campaign_requires_explicit_reevaluation(self):
        generator=load_generator_module()
        with tempfile.TemporaryDirectory() as tmp:
            scan=Path(tmp)/"scan.dat"
            metadata=generator.build_scan_metadata(generator.cli_args,"scan",scan)
            metadata.pop("physics_version")
            generator.write_scan_metadata_file(scan.with_suffix(".metadata.json"),metadata)
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                generator.parse_args(["--resume-from",str(scan)])

    def test_missing_executable_fails_before_output_creation(self):
        generator = load_generator_module(["--micromegas-main", "/nonexistent/TRSM/main"])
        generator.reset_output = Mock()
        with self.assertRaisesRegex(generator.CampaignStateError, "not found or not executable"):
            generator.main()
        generator.reset_output.assert_not_called()

    def test_resume_rejects_a_changed_recorded_backend_before_writing(self):
        generator = load_generator_module(["--micromegas-version", "7"])
        generator.cli_args.resume_from = Path("scan.dat")
        generator.cli_args.resume_metadata = {
            "micromegas": {"version": "6.1.15", "executable": str(default_micromegas_main())},
        }
        generator.reset_output = Mock()
        with self.assertRaisesRegex(generator.CampaignStateError, "differs from the saved campaign"):
            generator.main()
        generator.reset_output.assert_not_called()


if __name__ == "__main__":
    unittest.main()
