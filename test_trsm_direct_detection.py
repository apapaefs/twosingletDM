import contextlib
import copy
import io
import json
import math
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock

from test_generate_trsm_points_runner import load_generator_module
from test_trsm_DM import DMPoint, MicromegasResult, summarize_dm_result, test_dm
from trsm_direct_detection import load_si_limit_table


def synthetic_table():
    """Artificial numbers for software tests only; not an experimental limit."""
    return {
        "schema": "trsm_si_upper_limit_v1",
        "label": "synthetic-test-only",
        "source": "Artificial software-test fixture, not LZ data",
        "confidence_level": 0.9,
        "interaction": "elastic_isoscalar_si",
        "limit_kind": "observed_upper",
        "cross_section": "per_nucleon",
        "cross_section_unit": "cm2",
        "points": [
            {"mass_GeV": 10.0, "upper_limit": 1e-45},
            {"mass_GeV": 100.0, "upper_limit": 1e-47},
            {"mass_GeV": 1000.0, "upper_limit": 1e-46},
        ],
    }


class TestSILimitTable(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / "synthetic.json"
        self.write_table(synthetic_table())

    def write_table(self, data):
        self.path.write_text(json.dumps(data), encoding="utf-8")

    def test_units_anchors_and_logarithmic_interpolation(self):
        table = load_si_limit_table(self.path)
        for mass, limit in ((10, 1e-9), (100, 1e-11), (1000, 1e-10)):
            self.assertAlmostEqual(table.upper_limit_pb(mass) / limit, 1)
        self.assertAlmostEqual(table.upper_limit_pb(math.sqrt(10 * 100)) / 1e-10, 1)
        data = synthetic_table()
        data["cross_section_unit"] = "pb"
        for point in data["points"]:
            point["upper_limit"] *= 1e36
        self.write_table(data)
        self.assertEqual(load_si_limit_table(self.path).limits_pb, table.limits_pb)

    def test_rejects_unphysical_masses_and_extrapolation(self):
        table = load_si_limit_table(self.path)
        for mass in (0, -1, math.nan, math.inf, 9.99, 1000.01):
            with self.subTest(mass=mass), self.assertRaises(ValueError):
                table.upper_limit_pb(mass)

    def test_rejects_wrong_physics_or_confidence_interval(self):
        for field, value in (
            ("schema", "hepdata"),
            ("confidence_level", 0.95),
            ("interaction", "inelastic_isoscalar_si"),
            ("limit_kind", "median_expected"),
            ("limit_kind", "observed_lower"),
            ("cross_section", "per_nucleus"),
            ("cross_section_unit", "GeV^-2"),
            ("source", ""),
        ):
            data = synthetic_table()
            data[field] = value
            self.write_table(data)
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                load_si_limit_table(self.path)

    def test_rejects_bad_or_nonmonotonic_table(self):
        original = synthetic_table()
        variants = [None, [], original["points"][:1], list(reversed(original["points"]))]
        for field, value in (
            ("mass_GeV", 10), ("mass_GeV", math.nan),
            ("upper_limit", 0), ("upper_limit", -1),
            ("upper_limit", math.inf), ("upper_limit", "1e-47"),
        ):
            points = copy.deepcopy(original["points"])
            points[1][field] = value
            variants.append(points)
        for points in variants:
            data = dict(original, points=points)
            self.write_table(data)
            with self.subTest(points=points), self.assertRaises(ValueError):
                load_si_limit_table(self.path)

    def test_rescaling_and_exclusion_boundary_use_observed_upper_limit(self):
        table = load_si_limit_table(self.path)
        point = DMPoint(0.1, 0.01, 0.01, 100, 200, -0.15, 380)
        base = table.upper_limit_pb(100)
        for fraction in (0, 0.25, 1):
            for factor in (0.999, 1.0, 1.001):
                # A finite cross section cannot exclude a zero-abundance candidate.
                cross_section = factor * base / fraction if fraction else 1
                result = MicromegasResult(100, 0.121 * fraction, cross_section)
                summary = summarize_dm_result(point, result, limit_table=table)
                self.assertEqual(summary.direct_detection_excluded, bool(fraction and factor > 1))
                self.assertEqual(summary.lux_base_limit, base)
        result = MicromegasResult(100, 0.121 / 4, 2 * base)
        self.assertFalse(summarize_dm_result(point, result, limit_table=table).direct_detection_excluded)
        self.assertTrue(summarize_dm_result(point, result, limit_table=table, rescale=False).direct_detection_excluded)

    def test_provider_records_table_identity_and_retains_gamma_constraint(self):
        from test_trsm_DM_unit import MICROMEGAS_OUTPUT_WITH_INDIRECT

        passed, info, diagnostics = test_dm(
            0.2, 0.1, 0.5, 1000, 500, math.asin(0.3), 500,
            raw_output=MICROMEGAS_OUTPUT_WITH_INDIRECT, limit_table=self.path,
        )
        self.assertFalse(passed)
        self.assertIn("Fermi-LAT gamma-line flux above limit", info)
        self.assertTrue(diagnostics["dm_indirect_detection_excluded"])
        self.assertEqual(diagnostics["dm_limit_model"], load_si_limit_table(self.path).model_id)
        self.assertAlmostEqual(diagnostics["dm_lux_base_limit"] / 1e-10, 1)
        self.assertAlmostEqual(diagnostics["dm_dir_det_limit"] / (1e-10 * 0.121 / 0.049), 1)

    def test_scan_metadata_run_name_and_legacy_fingerprint(self):
        baseline = load_generator_module(["123"])
        before = baseline.immutable_scan_configuration(baseline.cli_args)
        del baseline.cli_args.dm_limit_table
        del baseline.cli_args._dm_limit_table
        self.assertEqual(before, baseline.immutable_scan_configuration(baseline.cli_args))
        generator = load_generator_module(["123", "--dm-limit-table", str(self.path)])
        table = generator.cli_args._dm_limit_table
        self.assertIn("-dd-synthetic-test-only-" + table.sha256[:12], generator.RunTag)
        metadata = generator.build_scan_metadata(generator.cli_args, generator.RunTag, Path("scan.dat"))
        self.assertEqual(metadata["direct_detection"]["sha256"], table.sha256)
        self.assertEqual(metadata["direct_detection"]["comparison_unit"], "pb")
        self.assertEqual(metadata["direct_detection"]["limit_kind"], "observed_upper")
        self.assertNotIn("_dm_limit_table", metadata["options"])
        json.dumps(metadata)
        self.assertNotEqual(before, generator.immutable_scan_configuration(generator.cli_args))

    def test_resume_restores_table_and_rejects_changed_contents_before_output(self):
        generator = load_generator_module(["123", "--dm-limit-table", str(self.path)])
        args = generator.cli_args
        before = generator.immutable_scan_configuration(args)
        scan = self.path.with_name("scan.dat")
        generator.write_scan_metadata_file(
            scan.with_suffix(".metadata.json"),
            generator.build_scan_metadata(args, generator.RunTag, scan),
        )
        resumed = generator.parse_args(["--resume-from", str(scan)])
        self.assertEqual(before, generator.immutable_scan_configuration(resumed))
        self.assertEqual(resumed._dm_limit_table.sha256, args._dm_limit_table.sha256)
        data = synthetic_table()
        data["points"][1]["upper_limit"] *= 2
        self.write_table(data)
        # The ongoing scan keeps the original values loaded at startup.
        self.assertAlmostEqual(args._dm_limit_table.upper_limit_pb(100) / 1e-11, 1)
        generator.cli_args = generator.parse_args(["--resume-from", str(scan)])
        self.assertNotEqual(before, generator.immutable_scan_configuration(generator.cli_args))
        generator.reset_output = Mock()
        with self.assertRaisesRegex(generator.CampaignStateError, "contents differ"):
            generator.main()
        generator.reset_output.assert_not_called()

    def test_missing_or_malformed_table_fails_at_argument_parsing(self):
        generator = load_generator_module()
        for content in ("missing", "not JSON"):
            if content == "missing":
                self.path.unlink()
            else:
                self.path.write_text(content)
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                generator.parse_args(["--dm-limit-table", str(self.path)])

    def test_scan_range_must_be_covered_before_any_output(self):
        for options in (
            ["--independent-m3", "--m3-min", "4"],
            ["--m2", "380", "--m3", "5", "--vs", "200", "--a12", "-0.15",
             "--lx", "0.1", "--lphix", "0.01", "--lsx", "0.01"],
        ):
            generator = load_generator_module(["123", "--dm-limit-table", str(self.path), *options])
            generator.reset_output = Mock()
            with self.assertRaisesRegex(generator.CampaignStateError, "does not cover"):
                generator.main()
            generator.reset_output.assert_not_called()


class TestLZ2026TemporaryLimit(unittest.TestCase):
    path = Path(__file__).resolve().parent / "DM/data/lz2026/lz2026-figs7-highmass-approx.json"

    def test_paper_anchor_units_and_requested_high_mass_scaling(self):
        table = load_si_limit_table(self.path)
        self.assertEqual(table.masses_gev, (400, 1000, 4000))
        for mass, expected_pb in (
            (400, 1.824e-11), (500, 2.28e-11), (750, 3.42e-11),
            (1000, 4.56e-11), (2000, 9.12e-11), (4000, 1.824e-10),
        ):
            self.assertAlmostEqual(table.upper_limit_pb(mass) / expected_pb, 1, places=12)
        for mass in (62, 399.99, 4000.01):
            with self.assertRaisesRegex(ValueError, "outside"):
                table.upper_limit_pb(mass)

    def test_provenance_records_approximation_and_does_not_mutate_loaded_table(self):
        table = load_si_limit_table(self.path)
        metadata = table.metadata()
        provenance = metadata["provenance"]
        self.assertEqual(provenance["status"], "temporary_high_mass_approximation")
        self.assertFalse(provenance["official_LZ_mass_table"])
        self.assertEqual(provenance["published_anchor_delta_keV"], 0)
        self.assertEqual(provenance["published_anchor_mass_GeV"], 1000)
        provenance["official_LZ_mass_table"] = True
        self.assertFalse(table.metadata()["provenance"]["official_LZ_mass_table"])

    def test_inelastic_audit_table_cannot_be_used_as_elastic_mass_limits(self):
        path = self.path.with_name("figure-s7-o1-1tev-digitized.json")
        with self.assertRaisesRegex(ValueError, "schema"):
            load_si_limit_table(path)
        report = json.loads(path.read_text())
        self.assertEqual(report["DM_mass_GeV"], 1000)
        self.assertEqual([point["delta_keV"] for point in report["points"]], list(range(0, 351, 50)))
        self.assertAlmostEqual(report["points"][0]["upper_limit_cm2"] / 4.55555e-47, 1, places=6)
        self.assertLess(report["normalization_check"]["max_relative_difference"], 0.005)

    def test_scan_preserves_approximation_provenance_for_both_backends(self):
        for version in ("6", "7"):
            generator = load_generator_module([
                "123", "--micromegas-version", version, "--dm-limit-table", str(self.path),
                "--m3-min", "400", "--m3-max", "4000",
            ])
            metadata = generator.build_scan_metadata(generator.cli_args, generator.RunTag, Path("scan.dat"))
            self.assertIn("highmass-approx", generator.RunTag)
            self.assertFalse(metadata["direct_detection"]["provenance"]["official_LZ_mass_table"])
            self.assertEqual(metadata["direct_detection"]["mass_range_gev"], [400, 4000])

    def test_approximate_limit_retains_abundance_rescaling(self):
        from test_trsm_DM_unit import MICROMEGAS_OUTPUT

        _passed, _info, diagnostics = test_dm(
            0.2, 0.1, 0.5, 1000, 500, math.asin(0.3), 500,
            raw_output=MICROMEGAS_OUTPUT, limit_table=self.path,
        )
        self.assertAlmostEqual(diagnostics["dm_lux_base_limit"] / 4.56e-11, 1)
        self.assertAlmostEqual(diagnostics["dm_dir_det_limit"] / (4.56e-11 * 0.121 / 0.049), 1)
        self.assertIn("highmass-approx", diagnostics["dm_limit_model"])


if __name__ == "__main__":
    unittest.main()
