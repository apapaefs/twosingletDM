import csv
import importlib.util
import math
import json
import sys
import tempfile
import unittest
from pathlib import Path
from trsm_cmb import CMB_COLUMNS, CMBSignal, assess_cmb_limit, cmb_diagnostics
from trsm_direct_detection import load_si_limit_table


SCRIPT_PATH = Path(__file__).resolve().parent / "reevaluate_trsm_dm_higgs.py"


def load_module():
    spec = importlib.util.spec_from_file_location(
        "reevaluate_trsm_dm_higgs", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


INPUT_HEADER = [
    "M2",
    "M3",
    "vs",
    "vx",
    "a12",
    "lX",
    "lPhiX",
    "lSX",
    "w1",
    "hb",
    "dm_omega",
    "unknown_column",
]


def write_input(path, count=3):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(INPUT_HEADER)
        for index in range(count):
            writer.writerow(
                [
                    300 + index,
                    50 + index,
                    500,
                    0,
                    0.1,
                    0.2,
                    0.3,
                    0.4,
                    999,
                    False,
                    999,
                    f"keep-{index}",
                ]
            )


def valid_updates(module, point_index):
    updates = {
        "K133": 1.0 + point_index,
        "K233": 2.0 + point_index,
        "w1": 0.005,
        "w2": 1.0 + point_index,
        "w3": 0.0,
        "hb": point_index % 2 == 1,
        "hs": True,
        "dm": False,
        "dm_mdm": 50.0 + point_index,
        "dm_omega": 0.2,
        "dm_xf": 25.0,
        "dm_freezeout_temperature_GeV": (50.0 + point_index) / 25.0,
        "dm_relic_upper_limit": 0.121,
        "dm_dir_det": 2.0e-9,
        "dm_dir_det_limit": 1.0e-9,
        "dm_lux_base_limit": 1.0e-11,
        "dm_relic_excluded": True,
        "dm_direct_detection_excluded": True,
        "dm_indirect_available": False,
        "dm_indirect_channels_seen": 0,
        "dm_indirect_channels_used": 0,
        "dm_indirect_energy": math.nan,
        "dm_indirect_flux": math.nan,
        "dm_indirect_limit": math.nan,
        "dm_indirect_ratio": 0.0,
        "dm_indirect_detection_excluded": False,
        "dm_limit_model": "lz2025-source",
        "dm_rescale": True,
        "higgstools_hb_selected_limits": "{}",
        "higgstools_hb_top_obs": "[]",
        "higgstools_hs_chi2": 150.0,
        "higgstools_hs_delta_chi2": 1.0,
        "higgstools_hs_top_chi2": "[]",
        "h1_h3h3_width": 0.001,
        "h1_h3h3_br": 0.2,
        "h1_h2h2_width": 0.0,
        "h1_h2h2_br": 0.0,
        "h2_h3h3_width": 0.1,
        "h2_h3h3_br": 0.1 / (1.0 + point_index),
        "higgs_invisible_widths_included": True,
        "portal_convention": module.EXPECTED_CONVENTION_ID,
        "micromegas_model_convention": module.EXPECTED_CONVENTION_ID,
    }
    return updates


class RecordingEvaluator:
    def __init__(self, module, fail_at=None):
        self.module = module
        self.fail_at = fail_at
        self.calls = []

    def __call__(self, row, point_index):
        self.calls.append(point_index)
        if point_index == self.fail_at:
            raise RuntimeError("injected provider failure")
        return valid_updates(self.module, point_index)


class TestReevaluateTRSMDMHiggs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = load_module()

    def test_cmb_rescaling_and_unavailable_aggregate_validation(self):
        updates = valid_updates(self.module, 1)
        updates.update(dm_omega=.06, dm_relic_excluded=False,
                       dm_direct_detection_excluded=False)
        updates.update(cmb_diagnostics(assess_cmb_limit(CMBSignal(True, 8, "ok", ""), .06)))
        self.module.validate_updates(updates, 2)
        updates["dm"] = True
        with self.assertRaisesRegex(self.module.ReEvaluationError, "aggregate DM"):
            self.module.validate_updates(updates, 2)
        updates["dm"] = False
        updates.update(cmb_diagnostics(assess_cmb_limit(CMBSignal(), .06)))
        self.module.validate_updates(updates, 2)
        updates.update(dm_omega=0, dm_dir_det_limit=math.inf, dm=True)
        updates.update(cmb_diagnostics(assess_cmb_limit(CMBSignal(True, 8, "ok", ""), 0)))
        self.module.validate_updates(updates, 2)
        updates["dm_cmb_abundance_fraction"] = .5
        with self.assertRaises(self.module.ReEvaluationError):
            self.module.validate_updates(updates, 2)

    def test_disabled_reevaluation_clears_old_cmb_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, output = Path(tmp) / "old.dat", Path(tmp) / "new.dat"
            write_input(source, count=1)
            with source.open() as stream:
                rows = list(csv.reader(stream, delimiter="\t"))
            stale = cmb_diagnostics(assess_cmb_limit(CMBSignal(True, 100, "ok", ""), .12))
            with source.open("w") as stream:
                writer = csv.writer(stream, delimiter="\t")
                writer.writerow(rows[0] + list(CMB_COLUMNS))
                writer.writerow(rows[1] + list(stale.values()))
            original = source.read_bytes()
            self.module.reevaluate(source, output, RecordingEvaluator(self.module))
            self.assertEqual(source.read_bytes(), original)
            with output.open() as stream:
                row = next(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual(row["dm_cmb_status"], "disabled")
            self.assertEqual(row["dm_cmb_enabled"], "False")
            self.assertNotEqual(row["dm_cmb_ratio_raw"], "100")

    def test_new_checkpoint_rejects_changed_physics_configuration(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, output = Path(tmp) / "old.dat", Path(tmp) / "new.dat"
            write_input(source, count=2)
            config = {"micromegas": {"version": "7.1.4", "executable": "/test/main"},
                      "direct_detection": {"model": "lz2025-source"},
                      "planck_cmb": {"enabled": True, "method": "test"}}
            def evaluator(row, index):
                if index == 2:
                    raise RuntimeError("stop for resume")
                updates = valid_updates(self.module, index)
                updates.update(cmb_diagnostics(assess_cmb_limit(CMBSignal(True, 0, "ok", ""), .2)))
                return updates
            with self.assertRaisesRegex(self.module.ReEvaluationError, "stop for resume"):
                self.module.reevaluate(source, output, evaluator, checkpoint_every=1, evaluation_configuration=config)
            for key, change in (("planck_cmb", {"enabled": False}),
                                ("direct_detection", {"model": "legacy-output"}),
                                ("micromegas", {"version": "6.1.15", "executable": "/test/main"})):
                with self.assertRaises(self.module.ReEvaluationError):
                    self.module.reevaluate(source, output, evaluator, resume=True,
                                           evaluation_configuration={**config, key: change})
            def resumed(row, index):
                self.assertEqual(index, 2)
                updates = valid_updates(self.module, index)
                updates.update(cmb_diagnostics(assess_cmb_limit(CMBSignal(True, 0, "ok", ""), .2)))
                return updates
            self.module.reevaluate(source, output, resumed, resume=True, evaluation_configuration=config)
            self.assertEqual(json.loads(output.with_suffix(".metadata.json").read_text())["evaluation_configuration"], config)

    def test_defaults_to_v2_si_table_and_preserves_explicit_replacement(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, output = Path(tmp) / "old.dat", Path(tmp) / "new.dat"
            write_input(source, count=1)
            table_path = Path(tmp) / "limits.json"
            reference = Path(__file__).parent / "DM/data/lz2026/lz2026-figs7-highmass-approx.json"
            table_path.write_bytes(reference.read_bytes())
            table = load_si_limit_table(table_path)
            source.write_text(source.read_text().replace("300\t50\t", "300\t500\t"))
            sidecar = source.with_suffix(".metadata.json")
            sidecar.write_text(json.dumps({"direct_detection": table.metadata()}))
            from test_trsm_cmb import capable_driver
            driver=capable_driver(Path(tmp)/"main")
            args = self.module.parse_args([str(source), "--output", str(output), "--micromegas-main", str(driver)])
            config, recovered, _ = self.module.resolve_evaluation_configuration(args)
            self.assertNotEqual(recovered.sha256, table.sha256)
            self.assertTrue(recovered.model_id.startswith("table:lz-ws2024-observed-v2:"))
            self.assertTrue(config["planck_cmb"]["enabled"])
            table_path.unlink()
            self.module.resolve_evaluation_configuration(args)
            args.dm_limit_table = reference
            config, recovered, _ = self.module.resolve_evaluation_configuration(args)
            self.assertEqual(recovered.sha256, table.sha256)

    def test_legacy_checkpoint_adoption_freezes_configuration(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, output = Path(tmp) / "old.dat", Path(tmp) / "new.dat"
            write_input(source, count=2)
            evaluator = RecordingEvaluator(self.module, fail_at=2)
            with self.assertRaises(self.module.ReEvaluationError):
                self.module.reevaluate(source, output, evaluator, checkpoint_every=1)
            config = {"micromegas": {"version": "6.1.15", "executable": "/test/main"},
                      "direct_detection": {"model": "lz2025-source"},
                      "planck_cmb": {"enabled": False}}
            with self.assertRaises(self.module.ReEvaluationError):
                self.module.reevaluate(source, output, evaluator, resume=True, evaluation_configuration=config)
            saved = self.module.read_checkpoint_metadata(self.module.checkpoint_path(output))
            self.assertNotIn("evaluation_configuration", saved)
            changed = {**config, "micromegas": {"version": "7.1.4", "executable": "/test/main"}}
            with self.assertRaisesRegex(self.module.ReEvaluationError, "configuration"):
                self.module.reevaluate(source, output, evaluator, resume=True, evaluation_configuration=changed)

    def test_reevaluation_defaults_overrides_and_outdated_driver_preflight(self):
        from test_trsm_cmb import capable_driver
        with tempfile.TemporaryDirectory() as tmp:
            source, output = Path(tmp) / "old.dat", Path(tmp) / "new.dat"
            write_input(source, count=1)
            driver = capable_driver(Path(tmp) / "main")
            for options, enabled in (([], True), (["--micromegas-version", "7"], True),
                                     (["--micromegas-version", "7", "--no-planck-cmb"], False),
                                     (["--planck-cmb"], True)):
                args = self.module.parse_args([str(source), "--output", str(output), "--micromegas-main", str(driver), *options])
                config, _, _ = self.module.resolve_evaluation_configuration(args)
                self.assertIs(config["planck_cmb"]["enabled"], enabled)
            driver.write_text(f"#!{sys.executable}\nprint('old driver')\n")
            with self.assertRaisesRegex(ValueError, "rebuilt TRSM driver"):
                self.module.run([str(source), "--output", str(output), "--micromegas-version", "7", "--micromegas-main", str(driver)])
            self.assertFalse(output.exists())
            self.assertFalse(self.module.checkpoint_path(output).exists())

    def test_stale_or_incomplete_cmb_results_cannot_validate(self):
        updates = valid_updates(self.module, 1)
        updates["dm_cmb_enabled"] = True
        with self.assertRaisesRegex(self.module.ReEvaluationError, "incomplete CMB"):
            self.module.validate_updates(updates, 2)
        updates.update(cmb_diagnostics())
        updates["dm_cmb_ratio_raw"] = 99
        with self.assertRaisesRegex(self.module.ReEvaluationError, "unavailable CMB"):
            self.module.validate_updates(updates, 2)

    def test_preserves_unknown_columns_and_replaces_recomputed_values(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            source = tmpdir / "legacy.dat"
            output = tmpdir / "canonical.dat"
            write_input(source, count=2)
            evaluator = RecordingEvaluator(self.module)

            result = self.module.reevaluate(
                source, output, evaluator, checkpoint_every=1
            )

            self.assertEqual(result, output.resolve())
            self.assertEqual(evaluator.calls, [1, 2])
            self.assertFalse(self.module.checkpoint_path(output).exists())
            with output.open(encoding="utf-8", newline="") as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))

            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["unknown_column"], "keep-0")
            self.assertNotEqual(rows[0]["w1"], "999")
            self.assertEqual(rows[0]["dm_omega"], "0.2")
            self.assertEqual(rows[0]["higgs_invisible_widths_included"], "True")
            self.assertEqual(
                rows[0]["portal_convention"], self.module.EXPECTED_CONVENTION_ID
            )
            self.assertEqual(len(rows[0]), len(set(rows[0])))

    def test_dm_reevaluation_clears_stale_freezeout_phase_flags(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "with_ewpt.dat"
            output = root / "updated.dat"
            write_input(source, count=1)
            with source.open(encoding="utf-8", newline="") as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            rows[0].update({
                "ewpt_x_broken_min_T_GeV": "51.0",
                "ewpt_x_phase_at_freezeout": "broken",
                "ewpt_x_broken_at_or_after_freezeout": "True",
                "dm_relic_z2_freezeout_compatible": "False",
            })
            rows[0].update({column: "stale" for column in self.module.THERMAL_VEV_COLUMNS})
            with source.open("w", encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=rows[0], delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)

            self.module.reevaluate(source, output, RecordingEvaluator(self.module))
            with output.open(encoding="utf-8", newline="") as stream:
                updated = next(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual(updated["ewpt_x_broken_min_T_GeV"], "51.0")
            self.assertEqual(updated["dm_resonance_nearest_mediator"], "h1")
            for column in self.module.FREEZEOUT_DEPENDENT_COLUMNS:
                self.assertEqual(updated[column], "nan")

    def test_failure_leaves_no_final_output_and_resume_skips_committed_rows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            source = tmpdir / "legacy.dat"
            output = tmpdir / "canonical.dat"
            write_input(source, count=4)
            first = RecordingEvaluator(self.module, fail_at=3)

            with self.assertRaisesRegex(
                self.module.ReEvaluationError, "injected provider failure"
            ):
                self.module.reevaluate(
                    source, output, first, checkpoint_every=2
                )

            self.assertEqual(first.calls, [1, 2, 3])
            self.assertFalse(output.exists())
            self.assertTrue(self.module.checkpoint_path(output).exists())

            second = RecordingEvaluator(self.module)
            self.module.reevaluate(
                source,
                output,
                second,
                resume=True,
                checkpoint_every=2,
            )
            self.assertEqual(second.calls, [3, 4])
            with output.open(encoding="utf-8", newline="") as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual(len(rows), 4)
            self.assertEqual([row["unknown_column"] for row in rows], [
                "keep-0",
                "keep-1",
                "keep-2",
                "keep-3",
            ])

    def test_resume_rejects_changed_input_identity(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            source = tmpdir / "legacy.dat"
            output = tmpdir / "canonical.dat"
            write_input(source, count=2)
            first = RecordingEvaluator(self.module, fail_at=2)
            with self.assertRaises(self.module.ReEvaluationError):
                self.module.reevaluate(
                    source, output, first, checkpoint_every=1
                )

            with source.open("a", encoding="utf-8") as stream:
                stream.write("\t".join(["301", "51", "500", "0", "0.1", "0.2", "0.3", "0.4", "999", "False", "999", "changed"]) + "\n")

            with self.assertRaisesRegex(
                self.module.ReEvaluationError, "checkpoint does not match"
            ):
                self.module.reevaluate(
                    source,
                    output,
                    RecordingEvaluator(self.module),
                    resume=True,
                    checkpoint_every=1,
                )
            self.assertFalse(output.exists())

    def test_refuses_overwrite_and_requires_explicit_resume(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            source = tmpdir / "legacy.dat"
            output = tmpdir / "canonical.dat"
            write_input(source, count=1)
            output.write_text("do not replace\n", encoding="utf-8")
            with self.assertRaisesRegex(FileExistsError, "refusing to overwrite"):
                self.module.reevaluate(
                    source, output, RecordingEvaluator(self.module)
                )

            output.unlink()
            partial = self.module.checkpoint_path(output)
            partial.write_bytes(b"existing checkpoint")
            with self.assertRaisesRegex(FileExistsError, "use --resume"):
                self.module.reevaluate(
                    source, output, RecordingEvaluator(self.module)
                )

    def test_validation_rejects_unusable_core_results(self):
        updates = valid_updates(self.module, 1)
        updates["dm_omega"] = None
        with self.assertRaisesRegex(
            self.module.ReEvaluationError, "dm_omega.*not finite"
        ):
            self.module.validate_updates(updates, 2)

        updates = valid_updates(self.module, 1)
        updates["dm_indirect_available"] = True
        with self.assertRaisesRegex(
            self.module.ReEvaluationError, "available indirect result"
        ):
            self.module.validate_updates(updates, 2)

    def test_unavailable_dm_is_recorded_as_failed_without_losing_higgs_results(self):
        updates = valid_updates(self.module, 1)
        for column in self.module.DM_COLUMNS:
            updates[column] = None
        self.module.validate_updates(updates, 2)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            source = tmpdir / "legacy.dat"
            output = tmpdir / "canonical.dat"
            write_input(source, count=1)
            self.module.reevaluate(
                source, output, lambda _row, _index: updates, checkpoint_every=1
            )
            with output.open(encoding="utf-8", newline="") as stream:
                row = next(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual(row["dm"], "False")
            self.assertEqual(row["dm_omega"], "nan")
            self.assertEqual(row["dm_relic_excluded"], "nan")
            self.assertEqual(row["hs"], "True")
            self.assertEqual(row["higgstools_hs_chi2"], "150.0")

        updates["dm"] = True
        with self.assertRaisesRegex(
            self.module.ReEvaluationError, "unavailable DM details but is marked DM-passing"
        ):
            self.module.validate_updates(updates, 2)

    def test_cli_requires_versioned_output_and_positive_checkpoint(self):
        args = self.module.parse_args(
            ["legacy.dat", "--output", "canonical.dat", "--resume"]
        )
        self.assertEqual(args.output, Path("canonical.dat"))
        self.assertTrue(args.resume)
        self.assertEqual(args.checkpoint_every, 25)

        with self.assertRaises(SystemExit):
            self.module.parse_args(
                [
                    "legacy.dat",
                    "--output",
                    "canonical.dat",
                    "--checkpoint-every",
                    "0",
                ]
            )


if __name__ == "__main__":
    unittest.main()
