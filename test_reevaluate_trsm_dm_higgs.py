import csv
import importlib.util
import math
import sys
import tempfile
import unittest
from pathlib import Path


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
