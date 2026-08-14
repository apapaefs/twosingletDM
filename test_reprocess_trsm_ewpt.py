import csv
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import reprocess_trsm_ewpt as reprocessor


class FakePoint:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class FakeConfig:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class FakeEWPTModule:
    TRSMEWPTPoint = FakePoint
    EWPTConfig = FakeConfig

    def __init__(self):
        self.calls = []
        self.fail_indexes = set()
        self.interrupt_indexes = set()
        self.eq418_failed_indexes = set()

    def check_eq_4_18(self, row):
        index = int(round(float(row["m2"])))
        return SimpleNamespace(satisfied=index not in self.eq418_failed_indexes)

    def run_trsm_ewpt(self, point, *, config, workdir, keep_files):
        self.calls.append(
            {
                "index": point.index,
                "dm_independent_point": (point.m2, point.m3),
                "config": config,
                "workdir": Path(workdir),
                "keep_files": keep_files,
            }
        )
        if point.index in self.interrupt_indexes:
            raise KeyboardInterrupt("simulated interruption")
        if point.index in self.fail_indexes:
            raise RuntimeError(f"simulated BSMPT failure for point {point.index}")
        return SimpleNamespace(point=point)

    @staticmethod
    def summarize_result(result):
        return f"summary for point {result.point.index}"

    @staticmethod
    def result_to_json(result):
        index = result.point.index
        labels = ["SYM", "X_BROKEN", "EW"] if index == 2 else ["SYM", "EW"]
        return {
            "transition_strengths": [
                {
                    "temperature_kind": "crit",
                    "ew_true_over_T": 10.0 + index,
                },
                {
                    "temperature_kind": "nucl",
                    "ew_true_over_T": 1.0 + index / 10.0,
                },
            ],
            "minimatracer": {
                "global_phase_path": labels,
                "ew_step_index": index - 1,
            },
        }


BASE_HEADER = [
    *reprocessor.REQUIRED_INPUT_COLUMNS,
    *reprocessor.EWPT_COLUMNS,
    "mg5_xsec_gg_heta0_pb",
]


def scan_row(index, *, dm=True, **updates):
    row = {
        "M2": str(index),
        "M3": "5.0",
        "vs": "100.0",
        "vx": "0.0",
        "a12": "0.1",
        "lX": "0.2",
        "lPhiX": "0.01",
        "lSX": "0.02",
        "evo": "True",
        "thc": "True",
        "hb": "True",
        "hs": "True",
        "ewpo": "True",
        "wmass": "True",
        "dm": ("True" if dm else "False") if type(dm) is bool else str(dm),
        "ewpt_ew_true_over_T": "nan",
        "ewpt_global_phase_path": "nan",
        "ewpt_has_x_broken": "nan",
        "ewpt_ew_step_index": "nan",
        "ewpt_status": "nan",
        "ewpt_error": "nan",
        "mg5_xsec_gg_heta0_pb": str(0.01 * index),
    }
    row.update({key: str(value) for key, value in updates.items()})
    return row


def write_scan(path, rows, *, metadata=False):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=BASE_HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if metadata:
        path.with_suffix(".metadata.json").write_text(
            json.dumps(
                {
                    "schema": "trsm_scan_metadata_v1",
                    "scan_file": path.name,
                    "postprocessing_history": [],
                }
            ),
            encoding="utf-8",
        )


def read_scan(path):
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


class TestReprocessTRSMEWPT(unittest.TestCase):
    def test_non_dm_selection_preserves_rows_and_skips_existing_results(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "scan.dat"
            output = root / "scan_ewpt.dat"
            workdir = root / "ewpt"
            rows = [
                scan_row(1, dm=True),
                scan_row(2, dm=False),
                scan_row(3, dm=False, hb="False"),
                scan_row(
                    4,
                    dm=True,
                    ewpt_ew_true_over_T="2.4",
                    ewpt_status="success",
                ),
            ]
            write_scan(source, rows, metadata=True)
            fake = FakeEWPTModule()

            result = reprocessor.run(
                [
                    str(source),
                    "--output",
                    str(output),
                    "--ewpt-workdir",
                    str(workdir),
                    "--ewpt-thigh",
                    "1000",
                ],
                ewpt_module=fake,
            )

            self.assertEqual([call["index"] for call in fake.calls], [1, 2])
            self.assertEqual(fake.calls[0]["config"].thigh, 1000.0)
            written = read_scan(output)
            self.assertEqual(len(written), 4)
            self.assertEqual(written[0]["ewpt_status"], "success")
            self.assertEqual(float(written[0]["ewpt_ew_true_over_T"]), 1.1)
            self.assertEqual(written[1]["ewpt_global_phase_path"], "SYM -> X_BROKEN -> EW")
            self.assertEqual(written[1]["ewpt_has_x_broken"], "True")
            self.assertEqual(written[1]["dm"], "False")
            self.assertEqual(written[2]["ewpt_status"], "nan")
            self.assertEqual(written[3]["ewpt_ew_true_over_T"], "2.4")
            self.assertEqual(written[1]["mg5_xsec_gg_heta0_pb"], "0.02")
            self.assertEqual(result.counts["success"], 2)
            self.assertEqual(result.counts["existing"], 1)
            self.assertEqual(result.counts["ineligible"], 1)
            self.assertEqual(result.counts["eligible_dm_pass"], 2)
            self.assertEqual(result.counts["eligible_dm_fail"], 1)
            self.assertFalse(reprocessor.checkpoint_path(output).exists())
            self.assertTrue(result.provenance.is_file())
            metadata = json.loads(result.metadata.read_text(encoding="utf-8"))
            self.assertEqual(metadata["scan_file"], output.name)
            self.assertEqual(len(metadata["postprocessing_history"]), 1)

    def test_bsmpt_failure_is_recorded_and_does_not_stop_later_rows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "scan.dat"
            output = root / "scan_ewpt.dat"
            write_scan(source, [scan_row(1), scan_row(2, dm=False), scan_row(3)])
            fake = FakeEWPTModule()
            fake.fail_indexes.add(2)

            result = reprocessor.run(
                [str(source), "--output", str(output)],
                ewpt_module=fake,
            )

            written = read_scan(output)
            self.assertEqual([call["index"] for call in fake.calls], [1, 2, 3])
            self.assertEqual(written[1]["ewpt_status"], "failed")
            self.assertIn("simulated BSMPT failure", written[1]["ewpt_error"])
            self.assertEqual(written[1]["ewpt_ew_true_over_T"], "nan")
            self.assertEqual(result.counts["success"], 2)
            self.assertEqual(result.counts["failed"], 1)

    def test_rerun_existing_ewpt_replaces_the_stored_result(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "scan.dat"
            output = root / "scan_ewpt.dat"
            write_scan(
                source,
                [
                    scan_row(
                        1,
                        ewpt_ew_true_over_T="9.9",
                        ewpt_status="success",
                    )
                ],
            )
            fake = FakeEWPTModule()

            result = reprocessor.run(
                [
                    str(source),
                    "--output",
                    str(output),
                    "--rerun-existing-ewpt",
                ],
                ewpt_module=fake,
            )

            self.assertEqual([call["index"] for call in fake.calls], [1])
            self.assertEqual(read_scan(output)[0]["ewpt_ew_true_over_T"], "1.1")
            self.assertEqual(result.counts["success"], 1)
            self.assertEqual(result.counts["existing"], 0)

    def test_eq418_prefilter_skips_only_the_failing_eligible_row(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "scan.dat"
            output = root / "scan_ewpt.dat"
            write_scan(source, [scan_row(1), scan_row(2, dm=False)])
            fake = FakeEWPTModule()
            fake.eq418_failed_indexes.add(1)

            result = reprocessor.run(
                [
                    str(source),
                    "--output",
                    str(output),
                    "--ewpt-require-eq418",
                ],
                ewpt_module=fake,
            )

            self.assertEqual([call["index"] for call in fake.calls], [2])
            self.assertEqual(result.counts["skipped_eq418"], 1)
            self.assertEqual(read_scan(output)[0]["ewpt_status"], "nan")

    def test_resume_continues_after_the_last_transactional_checkpoint(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "scan.dat"
            output = root / "scan_ewpt.dat"
            workdir = root / "work"
            write_scan(source, [scan_row(1), scan_row(2), scan_row(3)])
            fake = FakeEWPTModule()
            fake.interrupt_indexes.add(2)
            argv = [
                str(source),
                "--output",
                str(output),
                "--ewpt-workdir",
                str(workdir),
            ]

            with self.assertRaises(KeyboardInterrupt):
                reprocessor.run(argv, ewpt_module=fake)

            self.assertFalse(output.exists())
            self.assertTrue(reprocessor.checkpoint_path(output).is_file())
            fake.interrupt_indexes.clear()
            result = reprocessor.run([*argv, "--resume"], ewpt_module=fake)

            self.assertEqual([call["index"] for call in fake.calls], [1, 2, 2, 3])
            self.assertEqual(result.counts["success"], 3)
            self.assertTrue(output.is_file())
            self.assertFalse(reprocessor.checkpoint_path(output).exists())

    def test_invalid_stored_constraint_is_not_silently_treated_as_passing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "scan.dat"
            output = root / "scan_ewpt.dat"
            write_scan(source, [scan_row(1, dm="nan")])

            with self.assertRaises(reprocessor.EWPTReprocessingError):
                reprocessor.run(
                    [str(source), "--output", str(output)],
                    ewpt_module=FakeEWPTModule(),
                )

    def test_invalid_source_metadata_is_rejected_before_bsmpt_runs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "scan.dat"
            output = root / "scan_ewpt.dat"
            write_scan(source, [scan_row(1)])
            source.with_suffix(".metadata.json").write_text(
                json.dumps({"postprocessing_history": "invalid"}),
                encoding="utf-8",
            )
            fake = FakeEWPTModule()

            with self.assertRaises(reprocessor.EWPTReprocessingError):
                reprocessor.run(
                    [str(source), "--output", str(output)],
                    ewpt_module=fake,
                )

            self.assertEqual(fake.calls, [])
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
