import json
import math
import os
import random
import socket
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from test_generate_trsm_points_runner import load_generator_module
from trsm_scan_campaign import (
    CampaignLock,
    CampaignStateError,
    advance_output_records,
    atomic_write_json,
    decode_rng_state,
    encode_rng_state,
    inspect_tsv,
    quarantine_uncheckpointed_point_dirs,
    reconcile_outputs,
    repair_incomplete_tsv_tail,
)


def round_sig(value, digits=2):
    if value == 0:
        return 0.0
    return round(
        value,
        digits - int(math.floor(math.log10(abs(value)))) - 1,
    )


def configure_numeric_stubs(generator):
    generator.np = SimpleNamespace(arccos=math.acos, sqrt=math.sqrt)
    generator.round_sig = round_sig


def scan_args(generator, seed, target):
    return generator.parse_args(
        [
            str(seed),
            "--nrandom",
            str(target),
            "--nrandom-count-evo-thc",
            "--write-evo-thc-points",
            "--scan-k133-k233-log",
            "--independent-m3",
        ]
    )


class TestCampaignStateUtilities(unittest.TestCase):
    def test_rng_state_has_safe_json_roundtrip(self):
        rng = random.Random(879)
        for _ in range(17):
            rng.random()
        encoded = encode_rng_state(rng.getstate())
        json.dumps(encoded)
        restored = random.Random()
        restored.setstate(decode_rng_state(encoded))
        self.assertEqual(
            [rng.random() for _ in range(10)],
            [restored.random() for _ in range(10)],
        )

    def test_output_ahead_of_checkpoint_is_backed_up_and_truncated(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "scan.dat"
            committed = b"x\n1\n"
            path.write_bytes(committed)
            expected = {"main": inspect_tsv(path)}
            path.write_bytes(committed + b"2\n")

            backups = reconcile_outputs({"main": path}, expected)

            self.assertEqual(path.read_bytes(), committed)
            self.assertEqual(len(backups), 1)
            self.assertEqual(backups[0].read_bytes(), committed + b"2\n")

    def test_output_shorter_than_checkpoint_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "scan.dat"
            path.write_text("x\n1\n2\n", encoding="ascii")
            expected = {"main": inspect_tsv(path)}
            path.write_text("x\n1\n", encoding="ascii")
            with self.assertRaises(CampaignStateError):
                reconcile_outputs({"main": path}, expected)

    def test_same_size_output_corruption_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "scan.dat"
            path.write_bytes(b"x\n1\n")
            expected = {"main": inspect_tsv(path)}
            path.write_bytes(b"x\n2\n")
            with self.assertRaises(CampaignStateError):
                reconcile_outputs({"main": path}, expected)

    def test_incomplete_legacy_tail_is_backed_up_and_removed(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "scan.dat"
            path.write_bytes(b"x\n1\npartial")

            backup = repair_incomplete_tsv_tail(path)

            self.assertIsNotNone(backup)
            self.assertEqual(path.read_bytes(), b"x\n1\n")
            self.assertEqual(backup.read_bytes(), b"x\n1\npartial")

    def test_output_records_advance_without_recounting_prefix(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "scan.dat"
            path.write_bytes(b"")
            previous = {"main": inspect_tsv(path)}
            path.write_bytes(b"x\n1\n")
            current = advance_output_records({"main": path}, previous)
            self.assertEqual(current["main"]["data_rows"], 1)
            path.write_bytes(b"x\n1\n2\n")
            current = advance_output_records({"main": path}, current)
            self.assertEqual(current["main"]["data_rows"], 2)

    def test_uncheckpointed_ewpt_directories_are_quarantined(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            workdir = Path(tmpdir)
            (workdir / "point_000002").mkdir()
            (workdir / "point_000003").mkdir()
            (workdir / "notes").mkdir()

            moved = quarantine_uncheckpointed_point_dirs(workdir, 2)

            self.assertTrue((workdir / "point_000002").is_dir())
            self.assertFalse((workdir / "point_000003").exists())
            self.assertEqual(len(moved), 1)
            self.assertTrue(moved[0].name.startswith("point_000003.interrupted-"))
            self.assertTrue((workdir / "notes").is_dir())

    def test_live_and_stale_campaign_locks(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scan = Path(tmpdir) / "scan.dat"
            first = CampaignLock(scan, ["test"])
            first.acquire()
            try:
                with self.assertRaises(CampaignStateError):
                    CampaignLock(scan, ["second"]).acquire()
            finally:
                first.release()

            lock = scan.with_suffix(".lock")
            atomic_write_json(
                lock,
                {
                    "schema": "trsm_scan_lock_v1",
                    "hostname": socket.gethostname(),
                    "pid": 99_999_999,
                    "token": "stale",
                },
            )
            replacement = CampaignLock(scan, ["replacement"])
            replacement.acquire()
            try:
                self.assertIsNotNone(replacement.archived_stale_lock)
                self.assertTrue(replacement.archived_stale_lock.exists())
            finally:
                replacement.release()


class TestGeneratorResume(unittest.TestCase):
    def test_resume_cli_hydrates_saved_settings_and_total_target(self):
        generator = load_generator_module()
        configure_numeric_stubs(generator)
        with tempfile.TemporaryDirectory() as tmpdir:
            scan = Path(tmpdir) / "trsm_points_test.dat"
            args = scan_args(generator, 879, 100_000)
            metadata = generator.build_scan_metadata(
                args,
                "test",
                scan,
            )
            generator.write_scan_metadata_file(
                scan.with_suffix(".metadata.json"),
                metadata,
            )

            resumed = generator.parse_args(
                [
                    "--resume-from",
                    str(scan),
                    "--nrandom",
                    "10000",
                    "--no-print-info",
                ]
            )

            self.assertEqual(resumed.seed, 879)
            self.assertEqual(resumed.nrandom, 10_000)
            self.assertTrue(resumed.nrandom_count_evo_thc)
            self.assertTrue(resumed.write_evo_thc_points)
            self.assertTrue(resumed.scan_k133_k233_log)
            self.assertTrue(resumed.independent_m3)
            self.assertFalse(resumed.print_info)

    def test_resume_cli_rejects_immutable_override(self):
        generator = load_generator_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            scan = Path(tmpdir) / "trsm_points_test.dat"
            args = scan_args(generator, 879, 10)
            metadata = generator.build_scan_metadata(args, "test", scan)
            generator.write_scan_metadata_file(
                scan.with_suffix(".metadata.json"),
                metadata,
            )
            with self.assertRaises(SystemExit):
                generator.parse_args(
                    [
                        "--resume-from",
                        str(scan),
                        "--nrandom",
                        "20",
                        "--independent-m3",
                    ]
                )

    def test_checkpointless_replay_recovers_exact_next_candidate(self):
        generator = load_generator_module()
        configure_numeric_stubs(generator)
        args = scan_args(generator, 879, 100)
        reference_rng = random.Random(879)
        selected = []
        candidate = None
        for draw in range(1, 9):
            candidate = generator.draw_random_vxzero_candidate(
                args,
                rng=reference_rng,
            )
            if draw in (2, 5, 8):
                selected.append(candidate)
        expected_next = generator.draw_random_vxzero_candidate(
            args,
            rng=reference_rng,
        )

        columns = [
            "M2",
            "M3",
            "vs",
            "a12",
            "lX",
            "lPhiX",
            "lSX",
            "K133",
            "K233",
            "evo",
            "thc",
            "hb",
            "hs",
            "ewpo",
            "wmass",
            "dm",
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            scan = Path(tmpdir) / "scan.dat"
            rows = []
            for index, item in enumerate(selected):
                row = {
                    "M2": item["m2"],
                    "M3": item["m3"],
                    "vs": item["vs"],
                    "a12": item["a12"],
                    "lX": item["lX"],
                    "lPhiX": item["lPhiX"],
                    "lSX": item["lSX"],
                    "K133": item["K133"],
                    "K233": item["K233"],
                    "evo": True,
                    "thc": True,
                    "hb": index == 0,
                    "hs": index == 0,
                    "ewpo": index == 0,
                    "wmass": index == 0,
                    "dm": index == 0,
                }
                rows.append(row)
            with scan.open("w", encoding="ascii") as stream:
                stream.write("\t".join(columns) + "\n")
                for row in rows:
                    stream.write(
                        "\t".join(str(row[column]) for column in columns) + "\n"
                    )
            args.resume_metadata = generator.build_scan_metadata(
                args,
                "test",
                scan,
            )
            recovered_rng = random.Random()

            recovered = generator.recover_legacy_campaign(
                scan,
                args,
                recovered_rng,
                columns,
            )
            actual_next = generator.draw_random_vxzero_candidate(
                args,
                rng=recovered_rng,
            )

            self.assertEqual(recovered["draw_count"], 8)
            self.assertEqual(recovered["evo_thc_count"], 3)
            self.assertEqual(recovered["viable_count"], 1)
            for key in generator.LEGACY_CANDIDATE_COLUMNS.values():
                self.assertEqual(actual_next[key], expected_next[key])

    def _run_fake_campaign(self, generator, directory, target):
        args = scan_args(generator, 73, target)
        args.checkpoint_every = 1
        args.print_info = False
        generator.cli_args = args
        generator.ini_seed = args.seed
        generator.nrandom = args.nrandom
        generator.OutputDir = str(directory) + os.sep
        generator.RunTag = "resume-test"
        runtime = generator.prepare_campaign_runtime(generator.RunTag, args)
        path = Path(generator.output_path(generator.RunTag))

        def fake_evaluate(
            seed,
            m2,
            m3,
            vs,
            a12,
            lx,
            lphix,
            lsx,
            **kwargs,
        ):
            if path.stat().st_size == 0:
                with path.open("a", encoding="ascii") as stream:
                    stream.write("M2\tM3\tvs\ta12\tlX\tlPhiX\tlSX\n")
            with path.open("a", encoding="ascii") as stream:
                stream.write(
                    "\t".join(
                        str(value)
                        for value in (m2, m3, vs, a12, lx, lphix, lsx)
                    )
                    + "\n"
                )
            return {"viable": 1, "evo_thc": True}

        generator.evaluate_trsm_point_vxzero = fake_evaluate
        generator.run_random_vxzero_scan(runtime)
        return path

    def test_interrupted_prefix_then_resume_matches_uninterrupted_output(self):
        generator = load_generator_module()
        configure_numeric_stubs(generator)
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            uninterrupted = root / "uninterrupted"
            split = root / "split"
            uninterrupted.mkdir()
            split.mkdir()

            full_path = self._run_fake_campaign(generator, uninterrupted, 7)
            full_bytes = full_path.read_bytes()
            split_path = self._run_fake_campaign(generator, split, 3)

            metadata = json.loads(
                split_path.with_suffix(".metadata.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(metadata["campaign_progress"]["evo_thc_count"], 3)
            resume_args = generator.parse_args(
                [
                    "--resume-from",
                    str(split_path),
                    "--nrandom",
                    "7",
                    "--no-print-info",
                ]
            )
            generator.cli_args = resume_args
            generator.ini_seed = resume_args.seed
            generator.nrandom = resume_args.nrandom
            generator.OutputDir = str(split) + os.sep
            generator.RunTag = resume_args.resume_metadata["run_tag"]
            runtime = generator.prepare_campaign_runtime(
                generator.RunTag,
                resume_args,
            )

            def append_evaluated(seed, m2, m3, vs, a12, lx, lphix, lsx, **kwargs):
                with split_path.open("a", encoding="ascii") as stream:
                    stream.write(
                        "\t".join(
                            str(value)
                            for value in (m2, m3, vs, a12, lx, lphix, lsx)
                        )
                        + "\n"
                    )
                return {"viable": 1, "evo_thc": True}

            generator.evaluate_trsm_point_vxzero = append_evaluated
            generator.run_random_vxzero_scan(runtime)

            self.assertEqual(split_path.read_bytes(), full_bytes)
            final_metadata = json.loads(
                split_path.with_suffix(".metadata.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(final_metadata["initial_requested_points"], 3)
            self.assertEqual(final_metadata["requested_points"], 7)
            self.assertEqual(final_metadata["campaign_status"], "complete")
            self.assertEqual(len(final_metadata["resume_history"]), 1)

    def test_resume_rejects_target_below_completed_total(self):
        generator = load_generator_module()
        configure_numeric_stubs(generator)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self._run_fake_campaign(generator, Path(tmpdir), 3)
            args = generator.parse_args(
                [
                    "--resume-from",
                    str(path),
                    "--nrandom",
                    "2",
                ]
            )
            generator.cli_args = args
            generator.OutputDir = str(Path(tmpdir)) + os.sep
            generator.RunTag = args.resume_metadata["run_tag"]
            with self.assertRaises(CampaignStateError):
                generator.prepare_campaign_runtime(generator.RunTag, args)

    def test_equal_target_resume_is_a_complete_noop(self):
        generator = load_generator_module()
        configure_numeric_stubs(generator)
        with tempfile.TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            path = self._run_fake_campaign(generator, directory, 3)
            original = path.read_bytes()
            args = generator.parse_args(
                ["--resume-from", str(path), "--nrandom", "3"]
            )
            generator.cli_args = args
            generator.OutputDir = str(directory) + os.sep
            generator.RunTag = args.resume_metadata["run_tag"]
            runtime = generator.prepare_campaign_runtime(generator.RunTag, args)
            generator.evaluate_trsm_point_vxzero = lambda *a, **k: self.fail(
                "equal-target resume evaluated a new point"
            )

            result = generator.run_random_vxzero_scan(runtime)

            self.assertEqual(result, 3)
            self.assertEqual(path.read_bytes(), original)
            metadata = json.loads(
                path.with_suffix(".metadata.json").read_text(encoding="utf-8")
            )
            self.assertEqual(metadata["campaign_status"], "complete")

    def test_corrupt_checkpoint_is_rejected(self):
        generator = load_generator_module()
        configure_numeric_stubs(generator)
        with tempfile.TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            path = self._run_fake_campaign(generator, directory, 2)
            path.with_suffix(".checkpoint.json").write_text(
                "{not-json",
                encoding="utf-8",
            )
            args = generator.parse_args(
                ["--resume-from", str(path), "--nrandom", "3"]
            )
            generator.cli_args = args
            generator.OutputDir = str(directory) + os.sep
            generator.RunTag = args.resume_metadata["run_tag"]
            with self.assertRaises(CampaignStateError):
                generator.prepare_campaign_runtime(generator.RunTag, args)

    def test_configuration_fingerprint_mismatch_is_rejected(self):
        generator = load_generator_module()
        configure_numeric_stubs(generator)
        with tempfile.TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            path = self._run_fake_campaign(generator, directory, 2)
            args = generator.parse_args(
                ["--resume-from", str(path), "--nrandom", "3"]
            )
            generator.cli_args = args
            generator.OutputDir = str(directory) + os.sep
            generator.RunTag = args.resume_metadata["run_tag"]
            original_m3_max = generator.m3_max
            generator.m3_max = original_m3_max + 1
            try:
                with self.assertRaises(CampaignStateError):
                    generator.prepare_campaign_runtime(generator.RunTag, args)
            finally:
                generator.m3_max = original_m3_max

    def test_fresh_run_refuses_to_overwrite_campaign_state(self):
        generator = load_generator_module()
        configure_numeric_stubs(generator)
        with tempfile.TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            self._run_fake_campaign(generator, directory, 2)
            args = scan_args(generator, 73, 2)
            generator.cli_args = args
            generator.OutputDir = str(directory) + os.sep
            generator.RunTag = "resume-test"
            with self.assertRaises(CampaignStateError):
                generator.prepare_campaign_runtime(generator.RunTag, args)


if __name__ == "__main__":
    unittest.main()
