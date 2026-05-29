import csv
import contextlib
import importlib.util
import io
import json
import sys
import tempfile
import textwrap
import unittest
from concurrent.futures import Future
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().parent / "run_trsm_seed_campaign.py"


def load_module():
    spec = importlib.util.spec_from_file_location("run_trsm_seed_campaign", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_fake_generator(path):
    path.write_text(
        textwrap.dedent(
            r'''
            import argparse
            import csv
            import sys
            from datetime import date
            from pathlib import Path

            parser = argparse.ArgumentParser()
            parser.add_argument("seed", type=int)
            parser.add_argument("--nrandom", type=int)
            parser.add_argument("--sleep", type=float, default=0.0)
            parser.add_argument("--run-ewpt", action="store_true")
            parser.add_argument("--run-ewpt-on-dm-failed", action="store_true")
            parser.add_argument("--ewpt-workdir", type=Path)
            parser.add_argument("--ewpt-require-eq418", action="store_true")
            parser.add_argument("--ewpt-thigh")
            parser.add_argument("--ewpt-plot-phases", action="store_true")
            parser.add_argument("--ewpt-plot-format")
            parser.add_argument("--ewpt-sym-threshold")
            parser.add_argument("--ewpt-w1-threshold")
            parser.add_argument("--ewpt-wx-threshold")
            parser.add_argument("--ewpt-ws-threshold")
            args = parser.parse_args()

            print(f"fake generator seed={args.seed} nrandom={args.nrandom}")
            sys.stdout.flush()
            if args.sleep:
                import time
                time.sleep(args.sleep)
            if args.seed == 2:
                print("intentional fake failure", file=sys.stderr)
                sys.exit(7)

            tag = f"13.6-{date.today().strftime('%Y%m%d')}-{args.seed}-False_vxzero"
            output_path = Path.cwd() / "output" / f"trsm_points_{tag}.dat"
            output_path.parent.mkdir(parents=True, exist_ok=True)
            rows_by_seed = {
                1: [
                    ["160", "68", "50", "0", "-0.18", "0", "0", "1.05", "0.34", "0.36", "True"],
                    ["380", "500", "200", "0", "-0.15", "0", "0", "0.10", "0.044", "0.15", "True"],
                ],
                3: [
                    ["340", "500", "120", "0", "-0.20", "0", "0", "0.15", "0.035", "0.10", "True"],
                ],
            }
            header = ["M2", "M3", "vs", "vx", "a12", "a13", "a23", "lX", "lPhiX", "lSX", "dm"]
            with output_path.open("w", encoding="ascii", newline="") as stream:
                writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
                writer.writerow(header)
                writer.writerows(rows_by_seed.get(args.seed, []))

            if args.run_ewpt and args.ewpt_workdir is not None:
                args.ewpt_workdir.mkdir(parents=True, exist_ok=True)
                strengths = {
                    1: [
                        ("point_000001", "fopt_strength_crit_0: T=150, ew_jump/T=0.9, ew_true/T=0.9, field_jump/T=1.0\nfopt_strength_nucl_0: T=150, ew_jump/T=0.4, ew_true/T=0.4, field_jump/T=0.5\n", ["160", "68", "50", "-0.18", "1.05", "0.34", "0.36"]),
                        ("point_000002", "fopt_strength_perc_0: T=150, ew_jump/T=0.6, ew_true/T=0.6, field_jump/T=0.7\n", ["380", "500", "200", "-0.15", "0.10", "0.044", "0.15"]),
                    ],
                    3: [
                        ("point_000001", "fopt_strength_nucl_0: T=140, ew_jump/T=0.8, ew_true/T=0.8, field_jump/T=0.9\n", ["340", "500", "120", "-0.20", "0.15", "0.035", "0.10"]),
                    ],
                }
                for directory, summary, values in strengths.get(args.seed, []):
                    point_dir = args.ewpt_workdir / directory
                    point_dir.mkdir(parents=True, exist_ok=True)
                    (point_dir / "ewpt_summary.txt").write_text(summary, encoding="ascii")
                    with (point_dir / "TRSM_Input.tsv").open("w", encoding="ascii", newline="") as stream:
                        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
                        writer.writerow(["", "m1", "m2", "m3", "vs", "a12", "lx", "lphix", "lsx"])
                        writer.writerow([directory.split("_")[-1], "125.09"] + values)
            '''
        ).lstrip(),
        encoding="ascii",
    )


class RecordingExecutor:
    max_workers_seen = []

    def __init__(self, max_workers=None):
        self.max_workers_seen.append(max_workers)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def submit(self, fn, *args, **kwargs):
        future = Future()
        try:
            future.set_result(fn(*args, **kwargs))
        except BaseException as error:
            future.set_exception(error)
        return future


class TestTRSMSeedCampaign(unittest.TestCase):
    def test_seed_range_expansion(self):
        campaign = load_module()

        self.assertEqual(campaign.expand_seeds(10, 4), [10, 11, 12, 13])

    def test_builds_generator_command_with_forwarded_options(self):
        campaign = load_module()
        args = campaign.parse_args(
            [
                "--seed-start",
                "5",
                "--nseeds",
                "1",
                "--nrandom",
                "250",
                "--campaign-dir",
                "campaign",
                "--run-ewpt",
                "--write-dm-failed",
                "--run-ewpt-on-dm-failed",
                "--ewpt-require-eq418",
                "--ewpt-thigh",
                "1000",
                "--ewpt-plot-phases",
                "--ewpt-wx-threshold",
                "2.0",
            ]
        )

        command = campaign.build_generator_command(
            args,
            seed=5,
            generator_script=Path("/tmp/fake_generator.py"),
            ewpt_workdir=Path("/tmp/ewpt/seed_5"),
        )

        self.assertEqual(command[:4], [sys.executable, "-u", "/tmp/fake_generator.py", "5"])
        self.assertIn("--nrandom", command)
        self.assertIn("250", command)
        self.assertIn("--run-ewpt", command)
        self.assertIn("--write-dm-failed", command)
        self.assertIn("--run-ewpt-on-dm-failed", command)
        self.assertIn("--ewpt-require-eq418", command)
        self.assertIn("--ewpt-workdir", command)
        self.assertIn("/tmp/ewpt/seed_5", command)
        self.assertIn("--ewpt-wx-threshold", command)
        self.assertIn("2.0", command)

    def test_python_executable_can_be_overridden(self):
        campaign = load_module()
        args = campaign.parse_args(
            [
                "--seed-start",
                "5",
                "--nseeds",
                "1",
                "--campaign-dir",
                "campaign",
                "--python-executable",
                "/tmp/custom-python",
            ]
        )

        command = campaign.build_generator_command(
            args,
            seed=5,
            generator_script=Path("/tmp/fake_generator.py"),
            ewpt_workdir=Path("/tmp/ewpt/seed_5"),
        )

        self.assertEqual(command[0], "/tmp/custom-python")

    def test_dm_failed_ewpt_only_still_forwards_campaign_workdir(self):
        campaign = load_module()
        args = campaign.parse_args(
            [
                "--seed-start",
                "5",
                "--nseeds",
                "1",
                "--campaign-dir",
                "campaign",
                "--run-ewpt-on-dm-failed",
            ]
        )

        command = campaign.build_generator_command(
            args,
            seed=5,
            generator_script=Path("/tmp/fake_generator.py"),
            ewpt_workdir=Path("/tmp/ewpt/seed_5"),
        )

        self.assertIn("--run-ewpt-on-dm-failed", command)
        self.assertIn("--ewpt-workdir", command)
        self.assertIn("/tmp/ewpt/seed_5", command)
        self.assertNotIn("--run-ewpt", command)

    def test_default_python_executable_prefers_active_virtualenv(self):
        campaign = load_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            python_path = Path(tmpdir) / "bin" / "python"
            python_path.parent.mkdir()
            python_path.write_text("", encoding="ascii")

            self.assertEqual(
                campaign.default_python_executable({"VIRTUAL_ENV": tmpdir}),
                python_path,
            )

    def test_ewpt_strength_ranking_prefers_nucleation_over_fallbacks(self):
        campaign = load_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            summary = Path(tmpdir) / "ewpt_summary.txt"
            summary.write_text(
                "\n".join(
                    [
                        "fopt_strength_crit_0: T=150, ew_jump/T=0.9, ew_true/T=0.9, field_jump/T=0.9",
                        "fopt_strength_nucl_0: T=150, ew_jump/T=0.4, ew_true/T=0.4, field_jump/T=0.4",
                    ]
                ),
                encoding="ascii",
            )

            best = campaign.best_strength_from_summary(summary)

        self.assertEqual(best.kind, "nucl")
        self.assertEqual(best.value, 0.4)

    def test_parallel_job_limit_wiring(self):
        campaign = load_module()
        args = campaign.parse_args(
            ["--seed-start", "1", "--nseeds", "2", "--campaign-dir", "campaign", "--jobs", "3"]
        )

        list(
            campaign.run_seed_jobs(
                args,
                [1, 2],
                run_seed_func=lambda seed, args: campaign.SeedResult(seed=seed, returncode=0),
                executor_cls=RecordingExecutor,
            )
        )

        self.assertEqual(RecordingExecutor.max_workers_seen[-1], 3)

    def test_campaign_runs_fake_generator_and_writes_aggregate_outputs(self):
        campaign = load_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            fake_generator = tmpdir / "fake_generate_trsm_points.py"
            write_fake_generator(fake_generator)
            campaign_dir = tmpdir / "campaign"
            run_cwd = tmpdir / "run"
            run_cwd.mkdir()

            args = campaign.parse_args(
                [
                    "--seed-start",
                    "1",
                    "--nseeds",
                    "3",
                    "--nrandom",
                    "20",
                    "--jobs",
                    "2",
                    "--campaign-dir",
                    str(campaign_dir),
                    "--generator-script",
                    str(fake_generator),
                    "--run-cwd",
                    str(run_cwd),
                    "--run-ewpt",
                    "--ewpt-require-eq418",
                    "--ewpt-thigh",
                    "1000",
                ]
            )

            result = campaign.run_campaign(args)

            self.assertEqual(len(result.seed_results), 3)
            self.assertEqual({seed.seed for seed in result.seed_results}, {1, 2, 3})
            failed = [seed for seed in result.seed_results if seed.seed == 2][0]
            self.assertEqual(failed.returncode, 7)

            summary_path = campaign_dir / "campaign_summary.tsv"
            combined_path = campaign_dir / "combined_points.tsv"
            best_path = campaign_dir / "best_points.tsv"
            json_path = campaign_dir / "campaign_summary.json"
            for path in [summary_path, combined_path, best_path, json_path]:
                self.assertTrue(path.exists(), path)

            with summary_path.open(encoding="ascii", newline="") as stream:
                summary_rows = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual(len(summary_rows), 3)
            self.assertEqual(summary_rows[0]["viable_count"], "2")
            self.assertEqual(summary_rows[1]["returncode"], "7")
            self.assertEqual(summary_rows[2]["best_ew_jump_over_T"], "0.8")

            with combined_path.open(encoding="ascii", newline="") as stream:
                combined_rows = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual(len(combined_rows), 3)
            self.assertEqual(combined_rows[0]["seed"], "1")
            self.assertEqual(combined_rows[-1]["M2"], "340")

            with best_path.open(encoding="ascii", newline="") as stream:
                best_rows = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual(best_rows[0]["seed"], "3")
            self.assertEqual(best_rows[0]["best_ew_jump_over_T"], "0.8")
            self.assertEqual(best_rows[1]["strength_kind"], "nucl")
            self.assertEqual(best_rows[1]["best_ew_jump_over_T"], "0.4")

            payload = json.loads(json_path.read_text(encoding="ascii"))
            self.assertEqual(payload["metadata"]["seed_start"], 1)
            self.assertEqual(payload["metadata"]["nseeds"], 3)
            self.assertEqual(payload["totals"]["viable_count"], 3)
            self.assertEqual(payload["totals"]["failed_seeds"], 1)
            self.assertTrue((campaign_dir / "logs" / "seed_1.log").exists())

    def test_run_seed_streams_subprocess_output_to_log_before_completion(self):
        campaign = load_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            fake_generator = tmpdir / "fake_generate_trsm_points.py"
            write_fake_generator(fake_generator)
            campaign_dir = tmpdir / "campaign"
            run_cwd = tmpdir / "run"
            run_cwd.mkdir()
            args = campaign.parse_args(
                [
                    "--seed-start",
                    "4",
                    "--nseeds",
                    "1",
                    "--nrandom",
                    "0",
                    "--campaign-dir",
                    str(campaign_dir),
                    "--generator-script",
                    str(fake_generator),
                    "--run-cwd",
                    str(run_cwd),
                    "--generator-extra-arg=--sleep",
                    "--generator-extra-arg",
                    "0.5",
                ]
            )

            result = campaign.run_seed(4, args)

            log_text = Path(result.log_path).read_text(encoding="utf-8")
            self.assertIn("command:", log_text)
            self.assertIn("fake generator seed=4 nrandom=0", log_text)
            self.assertIn("returncode: 0", log_text)

    def test_campaign_prints_immediate_launch_and_heartbeat_messages(self):
        campaign = load_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            fake_generator = tmpdir / "fake_generate_trsm_points.py"
            write_fake_generator(fake_generator)
            campaign_dir = tmpdir / "campaign"
            run_cwd = tmpdir / "run"
            run_cwd.mkdir()
            args = campaign.parse_args(
                [
                    "--seed-start",
                    "4",
                    "--nseeds",
                    "1",
                    "--nrandom",
                    "0",
                    "--jobs",
                    "1",
                    "--campaign-dir",
                    str(campaign_dir),
                    "--generator-script",
                    str(fake_generator),
                    "--run-cwd",
                    str(run_cwd),
                    "--heartbeat-seconds",
                    "0.1",
                    "--generator-extra-arg=--sleep",
                    "--generator-extra-arg",
                    "0.3",
                ]
            )
            stdout = io.StringIO()

            with contextlib.redirect_stdout(stdout):
                campaign.run_campaign(args)

            output = stdout.getvalue()
            self.assertIn("Starting TRSM seed campaign", output)
            self.assertIn("[launch] seed=4", output)
            self.assertIn("[running]", output)
            self.assertIn("[1/1 done]", output)


if __name__ == "__main__":
    unittest.main()
