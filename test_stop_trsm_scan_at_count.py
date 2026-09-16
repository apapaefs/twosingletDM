import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import stop_trsm_scan_at_count as monitor


HEADER = "M2\tevo\tthc\n"


def row(m2="300", evo="True", thc="True"):
    return f"{m2}\t{evo}\t{thc}\n"


class EvoThcFileCounterTests(unittest.TestCase):
    def test_counts_existing_and_incremental_complete_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "scan.dat"
            path.write_text(HEADER + row("300") + row("301"), encoding="utf-8")
            counter = monitor.EvoThcFileCounter(path)
            self.assertEqual(counter.count, 2)

            with path.open("ab") as stream:
                stream.write(b"302\tTrue\t")
                stream.flush()
                os.fsync(stream.fileno())
            self.assertEqual(counter.update(), 2)
            self.assertGreater(counter.pending_bytes, 0)

            with path.open("ab") as stream:
                stream.write(b"True\n")
                stream.flush()
                os.fsync(stream.fileno())
            self.assertEqual(counter.update(), 3)
            self.assertEqual(counter.pending_bytes, 0)

    def test_rejects_a_row_that_does_not_pass_both_flags(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "scan.dat"
            path.write_text(HEADER + row(evo="False"), encoding="utf-8")
            with self.assertRaisesRegex(
                monitor.StopMonitorError, "not an evo/thc-passing point"
            ):
                monitor.EvoThcFileCounter(path)

    def test_rejects_a_malformed_row(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "scan.dat"
            path.write_text(HEADER + "300\tTrue\n", encoding="utf-8")
            with self.assertRaisesRegex(monitor.StopMonitorError, "Malformed TSV"):
                monitor.EvoThcFileCounter(path)

    def test_rejects_output_truncation(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "scan.dat"
            path.write_text(HEADER + row(), encoding="utf-8")
            counter = monitor.EvoThcFileCounter(path)
            path.write_text(HEADER, encoding="utf-8")
            with self.assertRaisesRegex(monitor.StopMonitorError, "shrank"):
                counter.update()


class MetadataAndProcessTests(unittest.TestCase):
    def test_metadata_requires_counted_written_evo_thc_points(self):
        with tempfile.TemporaryDirectory() as directory:
            scan = Path(directory) / "scan.dat"
            scan.write_text(HEADER, encoding="utf-8")
            metadata = {
                "seed": 879,
                "scan_file": scan.name,
                "options": {
                    "nrandom_count_evo_thc": True,
                    "write_evo_thc_points": True,
                },
            }
            monitor.metadata_path(scan).write_text(
                json.dumps(metadata), encoding="utf-8"
            )
            _, loaded, seed = monitor.load_campaign_metadata(scan)
            self.assertEqual(seed, 879)
            self.assertEqual(loaded, metadata)

            metadata["options"]["write_evo_thc_points"] = False
            monitor.metadata_path(scan).write_text(
                json.dumps(metadata), encoding="utf-8"
            )
            with self.assertRaisesRegex(
                monitor.StopMonitorError, "--write-evo-thc-points"
            ):
                monitor.load_campaign_metadata(scan)

    def test_parse_and_validate_generator_process(self):
        record = (
            " 3616  3616 Thu Jul 23 13:02:18 2026     "
            "/usr/bin/python3 generate_trsm_points.py 879 --nrandom 100000 "
            "--nrandom-count-evo-thc --write-evo-thc-points"
        )
        process = monitor.parse_ps_line(record)
        self.assertEqual(process.pid, 3616)
        self.assertEqual(process.pgid, 3616)
        self.assertEqual(process.start_time, "Thu Jul 23 13:02:18 2026")
        self.assertIsNone(monitor.scan_process_mismatch(process, 879))
        self.assertIn("seed", monitor.scan_process_mismatch(process, 888))

    def test_select_pid_rejects_an_unrelated_process(self):
        process = monitor.ProcessIdentity(
            12,
            12,
            "Thu Jul 23 13:02:18 2026",
            "/usr/bin/python3 something_else.py",
        )
        with mock.patch.object(monitor, "get_process_identity", return_value=process):
            with self.assertRaisesRegex(monitor.StopMonitorError, "unsafe to signal"):
                monitor.select_scan_process(879, 12)


class MonitorLockTests(unittest.TestCase):
    def test_live_monitor_lock_is_not_replaced(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "monitor.lock"
            path.write_text(
                json.dumps(
                    {
                        "schema": monitor.LOCK_SCHEMA,
                        "hostname": monitor.socket.gethostname(),
                        "pid": os.getpid(),
                        "token": "other",
                    }
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                monitor.StopMonitorError, "live stop monitor"
            ):
                monitor.MonitorLock(path).acquire()


if __name__ == "__main__":
    unittest.main()
