"""CalcGW parsing, point identity, and transactional campaign regression checks."""

import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

from reprocess_trsm_gw import reprocess
from trsm_gw import GW_COLUMNS, parse_gw_row, run_calcgw


def native_row():
    row = dict(index=7, status_nlo_stability="success", status_tracing="success", status_coex_pairs="success")
    for index, alpha, beta in ((0, .04, 200), (3, .7, 30)):
        row.update({f"status_gw_{index}": "success", f"status_perc_{index}": "success",
                    f"T_star_{index}": 100 - index, f"alpha_PT_{index}": alpha, f"beta/H_{index}": beta,
                    f"f_2_sw_{index}": .001, f"h2Omega_2_sw_{index}": 1e-12,
                    f"SNR(LISA-3yrs)_{index}": 10, f"kappa_col_{index}": 0})
        row.update({f"v_wall_{index}": .95, f"eps_turb_{index}": .1})
    return row


def scan_row():
    return dict(M2="326.5", M3="23.93", vs="224.76", vx="0", a12="0.1116",
                lX="7.66", lPhiX="0.001689", lSX="-0.0001516", ewpt_eligible="True",
                mg5_xsec_hh_pb="0.01617", original_json='{"note":"preserve cells"}')


def write_scan(path):
    rows = [scan_row() for _ in range(4)]
    rows[1]["ewpt_eligible"] = "False"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


class CalcGWTests(unittest.TestCase):
    def test_pairing_native_statuses_and_no_zero_substitution(self):
        result = parse_gw_row(native_row())
        transitions = result["gw_transitions"]
        self.assertEqual([(t["transition_index"], t["alpha"], t["beta_over_H"]) for t in transitions],
                         [(0, .04, 200), (3, .7, 30)])
        self.assertEqual(result["gw_successful_transition_count"], 2)
        self.assertEqual(transitions[0]["kappa_col"], 0)
        self.assertIsNone(transitions[0]["h2omega_col"])
        for key, value in (("status_gw_0", "failure"), ("status_perc_0", "failure"),
                           ("alpha_PT_0", "nan"), ("alpha_PT_0", 0), ("beta/H_0", .9),
                           ("T_star_0", -1)):
            row = native_row()
            row[key] = value
            result = parse_gw_row(row)
            self.assertEqual(result["gw_successful_transition_count"], 1)
            self.assertFalse(result["gw_transitions"][0]["valid"])
        row = native_row()
        row["status_nlo_stability"] = "failure"
        self.assertEqual(parse_gw_row(row)["gw_successful_transition_count"], 0)
        self.assertEqual(parse_gw_row(row)["gw_status"], "native_failure")

    def test_missing_native_output_is_not_a_success(self):
        with self.assertRaisesRegex(ValueError, "statuses"):
            parse_gw_row({})
        row = dict(status_nlo_stability="success", status_tracing="success", status_coex_pairs="no_coex_pair",
                   transition_history=float("nan"))
        self.assertEqual(parse_gw_row(row)["gw_status"], "no_successful_gw")
        self.assertIsNone(json.loads(json.dumps(parse_gw_row(row), allow_nan=False))["gw_native_summary"]["transition_history"])
        row["status_coex_pairs"] = "success"
        with self.assertRaisesRegex(ValueError, "no GW fields"):
            parse_gw_row(row)

    def test_native_adapter_parameters_identity_and_timeout(self):
        options = dict(executable="/fake/CalcGW", thigh=1000., trans_temp="perc", vwall=.95, epsturb=.1, timeout=2.)
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)

            def native(command, **kwargs):
                self.assertIn("--thigh=1000.0", command)
                self.assertIn("--usemultithreading=false", command)
                self.assertIn("--trans_temp=perc", command)
                self.assertEqual(kwargs["env"]["OMP_NUM_THREADS"], "1")
                row = native_row()
                with (directory / "CalcGW.tsv").open("w") as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(row), delimiter="\t")
                    writer.writeheader()
                    writer.writerow(row)
                return mock.Mock(returncode=0)

            with mock.patch("trsm_gw.subprocess.run", side_effect=native):
                result = run_calcgw(scan_row(), 7, directory, options)
            self.assertEqual(result["gw_status"], "success")
            fields = (directory / "TRSM_Input.tsv").read_text().splitlines()[1].split("\t")
            self.assertEqual(float(fields[0]), 7)
            self.assertEqual(float(fields[-1]), float(scan_row()["lSX"]))
            # A stale output cannot turn a new empty run into a successful calculation.
            with mock.patch("trsm_gw.subprocess.run", return_value=mock.Mock(returncode=0)):
                self.assertEqual(run_calcgw(scan_row(), 7, directory, options)["gw_status"], "invalid_output")
            import subprocess
            with mock.patch("trsm_gw.subprocess.run", side_effect=subprocess.TimeoutExpired("CalcGW", 2)):
                self.assertEqual(run_calcgw(scan_row(), 7, directory, options)["gw_status"], "timeout")

    def test_pilot_resume_parallel_matches_serial_and_preserves_cells(self):
        with tempfile.TemporaryDirectory() as temp:
            source, output, serial = (Path(temp) / name for name in ("input.tsv", "gw.tsv", "serial.tsv"))
            original = write_scan(source)
            before = source.read_bytes()
            source.with_suffix(".metadata.json").write_text('{"campaign":"100k"}')
            calls = []

            def runner(row, index, directory, options):
                calls.append(index)
                return parse_gw_row(native_row())

            kwargs = dict(executable=Path(sys.executable), runner=runner)
            counts = reprocess(source, output, jobs=2, max_points=2, **kwargs)
            self.assertFalse(counts["complete"])
            self.assertFalse(output.exists())
            self.assertEqual(set(calls), {1, 3})
            counts = reprocess(source, output, jobs=1, resume=True, **kwargs)
            self.assertTrue(counts["complete"])
            self.assertEqual(calls.count(1), 1)
            self.assertEqual(calls.count(3), 1)
            self.assertEqual(calls.count(4), 1)
            self.assertEqual(counts["successful_transitions"], 6)
            with output.open() as stream:
                records = list(csv.DictReader(stream, delimiter="\t"))
            for actual, expected in zip(records, original):
                self.assertEqual({key: actual[key] for key in expected}, expected)
            self.assertEqual(records[1]["gw_status"], "not_selected")
            self.assertEqual(source.read_bytes(), before)
            self.assertEqual(json.loads(output.with_suffix(".metadata.json").read_text())["campaign"], "100k")
            reprocess(source, serial, jobs=1, **kwargs)
            self.assertEqual(output.read_bytes(), serial.read_bytes())
            calls.clear()
            reprocess(source, output, resume=True, **kwargs)
            self.assertFalse(calls)
            with self.assertRaises(FileExistsError):
                reprocess(source, output, **kwargs)

    def test_resume_rejects_changed_configuration_and_source(self):
        with tempfile.TemporaryDirectory() as temp:
            source, output = Path(temp) / "scan.tsv", Path(temp) / "gw.tsv"
            write_scan(source)
            kwargs = dict(executable=Path(sys.executable), runner=lambda *args: parse_gw_row(native_row()))
            reprocess(source, output, max_points=1, **kwargs)
            with self.assertRaisesRegex(ValueError, "checkpoint does not match"):
                reprocess(source, output, resume=True, vwall=.8, **kwargs)
            source.write_text(source.read_text().replace("326.5", "326.6"))
            with self.assertRaisesRegex(ValueError, "checkpoint does not match"):
                reprocess(source, output, resume=True, **kwargs)
            with self.assertRaisesRegex(ValueError, "overlap"):
                reprocess(source, source, **kwargs)


if __name__ == "__main__":
    unittest.main()
