import csv
import fcntl
import json
import math
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import reprocess_trsm_mg5 as mg
from scan_output import output_columns, output_row
from trsm_mg5_rates import derive_mg5_rates, stored_full_viability, ADDITIONAL_RATE_COLUMNS


class TestRates(unittest.TestCase):
    def test_components_and_incomplete_total(self):
        row = {"h2_h3h3_br": .25, "mg5_xsec_gg_heta0_pb": 4,
               "mg5_xsec_pp_eta0Z_pb": 2, "mg5_xsec_gg_eta0Z_pb": .5}
        result = derive_mg5_rates(row)
        self.assertEqual(result["mono_z_xsec_pb"], .5)
        self.assertEqual(result["mono_z_gg_xsec_pb"], .125)
        self.assertEqual(result["mono_z_total_xsec_pb"], .625)
        self.assertEqual(result["mg5_xsec_pp_eta0Z_total_pb"], 2.5)
        for invalid in (None, "nan", "inf", -1):
            row["mg5_xsec_gg_eta0Z_pb"] = invalid
            self.assertNotIn("mono_z_total_xsec_pb", derive_mg5_rates(row))
        row["mg5_xsec_gg_eta0Z_pb"] = 0
        self.assertEqual(derive_mg5_rates(row)["mono_z_gg_xsec_pb"], 0)
        row["h2_h3h3_br"] = 1.1
        self.assertNotIn("mono_z_total_xsec_pb", derive_mg5_rates(row))

    def test_selection_matches_versioned_plot_selection(self):
        row = {key: "True" for key in ("evo", "thc", "hb", "hs", "ewpo", "wmass", "flavour", "dm")}
        self.assertTrue(stored_full_viability(row))
        row["evo"] = "False"
        self.assertFalse(stored_full_viability(row))
        row.update(constraint_version="trsm_constraints_v2", experimental_subset="True")
        self.assertTrue(stored_full_viability(row))
        row["flavour"] = "nan"
        self.assertFalse(stored_full_viability(row))

    def test_existing_scan_schema_is_unchanged_without_new_process(self):
        old_processes = {"gg_heta0": math.nan, "pp_eta0Z": math.nan}
        old = output_columns(old_processes)
        self.assertTrue(all(column not in old for column in ADDITIONAL_RATE_COLUMNS))
        new_processes = {**old_processes, "gg_eta0Z": math.nan}
        new = output_columns(new_processes)
        self.assertEqual(new[-3:], list(ADDITIONAL_RATE_COLUMNS))
        self.assertEqual(len(new), len(output_row({}, new_processes).split("\t")))


class TestEnrichment(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.source, self.output = self.root / "source.tsv", self.root / "result.tsv"
        self.rows = [self.row(), self.row(dm="False"), self.row(point_index="1")]
        self.write_source()
        self.calls = []

    def row(self, **updates):
        row = {key: "1.0" for key in mg.LAMBDA_COLUMNS}
        row.update({key: "True" for key in ("thc", "experimental_subset", "dm", "flavour")})
        row.update(constraint_version="trsm_constraints_v2", M2="250", M3="30",
                   w1=".004", w2=".2", w3="0", k1=".99", k2=".1", k3="0",
                   h2_h3h3_br="0.25", mg5_xsec_gg_heta0_pb="4.000000",
                   mg5_xsec_pp_eta0Z_pb="2.000000", mono_z_xsec_pb="0.500000",
                   point_index="1", untouched="an unchanged cell")
        row.update(updates)
        return row

    def write_source(self):
        with self.source.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, list(self.rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(self.rows)

    def rate(self, process, name, *args, **kwargs):
        self.calls.append((process, name, args, kwargs))
        return .5

    def run_enrichment(self, **kwargs):
        return mg.reprocess(self.source, self.output, get_xsec=kwargs.pop("get_xsec", self.rate),
                            runtime_receipt=kwargs.pop("runtime_receipt", lambda processes: {"processes": list(processes)}),
                            **kwargs)

    def read_output(self):
        with self.output.open() as stream:
            return list(csv.DictReader(stream, delimiter="\t"))

    def test_streamed_copy_preserves_existing_cells_and_adds_only_missing(self):
        original = self.source.read_bytes()
        self.source.with_suffix(".metadata.json").write_text(json.dumps({"variable_ranges": ["kept"]}))
        counts = self.run_enrichment()
        self.assertEqual(counts["eligible"], 2)
        self.assertEqual(len(self.calls), 2)
        self.assertEqual({call[0] for call in self.calls}, {"gg_eta0Z"})
        self.assertNotEqual(self.calls[0][1], self.calls[1][1])
        self.assertEqual(self.calls[0][3]["ecm"], 13.6)
        self.assertEqual(self.calls[0][3]["w1"], .004)
        result = self.read_output()
        for original_row, new_row in zip(self.rows, result):
            self.assertEqual({key: new_row[key] for key in original_row}, original_row)
        self.assertEqual(result[0]["mono_z_total_xsec_pb"], "0.625")
        self.assertEqual(result[1]["mg5_xsec_gg_eta0Z_pb"], "nan")
        self.assertEqual(result[1]["mono_z_total_xsec_pb"], "nan")
        self.assertEqual(self.source.read_bytes(), original)
        metadata = json.loads(self.output.with_suffix(".metadata.json").read_text())
        self.assertEqual(metadata["variable_ranges"], ["kept"])
        self.assertEqual(metadata["mg5_enrichment"]["configuration"]["source_energy"]["status"],
                         "assumed_from_requested_energy")
        self.run_enrichment(resume=True)
        self.assertEqual(len(self.calls), 2)
        with self.assertRaises(FileExistsError):
            self.run_enrichment()

    def test_resume_per_process_and_reject_mismatched_state(self):
        self.rows[0]["mg5_xsec_gg_heta0_pb"] = "nan"
        self.write_source()
        def interrupted(process, *args, **kwargs):
            if process == "gg_eta0Z":
                raise RuntimeError("native failure")
            return self.rate(process, *args, **kwargs)
        with self.assertRaisesRegex(RuntimeError, "data row 1, process gg_eta0Z"):
            self.run_enrichment(get_xsec=interrupted)
        self.assertFalse(self.output.exists())
        self.assertTrue(Path(str(self.output) + ".partial").exists())
        with self.assertRaisesRegex(ValueError, "does not match"):
            self.run_enrichment(resume=True, energy=14)
        with self.assertRaisesRegex(ValueError, "does not match"):
            self.run_enrichment(resume=True, runtime_receipt=lambda _: {"changed": True})
        self.source.write_text(self.source.read_text() + "\n")
        with self.assertRaisesRegex(ValueError, "does not match"):
            self.run_enrichment(resume=True)
        self.write_source()
        self.run_enrichment(resume=True)
        self.assertEqual([call[0] for call in self.calls].count("gg_heta0"), 1)
        self.assertEqual(len(self.calls), 3)
        self.assertFalse(Path(str(self.output) + ".partial").exists())

    def test_dry_run_and_invalid_native_result(self):
        counts = self.run_enrichment(dry_run=True)
        self.assertEqual(counts["requested_missing_rates"], 2)
        self.assertFalse(self.output.exists())
        self.assertFalse(Path(str(self.output) + ".partial").exists())
        self.assertFalse(self.calls)
        with self.assertRaisesRegex(RuntimeError, "invalid cross section"):
            self.run_enrichment(get_xsec=lambda *args, **kwargs: math.nan)

    def test_refuse_same_file_and_changed_completed_output(self):
        with self.assertRaisesRegex(ValueError, "paths must differ"):
            mg.reprocess(self.source, self.source)
        self.run_enrichment()
        self.output.write_text(self.output.read_text() + "\n")
        with self.assertRaises(FileExistsError):
            self.run_enrichment(resume=True)

    def test_preserve_source_sidecar_and_refuse_concurrent_writer(self):
        source_metadata = self.source.with_suffix(".metadata.json")
        source_metadata.write_text('{"kept": true}')
        for collision in (self.source.with_suffix(".dat"), source_metadata):
            with self.assertRaisesRegex(ValueError, "collides"):
                mg.reprocess(self.source, collision)
        self.assertEqual(source_metadata.read_text(), '{"kept": true}')
        alias = self.root / "alias.tsv"
        alias.with_suffix(".metadata.json").symlink_to(source_metadata)
        with self.assertRaisesRegex(ValueError, "collides"):
            mg.reprocess(self.source, alias)
        with Path(str(self.output) + ".lock").open("a+") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            with self.assertRaisesRegex(RuntimeError, "another MG5 enrichment"):
                self.run_enrichment()
        self.assertFalse(self.output.exists())
        self.assertFalse(self.calls)

    def test_source_energy_mismatch_rejected_before_native_call(self):
        variants = (
            {"fixed_parameters": [{"variable": "sqrt(s)", "value": 13, "unit": "TeV"}]},
            {"fixed": {"Energy": 13}},
            {"configuration": {"fixed": {"Energy": 13}}},
            {"mg5_enrichment": {"configuration": {"energy_TeV": 13}}},
        )
        for metadata in variants:
            with self.subTest(metadata=metadata):
                self.source.with_suffix(".metadata.json").write_text(json.dumps(metadata))
                with self.assertRaisesRegex(ValueError, "would mix existing and new rates"):
                    self.run_enrichment()
                self.assertFalse(self.calls)
                self.assertFalse(Path(str(self.output) + ".partial").exists())
        valid = {"fixed_parameters": [{"variable": "sqrt(s)", "value": 13600, "unit": "GeV"}]}
        self.source.with_suffix(".metadata.json").write_text(json.dumps(valid))
        self.run_enrichment()
        metadata = json.loads(self.output.with_suffix(".metadata.json").read_text())
        self.assertEqual(metadata["mg5_enrichment"]["configuration"]["source_energy"]["status"],
                         "verified_source_metadata")

    def test_cli_defaults_one_core_and_preserves_explicit_core_setting(self):
        for explicit, expected in ((None, "1"), ("4", "4")):
            with self.subTest(explicit=explicit), mock.patch.dict(os.environ, {}, clear=True):
                if explicit is not None:
                    os.environ["TRSM_MG5_CORES"] = explicit
                def receipt(processes):
                    return {"cores": os.environ.get("TRSM_MG5_CORES")}
                real_reprocess = mg.reprocess
                def run(*args, **kwargs):
                    return real_reprocess(*args, **kwargs, get_xsec=self.rate, runtime_receipt=receipt)
                self.output = self.root / f"result-{expected}.tsv"
                with mock.patch.object(mg, "reprocess", side_effect=run):
                    self.assertEqual(mg.main([str(self.source), "--output", str(self.output)]), 0)
                metadata = json.loads(self.output.with_suffix(".metadata.json").read_text())
                self.assertEqual(metadata["mg5_enrichment"]["configuration"]["runtime"]["cores"], expected)


if __name__ == "__main__":
    unittest.main()
