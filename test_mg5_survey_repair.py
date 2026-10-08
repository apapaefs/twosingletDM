import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock

from generate_mg5_trsm_xsecs import _madevent_lock
from tools.repair_mg5_survey import exact_helicity_card, guarded_source, pending_updates, repair_survey


SOURCE = """class Survey:
    def write_parameter(self, parralelization, Pdirs=None):
        options = {'maxiter': self.iterations, 'miniter': self.min_iterations}
        if parralelization:
            options['maxiter'] = 1
            options['miniter'] = 1
        if not Pdirs:
            Pdirs = self.subproc
        return options
"""


class SurveyRepairTests(unittest.TestCase):
    def make_process(self, process, *, loop=False):
        source = process / "bin/internal/gen_ximprove.py"
        source.parent.mkdir(parents=True)
        source.write_text(SOURCE)
        characteristics = process / "SubProcesses/proc_characteristics"
        characteristics.parent.mkdir()
        characteristics.write_text(f"loop_induced = {loop}\n")
        (process / "Cards").mkdir()
        for name in ("run_card.dat", "run_card_default.dat"):
            (process / "Cards" / name).write_text("  1 = nhel ! sample\n 6800 = ebeam1\n")
        return source

    def test_invalid_unsplit_requests_are_raised_to_minimum_and_valid_chunks_stay_one(self):
        patched = guarded_source(SOURCE)
        namespace = {"logger": Mock()}
        exec(compile(patched, "guard_test", "exec"), namespace)
        survey = namespace["Survey"]()
        survey.min_iterations = 3
        survey.subproc = []
        for iterations in (1, 2, 3, 7):
            survey.iterations = iterations
            self.assertEqual(survey.write_parameter(False),
                             {"maxiter": max(iterations, 3), "miniter": 3})
            self.assertEqual(survey.write_parameter(True), {"maxiter": 1, "miniter": 1})
        self.assertEqual(namespace["logger"].warning.call_count, 2)
        self.assertEqual(guarded_source(patched), patched)

    def test_unknown_source_is_refused(self):
        for source in ("class Unrelated: pass\n", SOURCE.replace("options['miniter'] = 1", "pass")):
            with self.assertRaisesRegex(ValueError, "Unsupported MG5"):
                guarded_source(source)

    def test_repair_retains_backup_receipt_and_is_idempotent(self):
        with tempfile.TemporaryDirectory() as tmp:
            process = Path(tmp).resolve()
            source = self.make_process(process, loop=True)
            event = process / "saved_event.lhe"
            event.write_text("preserve")
            receipt_path = repair_survey(process)
            receipt = json.loads(receipt_path.read_text())
            self.assertEqual(len(receipt["files"]), 3)
            for item in receipt["files"]:
                self.assertEqual(item["old_sha256"], hashlib.sha256(Path(item["backup"]).read_bytes()).hexdigest())
                self.assertEqual(item["new_sha256"], hashlib.sha256(Path(item["source"]).read_bytes()).hexdigest())
            self.assertEqual(Path(receipt["files"][0]["backup"]).read_text(), SOURCE)
            self.assertEqual((process / "Cards/run_card.dat").read_text(), "  0 = nhel ! sample\n 6800 = ebeam1\n")
            self.assertEqual(event.read_text(), "preserve")
            mtime = source.stat().st_mtime_ns
            self.assertIsNone(repair_survey(process))
            self.assertEqual(source.stat().st_mtime_ns, mtime)
            self.assertEqual(len(list((process / ".trsm-maintenance").iterdir())), 1)
            self.assertEqual(pending_updates(process), {})

    def test_tree_repair_preserves_helicity_card(self):
        with tempfile.TemporaryDirectory() as tmp:
            process = Path(tmp)
            self.make_process(process)
            repair_survey(process)
            self.assertEqual((process / "Cards/run_card.dat").read_text(), "  1 = nhel ! sample\n 6800 = ebeam1\n")

    def test_exact_helicity_card_changes_only_sampling_parameters(self):
        source = ("# 1 = nhel in a comment\n  1 = nhel ! sampling\n\t1 = nhel_survey\n"
                  " 0 = nhel_refine\n 21 = iseed\n 6800 = ebeam1\n nn23lo1 = pdlabel\n")
        expected = source.replace("  1 = nhel", "  0 = nhel").replace("\t1 = nhel", "\t0 = nhel")
        self.assertEqual(exact_helicity_card(source), expected)
        self.assertEqual(exact_helicity_card(expected), expected)
        for invalid in ("", " 2 = nhel\n", " 1 = nhel\n 0 = nhel\n"):
            with self.assertRaisesRegex(ValueError, "Unsupported MG5 run card"):
                exact_helicity_card(invalid)

    def test_invalid_card_aborts_before_changing_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            process = Path(tmp)
            source = self.make_process(process, loop=True)
            (process / "Cards/run_card.dat").write_text("unrecognized card\n")
            with self.assertRaises(ValueError):
                repair_survey(process)
            self.assertEqual(source.read_text(), SOURCE)
            self.assertFalse((process / ".trsm-maintenance").exists())

    def test_active_campaign_writer_prevents_repair(self):
        with tempfile.TemporaryDirectory() as tmp:
            process = Path(tmp).resolve()
            source = self.make_process(process)
            with _madevent_lock(process):
                with self.assertRaisesRegex(ValueError, "in use"):
                    repair_survey(process)
            self.assertEqual(source.read_text(), SOURCE)
            self.assertFalse((process / ".trsm-maintenance").exists())

    def test_external_generator_source_link_is_not_modified(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            external = root / "external.py"
            external.write_text(SOURCE)
            source = self.make_process(root / "process")
            source.unlink()
            source.symlink_to(external)
            with self.assertRaisesRegex(ValueError, "outside the selected process"):
                repair_survey(root / "process")
            self.assertEqual(external.read_text(), SOURCE)


if __name__ == "__main__":
    unittest.main()
