import json
import copy
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import textwrap
import unittest
from unittest.mock import patch

from generate_mg5_trsm_xsecs import mg5_runtime_receipt
from tools.setup_mg5_pool import setup_pool
from tools.migrate_mg5_pool_campaign import validate_transition
from trsm_mg5_pool import process_workspace


class TestMG5Pool(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.source = Path(self.tmp.name).resolve() / "source"
        self.pool = Path(self.tmp.name).resolve() / "pool"
        self.process = self.source / "pp_eta0Z"
        (self.process / "bin").mkdir(parents=True)
        (self.source / "VERSION").write_text("test version")
        sub = self.process / "SubProcesses"
        (sub / "P1").mkdir(parents=True)
        (sub / "subproc.mg").write_text("P1\n")
        (sub / "proc_characteristics").write_text("loop_induced = False\n")
        (sub / "P1/madevent").write_text("binary fixture")
        (sub / "P1/madevent").chmod(0o755)
        (self.process / "Cards").mkdir()
        (self.process / "Cards/proc_card_mg5.dat").write_text("unchanged physics")
        madevent = self.process / "bin/madevent"
        madevent.write_text(textwrap.dedent('''\
            #!/usr/bin/env python3
            import gzip,json,sys,time
            from pathlib import Path
            text=Path(sys.argv[1]).read_text()
            run=next(x.split()[1] for x in text.splitlines() if x.startswith('generate_events '))
            mass=next(x.split()[2] for x in text.splitlines() if x.startswith('set Meta '))
            guard=Path('active')
            with guard.open('x') as f:f.write(run)
            card=Path('Cards/param_card.dat');card.write_text(mass)
            with Path('timeline.jsonl').open('a') as f:f.write(json.dumps(['start',run,time.time()])+'\\n')
            time.sleep(.4)
            assert card.read_text()==mass
            lhe=Path('Events')/run/'unweighted_events.lhe.gz';lhe.parent.mkdir(parents=True)
            with gzip.open(lhe,'wt') as f:f.write('Integrated weight : '+mass+'\\n')
            with Path('timeline.jsonl').open('a') as f:f.write(json.dumps(['end',run,time.time()])+'\\n')
            guard.unlink()
            '''))
        madevent.chmod(0o755)

    def setup_pool(self, slots=2):
        return setup_pool(self.source, self.pool, ["pp_eta0Z"], slots)

    def test_copy_is_independent_and_omits_old_integration_output(self):
        (self.process / "SubProcesses/P1/G1").mkdir()
        (self.process / "SubProcesses/P1/G1/large-old-output").write_text("old")
        (self.process / "Cards/relative").symlink_to("proc_card_mg5.dat")
        (self.process / "Cards/absolute").symlink_to(self.process / "Cards/proc_card_mg5.dat")
        self.setup_pool()
        first = self.pool / "slots/000/pp_eta0Z"
        second = self.pool / "slots/001/pp_eta0Z"
        self.assertFalse((first / "SubProcesses/P1/G1").exists())
        self.assertEqual((first / "Cards/absolute").resolve(), first / "Cards/proc_card_mg5.dat")
        (first / "Cards/proc_card_mg5.dat").write_text("private change")
        self.assertEqual((second / "Cards/proc_card_mg5.dat").read_text(), "unchanged physics")
        self.assertEqual((self.process / "Cards/proc_card_mg5.dat").read_text(), "unchanged physics")

    def test_preflight_checks_every_copy_and_source(self):
        self.setup_pool()
        with patch.dict(os.environ, TRSM_MG5_POOL_ROOT=str(self.pool)):
            receipt = mg5_runtime_receipt(["pp_eta0Z"], self.source)
            self.assertEqual(receipt["pool"]["slots_per_process"], 2)
            (self.pool / "slots/001/pp_eta0Z/Cards/proc_card_mg5.dat").write_text("wrong model")
            with self.assertRaisesRegex(ValueError, "mismatch"):
                mg5_runtime_receipt(["pp_eta0Z"], self.source)

    def test_no_pool_preserves_reference_location(self):
        with patch.dict(os.environ, {"TRSM_MG5_POOL_ROOT": ""}):
            with process_workspace("pp_eta0Z", str(self.source), {"pp_eta0Z": "pp_eta0Z"}) as path:
                self.assertEqual(path, str(self.source))

    def test_concurrent_jobs_use_both_copies_without_mixing_cards_or_results(self):
        self.setup_pool()
        code = ('from generate_mg5_trsm_xsecs import get_mg5_xsec; import sys; '
                'print("RESULT",get_mg5_xsec("pp_eta0Z",sys.argv[1],list(range(1,11)),'
                '.99,.1,0,float(sys.argv[2]),1.,50.,0.,ecm=13.6,w1=.004,k233=2.))')
        env = dict(os.environ, TRSM_MG5_LOCATION=str(self.source), TRSM_MG5_POOL_ROOT=str(self.pool),
                   TRSM_MG5_CORES="1", PYTHONPATH=str(Path(__file__).resolve().parent))
        children = [subprocess.Popen([sys.executable, "-B", "-c", code, f"test{i}", str(300+i)],
                                     env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                    for i in range(6)]
        try:
            for i, child in enumerate(children):
                output, _ = child.communicate(timeout=20)
                self.assertEqual(child.returncode, 0, output)
                self.assertIn(f"RESULT {300+i}.0", output)
        finally:
            for child in children:
                if child.poll() is None:child.kill()
                child.wait()
        spans = []
        runs = []
        for slot in ("000", "001"):
            path = self.pool / f"slots/{slot}/pp_eta0Z/timeline.jsonl"
            self.assertTrue(path.exists(), "Every copy should be usable")
            rows = [json.loads(x) for x in path.read_text().splitlines()]
            self.assertEqual([r[0] for r in rows], ["start", "end"]*(len(rows)//2))
            for start, end in zip(rows[::2], rows[1::2]):
                self.assertEqual(start[1], end[1]);runs.append(start[1])
                spans.append((slot, start[2], end[2]))
        self.assertEqual(len(runs), 6)
        self.assertEqual(len(set(runs)), 6)
        self.assertTrue(any(a[0]!=b[0] and max(a[1],b[1])<min(a[2],b[2])
                            for a in spans for b in spans))

    def test_exception_releases_slot(self):
        self.setup_pool(slots=1)
        with patch.dict(os.environ, TRSM_MG5_POOL_ROOT=str(self.pool)):
            with self.assertRaisesRegex(RuntimeError, "injected"):
                with process_workspace("pp_eta0Z", self.source, {"pp_eta0Z": "pp_eta0Z"}):
                    raise RuntimeError("injected")
            with process_workspace("pp_eta0Z", self.source, {"pp_eta0Z": "pp_eta0Z"}) as path:
                self.assertEqual(path, self.pool / "slots/000")


class TestMG5PoolTransition(unittest.TestCase):
    def fixtures(self):
        old = {"physics_manifest": {"sources_sha256": {"generate_mg5_trsm_xsecs.py": "old", "physics.py": "fixed"},
                                    "sha256": "oldmanifest", "settings": {"energy": 13.6}},
               "mg5_runtime": {"sources_sha256": {"model": "fixed"}, "cores": "1"},
               "generator_sha256": "fixed-generator"}
        new = copy.deepcopy(old)
        new["physics_manifest"]["sources_sha256"]["generate_mg5_trsm_xsecs.py"] = "new"
        new["physics_manifest"]["sha256"] = "newmanifest"
        new["mg5_runtime"]["pool"] = {"root": "/pool", "slots_per_process": 16}
        changes = {"generate_mg5_trsm_xsecs.py": {"before": "old", "after": "new"}}
        return old, new, changes

    def test_exact_operational_change_is_accepted(self):
        old, new, changes = self.fixtures()
        self.assertEqual(validate_transition(old, new, changes)["slots_per_process"], 16)

    def test_physics_or_unlisted_source_change_is_rejected(self):
        for kind in ("physics", "source", "runtime"):
            old, new, changes = self.fixtures()
            if kind == "physics":new["physics_manifest"]["settings"]["energy"] = 14
            if kind == "source":new["physics_manifest"]["sources_sha256"]["physics.py"] = "changed"
            if kind == "runtime":new["mg5_runtime"]["sources_sha256"]["model"] = "changed"
            with self.assertRaises(ValueError):validate_transition(old, new, changes)


if __name__ == "__main__":
    unittest.main()
