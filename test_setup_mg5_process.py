import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from tools.bootstrap_runtime import SetupError, prepare_madloop_rpaths
from tools.setup_mg5_process import (
    ROOT, generation_card, model_hashes, setup_process, validate_definition,
)


class SetupMG5ProcessTests(unittest.TestCase):
    def setUp(self):
        self.spec = json.loads((ROOT / "config/runtime-sources-v2.json").read_text())["mg5_processes"]["gg_eta0Z"]

    def test_generation_card_does_not_save_global_settings(self):
        card = generation_card(self.spec, jobs=1, fc="gfortran", cxx="c++", collier="/runtime/collier")
        self.assertIn("generate g g > eta0 z [noborn=QCD]\n", card)
        self.assertIn("output gg_eta0Z\n", card)
        self.assertTrue(all(line.endswith(" --no_save") for line in card.splitlines() if line.startswith("set ")))

    def fixture(self, root):
        (root / "bin").mkdir()
        (root / "bin/mg5_aMC").touch()
        (root / "VERSION").write_text("3.5.15")
        model = root / "models/loop_sm_twoscalar_generic"
        model.mkdir(parents=True)
        (model / "parameters.py").write_text("# UFO source\n")
        process = root / "gg_eta0Z"
        (process / "Cards").mkdir(parents=True)
        (process / "Cards/proc_card_mg5.dat").write_text(
            "import model loop_sm_twoscalar_generic\n"
            "generate g g > eta0 z [noborn=QCD]\noutput gg_eta0Z\n")
        exported = process / "bin/internal/ufomodel"
        exported.mkdir(parents=True)
        (exported / "parameters.py").write_text((model / "parameters.py").read_text())
        return process, model

    def test_matching_existing_process_is_validated_without_generation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            process, model = self.fixture(root)
            original = (process / "Cards/proc_card_mg5.dat").read_bytes()
            with patch("tools.setup_mg5_process.mg5_runtime_receipt", return_value={"checked": True}) as receipt, \
                    patch("tools.setup_mg5_process.subprocess.run") as native:
                self.assertEqual(setup_process(root, "gg_eta0Z"), {"checked": True})
                native.assert_not_called()
                receipt.assert_called_once_with(["gg_eta0Z"], mgloc=root.resolve())
            self.assertEqual((process / "Cards/proc_card_mg5.dat").read_bytes(), original)

    def test_conflicting_process_and_model_are_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            process, model = self.fixture(root)
            card = process / "Cards/proc_card_mg5.dat"
            original = card.read_text()
            card.write_text(original + "add process g g > h h [noborn=QCD]\n")
            with self.assertRaisesRegex(SetupError, "different process/model"):
                setup_process(root, "gg_eta0Z")
            card.write_text(original)
            (model / "parameters.py").write_text("# changed model\n")
            with self.assertRaisesRegex(SetupError, "UFO differs"):
                setup_process(root, "gg_eta0Z")
            self.assertEqual(card.read_text(), original)

    def test_recorded_identity_change_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            process, model = self.fixture(root)
            (process / "trsm-process-setup.json").write_text(json.dumps({"identity": {}}))
            with self.assertRaisesRegex(SetupError, "source identity differs"):
                setup_process(root, "gg_eta0Z")

    def test_model_restrictions_are_part_of_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "parameters.py").write_text("# model")
            before = model_hashes(root)
            (root / "restrict_default.dat").write_text("# restriction")
            self.assertNotEqual(before, model_hashes(root))

    def test_madloop_check_keeps_runtime_library_paths_in_fresh_export(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            original = root / "template"
            original.write_text("LINKLIBS = -lmodel $(LINK_LOOP_LIBS) $(LDFLAGS)\n")
            (root / "SubProcesses").mkdir()
            makefile = root / "SubProcesses/makefile_MadLoop"
            makefile.symlink_to(original)
            prepare_madloop_rpaths(root)
            prepare_madloop_rpaths(root)
            self.assertFalse(makefile.is_symlink())
            self.assertEqual(makefile.read_text().count("$(RPATH_LIBS)"), 1)
            self.assertNotIn("RPATH_LIBS", original.read_text())


if __name__ == "__main__":
    unittest.main()
