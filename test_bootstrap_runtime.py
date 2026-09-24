import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

from tools.bootstrap_runtime import (Installer, SetupError, activation_text,
    check_bsmpt_startup, component_selection, download_archive, extract_source, madgraph_subprocesses,
    main, parse_args, process_card)
from trsm_paths import higgs_dataset_path


class TestSources(unittest.TestCase):
    def test_download_and_cache_require_matching_digest(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "upstream.tgz"
            source.write_bytes(b"test source release")
            spec = {"filename": "release.tgz", "url": source.as_uri(),
                    "sha256": hashlib.sha256(source.read_bytes()).hexdigest()}
            result = download_archive(spec, root / "cache")
            source.unlink()
            self.assertEqual(download_archive(spec, root / "cache"), result)
            result.write_bytes(b"corrupted cache")
            with self.assertRaisesRegex(SetupError, "Cached archive checksum mismatch"):
                download_archive(spec, root / "cache")

    def test_bad_download_is_never_published(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "upstream.tgz"
            source.write_bytes(b"different release")
            spec = {"filename": "release.tgz", "url": source.as_uri(), "sha256": "0" * 64}
            with self.assertRaisesRegex(SetupError, "Checksum mismatch"):
                download_archive(spec, root / "cache")
            self.assertEqual(list((root / "cache").iterdir()), [])

    def archive(self, root, entries):
        archive = root / "source.tgz"
        with tarfile.open(archive, "w:gz") as stream:
            for name, content in entries:
                member = tarfile.TarInfo(name)
                member.size = len(content)
                stream.addfile(member, io.BytesIO(content))
        return archive

    def test_ufo_extraction_omits_machine_caches_and_preserves_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive = self.archive(root, [("model/parameters.py", b"source"),
                ("model/model.pkl", b"cache"), ("model/._parameters.py", b"macOS metadata"),
                ("model/__pycache__/parameters.pyc", b"cache")])
            extract_source(archive, root / "model", "model", ufo=True)
            self.assertEqual((root / "model/parameters.py").read_bytes(), b"source")
            self.assertEqual([p.name for p in (root / "model").iterdir()], ["parameters.py"])
            with self.assertRaisesRegex(SetupError, "Refusing to extract over"):
                extract_source(archive, root / "model", "model")

    def test_archive_traversal_cannot_escape_destination(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive = self.archive(root, [("../../escaped", b"bad")])
            with self.assertRaises(tarfile.FilterError):
                extract_source(archive, root / "model", "model")
            self.assertFalse((root / "model").exists())
            self.assertFalse((root / "escaped").exists())


class TestInstaller(unittest.TestCase):
    def installer(self, prefix, *options):
        return Installer(parse_args(["--prefix", str(prefix), "--components", "python", *options]))

    def test_dry_run_never_creates_prefix(self):
        with tempfile.TemporaryDirectory() as temporary, patch("sys.stdout", new_callable=io.StringIO):
            prefix = Path(temporary) / "new"
            self.assertEqual(main(["--prefix", str(prefix), "--dry-run"]), 0)
            self.assertFalse(prefix.exists())

    def test_existing_installation_is_untouched(self):
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary)
            existing = prefix / "main"
            existing.write_text("active runtime")
            with self.assertRaisesRegex(SetupError, "unmanaged runtime"):
                with self.installer(prefix).workspace():
                    self.fail("Must not enter existing installation")
            self.assertEqual(existing.read_text(), "active runtime")
            self.assertEqual(list(prefix.iterdir()), [existing])

    def test_successful_step_is_not_repeated_and_changed_artifact_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary) / "runtime"
            first = self.installer(prefix)
            calls = []
            def build():
                calls.append("built")
                artifact = prefix / "executable"
                artifact.write_text("correct binary")
                return [artifact]
            with first.workspace():
                first.step("native", build)
            second = self.installer(prefix, "--resume")
            with second.workspace():
                second.step("native", build)
            self.assertEqual(calls, ["built"])
            (prefix / "executable").write_text("changed binary")
            third = self.installer(prefix, "--resume")
            with third.workspace(), self.assertRaisesRegex(SetupError, "artifact changed"):
                third.step("native", build)

    def test_failed_build_is_preserved_on_resume(self):
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary) / "runtime"
            first = self.installer(prefix)
            def failing():
                path = first.fresh(prefix / "build/model")
                path.mkdir()
                (path / "diagnostic").write_text("partial build evidence")
                raise SetupError("compiler failed")
            with first.workspace(), self.assertRaisesRegex(SetupError, "compiler failed"):
                first.step("native", failing)
            second = self.installer(prefix, "--resume")
            def succeed():
                path = second.fresh(prefix / "build/model")
                path.mkdir()
                (path / "main").write_text("new executable")
                return [path / "main"]
            with second.workspace():
                second.step("native", succeed)
            old = list((prefix / "failed-builds").glob("model-*/diagnostic"))
            self.assertEqual(len(old), 1)
            self.assertEqual(old[0].read_text(), "partial build evidence")
            self.assertEqual(json.loads((prefix / "bootstrap-state.json").read_text())["steps"]["native"]["status"], "done")

    def test_changed_build_configuration_rejects_resume(self):
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary) / "runtime"
            with self.installer(prefix).workspace():
                pass
            changed = self.installer(prefix, "--resume")
            changed.compilers["FC"] = "different-gfortran"
            with self.assertRaisesRegex(SetupError, "Setup inputs changed"):
                with changed.workspace():
                    pass

    def test_download_only_can_resume_into_build(self):
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary) / "runtime"
            first = self.installer(prefix, "--download-only")
            with patch.object(first, "preflight"), patch("sys.stdout", new_callable=io.StringIO):
                first.install()
            self.assertEqual(first.state["steps"], {})
            resumed = self.installer(prefix, "--resume")
            def built():
                artifact = prefix / "installed"
                artifact.write_text("built after download-only")
                return [artifact]
            with patch.object(resumed, "preflight"), patch.object(resumed, "python_environment", side_effect=built), patch.object(resumed, "verify", return_value=[]):
                resumed.install()
            self.assertEqual(resumed.state["steps"]["python"]["status"], "done")
            self.assertEqual(resumed.state["steps"]["verify"]["status"], "done")

    def test_concurrent_install_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary) / "runtime"
            with self.installer(prefix).workspace():
                with self.assertRaisesRegex(SetupError, "Another installer"):
                    with self.installer(prefix, "--resume").workspace():
                        pass

    def test_process_plan_matches_scan_interface(self):
        from generate_mg5_trsm_xsecs import ProcLocation
        with tempfile.TemporaryDirectory() as temporary:
            installer = self.installer(Path(temporary) / "runtime")
            for key, spec in installer.config["mg5_processes"].items():
                self.assertEqual(spec["directory"], ProcLocation[key].rstrip("/"))
                card = process_card(spec, jobs=2, fc="gfortran", cxx="c++", collier="/runtime/collier")
                self.assertIn("generate " + spec["generate"], card)
                self.assertNotIn("launch", card)
                self.assertNotIn("generate_events", card)

    def test_zero_exit_without_generated_process_is_failure(self):
        with tempfile.TemporaryDirectory() as temporary:
            installer = self.installer(Path(temporary) / "runtime")
            (installer.prefix / "logs").mkdir(parents=True)
            with patch.object(installer, "run"):
                with self.assertRaisesRegex(SetupError, "did not generate"):
                    installer.process("pp_eta0Z")

    def test_madgraph_staging_preserves_deferred_template_links(self):
        with tempfile.TemporaryDirectory() as temporary:
            installer = self.installer(Path(temporary) / "runtime")
            source = installer.prefix / "sources/madgraph"
            (source / "bin").mkdir(parents=True)
            (source / "input").mkdir()
            (source / "bin/mg5_aMC").write_text("launcher")
            (source / "bin/create_release.py").write_text("upstream packaging utility")
            (source / "VERSION").write_text("3.5.15")
            (source / "input/mg5_configuration.txt").touch()
            (source / "template.inc").symlink_to("generated_later.inc")
            installer.madgraph()
            link = installer.prefix / "MG5_aMC_v3_5_15/template.inc"
            self.assertTrue(link.is_symlink())
            self.assertEqual(os.readlink(link), "generated_later.inc")
            self.assertFalse((installer.prefix / "MG5_aMC_v3_5_15/bin/create_release.py").exists())
            self.assertTrue((source / "bin/create_release.py").exists())

    def test_only_integration_subprocesses_receive_madevent_target(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name in ("P0_gg_heta0", "PV0_0_1_gg_heta0"):
                path = root / "SubProcesses" / name
                path.mkdir(parents=True)
                (path / "makefile").touch()
            listing = root / "SubProcesses/subproc.mg"
            listing.write_text("P0_gg_heta0\n")
            self.assertEqual(madgraph_subprocesses(root), [root / "SubProcesses/P0_gg_heta0"])
            listing.write_text("../outside\n")
            with self.assertRaisesRegex(SetupError, "Invalid generated subprocess"):
                madgraph_subprocesses(root)

    def test_collier_installation_has_madgraph_module_layout(self):
        with tempfile.TemporaryDirectory() as temporary:
            installer = self.installer(Path(temporary) / "runtime")
            def extracted(archive, target, directory):
                (target / "modules").mkdir(parents=True)
                (target / "modules/collier.mod").write_bytes(b"Fortran module")
                (target / "libcollier.a").write_bytes(b"static library")
            with patch("tools.bootstrap_runtime.extract_source", side_effect=extracted), patch.object(installer, "run"):
                artifacts = installer.collier()
            self.assertIn(installer.prefix / "collier/modules/collier.mod", artifacts)
            self.assertEqual((installer.prefix / "collier/libcollier.a").read_bytes(), b"static library")

    def test_bsmpt_help_exit_one_requires_expected_help_and_model(self):
        help_text = "CalcTemps calculates characteristic temperatures\nThe implemented models are\ntrsm\nsm\n"
        for status, output, valid in (
                (0, help_text, True),
                (1, help_text + "Not all required parameters are set.\n", True),
                (1, help_text + "Unexpected initialization failure\n", False),
                (0, help_text.replace("trsm\n", ""), False),
                (-11, help_text, False)):
            with self.subTest(status=status, output=output), patch("tools.bootstrap_runtime.subprocess.run",
                    return_value=subprocess.CompletedProcess([], status, output)):
                if valid:
                    check_bsmpt_startup("CalcTemps", {}, io.StringIO())
                else:
                    with self.assertRaisesRegex(SetupError, "BSMPT startup check failed"):
                        check_bsmpt_startup("CalcTemps", {}, io.StringIO())

    def test_component_dependencies_and_activation_paths(self):
        self.assertEqual(component_selection(["higgstools"]), ("python", "datasets", "higgstools"))
        text = activation_text(Path("/tmp/test runtime"), ("python", "datasets", "mg5"))
        self.assertIn("export TRSM_RUNTIME_ROOT='/tmp/test runtime'", text)
        self.assertIn("export TRSM_HB_DATASET=", text)
        self.assertIn("export TRSM_MG5_LOCATION=", text)
        self.assertNotIn("export HOME=", text)

    def test_serial_builds_ignore_inherited_make_parallelism(self):
        with tempfile.TemporaryDirectory() as temporary, patch.dict(os.environ, {"MAKEFLAGS": "-j8", "MFLAGS": "-j8", "MAKELEVEL": "2"}):
            installer = self.installer(Path(temporary) / "runtime")
            for key in ("MAKEFLAGS", "MFLAGS", "MAKELEVEL"):
                self.assertNotIn(key, installer.env)

    def test_higgs_dataset_override_matches_provider_and_provenance(self):
        with tempfile.TemporaryDirectory() as temporary, patch.dict(os.environ, {"TRSM_HB_DATASET": temporary}):
            self.assertEqual(higgs_dataset_path("hbdataset"), Path(temporary).resolve())
        with patch.dict(os.environ, {}, clear=True):
            self.assertEqual(higgs_dataset_path("hsdataset").name, "hsdataset")

    @unittest.skipUnless(shutil.which("rsync"), "BSMPT staging requires rsync")
    def test_bsmpt_staging_uses_upstream_filename_case_and_shared_inputs(self):
        root = Path(__file__).resolve().parent
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "upstream"
            (source / "src/models").mkdir(parents=True)
            (source / "src/prog").mkdir()
            (source / "CMakeLists.txt").write_text("set(BSMPT_VERSION 3.2.1)\n")
            (source / "src/CMakeLists.txt").touch()
            fields = {"C_MassW": "MW_GeV", "C_MassZ": "MZ_GeV", "C_MassSMHiggs": "M1_GeV", "C_GF": "GF_GeV^-2"}
            (source / "src/models/SMParam.cpp").write_text("\n".join(f"SM.{name} = 0;" for name in fields))
            for name in ("CalcTemps", "MinimaTracer"):
                (source / f"src/prog/{name}.cpp").write_text("std::ofstream outfile(filename);\n")
            target = Path(temporary) / "staged"
            subprocess.run([sys.executable, str(root / "tools/stage_bsmpt_v2.py"), str(source), str(target)],
                           check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            inputs = json.loads((root / "config/sm-inputs-v2.json").read_text())
            values = (target / "src/models/SMParam.cpp").read_text()
            for name, key in fields.items():
                self.assertIn(f"SM.{name} = {inputs[key]:.17g};", values)
            self.assertIn("std::setprecision(17)", (target / "src/prog/CalcTemps.cpp").read_text())


if __name__ == "__main__":
    unittest.main()
