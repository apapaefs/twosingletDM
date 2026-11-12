import contextlib
import io
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import textwrap
import unittest
from unittest.mock import patch

from generate_mg5_trsm_xsecs import _madevent_lock
from tools.madloop_build import (MadLoopBuildError, build_madloop_checks, main,
                                polynomial_modules, repair_madloop, validate_madloop)


class MadLoopBuildTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.process = Path(self.temporary.name).resolve() / "process"
        self.root = self.process / "SubProcesses"
        self.directory = self.root / "PV7_2_9_test"
        self.directory.mkdir(parents=True)
        (self.root / "proc_characteristics").write_text("loop_induced = True\n")
        self.module = self.directory / "ml5_7_2_9_polynomial_constants.mod"
        (self.directory / "polynomial.f").write_text(textwrap.dedent("""\
                  MODULE ML5_7_2_9_POLYNOMIAL_CONSTANTS
                  INTEGER, PARAMETER :: VALUE = 7
                  END MODULE ML5_7_2_9_POLYNOMIAL_CONSTANTS
            """))
        # The actual MG5 export is fixed-form Fortran.
        source = self.directory / "polynomial.f"
        source.write_text("".join("      " + line + "\n" for line in source.read_text().splitlines()))
        (self.directory / "loop_matrix.f").write_text(
            "      SUBROUTINE MATRIX(ANSWER)\n"
            "      USE ML5_7_2_9_POLYNOMIAL_CONSTANTS\n"
            "      INTEGER ANSWER\n"
            "      ANSWER=VALUE\n"
            "      END\n")
        (self.directory / "check_sa.f").write_text(
            "      PROGRAM CHECK\n"
            "      INTEGER ANSWER\n"
            "      CALL MATRIX(ANSWER)\n"
            "      IF (ANSWER.NE.7) STOP 1\n"
            "      END\n")
        (self.directory / "makefile").write_text(textwrap.dedent("""\
            FC = gfortran
            check: check_sa.o polynomial.o loop_matrix.o
            \t$(FC) -o $@ $^
            polynomial.o: polynomial.f
            \t$(FC) -c $< -o $@
            loop_matrix.o: loop_matrix.f polynomial.o
            \t$(FC) -c $< -o $@
            check_sa.o: check_sa.f
            \t$(FC) -c $< -o $@
            """))

    def ready(self):
        check = self.directory / "check"
        check.write_text("compiled executable")
        check.chmod(0o755)
        self.module.write_text("compiled module")

    def test_preflight_rejects_parent_only_module_without_writing(self):
        self.ready()
        self.module.rename(self.root / self.module.name)
        before = set(self.process.rglob("*"))
        with self.assertRaisesRegex(MadLoopBuildError, r"Missing MadLoop.*--repair"):
            validate_madloop(self.process)
        self.assertEqual(set(self.process.rglob("*")), before)

    def test_empty_module_and_nonexecutable_check_are_rejected(self):
        self.ready()
        self.module.write_bytes(b"")
        with self.assertRaisesRegex(MadLoopBuildError, self.module.name):
            validate_madloop(self.process)
        self.ready()
        (self.directory / "check").chmod(0o644)
        with self.assertRaisesRegex(MadLoopBuildError, "/check"):
            validate_madloop(self.process)

    def test_module_names_come_from_source(self):
        self.assertEqual(polynomial_modules(self.directory), [self.module])
        with (self.directory / "polynomial.f").open("a") as stream:
            stream.write("      module Different_Prefix\n      end module Different_Prefix\n")
        self.assertEqual(polynomial_modules(self.directory),
                         [self.module, self.directory / "different_prefix.mod"])

    def test_unrecognized_polynomial_source_is_an_error(self):
        (self.directory / "polynomial.f").write_text("C incomplete export\n")
        with self.assertRaisesRegex(MadLoopBuildError, "No Fortran module"):
            validate_madloop(self.process)

    def test_nonoptimized_and_tree_exports_remain_supported(self):
        self.ready()
        (self.directory / "polynomial.f").unlink()
        self.assertEqual(validate_madloop(self.process), [self.directory / "check"])
        self.directory.rename(self.root / "unused")
        (self.root / "proc_characteristics").write_text("loop_induced = False\n")
        self.assertEqual(validate_madloop(self.process), [])
        self.assertEqual(repair_madloop(self.process), [])
        self.assertFalse((self.process / ".trsm_madevent.lock").exists())

    def test_missing_loop_directory_is_not_mistaken_for_tree_export(self):
        self.directory.rename(self.root / "unused")
        with self.assertRaisesRegex(MadLoopBuildError, "no MadLoop subprocesses"):
            validate_madloop(self.process)

    def test_repair_rejects_out_of_process_subdirectory(self):
        outside = Path(self.temporary.name) / "outside"
        self.directory.rename(outside)
        self.directory.symlink_to(outside, target_is_directory=True)
        with self.assertRaisesRegex(MadLoopBuildError, "leaves the selected process"):
            repair_madloop(self.process)
        self.assertFalse((self.process / ".trsm_madevent.lock").exists())

    def test_repair_uses_campaign_writer_lock_and_refuses_active_process(self):
        with _madevent_lock(self.process), patch("tools.madloop_build.subprocess.run") as run:
            with self.assertRaisesRegex(MadLoopBuildError, "in use"):
                repair_madloop(self.process)
            run.assert_not_called()
        self.assertEqual(list(self.process.glob("madloop-build-*.log")), [])

    def test_failed_build_preserves_log_and_unrelated_data(self):
        event = self.process / "saved_event.lhe"
        event.write_text("keep")
        with contextlib.redirect_stdout(io.StringIO()), patch("tools.madloop_build.subprocess.run",
                return_value=subprocess.CompletedProcess([], 2)):
            with self.assertRaisesRegex(MadLoopBuildError, "compilation failed.*exit 2"):
                repair_madloop(self.process)
        log, = self.process.glob("madloop-build-*.log")
        self.assertIn("-W polynomial.f check", log.read_text())
        self.assertEqual(event.read_text(), "keep")
        self.assertFalse(self.module.exists())

    def test_zero_exit_without_artifacts_is_not_success(self):
        with self.assertRaisesRegex(MadLoopBuildError, "Missing MadLoop"):
            build_madloop_checks(self.process, lambda command: None)

    def test_cli_check_reports_failure_and_success(self):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(main([str(self.process)]), 1)
            self.ready()
            self.assertEqual(main([str(self.process)]), 0)

    @unittest.skipUnless(shutil.which("gfortran") and shutil.which("make"), "requires native Fortran and make")
    def test_native_parent_build_repair_and_repeated_madloop_initialization(self):
        sources = {p: p.read_bytes() for p in self.directory.iterdir()}
        # Match OLP_static: the object goes in PV*/, the module in its parent.
        subprocess.run(["gfortran", "-c", str(self.directory / "polynomial.f"),
                        "-o", str(self.directory / "polynomial.o")], cwd=self.root,
                       check=True, capture_output=True)
        self.assertTrue((self.root / self.module.name).is_file())
        broken = subprocess.run(["make", "-j1", "check"], cwd=self.directory,
                                text=True, capture_output=True)
        self.assertNotEqual(broken.returncode, 0)
        self.assertIn(self.module.name, broken.stderr)

        with contextlib.redirect_stdout(io.StringIO()):
            artifacts = repair_madloop(self.process)
        self.assertIn(self.module, artifacts)
        subprocess.run([str(self.directory / "check")], check=True)

        # MadLoop removes these files and recompiles at each initialization.
        for name in ("check", "check_sa.o", "loop_matrix.o"):
            (self.directory / name).unlink()
        subprocess.run(["make", "-j1", "check"], cwd=self.directory, check=True, capture_output=True)
        subprocess.run([str(self.directory / "check")], check=True)

        mtimes = {p: p.stat().st_mtime_ns for p in artifacts}
        with contextlib.redirect_stdout(io.StringIO()), patch.dict(os.environ,
                {"MAKEFLAGS": "-f nonexistent -j8", "GNUMAKEFLAGS": "-f nonexistent"}):
            repair_madloop(self.process)
        self.assertEqual(mtimes, {p: p.stat().st_mtime_ns for p in artifacts})
        # Recover even when check and polynomial.o exist but a module is lost.
        self.module.unlink()
        with contextlib.redirect_stdout(io.StringIO()):
            repair_madloop(self.process)
        subprocess.run([str(self.directory / "check")], check=True)
        self.assertTrue(self.module.is_file())
        self.assertEqual(sources, {p: p.read_bytes() for p in sources})


if __name__ == "__main__":
    unittest.main()
