import unittest
from contextlib import redirect_stdout
from io import StringIO

from mg5_process_runner import run_mg5_processes, selected_mg5_processes


class TestMG5ProcessRunner(unittest.TestCase):
    def test_selected_processes_skips_unavailable_with_warning(self):
        buffer = StringIO()

        with redirect_stdout(buffer):
            selected = selected_mg5_processes(
                ["hh", "missing", "hhh"],
                proc_location={"hh": "gg_hh_twoscalar/", "hhh": "gg_hhh_twoscalar/"},
            )

        self.assertEqual(selected, ["hh", "hhh"])
        self.assertIn("Warning: MG5 process 'missing' is not available", buffer.getvalue())

    def test_run_mg5_processes_returns_cross_section_dictionary(self):
        calls = []

        def fake_get_xsec(process, run_tag, lambdas, k1, k2, k3, m2, w2, m3, w3, ecm):
            calls.append((process, run_tag, ecm))
            return {"hh": 1.25, "hhh": 0.031}[process]

        buffer = StringIO()
        with redirect_stdout(buffer):
            xsecs = run_mg5_processes(
                ["hh", "missing", "hhh"],
                "SCAN13.6",
                [1, 2, 3],
                0.9,
                0.1,
                0.0,
                300.0,
                1.0,
                500.0,
                2.0,
                13.6,
                get_xsec=fake_get_xsec,
                proc_location={"hh": "gg_hh_twoscalar/", "hhh": "gg_hhh_twoscalar/"},
            )

        self.assertEqual(xsecs, {"hh": 1.25, "hhh": 0.031})
        self.assertEqual(calls, [("hh", "SCAN13.6", 13.6), ("hhh", "SCAN13.6", 13.6)])
        self.assertIn("MG5 hh xsec [pb] = 1.25", buffer.getvalue())

    def test_run_mg5_processes_warns_when_selection_is_empty(self):
        buffer = StringIO()

        with redirect_stdout(buffer):
            xsecs = run_mg5_processes(
                [],
                "SCAN13.6",
                [1, 2, 3],
                0.9,
                0.1,
                0.0,
                300.0,
                1.0,
                500.0,
                2.0,
                13.6,
                get_xsec=lambda *args, **kwargs: self.fail("MG5 should not run"),
            )

        self.assertEqual(xsecs, {})
        self.assertIn("MG5ProcessesToRun is empty", buffer.getvalue())


if __name__ == "__main__":
    unittest.main()
