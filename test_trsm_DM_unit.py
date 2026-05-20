import math
import unittest
from contextlib import redirect_stdout
from io import StringIO

from test_trsm_DM import parse_micromegas_output, print_dm_info, test_dm


MICROMEGAS_OUTPUT = """
Dark matter candidate is '~X' with spin=0/2 mass=1.00E+03

Masses of odd sector Particles:
~X       : MX      = 1000.000 ||

==== Calculation of relic density =====
Xf=2.61e+01 Omega=4.90e-02

~X[~X]-nucleon cross sections[pb]:
 proton  SI 1.194E-09 [1.194E-09] SD 0.000E+00 [0.000E+00]
 neutron SI 1.218E-09 [1.218E-09] SD 0.000E+00 [0.000E+00]
"""


class TestTrsmDM(unittest.TestCase):
    def test_parse_micromegas_output(self):
        result = parse_micromegas_output(MICROMEGAS_OUTPUT)

        self.assertTrue(math.isclose(result.mdm, 1000.0))
        self.assertTrue(math.isclose(result.omega, 0.049))
        self.assertTrue(math.isclose(result.dir_det, 1.218e-9))

    def test_test_dm_returns_boolean_and_debug_info(self):
        passed, info, dm_exclusion_info = test_dm(
            lX=0.2,
            lPhiX=0.1,
            lSX=0.5,
            M3=1000.0,
            vs=500.0,
            a12=math.asin(0.3),
            M2=500.0,
            raw_output=MICROMEGAS_OUTPUT,
        )

        self.assertIs(passed, True)
        self.assertIn("DM check: Pass", info)
        self.assertIn("Omega=0.049", info)
        self.assertIn("DirDetLimit=2.05907e-09", info)
        self.assertTrue(math.isclose(dm_exclusion_info["dm_mdm"], 1000.0))
        self.assertTrue(math.isclose(dm_exclusion_info["dm_omega"], 0.049))
        self.assertTrue(math.isclose(dm_exclusion_info["dm_relic_upper_limit"], 0.1224))
        self.assertTrue(math.isclose(dm_exclusion_info["dm_dir_det"], 1.218e-9))
        self.assertTrue(math.isclose(dm_exclusion_info["dm_dir_det_limit"], 2.05907e-9, rel_tol=1e-5))
        self.assertFalse(dm_exclusion_info["dm_relic_excluded"])
        self.assertFalse(dm_exclusion_info["dm_direct_detection_excluded"])

    def test_print_dm_info_formats_debug_info(self):
        buffer = StringIO()

        with redirect_stdout(buffer):
            print_dm_info("DM check: Pass\n  LX=0.2 LHX=0.1\n  Reason: all DM checks passed")

        output = buffer.getvalue()
        self.assertIn("Dark Matter Check", output)
        self.assertIn("DM check: Pass", output)
        self.assertIn("LX=0.2 LHX=0.1", output)
        self.assertIn("+", output)


if __name__ == "__main__":
    unittest.main()
