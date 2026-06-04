import math
import unittest
from contextlib import redirect_stdout
from io import StringIO

from test_trsm_DM import (
    fermi_lat_r16_line_limit,
    parse_micromegas_output,
    print_dm_info,
    test_dm,
)


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


MICROMEGAS_OUTPUT_WITH_INDIRECT = (
    MICROMEGAS_OUTPUT
    + """
Fermi-LAT gamma-line flux estimate using R16 ROI approximation
FermiLAT_line_channel A A: E_gamma=5.000000E+01[GeV], sigmaV=2.098270E-18[cm^3 s^-1], N_gamma*sigmaV=4.196541E-18[cm^3 s^-1], Phi_R16=5.148348E+00[cm^-2 s^-1]
"""
)


class TestTrsmDM(unittest.TestCase):
    def test_parse_micromegas_output(self):
        result = parse_micromegas_output(MICROMEGAS_OUTPUT)

        self.assertTrue(math.isclose(result.mdm, 1000.0))
        self.assertTrue(math.isclose(result.omega, 0.049))
        self.assertTrue(math.isclose(result.dir_det, 1.218e-9))
        self.assertEqual(result.indirect_line_channels, ())

    def test_parse_micromegas_output_extracts_indirect_line_channels(self):
        result = parse_micromegas_output(MICROMEGAS_OUTPUT_WITH_INDIRECT)

        self.assertEqual(len(result.indirect_line_channels), 1)
        channel = result.indirect_line_channels[0]
        self.assertTrue(math.isclose(channel.energy_gev, 50.0))
        self.assertTrue(math.isclose(channel.flux_cm2_s, 5.148348))

    def test_fermi_lat_line_limit_interpolates_known_table(self):
        self.assertTrue(
            math.isclose(
                fermi_lat_r16_line_limit(49.1),
                6.40e-10,
                rel_tol=1e-12,
            )
        )

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
        self.assertFalse(dm_exclusion_info["dm_indirect_detection_excluded"])

    def test_indirect_detection_failure_fails_dm_check(self):
        passed, info, dm_exclusion_info = test_dm(
            lX=0.2,
            lPhiX=0.1,
            lSX=0.5,
            M3=1000.0,
            vs=500.0,
            a12=math.asin(0.3),
            M2=500.0,
            raw_output=MICROMEGAS_OUTPUT_WITH_INDIRECT,
        )

        self.assertIs(passed, False)
        self.assertIn("DM check: Fail", info)
        self.assertIn("Fermi-LAT gamma-line flux above limit", info)
        self.assertTrue(dm_exclusion_info["dm_indirect_available"])
        self.assertTrue(dm_exclusion_info["dm_indirect_detection_excluded"])
        self.assertTrue(dm_exclusion_info["dm_indirect_ratio"] > 1.0)
        self.assertTrue(math.isclose(dm_exclusion_info["dm_indirect_energy"], 50.0))

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
