import importlib.util
import math
from pathlib import Path
import unittest


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "run_single_point.py"


def load_module():
    spec = importlib.util.spec_from_file_location("run_single_point", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


MICROMEGAS_OUTPUT = """
Dark matter candidate is '~X' with spin=0/2 mass=1.00E+03

Masses of odd sector Particles:
~X       : MX      = 1000.000 ||

==== Calculation of relic density =====
Xf=2.61e+01 Omega=4.90e-02

==== Calculation of CDM-nucleons amplitudes  =====
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


class RunSinglePointTests(unittest.TestCase):
    def test_parses_micromegas_output_values(self):
        single = load_module()

        result = single.parse_micromegas_output(MICROMEGAS_OUTPUT)

        self.assertTrue(math.isclose(result.mdm, 1000.0))
        self.assertTrue(math.isclose(result.omega, 0.049))
        self.assertTrue(math.isclose(result.dir_det, 1.218e-9))
        self.assertEqual(result.indirect_line_channels, ())

    def test_parses_indirect_detection_channels(self):
        single = load_module()

        result = single.parse_micromegas_output(MICROMEGAS_OUTPUT_WITH_INDIRECT)

        self.assertEqual(len(result.indirect_line_channels), 1)
        channel = result.indirect_line_channels[0]
        self.assertTrue(math.isclose(channel.energy_gev, 50.0))
        self.assertTrue(math.isclose(channel.flux_cm2_s, 5.148348))

    def test_reads_single_card(self):
        single = load_module()
        import tempfile

        with tempfile.TemporaryDirectory() as tmpdir:
            card = Path(tmpdir) / "MO_inp100.dat"
            card.write_text(
                "\n".join(
                    [
                        "LX\t\t0.2",
                        "LHX\t\t0.1",
                        "LSX\t\t0.5",
                        "MX\t\t1000",
                        "vevs\t500",
                        "SinT\t0.3",
                        "Mh2\t\t500",
                    ]
                )
                + "\n",
                encoding="ascii",
            )

            point = single.read_card(card, index=100)

        self.assertEqual(point.index, 100)
        self.assertTrue(math.isclose(point.lx, 0.2))
        self.assertTrue(math.isclose(point.lsx, 0.5))
        self.assertTrue(math.isclose(point.mx, 1000.0))

    def test_summary_matches_cpp_rescaled_direct_detection_limit(self):
        single = load_module()
        point = single.PointInput(
            index=100,
            lx=0.2,
            lhx=0.1,
            lsx=0.5,
            mx=1000.0,
            vevs=500.0,
            sint=0.3,
            mh2=500.0,
        )
        result = single.MicromegasResult(mdm=1000.0, omega=0.049, dir_det=1.218e-9)

        summary = single.build_summary(point, result)

        self.assertTrue(math.isclose(summary.lux_base_limit, 8.243e-10, rel_tol=1e-5))
        self.assertTrue(math.isclose(summary.dir_det_limit, 2.05907e-9, rel_tol=1e-5))
        self.assertFalse(summary.relic_excluded)
        self.assertFalse(summary.direct_detection_excluded)
        self.assertFalse(summary.dm_excluded)

    def test_summary_includes_indirect_detection_exclusion(self):
        single = load_module()
        point = single.PointInput(
            index=100,
            lx=0.2,
            lhx=0.1,
            lsx=0.5,
            mx=1000.0,
            vevs=500.0,
            sint=0.3,
            mh2=500.0,
        )
        result = single.parse_micromegas_output(MICROMEGAS_OUTPUT_WITH_INDIRECT)

        summary = single.build_summary(point, result)

        self.assertFalse(summary.relic_excluded)
        self.assertFalse(summary.direct_detection_excluded)
        self.assertTrue(summary.indirect_limit.excluded)
        self.assertTrue(summary.dm_excluded)
        self.assertTrue(summary.indirect_limit.max_ratio > 1.0)
        self.assertIn("5.14835", single.summary_line(summary))

    def test_nonpositive_dark_matter_mass_is_rejected(self):
        single = load_module()

        with self.assertRaisesRegex(ValueError, "positive"):
            single.direct_detection_base_limit(0.0)


if __name__ == "__main__":
    unittest.main()
