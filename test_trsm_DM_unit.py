import math
import unittest
from contextlib import redirect_stdout
from io import StringIO

from test_trsm_DM import (
    INDIRECT_LINE_PRESCRIPTION,
    IndirectLineChannel,
    assess_indirect_limit,
    fermi_lat_energy_resolution,
    fermi_lat_r16_line_limit,
    indirect_detection_configuration,
    line_template_overlap,
    parse_micromegas_output,
    print_dm_info,
    test_dm,
)


MICROMEGAS_OUTPUT = """
Dark matter candidate is '~X' with spin=0/2 mass=1.00E+03

Masses of odd sector Particles:
~X       : MX      = 1000.000 ||

==== Calculation of relic density =====
darkOmega_error=0
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


MICROMEGAS_OUTPUT_WITH_RESCALED_INDIRECT = """
Dark matter candidate is '~X' with spin=0/2 mass=1.00E+03

Masses of odd sector Particles:
~X       : MX      = 1000.000 ||

==== Calculation of relic density =====
darkOmega_error=0
Xf=2.61e+01 Omega=4.90e-02

~X[~X]-nucleon cross sections[pb]:
 proton  SI 1.000E-12 [1.000E-12] SD 0.000E+00 [0.000E+00]
neutron SI 1.000E-12 [1.000E-12] SD 0.000E+00 [0.000E+00]

Fermi-LAT gamma-line flux estimate using R16 ROI approximation
FermiLAT_line_channel A A: E_gamma=5.000000E+01[GeV], sigmaV=1.000000E-27[cm^3 s^-1], N_gamma*sigmaV=1.000000E-27[cm^3 s^-1], Phi_R16=1.000000E-09[cm^-2 s^-1]
"""


MICROMEGAS_OUTPUT_WITH_SIGNED_AMPLITUDES = """
Dark matter candidate is '~X' with spin=0/2 mass=7.303E+02
~X       : MX      = 730.3 ||
darkOmega_error=0
Xf=2.61e+01 Omega=1.79e-04

==== Calculation of CDM-nucleons amplitudes  =====
~X[~X]-nucleon micrOMEGAs amplitudes
proton SI -1.965E-09 [-1.965E-09] SD 0.000E+00 [0.000E+00]
neutron SI -1.984E-09 [-1.984E-09] SD 0.000E+00 [0.000E+00]
~X[~X]-nucleon cross sections[pb]:
 proton  SI 1.683E-09 [1.683E-09] SD 0.000E+00 [0.000E+00]
 neutron SI 1.717E-09 [1.717E-09] SD 0.000E+00 [0.000E+00]
"""


MICROMEGAS_OUTPUT_WITH_AMPLITUDE_ONLY = """
~X       : MX      = 730.3 ||
darkOmega_error=0
Xf=2.61e+01 Omega=1.79e-04
~X[~X]-nucleon micrOMEGAs amplitudes
neutron SI -1.984E-09 [-1.984E-09] SD 0.000E+00 [0.000E+00]
"""


class TestTrsmDM(unittest.TestCase):
    def test_parse_micromegas_output(self):
        result = parse_micromegas_output(MICROMEGAS_OUTPUT)

        self.assertTrue(math.isclose(result.mdm, 1000.0))
        self.assertTrue(math.isclose(result.omega, 0.049))
        self.assertTrue(math.isclose(result.xf, 26.1))
        self.assertTrue(math.isclose(result.dir_det, 1.218e-9))
        self.assertEqual(result.indirect_line_channels, ())

    def test_parse_micromegas_output_extracts_indirect_line_channels(self):
        result = parse_micromegas_output(MICROMEGAS_OUTPUT_WITH_INDIRECT)

        self.assertEqual(len(result.indirect_line_channels), 1)
        channel = result.indirect_line_channels[0]
        self.assertTrue(math.isclose(channel.energy_gev, 50.0))
        self.assertTrue(math.isclose(channel.flux_cm2_s, 5.148348))

    def test_parse_uses_cross_section_not_earlier_signed_amplitude(self):
        result = parse_micromegas_output(MICROMEGAS_OUTPUT_WITH_SIGNED_AMPLITUDES)

        self.assertTrue(math.isclose(result.dir_det, 1.717e-9))
        self.assertGreater(result.dir_det, 0.0)

    def test_parse_rejects_amplitude_without_cross_section_block(self):
        with self.assertRaisesRegex(ValueError, r"cross-sections\[pb\] section"):
            parse_micromegas_output(MICROMEGAS_OUTPUT_WITH_AMPLITUDE_ONLY)

    def test_fermi_lat_line_limit_interpolates_known_table(self):
        self.assertTrue(
            math.isclose(
                fermi_lat_r16_line_limit(49.1),
                6.40e-10,
                rel_tol=1e-12,
            )
        )

    def test_test_dm_uses_updated_direct_detection_limit(self):
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

        self.assertIs(passed, False)
        self.assertIn("DM check: Fail", info)
        self.assertIn("Omega=", info)
        self.assertIn("DirDetLimit=", info)
        self.assertIn("DirDet above rescaled direct-detection limit", info)
        self.assertTrue(math.isclose(dm_exclusion_info["dm_mdm"], 1000.0))
        self.assertTrue(math.isclose(dm_exclusion_info["dm_omega"], 0.049))
        self.assertTrue(math.isclose(dm_exclusion_info["dm_xf"], 26.1))
        self.assertTrue(math.isclose(dm_exclusion_info["dm_freezeout_temperature_GeV"], 1000.0 / 26.1))
        self.assertTrue(math.isclose(dm_exclusion_info["dm_relic_upper_limit"], 0.121))
        self.assertTrue(math.isclose(dm_exclusion_info["dm_dir_det"], 1.218e-9))
        self.assertTrue(math.isclose(dm_exclusion_info["dm_dir_det_limit"], __import__("test_trsm_DM").direct_detection_base_limit(1000)*.12/.049, rel_tol=1e-5))
        self.assertTrue(math.isclose(dm_exclusion_info["dm_lux_base_limit"], __import__("test_trsm_DM").direct_detection_base_limit(1000), rel_tol=1e-5))
        self.assertFalse(dm_exclusion_info["dm_relic_excluded"])
        self.assertTrue(dm_exclusion_info["dm_direct_detection_excluded"])
        self.assertFalse(dm_exclusion_info["dm_indirect_detection_excluded"])
        self.assertEqual(dm_exclusion_info["dm_limit_model"], "lz-ws2024-observed-v2")
        self.assertEqual(dm_exclusion_info["dm_indirect_channels_seen"], 0)
        self.assertEqual(dm_exclusion_info["dm_indirect_channels_used"], 0)

    def test_indirect_detection_limit_is_relic_rescaled(self):
        passed, info, dm_exclusion_info = test_dm(
            lX=0.2,
            lPhiX=0.1,
            lSX=0.5,
            M3=1000.0,
            vs=500.0,
            a12=math.asin(0.3),
            M2=500.0,
            raw_output=MICROMEGAS_OUTPUT_WITH_RESCALED_INDIRECT,
        )

        self.assertIs(passed, True)
        self.assertIn("DM check: Pass", info)
        self.assertFalse(dm_exclusion_info["dm_direct_detection_excluded"])
        self.assertFalse(dm_exclusion_info["dm_indirect_detection_excluded"])
        expected_limit = fermi_lat_r16_line_limit(50.0) * (0.12 / 0.049) ** 2
        self.assertTrue(
            math.isclose(
                dm_exclusion_info["dm_indirect_limit"],
                expected_limit,
                rel_tol=1e-5,
            )
        )
        self.assertLess(dm_exclusion_info["dm_indirect_ratio"], 1.0)

    def test_loop_derived_z_gamma_line_is_combined_with_gamma_gamma(self):
        # DM/main.c derives XX -> Z A from XX -> A A (trsm_loop_zgamma_ratio)
        # and prints it in the standard line format at its own energy.
        output = MICROMEGAS_OUTPUT_WITH_RESCALED_INDIRECT.replace(
            "FermiLAT_line_channel A A: E_gamma=5.000000E+01[GeV], sigmaV=1.000000E-27[cm^3 s^-1], N_gamma*sigmaV=1.000000E-27[cm^3 s^-1], Phi_R16=1.000000E-09[cm^-2 s^-1]",
            "FermiLAT_line_channel A A: E_gamma=3.000000E+02[GeV], sigmaV=1.0E-30[cm^3 s^-1], N_gamma*sigmaV=2.0E-30[cm^3 s^-1], Phi_R16=1.000000E-10[cm^-2 s^-1]\n"
            'TRSM_zgamma_line_v1 {"ratio_ZA_over_AA":40,"E_gamma_GeV":293.071,"sigma_v_cm3_s":4.0e-29}\n'
            "FermiLAT_line_channel Z A: E_gamma=2.930710E+02[GeV], sigmaV=4.0E-29[cm^3 s^-1], N_gamma*sigmaV=4.0E-29[cm^3 s^-1], Phi_R16=2.000000E-09[cm^-2 s^-1]",
        )
        result = parse_micromegas_output(output)
        self.assertEqual(
            [channel.energy_gev for channel in result.indirect_line_channels],
            [300.0, 293.071],
        )

        passed, _, dm_exclusion_info = test_dm(
            lX=0.2, lPhiX=0.1, lSX=0.5, M3=1000.0, vs=500.0,
            a12=math.asin(0.3), M2=500.0, raw_output=output,
        )
        # Alone, the A A line has ratio ~0.12; the Z A line (~2.6) excludes.
        # The lines are 0.26 sigma_E apart, so the A A flux adds to the Z A test
        # with weight exp(-(dE)^2/(4 sigma_E^2)) ~ 0.98.
        self.assertIs(passed, False)
        self.assertTrue(dm_exclusion_info["dm_indirect_detection_excluded"])
        self.assertEqual(dm_exclusion_info["dm_indirect_channels_used"], 2)
        self.assertTrue(math.isclose(dm_exclusion_info["dm_indirect_energy"], 293.071))
        sigma = fermi_lat_energy_resolution(293.071) * 293.071
        weight = math.exp(-((300.0 - 293.071) ** 2) / (4 * sigma**2))
        self.assertTrue(0.98 < weight < 0.99)
        expected = (2.0e-9 + weight * 1.0e-10) / (fermi_lat_r16_line_limit(293.071) * (0.12 / 0.049) ** 2)
        self.assertTrue(math.isclose(dm_exclusion_info["dm_indirect_flux"], 2.0e-9 + weight * 1.0e-10, rel_tol=1e-12))
        self.assertTrue(math.isclose(dm_exclusion_info["dm_indirect_ratio"], expected, rel_tol=1e-9))

    def test_negative_relic_density_fails_dm_check(self):
        invalid_output = MICROMEGAS_OUTPUT.replace(
            "Omega=4.90e-02",
            "Omega=-1.00e+00",
        )

        passed, info, dm_exclusion_info = test_dm(
            lX=0.2,
            lPhiX=0.1,
            lSX=0.5,
            M3=1000.0,
            vs=500.0,
            a12=math.asin(0.3),
            M2=500.0,
            raw_output=invalid_output,
        )

        self.assertIsNone(passed)
        self.assertIn("Relic density Omega must be finite and non-negative", info)
        self.assertIsNone(dm_exclusion_info["dm_omega"])
        self.assertIsNone(dm_exclusion_info["dm_freezeout_temperature_GeV"])

    def test_missing_or_invalid_xf_does_not_invent_a_freezeout_temperature(self):
        for raw_output in (
            MICROMEGAS_OUTPUT.replace("Xf=2.61e+01 ", ""),
            MICROMEGAS_OUTPUT.replace("Xf=2.61e+01", "Xf=0"),
        ):
            with self.subTest(raw_output=raw_output):
                passed, _, diagnostics = test_dm(
                    .2, .1, .5, 1000, 500, .3, 500, raw_output=raw_output,
                )
                self.assertFalse(passed)  # This point fails direct detection.
                self.assertIsNone(diagnostics["dm_xf"])
                self.assertIsNone(diagnostics["dm_freezeout_temperature_GeV"])

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
        self.assertEqual(dm_exclusion_info["dm_indirect_channels_seen"], 1)
        self.assertEqual(dm_exclusion_info["dm_indirect_channels_used"], 1)

    def test_print_dm_info_formats_debug_info(self):
        buffer = StringIO()

        with redirect_stdout(buffer):
            print_dm_info("DM check: Pass\n  LX=0.2 LHX=0.1\n  Reason: all DM checks passed")

        output = buffer.getvalue()
        self.assertIn("Dark Matter Check", output)
        self.assertIn("DM check: Pass", output)
        self.assertIn("LX=0.2 LHX=0.1", output)
        self.assertIn("+", output)


class TestGammaLineCombination(unittest.TestCase):
    """Lines within the Fermi-LAT energy resolution are one line to the LAT."""

    def ratio(self, channels, combine=True):
        return assess_indirect_limit(channels, omega=0.12, combine_lines=combine)

    def test_energy_resolution_table_and_interpolation(self):
        self.assertTrue(math.isclose(fermi_lat_energy_resolution(9.7167421875), 0.0632673773133))
        between = fermi_lat_energy_resolution(math.sqrt(170.734375 * 302.88728125))
        self.assertTrue(math.isclose(between, (0.0804665341422 + 0.0888698056022) / 2))
        for energy in (0.214, 1.0, 62.5, 462.0):
            self.assertTrue(0.06 < fermi_lat_energy_resolution(energy) < 0.16)
        with self.assertRaises(ValueError):
            fermi_lat_energy_resolution(2000.0)

    def test_overlap_weight(self):
        self.assertEqual(line_template_overlap(100.0, 100.0), 1.0)
        sigma = fermi_lat_energy_resolution(100.0) * 100.0
        self.assertTrue(math.isclose(line_template_overlap(100.0, 100.0 + sigma), math.exp(-0.25)))
        self.assertLess(line_template_overlap(100.0, 100.0 + 8 * sigma), 1e-6)

    def test_isolated_line_is_unchanged(self):
        line = (IndirectLineChannel(62.5, 1.0e-10),)
        self.assertEqual(self.ratio(line), self.ratio(line, combine=False))
        expected = 1.0e-10 / fermi_lat_r16_line_limit(62.5)
        self.assertTrue(math.isclose(self.ratio(line).max_ratio, expected, rel_tol=1e-15))

    def test_coincident_lines_add(self):
        # Each line alone is at 60% of the limit; together they exclude.
        flux = 0.6 * fermi_lat_r16_line_limit(200.0)
        lines = (IndirectLineChannel(200.0, flux), IndirectLineChannel(200.0, flux))
        self.assertFalse(self.ratio(lines, combine=False).excluded)
        combined = self.ratio(lines)
        self.assertTrue(combined.excluded)
        self.assertTrue(math.isclose(combined.max_ratio, 1.2, rel_tol=1e-12))

    def test_lines_one_sigma_apart(self):
        sigma = fermi_lat_energy_resolution(200.0) * 200.0
        lines = (IndirectLineChannel(200.0, 3.0e-11), IndirectLineChannel(200.0 - sigma, 1.0e-11))
        result = self.ratio(lines)
        # Each line energy is tested with the other line's flux added. Near
        # 200 GeV the weight is exp(-1/4) (sigma_E is evaluated at the test energy).
        low = 200.0 - sigma
        at_high = (3.0e-11 + math.exp(-0.25) * 1.0e-11) / fermi_lat_r16_line_limit(200.0)
        at_low = (1.0e-11 + line_template_overlap(low, 200.0) * 3.0e-11) / fermi_lat_r16_line_limit(low)
        self.assertTrue(math.isclose(result.max_ratio, max(at_high, at_low), rel_tol=1e-12))
        self.assertEqual(result.channels_used, 2)

    def test_well_separated_lines_do_not_combine(self):
        # gamma-gamma and Z-gamma lines at M3 = 62.5 GeV are about 8 sigma_E apart.
        lines = (IndirectLineChannel(62.5, 1.0e-10), IndirectLineChannel(29.24, 3.0e-11))
        separate, combined = self.ratio(lines, combine=False), self.ratio(lines)
        self.assertEqual(combined.energy_gev, separate.energy_gev)
        self.assertTrue(math.isclose(combined.max_ratio, separate.max_ratio, rel_tol=1e-6))

    def test_combination_respects_abundance_rescaling(self):
        flux = 0.6 * fermi_lat_r16_line_limit(200.0)
        lines = (IndirectLineChannel(200.0, flux), IndirectLineChannel(200.0, flux))
        half = assess_indirect_limit(lines, omega=0.06)
        self.assertTrue(math.isclose(half.max_ratio, 1.2 * 0.25, rel_tol=1e-12))
        self.assertFalse(assess_indirect_limit(lines, omega=0.0).available)

    def test_configuration_records_prescription(self):
        configuration = indirect_detection_configuration()
        self.assertEqual(configuration["line_combination"], INDIRECT_LINE_PRESCRIPTION)
        self.assertEqual(configuration["energy_range_GeV"], [0.214, 462.0])


if __name__ == "__main__":
    unittest.main()
