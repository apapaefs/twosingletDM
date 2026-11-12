"""Tests for trsm_sidm.py (H3 H3 -> H3 H3 self-scattering).

Run with ``python -m unittest test_trsm_sidm.py``.  The CalcHEP comparison
runs only when TRSM_SIDM_MAIN points to a compiled ``TRSM/main_sidm``; two
solver cross-checks run only when SciPy is installed.
"""

import csv
import gzip
import json
import math
import os
import tempfile
import unittest
from pathlib import Path

import numpy as np

import trsm_sidm as S

try:
    from scipy import integrate as _integrate, special as _special
except ImportError:  # the module itself needs only NumPy
    _integrate = _special = None


def contact_point(lX=0.5, M3=50.0):
    return S.SIDMPoint(M3=M3, M2=300.0, K133=0.0, K233=0.0, lX=lX)


def dark_h2_point(M3=50.0, dM=1.25e-5, K233=1.0, lX=0.1, gamma2_other=0.0):
    return S.SIDMPoint(M3=M3, M2=2.0 * M3 + dM, K133=0.0, K233=K233, lX=lX,
                       gamma1_other=4.07e-3, gamma2_other=gamma2_other)


def yukawa_point(b, M3=60.0, M2=20.0, lX=0.0):
    """Single light mediator with Born parameter b = alpha M3/M2."""
    K233 = math.sqrt(4.0 * math.pi * b * M2 * M3)
    return S.SIDMPoint(M3=M3, M2=M2, K133=0.0, K233=K233, lX=lX)


class ConventionTests(unittest.TestCase):
    def test_portal_couplings_match_rotation_matrix(self):
        lPhiX, lSX, vs, a12, v = 0.13, -0.42, 377.0, -0.31, 246.0
        R11, R12 = math.cos(a12), -math.sin(a12)  # generate_trsm_info.Rmatrix
        R21, R22 = math.sin(a12), math.cos(a12)
        K133, K233 = S.vxzero_portal_couplings(lPhiX, lSX, vs, a12, v)
        self.assertAlmostEqual(K133, 0.5 * (lPhiX * R11 * v + lSX * R12 * vs), places=12)
        self.assertAlmostEqual(K233, 0.5 * (lPhiX * R21 * v + lSX * R22 * vs), places=12)

    def test_identical_scalar_width_normalisation(self):
        # Feynman rule 2K: Gamma = (2K)^2 beta / (32 pi M) = K^2 beta / (8 pi M)
        width = S.scalar_to_identical_scalar_width(10.0, 100.0, 3.0)
        beta = math.sqrt(1.0 - 4.0 * 100.0 / 1.0e4)
        self.assertAlmostEqual(width, 4.0 * 9.0 * beta / (32.0 * math.pi * 100.0), places=15)
        self.assertEqual(S.scalar_to_identical_scalar_width(50.0, 100.0, 3.0), 0.0)

    def test_unit_conversion(self):
        self.assertAlmostEqual(float(S.per_mass_cm2_g(1.0, 1.0)), 0.3893793721e-27 / 1.78266192e-24, places=18)

    def test_born_parameter_and_regime(self):
        p = S.SIDMPoint(M3=10.0, M2=20.0, K133=1.0, K233=125.0, lX=0.0)
        self.assertAlmostEqual(p.alpha(2), 125.0**2 / (4.0 * math.pi * 100.0))
        self.assertAlmostEqual(p.born_parameter(2), p.alpha(2) * 10.0 / 20.0)
        self.assertEqual(p.regime(), "strong")
        self.assertEqual(yukawa_point(0.2).regime(), "nonperturbative")
        self.assertEqual(yukawa_point(1.0e-4).regime(), "perturbative")

    def test_s_minus_mass_squared_is_accurate_near_threshold(self):
        M3, eps = 50.0, 1.0e-12
        Mi = 2.0 * M3 * (1.0 + eps)
        vr = S.resonance_velocity(M3, Mi)
        v = np.array([0.5 * vr, vr, 2.0 * vr])
        exact = M3**2 * (v**2 - vr**2) / ((1 - v**2 / 4) * (1 - vr**2 / 4))
        np.testing.assert_allclose(S._s_minus_mass2(M3, v, Mi), exact, rtol=1e-12, atol=0.0)

    def test_mean_relative_velocity(self):
        # <v> = 4 sigma_1D / sqrt(pi) for f(v) ~ v^2 exp(-v^2/(4 sigma^2))
        sigma = S.mean_relative_velocity_to_dispersion(1900.0)
        v = np.linspace(0.0, 30.0 * sigma, 200001)
        f = v**2 * np.exp(-(v**2) / (4.0 * sigma**2))
        self.assertAlmostEqual(np.sum(v * f) / np.sum(f), 1900.0, places=3)


class SelfEnergyTests(unittest.TestCase):
    def test_loop_function_is_continuous_at_threshold(self):
        self.assertAlmostEqual(float(S.loop_real_part(1.0e-9)), 0.0, places=15)
        self.assertAlmostEqual(S.loop_real_part_at_mass(50.0, 100.0), 0.0)
        self.assertAlmostEqual(S.loop_real_part_at_mass(50.0, 100.0 + 1e-9), 0.0, places=8)
        self.assertAlmostEqual(S.loop_real_part_at_mass(50.0, 100.0 - 1e-9), 0.0, places=4)

    def test_cusp_below_threshold(self):
        # Re Pi(s) = K^2/(8 pi^2) (-pi betabar + ...) just below threshold
        M3, Mi = 50.0, 99.999
        betabar = math.sqrt(4.0 * M3**2 / Mi**2 - 1.0)
        self.assertAlmostEqual(S.loop_real_part_at_mass(M3, Mi) / (-math.pi * betabar), 1.0, places=2)

    def test_dispersion_relation(self):
        # Re[B(s) - B(s0)] = (s - s0)/pi PV int ds' Im B(s') / ((s'-s)(s'-s0)),
        # Im B = pi beta(s').  Checked with s, s0 below threshold (no PV needed).
        m = 1.0

        def real_part(s):
            bb = math.sqrt(4.0 * m * m / s - 1.0)
            return -2.0 * bb * math.atan2(1.0, bb)

        s, s0 = 3.2, 1.1
        u = np.linspace(0.0, 1.0, 400001)[1:-1]       # s' = 4 m^2 / u^2
        sp = 4.0 * m * m / u**2
        integrand = np.sqrt(1.0 - 4.0 * m * m / sp) / ((sp - s) * (sp - s0)) * 8.0 * m * m / u**3
        integral = np.sum(integrand) * (u[1] - u[0])
        self.assertAlmostEqual((s - s0) * integral / (real_part(s) - real_part(s0)), 1.0, places=3)

    def test_subtraction_vanishes_on_shell_above_threshold(self):
        p = dark_h2_point(dM=0.3, K233=20.0)
        vr, _ = p.resonance(2)
        self.assertAlmostEqual(float(S.self_energy_real(p, np.array([vr]), 2)[0]), 0.0, places=9)

    def test_pole_below_threshold_stays_at_the_physical_mass(self):
        # D(s) = s - M^2 + Pi(s) - Re Pi(M^2).  Near threshold, with
        # 4 kb^2 = 4 M3^2 - M^2 and g = K^2/(8 pi M3):
        # D = (kb + i k) [4 (kb - i k) + g], which vanishes at k = i kb.
        M3, dM, K = 10.0, -1.0e-6, 5.0
        p = S.SIDMPoint(M3=M3, M2=2.0 * M3 + dM, K133=0.0, K233=K, lX=0.0)
        kb = 0.5 * math.sqrt((2.0 * M3 - p.M2) * (2.0 * M3 + p.M2))
        g = K * K / (8.0 * math.pi * M3)
        v = np.array([1e-7, 1e-5, 1e-4, 1e-3])
        k = 0.5 * M3 * v
        sigma, _ = S.cross_sections(p, v, yukawa="off")
        flatte = 8.0 * math.pi / (4.0 * M3**2) * (K * K / (4.0 * math.pi)) ** 2 \
            / np.abs((kb + 1j * k) * (4.0 * (kb - 1j * k) + g)) ** 2
        np.testing.assert_allclose(sigma, flatte, rtol=2e-3)

    def test_strongly_coupled_shallow_pole_gives_universal_cross_section(self):
        # g >> kb: sigma -> 8 pi / (kb^2 + k^2), independent of the coupling.
        M3, dM = 10.0, -1.0e-10
        for K in (20.0, 40.0):
            p = S.SIDMPoint(M3=M3, M2=2.0 * M3 + dM, K133=0.0, K233=K, lX=0.0)
            kb = 0.5 * math.sqrt((2.0 * M3 - p.M2) * (2.0 * M3 + p.M2))
            v = np.array([1e-6, 1e-5])
            k = 0.5 * M3 * v
            sigma, _ = S.cross_sections(p, v, yukawa="off")
            np.testing.assert_allclose(sigma, 8.0 * math.pi / (kb**2 + k**2), rtol=2e-2)


class AmplitudeTests(unittest.TestCase):
    def test_contact_term_cross_section(self):
        p = contact_point(lX=0.5)
        v = np.array([1e-4, 1e-2, 0.3])
        s, _, _ = S.kinematics(p.M3, v)
        expected = 9.0 * 0.5**2 / (8.0 * math.pi * s)  # 36 lX^2/(32 pi s)
        for mode in ("tree-fixed-width", "tree-running-width"):
            sigma, sigma_v = S.cross_sections(p, v, mode=mode)
            np.testing.assert_allclose(sigma, expected, rtol=1e-13)
            np.testing.assert_allclose(sigma_v, sigma, rtol=1e-13)  # isotropic

    def test_exchange_integral_with_narrow_forward_peak(self):
        # int_-1^1 [1/(t-m^2) + 1/(u-m^2)]^2 dx in closed form, for k up to 1000 m.
        p = S.SIDMPoint(M3=1000.0, M2=0.05, K133=0.0, K233=1.0, lX=0.0)
        for v in (1e-4, 1e-3, 1e-2, 0.1):
            s, k2, _ = S.kinematics(p.M3, np.array([v]))
            A, B = p.M2**2 + 2.0 * k2[0], 2.0 * k2[0]
            integral = 4.0 * A * A * (1.0 / (A * A * (A * A - B * B)) + math.log((A + B) / (A - B)) / (2.0 * A**3 * B))
            exact = 16.0 * p.K233**4 * integral / (64.0 * math.pi * s[0])
            numeric, _ = S.cross_sections(p, [v], mode="tree-fixed-width")
            self.assertAlmostEqual(numeric[0] / exact, 1.0, places=6)

    def test_amplitude_symmetric_under_t_u_exchange(self):
        p = S.SIDMPoint(M3=300.0, M2=20.0, K133=-1.2, K233=60.0, lX=0.1, gamma2_other=1e-7)
        x = np.linspace(-0.9, 0.9, 7)
        for mode in S.AMPLITUDE_MODES:
            amp = S.scattering_amplitude(p, [0.3], x, mode=mode)[0]
            np.testing.assert_allclose(amp, amp[::-1], rtol=1e-12)

    def test_kmatrix_satisfies_elastic_unitarity(self):
        for yukawa in ("off", "on"):
            p = S.SIDMPoint(M3=40.0, M2=90.0, K133=30.0, K233=70.0, lX=-0.2)
            v = np.geomspace(1e-5, 1.0, 40)
            _, _, a = S.partial_waves(p, v, yukawa=yukawa)
            _, _, beta = S.kinematics(p.M3, v)
            lhs, rhs = a.imag, 0.5 * beta[:, None] * np.abs(a) ** 2
            np.testing.assert_allclose(lhs, rhs, rtol=1e-9, atol=1e-30)

    def test_kmatrix_with_inelastic_width_is_below_unitarity_circle(self):
        p = S.SIDMPoint(M3=40.0, M2=90.0, K133=30.0, K233=70.0, lX=-0.2,
                        gamma1_other=0.2, gamma2_other=0.8)
        v = np.geomspace(1e-5, 1.0, 40)
        _, _, a = S.partial_waves(p, v, yukawa="off")
        _, _, beta = S.kinematics(p.M3, v)
        self.assertTrue(np.all(a[:, 0].imag >= 0.5 * beta * np.abs(a[:, 0]) ** 2 * (1 - 1e-12)))

    def test_resonance_saturates_swave_unitarity(self):
        p = dark_h2_point()
        vr, _ = p.resonance(2)
        sigma, sigma_v = S.cross_sections(p, [vr])
        _, k2, _ = S.kinematics(p.M3, np.array([vr]))
        self.assertAlmostEqual(sigma[0] * k2[0] / (8.0 * math.pi), 1.0, places=5)
        self.assertAlmostEqual(sigma_v[0] / sigma[0], 1.0, places=6)   # pure s-wave
        self.assertAlmostEqual(float(S.unitarity_limit_cm2_g(p.M3, vr * S.C_KMS))
                               / float(S.per_mass_cm2_g(sigma[0], p.M3)), 1.0, places=5)

    def test_unitarised_partial_waves_respect_bound_everywhere(self):
        p = dark_h2_point(M3=10.0, dM=1e-6, K233=125.0)
        v = np.geomspace(1e-6, 0.1, 60)
        _, _, a = S.partial_waves(p, v)
        _, _, beta = S.kinematics(p.M3, v)
        self.assertTrue(np.all(0.5 * beta[:, None] * np.abs(a) <= 1.0 + 1e-12))
        # whereas the CalcHEP-like fixed-width tree amplitude violates it
        _, K, _ = S.partial_waves(p, v, mode="tree-fixed-width")
        self.assertGreater(np.max(0.5 * beta * np.abs(K[:, 0])), 10.0)

    def test_kmatrix_matches_dyson_propagator_when_s_channel_dominates(self):
        p = dark_h2_point(K233=1.0, lX=0.0)
        vr, dv = p.resonance(2)
        v = vr + dv * np.linspace(-5, 5, 11)
        a, _ = S.cross_sections(p, v, mode="kmatrix")
        b, _ = S.cross_sections(p, v, mode="tree-running-width")
        np.testing.assert_allclose(a, b, rtol=1e-4)

    def test_angular_integral_matches_partial_wave_sum(self):
        p = S.SIDMPoint(M3=200.0, M2=500.0, K133=8.7, K233=280.0, lX=0.5, gamma2_other=1.0)
        v = np.array([1e-3, 0.3, 1.0])
        sigma, _ = S.cross_sections(p, v)
        ls, _, a = S.partial_waves(p, v)
        s, _, _ = S.kinematics(p.M3, v)
        pw = 8.0 * math.pi / s * np.sum((2 * ls + 1) * np.abs(a) ** 2, axis=1)
        np.testing.assert_allclose(sigma, pw, rtol=1e-8)

    def test_viscosity_matches_phase_shift_formula(self):
        # Colquhoun et al. (2011.04679) Eq. 28 for identical bosons:
        # sigma_V^even = 8 pi/k^2 sum_n (2n+1)(2n+2)/(4n+3) sin^2(delta_{2n+2} - delta_{2n});
        # this module's sigma_V carries the extra factor 3/2.
        rng = np.random.default_rng(7)
        delta = np.concatenate([rng.uniform(-1.5, 1.5, 5), [0.0]])   # l = 0, 2, ..., 8, then 10
        beta, k = 0.01, 0.37
        s = 4.0 * k * k / beta**2
        ls = np.arange(0, 10, 2)
        a = 2.0 / beta * np.exp(1j * delta[:5]) * np.sin(delta[:5])
        x, w = np.polynomial.legendre.leggauss(64)
        amp = 16.0 * math.pi * (S._even_legendre(x, 8) * ((2 * ls + 1) * a)[None, :]).sum(axis=1)
        mine = 3.0 / (128.0 * math.pi * s) * np.sum(w * (1.0 - x * x) * np.abs(amp) ** 2)
        n = np.arange(5)
        colquhoun = 8.0 * math.pi / k**2 * np.sum((2 * n + 1) * (2 * n + 2) / (4 * n + 3)
                                                 * np.sin(delta[1:] - delta[:-1]) ** 2)
        self.assertAlmostEqual(mine / (1.5 * colquhoun), 1.0, places=10)

    def test_forward_peaked_scattering_has_smaller_viscosity(self):
        # light mediator at high velocity: t/u-channel peaks at x = +-1
        p = S.SIDMPoint(M3=1000.0, M2=4.0, K133=0.0, K233=5.0, lX=0.0)
        sigma, sigma_v = S.cross_sections(p, [0.05])
        self.assertLess(sigma_v[0], 0.9 * sigma[0])


class SchroedingerTests(unittest.TestCase):
    def test_riccati_bessel_wronskian_and_closed_forms(self):
        z = np.concatenate([np.geomspace(1e-9, 1.0, 30), np.linspace(1.0, 150.0, 200)])
        for l in (0, 2, 4, 6, 8):
            j, dj, n, dn = S.riccati_bessel(l, z)
            np.testing.assert_allclose(j * dn - dj * n, 1.0, rtol=1e-9)
        j, dj, n, dn = S.riccati_bessel(2, z)
        big = z > 0.5
        np.testing.assert_allclose(j[big], ((3 / z**2 - 1) * np.sin(z) - 3 / z * np.cos(z))[big], rtol=1e-9, atol=1e-13)
        np.testing.assert_allclose(n, -(3 / z**2 - 1) * np.cos(z) - 3 / z * np.sin(z), rtol=1e-12)

    @unittest.skipIf(_special is None, "SciPy not installed")
    def test_riccati_bessel_against_scipy(self):
        z = np.concatenate([np.geomspace(1e-9, 1.0, 30), np.linspace(1.0, 150.0, 300)])
        for l in (0, 2, 4, 6, 8):
            j, _, n, _ = S.riccati_bessel(l, z)
            np.testing.assert_allclose(j, z * _special.spherical_jn(l, z), rtol=1e-11)
            np.testing.assert_allclose(n, z * _special.spherical_yn(l, z), rtol=1e-11)

    def test_weak_coupling_scattering_length(self):
        # a = (b/m)(1 + b/2 + ...): first and second Born terms of the Yukawa potential
        M3, m = 50.0, 5.0
        for b in (1e-3, 1e-2):
            k = np.array([1e-6 * m])
            N, D = S.yukawa_tan_delta(M3, [b * m / M3], [m], k, [0])
            self.assertAlmostEqual((N / D)[0, 0] / k[0] / (b / m) / (1.0 + 0.5 * b), 1.0, places=4)

    def test_bound_state_thresholds(self):
        # A zero-energy s-wave state appears at alpha M3/m = 1.6798, a d-wave one at 21.895.
        M3, m = 50.0, 5.0
        k = np.array([1e-7 * m])

        def below_resonance(b, l):
            _, D = S.yukawa_tan_delta(M3, [b * m / M3], [m], k, [l])
            return D[0, 0]

        self.assertLess(below_resonance(1.6790, 0) * below_resonance(1.6806, 0), 0.0)
        self.assertGreater(below_resonance(1.60, 0) * below_resonance(1.6790, 0), 0.0)
        self.assertLess(below_resonance(21.88, 2) * below_resonance(21.91, 2), 0.0)

    def test_step_size_convergence(self):
        M3, m = 50.0, 5.0
        k = np.array([1e-7, 1e-3, 0.05, 0.3, 1.0, 5.0, 20.0]) * m
        for b in (0.3, 5.0, 60.0):
            for l in (0, 4, 8):
                N, D = S.yukawa_tan_delta(M3, [b * m / M3], [m], k, [l])
                N2, D2 = S.yukawa_tan_delta(M3, [b * m / M3], [m], k, [l], eta=S._SOLVER_ETA / 4)
                error = np.abs(np.sin(np.arctan2(N, D) - np.arctan2(N2, D2)))
                self.assertLess(error.max(), 2e-7)

    @unittest.skipIf(_integrate is None, "SciPy not installed")
    def test_against_direct_integration(self):
        M3, alphas, masses = 60.0, [0.3, 0.05], [20.0, 125.09]

        def reference(k, l):
            r0, rm = 1e-9, 40.0 / min(masses)

            def rhs(r, y):
                U = M3 * sum(a * math.exp(-m * r) for a, m in zip(alphas, masses)) / r
                return [y[1], (l * (l + 1) / r**2 - k * k - U) * y[0]]

            sol = _integrate.solve_ivp(rhs, [r0, rm], [1.0, (l + 1) / r0], method="DOP853",
                                       rtol=1e-12, atol=1e-300)
            L, z = rm * sol.y[1, -1] / sol.y[0, -1], k * rm
            j, dj, n, dn = S.riccati_bessel(l, np.array([z]))
            return math.atan2((z * dj - L * j)[0], (z * dn - L * n)[0])

        for k in (0.01, 3.0, 40.0):
            N, D = S.yukawa_tan_delta(M3, alphas, masses, np.array([k]), [0, 2])
            for column, l in enumerate((0, 2)):
                self.assertLess(abs(math.sin(math.atan2(N[0, column], D[0, column]) - reference(k, l))), 1e-6)

    def test_solver_agrees_with_born_at_weak_coupling(self):
        v = np.array([30.0, 1000.0, 30000.0]) / S.C_KMS
        p = yukawa_point(1.0e-3)
        born, born_v = S.cross_sections(p, v, yukawa="off")
        exact, exact_v = S.cross_sections(p, v, yukawa="on")
        # the exchange dominates here, so sigma changes by about b = 1e-3
        np.testing.assert_allclose(exact, born, rtol=4e-3)
        np.testing.assert_allclose(exact_v, born_v, rtol=4e-3)

    def test_zero_energy_resonance_saturates_unitarity(self):
        p = yukawa_point(1.6798077622)
        v = np.array([3.0e-4])
        sigma, _ = S.cross_sections(p, v)
        _, k2, _ = S.kinematics(p.M3, v)
        self.assertGreater(sigma[0] * k2[0] / (8.0 * math.pi), 0.98)
        born, _ = S.cross_sections(p, v, yukawa="off")
        self.assertLess(born[0], 1e-3 * sigma[0])        # unitarised Born misses it

    def test_auto_switch(self):
        v = np.array([1e-3])
        weak, strong = yukawa_point(1e-3), yukawa_point(0.5)
        np.testing.assert_array_equal(S.cross_sections(weak, v)[0], S.cross_sections(weak, v, yukawa="off")[0])
        np.testing.assert_array_equal(S.cross_sections(strong, v)[0], S.cross_sections(strong, v, yukawa="on")[0])
        self.assertNotEqual(S.cross_sections(strong, v)[0][0], S.cross_sections(strong, v, yukawa="off")[0][0])


class AverageTests(unittest.TestCase):
    def test_constant_cross_section_average(self):
        p = contact_point()
        avg = S.maxwell_averages(p, [30.0], [0, 1, 5])
        sigma0, _ = S.cross_sections(p, [1e-6])
        for value in avg.values():
            self.assertAlmostEqual(value / sigma0[0], 1.0, places=6)

    def test_maxwell_weights_are_normalised(self):
        nodes, weights = S.velocity_quadrature(contact_point(), 1e-9, 0.01)
        disp = 100.0 / S.C_KMS
        for n in (0, 1, 5):
            numeric = np.sum(weights * nodes ** (2 + n) * np.exp(-nodes**2 / (4 * disp**2)))
            # int_0^inf v^(2+n) exp(-v^2/(4 d^2)) dv = (2d)^(3+n) Gamma((3+n)/2)/2
            exact = (2 * disp) ** (3 + n) * math.gamma((3 + n) / 2) / 2
            self.assertAlmostEqual(numeric / exact, 1.0, places=10)

    def test_resonant_average_against_brute_force(self):
        p = dark_h2_point(M3=50.0, dM=1.25e-5, K233=1.0, gamma2_other=1e-9)
        vr, dv = p.resonance(2)
        disp_kms = 200.0
        avg, err = S.maxwell_averages(p, [disp_kms], [1], n_angle=16, return_error=True)
        self.assertLess(err, 1e-6)
        disp = disp_kms / S.C_KMS
        v = np.unique(np.concatenate([
            np.geomspace(1e-9, 14 * disp, 20000),
            vr + dv * np.sinh(np.linspace(-14, 14, 60001)),
        ]))
        v = v[(v > 0) & (v < 14 * disp)]
        _, sv = S.cross_sections(p, v, n_angle=16)
        f = v**3 * np.exp(-v**2 / (4 * disp**2))
        trap = getattr(np, "trapezoid", None) or np.trapz  # NumPy < 2 has only trapz
        brute = trap(sv * f, v) / trap(f, v)
        self.assertAlmostEqual(avg[(disp_kms, 1)] / brute, 1.0, places=4)


class ScanInterfaceTests(unittest.TestCase):
    ROW = {"M2": "100.0001", "M3": "50", "vs": "500", "vx": "0", "a12": "0.001",
           "lX": "0.1", "lPhiX": "0.0", "lSX": "0.004", "K133": "nan", "K233": "nan",
           "w1": "0.0041", "w2": "1e-6", "dm": "True"}

    def test_from_scan_row_subtracts_invisible_width(self):
        p = S.SIDMPoint.from_scan_row(self.ROW)
        K133, K233 = S.vxzero_portal_couplings(0.0, 0.004, 500.0, 0.001)
        self.assertAlmostEqual(p.K233, K233)
        inv = S.scalar_to_identical_scalar_width(50.0, 100.0001, K233)
        self.assertAlmostEqual(p.gamma2_other, 1e-6 - inv, places=15)

    def test_from_scan_row_rejects_vx_branch(self):
        row = dict(self.ROW, vx="10")
        with self.assertRaises(ValueError):
            S.SIDMPoint.from_scan_row(row)

    def _write(self, path, rows, opener=open):
        with opener(path, "wt", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(self.ROW), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def test_scan_command_round_trip(self):
        with tempfile.TemporaryDirectory() as tmp:
            src, dst = Path(tmp) / "scan.dat", Path(tmp) / "scan_sidm.dat"
            self._write(src, [self.ROW, dict(self.ROW, w2="nan")])
            S.main(["scan", str(src), "--output", str(dst), "--velocities", "30", "1000",
                    "--dispersions", "50", "--powers", "1"])
            with dst.open() as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual([r["sidm_status"] for r in rows], ["ok", "failed"])
            self.assertEqual(rows[0]["dm"], "True")
            self.assertEqual(rows[0]["sidm_regime"], "perturbative")
            self.assertEqual(rows[0]["sidm_reliable"], "True")
            self.assertGreater(float(rows[0]["sidm_sigmaV_over_m_v30kms"]), 0.0)
            self.assertIn("sidm_sigmaV_over_m_avg_n1_disp50kms", rows[0])
            self.assertEqual(rows[1]["sidm_sigma_over_m_v30kms"], "nan")
            meta = json.loads((Path(tmp) / "scan_sidm.dat.sidm.json").read_text())
            self.assertEqual(meta["rows_ok"], 1)
            self.assertEqual(meta["method"], S.SIDM_METHOD)

    def test_gzip_input_and_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            src, dst = Path(tmp) / "scan.tsv.gz", Path(tmp) / "out.tsv.gz"
            self._write(src, [self.ROW], opener=gzip.open)
            S.main(["scan", str(src), "--output", str(dst), "--velocities", "30", "--dispersions"])
            header, rows = S.read_scan(dst)
            self.assertIn("sidm_sigmaV_over_m_v30kms", header)
            self.assertEqual(rows[0]["sidm_status"], "ok")


@unittest.skipUnless(os.environ.get("TRSM_SIDM_MAIN"), "set TRSM_SIDM_MAIN to a compiled TRSM/main_sidm")
class CalcHEPComparisonTest(unittest.TestCase):
    def test_against_calchep(self):
        import trsm_sidm_calchep_check
        self.assertEqual(trsm_sidm_calchep_check.main(["--main-sidm", os.environ["TRSM_SIDM_MAIN"]]), 0)


if __name__ == "__main__":
    unittest.main()
