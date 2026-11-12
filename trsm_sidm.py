"""Dark-matter self-scattering H3 H3 -> H3 H3 in the vx=0 TRSM.

This module is standalone: it needs only NumPy.  It computes, for one TRSM
point or for every row of a scan TSV,

* the total self-scattering cross section sigma(v), and
* the viscosity cross section sigma_V(v) used in astrophysical SIDM studies,

at fixed relative velocities and as Maxwell-Boltzmann averages over a halo
velocity dispersion.  It applies no astrophysical bound; that is done by
``trsm_sidm_bounds.py``.

Conventions
-----------
Potential coefficients follow ``trsm_vxzero_canonical_v1`` (see
``generate_trsm_info.vxzero_portal_couplings``)::

    V  >  K133 h1 h3 h3 + K233 h2 h3 h3 + (lX/4) h3^4,     h3 = X.

The Feynman rules are ``-2 i K_i33`` and ``-6 i lX``, as in
``DM/models/h4GOn/lgrng1.mdl``.  The tree-level amplitude is

    M = -6 lX - 4 sum_{i=1,2} K_i33^2 [ 1/D_i(s) + 1/(t - M_i^2) + 1/(u - M_i^2) ].

Velocities are relative velocities ``v = 2 beta`` in units of c (km/s in the
CLI and in column names): ``s = 4 M3^2/(1 - v^2/4)``, ``k ~ M3 v/2``.
Cross sections for identical particles::

    dsigma/dOmega = |M|^2 / (64 pi^2 s)
    sigma   = (1/2) int dOmega dsigma/dOmega
    sigma_V = (3/2) (1/2) int dOmega sin^2(theta) dsigma/dOmega

so sigma_V = sigma for isotropic scattering.  This is the normalisation of
Yang & Yu (arXiv:2205.03392, Eq. 2.5); the sigma_V of Colquhoun et al.
(arXiv:2011.04679, Eq. 28) is 2/3 of it.  Halo averages use the relative
velocity distribution ``f(v) ~ v^2 exp(-v^2/(4 sigma_1D^2))``, for which
``<v> = 4 sigma_1D/sqrt(pi)``, and report ``<sigma_V v^n>/<v^n>``.  n = 1 is
the ``<sigma v>/<v>`` of Sagunski et al. (arXiv:2006.12515) and n = 5 the
effective cross section of Yang & Yu.

What the default mode (``kmatrix``) does
----------------------------------------
1. s-channel propagators.  ``D_i(s) = s - M_i^2 + i M_i Gamma_i,other
   + Pi_i(s) - Re Pi_i(M_i^2)`` with the H3-loop self-energy
   ``Pi_i(s) = K_i33^2/(8 pi^2) [-beta ln((1+beta)/(1-beta)) + i pi beta]``.
   Its imaginary part is the energy-dependent H_i -> H3 H3 width
   ``K_i33^2 beta(s)/(8 pi)``.  Its real part is subtracted on shell, so
   ``M_i`` stays the physical mass.  For ``M_i < 2 M3`` the subtraction
   leaves ``+K_i33^2 betabar_i/(8 pi)`` at threshold, with
   ``betabar_i = sqrt(4 M3^2/M_i^2 - 1)``: the continuation of the same
   threshold cusp.  It is a leading effect near ``M_i = 2 M3`` and keeps the
   pole of the amplitude at ``s = M_i^2``.
2. Unitarity.  Every even partial wave ``a_l = (1/32 pi) int dx P_l(x) M``
   is unitarised, ``a_l = K_l/(1 - i (beta/2) K_l)``, so that
   ``sigma_l <= 8 pi (2l+1)/k^2``.
3. Non-perturbative exchange.  t/u-channel H_i exchange between slow H3's
   is an attractive Yukawa potential ``V = -alpha_i exp(-M_i r)/r`` with
   ``alpha_i = K_i33^2/(4 pi M3^2)``.  The Born approximation needs
   ``b_i = alpha_i M3/M_i << 1``; bound states appear from b = 1.68.  When
   ``max b_i`` exceeds ``yukawa_threshold`` the Born partial waves of the
   exchange are replaced by phase shifts from the radial Schroedinger
   equation: ``K_l = K_l^tree - K_l^Born-exchange + (2/beta) tan(delta_l)``.
   The contact and s-channel pieces are kept additively in the K-matrix;
   their distortion by the potential (Sommerfeld-type factors) is not
   included.

Validity (column ``sidm_regime``)
---------------------------------
``perturbative``     Born exchange is adequate (max b_i < yukawa_threshold).
``nonperturbative``  Schroedinger phase shifts used, alpha_i < alpha_strong.
``strong``           alpha_i >= alpha_strong (default 0.5, K_i33 = 2.5 M3).
                     The non-relativistic potential picture and the tree
                     couplings are then not under control; numbers are
                     written but ``sidm_reliable`` is False.

Two tree-level modes exist for cross-checks: ``tree-running-width`` (Dyson
propagator D_i above, no unitarisation, no Schroedinger solver) and
``tree-fixed-width`` (``D_i = s - M_i^2 + i M_i Gamma_i,total``, as CalcHEP).

Inputs from a scan row: M2, M3, lX, K133, K233 (or lPhiX, lSX, vs, a12), and
the physical widths w1, w2 with h1_h3h3_width/h2_h3h3_width (recomputed from
K_i33 when absent), from which Gamma_i,other = w_i - Gamma(h_i -> h3 h3).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import sys
import tempfile
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from numpy.polynomial import legendre as _legendre


SIDM_METHOD = "trsm_sidm_v2"
SIDM_CONVENTION = "trsm_vxzero_canonical_v1"

GEV2_TO_CM2 = 0.3893793721e-27  # (hbar c)^2 [GeV^2 cm^2]
GEV_TO_GRAM = 1.78266192e-24
C_KMS = 299792.458

DEFAULT_M1 = 125.09
DEFAULT_SM_VEV = 246.0  # generate_trsm_info.TRSM_SM_VEV
SM_HIGGS_WIDTH_125 = 4.07e-3  # GeV, only for the CLI default of Gamma_1,other

AMPLITUDE_MODES = ("kmatrix", "tree-running-width", "tree-fixed-width")
DEFAULT_MODE = "kmatrix"
DEFAULT_N_ANGLE = 48          # Gauss nodes on cos(theta) in [0, 1]
DEFAULT_LMAX = 8
DEFAULT_VELOCITIES_KMS = (10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0)
DEFAULT_DISPERSIONS_KMS = (10.0, 30.0, 100.0, 300.0, 1000.0)
DEFAULT_POWERS = (1, 5)
DEFAULT_YUKAWA_THRESHOLD = 0.01   # Born parameter above which the solver is used
DEFAULT_ALPHA_STRONG = 0.5
YUKAWA_CHOICES = ("auto", "off", "on")

# Maxwell-average quadrature.
_PANEL_NODES = 16
_LOG_PANEL_RATIO = 2.0
_RESONANCE_SCALES = (0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0)
_RESONANCE_REACH = 1000.0
_V_LOW_FACTOR = 1.0e-5
_V_HIGH_FACTOR = 14.0

# Radial Schroedinger solver.
_SOLVER_ETA = 0.05            # phase advance per RK4 step
_SOLVER_RANGE = 30.0          # match at r = _SOLVER_RANGE / (lightest mediator mass)
_SOLVER_START = 1.0e-12       # target error of the starting value
# Rough b = alpha M3/M_i at which the first bound state appears in partial wave l.
_BOUND_STATE_B = {0: 1.68, 2: 21.9, 4: 64.0, 6: 125.0, 8: 200.0}


# ---------------------------------------------------------------------------
# Point definition
# ---------------------------------------------------------------------------

def vxzero_portal_couplings(lPhiX, lSX, vs, a12, v=DEFAULT_SM_VEV):
    """Canonical-v1 K133, K233 (same formula as generate_trsm_info)."""
    c12 = math.cos(a12)
    s12 = math.sin(a12)
    k133 = 0.5 * (lPhiX * v * c12 - lSX * vs * s12)
    k233 = 0.5 * (lPhiX * v * s12 + lSX * vs * c12)
    return k133, k233


def scalar_to_identical_scalar_width(daughter_mass, parent_mass, coupling):
    """Gamma(parent -> daughter daughter) for V containing K parent*daughter^2."""
    if parent_mass <= 2.0 * daughter_mass:
        return 0.0
    beta = math.sqrt(max(0.0, 1.0 - 4.0 * daughter_mass**2 / parent_mass**2))
    return coupling**2 * beta / (8.0 * math.pi * parent_mass)


def resonance_velocity(M3, Mi):
    """Relative velocity v (units of c) at which sqrt(s) = Mi; nan below threshold."""
    if Mi <= 2.0 * M3:
        return math.nan
    return 2.0 * math.sqrt((Mi - 2.0 * M3) * (Mi + 2.0 * M3)) / Mi


@dataclass(frozen=True)
class SIDMPoint:
    """Inputs needed for H3 H3 -> H3 H3 (masses and widths in GeV)."""

    M3: float
    M2: float
    K133: float
    K233: float
    lX: float
    gamma1_other: float = 0.0
    gamma2_other: float = 0.0
    M1: float = DEFAULT_M1

    def __post_init__(self):
        for name, value in asdict(self).items():
            if not math.isfinite(float(value)):
                raise ValueError(f"{name} must be finite, got {value!r}")
        if self.M3 <= 0.0 or self.M2 <= 0.0 or self.M1 <= 0.0:
            raise ValueError("scalar masses must be positive")
        if self.gamma1_other < 0.0 or self.gamma2_other < 0.0:
            raise ValueError("widths must be non-negative")

    @property
    def mediators(self):
        """(mass, K_i33, Gamma_i,other) for i = 1, 2."""
        return (
            (self.M1, self.K133, self.gamma1_other),
            (self.M2, self.K233, self.gamma2_other),
        )

    def h3h3_width(self, index):
        Mi, Ki, _ = self.mediators[index - 1]
        return scalar_to_identical_scalar_width(self.M3, Mi, Ki)

    def total_width(self, index):
        return self.mediators[index - 1][2] + self.h3h3_width(index)

    def alpha(self, index):
        """Yukawa strength alpha_i = K_i33^2/(4 pi M3^2) of t-channel H_i exchange."""
        Ki = self.mediators[index - 1][1]
        return Ki**2 / (4.0 * math.pi * self.M3**2)

    def born_parameter(self, index):
        """b_i = alpha_i M3 / M_i; the Born approximation needs this << 1."""
        return self.alpha(index) * self.M3 / self.mediators[index - 1][0]

    def regime(self, yukawa_threshold=DEFAULT_YUKAWA_THRESHOLD, alpha_strong=DEFAULT_ALPHA_STRONG):
        """'perturbative', 'nonperturbative' or 'strong' (see module docstring)."""
        if max(self.alpha(1), self.alpha(2)) >= alpha_strong:
            return "strong"
        if max(self.born_parameter(1), self.born_parameter(2)) >= yukawa_threshold:
            return "nonperturbative"
        return "perturbative"

    def resonance(self, index):
        """(v_res, dv_res) in units of c for the H_i s-channel pole (nan if closed)."""
        Mi, Ki, gamma_other = self.mediators[index - 1]
        vr = resonance_velocity(self.M3, Mi)
        if not math.isfinite(vr):
            return math.nan, math.nan
        sigma_eff = Mi * gamma_other + Ki**2 * (0.5 * vr) / (8.0 * math.pi)
        ds_dv = 2.0 * self.M3**2 * vr / (1.0 - 0.25 * vr * vr) ** 2
        return vr, sigma_eff / ds_dv

    @classmethod
    def from_potential(
        cls, M2, M3, vs, a12, lX, lPhiX, lSX, *,
        w1=None, w2=None, gamma1_other=None, gamma2_other=None,
        M1=DEFAULT_M1, v=DEFAULT_SM_VEV,
    ):
        """Build a point from the scan basis (M2, M3, vs, a12, lX, lPhiX, lSX).

        Either physical total widths ``w1``/``w2`` (which include H_i -> H3 H3)
        or the non-H3H3 widths ``gamma*_other`` may be given.  Missing h1
        information defaults to 4.07 MeV * cos^2(a12); missing h2 information
        defaults to zero (the dark-H2 limit).
        """
        K133, K233 = vxzero_portal_couplings(lPhiX, lSX, vs, a12, v=v)
        g1 = _other_width(w1, gamma1_other, scalar_to_identical_scalar_width(M3, M1, K133))
        g2 = _other_width(w2, gamma2_other, scalar_to_identical_scalar_width(M3, M2, K233))
        if g1 is None:
            g1 = SM_HIGGS_WIDTH_125 * math.cos(a12) ** 2
        if g2 is None:
            g2 = 0.0
        return cls(M3=M3, M2=M2, K133=K133, K233=K233, lX=lX,
                   gamma1_other=g1, gamma2_other=g2, M1=M1)

    @classmethod
    def from_scan_row(cls, row, *, M1=DEFAULT_M1, v=DEFAULT_SM_VEV):
        """Build a point from a ``generate_trsm_points.py`` TSV row (dict of str)."""
        get = lambda key: _row_float(row, key)
        vx = get("vx")
        if math.isfinite(vx) and vx != 0.0:
            raise ValueError("only the vx=0 branch (stable H3) is supported")
        M2, M3, lX = get("M2"), get("M3"), get("lX")
        K133, K233 = get("K133"), get("K233")
        if not (math.isfinite(K133) and math.isfinite(K233)):
            K133, K233 = vxzero_portal_couplings(get("lPhiX"), get("lSX"), get("vs"), get("a12"), v=v)
        w1, w2 = get("w1"), get("w2")
        if not (math.isfinite(w1) and math.isfinite(w2)):
            raise ValueError("scan row needs finite physical widths w1 and w2")
        inv1, inv2 = get("h1_h3h3_width"), get("h2_h3h3_width")
        if not math.isfinite(inv1):
            inv1 = scalar_to_identical_scalar_width(M3, M1, K133)
        if not math.isfinite(inv2):
            inv2 = scalar_to_identical_scalar_width(M3, M2, K233)
        return cls(M3=M3, M2=M2, K133=K133, K233=K233, lX=lX,
                   gamma1_other=_subtract_width(w1, inv1, "w1"),
                   gamma2_other=_subtract_width(w2, inv2, "w2"), M1=M1)


def _other_width(total, other, invisible):
    if other is not None:
        return float(other)
    if total is not None:
        return _subtract_width(float(total), invisible, "total width")
    return None


def _subtract_width(total, invisible, label):
    other = total - invisible
    if other < -1.0e-6 * max(abs(total), 1.0e-300):
        raise ValueError(f"{label}={total:g} is smaller than its H3H3 part {invisible:g}")
    return max(0.0, other)


def _row_float(row, key):
    value = row.get(key)
    if value is None:
        return math.nan
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


# ---------------------------------------------------------------------------
# Kinematics and the tree amplitude
# ---------------------------------------------------------------------------

def kinematics(M3, v):
    """Return s, k^2 and beta for relative velocity v (array, units of c)."""
    v = np.asarray(v, dtype=float)
    if np.any(v < 0.0) or np.any(v >= 2.0):
        raise ValueError("relative velocity must satisfy 0 <= v < 2")
    q = 1.0 - 0.25 * v * v
    s = 4.0 * M3**2 / q
    k2 = M3**2 * 0.25 * v * v / q
    return s, k2, 0.5 * v


def _s_minus_mass2(M3, v, Mi):
    """s - Mi^2 without cancellation near the pole."""
    vr = resonance_velocity(M3, Mi)
    q = 1.0 - 0.25 * v * v
    if math.isfinite(vr):
        return M3**2 * (v - vr) * (v + vr) / (q * (1.0 - 0.25 * vr * vr))
    return (2.0 * M3 - Mi) * (2.0 * M3 + Mi) + M3**2 * v * v / q


def loop_real_part(beta):
    """Re of the H3 bubble above threshold, -beta ln((1+beta)/(1-beta))."""
    beta = np.asarray(beta, dtype=float)
    return -2.0 * beta * np.arctanh(beta)


def loop_real_part_at_mass(M3, Mi):
    """The same function at s = Mi^2, continued below threshold if Mi < 2 M3."""
    if Mi > 2.0 * M3:
        beta = math.sqrt((Mi - 2.0 * M3) * (Mi + 2.0 * M3)) / Mi
        return -2.0 * beta * math.atanh(beta)
    if Mi < 2.0 * M3:
        betabar = math.sqrt((2.0 * M3 - Mi) * (2.0 * M3 + Mi)) / Mi
        return -2.0 * betabar * math.atan2(1.0, betabar)
    return 0.0


def self_energy_real(point, v, index):
    """Re[Pi_i(s) - Pi_i(M_i^2)] of the H3 loop, in GeV^2 (on-shell subtracted)."""
    Mi, Ki, _ = point.mediators[index - 1]
    _, _, beta = kinematics(point.M3, v)
    return Ki**2 / (8.0 * math.pi**2) * (loop_real_part(beta) - loop_real_part_at_mass(point.M3, Mi))


def _s_channel_amplitude(point, v, mode):
    """s-channel H1 + H2 exchange (complex), shape (len(v),)."""
    v = np.asarray(v, dtype=float)
    _, _, beta = kinematics(point.M3, v)
    amp = np.zeros(v.shape, dtype=complex)
    for index, (Mi, Ki, gamma_other) in enumerate(point.mediators, start=1):
        if Ki == 0.0:
            continue
        real = _s_minus_mass2(point.M3, v, Mi)
        if mode == "kmatrix":
            real = real + self_energy_real(point, v, index)
            imag = Mi * gamma_other * np.ones_like(v)
        elif mode == "tree-running-width":
            real = real + self_energy_real(point, v, index)
            imag = Mi * gamma_other + Ki**2 * beta / (8.0 * math.pi)
        elif mode == "tree-fixed-width":
            imag = Mi * point.total_width(index) * np.ones_like(v)
        else:
            raise ValueError(f"unknown amplitude mode {mode!r}; use one of {AMPLITUDE_MODES}")
        denom = real + 1j * imag
        # An exact zero only happens for Gamma_other = 0 on the pole; a tiny
        # offset keeps K finite so that the K-matrix returns 2i/beta.
        denom = np.where(denom == 0.0, 1.0e-280, denom)
        amp = amp - 4.0 * Ki**2 / denom
    return amp


@lru_cache(maxsize=8)
def _unit_rule(n):
    t, w = _legendre.leggauss(int(n))
    t, w = 0.5 * (t + 1.0), 0.5 * w
    t.setflags(write=False)
    w.setflags(write=False)
    return t, w


def _angular_nodes(point, k2, n_angle):
    """Quadrature in cos(theta) for integrands even in cos(theta).

    Returns ``omx = 1 - cos(theta)`` and weights, both of shape (nv, n_angle),
    such that ``sum(w * g) = int_{-1}^{1} g(x) dx``.  The nodes cover
    x in (0, 1) and cluster logarithmically towards x = 1 when the lightest
    exchanged mass is small compared with k, where t-channel exchange peaks.
    """
    tau, wt = _unit_rule(n_angle)
    masses = [Mi for Mi, Ki, _ in point.mediators if Ki != 0.0]
    k2 = np.asarray(k2, dtype=float)
    if masses:
        with np.errstate(divide="ignore"):
            c = np.minimum(min(masses) ** 2 / (2.0 * k2), 1.0e30)
    else:
        c = np.full(k2.shape, 1.0e30)
    rho = np.log1p(1.0 / c)
    omx = np.minimum(c[:, None] * np.expm1(rho[:, None] * tau[None, :]), 1.0)
    weights = 2.0 * wt[None, :] * (c * rho)[:, None] * np.exp(rho[:, None] * tau[None, :])
    return omx, weights


def _exchange_amplitude(point, k2, omx, index):
    """t + u exchange of mediator ``index`` (real), for omx = 1 - cos(theta)."""
    Mi, Ki, _ = point.mediators[index - 1]
    k2 = np.asarray(k2, dtype=float)[:, None]
    t = -2.0 * k2 * omx
    u = -2.0 * k2 * (2.0 - omx)
    return -4.0 * Ki**2 * (1.0 / (t - Mi**2) + 1.0 / (u - Mi**2))


def _even_ls(lmax):
    return np.arange(0, int(lmax) + 1, 2)


def _even_legendre(x, lmax):
    """P_l(x) for even l <= lmax; appends an axis of length n_even."""
    return _legendre.legvander(np.asarray(x, dtype=float), int(lmax))[..., 0::2]


def _project(values, weights, legendre):
    """a_l = (1/32 pi) int dx P_l(x) values(x) for even l, shape (nv, n_even)."""
    return np.einsum("vj,vjl->vl", values * weights, legendre) / (32.0 * math.pi)


# ---------------------------------------------------------------------------
# Radial Schroedinger equation for the exchange (Yukawa) potential
# ---------------------------------------------------------------------------

def _double_factorial(n):
    out = 1.0
    while n > 1:
        out *= n
        n -= 2
    return out


def _spherical_j(l, z):
    z = np.asarray(z, dtype=float)
    out = np.empty_like(z)
    small = z < l + 2.0
    if np.any(small):
        zs = z[small]
        q = -0.5 * zs * zs
        term = np.ones_like(zs)
        total = term.copy()
        for n in range(1, 60):
            term = term * q / (n * (2 * l + 2 * n + 1))
            total += term
        out[small] = zs**l / _double_factorial(2 * l + 1) * total
    if np.any(~small):
        zl = z[~small]
        jm = np.sin(zl) / zl
        if l > 0:
            jc = np.sin(zl) / zl**2 - np.cos(zl) / zl
            for n in range(1, l):
                jm, jc = jc, (2 * n + 1) / zl * jc - jm
            jm = jc
        out[~small] = jm
    return out


def _spherical_y(l, z):
    z = np.asarray(z, dtype=float)
    ym = -np.cos(z) / z
    if l == 0:
        return ym
    yc = -np.cos(z) / z**2 - np.sin(z) / z
    for n in range(1, l):
        ym, yc = yc, (2 * n + 1) / z * yc - ym
    return yc


def riccati_bessel(l, z):
    """Return (j, j', n, n') for the Riccati-Bessel functions z j_l(z), z y_l(z)."""
    z = np.asarray(z, dtype=float)
    jl, yl = _spherical_j(l, z), _spherical_y(l, z)
    if l == 0:
        return z * jl, np.cos(z), z * yl, np.sin(z)
    return (z * jl, z * _spherical_j(l - 1, z) - l * jl,
            z * yl, z * _spherical_y(l - 1, z) - l * yl)


def yukawa_tan_delta(M3, alphas, masses, k, ls, *, eta=_SOLVER_ETA,
                     range_factor=_SOLVER_RANGE, start_tolerance=_SOLVER_START):
    """Phase shifts of ``V(r) = -sum_i alpha_i exp(-m_i r)/r`` for reduced mass M3/2.

    Returns (N, D) of shape (len(k), len(ls)) with ``tan(delta_l) = N/D``.
    The radial equation ``u'' + [k^2 - l(l+1)/r^2 + M3 sum_i alpha_i
    exp(-m_i r)/r] u = 0`` is written for ``w = u/sqrt(r)`` in ``t = ln r``,
    ``w'' = -Q w`` with ``Q = r^2 (k^2 + M3 sum ...) - (l + 1/2)^2``.  It is
    integrated with RK4 in a Pruefer phase, ``w'/w = kappa cot(phi)``, whose
    scale ``kappa`` is constant within a step and reset at every step to the
    local frequency ``sqrt(max(|Q|, |dQ/dt|^(2/3), 1))``.  The phase stays
    finite through nodes, turning points and the centrifugal barrier, and
    small k costs nothing.  The solution is matched to Riccati-Bessel
    functions at ``r = range_factor / min(m_i)``.
    """
    k = np.atleast_1d(np.asarray(k, dtype=float))
    alphas = np.asarray(alphas, dtype=float)
    masses = np.asarray(masses, dtype=float)
    if np.any(k <= 0.0):
        raise ValueError("the Schroedinger solver needs k > 0")
    if alphas.size == 0 or np.any(alphas <= 0.0) or np.any(masses <= 0.0):
        raise ValueError("need positive couplings and mediator masses")
    r_match = range_factor / masses.min()
    t_end = math.log(r_match)
    k2 = k * k
    k2_max, k2_min = float(k2.max()), float(k2.min())
    third = 1.0 / 3.0

    def potential(t):
        r = math.exp(t)
        e = alphas * np.exp(-masses * r)
        return r * r, M3 * r * float(e.sum()), M3 * float(((r - masses * r * r) * e).sum())

    N = np.empty((k.size, len(ls)))
    D = np.empty_like(N)
    for column, l in enumerate(ls):
        nu2 = (l + 0.5) ** 2
        # Start deep inside the centrifugal barrier, where w ~ r^(l+1/2):
        # r^2 (k^2 + U) <= eps nu^2.  The regular solution is an attractor, so
        # the error of the starting value decays like eps^(2l+2).
        eps = start_tolerance ** (1.0 / (2 * l + 2))
        t_start = math.log(min(math.sqrt(0.5 * eps * nu2 / k2_max),
                               0.5 * eps * nu2 / (M3 * alphas.sum()), 1.0e-2 * r_match))

        # Step grid: resolve the fastest oscillation, the turning points (Airy
        # scale |dQ/dt|^(-1/3)) and the relaxation rate in the barrier.
        grid = [t_start]
        while grid[-1] < t_end:
            r2, U, dU = potential(grid[-1])
            oscillation = max(r2 * k2_max + U - nu2, 0.0)
            barrier = max(-(r2 * k2_min + U - nu2), 0.0)
            slope = max(abs(2.0 * r2 * k2_max + dU), abs(2.0 * r2 * k2_min + dU))
            rate = max(math.sqrt(oscillation), math.sqrt(barrier), slope**third, 1.0)
            grid.append(min(grid[-1] + eta / rate, t_end))

        pot_a = potential(grid[0])
        Q = pot_a[0] * k2 + (pot_a[1] - nu2)
        kappa = np.sqrt(np.maximum(-Q, 1.0e-12))
        phi = np.full(k.size, 0.25 * math.pi)            # w'/w = sqrt(-Q) in the barrier
        for ta, tb in zip(grid[:-1], grid[1:]):
            h = tb - ta
            pot_m, pot_b = potential(ta + 0.5 * h), potential(tb)
            Qa = pot_a[0] * k2 + (pot_a[1] - nu2)
            dQa = 2.0 * pot_a[0] * k2 + pot_a[2]
            scale = np.sqrt(np.maximum(np.maximum(np.abs(Qa), np.abs(dQa) ** (2.0 * third)), 1.0))
            phi = np.arctan2(scale * np.sin(phi), kappa * np.cos(phi))   # keeps w'/w
            kappa = scale
            qa, qm, qb = Qa / kappa, (pot_m[0] * k2 + (pot_m[1] - nu2)) / kappa, \
                (pot_b[0] * k2 + (pot_b[1] - nu2)) / kappa

            def rhs(q, angle):
                c = np.cos(angle)
                return kappa * c * c + q * (1.0 - c * c)

            f1 = rhs(qa, phi)
            f2 = rhs(qm, phi + 0.5 * h * f1)
            f3 = rhs(qm, phi + 0.5 * h * f2)
            f4 = rhs(qb, phi + h * f3)
            phi = phi + (h / 6.0) * (f1 + 2.0 * f2 + 2.0 * f3 + f4)
            pot_a = pot_b

        z = k * r_match
        j, dj, n, dn = riccati_bessel(l, z)
        s, c = np.sin(phi), np.cos(phi)
        # r u'/u = 1/2 + kappa cot(phi); multiplied through by sin(phi).
        N[:, column] = (z * dj - 0.5 * j) * s - kappa * j * c
        D[:, column] = (z * dn - 0.5 * n) * s - kappa * n * c
    return N, D


def _solver_mediators(point):
    """Indices of mediators kept in the potential (the rest stay in Born)."""
    born = [point.born_parameter(1), point.born_parameter(2)]
    top = max(born)
    return [i + 1 for i, b in enumerate(born) if b > 0.0 and b >= 1.0e-6 * top]


def _use_solver(point, mode, yukawa, yukawa_threshold):
    if mode != "kmatrix" or yukawa == "off":
        return False
    if yukawa not in YUKAWA_CHOICES:
        raise ValueError(f"yukawa must be one of {YUKAWA_CHOICES}")
    top = max(point.born_parameter(1), point.born_parameter(2))
    return top > 0.0 and (yukawa == "on" or top >= yukawa_threshold)


# ---------------------------------------------------------------------------
# Partial waves and cross sections
# ---------------------------------------------------------------------------

def _waves(point, v, mode, lmax, n_angle, yukawa, yukawa_threshold):
    """Shared core: angular nodes, Born projections and final partial waves."""
    v = np.atleast_1d(np.asarray(v, dtype=float))
    s, k2, beta = kinematics(point.M3, v)
    omx, weights = _angular_nodes(point, k2, n_angle)
    legendre = _even_legendre(1.0 - omx, lmax)
    ls = _even_ls(lmax)

    exchange = {i: _exchange_amplitude(point, k2, omx, i)
                for i in (1, 2) if point.mediators[i - 1][1] != 0.0}
    background = np.full(omx.shape, -6.0 * point.lX)
    for values in exchange.values():
        background = background + values
    born = _project(background, weights, legendre)          # contact + exchange, real

    K = born.astype(complex)
    K[:, 0] += _s_channel_amplitude(point, v, mode) / (16.0 * math.pi)
    exact = np.zeros(ls.size, dtype=bool)
    if mode != "kmatrix":
        a = K
    elif not _use_solver(point, mode, yukawa, yukawa_threshold):
        a = K / (1.0 - 0.5j * beta[:, None] * K)
    else:
        if np.any(v <= 0.0):
            raise ValueError("the Schroedinger solver needs v > 0")
        kept = _solver_mediators(point)
        born_exchange = _project(sum(exchange[i] for i in kept), weights, legendre)
        top = max(point.born_parameter(i) for i in kept)
        size = 0.5 * np.max(beta[:, None] * np.abs(born_exchange), axis=0)
        for column, l in enumerate(ls):
            exact[column] = (l == 0 or top >= 0.03 * _BOUND_STATE_B.get(int(l), 2.4 * (l + 1) ** 2)
                             or size[column] > 1.0e-3)
        a = K / (1.0 - 0.5j * beta[:, None] * K)
        cols = np.flatnonzero(exact)
        N, D = yukawa_tan_delta(
            point.M3, [point.alpha(i) for i in kept], [point.mediators[i - 1][0] for i in kept],
            np.sqrt(k2), ls[cols])
        rest = K[:, cols] - born_exchange[:, cols]          # contact, s-channel, dropped mediators
        top_part = rest * D + (2.0 / beta)[:, None] * N
        a[:, cols] = top_part / (D - 1j * N - 0.5j * beta[:, None] * rest * D)
    return {"v": v, "s": s, "k2": k2, "beta": beta, "omx": omx, "weights": weights,
            "legendre": legendre, "ls": ls, "background": background, "born": born,
            "K": K, "a": a, "exact": exact}


def partial_waves(point, v, *, mode=DEFAULT_MODE, lmax=DEFAULT_LMAX, n_angle=DEFAULT_N_ANGLE,
                  yukawa="auto", yukawa_threshold=DEFAULT_YUKAWA_THRESHOLD):
    """Return (ls, K_l, a_l) for even l, arrays of shape (len(v), n_even).

    ``K_l`` are the projections of the tree amplitude in the chosen mode and
    ``a_l`` the amplitudes actually used.  ``sigma = 8 pi/s sum (2l+1)|a_l|^2``
    for l <= lmax; higher partial waves are kept at Born level.
    """
    w = _waves(point, v, mode, lmax, n_angle, yukawa, yukawa_threshold)
    return w["ls"], w["K"], w["a"]


def _amplitude_on_nodes(w):
    """Full amplitude on the internal angular nodes, shape (nv, n_angle)."""
    correction = 16.0 * math.pi * (2 * w["ls"] + 1)[None, :] * (w["a"] - w["born"])
    return w["background"] + np.einsum("vl,vjl->vj", correction, w["legendre"])


def scattering_amplitude(point, v, x, *, mode=DEFAULT_MODE, lmax=DEFAULT_LMAX,
                         n_angle=DEFAULT_N_ANGLE, yukawa="auto",
                         yukawa_threshold=DEFAULT_YUKAWA_THRESHOLD):
    """Full amplitude M(v, x) with x = cos(theta) in the CM frame, shape (nv, nx)."""
    w = _waves(point, v, mode, lmax, n_angle, yukawa, yukawa_threshold)
    x = np.atleast_1d(np.asarray(x, dtype=float))
    omx = np.broadcast_to(1.0 - x, (w["v"].size, x.size))
    background = np.full(omx.shape, -6.0 * point.lX)
    for index in (1, 2):
        if point.mediators[index - 1][1] != 0.0:
            background = background + _exchange_amplitude(point, w["k2"], omx, index)
    correction = 16.0 * math.pi * (2 * w["ls"] + 1)[None, :] * (w["a"] - w["born"])
    return background + correction @ _even_legendre(x, lmax).T


def cross_sections(point, v, *, mode=DEFAULT_MODE, lmax=DEFAULT_LMAX, n_angle=DEFAULT_N_ANGLE,
                   yukawa="auto", yukawa_threshold=DEFAULT_YUKAWA_THRESHOLD):
    """Return (sigma, sigma_V) in GeV^-2 at relative velocities v (units of c)."""
    w = _waves(point, v, mode, lmax, n_angle, yukawa, yukawa_threshold)
    amp2 = np.abs(_amplitude_on_nodes(w)) ** 2
    omx, weights = w["omx"], w["weights"]
    sigma = np.sum(amp2 * weights, axis=1) / (64.0 * math.pi * w["s"])
    sin2 = omx * (2.0 - omx)
    sigma_v = np.sum(amp2 * weights * sin2, axis=1) * 3.0 / (128.0 * math.pi * w["s"])
    return sigma, sigma_v


def differential_cross_section(point, v, x, **options):
    """dsigma/dOmega = |M|^2/(64 pi^2 s) in GeV^-2 (sigma is half its 4pi integral)."""
    v = np.atleast_1d(np.asarray(v, dtype=float))
    s, _, _ = kinematics(point.M3, v)
    amp2 = np.abs(scattering_amplitude(point, v, x, **options)) ** 2
    return amp2 / (64.0 * math.pi**2 * s[:, None])


def tree_unitarity_ratio(point, v, *, lmax=DEFAULT_LMAX, n_angle=DEFAULT_N_ANGLE):
    """max_l (beta/2)|K_l| of the tree amplitude; > 1 violates elastic unitarity."""
    v = np.atleast_1d(np.asarray(v, dtype=float))
    _, K, _ = partial_waves(point, v, mode="tree-running-width", lmax=lmax, n_angle=n_angle)
    _, _, beta = kinematics(point.M3, v)
    return np.max(0.5 * beta[:, None] * np.abs(K), axis=1)


def per_mass_cm2_g(sigma_gev2, M3):
    """Convert sigma [GeV^-2] to sigma/m [cm^2/g] for DM mass M3 [GeV]."""
    return np.asarray(sigma_gev2) * GEV2_TO_CM2 / (M3 * GEV_TO_GRAM)


def unitarity_limit_cm2_g(M3, v_kms, l=0):
    """Largest sigma/m of one partial wave, 8 pi (2l+1)/(k^2 M3), in cm^2/g."""
    _, k2, _ = kinematics(M3, np.asarray(v_kms, dtype=float) / C_KMS)
    return per_mass_cm2_g(8.0 * math.pi * (2 * l + 1) / k2, M3)


# ---------------------------------------------------------------------------
# Halo (Maxwell-Boltzmann) averages
# ---------------------------------------------------------------------------

def _panel(a, b, n, mapping, centre=0.0, width=1.0):
    t, w = _legendre.leggauss(n)
    if mapping == "linear":
        nodes = a + 0.5 * (b - a) * (t + 1.0)
        return nodes, 0.5 * (b - a) * w
    if mapping == "log":
        la, lb = math.log(a), math.log(b)
        nodes = np.exp(la + 0.5 * (lb - la) * (t + 1.0))
        return nodes, 0.5 * (lb - la) * w * nodes
    if mapping == "tan":
        ta, tb = math.atan((a - centre) / width), math.atan((b - centre) / width)
        theta = ta + 0.5 * (tb - ta) * (t + 1.0)
        nodes = centre + width * np.tan(theta)
        return nodes, 0.5 * (tb - ta) * w * width / np.cos(theta) ** 2
    raise ValueError(mapping)


def velocity_quadrature(point, v_low, v_high, *, panel_nodes=_PANEL_NODES):
    """Nodes/weights on [0, v_high] resolving both s-channel resonances."""
    breaks = {0.0, v_low, v_high}
    edge = v_low
    while edge < v_high:
        breaks.add(edge)
        edge *= _LOG_PANEL_RATIO
    resonances = []
    for index in (1, 2):
        vr, dv = point.resonance(index)
        if math.isfinite(vr) and dv > 0.0 and vr < v_high:
            resonances.append((vr, dv))
            breaks.add(vr)
            for scale in _RESONANCE_SCALES:
                for edge in (vr - scale * dv, vr + scale * dv):
                    if 0.0 < edge < v_high:
                        breaks.add(edge)
    edges = np.array(sorted(breaks))
    nodes, weights = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        if b <= a:
            continue
        mid = 0.5 * (a + b)
        near = [(abs(mid - vr) / dv, vr, dv) for vr, dv in resonances
                if abs(mid - vr) <= _RESONANCE_REACH * dv]
        if near:
            _, vr, dv = min(near)
            n, w = _panel(a, b, panel_nodes, "tan", vr, dv)
        elif a == 0.0:
            n, w = _panel(a, b, panel_nodes, "linear")
        else:
            n, w = _panel(a, b, panel_nodes, "log")
        nodes.append(n)
        weights.append(w)
    return np.concatenate(nodes), np.concatenate(weights)


def maxwell_averages(
    point, dispersions_kms=DEFAULT_DISPERSIONS_KMS, powers=DEFAULT_POWERS, *,
    panel_nodes=_PANEL_NODES, return_error=False, **options,
):
    """<sigma_V v^n>/<v^n> [GeV^-2] for Maxwellian halos.

    ``dispersions_kms`` are 1D velocity dispersions sigma_1D; the mean
    relative velocity of the distribution is 4 sigma_1D/sqrt(pi).  Returns a
    dict keyed by (dispersion_kms, n).  With ``return_error`` the quadrature
    is repeated with twice the nodes per panel and the largest relative
    difference is returned as well.
    """
    dispersions_kms = [float(d) for d in dispersions_kms]
    disps = np.asarray(dispersions_kms) / C_KMS
    if np.any(disps <= 0.0):
        raise ValueError("velocity dispersions must be positive")
    v_low = _V_LOW_FACTOR * float(np.min(disps))
    v_high = min(_V_HIGH_FACTOR * float(np.max(disps)), 1.9)

    def evaluate(n_nodes):
        nodes, weights = velocity_quadrature(point, v_low, v_high, panel_nodes=n_nodes)
        _, sigma_v = cross_sections(point, nodes, **options)
        out = {}
        for disp_kms, disp in zip(dispersions_kms, disps):
            f = nodes**2 * np.exp(-(nodes**2) / (4.0 * disp**2))
            for n in powers:
                wv = weights * f * nodes**n
                out[(disp_kms, n)] = float(np.dot(wv, sigma_v) / np.sum(wv))
        return out

    result = evaluate(panel_nodes)
    if not return_error:
        return result
    refined = evaluate(2 * panel_nodes)
    error = max(
        (abs(refined[key] - value) / abs(refined[key]) if refined[key] != 0.0 else abs(value))
        for key, value in result.items()
    )
    return refined, error


def mean_relative_velocity_to_dispersion(mean_v):
    """sigma_1D of the Maxwellian whose mean relative velocity is ``mean_v``."""
    return mean_v * math.sqrt(math.pi) / 4.0


# ---------------------------------------------------------------------------
# Per-point analysis and output columns
# ---------------------------------------------------------------------------

def _tag(value):
    return f"{value:g}".replace(".", "p").replace("+", "")


def velocity_columns(velocities_kms):
    cols = []
    for v in velocities_kms:
        cols += [f"sidm_sigma_over_m_v{_tag(v)}kms", f"sidm_sigmaV_over_m_v{_tag(v)}kms"]
    return cols


def average_columns(dispersions_kms, powers):
    return [f"sidm_sigmaV_over_m_avg_n{n}_disp{_tag(d)}kms" for d in dispersions_kms for n in powers]


DIAGNOSTIC_COLUMNS = (
    "sidm_status", "sidm_error", "sidm_method", "sidm_mode", "sidm_regime", "sidm_reliable",
    "sidm_alpha_h1", "sidm_alpha_h2", "sidm_born_h1", "sidm_born_h2",
    "sidm_v_res_h1_kms", "sidm_v_res_h2_kms", "sidm_dv_res_h1_kms", "sidm_dv_res_h2_kms",
    "sidm_gamma1_other_GeV", "sidm_gamma2_other_GeV",
    "sidm_tree_unitarity_max", "sidm_avg_quadrature_relerr",
)


def sidm_columns(velocities_kms=DEFAULT_VELOCITIES_KMS, dispersions_kms=DEFAULT_DISPERSIONS_KMS,
                 powers=DEFAULT_POWERS):
    return list(DIAGNOSTIC_COLUMNS) + velocity_columns(velocities_kms) + average_columns(dispersions_kms, powers)


def point_diagnostics(point, *, mode=DEFAULT_MODE, yukawa_threshold=DEFAULT_YUKAWA_THRESHOLD,
                      alpha_strong=DEFAULT_ALPHA_STRONG):
    """Regime and coupling diagnostics (no cross sections)."""
    regime = point.regime(yukawa_threshold, alpha_strong)
    out = {"sidm_status": "ok", "sidm_error": "", "sidm_method": SIDM_METHOD, "sidm_mode": mode,
           "sidm_regime": regime, "sidm_reliable": regime != "strong",
           "sidm_alpha_h1": point.alpha(1), "sidm_alpha_h2": point.alpha(2),
           "sidm_born_h1": point.born_parameter(1), "sidm_born_h2": point.born_parameter(2)}
    for index in (1, 2):
        vr, dv = point.resonance(index)
        out[f"sidm_v_res_h{index}_kms"] = vr * C_KMS
        out[f"sidm_dv_res_h{index}_kms"] = dv * C_KMS
    out["sidm_gamma1_other_GeV"] = point.gamma1_other
    out["sidm_gamma2_other_GeV"] = point.gamma2_other
    return out


def analyse_point(
    point, *, velocities_kms=DEFAULT_VELOCITIES_KMS, dispersions_kms=DEFAULT_DISPERSIONS_KMS,
    powers=DEFAULT_POWERS, mode=DEFAULT_MODE, lmax=DEFAULT_LMAX, n_angle=DEFAULT_N_ANGLE,
    yukawa="auto", yukawa_threshold=DEFAULT_YUKAWA_THRESHOLD, alpha_strong=DEFAULT_ALPHA_STRONG,
):
    """Return an ordered dict of SIDM output columns for one point (sigma/m in cm^2/g)."""
    options = {"mode": mode, "lmax": lmax, "n_angle": n_angle, "yukawa": yukawa,
               "yukawa_threshold": yukawa_threshold}
    out = point_diagnostics(point, mode=mode, yukawa_threshold=yukawa_threshold, alpha_strong=alpha_strong)
    vel = np.asarray(velocities_kms, dtype=float) / C_KMS
    if vel.size:
        sigma, sigma_v = cross_sections(point, vel, **options)
        out["sidm_tree_unitarity_max"] = float(np.max(tree_unitarity_ratio(point, vel, lmax=lmax, n_angle=n_angle)))
    else:
        sigma = sigma_v = np.array([])
        out["sidm_tree_unitarity_max"] = math.nan
    if len(dispersions_kms) and len(powers):
        averages, error = maxwell_averages(point, dispersions_kms, powers, return_error=True, **options)
    else:
        averages, error = {}, math.nan
    out["sidm_avg_quadrature_relerr"] = error

    names = velocity_columns(velocities_kms)
    for i in range(len(vel)):
        out[names[2 * i]] = float(per_mass_cm2_g(sigma[i], point.M3))
        out[names[2 * i + 1]] = float(per_mass_cm2_g(sigma_v[i], point.M3))
    keys = [(float(d), n) for d in dispersions_kms for n in powers]
    for name, key in zip(average_columns(dispersions_kms, powers), keys):
        out[name] = float(per_mass_cm2_g(averages[key], point.M3))
    return out


def failed_columns(error, columns, mode=DEFAULT_MODE):
    out = {name: math.nan for name in columns}
    out.update({"sidm_status": "failed", "sidm_error": str(error).replace("\t", " ").replace("\n", " "),
                "sidm_method": SIDM_METHOD, "sidm_mode": mode, "sidm_regime": "", "sidm_reliable": ""})
    return out


# ---------------------------------------------------------------------------
# Scan-file helpers (shared with trsm_sidm_bounds.py)
# ---------------------------------------------------------------------------

def format_value(value):
    if isinstance(value, (bool, np.bool_)):
        return "True" if value else "False"
    if isinstance(value, float):
        return "nan" if math.isnan(value) else repr(value)
    return str(value)


def open_text(path):
    """Open a plain or gzip-compressed text file for reading."""
    opener = gzip.open if Path(path).suffix == ".gz" else open
    return opener(path, "rt", encoding="utf-8", newline="")


def read_scan(path):
    with open_text(path) as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_scan(destination, header, rows, new_columns, evaluate, *, progress=0, label="sidm"):
    """Write ``rows`` plus ``new_columns`` computed by ``evaluate(row) -> dict``.

    The source is never modified; the output appears only when complete.
    Returns the number of rows for which ``evaluate`` did not raise.
    """
    destination = Path(destination)
    out_header = header + [c for c in new_columns if c not in header]
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=destination.name, dir=destination.parent)
    os.close(fd)
    n_ok = 0
    try:
        opener = gzip.open if destination.suffix == ".gz" else open
        with opener(tmp, "wt", encoding="utf-8", newline="") as stream:
            stream.write("\t".join(out_header) + "\n")
            for number, row in enumerate(rows, start=1):
                result, ok = evaluate(row)
                n_ok += bool(ok)
                merged = dict(row)
                merged.update({k: format_value(v) for k, v in result.items()})
                stream.write("\t".join(merged.get(c, "nan") or "nan" for c in out_header) + "\n")
                if progress and number % progress == 0:
                    print(f"[{label}] {number}/{len(rows)} rows", file=sys.stderr)
        os.replace(tmp, destination)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise
    return n_ok


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def configuration(args):
    return {
        "method": SIDM_METHOD, "convention": SIDM_CONVENTION, "mode": args.mode,
        "n_angle": args.n_angle, "lmax": args.lmax, "M1_GeV": args.m1,
        "sm_vev_for_missing_K_GeV": args.vev, "yukawa": args.yukawa,
        "yukawa_threshold": args.yukawa_threshold, "alpha_strong": args.alpha_strong,
        "velocities_kms": list(args.velocities), "dispersions_kms": list(args.dispersions),
        "average_powers": list(args.powers),
        "velocity_definition": "relative velocity v = 2 beta_cm; s = 4 M3^2/(1 - v^2/4)",
        "maxwellian": "f(v) ~ v^2 exp(-v^2/(4 sigma_1D^2)); average <sigma_V v^n>/<v^n>",
        "sigma_normalisation": "sigma = (1/2) int dOmega |M|^2/(64 pi^2 s) (identical particles)",
        "sigmaV_normalisation": "sigma_V = (3/2)(1/2) int dOmega sin^2 dsigma/dOmega; equals sigma if isotropic",
        "units": "sigma/m in cm^2/g; velocities in km/s; widths in GeV",
        "widths": "Gamma_i,other = w_i - Gamma(h_i -> h3 h3); H3 loop self-energy added on-shell subtracted",
    }


def _options(args):
    return {"mode": args.mode, "lmax": args.lmax, "n_angle": args.n_angle, "yukawa": args.yukawa,
            "yukawa_threshold": args.yukawa_threshold}


def run_scan(args):
    source = Path(args.input).resolve()
    destination = Path(args.output).resolve()
    if destination == source:
        raise SystemExit("--output must differ from the input scan")
    columns = sidm_columns(args.velocities, args.dispersions, args.powers)
    header, rows = read_scan(source)

    def evaluate(row):
        try:
            point = SIDMPoint.from_scan_row(row, M1=args.m1, v=args.vev)
            return analyse_point(point, velocities_kms=args.velocities, dispersions_kms=args.dispersions,
                                 powers=args.powers, alpha_strong=args.alpha_strong, **_options(args)), True
        except Exception as error:  # record and continue
            return failed_columns(error, columns, args.mode), False

    n_ok = write_scan(destination, header, rows, columns, evaluate, progress=args.progress)
    meta = configuration(args)
    meta.update({"input": str(source), "input_sha256": sha256(source), "output": str(destination),
                 "rows": len(rows), "rows_ok": n_ok, "columns": columns})
    meta_path = destination.with_name(destination.name + ".sidm.json")
    meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {destination} ({n_ok}/{len(rows)} rows ok) and {meta_path.name}")


def point_from_args(args):
    return SIDMPoint.from_potential(
        args.m2, args.m3, args.vs, args.a12, args.lx, args.lphix, args.lsx,
        w1=args.w1, w2=args.w2, gamma1_other=args.gamma1_other,
        gamma2_other=args.gamma2_other, M1=args.m1, v=args.vev)


def run_point(args):
    point = point_from_args(args)
    result = analyse_point(point, velocities_kms=args.velocities, dispersions_kms=args.dispersions,
                           powers=args.powers, alpha_strong=args.alpha_strong, **_options(args))
    if args.json:
        print(json.dumps({"point": asdict(point), "result": result}, indent=2, default=str))
        return
    print("point:", ", ".join(f"{k}={v:.6g}" for k, v in asdict(point).items()))
    for key, value in result.items():
        print(f"{key:44s} {format_value(value)}")


def run_curve(args):
    point = point_from_args(args)
    v = np.geomspace(args.vmin, args.vmax, args.n) / C_KMS
    extra = [point.resonance(i)[0] for i in (1, 2)]
    v = np.unique(np.concatenate([v, [e for e in extra if math.isfinite(e) and v[0] <= e <= v[-1]]]))
    sigma, sigma_v = cross_sections(point, v, **_options(args))
    print("v_kms\tsigma_over_m_cm2_g\tsigmaV_over_m_cm2_g")
    for vi, a, b in zip(v, per_mass_cm2_g(sigma, point.M3), per_mass_cm2_g(sigma_v, point.M3)):
        print(f"{vi * C_KMS:.10g}\t{a:.10g}\t{b:.10g}")


def add_physics_arguments(parser):
    parser.add_argument("--mode", choices=AMPLITUDE_MODES, default=DEFAULT_MODE)
    parser.add_argument("--n-angle", type=int, default=DEFAULT_N_ANGLE)
    parser.add_argument("--lmax", type=int, default=DEFAULT_LMAX)
    parser.add_argument("--m1", type=float, default=DEFAULT_M1)
    parser.add_argument("--vev", type=float, default=DEFAULT_SM_VEV,
                        help="v used only when K133/K233 must be rebuilt from lPhiX, lSX")
    parser.add_argument("--yukawa", choices=YUKAWA_CHOICES, default="auto",
                        help="Schroedinger phase shifts for t/u exchange (auto: above --yukawa-threshold)")
    parser.add_argument("--yukawa-threshold", type=float, default=DEFAULT_YUKAWA_THRESHOLD,
                        help="Born parameter alpha M3/M_i above which the solver is used")
    parser.add_argument("--alpha-strong", type=float, default=DEFAULT_ALPHA_STRONG,
                        help="alpha_i above which the point is flagged unreliable")


def _add_common(parser):
    add_physics_arguments(parser)
    parser.add_argument("--velocities", type=float, nargs="*", default=list(DEFAULT_VELOCITIES_KMS),
                        help="fixed relative velocities [km/s]")
    parser.add_argument("--dispersions", type=float, nargs="*", default=list(DEFAULT_DISPERSIONS_KMS),
                        help="halo 1D velocity dispersions [km/s] for Maxwell averages")
    parser.add_argument("--powers", type=int, nargs="*", default=list(DEFAULT_POWERS),
                        help="n in <sigma_V v^n>/<v^n>")


def add_point_arguments(parser):
    for name in ("m2", "m3", "vs", "a12", "lx", "lphix", "lsx"):
        parser.add_argument(f"--{name}", type=float, required=True)
    parser.add_argument("--w1", type=float, help="physical h1 total width [GeV]")
    parser.add_argument("--w2", type=float, help="physical h2 total width [GeV]")
    parser.add_argument("--gamma1-other", type=float, help="h1 width without h1->h3h3 [GeV]")
    parser.add_argument("--gamma2-other", type=float, help="h2 width without h2->h3h3 [GeV]")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    scan = sub.add_parser("scan", help="add SIDM columns to a scan TSV (written to a new file)")
    scan.add_argument("input")
    scan.add_argument("--output", required=True)
    scan.add_argument("--progress", type=int, default=0, help="report every N rows")
    _add_common(scan)
    point = sub.add_parser("point", help="evaluate one point")
    add_point_arguments(point)
    _add_common(point)
    point.add_argument("--json", action="store_true")
    curve = sub.add_parser("curve", help="print sigma/m and sigma_V/m versus v for one point")
    add_point_arguments(curve)
    _add_common(curve)
    curve.add_argument("--vmin", type=float, default=1.0, help="km/s")
    curve.add_argument("--vmax", type=float, default=1.0e4, help="km/s")
    curve.add_argument("--n", type=int, default=400)
    args = parser.parse_args(argv)
    {"scan": run_scan, "point": run_point, "curve": run_curve}[args.command](args)


if __name__ == "__main__":
    main()
