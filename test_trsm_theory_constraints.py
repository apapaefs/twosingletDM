"""Tree-level boundedness and perturbative-unitarity checks for the TRSM."""

import math

import numpy as np


__all__ = ["theory_constraints", "theory_constraints_vxzero"]


VH = 246.0
M1 = 125.09
COPOSITIVITY_TOLERANCE = 1.0e-12
ROOT_IMAGINARY_TOLERANCE = 1.0e-10


def _all_finite(*values):
    try:
        return all(math.isfinite(float(value)) for value in values)
    except (TypeError, ValueError, OverflowError):
        return False


def _rotation_matrix(a12, a13, a23):
    c1, c2, c3 = math.cos(a12), math.cos(a13), math.cos(a23)
    s1, s2, s3 = math.sin(a12), math.sin(a13), math.sin(a23)
    return (
        (c1 * c2, -s1 * c2, -s2),
        (s1 * c3 - c1 * s2 * s3, c1 * c3 + s1 * s2 * s3, -c2 * s3),
        (c1 * s2 * c3 + s1 * s3, c1 * s3 - s1 * s2 * c3, c2 * c3),
    )


def _copositive_quartic_matrix(lphi, ls, lx, lphis, lphix, lsx):
    """Check boundedness of the quartic potential in all field directions.

    For V4 = (lambda_i phi_i^4 + lambda_ij phi_i^2 phi_j^2) / 4,
    the matrix acting on the non-negative variables phi_i^2 has diagonal
    entries lambda_i and off-diagonal entries lambda_ij / 2.
    """

    diagonal = (lphi, ls, lx)
    if any(value < -COPOSITIVITY_TOLERANCE for value in diagonal):
        return False
    a11, a22, a33 = (max(0.0, value) for value in diagonal)
    a12, a13, a23 = lphis / 2.0, lphix / 2.0, lsx / 2.0

    bar12 = a12 + math.sqrt(a11 * a22)
    bar13 = a13 + math.sqrt(a11 * a33)
    bar23 = a23 + math.sqrt(a22 * a33)
    if any(
        value < -COPOSITIVITY_TOLERANCE
        for value in (bar12, bar13, bar23)
    ):
        return False

    final_condition = (
        math.sqrt(a11 * a22 * a33)
        + a12 * math.sqrt(a33)
        + a13 * math.sqrt(a22)
        + a23 * math.sqrt(a11)
        + math.sqrt(2.0 * max(0.0, bar12 * bar13 * bar23))
    )
    return final_condition >= -COPOSITIVITY_TOLERANCE


def _unitarity_eigenvalues(lphi, ls, lx, lphis, lphix, lsx):
    coefficients = [
        1.0,
        -12.0 * lphi - 6.0 * ls - 6.0 * lx,
        72.0 * lphi * (ls + lx)
        - 4.0 * (lphis**2 + lphix**2)
        + 36.0 * ls * lx
        - lsx**2,
        12.0 * lphi * lsx**2
        + 24.0 * lphis**2 * lx
        + 24.0 * lphix**2 * ls
        - 8.0 * lphis * lphix * lsx
        - 432.0 * lphi * ls * lx,
    ]
    if not _all_finite(*coefficients):
        return None
    roots = np.roots(coefficients)
    if len(roots) != 3:
        return None
    return roots


def _passes_tree_level_constraints(lphi, ls, lx, lphis, lphix, lsx):
    quartics = (lphi, ls, lx, lphis, lphix, lsx)
    if not _all_finite(*quartics):
        return False
    if not _copositive_quartic_matrix(*quartics):
        return False

    if any(abs(value) >= 4.0 * math.pi for value in (lphi, ls, lx)):
        return False
    if any(abs(value) >= 8.0 * math.pi for value in (lphis, lphix, lsx)):
        return False

    roots = _unitarity_eigenvalues(*quartics)
    if roots is None:
        return False
    for root in roots:
        real = float(np.real(root))
        imaginary = float(np.imag(root))
        if not _all_finite(real, imaginary):
            return False
        if abs(imaginary) > ROOT_IMAGINARY_TOLERANCE * (1.0 + abs(real)):
            return False
        if abs(real) >= 16.0 * math.pi:
            return False
    return True


def theory_constraints(vs, vx, M2, M3, a12, a13, a23):
    """Check the general two-singlet model at tree level."""

    inputs = (vs, vx, M2, M3, a12, a13, a23)
    if not _all_finite(*inputs) or vs == 0.0 or vx == 0.0:
        return False

    try:
        rotation = _rotation_matrix(a12, a13, a23)
        masses = (M1, M2, M3)
        lphi = sum(mass**2 * rotation[i][0] ** 2 for i, mass in enumerate(masses)) / (
            2.0 * VH**2
        )
        ls = sum(mass**2 * rotation[i][1] ** 2 for i, mass in enumerate(masses)) / (
            2.0 * vs**2
        )
        lx = sum(mass**2 * rotation[i][2] ** 2 for i, mass in enumerate(masses)) / (
            2.0 * vx**2
        )
        lphis = sum(
            mass**2 * rotation[i][0] * rotation[i][1]
            for i, mass in enumerate(masses)
        ) / (VH * vs)
        lphix = sum(
            mass**2 * rotation[i][0] * rotation[i][2]
            for i, mass in enumerate(masses)
        ) / (VH * vx)
        lsx = sum(
            mass**2 * rotation[i][1] * rotation[i][2]
            for i, mass in enumerate(masses)
        ) / (vs * vx)
    except (OverflowError, TypeError, ValueError, ZeroDivisionError):
        return False

    return _passes_tree_level_constraints(lphi, ls, lx, lphis, lphix, lsx)


def theory_constraints_vxzero(vs, M2, M3, a12, lX, lPhiX, lSX):
    """Check the Z2-symmetric vx=0 dark-matter branch at tree level."""

    inputs = (vs, M2, M3, a12, lX, lPhiX, lSX)
    if not _all_finite(*inputs) or vs == 0.0:
        return False

    try:
        cosine = math.cos(a12)
        sine = math.sin(a12)
        lphi = (M1**2 * cosine**2 + M2**2 * sine**2) / (2.0 * VH**2)
        ls = (M1**2 * sine**2 + M2**2 * cosine**2) / (2.0 * vs**2)
        lphis = (M2**2 - M1**2) * sine * cosine / (VH * vs)
    except (OverflowError, TypeError, ValueError, ZeroDivisionError):
        return False

    return _passes_tree_level_constraints(lphi, ls, lX, lphis, lPhiX, lSX)
