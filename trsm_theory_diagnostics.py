"""Vacuum and running-coupling screens for the canonical vx=0 potential.

V = sum(mu_i^2 phi_i^2)/2 + sum(lambda_i phi_i^4)/4
    + sum(lambda_ij phi_i^2 phi_j^2)/4, in (h,s,x) order.
These screens never determine whether an exploratory point is retained.
"""

import itertools
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import bisect, linprog
from trsm_inputs import M1, VEV, RG_BOUNDARY, nullable_and
from test_trsm_theory_constraints import _copositive_quartic_matrix, _unitarity_eigenvalues
from test_trsm_evolution import rhs

THEORY_COLUMNS = (
    "vacuum_tree_global", "vacuum_tree_status", "vacuum_tree_depth_gap_GeV4",
    "vacuum_tree_competitor_h_GeV", "vacuum_tree_competitor_s_GeV", "vacuum_tree_competitor_x_GeV",
    "rg_integration_success", "rg_bfb", "rg_unitarity", "rg_status",
    "rg_first_bfb_failure_GeV", "rg_first_unitarity_failure_GeV",
    "rg_reached_scale_GeV", "theory_strict_subset",
)


def potential_parameters(vs, m2, m3, angle, lx, lphix, lsx):
    args = (vs, m2, m3, angle, lx, lphix, lsx)
    if not all(math.isfinite(float(x)) for x in args) or vs == 0 or min(m2, m3) <= 0:
        raise ValueError("finite inputs, nonzero vs and positive scalar masses required")
    c, s = math.cos(angle), math.sin(angle)
    lh = (M1*M1*c*c + m2*m2*s*s)/(2*VEV**2)
    ls = (M1*M1*s*s + m2*m2*c*c)/(2*vs*vs)
    lhs = (m2*m2-M1*M1)*s*c/(VEV*vs)
    matrix = np.array([[lh, lhs/2, lphix/2], [lhs/2, ls, lsx/2], [lphix/2, lsx/2, lx]])
    mu = np.array([-lh*VEV**2-lhs*vs**2/2, -ls*vs**2-lhs*VEV**2/2,
                   m3*m3-lphix*VEV**2/2-lsx*vs**2/2])
    return np.array([lh, ls, lx, lhs, lphix, lsx]), matrix, mu


def assess_vacuum(vs, m2, m3, angle, lx, lphix, lsx):
    quartics, matrix, mu = potential_parameters(vs, m2, m3, angle, lx, lphix, lsx)
    target = np.array([VEV**2, vs**2, 0.0])
    def depth(y):
        return float(mu @ y/2 + y @ matrix @ y/4)
    def hessian(y):
        fields = np.sqrt(np.maximum(y, 0))
        return np.diag(mu + matrix @ y) + 2*matrix*np.outer(fields, fields)
    target_depth = depth(target)
    result = {"vacuum_tree_global": None, "vacuum_tree_status": "unresolved",
              "vacuum_tree_depth_gap_GeV4": None, "stationary_points": []}
    if not _copositive_quartic_matrix(*quartics):
        result.update(vacuum_tree_global=False, vacuum_tree_status="unbounded_quartic")
        return result
    unresolved = False
    unbounded_flat = False
    minima = []
    for size in range(4):
        for support in itertools.combinations(range(3), size):
            y = np.zeros(3)
            if support:
                ids = list(support)
                sub = matrix[np.ix_(ids, ids)]
                b = -mu[ids]
                if np.linalg.matrix_rank(sub, tol=1e-12*max(1.0, np.linalg.norm(sub))) < size:
                    # A nonnegative quartic-null direction with negative mass
                    # term is unbounded even though quartic copositivity passes.
                    lp = linprog(mu[ids], A_eq=np.vstack([sub, np.ones(size)]),
                                 b_eq=np.r_[np.zeros(size), 1.0], bounds=(0, None), method="highs")
                    if lp.success and lp.fun < -1e-9*max(1.0, np.linalg.norm(mu)):
                        unbounded_flat = True
                    feasible = linprog(np.zeros(size), A_eq=sub, b_eq=b,
                                       bounds=(0, None), method="highs")
                    if not feasible.success:
                        continue
                    unresolved = True
                    solution = feasible.x
                else:
                    solution = np.linalg.solve(sub, b)
                tolerance = 1e-10*max(1.0, np.max(np.abs(solution)))
                if np.any(solution < -tolerance):
                    continue
                y[ids] = np.maximum(solution, 0)
            eigen = np.linalg.eigvalsh(hessian(y))
            local = bool(np.min(eigen) >= -1e-8*max(1.0, np.max(np.abs(eigen))))
            record = {"fields_GeV": np.sqrt(y).tolist(), "potential_GeV4": depth(y),
                      "hessian_eigenvalues_GeV2": eigen.tolist(), "local_minimum": local}
            result["stationary_points"].append(record)
            if local:
                minima.append(record)
    if unbounded_flat:
        result.update(vacuum_tree_global=False, vacuum_tree_status="unbounded_flat_direction")
        return result
    if not minima:
        return result
    lowest = min(minima, key=lambda p: p["potential_GeV4"])
    gap = target_depth-lowest["potential_GeV4"]
    tolerance = 1e-8*max(1.0, abs(target_depth), abs(lowest["potential_GeV4"]))
    result["vacuum_tree_depth_gap_GeV4"] = gap
    for name, value in zip(("h", "s", "x"), lowest["fields_GeV"]):
        result[f"vacuum_tree_competitor_{name}_GeV"] = value
    target_local = np.min(np.linalg.eigvalsh(hessian(target))) >= -1e-7
    if not target_local or gap > tolerance:
        result.update(vacuum_tree_global=False, vacuum_tree_status="nonglobal" if target_local else "not_local_minimum")
    elif unresolved:
        result["vacuum_tree_status"] = "flat_stationary_family"
    else:
        degenerate = any(abs(p["potential_GeV4"]-target_depth) <= tolerance and
                         not np.allclose(np.square(p["fields_GeV"]), target, rtol=1e-6, atol=1e-6)
                         for p in minima)
        result.update(vacuum_tree_global=True, vacuum_tree_status="degenerate_global" if degenerate else "global")
    return result


def unitarity_margin(q):
    roots = _unitarity_eigenvalues(*q)
    if roots is None:
        return -math.inf
    return min(1-max(abs(q[:3]))/(4*math.pi), 1-max(abs(q[3:]))/(8*math.pi),
               1-max(abs(roots))/(16*math.pi))


def assess_running(vs, m2, m3, angle, lx, lphix, lsx, *, max_scale=1000.0):
    q, _, mu = potential_parameters(vs, m2, m3, angle, lx, lphix, lsx)
    lo = RG_BOUNDARY["scale_GeV"]
    if not math.isfinite(max_scale) or max_scale <= lo:
        raise ValueError("RG maximum scale must exceed the boundary scale")
    initial = np.r_[[RG_BOUNDARY[k] for k in ("g3", "g2", "g1", "yt")], q, mu]
    result = {"rg_integration_success": False, "rg_bfb": None, "rg_unitarity": None,
              "rg_status": "integration_failed", "rg_reached_scale_GeV": lo,
              "rg_first_bfb_failure_GeV": None, "rg_first_unitarity_failure_GeV": None}
    try:
        sol = solve_ivp(rhs, (lo, max_scale), initial, rtol=1e-8, atol=1e-10,
                        max_step=(max_scale-lo)/256, dense_output=True)
    except (ValueError, FloatingPointError, OverflowError, RuntimeError):
        return result
    result["rg_reached_scale_GeV"] = float(sol.t[-1])
    success = bool(sol.success and np.all(np.isfinite(sol.y)))
    result["rg_integration_success"] = success
    grid = np.unique(np.r_[sol.t, np.geomspace(lo, sol.t[-1], 513)])
    for name, predicate in (("bfb", lambda q: _copositive_quartic_matrix(*q)),
                            ("unitarity", lambda q: unitarity_margin(q) > 0)):
        previous = lo
        failure = None
        for scale in grid:
            if not predicate(sol.sol(scale)[4:10]):
                if scale == lo:
                    failure = lo
                else:
                    left, right = previous, float(scale)
                    for _ in range(50):
                        if right-left < 1e-5:
                            break
                        mid = (left+right)/2
                        if predicate(sol.sol(mid)[4:10]): left = mid
                        else: right = mid
                    failure = right
                break
            previous = float(scale)
        result[f"rg_{name}"] = False if failure is not None else (True if success else None)
        result[f"rg_first_{name}_failure_GeV"] = failure
    result["rg_status"] = "assessed" if success else "integration_failed"
    return result


def theory_diagnostics(vs, m2, m3, angle, lx, lphix, lsx):
    args = (vs, m2, m3, angle, lx, lphix, lsx)
    vacuum = assess_vacuum(*args)
    result = {key: value for key, value in vacuum.items() if key != "stationary_points"}
    result.update(assess_running(*args))
    result["theory_strict_subset"] = nullable_and(result[k] for k in
        ("vacuum_tree_global", "rg_integration_success", "rg_bfb", "rg_unitarity"))
    return result
