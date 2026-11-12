"""Read-only numerical probes for the constraint audit at b7b0d93.

Run from the repository root with its SciPy/HiggsTools Python environment:
    python audits/2026-09-23/reproduce_constraints.py
No micrOMEGAs/BSMPT jobs are launched and no scan files are modified.
"""

import contextlib
import io
import json
import math
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path.cwd()))
import numpy as np
from scipy.integrate import solve_ivp
import test_trsm_evolution as evolution
import test_trsm_theory_constraints as theory
import singlet_EWPO as ewpo
from test_trsm_DM import lz2025_source_direct_detection_base_limit
from ewpt_entry_criterion import ew_entry_updates


def quartics(vs, m2, angle, lx, lphix, lsx):
    c, s = math.cos(angle), math.sin(angle)
    return np.array([
        (125.09**2 * c*c + m2*m2 * s*s) / (2 * 246**2),
        (125.09**2 * s*s + m2*m2 * c*c) / (2 * vs*vs),
        lx, (m2*m2 - 125.09**2) * s*c / (246 * vs), lphix, lsx,
    ])


def initial_conditions(q, vs, m3):
    lh, ls, lx, lhs, lhx, lsx = q
    # Reproduce the code's current mass initialization (including its missing m3^2).
    mu = [-lh*246**2 - lhs*vs*vs/2,
          -ls*vs*vs - lhs*246**2/2,
          -lhx*246**2/2 - lsx*vs*vs/2]
    return [1.221, math.sqrt(.424), math.sqrt(.1273), .96738, *q, *mu]


report = {'SM_input_conventions': {
    'scan_higgs_mass_GeV': 125.09, 'micromegas_higgs_mass_GeV': 125.,
    'scan_v_GeV': 246., 'bsmpt_v_GeV': math.sqrt(1/math.sqrt(2)/1.1663787e-5),
    'micromegas_v_GeV': 2*80.385*.474/.31333,
    'micromegas_derived_GF': .31333**2/(2*.474*80.385)**2/math.sqrt(2),
}}
params = dict(vs=300., M2=200., M3=50., a12=.1, lX=.1, lPhiX=.1, lSX=1.)
q = quartics(params['vs'], params['M2'], params['a12'], params['lX'], params['lPhiX'], params['lSX'])
lh, ls, lx, lhs, lhx, lsx = q
A = np.array([[lh, lhs/2, lhx/2], [lhs/2, ls, lsx/2], [lhx/2, lsx/2, lx]])
mu = np.array([-lh*246**2-lhs*300**2/2, -ls*300**2-lhs*246**2/2,
               50**2-lhx*246**2/2-lsx*300**2/2])
y_wanted = np.array([246.**2, 300.**2, 0.])
y_x = np.array([0., 0., -mu[2]/lx])
def stationary_record(y):
    fields = np.sqrt(y)
    hessian = np.diag(mu + A @ y) + 2 * np.outer(fields, fields) * A
    return {'fields_h_s_x_GeV': fields.tolist(),
            'potential_GeV4': float(mu @ y / 2 + y @ A @ y / 4),
            'radial_hessian_eigenvalues_GeV2': np.linalg.eigvalsh(hessian).tolist()}
report['deeper_x_vacuum'] = {'input': params, 'thc': theory.theory_constraints_vxzero(**params),
    'evo': bool(evolution.solv_eq(91, 1000, .0001, initial_conditions(q, 300, 50), 900)),
    'wanted_vacuum': stationary_record(y_wanted), 'x_vacuum': stationary_record(y_x)}

q_run = quartics(300, 200, 0, 3., 0, 0)
ini = initial_conditions(q_run, 300, 50)
run = solve_ivp(evolution.rhs, (91, 1000), ini, t_eval=[900, 1000], rtol=1e-9, atol=1e-10)
report['finite_but_nonperturbative_running'] = {
    'input': {'vs':300, 'M2':200, 'M3':50, 'a12':0, 'lX':3, 'lPhiX':0, 'lSX':0},
    'thc_at_input': theory._passes_tree_level_constraints(*q_run),
    'evo_flag': bool(evolution.solv_eq(91, 1000, .0001, ini, 900)),
    'solver_success': bool(run.success), 'lambdaX_900_1000_GeV': run.y[6].tolist(),
    'thc_at_900_1000_GeV': [theory._passes_tree_level_constraints(*row) for row in run.y[4:10].T],
}

report['degenerate_unitarity_roots'] = []
for q_degenerate in [(0.5, 1, 1, 0, 0, 0), (125.09**2/(2*246**2), 1, 1, 0, 0, 0)]:
    roots = theory._unitarity_eigenvalues(*q_degenerate)
    report['degenerate_unitarity_roots'].append({'quartics':q_degenerate,
        'roots_real_imag':[[float(z.real), float(z.imag)] for z in roots],
        'thc':theory._passes_tree_level_constraints(*q_degenerate)})

report['ewpo_equal_mass_probes'] = []
for m2 in (80.385, 91.1876, 91.1876*(1+1e-8), 91.1876*(1-1e-8)):
    try:
        values = [ewpo.Delta_S(125.09, m2, .2, 91.1876, 80.385),
                  ewpo.Delta_T(125.09, m2, .2, 91.1876, 80.385),
                  ewpo.Delta_U(125.09, m2, .2, 91.1876, 80.385)]
        report['ewpo_equal_mass_probes'].append({'M2':m2,'STU':values})
    except Exception as error:
        report['ewpo_equal_mass_probes'].append({'M2':m2,'error':f'{type(error).__name__}: {error}'})
with warnings.catch_warnings(record=True) as recorded:
    report['wmass_out_of_range'] = {'pass_at_M2_20_sin_05': ewpo.check_wmass_tania(20, .5),
        'warnings':[str(w.message) for w in recorded]}
report['dd_piecewise_boundaries'] = [{'mass':m, 'below_pb':lz2025_source_direct_detection_base_limit(m-1e-8),
    'at_pb':lz2025_source_direct_detection_base_limit(m)} for m in (10,40,60,80,100)]
report['empty_ewpt_payload_flags'] = {k:v for k,v in ew_entry_updates({
    'calctemps':{'status_tracing':'no_coverage', 'status_coex_pairs':'not_set'},
    'transition_strengths':[]}).items() if k.endswith('_candidate')}

try:
    with contextlib.redirect_stdout(io.StringIO()):
        import test_trsm_higgstools as ht
        from generate_trsm_info import generate_lams
        info = generate_lams(1, 200., 300., 300., 0., 0., 0., 0., False, lX=.1, lPhiX=0., lSX=0.)
        hb, hs, details = ht.analyze_parampoint(ht.pred, ht.H1, ht.H2, ht.H3,
            125.09, 200., 300., 1., 0., 0., info[22], info[23], info[24], return_details=True)
    report['higgssignals_decoupling_limit'] = {'hb':hb, 'hs':hs,
        **{k:v for k,v in details['higgssignals'].items() if k!='top_chi2_contributors'}}
except Exception as error:
    report['higgssignals_decoupling_limit'] = {'error':f'{type(error).__name__}: {error}'}

print(json.dumps(report, indent=2, allow_nan=False))
