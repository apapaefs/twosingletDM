#!/usr/bin/env python3
"""Evaluate the audit's numerical counterexamples under the v2 convention."""
import sys,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from trsm_theory_diagnostics import assess_vacuum,assess_running
from test_trsm_theory_constraints import _unitarity_eigenvalues
from test_trsm_DM import direct_detection_base_limit,test_dm
from trsm_inputs import M1,MW,MZ,sm_inputs,json_safe
from ewpt_entry_criterion import ew_entry_updates
from trsm_constraint_profile import profile_updates
import singlet_EWPO as ew
import test_trsm_higgstools as ht
old=json.loads((ROOT/'benchmarks/v2/audit-baseline.json').read_text())
v=assess_vacuum(300,200,50,.1,.1,.1,1);r=assess_running(300,200,50,0,3,0,0)
assert v['vacuum_tree_global'] is False and r['rg_unitarity'] is False and r['rg_integration_success'] is True
report={'tree_vacuum':{'old':old['deeper_x_vacuum'],'new':v},'running':{'old':old['finite_but_nonperturbative_running'],'new':r},
        'unitarity':[{'old':case,'new_eigenvalues':_unitarity_eigenvalues(*case['quartics']).tolist()} for case in old['degenerate_unitarity_roots']],
        'EWPO':[{'M2':mass,'new_STU':[fn(M1,mass,.2,MZ,MW) for fn in (ew.Delta_S,ew.Delta_T,ew.Delta_U)]} for mass in (MW,MZ,MZ*(1+1e-8),MZ*(1-1e-8))],
        'DD_boundaries':[{'old':case,'new_below_pb':direct_detection_base_limit(case['mass']*(1-1e-10)),'new_at_pb':direct_detection_base_limit(case['mass'])} for case in old['dd_piecewise_boundaries']],
        'empty_EWPT':{'old':old['empty_ewpt_payload_flags'],'new':{k:ew_entry_updates({})[k] for k in ('ewpt_baryo_candidate','ewpt_gw_candidate')}},
        'HiggsSignals_SM':{'old':old['higgssignals_decoupling_limit'],'new_chi2_reference':ht.ress_SM,'new_delta_chi2':ht._sm_reference_chi2()-ht.ress_SM},
        'SM_inputs':{'old':old['SM_input_conventions'],'new':sm_inputs()},
        'W_mass_coverage':profile_updates(dict(M2=20,wmass=True,hb=True,hs=True,ewpo=True,thc=True,dm=False))}
assert all(math.isfinite(x) for p in report['EWPO'] for x in p['new_STU'])
assert report['HiggsSignals_SM']['new_delta_chi2']==0
for p in report['DD_boundaries']:assert math.isclose(p['new_below_pb'],p['new_at_pb'],rel_tol=1e-8)
report['status']='passed'
out=Path(sys.argv[1]);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(json_safe(report),indent=2,allow_nan=False)+'\n');print(out)
