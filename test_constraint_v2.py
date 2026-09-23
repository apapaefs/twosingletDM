"""Scientific and status regressions for the active vx=0 release."""
import csv,json,math,tempfile,unittest
from pathlib import Path
import numpy as np
from trsm_inputs import M1,VEV,GF,MW,MZ,EE,SW,PHYSICS_VERSION
from trsm_theory_diagnostics import assess_vacuum,assess_running,potential_parameters
from test_trsm_theory_constraints import _unitarity_eigenvalues
from test_trsm_DM import DMPoint,write_micromegas_card,test_dm,direct_detection_base_limit
from trsm_direct_detection import load_si_limit_table,DEFAULT_LIMIT_TABLE
from ewpt_entry_criterion import ew_entry_updates
from test_trsm_ewpt import calculate_fopt_strengths
from dm_thermal_relic_diagnostic import coupling_sensitivity_updates,resonance_proximity_updates
from ewpt_x_history import x_history_updates
from trsm_constraint_profile import profile_updates
from test_ewpt_entry_criterion import transition
from scan_output import output_columns,output_row

ROOT=Path(__file__).resolve().parent

def assessed(*strengths):
    row={'status_nlo_stability':'success','status_tracing':'success','status_coex_pairs':'success'}
    for s in strengths:
        index=s['transition_index'];row[f'status_{s["temperature_kind"]}_{index}']='success'
        for k in ('crit','nucl','perc','compl'):row.setdefault(f'status_{k}_{index}','not_met')
    return {'calctemps':row,'transition_strengths':list(strengths)}

class ConstraintV2Tests(unittest.TestCase):
    def test_shared_SM_and_exact_card_pole(self):
        self.assertAlmostEqual(VEV,(math.sqrt(2)*GF)**-.5)
        self.assertAlmostEqual(SW**2,1-MW**2/MZ**2)
        self.assertAlmostEqual(EE,2*MW*SW/VEV)
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'card';write_micromegas_card(DMPoint(.1,.02,.03,M1/2,300,.1,200),path)
            p={k:float(v) for k,v in (line.split() for line in path.read_text().splitlines())}
            self.assertEqual(2*p['MX'],p['Mh']);self.assertEqual(p['EE'],EE)
    def test_repeated_eigenvalues(self):
        eig=_unitarity_eigenvalues(M1*M1/(2*VEV*VEV),1,1,0,0,0)
        np.testing.assert_allclose(eig,[12*M1*M1/(2*VEV*VEV),6,6],atol=1e-12)
        self.assertTrue(np.isrealobj(eig))
    def test_deeper_X_vacuum(self):
        result=assess_vacuum(300,200,50,.1,.1,.1,1)
        self.assertIs(result['vacuum_tree_global'],False)
        self.assertGreater(result['vacuum_tree_depth_gap_GeV4'],4.5e9)
        self.assertAlmostEqual(result['vacuum_tree_competitor_x_GeV'],674.76815,places=3)
        self.assertGreaterEqual(len(result['stationary_points']),4)
    def test_degenerate_global_and_flat_branches(self):
        lx=(50**2-300**2/2)**2/((200**2/(2*300**2))*300**4)
        result=assess_vacuum(300,200,50,0,lx,0,1)
        self.assertEqual(result['vacuum_tree_status'],'degenerate_global')
        self.assertIs(result['vacuum_tree_global'],True)
        flat=assess_vacuum(300,200,VEV/math.sqrt(2),0,200**2/(2*300**2),1,200**2/300**2)
        self.assertEqual(flat['vacuum_tree_status'],'flat_stationary_family')
        self.assertIsNone(flat['vacuum_tree_global'])
        self.assertEqual(assess_vacuum(300,200,50,0,0,.1,1)['vacuum_tree_status'],'unbounded_flat_direction')
    def test_RG_counterexample_and_mass_initialization(self):
        result=assess_running(300,200,50,0,3,0,0)
        self.assertIs(result['rg_integration_success'],True)
        self.assertIs(result['rg_bfb'],True);self.assertIs(result['rg_unitarity'],False)
        self.assertTrue(594.65<result['rg_first_unitarity_failure_GeV']<594.67)
        self.assertEqual(result['rg_reached_scale_GeV'],1000)
        self.assertEqual(potential_parameters(300,200,50,0,3,0,0)[2][2],2500)
    def test_DD_units_boundaries_continuity(self):
        table=load_si_limit_table(DEFAULT_LIMIT_TABLE)
        self.assertAlmostEqual(table.upper_limit_pb(40)/2.1816833824484827e-12,1)
        self.assertTrue(math.isnan(direct_detection_base_limit(8.9999)))
        self.assertTrue(math.isfinite(direct_detection_base_limit(9)))
        self.assertTrue(math.isfinite(direct_detection_base_limit(10000)))
        self.assertTrue(math.isnan(direct_detection_base_limit(10000.01)))
        self.assertAlmostEqual(direct_detection_base_limit(10-1e-7)/direct_detection_base_limit(10+1e-7),1,places=6)
    def test_relic_boundary_solver_error_and_abundance(self):
        template='MX = 40\nXf=20 Omega={omega:.17g}\ndarkOmega_error={error}\n~X[~X]-nucleon cross sections[pb]:\nneutron SI 1e-30\n'
        for omega,excluded in [(.121-1e-12,False),(.121,False),(.121+1e-12,True)]:
            passed,_,r=test_dm(.1,.1,.1,40,300,.1,200,raw_output=template.format(omega=omega,error=0))
            self.assertIs(r['dm_relic_excluded'],excluded);self.assertIs(passed,not excluded)
        passed,_,r=test_dm(.1,.1,.1,40,300,.1,200,raw_output=template.format(omega=.06,error=1))
        self.assertIsNone(passed);self.assertEqual(r['dm_solver_error'],1);self.assertIsNone(r['dm_relic_excluded'])
        _,_,r=test_dm(.1,.1,.1,40,300,.1,200,raw_output=template.format(omega=.06,error=0))
        self.assertAlmostEqual(r['dm_dir_det_limit']/r['dm_lux_base_limit'],2)
    def test_EWPO_removable_limits(self):
        from singlet_EWPO import H_S,H_T
        for function,x in [(H_S,1),(lambda x:H_T(x,MW/MZ),1),(lambda x:H_T(x,MW/MZ),(MW/MZ)**2),(H_S,4)]:
            center=function(x);self.assertTrue(math.isfinite(center))
            self.assertAlmostEqual(function(x*(1+1e-7)),center,places=5)
            self.assertAlmostEqual(function(x*(1-1e-7)),center,places=5)
        self.assertEqual(H_S(0),0);self.assertEqual(H_T(0,MW/MZ),0)
        from singlet_EWPO import check_wmass_tania
        self.assertIsNone(check_wmass_tania(float('nan'),.1))
    def test_failed_solver_preserves_evidence(self):
        text='TRSM_inputs_v2 {"MX":40,"width_h1":0.004,"width_h2":0.02}\nMX = 40\nXf=nan Omega=nan\ndarkOmega_error=4\n'
        passed,_,r=test_dm(.1,.1,.1,40,300,.1,200,raw_output=text)
        self.assertIsNone(passed);self.assertEqual(r['dm_solver_error'],4)
        self.assertEqual(r['dm_h1_width_GeV'],.004);self.assertIsNone(r['dm_xf'])
        self.assertEqual(json.loads(r['dm_actual_inputs'])['raw_Omega'],'nan')
    def test_strict_JSON_null_and_partial_executable_status(self):
        from trsm_inputs import json_safe
        from test_trsm_ewpt import run_trsm_ewpt,EWPTConfig,TRSMEWPTPoint,result_to_json
        self.assertEqual(json.dumps(json_safe({'x':float('nan')}),allow_nan=False),'{"x": null}')
        with tempfile.TemporaryDirectory() as directory:
            result=run_trsm_ewpt(TRSMEWPTPoint(200,50,300,.1,.1,.01,.02),
                config=EWPTConfig(executable=Path(directory)/'missing',minima_executable=Path(directory)/'also-missing'),workdir=directory)
            payload=result_to_json(result)
            self.assertEqual(payload['execution']['calctemps_returncode'],127)
            self.assertEqual(payload['execution']['minimatracer_returncode'],127)
            self.assertIsNone(ew_entry_updates(payload)['ewpt_gw_candidate'])
    def test_reevaluation_preserves_unassessed_Higgs_and_DM_diagnostics(self):
        from unittest.mock import patch
        from reevaluate_trsm_dm_higgs import CoreEvaluator
        raw='MX = 40\nXf=20 Omega=0.06\ndarkOmega_error=0\n~X[~X]-nucleon cross sections[pb]:\nneutron SI 1e-30\n'
        dm=test_dm(.1,.01,.02,40,300,.1,200,raw_output=raw)
        with tempfile.TemporaryDirectory() as directory:
            executable=Path(directory)/'main';executable.touch()
            evaluator=CoreEvaluator(executable)
            with patch.object(evaluator.higgs,'analyze_parampoint',side_effect=RuntimeError('fixture Higgs numerical failure')),patch.object(evaluator.dm_provider,'test_dm',return_value=dm):
                result=evaluator(dict(M2='200',M3='40',vs='300',vx='0',a12='.1',lX='.1',lPhiX='.01',lSX='.02'),1)
        self.assertIsNone(result['experimental_subset']);self.assertIsNone(result['hb']);self.assertIsNone(result['hs'])
        self.assertIs(result['dm'],True);self.assertIs(result['ewpt_eligible'],True)
        self.assertEqual(result['point_assessment_status'],'partially_assessed')
        self.assertIn('fixture Higgs numerical failure',result['point_assessment_reason'])
    def test_nullable_candidates_and_strict_cut(self):
        for jump,expected in [(100,False),(100.000001,True)]:
            payload=assessed(transition(0,'crit',100,(0,0,0),(jump,0,0)))
            r=ew_entry_updates(payload);self.assertIs(r['ewpt_baryo_candidate'],expected);self.assertIs(r['ewpt_gw_candidate'],expected)
        payload=assessed(transition(0,'crit',100,(0,0,0),(60,0,0)),transition(0,'nucl',40,(0,0,0),(70,0,0)))
        r=ew_entry_updates(payload);self.assertIs(r['ewpt_baryo_candidate'],False);self.assertIs(r['ewpt_gw_candidate'],True)
        r=ew_entry_updates(assessed(transition(0,'crit',100,(0,0,0),(0,0,110))))
        self.assertIs(r['ewpt_baryo_candidate'],False);self.assertIs(r['ewpt_gw_candidate'],True)
        payload['calctemps']['status_tracing']='failed';r=ew_entry_updates(payload)
        self.assertIsNone(r['ewpt_baryo_candidate']);self.assertIs(r['ewpt_gw_candidate'],True)
        self.assertIsNone(ew_entry_updates({})['ewpt_gw_candidate'])
    def test_symmetry_copies_do_not_create_a_GW_jump(self):
        r=ew_entry_updates(assessed(transition(0,'crit',100,(240,0,-300),(-240,0,300))))
        self.assertIs(r['ewpt_gw_candidate'],False)
    def test_point223615_classification(self):
        p=json.loads((ROOT/'benchmarks/v2/point_223615/legacy-calctemps.json').read_text())
        p['transition_strengths']=[x.to_dict() for x in calculate_fopt_strengths(p['calctemps'])]
        r=ew_entry_updates(p);self.assertIs(r['ewpt_baryo_candidate'],False);self.assertIs(r['ewpt_gw_candidate'],True)
        self.assertAlmostEqual(r['ewpt_gw_crit_field_jump_over_T'],6.5908,places=4)
    def test_unknown_equilibrium_never_passes_Z2(self):
        payload={'minimatracer':{'global_branch':[{'temp':0,'label':'EW'},{'temp':2,'label':'UNRESOLVED'},{'temp':4,'label':'EW'}]}}
        self.assertIsNone(x_history_updates(payload,3)['dm_relic_z2_freezeout_compatible'])
    def test_coupling_cancellation_and_width_gap(self):
        payload={'minimatracer':{'global_branch':[{'temp':t,'w1':VEV-t,'wx':0,'ws':300,'phase_index':0} for t in (0,1,2,4,5)]}}
        point={'dm_freezeout_temperature_GeV':2,'lPhiX':1,'lSX':VEV/300,'a12':math.pi/4,'vs':300}
        r=coupling_sensitivity_updates(payload,point)
        self.assertIsNone(r['dm_thermal_K133_max_fractional_change'])
        self.assertGreater(r['dm_thermal_K133_max_abs_change_GeV'],0)
        self.assertIs(r['dm_thermal_K133_cancellation_sensitive'],True)
        r=resonance_proximity_updates(200,60,3,widths=(.005,0))
        self.assertAlmostEqual(r['dm_resonance_h1_abs_gap_over_width'],5.09/.005)
        self.assertIsNone(r['dm_resonance_h2_abs_gap_over_width'])
    def test_EWPT_selection_independent_of_DM_vacuum_RG(self):
        r=profile_updates(dict(M2=200,thc=True,hb=True,hs=True,ewpo=True,wmass=True,dm=False,vacuum_tree_global=False,rg_bfb=False))
        self.assertIs(r['ewpt_eligible'],True)
        r=profile_updates(dict(M2=200,thc=True,hb=None,hs=True,ewpo=True,wmass=True,dm=None))
        self.assertIsNone(r['experimental_subset']);self.assertIs(r['ewpt_eligible'],True)
    def test_column_alignment_with_nullable_values(self):
        for cmb in (False,True):
            columns=output_columns({'test':1},planck_cmb=cmb)
            self.assertEqual(len(columns),len(set(columns)))
            payload={key:None for key in columns};payload['constraint_version']=PHYSICS_VERSION;payload['dm']=False
            values=output_row(payload,{'test':1},planck_cmb=cmb).split('\t')
            self.assertEqual(len(values),len(columns));r=dict(zip(columns,values))
            self.assertEqual(r['constraint_version'],PHYSICS_VERSION);self.assertEqual(r['dm'],'False');self.assertEqual(r['ewpt_gw_candidate'],'nan')

if __name__=='__main__':unittest.main()
