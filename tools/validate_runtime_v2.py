#!/usr/bin/env python3
"""Matched backend/virtual-boson benchmarks, recording solver and hook evidence."""
import sys,os,json,re,subprocess,tempfile,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from test_trsm_DM import DMPoint,write_micromegas_card,parse_micromegas_output
from trsm_micromegas import default_micromegas_main
from trsm_inputs import M1,MW,MZ,json_safe

def run(destination):
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    cases=[('h1_resonance',M1/2,200,.1,.01,.02),('h2_resonance',100,200,.15,.02,.03),
           ('W_below',MW-3,230,.2,.03,.01),('W_above',MW+3,230,.2,.03,.01),
           ('Z_below',MZ-3,250,.2,.03,.01),('Z_above',MZ+3,250,.2,.03,.01),
           ('heavy',200,420,.1,.02,.03)]
    records=[]
    for version in ('7.1.4','6.1.15'):
      with tempfile.TemporaryDirectory(prefix='trsm-runtime-validation-') as tmp:
       for name,mx,m2,angle,lhx,lsx in cases:
        for virtual in (True,False):
          card=Path(tmp)/'point.dat';write_micromegas_card(DMPoint(.1,lhx,lsx,mx,300,angle,m2),card)
          env=dict(os.environ,TRSM_RUNTIME_DIR=tmp)
          env.pop('TRSM_LEGACY_VIRTUAL_OFF',None)
          if not virtual:env['TRSM_LEGACY_VIRTUAL_OFF']='1'
          done=subprocess.run([str(default_micromegas_main(version)),str(card),'--planck-cmb'],cwd=tmp,env=env,capture_output=True,text=True,timeout=240,check=True)
          label=f'{version}-{name}-virtual-{int(virtual)}';(destination/(label+'.log')).write_text(done.stdout+done.stderr)
          result=parse_micromegas_output(done.stdout)
          assert result.solver_error==0,(label,result.solver_error)
          assert result.hook_calls and min(result.hook_calls)>0,(label,result.hook_calls)
          assert result.inputs['Mh']==M1 and result.mdm==mx
          records.append({'backend':version,'case':name,'virtual_WZ':virtual,'omega':result.omega,'Xf':result.xf,'Tf_GeV':result.mdm/result.xf,
                          'SI_pb':result.dir_det,'hook_calls':result.hook_calls,'inputs':result.inputs})
          print(label,format(result.omega,'.12g'),result.hook_calls,flush=True)
          (destination/'runtime-comparison.json').write_text(json.dumps(json_safe(records),indent=2,allow_nan=False)+'\n')
    for record in records:
        counterpart=next(r for r in records if r['case']==record['case'] and r['backend']==record['backend'] and r['virtual_WZ']!=record['virtual_WZ'])
        record['omega_virtual_on_over_off']=(record['omega']/counterpart['omega'] if record['virtual_WZ'] else counterpart['omega']/record['omega'])
    (destination/'runtime-comparison.json').write_text(json.dumps(records,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':run(sys.argv[1])
