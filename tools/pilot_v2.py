#!/usr/bin/env python3
"""Real interrupted/resumed scan, fingerprint rejection, DM-failing EWPT and plots."""
import sys,os,json,csv,time,signal,subprocess,hashlib,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from plot_trsm_v2 import render_suite

def read(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def run(destination):
    dest=Path(destination).resolve();dest.mkdir(parents=True,exist_ok=False)
    initial=dest/'interrupted';reference=dest/'reference';explicit=dest/'explicit'
    for p in (initial,reference,explicit):p.mkdir()
    base=[sys.executable,str(ROOT/'generate_trsm_points.py'),'917319','--independent-m3','--m2-min','180','--m2-max','220',
          '--m3-min','50','--m3-max','55','--vs-min','290','--vs-max','310','--k1-min','.99','--k1-max','1',
          '--lx-min','.05','--lx-max','.2','--checkpoint-every','1','--no-print-info']
    manifest=initial/'outputs.json';command=base+['--nrandom','20','--output-manifest',str(manifest)]
    with (dest/'interrupted.log').open('w') as log:
        process=subprocess.Popen(command,cwd=initial,stdout=log,stderr=subprocess.STDOUT)
        deadline=time.monotonic()+300
        while time.monotonic()<deadline:
            checkpoints=list((initial/'output').glob('*.checkpoint.json'))
            if checkpoints:
                checkpoint=checkpoints[0];state=json.loads(checkpoint.read_text())
                if state['draw_count']>=1:break
            if process.poll() is not None:raise RuntimeError('Generator exited before interruption; inspect interrupted.log')
            time.sleep(.1)
        else:process.terminate();process.wait();raise RuntimeError('Pilot did not checkpoint within five minutes')
        process.send_signal(signal.SIGTERM);process.wait(timeout=120)
    state=json.loads(checkpoint.read_text());assert state['status']=='interrupted',state
    count=state['draw_count']+2
    output=Path(json.loads(manifest.read_text())['outputs']['main']['path'])
    def execute(label,command,cwd,env=None,expect=0):
        with (dest/(label+'.log')).open('w') as log:
            done=subprocess.run(command,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=1800)
        if expect==0 and done.returncode!=0:raise RuntimeError(label+' failed; see log')
        if expect!=0 and done.returncode==0:raise RuntimeError(label+' was unexpectedly accepted')
        return done.returncode
    execute('resume',[sys.executable,str(ROOT/'generate_trsm_points.py'),'--resume-from',str(output),'--nrandom',str(count)],initial)
    execute('reference',base+['--nrandom',str(count),'--output-manifest',str(reference/'outputs.json')],reference)
    reference_output=Path(json.loads((reference/'outputs.json').read_text())['outputs']['main']['path'])
    a,b=read(output),read(reference_output);assert len(a)==len(b)==count
    keys=('M2','M3','vs','a12','lX','lPhiX','lSX','thc','experimental_subset','dm','vacuum_tree_global','rg_bfb','rg_unitarity')
    for x,y in zip(a,b,strict=True):
        for key in keys:assert x[key]==y[key],(key,x[key],y[key])
        for key in ('dm_omega','dm_xf','dm_dir_det'):
            assert math.isclose(float(x[key]),float(y[key]),rel_tol=2e-6,abs_tol=1e-30),(key,x[key],y[key])
    before=hashlib.sha256(output.read_bytes()).hexdigest()
    rejected=execute('incompatible-resume',[sys.executable,str(ROOT/'generate_trsm_points.py'),'--resume-from',str(output),'--nrandom',str(count+1)],initial,dict(os.environ,TRSM_LEGACY_VIRTUAL_OFF='1'),expect=1)
    assert hashlib.sha256(output.read_bytes()).hexdigest()==before
    message=(dest/'incompatible-resume.log').read_text().lower()
    assert any(word in message for word in ('configuration','settings differ','fingerprint'))
    point=read(ROOT/'benchmarks/v2/point_223615/TRSM_Input.tsv')[0]
    command=[sys.executable,str(ROOT/'generate_trsm_points.py'),'917320','--run-ewpt','--ewpt-thigh','300','--no-print-info','--output-manifest',str(explicit/'outputs.json')]
    for key in ('m2','m3','vs','a12','lx','lphix','lsx'):command+=['--'+key,point[key]]
    execute('explicit-ewpt',command,explicit)
    point_output=Path(json.loads((explicit/'outputs.json').read_text())['outputs']['main']['path']);record=read(point_output)[0]
    assert record['dm']=='False' and record['ewpt_eligible']=='True'
    assert record['ewpt_baryo_candidate']=='False' and record['ewpt_gw_candidate']=='True'
    paths=render_suite(point_output,dest/'plots',ewpt_root=explicit/'output')
    assert len(paths)>=10
    report={'status':'passed','interrupted_after_draws':state['draw_count'],'resumed_draw_count':count,'matches_uninterrupted':True,
            'incompatible_resume_rejected':True,'output_unchanged_after_rejection':True,'point223615':{k:record[k] for k in ('dm','ewpt_eligible','ewpt_baryo_candidate','ewpt_gw_candidate')},
            'plot_files':[str(p) for p in paths]}
    (dest/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':run(sys.argv[1])
