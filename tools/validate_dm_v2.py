#!/usr/bin/env python3
"""Compare compiled loop amplitudes to independent high-precision formulae."""
import sys,os,json,math,re,subprocess,tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import mpmath as mp
from test_trsm_DM import DMPoint,write_micromegas_card
from trsm_micromegas import default_micromegas_main
from trsm_inputs import M1,VEV,EE,MW
mp.mp.dps=65

def independent_coefficient(q,masses,leptons,alpha_s,channel):
    q=mp.mpf(q)
    def f(t):
        if t<=1:return mp.asin(mp.sqrt(t))**2
        return -(mp.log((1+mp.sqrt(1-1/t))/(1-mp.sqrt(1-1/t)))-mp.j*mp.pi)**2/4
    def fermion(m):
        t=q*q/(4*mp.mpf(m)**2)
        return 2*(t+(t-1)*f(t))/t**2
    result=sum(fermion(m)*(3*c*c if channel==22 else 1) for m,c in zip(masses,[mp.mpf(2)/3,-mp.mpf(1)/3,-mp.mpf(1)/3,mp.mpf(2)/3,-mp.mpf(1)/3,mp.mpf(2)/3]) if m>0)
    if channel==21:return complex(-alpha_s*result/(16*mp.pi*VEV))
    t=q*q/(4*mp.mpf(MW)**2)
    result+=sum(fermion(m) for m in leptons)-(2*t*t+3*t+3*(2*t-1)*f(t))/t**2
    return complex(-EE*EE*result/(32*mp.pi**2*VEV))

def close(actual,expected,relative=2e-8,absolute=1e-22):
    if abs(actual-expected)>max(absolute,relative*abs(expected)):
        raise AssertionError(f'{actual} != {expected}; relative deviation {abs(actual-expected)/abs(expected) if expected else None}')


def run(output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    cases=[('calchep_normalization',60,200,0,.01,.02,math.sqrt(M1*M1/4-60**2)),('generic',60,200,.1,.01,.02,1),('zero_mixing',60,200,0,.01,.02,1),
           ('mixed',60,200,.4,.01,.02,1),('h1_pole',62,200,.25,.01,.03,math.sqrt(M1*M1/4-62**2)),
           ('h2_pole',99,200,.25,.01,.03,math.sqrt(100**2-99**2)),
           ('W_below',80,220,.2,.03,.01,.1),('W_above',81,220,.2,.03,.01,.1),
           ('top_below',172,410,.2,.03,.01,.1),('top_above',174,410,.2,.03,.01,.1),
           ('degenerate_cancellation',60,M1,math.pi/4,0,.03,1)]
    records=[]
    for version in ('7.1.4','6.1.15'):
      with tempfile.TemporaryDirectory(prefix='trsm-loop-validation-') as temporary:
       root=Path(temporary)
       for name,mx,m2,angle,lhx,lsx,pcm in cases:
        card=root/'card.dat';write_micromegas_card(DMPoint(.1,lhx,lsx,mx,300,angle,m2),card)
        env=dict(os.environ,TRSM_RUNTIME_DIR=temporary,TRSM_LOOP_PROBE=str(pcm))
        if name=='calchep_normalization':env['TRSM_CALCHEP_PROBE']='1'
        done=subprocess.run([str(default_micromegas_main(version)),str(card)],env=env,cwd=root,check=True,text=True,capture_output=True,timeout=180)
        (output/f'{version}-{name}.log').write_text(done.stdout+done.stderr)
        inputs=json.loads(re.search(r'^TRSM_inputs_v2 (.+)$',done.stdout,re.M)[1])
        params=[json.loads(s) for s in re.findall(r'^TRSM_loop_parameters (.+)$',done.stdout,re.M)]
        for p in params:
          for channel,prefix in [(22,'aa'),(21,'gg')]:
            reference=independent_coefficient(p['q'],p['quarks'],p['leptons'],p['alpha_s'],channel)
            close(complex(p[prefix+'_re'],p[prefix+'_im']),reference)
        widths=[float(x) for x in re.findall(r'TRSM_partial_width h\d aa=(\S+) gg=\S+',done.stdout)]
        gluon_widths=[float(x) for x in re.findall(r'TRSM_partial_width h\d aa=\S+ gg=(\S+)',done.stdout)]
        for i,(mass,proj) in enumerate([(M1,math.cos(angle)),(m2,math.sin(angle))]):
          p=params[i+1]
          for channel,actual in [(22,widths[i]),(21,gluon_widths[i])]:
            coeff=independent_coefficient(mass,p['quarks'],p['leptons'],p['alpha_s'],channel)*proj
            expected=(8 if channel==21 else 1)*abs(coeff)**2*mass**3/(4*math.pi)
            close(actual,expected,relative=2e-6) # findBr parses the backend's rounded BR text
        p=params[0];q=p['q'];s=q*q;c,st=math.cos(angle),math.sin(angle)
        g1=lhx*VEV*c-lsx*300*st;g2=lhx*VEV*st+lsx*300*c
        propagators=g1*c/complex(s-M1*M1,M1*inputs['width_h1'])+g2*st/complex(s-m2*m2,m2*inputs['width_h2'])
        probe=re.search(r'TRSM_loop_probe pcm=\S+ aa=(\S+) gg=(\S+)',done.stdout)
        for channel,actual in [(22,float(probe[1])),(21,float(probe[2]))]:
          coefficient=independent_coefficient(q,p['quarks'],p['leptons'],p['alpha_s'],channel)
          reference=(8 if channel==21 else 1)*s*abs(coefficient*propagators)**2/(4*math.pi*(2*pcm/q))
          close(actual,reference,absolute=1e-22)
        if name=='calchep_normalization':
          raw=float(re.search(r'TRSM_calchep_probe aa_pb=(\S+)',done.stdout)[1])
          close(raw,float(probe[1])*3.8937966e8)
          raw_gg=float(re.search(r'TRSM_calchep_probe gg_pb=(\S+)',done.stdout)[1])
          close(raw_gg,float(probe[2])*3.8937966e8)
        records.append({'backend':version,'case':name,'inputs':inputs,'aa_partial_widths':widths,'gg_partial_widths':gluon_widths,'status':'passed'})
        print(version,name,'passed',flush=True)
    (output/'independent-amplitude-report.json').write_text(json.dumps(records,indent=2)+'\n')

if __name__=='__main__':run(sys.argv[1] if len(sys.argv)>1 else 'validation/v2/dm-amplitudes')
