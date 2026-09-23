#!/usr/bin/env python3
"""Print or launch the versioned next-scan configuration in a fresh campaign directory."""
import argparse,json,sys,subprocess,shlex,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from trsm_inputs import PHYSICS_VERSION
p=argparse.ArgumentParser();p.add_argument('--config',type=Path,default=ROOT/'config/next-scan-v2.json');p.add_argument('--campaign-dir',type=Path,required=True)
p.add_argument('--pilot',action='store_true',help='Two workers, two draws each; same physics configuration')
p.add_argument('--run',action='store_true',help='Launch after printing the exact command');a=p.parse_args()
c=json.loads(a.config.read_text());assert c['constraint_version']==PHYSICS_VERSION
path=a.campaign_dir.expanduser().resolve()
command=[sys.executable,str(ROOT/'run_trsm_seed_campaign.py'),'--seed-start',str(c['seed_start']),
         '--nseeds',str(2 if a.pilot else c['nseeds']),'--nrandom',str(2 if a.pilot else c['draws_per_seed']),
         '--jobs',str(2 if a.pilot else c['jobs']),'--campaign-dir',str(path),'--run-cwd',str(path),
         '--python-executable',sys.executable,'--run-ewpt','--ewpt-thigh',str(c['ewpt_thigh_GeV'])]
command += ['--generator-extra-arg='+str(arg) for arg in c['generator_arguments']]
print(shlex.join(command),flush=True)
if a.run:
    if path.exists():raise SystemExit('Use a fresh campaign directory; resume existing seed outputs with generate_trsm_points.py --resume-from <manifest output path>')
    path.mkdir(parents=True)
    (path/'requested-configuration.json').write_text(json.dumps(c,indent=2)+'\n')
    raise SystemExit(subprocess.call(command,cwd=ROOT))
