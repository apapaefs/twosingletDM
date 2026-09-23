#!/usr/bin/env python3
"""Materialize the shared input convention in the LanHEP defaults and example card."""
import sys,re
from pathlib import Path
root=Path(__file__).resolve().parents[2];sys.path.insert(0,str(root))
from trsm_inputs import micromegas_sm_inputs
values=micromegas_sm_inputs()
p=root/'DM/models/lanhep_mdl/TRSM_mixed.mdl';s=p.read_text()
for name in ('EE','SW','MW','Mtp'):
    s,count=re.subn(r'\b'+name+r'\s*=\s*[-+0-9.eE]+',lambda m:name+' = '+format(values[name],'.17g'),s,count=1)
    if count!=1:raise ValueError('Missing LanHEP parameter '+name)
s,count=re.subn(r'mass Mh\s*=\s*[-+0-9.eE]+',lambda m:'mass Mh = '+format(values['Mh'],'.17g'),s,count=1)
if count!=1:raise ValueError('Missing LanHEP Higgs mass')
s=s.replace('Electromagnetic coupling constant (<->1/128)','Electromagnetic coupling: shared GF MW MZ scheme')
s=s.replace('sin of the Weinberg angle 0.474 - "on-shell",481 - "MS-bar" )','sin(thetaW): shared MW MZ on-shell scheme')
p.write_text(s)
p=root/'DM/data.par';card={k:v for k,v in (line.split() for line in p.read_text().splitlines() if line.strip())};card.update({k:format(v,'.17g') for k,v in values.items()})
p.write_text(''.join(f'{k} {v}\n' for k,v in card.items()))
