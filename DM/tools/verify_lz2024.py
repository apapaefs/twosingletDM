#!/usr/bin/env python3
"""Reproduce the normalized published observed SI table from the saved CC0 release."""
import hashlib,json
from pathlib import Path
import yaml
folder=Path(__file__).resolve().parents[1]/'data/lz2024'
raw=(folder/'hepdata-si-v2.yaml').read_bytes();data=yaml.safe_load(raw)
normalized=json.loads((folder/'observed-si-v2.json').read_text())
assert hashlib.sha256(raw).hexdigest()==normalized['provenance']['original_yaml_sha256']
x=next(v for v in data['independent_variables'] if v['header']['name']=='mass')['values']
y=next(v for v in data['dependent_variables'] if v['header']['name']=='limit')['values']
points=[{'mass_GeV':float(a['value']),'upper_limit':float(b['value'])} for a,b in zip(x,y,strict=True)]
assert points==normalized['points']
assert normalized['cross_section_unit']=='cm2' and normalized['cross_section']=='per_nucleon'
assert normalized['confidence_level']==.9 and len(points)==26
assert min(points,key=lambda p:p['upper_limit'])=={'mass_GeV':40.,'upper_limit':2.1816833824484827e-48}
print('Verified 26 observed 90% CL SI points, 9--10000 GeV; 1 pb = 1e-36 cm^2.')
