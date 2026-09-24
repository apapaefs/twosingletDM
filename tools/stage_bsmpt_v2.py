#!/usr/bin/env python3
"""Stage a source-only BSMPT copy, apply the archived model and v2 helpers.

No existing build is modified. Configure the destination with the site's own
Conan toolchain. SM values are generated from the same JSON as the Python code.
"""
import argparse,json,re,shutil,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('upstream',type=Path);p.add_argument('destination',type=Path);a=p.parse_args()
root=Path(__file__).resolve().parents[1]
if a.destination.exists(): raise SystemExit('Destination already exists; choose a fresh versioned source directory')
version=re.search(r'set\(BSMPT_VERSION\s+([^\)]+)\)',(a.upstream/'CMakeLists.txt').read_text())
if not version or version[1].strip()!='3.2.1':raise SystemExit('The validated v2 runtime requires BSMPT 3.2.1 with the archived TRSM model')
subprocess.run(['rsync','-a','--exclude=.git','--exclude=build','--exclude=__pycache__','--exclude=.venv*',str(a.upstream.resolve())+'/',str(a.destination.resolve())+'/'],check=True)
subprocess.run(['rsync','-a',str(root/'BSMPT')+'/',str(a.destination.resolve())+'/'],check=True)
inputs=json.loads((root/'config/sm-inputs-v2.json').read_text())
f=a.destination/'src/models/SMParam.cpp';s=f.read_text()
for name,key in [('C_MassW','MW_GeV'),('C_MassZ','MZ_GeV'),('C_MassSMHiggs','M1_GeV'),('C_GF','GF_GeV^-2')]:
 s,n=re.subn(r'(SM\.'+name+r'\s*=)[^;]+;',lambda m:m[1]+' '+format(inputs[key],'.17g')+';',s)
 if n!=1:raise RuntimeError('Cannot locate SM input '+name)
f.write_text(s)
f=a.destination/'src/CMakeLists.txt';s=f.read_text();s+='\nadd_executable(PhaseProbe prog/PhaseProbe.cpp)\ntarget_link_libraries(PhaseProbe Minimizer Models Utility MinimumTracer)\ntarget_compile_features(PhaseProbe PUBLIC cxx_std_17)\n';f.write_text(s)
# Default ostream precision is insufficient for reevaluation and mass-pole tests.
for name in ['CalcTemps','MinimaTracer']:
 f=a.destination/f'src/prog/{name}.cpp';s=f.read_text();s=s.replace('std::setprecision(16)','std::setprecision(17)')
 s=re.sub(r'(std::ofstream outfile\([^;]+;)',r'\1\n      outfile << std::setprecision(17);',s)
 f.write_text(s)
print(a.destination.resolve())
