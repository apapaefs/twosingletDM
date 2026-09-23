#!/usr/bin/env python3
"""Verify installed kernels against this source and write a host-specific receipt."""
import sys,json,os,subprocess,platform,argparse,importlib.metadata,hashlib,importlib.util
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from trsm_constraint_profile import sha256,git_revision
from trsm_inputs import sm_inputs,PHYSICS_VERSION
from trsm_micromegas import default_micromegas_main
def source_tree(directory):
    directory=Path(directory);files=[]
    for folder in ('src','include'):
        files += [p for p in (directory/folder).rglob('*') if p.is_file() and p.suffix in ('.cpp','.cc','.c','.h','.hpp','.in')]
    files += list(directory.glob('CMakeLists.txt'))
    if not files:return None
    data='\n'.join(str(p.relative_to(directory))+':'+sha256(p) for p in sorted(files))
    return hashlib.sha256(data.encode()).hexdigest()
p=argparse.ArgumentParser();p.add_argument('--bsmpt-upstream',type=Path,required=True);p.add_argument('--higgstools-upstream',type=Path,required=True);a=p.parse_args()
root=default_micromegas_main().parents[2]
installed={}
for version in ('7.1.4','6.1.15'):
    exe=default_micromegas_main(version);trsm=exe.parent
    assert (trsm/'main.c').read_bytes()==(ROOT/'DM/main.c').read_bytes(),'Installed driver source differs'
    loop=(trsm/'lib/trsm_loop.c').read_text();canonical=(ROOT/'DM/trsm_loop.c').read_text()
    assert loop==canonical or (version=='7.1.4' and loop.removeprefix('#define TRSM_MO7 1\n\n').removeprefix('#define TRSM_MO7 1\n')==canonical),'Installed loop differs'
    for source in (ROOT/'DM/models/h4GOn').glob('*.mdl'):
        value=(trsm/'work/models'/source.name).read_text()
        assert value.startswith(source.read_text()),f'Installed {version} {source.name} differs'
    installed['micromegas_'+version]={'executable':str(exe),'sha256':sha256(exe),'driver_sha256':sha256(trsm/'main.c'),'loop_sha256':sha256(trsm/'lib/trsm_loop.c'),
                                   'canonical_loop_sha256':sha256(ROOT/'DM/trsm_loop.c'),
                                   'model_sha256':{f.name:sha256(f) for f in sorted((trsm/'work/models').glob('*.mdl'))},
                                   'capabilities':json.loads(subprocess.check_output([str(exe),'--capabilities'],text=True))}
bs=root/'BSMPT-3.2.1'
assert (bs/'src/prog/PhaseProbe.cpp').read_bytes()==(ROOT/'BSMPT/src/prog/PhaseProbe.cpp').read_bytes()
installed['BSMPT']={'upstream_commit':git_revision(a.bsmpt_upstream),'version':'3.2.1','source':str(bs),
                    'source_tree_sha256':source_tree(bs),
                    'executables':{name:{'path':str(bs/'build/bin'/name),'sha256':sha256(bs/'build/bin'/name)} for name in ('CalcTemps','MinimaTracer','PhaseProbe','Test')},
                    'SMparam_sha256':sha256(bs/'src/models/SMparam.cpp'),'PhaseProbe_sha256':sha256(bs/'src/prog/PhaseProbe.cpp')}
payload={'physics_version':PHYSICS_VERSION,'release_source_commit':git_revision(ROOT),'host':platform.node(),'platform':platform.platform(),
         'python':sys.version,'packages':{name:importlib.metadata.version(name) for name in ('numpy','scipy','pandas','matplotlib','HiggsTools','mpmath','PyYAML')},
         'compiler':subprocess.check_output(['c++','--version'],text=True).splitlines()[0],
         'higgstools_source_commit':git_revision(a.higgstools_upstream),'higgstools_source_tree_sha256':source_tree(a.higgstools_upstream),
         'higgstools_extensions':{name:sha256(importlib.util.find_spec('Higgs.'+name).origin) for name in ('predictions','bounds','signals')},'sm_inputs':sm_inputs(),
         'micromegas_archives':{str(v):sha256(root/f'micromegas_{v}.tgz') or sha256(ROOT/f'DM/micromegas_{v}.tgz') for v in ('7.1.4','6.1.15')},'installed':installed}
path=root/'runtime-manifest.json';path.write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n');print(path)
