#!/usr/bin/env python3
"""Deterministic 100-point sample, with source identities and selection strata."""
import csv,json,hashlib,math
from pathlib import Path
root=Path(__file__).resolve().parents[1]
destination=root/'benchmarks/v2';destination.mkdir(parents=True,exist_ok=True)
keys=('M2','M3','vs','vx','a12','lX','lPhiX','lSX')
candidates={};sources={}
def add(row,source,line):
    try:values=tuple(float(row[k]) for k in keys)
    except (KeyError,ValueError):return
    if not all(math.isfinite(v) for v in values) or values[3]!=0:return
    digest=hashlib.sha256(repr(values).encode()).hexdigest()
    row=dict(row,source_file=source,source_line=str(line),selection_key=digest)
    candidates.setdefault(digest,row)
for path in sorted((root/'output').glob('*.dat')):
    if '20260923' in path.name:continue
    with path.open() as f:
        reader=csv.DictReader(f,delimiter='\t')
        if not set(keys)<=set(reader.fieldnames or []):continue
        sources[str(path.relative_to(root))]=hashlib.sha256(path.read_bytes()).hexdigest()
        for line,row in enumerate(reader,2):add(row,str(path.relative_to(root)),line)
for path in sorted((root/'plots/seed66666_xbroken_phases/sources').glob('point_*/ewpt_result.json')):
    p=json.loads(path.read_text());r=p['calctemps'];old=p['transition_strengths']
    row=dict(zip(keys,(r['m2'],r['m3'],r['vs'],0,r['a12'],r['lx'],r['lphix'],r['lsx'])))
    row['stored_transition_candidate']=True
    row['legacy_gw_critical_ratio']=max((s['field_jump_over_T'] for s in old if s['temperature_kind']=='crit'),default=0)
    row['original_point_index']=path.parent.name.removeprefix('point_')
    sources[str(path.relative_to(root))]=hashlib.sha256(path.read_bytes()).hexdigest()
    add(row,str(path.relative_to(root)),1)
rows=list(candidates.values()); selected=[];used=set()
def choose(name,rows,n):
    for row in rows:
        key=row['selection_key']
        if key in used:continue
        selected.append(dict(row,selection_stratum=name));used.add(key)
        n-=1
        if n==0:break
choose('stored_transition',sorted((r for r in rows if r.get('stored_transition_candidate')),key=lambda r:r['selection_key']),20)
choose('h1_resonance',sorted(rows,key=lambda r:(abs(125.09-2*float(r['M3'])),r['selection_key'])),20)
choose('h2_resonance',sorted(rows,key=lambda r:(abs(float(r['M2'])-2*float(r['M3'])),r['selection_key'])),20)
def dd_distance(r):
    try:return abs(math.log(float(r['dm_dir_det'])/float(r['dm_dir_det_limit'])))
    except (KeyError,ValueError,ZeroDivisionError):return math.inf
choose('DD_boundary',sorted(rows,key=lambda r:(dd_distance(r),r['selection_key'])),20)
choose('exclusions_and_general',sorted(rows,key=lambda r:(r.get('dm')!='False',r['selection_key'])),100-len(selected))
assert len(selected)==100
header=list(keys)+sorted(set().union(*(r.keys() for r in selected))-set(keys)-{None})
with (destination/'stored100.tsv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=header,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(selected)
(destination/'selection.json').write_text(json.dumps({'algorithm':'deterministic_stratified_v1','count':100,'source_files_sha256':sources,'output_sha256':hashlib.sha256((destination/'stored100.tsv').read_bytes()).hexdigest()},indent=2)+'\n')
print('Selected',len(selected),'from',len(rows),'stored points')
