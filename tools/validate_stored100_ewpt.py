#!/usr/bin/env python3
"""Bounded BSMPT validation; cache raw execution separately from interpretation."""
import sys,csv,json,subprocess,hashlib,math
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import test_trsm_ewpt as ew
from ewpt_assessment import status_updates
from ewpt_entry_criterion import ew_entry_updates
from ewpt_x_history import x_history_updates
from dm_thermal_relic_diagnostic import thermal_input_updates
from trsm_inputs import PHYSICS_VERSION,json_safe
from trsm_constraint_profile import sha256

def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True).encode()).hexdigest()
def finite(value):
    try:value=float(value)
    except (TypeError,ValueError):return None
    return value if math.isfinite(value) else None

def run(source,destination,workers=4):
    source=Path(source);destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    with source.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    executables={name:sha256(ew.DEFAULT_EXECUTABLE.with_name(name)) for name in ('CalcTemps','MinimaTracer','PhaseProbe')}
    old_fingerprint=executables['CalcTemps']+executables['PhaseProbe']
    def assess(index,row):
        directory=destination/f'point_{index:06d}';directory.mkdir(exist_ok=True)
        checkpoint=directory/'assessment.json';payload_path=directory/'ewpt_result.json'
        args={key:float(row[key]) for key in ('M2','M3','vs','a12','lX','lPhiX','lSX')}
        signature=digest({'parameters':args,'executables':executables,'thigh':300,'multithreading':False})
        payload=None
        if checkpoint.exists() and payload_path.exists():
            cached=json.loads(checkpoint.read_text())
            legacy=hashlib.sha256((json.dumps(row,sort_keys=True)+old_fingerprint).encode()).hexdigest()
            if cached.get('execution_signature')==signature or cached.get('signature')==legacy:
                payload=json.loads(payload_path.read_text())
        try:
            if payload is None:
                point=ew.TRSMEWPTPoint(args['M2'],args['M3'],args['vs'],args['a12'],args['lX'],args['lPhiX'],args['lSX'],index=index)
                result=ew.run_trsm_ewpt(point,config=ew.EWPTConfig(thigh=300,use_multithreading=False),workdir=directory,keep_files=True)
                payload=ew.result_to_json(result)
            # Recompute interpretation even when the expensive C++ outputs are cached.
            payload['transition_strengths']=[s.to_dict() for s in ew.calculate_fopt_strengths(payload.get('calctemps') or {})]
            payload_path.write_text(json.dumps(json_safe(payload),indent=2,allow_nan=False))
            updates=status_updates(payload);updates.update(ew_entry_updates(payload))
            updates.update(x_history_updates(payload,finite(row.get('dm_freezeout_temperature_GeV'))))
            updates.update(thermal_input_updates(payload,row));updates['ewpt_constraint_version']=PHYSICS_VERSION
        except (ValueError,RuntimeError,OSError,subprocess.SubprocessError) as error:
            updates={'ewpt_execution_status':'error','ewpt_status':'error','ewpt_error':str(error),'ewpt_baryo_candidate':None,'ewpt_gw_candidate':None}
        checkpoint.write_text(json.dumps(json_safe({'execution_signature':signature,'executables':executables,'updates':updates}),indent=2,allow_nan=False))
        return index,updates
    with ThreadPoolExecutor(max_workers=workers) as pool:
        jobs=[pool.submit(assess,i,row) for i,row in enumerate(rows,1) if row.get('ewpt_eligible')=='True']
        for future in as_completed(jobs):
            i,updates=future.result();rows[i-1].update(updates)
            print(i,updates.get('ewpt_status'),updates.get('ewpt_baryo_candidate'),updates.get('ewpt_gw_candidate'),flush=True)
    columns=list(rows[0]);columns+=sorted(set().union(*(r.keys() for r in rows))-set(columns))
    with (destination/'assessed100.tsv').open('w') as f:
        writer=csv.DictWriter(f,columns,delimiter='\t',lineterminator='\n');writer.writeheader()
        writer.writerows({k:('nan' if v is None else v) for k,v in r.items()} for r in rows)

if __name__=='__main__':run(sys.argv[1],sys.argv[2])
