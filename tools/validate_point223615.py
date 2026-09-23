#!/usr/bin/env python3
"""Fresh transition calculation and common-temperature check of archived point 223615."""
import sys,csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import test_trsm_ewpt as ew
from ewpt_equilibrium import PhaseProbe,refine_equilibrium
from ewpt_entry_criterion import ew_entry_updates
from ewpt_assessment import status_updates
from trsm_inputs import json_safe
from trsm_constraint_profile import sha256

def run(output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    fixture=ROOT/'benchmarks/v2/point_223615'
    with (fixture/'TRSM_Input.tsv').open() as f:row=next(csv.DictReader(f,delimiter='\t'))
    point=ew.TRSMEWPTPoint(**{k:float(row[k]) for k in ('m1','m2','m3','vs','a12','lx','lphix','lsx')},index=223615)
    fresh=ew.run_trsm_ewpt(point,config=ew.EWPTConfig(thigh=300,use_multithreading=False),workdir=output/'fresh',keep_files=True)
    payload=ew.result_to_json(fresh);(output/'fresh/ewpt_result.json').write_text(json.dumps(payload,indent=2,allow_nan=False))
    flags=ew_entry_updates(payload);assert flags['ewpt_baryo_candidate'] is False and flags['ewpt_gw_candidate'] is True
    phases=ew.parse_minimatracer_output(fixture/'minima_trace_1.tsv')
    point_file=output/'point.tsv';ew.write_trsm_input(point_file,point)
    with PhaseProbe(ew.DEFAULT_EXECUTABLE.with_name('PhaseProbe'),point_file,output/'phase-probe.log') as probe:
        archive=refine_equilibrium(phases,ew.PhaseClassificationThresholds(),probe)
    crossing=next(c for c in archive.crossings if 51.5<c['low_T_GeV']<51.6)
    assert crossing['status']=='resolved' and crossing['high_T_GeV']-crossing['low_T_GeV']<=.001
    assert crossing['low_T_GeV']<=51.5526<=crossing['high_T_GeV']
    report={'status':'passed','flags':flags,'calculation':status_updates(payload),'crossing':crossing,'equilibrium_status':archive.equilibrium_status,
            'executables':{name:sha256(ew.DEFAULT_EXECUTABLE.with_name(name)) for name in ('CalcTemps','MinimaTracer','PhaseProbe')}}
    (output/'refined_archive.json').write_text(json.dumps(json_safe(archive.to_dict()),indent=2,allow_nan=False))
    (output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':run(sys.argv[1])
