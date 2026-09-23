"""Common-temperature local-minimum refinement, independent of bubble history.

Only traced basins are compared. Resolution refers to their ordering, not a
proof that the tracer discovered every possible minimum. Unknown ordering is
preserved if refinement, local stability or numerical depth separation fails.
"""
import math
import select
import subprocess
from pathlib import Path

CROSSING_BRACKET_GEV = 1e-3
MINIMUM_MATCH_GEV = .05

class PhaseProbe:
    def __init__(self, executable, point_file, log_file):
        self.log = Path(log_file).open('w')
        self.process = subprocess.Popen([str(executable), str(point_file)], stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE, stderr=self.log, text=True, bufsize=1)
    def __enter__(self): return self
    def __exit__(self, *args):
        self.process.stdin.close()
        try: self.process.wait(timeout=5)
        except subprocess.TimeoutExpired: self.process.kill(); self.process.wait()
        self.process.stdout.close(); self.log.close()
    def refine(self, temperature, phase, sample):
        p=self.process
        p.stdin.write(f'{temperature:.17g} {phase} {sample.w1:.17g} {sample.wx:.17g} {sample.ws:.17g}\n');p.stdin.flush()
        if not select.select([p.stdout], [], [], 60)[0]: raise RuntimeError('PhaseProbe timed out')
        line=p.stdout.readline().split()
        if len(line)!=11 or line[0]!='PROBE': raise RuntimeError('Invalid PhaseProbe response: '+' '.join(line))
        values=list(map(float,line[4:]))
        return dict(temp=float(line[1]),phase_index=int(line[2]),status=line[3],
                    w1=values[0],wx=values[1],ws=values[2],veff=values[3],
                    uncertainty=values[4],gradient_norm=values[5],hessian_min=values[6])


def refine_equilibrium(phase_traces, thresholds, probe, bracket=CROSSING_BRACKET_GEV):
    from test_trsm_ewpt import (interpolate_phase_at,GlobalBranchPoint,MinimaTracerAnalysis,
                               classify_phase,compress_labels_for_cooling,find_ew_step_index)
    cache={}; potential_samples=[]
    def evaluate(t):
        if t in cache:return cache[t]
        minima=[]; uncertain=False
        for trace in phase_traces:
            seed=interpolate_phase_at(trace,t)
            if seed is None:continue
            result=probe.refine(t,trace.index,seed);potential_samples.append(result)
            if result['status']=='saddle':continue
            if result['status']!='minimum':uncertain=True;continue
            vector=[result[k] for k in ('w1','wx','ws')]
            # Remove copies related by independent field-sign symmetries.
            if any(math.dist(vector,[old[k] for k in ('w1','wx','ws')])<MINIMUM_MATCH_GEV for old in minima):continue
            minima.append(result)
        minima.sort(key=lambda x:x['veff'])
        best=minima[0] if minima else None
        separated=(len(minima)<2 or minima[1]['veff']-best['veff']>minima[1]['uncertainty']+best['uncertainty'])
        resolved=best is not None and not uncertain and separated
        if best:
            point=GlobalBranchPoint(t,best['phase_index'] if resolved else None,
                classify_phase(best['w1'],best['wx'],best['ws'],thresholds) if resolved else 'UNRESOLVED',
                best['w1'],best['wx'],best['ws'],best['veff'])
        else:point=GlobalBranchPoint(t,None,'UNRESOLVED',math.nan,math.nan,math.nan,math.nan)
        cache[t]=(point,minima);return cache[t]
    grid=sorted({s.temp for trace in phase_traces for s in trace.samples})
    for t in grid:evaluate(t)
    crossings=[]
    for low,high in zip(grid,grid[1:]):
        a=cache[low][0];b=cache[high][0]
        if a.phase_index==b.phase_index:continue
        # Continuous branch reconnections are still bracketed without claiming a FOPT.
        initial=(low,high); left=a.phase_index; right=b.phase_index
        status='resolved'
        while high-low>bracket:
            mid=(low+high)/2;pt,_=evaluate(mid)
            if pt.phase_index is None:
                status='potential_ordering_unresolved';break
            if pt.phase_index==left:low=mid
            elif pt.phase_index==right:high=mid
            else:status='multiple_or_unresolved_branches';break
        if left is None or right is None:status='endpoint_unresolved'
        crossings.append(dict(low_T_GeV=low,high_T_GeV=high,low_phase=left,high_phase=right,
                              initial_bracket_GeV=list(initial),status=status))
    grid=sorted(cache);branch=[cache[t][0] for t in grid];path=compress_labels_for_cooling(branch)
    return MinimaTracerAnalysis(phase_traces,grid,branch,path,find_ew_step_index(path),
        equilibrium_status='resolved' if all(p.phase_index is not None for p in branch) and all(c['status']=='resolved' for c in crossings) else 'partially_unresolved',
        crossings=crossings,potential_samples=potential_samples)
