#!/usr/bin/env python3
"""Plots of independent constraint subsets and explicitly nullable diagnostics."""
import argparse,csv,json,math,html
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

TITLES={
 'v2_candidate_subsets':'Qualitative candidates and independent DM / theory subsets',
 'v2_jump_temperatures':'Critical, nucleation and percolation jump ratios',
 'v2_freezeout_restoration':'Freeze-out and equilibrium X restoration',
 'v2_thermal_sensitivity':'VEV and fixed-mixing coupling sensitivity near resonances',
}

def number(value):
    try:return float(value)
    except (ValueError,TypeError):return math.nan

def verdict(value):
    return True if str(value)=='True' else False if str(value)=='False' else None

def render_suite(source,output,plot_format='both',dpi=150,ewpt_root=None):
    source=Path(source);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    with source.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    if not rows or not any(r.get('constraint_version')=='trsm_constraints_v2' for r in rows):return []
    rows=[r for r in rows if r.get('constraint_version')=='trsm_constraints_v2']
    def n(key):return np.array([number(r.get(key)) for r in rows])
    def b(key):return np.array([verdict(r.get(key)) for r in rows],dtype=object)
    paths=[]
    formats=('png','pdf') if plot_format=='both' else (plot_format,)
    def save(fig,stem):
        fig.tight_layout()
        for ext in formats:
            path=output/f'{stem}.{ext}';fig.savefig(path,dpi=dpi,bbox_inches='tight');paths.append(path)
        plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(13,8),sharex=True,sharey=True)
    x,y=n('M2'),n('M3')
    subsets=[('dm_subset','DM'),('vacuum_tree_global','Tree vacuum'),('rg_unitarity','RG')]
    for i,(flag,title) in enumerate([('ewpt_baryo_candidate','Baryogenesis'),('ewpt_gw_candidate','GW')]):
      candidates=b(flag);passed=candidates==True;unknown=candidates==None
      for j,(subset,label) in enumerate(subsets):
        ax=axes[i,j];state=b(subset)
        if label=='RG':
            state=np.array([False if False in (verdict(r.get('rg_bfb')),verdict(r.get('rg_unitarity')),verdict(r.get('rg_integration_success'))) else True if all(verdict(r.get(k)) is True for k in ('rg_bfb','rg_unitarity','rg_integration_success')) else None for r in rows],dtype=object)
        ax.scatter(x,y,s=9,c='.83',label='all assessed input points')
        ax.scatter(x[unknown],y[unknown],s=12,marker='x',c='#c7862a',label='candidate unassessed')
        ax.scatter(x[passed],y[passed],s=22,c='#347ca5',label=f'{title} candidate')
        exp_unknown=b('experimental_subset')==None
        ax.scatter(x[exp_unknown],y[exp_unknown],s=32,marker='s',facecolors='none',edgecolors='.25',label='experimental unassessed')
        mask=passed&(state==True);ax.scatter(x[mask],y[mask],s=58,facecolors='none',edgecolors='#187343',label=f'+ {label} subset')
        ax.set(xscale='log',yscale='log',title=f'{title}: {label} overlay',xlabel=r'$M_2$ [GeV]',ylabel=r'$M_3$ [GeV]')
        ax.text(.03,.97,f'Candidate unknown: {unknown.sum()}\nSubset unknown: {(state==None).sum()}',transform=ax.transAxes,va='top',fontsize=8)
        ax.legend(fontsize=7,loc='lower right')
    fig.suptitle('Qualitative candidate screens; nucleation failure does not veto positive evidence',fontsize=12)
    save(fig,'v2_candidate_subsets')

    fig,axes=plt.subplots(2,2,figsize=(10,8))
    def paired_strengths(ew_only,kind):
        pairs=[];missing=0;unassessed=0
        for row in rows:
            try:strengths=json.loads(row.get('ewpt_transition_strengths',''))
            except (ValueError,TypeError):unassessed+=1;continue
            by_index={}
            for s in strengths:by_index.setdefault(s['transition_index'],{})[s['temperature_kind']]=s
            for group in by_index.values():
                crit=group.get('crit')
                if not crit:continue
                if ew_only and not (abs(number(crit['false_vev']['w1']))<5<=abs(number(crit['true_vev']['w1']))):continue
                key='ew_jump_over_T' if ew_only else 'field_jump_over_T'
                x=number(crit.get(key));y=number(group.get(kind,{}).get(key))
                if math.isfinite(x) and math.isfinite(y):pairs.append((x,y))
                else:missing+=1
        return np.array(pairs).reshape(-1,2),missing,unassessed
    for i,(critical,prefix,title) in enumerate([('ewpt_ew_entry_jump_over_T','ewpt_ew_entry','EW entry'),('ewpt_gw_crit_field_jump_over_T','ewpt_gw','Any-field jump')]):
      for j,(kind,label) in enumerate([('nucl','Nucleation'),('perc','Percolation')]):
        yc=f'{prefix}_{kind}_'+('jump_over_T' if i==0 else 'field_jump_over_T')
        pairs,missing,unassessed=paired_strengths(i==0,kind);ax=axes[i,j]
        ax.scatter(pairs[:,0],pairs[:,1],s=20,alpha=.7);ax.axhline(1,c='.4',ls='--');ax.axvline(1,c='.4',ls='--')
        ax.set(xlabel=r'Critical $\Delta v/T_c$',ylabel=label+r' $\Delta v/T$',title=title)
        ax.set_xlim(left=0);ax.set_ylim(bottom=0,top=max(1.25,1.15*np.max(pairs[:,1]) if len(pairs) else 1.25))
        ax.text(.03,.97,f'Same transition; missing partner: {missing}\nPoints without transition records: {unassessed}',transform=ax.transAxes,va='top',fontsize=8,bbox=dict(facecolor='white',edgecolor='none',alpha=.9))
    save(fig,'v2_jump_temperatures')

    fig,ax=plt.subplots(figsize=(7,5));tf=n('dm_freezeout_temperature_GeV');lo=n('ewpt_x_final_restoration_low_T_GeV');hi=n('ewpt_x_final_restoration_high_T_GeV')
    valid=np.isfinite(tf)&np.isfinite(lo)&np.isfinite(hi)&(hi>=lo)&(lo>0)&(tf>0)
    mid=(lo+hi)/2;ax.errorbar(tf[valid],mid[valid],yerr=np.array([mid[valid]-lo[valid],hi[valid]-mid[valid]]),fmt='o',ms=4,alpha=.75,capsize=2)
    if valid.any():
        lower=min(tf[valid].min(),lo[valid].min());upper=max(tf[valid].max(),hi[valid].max());ax.plot([lower,upper],[lower,upper],'--',c='.5')
        ax.set(xscale='log',yscale='log')
    ax.set(xlabel=r'$T_f$ [GeV]',ylabel='Equilibrium X-restoration bracket [GeV]',title='Z₂ compatibility diagnostic')
    ax.text(.97,.03,f'Unresolved / no finite restoration bracket: {(~valid).sum()}',transform=ax.transAxes,va='bottom',ha='right',fontsize=9,bbox=dict(facecolor='white',edgecolor='none',alpha=.9))
    save(fig,'v2_freezeout_restoration')

    fig,axes=plt.subplots(2,2,figsize=(11,8))
    gaps=np.fmin(n('dm_resonance_h1_abs_gap_over_width'),n('dm_resonance_h2_abs_gap_over_width'))
    for ax,key,label in zip(axes.flat,['dm_relic_thermal_ew_vev_max_fractional_shift','dm_relic_thermal_s_vev_max_fractional_shift','dm_thermal_K133_max_fractional_change','dm_thermal_K233_max_fractional_change'],['EW VEV','S VEV','K133 fixed-mixing proxy','K233 fixed-mixing proxy']):
        y=n(key);valid=np.isfinite(gaps)&np.isfinite(y)&(gaps>=0)&(y>=0)
        ax.scatter(gaps[valid],y[valid],s=16,alpha=.7);ax.axhline(.1,c='.4',ls='--')
        ax.set_xscale('symlog',linthresh=1);ax.set_yscale('symlog',linthresh=.001)
        ax.set(xlim=(0,max(10,1.3*np.max(gaps[valid])) if valid.any() else 10),ylim=(0,max(.2,1.3*np.max(y[valid])) if valid.any() else .2),xlabel=r'Nearest $|M_i-2M_3|/\Gamma_i$',ylabel='Maximum fractional change',title=label)
        ax.text(.03,.97,f'Unavailable pair: {(~valid).sum()}',transform=ax.transAxes,va='top',fontsize=9,bbox=dict(facecolor='white',edgecolor='none',alpha=.9))
    save(fig,'v2_thermal_sensitivity')
    if ewpt_root:
      for payload_path in sorted(Path(ewpt_root).glob('**/ewpt_result.json')):
        payload=json.loads(payload_path.read_text());analysis=payload.get('minimatracer') or {}
        if not analysis.get('potential_samples'):continue
        stem='v2_trajectory_'+payload_path.parent.name
        if stem in [p.stem for p in paths]:stem+='_'+str(len(paths))
        render_trajectory(payload,output,stem,formats,dpi,paths)
    cards=''.join(f'<section><h2>{html.escape(TITLES.get(p.stem,p.stem))}</h2><a href="{p.name}"><img src="{p.name}" style="max-width:100%"></a></section>' for p in paths if p.suffix=='.png')
    (output/'v2-index.html').write_text('<!doctype html><html><meta charset="utf-8"><title>Constraint v2 diagnostics</title><body style="max-width:1200px;margin:auto;font-family:system-ui"><h1>Constraint v2 diagnostics</h1><p>Vacuum, running and DM results define separate subsets. Phase ordering is equilibrium ordering among traced basins. These plots do not establish a cosmological transition history.</p>'+cards+'</body></html>')
    return paths


def render_trajectory(payload,output,stem,formats,dpi,paths):
    analysis=payload['minimatracer'];fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True)
    colours=['#287aa4','#df8b22','#378047'];styles=['-','--',':'];keys=['w1','wx','ws']
    for trace in analysis['phase_traces']:
        for key,c,style in zip(keys,colours,styles):
            axes[0].plot([r['temp'] for r in trace['samples']],[r[key] for r in trace['samples']],c=c,ls=style,lw=.7,alpha=.35)
    branch=analysis['global_branch'];temps=[r['temp'] for r in branch]
    for key,c,style in zip(keys,colours,styles):
        axes[0].plot(temps,[r[key] if r['phase_index'] is not None else math.nan for r in branch],c=c,ls=style,lw=2,label=f'{key}: resolved equilibrium')
    by_temp={r['temp']:r for r in branch}
    for phase in sorted({r['phase_index'] for r in analysis['potential_samples']}):
        values=[r for r in analysis['potential_samples'] if r['phase_index']==phase and r['status']=='minimum' and by_temp.get(r['temp'],{}).get('phase_index') is not None]
        values.sort(key=lambda r:r['temp'])
        axes[1].plot([r['temp'] for r in values],[r['veff']-by_temp[r['temp']]['veff'] for r in values],label=f'local phase {phase}')
    for cross in analysis.get('crossings',[]):
        for ax in axes:ax.axvspan(cross['low_T_GeV'],cross['high_T_GeV'],color='#a23d32',alpha=.3)
    axes[0].set(ylabel='VEV [GeV]',title='Thin: local minima; thick: resolved equilibrium among traced minima')
    axes[1].set(xlabel='T [GeV]',ylabel=r'$V_i-V_{\rm eq}$ [GeV$^4$]',yscale='symlog')
    for ax in axes:ax.legend(fontsize=8);ax.grid(alpha=.15)
    fig.tight_layout()
    for ext in formats:
        path=Path(output)/f'{stem}.{ext}';fig.savefig(path,dpi=dpi,bbox_inches='tight');paths.append(path)
    plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('--output-dir',required=True);p.add_argument('--ewpt-root');p.add_argument('--format',default='both',choices=['png','pdf','both']);a=p.parse_args()
    render_suite(a.input,a.output_dir,a.format,ewpt_root=a.ewpt_root)
