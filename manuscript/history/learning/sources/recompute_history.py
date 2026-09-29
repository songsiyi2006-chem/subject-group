"""Deterministic post-hoc arithmetic; stdlib only; no model, QM or source execution.

Validates snapshot hashes, recomputes quantities used by Tables L2-L5 and writes
figure-ready rows. Source row numbers count the CSV header as row 1. JSON pointers
use zero-based indices. Time-series aggregate summaries are not raw observations.
"""
from pathlib import Path
from collections import Counter
import csv, hashlib, itertools, json, math, statistics

HERE=Path(__file__).resolve().parent
def load(f): return json.loads((HERE/f).read_text(encoding='utf8'))
def rows(f):
    with (HERE/f).open(encoding='utf8',newline='') as h: return list(csv.DictReader(h))
def dump(f,o): (HERE/f).write_text(json.dumps(o,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
def table(f,rs):
    with (HERE/f).open('w',encoding='utf8',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(rs[0]));w.writeheader();w.writerows(rs)
def close(a,b,tol=1e-11):
    assert abs(a-b)<=tol,(a,b)
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()

audit=json.loads((HERE.parent/'audit.json').read_text(encoding='utf8'))
for rec in audit['source_artifacts']:
    assert digest(HERE.parent/rec['snapshot'])==rec['sha256'],rec['id']
dedup=load('deduplication.json');assert all(r['identical'] for r in dedup['shared_files'])
assert not dedup['aqueous_csv_json_prediction_or_model_artifacts']

original={r['id']:r for r in load('phase01_original.json')['results']}
conformers=[]
for i,r in enumerate(load('phase01_audited.json')['results']):
    if r['status']!='ok': continue
    c=r['conformers']; old=original[r['id']]['conformers']
    for k in ['e_min','e_max','delta_e']: close(c[k],old[k])
    conformers.append({'id':r['id'],'name':r['name'],'force_field':c['ff'],
      'accepted_optimization_records':c['n_accepted'],'embedded':c['n_embedded'],
      'Emin_kcal_mol':c['e_min'],'Emax_kcal_mol':c['e_max'],
      'range_kcal_mol':c['delta_e'],'geometric_IMHB_count':c['imhb']['count'],
      'fragments':r['n_fragments'],'source_pointer':f'/results/{i}/conformers'})
table('derived_conformers.csv',conformers)

vmc=[];blockrows=[]
for name,s in load('phase11_results.json')['systems'].items():
    b=s['blocks'];e=statistics.mean(b);se=statistics.stdev(b)/math.sqrt(len(b))
    close(e,s['E_final']);close(se,s['E_err'])
    error=1000*(e-s['E_fci']);close(error,s['d_fci'])
    vmc.append({'system':name,'R_bohr':name.split('R')[-1] if name.startswith('H2') else '',
      'E_VMC_Hartree':e,'E_reference_Hartree':s['E_fci'],'signed_error_mEh':error,
      'block_SE_mEh':1000*se,'threshold_mEh':1.6,'point_gate_pass':abs(error)<1.6,
      'block_count':len(b),'epochs':s['epochs'],'walkers_per_system':s['n_walkers'],
      'variance_initial_Hartree2':s['var_first'],'variance_final_Hartree2':s['var_final'],
      'source_pointer':'/systems/'+name})
    for i,x in enumerate(b): blockrows.append({'system':name,'block_index':i,'energy_Hartree':x})
table('derived_vmc.csv',vmc);table('derived_vmc_blocks.csv',blockrows)

trace=[];windows=rows('phase28_window_robustness.csv')
for i,r in enumerate(rows('phase28_trace_metrics.csv')):
    n=int(r['n_source_rows']);u=int(r['unique_timestamps']);d=int(r['duplicate_time_count'])
    assert n-u==d
    w=[x for x in windows if x['trace_id']==r['trace_id']]
    late=[x for x in w if float(x['startup_exclusion_h'])>=10]
    trace.append({'trace_id':r['trace_id'],'source_DOI':r['source_DOI'],
      'n_rows':n,'unique_timestamps':u,'duplicate_rows':d,'duplicate_fraction':d/n,
      'duration_h':float(r['duration_h']),'unit':r['y_unit'],
      'timestamp_balanced_1h_change':float(r['timestamp_balanced_change']),
      'all_window_min':min(float(x['change']) for x in w),
      'all_window_max':max(float(x['change']) for x in w),
      'startup_ge10h_min':min(float(x['change']) for x in late),
      'startup_ge10h_max':max(float(x['change']) for x in late),
      'source_csv_row':i+2})
table('derived_trace_summary.csv',trace)

comparison=rows('phase29_comparison.csv');centered=[]
for i,r in enumerate(comparison):
    if r['basis']!='def2-svp':continue
    q=next(x for x in comparison if x['molecule_id']=='Q01' and x['basis']==r['basis'] and x['gfn']==r['gfn'])
    dft=float(r['PBE0_delta_E_eV']);gfn=float(r['GFN_delta_E_eV'])
    dd=dft-float(q['PBE0_delta_E_eV']);gd=gfn-float(q['GFN_delta_E_eV'])
    close(gfn-dft,float(r['GFN_minus_PBE0_eV']))
    centered.append({'molecule_id':r['molecule_id'],'matrix_id':r['matrix_id'],'basis':r['basis'],
      'gfn':int(r['gfn']),'geometry_sha256':r['geometry_sha256'],
      'PBE0_raw_delta_eV':dft,'GFN_raw_delta_eV':gfn,'raw_GFNmPBE0_eV':gfn-dft,
      'PBE0_Q01_centered_eV':dd,'GFN_Q01_centered_eV':gd,'double_difference_eV':gd-dd,
      'source_csv_row':i+2})
table('derived_reference_alignment.csv',centered)
max_dft_discrepancy=0.0
for r in rows('phase29_dft_differences.csv'):
    stem='phase29_'+r['molecule_id']+'_'+r['basis']
    cation=load(stem+'_cation.json');radical=load(stem+'_radical.json')
    assert cation['geometry_sha256']==radical['geometry_sha256']==r['geometry_sha256']
    assert cation['status']==radical['status']=='converged'
    assert cation['charge']==1 and radical['charge']==0
    assert cation['multiplicity']==1 and radical['multiplicity']==2
    value=(radical['energy_hartree']-cation['energy_hartree'])*27.211386245988
    error=abs(value-float(r['delta_E_radical_minus_cation_eV']))
    max_dft_discrepancy=max(max_dft_discrepancy,error)
    close(value,float(r['delta_E_radical_minus_cation_eV']),2e-10)

# Independently rebuild 72 xTB charge-state differences from the 144 energy ledger.
manifest=load('phase29_xtb_manifest.json');jobs=manifest['jobs'];state=rows('phase29_state_differences.csv')
hartree_to_eV=27.211386245988
max_state_discrepancy=0.0
for r in state:
    js=[j for j in jobs if j['molecule_id']==r['molecule_id'] and j['gfn']==int(r['gfn']) and j['environment']==r['environment']]
    assert len(js)==2
    assert len({j['geometry_sha256'] for j in js})==1
    jdict={j['state']:j for j in js}
    value=(jdict['radical']['energy_hartree']-jdict['cation']['energy_hartree'])*hartree_to_eV
    discrepancy=abs(value-float(r['delta_model_eV']));max_state_discrepancy=max(max_state_discrepancy,discrepancy)
    close(value,float(r['delta_model_eV']),2e-10)
rankrows=[]
settings=sorted({(r['gfn'],r['environment']) for r in state})
for a,b in itertools.combinations(settings,2):
    da={r['molecule_id']:float(r['delta_model_eV']) for r in state if (r['gfn'],r['environment'])==a}
    db={r['molecule_id']:float(r['delta_model_eV']) for r in state if (r['gfn'],r['environment'])==b}
    ids=sorted(set(da)&set(db)); rev=0;tie=0
    for x,y in itertools.combinations(ids,2):
        delta_a,delta_b=da[x]-da[y],db[x]-db[y]
        if abs(delta_a)<=1e-6 or abs(delta_b)<=1e-6:tie+=1
        elif delta_a*delta_b<0:rev+=1
    rankrows.append({'setting_a':'GFN'+a[0]+'_'+a[1],'setting_b':'GFN'+b[0]+'_'+b[1],
                     'common_molecules':len(ids),'pairs':len(ids)*(len(ids)-1)//2,'rank_reversals':rev,'ties':tie})
table('derived_rank_reversals.csv',rankrows)

meas=load('phase24_measurements.json');campaign=[]
for round_id in sorted({m['round'] for m in meas}):
    m=[r for r in meas if r['round']==round_id]
    assert all(x['provenance']=='simulated' for x in m)
    campaign.append({'round':round_id,'synthetic_measurements':len(m),
      'qc_passed':sum(x['qc_pass'] for x in m),'max_yield_pct':max(x['yield_pct'] for x in m),
      'max_abs_ee_pct':max(x['ee_abs_pct'] for x in m)})
table('derived_synthetic_campaign.csv',campaign)

q01svp=next(r for r in comparison if r['molecule_id']=='Q01' and r['basis']=='def2-svp')
q01tzvp=next(r for r in comparison if r['molecule_id']=='Q01' and r['basis']=='def2-tzvp')
summary={'no_new_quantum_or_training':True,'copied_snapshot_hashes_verified':len(audit['source_artifacts']),
 'deduplication':{'shared_ancestral_commit':dedup['common_root_commit'],'identical_code_and_figure_blobs':len(dedup['shared_files']),
                 'aqueous_recomputable_prediction_artifacts':0},
 'phase01':{'valid_structures':len(conformers),'failed_original_inputs':1,
            'accepted_optimization_records':sum(r['accepted_optimization_records'] for r in conformers),
            'forcefield_structure_counts':dict(Counter(r['force_field'] for r in conformers)),
            'old_and_audited_rounded_energies_equal':True,'independent_new_compounds_from_rerun':0},
 'phase11':{'systems':len(vmc),'production_blocks':len(blockrows),'epochs':sum(r['epochs'] for r in vmc),
            'point_gate_passed':sum(r['point_gate_pass'] for r in vmc),'point_gate_failed':sum(not r['point_gate_pass'] for r in vmc),
            'reference':'computed FCI CBS extrapolation; not error-free continuum truth'},
 'phase24':{'synthetic_measurements':len(meas),'qc_passed':sum(m['qc_pass'] for m in meas),'actual_experiments':0},
 'phase28':{'published_trace_rows':sum(r['n_rows'] for r in trace),
            'unique_timestamps_summed_within_traces':sum(r['unique_timestamps'] for r in trace),
            'duplicate_rows':sum(r['duplicate_rows'] for r in trace),'independent_electrode_count':None,
            'window_specifications':len(windows),'published_mean_sd_rows':len(rows('phase28_production_fe_inputs.csv'))},
 'phase29':{'matrix_xtb_calls':len(jobs),'matrix_unique_settings':len({(j['molecule_id'],j['gfn'],j['environment'],j['state']) for j in jobs}),
            'matrix_warning_calls':sum(bool(j['numerical_warnings']) for j in jobs),'matrix_state_pairs':len(state),
            'initial_xtb_calls_separate':12,'repeated_initial_settings_within_matrix':12,
            'cumulative_xtb_executions':156,'unique_xtb_settings_across_initial_and_matrix':144,
            'new_distinct_settings_in_extension':132,
            'matrix_energy_arithmetic_max_error_eV':max_state_discrepancy,
            'DFT_energy_arithmetic_max_error_eV':max_dft_discrepancy,
            'DFT_jobs':len(load('phase29_dft_manifest.json')['jobs']),
            'DFT_converged':sum(j['status']=='converged' for j in load('phase29_dft_manifest.json')['jobs']),
            'DFT_state_pairs':len(rows('phase29_dft_differences.csv')),'matched_DFT_GFNs_rows':len(comparison),
            'Q01_basis_increment_eV':float(q01tzvp['PBE0_delta_E_eV'])-float(q01svp['PBE0_delta_E_eV'])},
 'paclitaxel_conversion_check':{'reported_MW_g_mol':853.9,'reported_mass_range_ug_mL':[0.3,1.0],
      'converted_log10_mol_L':[math.log10(x/(853.9*1000)) for x in [.3,1]],
      'source_printed_range':[-6.1,-5.6],'interpretation':'Unit arithmetic only; literature mass range and identity not independently validated.'},
 'limits':['No historical model retrained or sampled.','No missing solubility predictions inferred from plots.',
           'Trace statistics recomputed from published aggregate snapshots, not independent electrodes.',
           'No bootstrapping or confidence intervals newly computed; saved block SE independently reproduced.',
           'Raw computational logs not independently re-parsed in this history audit; job ledgers and selected stored outputs used.']}
dump('posthoc_summary.json',summary)
audit['posthoc_summary']='sources/posthoc_summary.json'
audit['figure_ready_tables']=['sources/'+x for x in ['derived_conformers.csv','derived_vmc.csv','derived_vmc_blocks.csv',
   'derived_trace_summary.csv','phase28_window_robustness.csv','derived_reference_alignment.csv','derived_rank_reversals.csv','derived_synthetic_campaign.csv']]
audit['integrity_hashes']={'sources/'+p.name:digest(p) for p in sorted(HERE.iterdir()) if p.is_file()}
(HERE.parent/'audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
print(json.dumps(summary,ensure_ascii=True,indent=2))
