"""Read-only extraction of pinned historical sources; never executes source code.

Usage: python collect_history.py <directory containing the two isolated clones>
Only writes into this script's directory and ../audit.json. For arithmetic only,
run recompute_history.py without access to the clones.
"""
from pathlib import Path
import csv, hashlib, json, subprocess, sys

OUT = Path(__file__).resolve().parent
CLONES = Path(sys.argv[1])
REPOS = {
    'aqueous': ('aqueous-solubility-ml-benchmark', '5939547a936ed4baf7fa7ee9fac6fa6e90bfa8fe'),
    'complex': ('ai4chem-complex-scaffolds-benchmark', '5e19c519dbd11430926828fddf5375b919af67f0'),
}
def sha(b): return hashlib.sha256(b).hexdigest()
def git(key,*args):
    return subprocess.check_output(['git','-C',str(CLONES/REPOS[key][0]),*args])
def dump(name,obj): (OUT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
catalog=[]
def copy(key,path,name,kind,selectors=None):
    repo,commit=REPOS[key]
    b=git(key,'show',commit+':'+path)
    (OUT/name).write_bytes(b)
    catalog.append({'id':name,'repository':repo,'commit':commit,'original_path':path,
                    'sha256':sha(b),'bytes':len(b),'snapshot':'sources/'+name,
                    'evidence_type':kind,'selectors':selectors or [],
                    'url':f'https://github.com/songsiyi2006-chem/{repo}/blob/{commit}/{path}'})

for key in REPOS:
    assert git(key,'rev-parse','HEAD').decode().strip()==REPOS[key][1]
    copy(key,'LICENSE',key+'_LICENSE.txt','license')
    log=git(key,'log','--reverse','--format=%H%x09%aI%x09%s').decode('utf8')
    with (OUT/(key+'_history.csv')).open('w',encoding='utf8',newline='') as f:
        w=csv.writer(f);w.writerow(['commit','author_timestamp','subject','review_category'])
        for line in log.splitlines():
            h,t,s=line.split('\t',2)
            category=('audit_or_repair' if s.lower().startswith(('fix','audit','data(')) else
                      'organization_or_merge' if s.lower().startswith(('docs','refactor','merge')) else 'initial_or_stage_addition')
            w.writerow([h,t,s,category])
for f in ['README.md','aqsol_model.py','esol_model.py','paclitaxel.py','predict_solubility.py']:
    copy('aqueous',f,'aqueous_'+f,'report_only' if f.endswith('.md') else 'source_protocol')
for f in ['docs/EVIDENCE.md','docs/PHASE_INDEX.md','shared/audit/reports/PHASE1_19_CLOSEOUT.md',
          'shared/audit/reports/PHASE1_5_AUDIT_FIXES.md','shared/audit/reports/PHASE11_13_AUDIT_FIXES.md']:
    copy('complex',f,'complex_'+Path(f).name,'audit_document')
p11='shared/audit/evidence/phase1_19_rerun_20260910/phase11/20260910T142310/'
p01='shared/audit/evidence/phase1_19_rerun_20260910/phase01/20260910T135632/'
entries=[
 ('projects/phase01/results/benchmark_results.json','phase01_original.json','executed_forcefield', ['/results']),
 (p01+'bench_results/benchmark_results.json','phase01_audited.json','executed_forcefield', ['/results/*/conformers']),
 (p01+'acceptance_review.json','phase01_acceptance.json','acceptance_audit',[]),
 (p11+'results_phase11/phase11_results.json','phase11_results.json','executed_vmc', ['/systems/*/blocks','/systems/*/E_final','/self_tests']),
 (p11+'results_phase11/references.json','phase11_references.json','computed_reference',[]),
 (p11+'acceptance_review.json','phase11_acceptance.json','acceptance_audit',['/checks/chemical_accuracy','/checks/cusp']),
 ('projects/phase24/results/metrics.json','phase24_metrics.json','synthetic_response',[]),
 ('projects/phase24/results/measurements.json','phase24_measurements.json','synthetic_response',['/*/condition_id','/*/round','/*/yield_pct','/*/ee_abs_pct']),
 ('shared/phase25_27/results/acceptance.json','phase25_27_acceptance.json','acceptance_audit',['/criteria']),
 ('projects/phase28/results/run_20260920/published_trace_metrics.csv','phase28_trace_metrics.csv','published_measurement_reanalysis',['rows 2:4; trace_id,n_source_rows,unique_timestamps,timestamp_balanced_change']),
 ('projects/phase28/results/run_20260920/published_trace_summary.json','phase28_trace_summary.json','published_measurement_reanalysis',[]),
 ('projects/phase28/results/run_20260920/published_trace_hour_bins.csv','phase28_hour_bins.csv','published_measurement_reanalysis',[]),
 ('projects/phase28/results/additional_public_analysis/stability_window_robustness.csv','phase28_window_robustness.csv','published_measurement_reanalysis',['72 data rows; trace_id,startup_exclusion_h,window_h,change']),
 ('projects/phase28/results/additional_public_analysis/summary.json','phase28_additional_summary.json','published_measurement_reanalysis',[]),
 ('projects/phase28/data/literature/source_manifest.json','phase28_source_manifest.json','third_party_provenance',[]),
 ('projects/phase28/data/literature/extraction_manifest.json','phase28_extraction_manifest.json','third_party_provenance',[]),
 ('projects/phase28/data/additional_public/source_manifest.json','phase28_additional_source_manifest.json','third_party_provenance',[]),
 ('projects/phase28/data/additional_public/production_fe_summary_inputs.csv','phase28_production_fe_inputs.csv','published_summary_measurements',[]),
 ('projects/phase28/results/additional_public_analysis/replicate_mean_intervals.csv','phase28_replicate_mean_intervals.csv','published_measurement_reanalysis',[]),
 ('projects/phase28/results/additional_public_analysis/time_contrasts_correlation_sensitivity.csv','phase28_time_contrasts.csv','published_measurement_reanalysis',[]),
 ('projects/phase29/results/acceptance.json','phase29_acceptance.json','acceptance_audit',['/computed']),
 ('projects/phase29/results/xtb_matrix/manifest.json','phase29_xtb_manifest.json','executed_xtb_ledger',[]),
 ('projects/phase29/results/xtb_matrix/analysis.json','phase29_xtb_analysis.json','executed_xtb',[]),
 ('projects/phase29/results/xtb_matrix/state_differences.csv','phase29_state_differences.csv','executed_xtb',['72 data rows; molecule_id,gfn,environment']),
 ('projects/phase29/results/dft_crosscheck/manifest.json','phase29_dft_manifest.json','executed_dft_ledger',['/jobs']),
 ('projects/phase29/results/dft_crosscheck/summary.json','phase29_dft_summary.json','executed_dft',[]),
 ('projects/phase29/results/dft_crosscheck/energy_differences.csv','phase29_dft_differences.csv','executed_dft',['4 data rows']),
 ('projects/phase29/results/dft_crosscheck/comparison/comparison.csv','phase29_comparison.csv','matched_method_comparison',['8 data rows; geometry_sha256,basis,gfn,PBE0_delta_E_eV,GFN_delta_E_eV']),
 ('projects/phase29/results/dft_crosscheck/comparison/summary.json','phase29_comparison_summary.json','matched_method_comparison',[]),
 ('projects/phase30/results/acceptance.json','phase30_acceptance.json','uncalibrated_reduced_scenario',[]),
 ('projects/phase31/results/local_pilot/results/status.json','phase31_status.json','role_execution_ledger',[]),
 ('projects/phase31/results/local_pilot/results/pareto_top_leads.json','phase31_leads.json','acceptance_audit',[]),
]
for e in entries: copy('complex',*e)
for molecule,basis in [('Q01','def2-svp'),('Q01','def2-tzvp'),('Q02','def2-svp'),('Q03','def2-svp')]:
    for state in ['cation','radical']:
        job=f'{molecule}_{basis}_{state}'
        copy('complex',f'projects/phase29/results/dft_crosscheck/raw/{job}/result.json',
             'phase29_'+job+'.json','executed_dft_individual_result')

# Exact ancestry and shared blob checks, without interpreting repository names as replication.
root=git('complex','rev-list','--max-parents=0','HEAD').decode().strip()
assert root==REPOS['aqueous'][1]
shared=[]
for f in ['aqsol_model.py','esol_model.py','molecules.py','paclitaxel.py','predict_solubility.py']:
    a=git('aqueous','show',REPOS['aqueous'][1]+':'+f)
    cp='projects/legacy_solubility/code/'+f
    b=git('complex','show',REPOS['complex'][1]+':'+cp)
    shared.append({'aqueous_path':f,'complex_path':cp,'aqueous_sha256':sha(a),'complex_sha256':sha(b),'identical':a==b})
for f in ['aqsol_result.png','esol_result.png','paclitaxel.png','molecules.png']:
    a=git('aqueous','show',REPOS['aqueous'][1]+':figures/'+f)
    b=git('complex','show',REPOS['complex'][1]+':projects/phase01/figures/'+f)
    shared.append({'aqueous_path':'figures/'+f,'complex_path':'projects/phase01/figures/'+f,'aqueous_sha256':sha(a),'complex_sha256':sha(b),'identical':a==b})
tree=git('aqueous','ls-tree','-r','--name-only','HEAD').decode().splitlines()
dump('deduplication.json',{'common_root_commit':root,'shared_files':shared,
     'aqueous_tree':tree,'aqueous_csv_json_prediction_or_model_artifacts':
       [p for p in tree if p.endswith(('.csv','.json','.pkl','.pt','.joblib'))],
     'scope':'Exact committed bytes. No reconstruction of prediction points from raster images.'})

audit={'schema_version':'1.0','scope':'Read-only pinned history audit, no new scientific solver/training calls.',
 'repositories':[{'repository':r,'commit':c,'commit_count':int(git(k,'rev-list','--count','HEAD')),
                  'remote_verified_by':'git ls-remote HEAD with TLS validation; 2026-09-29',
                  'url':'https://github.com/songsiyi2006-chem/'+r} for k,(r,c) in REPOS.items()],
 'deduplication':'sources/deduplication.json','source_artifacts':catalog,
 'license':{'repository_material':'MIT; both original license notices preserved.',
            'third_party_data':'Phase28 source article/workbook rights remain with original authors/publisher. No article text, publisher figures or full workbooks redistributed. Selected numerical summaries and derived aggregates retained with DOI, extraction coordinates and original workbook hashes. No claim that repository MIT relicenses third-party material.',
            'source_license_verification_limit':'Publisher article metadata verified via primary web search; live rights text unavailable behind publisher identity redirect. Do not label third-party source material CC-BY without separate verification.',
            'aqsoldb':'No AqSolDB/ESOL raw measurements redistributed; only original MIT code and historical reported metrics.'},
 'historical_review_scope':{'all_commit_subjects_and_categories':True,'all_commit_diffs_or_all_scientific_outputs':False,
    'precedence':'Latest scoped acceptance and fresh-run ledgers supersede legacy figure/caption claims. Historical running statuses are not current job status.',
    'current_head_covers_phases':'01-31; phase32 not in this pinned remote history'},
 'new_work_accounting':{'quantum_energy_calls':0,'training_runs':0,'molecular_optimizations':0,'new_experiments':0,'analysis_kind':'deterministic arithmetic over frozen outputs'},
 'recompute_command':'python manuscript/history/learning/sources/recompute_history.py',
 'integrity_hashes':{}}
(OUT.parent/'audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
print(json.dumps({'copied_artifacts':len(catalog),'all_shared_files_identical':all(x['identical'] for x in shared),'commits':[x['commit_count'] for x in audit['repositories']]}))
