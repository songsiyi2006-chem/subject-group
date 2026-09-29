"""Validate expanded-paper provenance, independent arithmetic and final documents.

The checks establish consistency of archived records and publication artifacts,
not experimental validation, novelty, journal acceptance or universal accuracy.
Standard-library only. Final document requirements can be deferred with --data-only.
"""
import argparse,csv,hashlib,json,math,re,statistics,subprocess,sys
from collections import defaultdict
from pathlib import Path
BASE=Path(__file__).resolve().parents[1];ROOT=BASE.parents[1];HISTORY=BASE.parent/'history'
checks=[]
def check(name,ok):checks.append({'name':name,'passed':bool(ok)})
def load(p):return json.loads(p.read_text('utf-8-sig'))
def rows(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def close(a,b):return math.isclose(float(a),float(b),rel_tol=2e-11,abs_tol=2e-12)
def hashcheck(name,path,digest):check(name,path.is_file() and sha(path)==digest)
def main(data_only):
    inv=load(HISTORY/'repository_inventory.json')
    check('six distinct public repositories',len(inv['repositories'])==6 and len({r['name'] for r in inv['repositories']})==6)
    for r in inv['repositories']:
        check('full frozen commit '+r['name'],bool(re.fullmatch('[0-9a-f]{40}',r['frozen_commit'])))
    analysis=load(BASE/'results/multiscale_analysis.json')
    for path,digest in analysis['input_sha256'].items():hashcheck('archived multiscale input '+path,ROOT/path,digest)
    for r in analysis['figures']:hashcheck('multiscale figure '+r['path'],ROOT/r['path'],r['sha256'])
    hashcheck('multiscale generator',BASE/'scripts/analyze_multiscale.py',analysis['generator_sha256'])
    check('posthoc work not miscounted as new solvers',all(analysis[k]==0 for k in ['new_quantum_jobs','new_training_runs','new_PDE_solves','new_experiments']))
    # Independently group the saved fixed-budget records, without calling generator.
    raw=rows(ROOT/'electratwin/results/benchmark_extension/budget_prefixes.csv')
    final=defaultdict(dict)
    for r in raw:
        if int(r['evaluation'])==25:final[r['method']][int(r['seed'])]=float(r['full_pool_HV_fraction'])
    for r in rows(BASE/'results/optimization_final.csv'):
        values=list(final[r['method']].values())
        check('64 unique seeds '+r['method'],len(values)==64)
        for key,value in [('mean_HV_fraction',statistics.mean(values)),('sd_across_seeds',statistics.stdev(values)),('min_HV_fraction',min(values)),('max_HV_fraction',max(values))]:check('budget summary '+r['method']+' '+key,close(r[key],value))
    for r in analysis['paired']:
        diffs=[final[r['left']][s]-final[r['right']][s] for s in sorted(final[r['left']])]
        check('paired mean '+r['right'],close(r['mean_fraction_difference'],statistics.mean(diffs)))
        check('paired outcomes '+r['right'],[r['wins'],r['losses'],r['ties']]==[sum(x>1e-12 for x in diffs),sum(x< -1e-12 for x in diffs),sum(abs(x)<=1e-12 for x in diffs)])
    raw=rows(ROOT/'toolkit/results/audit/deconvolution_recovery.csv');groups=defaultdict(list)
    for r in raw:groups[(r['width_mismatch'],float(r['separation_min']))].append(r)
    check('80 synthetic chromatograms',len(raw)==80 and all(r['evidence_role']=='synthetic' for r in raw))
    for r in rows(BASE/'results/chromatography_grouped.csv'):
        group=groups[(r['width_mismatch'],float(r['separation_min']))]
        check('ten chromatograms per condition '+str(list(group[0].values())[1:3]),len(group)==10)
        for saved,original in [('mean_absolute_fitted_error_pct','relative_error_pct'),('mean_absolute_window_error_pct','window_error_pct')]:check('chromatography '+str((r['width_mismatch'],r['separation_min'],saved)),close(r[saved],statistics.mean(abs(float(x[original])) for x in group)))
    for r in rows(BASE/'results/thermostat_comparison.csv'):
        expected=float(r['mean_within_trajectory_temperature_variance_K2'])/float(r['canonical_reference_variance_temperature_K2'])
        check('canonical variance ratio '+r['thermostat'],close(r['variance_ratio_to_canonical'],expected))
    for r in rows(BASE/'results/kinetic_potential_comparison.csv'):
        err=float(r['mean_TOF_s_1'])/float(r['finite_window_CTMC_TOF_s_1'])-1
        check('finite-window relative error '+r['eta_V'],close(r['relative_mean_error_vs_finite_CTMC'],err))
        check('stationary occupations '+r['eta_V'],close(sum(float(r['steady_'+s]) for s in ['empty','substrate','radical','intermediate','product']),1))
    # Learning source bytes and selected arithmetic checks.
    a=load(HISTORY/'learning/audit.json')
    for r in a['source_artifacts']:hashcheck('learning snapshot '+r['snapshot'],HISTORY/'learning'/r['snapshot'],r['sha256'])
    for p,d in a.get('integrity_hashes',{}).items():hashcheck('learning derived artifact '+p,HISTORY/'learning'/p,d)
    for r in rows(HISTORY/'learning/sources/derived_reference_alignment.csv'):
        check('reference double difference '+r['molecule_id']+'/'+r['gfn'],close(r['double_difference_eV'],float(r['GFN_Q01_centered_eV'])-float(r['PBE0_Q01_centered_eV'])))
        check('raw reference offset '+r['molecule_id']+'/'+r['gfn'],close(r['raw_GFNmPBE0_eV'],float(r['GFN_raw_delta_eV'])-float(r['PBE0_raw_delta_eV'])))
    trace=rows(HISTORY/'learning/sources/derived_trace_summary.csv')
    check('trace timestamp deduplication',sum(int(r['n_rows']) for r in trace)==196496 and sum(int(r['unique_timestamps']) for r in trace)==39726)
    for r in trace:check('trace multiplicity '+r['trace_id'],int(r['n_rows'])-int(r['unique_timestamps'])==int(r['duplicate_rows']))
    # Public pharmaceutical source snapshots retain exact original bytes.
    pm=load(HISTORY/'pharmacology/sources/source_manifest.json')
    for r in pm['files']:hashcheck('pharmacology snapshot '+r['copied_path'],HISTORY/'pharmacology/sources'/r['copied_path'],r['sha256'])
    pf=load(HISTORY/'pharmacology/figures/manifest.json')
    check('sixteen pharmaceutical figure assets',len(pf['figures'])==16)
    for r in pf['figures']:
        hashcheck('pharmacology figure '+r['file'],HISTORY/'pharmacology'/r['file'],r['sha256'])
        for p,d in r['sources'].items():hashcheck('pharmacology plotted input '+p,HISTORY/'pharmacology'/p,d)
    for group in ['catalysis','pharmacology']:
        check('completed audit '+group,(HISTORY/group/'audit.json').exists())
    catalog=load(HISTORY/'catalysis/sources/catalog.json')
    for source_id,r in catalog['source_entries'].items():
        if r.get('local_excerpt'):hashcheck('catalysis excerpt '+source_id,HISTORY/'catalysis/sources'/r['local_excerpt'],r['excerpt_sha256'])
    native=load(HISTORY/'catalysis/sources/native_energy_checks.json')
    check('ten native catalytic energy checks',len(native)==10)
    for r in native:
        check('native energy identity '+r['name'],close(r['native_final_energy_hartree'],r['expected_summary_energy_hartree']) and r['energy_absolute_difference_hartree']==0)
    for r in rows(HISTORY/'catalysis/sources/posthoc_egnn_metrics.csv'):
        check('EGNN baseline difference '+r['split'],close(r['improvement_over_equal_eV'],float(r['equal_reference_MAE_eV'])-float(r['model_MAE_eV'])))
    jobs=rows(HISTORY/'catalysis/sources/posthoc_dft_jobs.csv')
    check('molecular preflight completed/timeout counts',len(jobs)==4 and sum(r['status']=='completed' for r in jobs)==2 and sum(r['status']=='timeout' for r in jobs)==2)
    check('molecular preflight no completed anion',not any(r['charge']=='-1' and r['status']=='completed' for r in jobs))
    for group,manifest,key,pathkey,inputkey,generator in [
        ('learning','figure_manifest.json','figures','path','inputs_sha256','script_sha256'),
        ('catalysis','manifest.json','files','file','source_sha256','generator_sha256')]:
        fm=load(HISTORY/group/'figures'/manifest)
        hashcheck(group+' plotting generator',HISTORY/group/'plotting.py',fm[generator])
        for r in fm[key]:
            hashcheck(group+' figure '+r[pathkey],HISTORY/group/r[pathkey],r['sha256'])
            for p,d in r[inputkey].items():hashcheck(group+' plotted input '+p,HISTORY/group/p,d)
    assembly=load(BASE/'results/assembly.json')
    for p,d in assembly['input_sha256'].items():hashcheck('assembly input '+p,ROOT/p,d)
    hashcheck('assembly script',BASE/'scripts/assemble_manuscripts.py',assembly['script_sha256'])
    for name,r in assembly['documents'].items():
        path=BASE/(name+'.md');text=path.read_text('utf-8')
        hashcheck('assembled manuscript '+name,path,r['source_sha256'])
        check('no unfinished placeholders '+name,not any(x in text for x in ['{{TABLE','{{HISTORY','TODO','待补充']))
        check('substantial source-linked figures '+name,r['figures']>=20)
        check('substantial tables '+name,r['tables']>=20)
        check('authorship '+name,('宋思毅' in text and '广西师范大学' in text) if 'chinese' in name else ('Siyi Song' in text and 'Guangxi Normal University' in text))
    en=assembly['documents']['manuscript_english'];cn=assembly['documents']['manuscript_chinese']
    check('bilingual figure/table inventory',en['figures']==cn['figures'] and en['tables']==cn['tables'])
    parity=load(BASE/'results/numeric_parity.json')
    check('bilingual table numeric tokens',parity['passed'])
    for lang,digest in parity['source_sha256'].items():hashcheck('numeric parity source '+lang,BASE/f'manuscript_{lang}.md',digest)
    run=subprocess.run([sys.executable,str(BASE/'scripts/finalize_references.py'),'--check'],capture_output=True,text=True,encoding='utf-8')
    check('verified references and bilingual citation parity',run.returncode==0)
    if run.returncode:print(run.stdout,run.stderr)
    if not data_only:
        qa=load(BASE/'results/document_qa.json')
        check('document machine checks',qa['passed'])
        visual=load(BASE/'results/visual_review.json')
        for p,d in visual['review_record_sha256'].items():hashcheck('visual review record '+p,BASE/'results'/p,d)
        repair=load(BASE/'results/visual_repair.json')
        before=load(BASE/'results/document_qa_before_kinetics_legend_repair.json')
        for old,current in zip(before['documents'],qa['documents']):
            changed=[b['page'] for a,b in zip(old['pages'],current['pages']) if a['image_sha256']!=b['image_sha256']]
            resolved=repair['documents'][current['name']]
            check('visual repair changed-page identity '+current['name'],changed==resolved['changed_pages'])
            check('visual repair final PDF '+current['name'],resolved['pdf_sha256']==current['pdf_sha256'])
            check('visual repair checked all changes '+current['name'],all(x['passed'] for x in resolved['changed_page_review']) and [x['page'] for x in resolved['changed_page_review']]==changed)
        for r in qa['documents']:
            name=r['name']
            check('30 pages before references '+name,r['reference_start_page'] is not None and r['reference_start_page']>=31)
            for ext in ['docx','pdf']:hashcheck('final '+name+'.'+ext,BASE/'documents'/(name+'.'+ext),r[ext+'_sha256'])
            reviewed=visual['documents'][name]
            check('all pages visually reviewed '+name,reviewed['pages_reviewed']==list(range(1,r['page_count']+1)))
            check('visual review matches PDF '+name,reviewed['pdf_sha256']==r['pdf_sha256'])
            check('no unresolved layout defects '+name,reviewed['unresolved_defects']==[])
    result={'passed':all(x['passed'] for x in checks),'checks':len(checks),'failed':[x for x in checks if not x['passed']],'scope':'Archived-data arithmetic, provenance and artifact integrity; not chemical validation or human peer review.','data_only':data_only,'items':checks}
    (BASE/'results/validation.json').write_bytes((json.dumps(result,ensure_ascii=False,indent=2)+'\n').encode())
    print(json.dumps({k:result[k] for k in ['passed','checks','failed']},ensure_ascii=False,indent=2))
    return 0 if result['passed'] else 1
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--data-only',action='store_true');args=ap.parse_args()
    raise SystemExit(main(args.data_only))
