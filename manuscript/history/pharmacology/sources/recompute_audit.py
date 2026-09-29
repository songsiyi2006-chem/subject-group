"""Read-only arithmetic audit of frozen results; no new MD, docking, QM or fits."""
from pathlib import Path
import csv, hashlib, json, math, statistics, subprocess
import numpy as np

BASE = Path(__file__).resolve().parents[1]
SRC = BASE/'sources'
REPO = 'ai4pharm-lead-developability-suite'
ROOT = SRC/REPO
MANIFEST = json.loads((SRC/'source_manifest.json').read_text(encoding='utf-8'))
PATHS = {x['path']:x for x in MANIFEST['files'] if x['repository']==REPO}
DER = SRC/'derived'
DER.mkdir(exist_ok=True)
FACTS=[]
CHECKS=[]
def read(path):
    p=ROOT/path
    return json.loads(p.read_text(encoding='utf-8-sig')) if p.suffix=='.json' else list(csv.DictReader(p.read_text(encoding='utf-8-sig').splitlines()))
def fact(key,value,unit,path,fields,derivation,evidence):
    FACTS.append(dict(id=key,value=value,unit=unit,source=PATHS[path],fields=fields,derivation=derivation,evidence_class=evidence))
def check(key,value,tolerance):
    CHECKS.append(dict(id=key,value=float(value),tolerance=tolerance,passed=bool(abs(value)<=tolerance)))
def write_csv(name,rows):
    p=DER/name
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def lo(x,a,b):return float(np.clip((b-x)/(b-a),0,1))

def main():
    p1='projects/task01_lead_developability/results_task1/developability_results.json'
    rows=read(p1); groups=[]; mpo=[];esol=[];oral=[]
    for x in rows:
        score=lo(x['clogp_rdkit'],3,5)+lo(x['logd74_approx'],2,4)+lo(x['mw_da'],360,500)+min(1-lo(x['tpsa_a2'],20,40),lo(x['tpsa_a2'],90,120))+lo(x['hbd'],.5,3.5)+(1 if x['basic_pka_assumed'] is None else lo(x['basic_pka_assumed'],8,10))
        mpo.append(score-x['cns_mpo_approx'])
        logs=.16-.63*x['clogp_rdkit']-.0062*x['mw_da']+.066*x['rotatable_bonds']-.74*x['aromatic_atom_fraction']
        esol.append(logs-x['esol_logs_mol_l'])
        oral.append(100*statistics.mean([x[k] for k in ['oral_d_solubility','oral_d_permeability','oral_d_herg_proxy','oral_d_polarity','oral_d_fsp3','oral_d_aromatic_rings']])-x['oral_developability_score_proxy'])
    check('mpo_six_terms_max_abs',max(map(abs,mpo)),1e-12)
    check('esol_formula_max_abs',max(map(abs,esol)),1e-12)
    check('oral_equal_weights_max_abs',max(map(abs,oral)),1e-12)
    for group in ['Oral Drugs','Toxic Dropouts','bRo5 Modalities']:
        g=[x for x in rows if x['archetype']==group]
        item=dict(archetype=group,n=len(g),mpo_median=statistics.median(x['cns_mpo_approx'] for x in g),oral_proxy_median=statistics.median(x['oral_developability_score_proxy'] for x in g),embedded=sum(x['conformers_embedded'] for x in g),converged=sum(x['conformers_converged'] for x in g),partial=sum(x['geometry_status']=='partial_convergence' for x in g))
        groups.append(item)
    fact('mpo_groups',groups,'mixed; score/conformer counts',p1,['archetype','cns_mpo_approx','oral_developability_score_proxy','conformers_embedded','conformers_converged','geometry_status'],'Group medians and sums of original full-precision records.','computed descriptors and uncalibrated heuristic scores')
    fact('absolute_chameleonic_missing',sum(x['absolute_chameleonic_hb_index'] is None for x in rows),'records',p1,['absolute_chameleonic_hb_index'],'Count null values; absence is not zero.','explicitly unidentified observable')
    write_csv('mpo_plot.csv',[{k:x[k] for k in ['name','archetype','cns_mpo_approx','custom_mpo_logp_2_4','mpo_pka_sensitivity_min','mpo_pka_sensitivity_max','geometry_status','conformers_embedded','conformers_converged','esol_logs_mol_l']} for x in rows])

    p2='projects/task02_tpd/data_task2/equilibrium_species.csv'
    eq=read(p2); maxima={k:0. for k in ['E','T','P','EP','PT','EPT']}
    plot=[]
    for x in eq:
        f={k:float(v) for k,v in x.items()}; e,t,p,ep,pt,ept=[f[k] for k in ['E_free_nM','T_free_nM','P_free_nM','EP_nM','PT_nM','EPT_nM']]
        residuals=dict(E=abs(e+ep+ept-100)/100,T=abs(t+pt+ept-100)/100,P=abs(p+ep+pt+ept-f['P_total_nM'])/f['P_total_nM'],EP=abs(ep-e*p/10)/max(ep,1e-300),PT=abs(pt-t*p/100)/max(pt,1e-300),EPT=abs(ept-f['alpha']*e*t*p/1000)/max(ept,1e-300))
        for key,value in residuals.items():maxima[key]=max(maxima[key],value)
        plot.append(dict(alpha=f['alpha'],P_total_nM=f['P_total_nM'],EPT_nM=ept,mass_balance_relative_max=max(residuals[k] for k in ['E','T','P'])))
    for key,value in maxima.items():check('tpd_'+key+'_relative',value,1e-9)
    write_csv('tpd_plot.csv',plot)
    p2m='projects/task02_tpd/data_task2/equilibrium_metrics.csv'
    peaks=[{k:float(v) for k,v in x.items()} for x in read(p2m)]
    fact('tpd_peaks',peaks,'nM except alpha',p2m,list(peaks[0]),'Original peak/window rows; source parameters E0=T0=100 nM, KD_E=10 nM, KD_T=100 nM.','analytic mass-action scenario')
    fact('tpd_supplied_peak_relative_excess_percent',100*(peaks[0]['supplied_formula_p_nm']/peaks[0]['exact_peak_total_p_nm']-1),'%',p2m,['supplied_formula_p_nm','exact_peak_total_p_nm'],'100*(supplied/exact-1), same parameter scenario.','analytic formula diagnostic')
    pl='projects/task02_tpd/data_task2/linker_conformers.csv'
    lr=read(pl); links=[]
    for name in sorted({x['linker'] for x in lr}):
        g=[x for x in lr if x['linker']==name];d=np.array([float(x['distance_A']) for x in g]);s=next(x for x in read('projects/task02_tpd/data_task2/linker_summary.csv') if x['linker']==name)
        check('linker_'+name+'_mean_error',d.mean()-float(s['mean_distance_A']),1e-12)
        check('linker_'+name+'_sd_error',d.std(ddof=1)-float(s['distance_sd_A']),1e-12)
        links.append(dict(linker=name,n=len(g),mean_A=float(d.mean()),sd_A=float(d.std(ddof=1)),compatible_fraction=float(((d>=8)&(d<=12)).mean()),converged=sum(int(x['mmff_status'])==0 for x in g)))
    fact('linker_statistics',links,'angstrom and fraction',pl,['linker','distance_A','mmff_status'],'Mean, sample SD (ddof=1), and inclusion in assumed 8–12 A interval; duplicate draws retained.','unweighted force-field sampling of capped linker proxies')

    p3='projects/task03_covalent_kinetics/data_task3/warhead_results.csv'
    war=read(p3);byid={x['id']:x for x in war};wp=[]
    for x in war:
        kd=float(x['koff_s'])/float(x['kon_M_inv_s']);ki=(float(x['koff_s'])+float(x['kinact_s']))/float(x['kon_M_inv_s'])
        check(x['id']+'_KD_relative',kd/float(x['KD_M'])-1,1e-14);check(x['id']+'_KI_relative',ki/float(x['KI_QSSA_M'])-1,1e-14)
        wp.append(dict(id=x['id'],rapid_equilibrium_efficiency=float(x['kinact_over_KD_M_inv_s']),qssa_efficiency=float(x['kinact_over_KI_QSSA_M_inv_s']),relative_overestimate_percent=100*(ki/kd-1),center_fplus_proxy=float(x['center_fplus_proxy']),center_rank=int(x['center_heavy_atom_rank']),LUMO_eV=float(x['LUMO_eV']),noncovalent_residence_s=float(x['noncovalent_residence_s'])))
    write_csv('covalent_plot.csv',wp)
    fact('covalent_efficiencies',wp,'M^-1 s^-1, %, eV and s',p3,['kon_M_inv_s','koff_s','kinact_s','kinact_over_KD_M_inv_s','kinact_over_KI_QSSA_M_inv_s','center_fplus_proxy','center_heavy_atom_rank','LUMO_eV','noncovalent_residence_s'],'KD=koff/kon; KI=(koff+kinact)/kon; percentage excess = 100*kinact/koff.','EHT descriptors plus hypothetical kinetic constants; no measured rates')
    pk='projects/task03_covalent_kinetics/data_task3/kobs_concentration.csv';kr=read(pk);errors=[];saved=[]
    for x in kr:
        w=byid[x['id']];a=float(w['kon_M_inv_s'])*float(x['concentration_M']);s=a+float(w['koff_s'])+float(w['kinact_s']);slow=2*a*float(w['kinact_s'])/(s+math.sqrt(s*s-4*a*float(w['kinact_s'])))
        errors.append(abs(float(x['kobs_ODE_fit_s'])/slow-1));saved.append(abs(float(x['kobs_exact_s'])/slow-1))
    check('kobs_saved_eigenvalue_max_relative',max(saved),1e-14);check('kobs_ODE_vs_eigenvalue_max_relative',max(errors),1e-6)
    fact('kobs_diagnostic',dict(rows=len(kr),max_ODE_vs_analytic_relative=max(errors)),'count and relative error',pk,['id','concentration_M','kobs_ODE_fit_s','kobs_exact_s'],'Stable smaller eigenvalue of the forward-only kinetic matrix. Parameters from warhead_results; reverse reaction is intentionally disabled in this diagnostic.','numerical consistency of a hypothetical kinetic scenario')

    p19='projects/task19_metadynamics/outputs/summary.json';summary=read(p19)
    c19=read('projects/task19_metadynamics/outputs/config.json');cv=read('projects/task19_metadynamics/outputs/atomistic_cv_trajectory.csv');hills=read('projects/task19_metadynamics/outputs/hill_history.csv')
    height_err=max(abs(float(h['height_kJ_mol'])-.25*math.exp(-float(h['bias_before_deposition_kJ_mol'])/(.008314462618*300*9))) for h in hills)
    check('hill_height_max_abs_kJ_mol',height_err,1e-12)
    check('MD_step_time_max_abs_ps',max(abs(float(x['step'])*.001-float(x['time_ps'])) for x in cv),1e-14)
    fact('md_pilot',dict(particles=summary['atomistic']['particles'],equilibration_ps=.001*c19['equilibration_steps'],biased_ps=.001*c19['metadynamics_steps'],recorded_frames=len(cv),hills=len(hills),distance_span_nm=max(float(x['distance_nm']) for x in cv)-min(float(x['distance_nm']) for x in cv),angle_span_rad=max(float(x['orientation_rad']) for x in cv)-min(float(x['orientation_rad']) for x in cv),pmf_converged=summary['atomistic']['pmf_converged'],cooperativity_alpha=summary['atomistic']['cooperativity_alpha']),'mixed; nm, rad, ps',p19,['atomistic'], 'Durations independently checked against config and trajectory steps; spans recomputed from archived CV rows.','short all-atom explicit-solvent pilot; no converged free energy')
    ppdb='projects/task19_metadynamics/inputs/5T35.pdb';atoms=[]
    for line in (ROOT/ppdb).read_text().splitlines():
        if line.startswith('ATOM  ') and line[16] in (' ','A') and line[21] in ('A','D') and line[76:78].strip()!='H':atoms.append((line[21],line[22:27],np.array([float(line[j:j+8]) for j in [30,38,46]])))
    va=[a for a in atoms if a[0]=='D'];ta=[a for a in atoms if a[0]=='A'];dist=np.linalg.norm(np.array([a[2] for a in va])[:,None,:]-np.array([a[2] for a in ta])[None,:,:],axis=2);pairs=np.argwhere(dist<=4.5);residues={(va[i][1],ta[j][1]) for i,j in pairs}
    contacts=dict(VHL_heavy_atoms=len(va),BRD4_heavy_atoms=len(ta),atom_pairs=len(pairs),residue_pairs=len(residues),minimum_distance_A=float(dist.min()))
    check('structure_atom_pair_count',len(pairs)-summary['structure']['PPI_atom_pairs_within_4_5_A'],0)
    check('structure_min_distance_A',dist.min()-summary['structure']['closest_PPI_distance_A'],1e-12)
    fact('crystal_contacts',contacts,'counts and angstrom',ppdb,['ATOM coordinates; author chains D/A; alternate locations blank/A'],'All heavy-atom distances; cutoff <=4.5 A. These are crystal contacts, not predicted affinities.','public experimental crystallographic structure; derived geometry')

    p20='projects/task20_denovo/outputs/summary.json';s20=read(p20);gen=read('projects/task20_denovo/outputs/generation_summary.csv');random=read('projects/task20_denovo/outputs/random_search_control.csv');candidate=read('projects/task20_denovo/outputs/evaluated_candidates.csv');dock=read('projects/task20_denovo/outputs/docking/candidate_vina_scores.csv')
    check('unique_candidate_count',len({x['id'] for x in candidate})-s20['unique_graphs_including_random_control'],0)
    check('random_unique_count',len({x['id'] for x in random})-s20['random_search_control']['unique_graph_budget'],0)
    check('random_best_geometry',max(float(x['pocket_geometry_score']) for x in random)-s20['random_search_control']['best_geometry_score'],1e-12)
    bestfinal=float(gen[-1]['best_pocket_score']);bestrand=max(float(x['pocket_geometry_score']) for x in random)
    stat=dict(evolution_budget=s20['unique_graphs_evaluated'],random_budget=len(random),union_unique=len(candidate),overlap=s20['unique_graphs_evaluated']+len(random)-len(candidate),population_rows=sum(int(x['population']) for x in gen),best_initial=float(gen[0]['best_pocket_score']),best_final=bestfinal,best_random=bestrand,best_gain_vs_initial_percent=100*(bestfinal/float(gen[0]['best_pocket_score'])-1),best_gain_vs_random_percent=100*(bestfinal/bestrand-1),selected_median_final=float(gen[-1]['median_pocket_score']),selected_median_random=s20['random_search_control']['selected_median_geometry_score'],candidate_docks=len(dock),vina_min=min(float(x['vina_score_kcal_mol']) for x in dock),vina_max=max(float(x['vina_score_kcal_mol']) for x in dock),redocking_rmsd_A=[x['symmetry_aware_heavy_atom_rmsd_A_without_superposition'] for x in s20['docking']['redocking_controls']])
    fact('optimization_comparison',stat,'count, arbitrary geometry score, %, nominal kcal/mol, angstrom',p20,['population','generations','history_rows','unique_graphs_evaluated','unique_graphs_including_random_control','random_search_control','best_geometry_score_initial','best_geometry_score_final','docking'],'Budget/count/extrema checked from generation, all-candidate, random-control and docking rows. One random seed only.','executed geometric search and empirical docking; no measured binding')
    write_csv('optimization_plot.csv',[dict(generation=int(x['generation']),population=int(x['population']),evaluated_unique=int(x['evaluated_unique']),best_geometry=float(x['best_pocket_score']),median_geometry=float(x['median_pocket_score']),random_best_geometry=bestrand,random_selected_median=s20['random_search_control']['selected_median_geometry_score']) for x in gen])

    histories=[]
    clone=BASE.parents[3]/'history_sources'/REPO
    oldcases=[('aa7c77f790a285250c8b635739f4eabdc2324fc8','results_task1/developability_results.csv','projects/task01_lead_developability/results_task1/developability_results.csv'),('12915ee798cdccf12130bc7c7d315d89b5bbbf71','data_task2/equilibrium_metrics.csv','projects/task02_tpd/data_task2/equilibrium_metrics.csv'),('541e7644b84ada2199abfcc90f3477a21b52470c','data_task3/warhead_results.csv','projects/task03_covalent_kinetics/data_task3/warhead_results.csv')]
    for commit,old,new in oldcases:
        ids=[subprocess.check_output(['git','rev-parse',rev+':'+path],cwd=clone).decode().strip() for rev,path in [(commit,old),('a7a2a4df501d5dad030dc488691b02c73fdb423f',new)]]
        histories.append(dict(original_commit=commit,original_path=old,current_path=new,original_blob_sha1=ids[0],current_blob_sha1=ids[1],identical_git_blob=ids[0]==ids[1],current_sha256=PATHS[new]['sha256']))
    duplicate=[]
    for rel in [oldcases[0][2],oldcases[1][2],oldcases[2][2]]:
        fresh='data_omnibus/recompute/'+rel
        duplicate.append(dict(archived=rel,recompute=fresh,identical_bytes=PATHS[rel]['sha256']==PATHS[fresh]['sha256'],archived_sha256=PATHS[rel]['sha256'],recompute_sha256=PATHS[fresh]['sha256'],interpretation='Repeated workflow snapshots, not independent chemical replicates.'))
    derived=[dict(path=p.relative_to(BASE).as_posix(),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),rows=len(list(csv.DictReader(p.read_text().splitlines())))) for p in sorted(DER.glob('*.csv'))]
    out=dict(status='pass' if all(x['passed'] for x in CHECKS) else 'needs_review',scope='Selected pharmacology evidence plus complete git history/tree inventories of two public repositories. No scientific reruns, new fits or experimental measurements.',repositories=MANIFEST['repositories'],source_manifest='sources/source_manifest.json',facts=FACTS,independent_arithmetic_checks=CHECKS,historical_identity=histories,recomputed_snapshots=duplicate,derived_csvs=derived,excluded_scope=['Clinical efficacy/dose claims','Synthetic assay readouts treated as measured labels','Literature companions treated as executed calculations','Tasks 4–18 not selected as mechanistic evidence'],limitations=['Structural identifiers are not measured activity labels.','Hills, finite trajectories and valid coordinates do not certify PMF convergence.','A single random control cannot estimate optimizer superiority.','Current source and historical executed source hashes are distinguished in original manifests.'])
    (BASE/'audit.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':out['status'],'checks':len(CHECKS),'facts':len(FACTS),'groups':groups,'md':next(x['value'] for x in FACTS if x['id']=='md_pilot'),'optimization':stat,'history':histories,'recompute':duplicate},ensure_ascii=False))

if __name__=='__main__':main()
