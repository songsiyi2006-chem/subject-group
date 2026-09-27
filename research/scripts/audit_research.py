"""Additional model audits and actual molecular calculations; no measured labels."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import ast,hashlib,json,platform,shutil,sys,time,warnings
from pathlib import Path
import numpy as np
import pandas as pd
import scipy,sklearn
from scipy.linalg import solve_banded
from scipy.optimize import brentq
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel,Matern
from sklearn.dummy import DummyRegressor
from sklearn.preprocessing import MinMaxScaler
from rdkit import Chem,rdBase
from rdkit.Chem import AllChem,Descriptors,rdMolDescriptors,Draw
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/audit'
OBJECTIVES=['Yield_pct','FE_pct','SEC_kWh_kg']
IDS=['melatonin','caffeine','phenylquinoline','tryptophol','thiophene_ethanol','ethyl_benzofuran_carboxylate','indoline','carbazole']

def substrates():
    tree=ast.parse((ROOT/'source/run_research_engine_original.py').read_text(encoding='utf-8'))
    for node in ast.walk(tree):
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='substrates' for t in node.targets):return ast.literal_eval(node.value)
    raise ValueError('Substrate source list absent')

def make_grid(area_cm2=1.):
    rows=[]
    for j in np.linspace(5,35,15):
        for cat in [2.,5.,8.]:
            for rpm in [300,600,900]:
                voltage=2.1+.045*np.log(j+1)+j*.001*area_cm2*12
                rate=j/20*(1-np.exp(-rpm/350))*(cat/5)**.4
                y=np.clip(94*rate/(1+rate)-.08*(j-18)**2,10,96.5)
                fe=np.clip(92-1.2*j+.015*(rpm/10),15,94)
                sec=2*96485.33*voltage/(3.6e6*.223*(fe/100))
                rows.append(dict(j_mA_cm2=float(j),catalyst_mol_pct=cat,stirring_rpm=rpm,U_cell_V=float(voltage),Yield_pct=float(y),FE_pct=float(fe),SEC_kWh_kg=float(sec)))
    return pd.DataFrame(rows)

def nondominated(values):
    costs=np.asarray(values)*np.array([-1,-1,1])
    # Row i dominates row j if it is no worse in all costs and better in at least one.
    dominates=(costs[:,None,:]<=costs[None,:,:]).all(2)&(costs[:,None,:]<costs[None,:,:]).any(2)
    return ~dominates.any(0)

def audit_pareto():
    raw=make_grid();rounded=raw.round({'j_mA_cm2':2,'U_cell_V':3,'Yield_pct':2,'FE_pct':2,'SEC_kWh_kg':3})
    for df,name in [(raw,'unrounded'),(rounded,'source_rounded')]:
        df['pareto']=nondominated(df[OBJECTIVES]);df['Green_Score']=df.Yield_pct*df.FE_pct/df.SEC_kWh_kg
        df.to_csv(OUT/f'conditions_{name}.csv',index_label='condition_id')
    front=rounded[rounded.pareto].sort_values('Yield_pct',ascending=False);front.to_csv(OUT/'pareto_complete.csv',index_label='condition_id')
    best=int(rounded.loc[rounded.pareto,'Green_Score'].idxmax());best_raw=int(raw.loc[raw.pareto,'Green_Score'].idxmax())
    weights=[]
    for a in range(21):
        for b in range(21-a):
            w=np.array([a,b,20-a-b])/20
            values=w[0]*np.log(raw.Yield_pct/100)+w[1]*np.log(raw.FE_pct/100)-w[2]*np.log(raw.SEC_kWh_kg)
            # Tie breaking is lowest condition id, retained explicitly.
            idx=int(values[raw.pareto].idxmax());weights.append(dict(w_yield=w[0],w_fe=w[1],w_sec=w[2],condition_id=idx))
    wdf=pd.DataFrame(weights);wdf.to_csv(OUT/'scalarization_weights.csv',index=False)
    choices=wdf.condition_id.value_counts().rename_axis('condition_id').reset_index(name='count').merge(raw.reset_index(names='condition_id'),on='condition_id');choices.to_csv(OUT/'scalarization_choices.csv',index=False)
    thresholds=[]
    for y,fe,sec in [(40,70,.85),(45,70,.85),(50,70,.85)]:
        valid=(raw.Yield_pct>=y)&(raw.FE_pct>=fe)&(raw.SEC_kWh_kg<=sec)
        thresholds.append(dict(min_yield=y,min_FE=fe,max_SEC=sec,feasible_count=int(valid.sum()),best_condition_id=int(raw.loc[valid,'Green_Score'].idxmax()) if valid.any() else None))
    pd.DataFrame(thresholds).to_csv(OUT/'constraint_scenarios.csv',index=False)
    area=[]
    for a in [1.,5.,10.]:
        df=make_grid(a);score=df.Yield_pct*df.FE_pct/df.SEC_kWh_kg;idx=int(score.idxmax())
        area.append(dict(area_cm2=a,condition_id=idx,**df.loc[idx].to_dict()))
    pd.DataFrame(area).to_csv(OUT/'area_sensitivity.csv',index=False)
    # A real, explicitly added multi-output GP audit; the original never fits a GP.
    x=raw[['j_mA_cm2','catalyst_mol_pct','stirring_rpm']].to_numpy();y=raw[OBJECTIVES].to_numpy()
    groups=np.repeat(np.arange(15)//3,9);metrics=[];predictions=[];warn=[]
    for fold in range(5):
        train=groups!=fold;test=~train;scaler=MinMaxScaler().fit(x[train]);xt=scaler.transform(x[train]);xv=scaler.transform(x[test])
        models=[('multioutput_gp',GaussianProcessRegressor(kernel=ConstantKernel(1,(.01,100))*Matern([1,1,1],(.02,100),nu=2.5),alpha=1e-6,normalize_y=True,n_restarts_optimizer=1,random_state=42)),('training_mean',DummyRegressor())]
        for name,model in models:
            with warnings.catch_warnings(record=True) as captured:
                warnings.simplefilter('always');model.fit(xt,y[train]);pred=model.predict(xv)
            warn.extend(dict(fold=fold,model=name,message=str(v.message)) for v in captured)
            for k,obj in enumerate(OBJECTIVES):metrics.append(dict(fold=fold,model=name,objective=obj,n_train=int(train.sum()),n_test=int(test.sum()),MAE=float(np.abs(pred[:,k]-y[test,k]).mean()),RMSE=float(np.sqrt(np.mean((pred[:,k]-y[test,k])**2)))))
            for i,idx in enumerate(np.where(test)[0]):predictions.append(dict(fold=fold,model=name,condition_id=int(idx),**{f'true_{o}':float(y[idx,k]) for k,o in enumerate(OBJECTIVES)},**{f'pred_{o}':float(pred[i,k]) for k,o in enumerate(OBJECTIVES)}))
    m=pd.DataFrame(metrics);p=pd.DataFrame(predictions);m.to_csv(OUT/'gp_group_metrics.csv',index=False);p.to_csv(OUT/'gp_predictions.csv',index=False)
    agg=m.groupby(['model','objective'])[['MAE','RMSE']].mean().reset_index();agg.to_csv(OUT/'gp_summary.csv',index=False)
    (OUT/'gp_warnings.json').write_text(json.dumps(warn,indent=2)+'\n',encoding='utf-8')
    gp=p[p.model=='multioutput_gp'].sort_values('condition_id');pred_front=set(gp.loc[nondominated(gp[['pred_'+o for o in OBJECTIVES]]),'condition_id']);true_front=set(raw.index[raw.pareto]);intersection=len(pred_front&true_front)
    source=ast.parse((ROOT/'source/run_research_engine_original.py').read_text(encoding='utf-8'))
    fit_calls=sum(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='fit' for n in ast.walk(source))
    return dict(evidence='synthetic response surface; added GP is retrospective surrogate validation, not Bayesian optimization',conditions=len(raw),rounded_front=len(front),unrounded_front=int(raw.pareto.sum()),front_membership_changed=sorted(set(rounded.index[rounded.pareto])^true_front),best_condition_id=best,unrounded_best_id=best_raw,best=rounded.loc[best].to_dict(),weight_cases=len(wdf),selected_conditions=len(choices),source_choice_weight_count=int((wdf.condition_id==best).sum()),source_fit_calls=fit_calls,gp_macro_metrics=agg.to_dict('records'),gp_out_of_fold_pareto_precision=intersection/len(pred_front),gp_out_of_fold_pareto_recall=intersection/len(true_front),gp_predicted_pareto_count=len(pred_front),gp_warnings=len(warn),constraint_scenarios=thresholds,area_sensitivity=area)

def audit_scope():
    original=json.loads((ROOT/'results/original/research_grade_results.json').read_text(encoding='utf-8'))['Module2_Substrate_Scope']['scope_summary_table']
    rows=[];sites=[];confrows=[];drawmols=[];xtbrows=[];charge_rows=[]
    sys.path.insert(0,str(ROOT.parent/'scripts'));from molecular_descriptors import run_xtb,HARTREE_EV
    executable=os.environ.get('XTB_EXE') or shutil.which('xtb')
    if not executable:raise RuntimeError('Existing xTB executable required; set XTB_EXE')
    for sid,sub,orig in zip(IDS,substrates(),original):
        print('Molecule '+sid,flush=True);folder=ROOT/'results/molecules'/sid;folder.mkdir(parents=True,exist_ok=True)
        base=Chem.MolFromSmiles(sub['smiles']);assert base is not None
        mol=Chem.AddHs(base);params=AllChem.ETKDGv3();params.randomSeed=42;params.numThreads=1
        embed=AllChem.EmbedMolecule(mol,params);assert embed==0
        status=AllChem.MMFFOptimizeMolecule(mol,maxIters=300);Chem.MolToMolFile(mol,str(folder/'source_geometry.mol'))
        AllChem.ComputeGasteigerCharges(mol);qs=[float(a.GetProp('_GasteigerCharge')) for a in mol.GetAtoms()]
        graph=Chem.AddHs(base);AllChem.ComputeGasteigerCharges(graph)
        qdiff=max(abs(float(a.GetProp('_GasteigerCharge'))-qs[a.GetIdx()]) for a in graph.GetAtoms())
        candidates=[]
        for a in mol.GetAtoms():
            hs=sum(n.GetAtomicNum()==1 for n in a.GetNeighbors())
            if a.GetAtomicNum()==6 and hs:
                entry=dict(molecule=sid,atom_index=a.GetIdx(),aromatic=a.GetIsAromatic(),hybridization=str(a.GetHybridization()),attached_H=hs,gasteiger_charge=qs[a.GetIdx()]);sites.append(entry);candidates.append(entry)
        primary=min(candidates,key=lambda a:a['gasteiger_charge']);tpsa=Descriptors.TPSA(mol);logp=Descriptors.MolLogP(mol);score=1.35+.12*logp-.008*tpsa+1.5*primary['gasteiger_charge'];clipped=float(np.clip(score,.75,2.3))
        row=dict(molecule=sid,name_en=sub['name'].split(' (')[0],name_source=sub['name'],input_smiles=sub['smiles'],canonical_smiles=Chem.MolToSmiles(base),formula=rdMolDescriptors.CalcMolFormula(base),charge=Chem.GetFormalCharge(base),multiplicity=1,heavy_atoms=base.GetNumAtoms(),atoms=mol.GetNumAtoms(),MW=Descriptors.MolWt(mol),TPSA=tpsa,logP=logp,embed_status=embed,source_mmff_status=status,primary_atom=primary['atom_index'],primary_aromatic=primary['aromatic'],primary_hybridization=primary['hybridization'],primary_q=primary['gasteiger_charge'],unclipped_proxy=score,proxy=clipped,no_geometry_charge_difference=qdiff,source_random_yield=orig['predicted_electrochemical_yield_pct'])
        assert round(clipped,3)==orig['estimated_E_ox_V_vs_SCE'];assert primary['atom_index']==orig['most_reactive_carbon_idx']
        ensemble=Chem.AddHs(base);ep=AllChem.ETKDGv3();ep.randomSeed=271828;ep.numThreads=1
        ids=list(AllChem.EmbedMultipleConfs(ensemble,numConfs=8,params=ep));optim=AllChem.MMFFOptimizeMoleculeConfs(ensemble,numThreads=1,maxIters=1000)
        for cid,(s,e) in zip(ids,optim):confrows.append(dict(molecule=sid,conformer_id=int(cid),status=s,MMFF_energy_kcal_mol=e))
        good=[(float(e),int(cid)) for cid,(s,e) in zip(ids,optim) if s==0];assert good
        energy,cid=min(good);writer=Chem.SDWriter(str(folder/'conformers.sdf'))
        for c in ids:writer.write(ensemble,confId=int(c))
        writer.close();row.update(conformers=len(ids),converged=len(good),selected_conformer=cid,selected_MMFF_energy=energy);rows.append(row)
        identity=dict(row,atom_order=[dict(index_0based=a.GetIdx(),element=a.GetSymbol(),neighbors=[n.GetIdx() for n in a.GetNeighbors()]) for a in ensemble.GetAtoms()])
        (folder/'identity.json').write_text(json.dumps(identity,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
        neutral=folder/'neutral';neutral.mkdir(exist_ok=True);(neutral/'input.xyz').write_text(Chem.MolToXYZBlock(ensemble,confId=cid),encoding='utf-8')
        nr=run_xtb(executable,neutral,['input.xyz','--gfn','2','--alpb','acetonitrile','--chrg','0','--uhf','0','--opt','tight','--json'],180)
        if not (nr['completed'] and nr.get('optimization_converged')):raise RuntimeError('Neutral job incomplete: '+sid)
        cation=folder/'cation';cation.mkdir(exist_ok=True);shutil.copyfile(neutral/'xtbopt.xyz',cation/'input.xyz')
        cr=run_xtb(executable,cation,['input.xyz','--gfn','2','--alpb','acetonitrile','--chrg','1','--uhf','1','--scc','--json'],180)
        if not cr['completed']:raise RuntimeError('Cation job incomplete: '+sid)
        q0=np.loadtxt(neutral/'charges');qp=np.loadtxt(cation/'charges');response=qp-q0
        for a in ensemble.GetAtoms():
            idx=a.GetIdx();charge_rows.append(dict(molecule=sid,atom_index=idx,element=a.GetSymbol(),neutral_charge=q0[idx],cation_charge=qp[idx],removal_response=response[idx]))
        xtbrows.append(dict(molecule=sid,neutral_Eh=nr['energy_Eh'],cation_Eh=cr['energy_Eh'],removal_difference_eV=(cr['energy_Eh']-nr['energy_Eh'])*HARTREE_EV,neutral_charge_sum=float(q0.sum()),cation_charge_sum=float(qp.sum()),source_site_removal_response=float(response[primary['atom_index']])))
        depiction=Chem.Mol(base);AllChem.Compute2DCoords(depiction)
        for a in depiction.GetAtoms():a.SetProp('atomNote',str(a.GetIdx()))
        drawmols.append(depiction)
        pd.DataFrame(xtbrows).to_csv(OUT/'scope_xtb.csv',index=False)
    pd.DataFrame(rows).to_csv(OUT/'scope_descriptors.csv',index=False);pd.DataFrame(sites).to_csv(OUT/'scope_ch_sites.csv',index=False);pd.DataFrame(confrows).to_csv(OUT/'scope_conformers.csv',index=False);pd.DataFrame(charge_rows).to_csv(OUT/'scope_xtb_atom_charges.csv',index=False)
    rngrows=[]
    for seed in range(1000):
        rng=np.random.RandomState(seed)
        for row in rows:
            mean,sd=(88,3) if row['proxy']<1.45 else ((72,4) if row['proxy']<1.85 else (42,6))
            rngrows.append(dict(seed=seed,molecule=row['molecule'],random_yield=round(float(rng.normal(mean,sd)),1)))
    random=pd.DataFrame(rngrows);random.to_csv(OUT/'scope_random_yield_seeds.csv',index=False)
    stats=random.groupby('molecule').random_yield.agg(['mean','std','min','max']).reset_index();stats.to_csv(OUT/'scope_random_yield_summary.csv',index=False)
    legends=[r['name_en']+' | atom '+str(r['primary_atom']) for r in rows]
    highlights=[[r['primary_atom']] for r in rows]
    Draw.MolsToGridImage(drawmols,molsPerRow=2,subImgSize=(600,360),legends=legends,highlightAtomLists=highlights).save(ROOT/'reports/figures/substrate_atom_maps.png')
    svg=Draw.MolsToGridImage(drawmols,molsPerRow=2,subImgSize=(600,360),legends=legends,highlightAtomLists=highlights,useSVG=True)
    (ROOT/'reports/figures/substrate_atom_maps.svg').write_text(svg,encoding='utf-8')
    return dict(evidence='actual geometry/descriptors and semiempirical jobs; source proxy and randomized yields remain unvalidated',molecules=len(rows),source_mmff_nonconverged=[r['molecule'] for r in rows if r['source_mmff_status']!=0],conformers=len(confrows),converged=sum(r['status']==0 for r in confrows),geometry_charge_max_difference=max(r['no_geometry_charge_difference'] for r in rows),nonaromatic_source_sites=[r['molecule'] for r in rows if not r['primary_aromatic']],clipped_proxy_count=sum(abs(r['proxy']-r['unclipped_proxy'])>1e-12 for r in rows),xtb_jobs=2*len(xtbrows),rows=rows,xtb=xtbrows,random_yield_seed_cases=1000,random_yield_summary=stats.to_dict('records'))

def eta_sphere(phi):
    x=np.asarray(phi,dtype=float);small=x<.01
    safe=np.maximum(x,.01);normal=3*(safe/np.tanh(safe)-1)/safe**2
    val=np.where(small,1-x*x/15+2*x**4/315-x**6/1575,normal)
    return float(val) if val.ndim==0 else val

def sphere_fvm(phi,n):
    faces=np.linspace(0,1,n+1);centers=(faces[:-1]+faces[1:])/2;dr=1/n
    volume=np.diff(faces**3)/3;internal=faces[1:-1]**2/dr
    diag=phi**2*volume;diag[:-1]+=internal;diag[1:]+=internal;surface=2/dr;diag[-1]+=surface
    band=np.zeros((3,n));band[1]=diag;band[0,1:]=-internal;band[2,:-1]=-internal
    rhs=np.zeros(n);rhs[-1]=surface;c=solve_banded((1,1),band,rhs)
    return 3*float(c@volume),centers,c

def audit_transport():
    kv=4.2e-4*850;data=[];grids=[];external=[];sensitivity=[]
    for name,d in [('Microporous_POP',1.5e-10),('Hierarchical_POP',8.5e-9)]:
        for radius in [10,25,50,100,200]:
            r=radius*1e-6;phi=r*np.sqrt(kv/d);eta=eta_sphere(phi);eta_n,x,c=sphere_fvm(phi,320)
            exact=np.sinh(phi*x)/(x*np.sinh(phi));center=phi/np.sinh(phi)
            data.append(dict(regime=name,radius_um=radius,D_eff=d,kv_s=kv,phi=phi,eta=eta,eta_fvm=eta_n,eta_abs_error=abs(eta_n-eta),center_C_over_Cs=center,surface_flux_dimensionless=phi/np.tanh(phi)-1))
            pd.DataFrame(dict(r_over_R=x,c_fvm=c,c_analytic=exact)).to_csv(OUT/f'profile_{name}_{radius}.csv',index=False)
            for n in [40,80,160,320]:
                value,_,_=sphere_fvm(phi,n);grids.append(dict(regime=name,radius_um=radius,n_cells=n,eta_fvm=value,eta_exact=eta,abs_error=abs(value-eta)))
            for bi in [.1,1,10,100,1000]:
                fraction=1/(1+eta*phi**2/(3*bi));external.append(dict(regime=name,radius_um=radius,Bi=bi,Cs_over_Cbulk=fraction,eta_overall=eta*fraction))
            for dfactor,kfactor in [(1,1),(.5,1),(2,1),(1,.5),(1,2),(1,.4)]:sensitivity.append(dict(regime=name,radius_um=radius,D_factor=dfactor,kv_factor=kfactor,eta=eta_sphere(r*np.sqrt(kv*kfactor/(d*dfactor)))))
    pd.DataFrame(data).to_csv(OUT/'transport_audit.csv',index=False);pd.DataFrame(grids).to_csv(OUT/'transport_grid_convergence.csv',index=False);pd.DataFrame(external).to_csv(OUT/'external_film_scenarios.csv',index=False);pd.DataFrame(sensitivity).to_csv(OUT/'transport_parameter_sensitivity.csv',index=False)
    thresholds=[]
    for target in [.4,.85]:
        phi=brentq(lambda x:eta_sphere(x)-target,.001,100)
        for name,d in [('Microporous_POP',1.5e-10),('Hierarchical_POP',8.5e-9)]:thresholds.append(dict(regime=name,target_eta=target,phi=phi,radius_um=phi*np.sqrt(d/kv)*1e6))
    pd.DataFrame(thresholds).to_csv(OUT/'transport_thresholds.csv',index=False)
    phis=np.logspace(-6,4,301);stable=eta_sphere(phis);source=np.clip(np.where(phis<1e-3,1,3/phis**2*(phis/np.tanh(phis)-1)),.01,1)
    pd.DataFrame(dict(phi=phis,eta_stable=stable,eta_source_clipped=source)).to_csv(OUT/'transport_asymptotes.csv',index=False)
    return dict(evidence='first-order spherical diffusion/reaction model with assigned coefficients',kv_s=kv,rows=data,thresholds=thresholds,max_fvm_error=max(r['eta_abs_error'] for r in data),grid_cases=len(grids),external_film_cases=len(external),parameter_cases=len(sensitivity),source_conclusion_false=True,source_floor_at_phi_10000=float(source[-1]),unclipped_eta_at_phi_10000=float(stable[-1]))

def main():
    OUT.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    summary={'environment':dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,sklearn=sklearn.__version__,rdkit=rdBase.rdkitVersion,threads=1)}
    for key,fn in [('multiobjective',audit_pareto),('scope',audit_scope),('transport',audit_transport)]:
        print('Starting '+key,flush=True);summary[key]=fn();(OUT/'audit_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    summary['elapsed_seconds']=time.perf_counter()-start
    (OUT/'audit_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf-8');print('Research audit complete',flush=True)
if __name__=='__main__':main()
