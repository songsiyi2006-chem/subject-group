"""Independent quantitative audit. Synthetic targets remain explicitly synthetic."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import argparse,importlib.util,json,sys,warnings,time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern,WhiteKernel,ConstantKernel
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.dummy import DummyRegressor
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error,mean_absolute_error,r2_score
from sklearn.exceptions import ConvergenceWarning
from rdkit import Chem,rdBase
from rdkit.Chem import AllChem,rdFreeSASA,rdMolDescriptors,Draw
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/audit';OUT.mkdir(exist_ok=True)
SMILES='COc1ccc2[nH]c(cc2c1)c3ccccc3'
CORE_LABELS={6:'N1',7:'C2',8:'C3',9:'C3a',10:'C4',2:'C5',3:'C6',4:'C7',5:'C7a'}
FEATURES=['current_mA_cm2','temp_K','feat_eps','feat_visc','feat_DN','feat_radius','feat_E_ox']

def latent_yield(df):
    return 88-.15*(df.current_mA_cm2.to_numpy()-14)**2-.02*(df.temp_K.to_numpy()-303.15)**2+.15*(.8*df.feat_eps.to_numpy()-12*df.feat_visc.to_numpy()+.3*df.feat_DN.to_numpy())-45*(df.feat_E_ox.to_numpy()<1.8)
def bootstrap(v):
    v=np.asarray(v);rng=np.random.default_rng(618)
    vals=v[rng.integers(0,len(v),(10000,len(v)))].mean(axis=1)
    return dict(mean=float(v.mean()),ci95=np.quantile(vals,[.025,.975]).tolist(),pairs=len(v),resamples=10000)
def audit_bo(seeds):
    pool=pd.read_csv(ROOT/'results/repaired/task_a_candidates.csv');latent=latent_yield(pool)
    raw=pool[FEATURES].to_numpy();physical=(raw-raw.min(axis=0))/np.ptp(raw,axis=0)
    onehot=np.column_stack([physical[:,:2],pd.get_dummies(pool[['solvent','electrolyte']],dtype=float).to_numpy()])
    rows=[];observations=[];warning_count=0
    for seed in range(seeds):
        initial=np.random.default_rng(9000+seed).choice(len(pool),6,replace=False)
        noise=np.random.default_rng(10000+seed).normal(0,1.2,14)
        random_order=np.random.default_rng(11000+seed).permutation(len(pool))
        for name in ['random','source_gp','scaled_physical_gp','scaled_onehot_gp']:
            chosen=list(initial);y=list(np.clip(latent[initial]+noise[:6],2,98))
            for step in range(8):
                if name=='random':idx=next(int(x) for x in random_order if x not in chosen)
                else:
                    if name=='source_gp':
                        x=raw;target=y;kernel=Matern(length_scale=np.ones(raw.shape[1]),nu=2.5)+WhiteKernel(noise_level=1.)
                        gp=GaussianProcessRegressor(kernel=kernel,n_restarts_optimizer=3,random_state=42)
                    else:
                        x=physical if name=='scaled_physical_gp' else onehot
                        scale=max(float(np.std(y)),1.);target=(np.asarray(y)-np.mean(y))/scale
                        kernel=ConstantKernel(1,(.01,100))*Matern(length_scale=np.ones(x.shape[1]),length_scale_bounds=(.03,30),nu=2.5)
                        gp=GaussianProcessRegressor(kernel=kernel,alpha=1.44/scale**2,n_restarts_optimizer=1,random_state=42)
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter('always',ConvergenceWarning);gp.fit(x[chosen],target)
                    warning_count+=sum(isinstance(w.message,ConvergenceWarning) for w in caught)
                    mu,sd=gp.predict(x,return_std=True);acq=mu+2*sd;acq[chosen]=-np.inf;idx=int(acq.argmax())
                chosen.append(idx);y.append(float(np.clip(latent[idx]+noise[6+step],2,98)))
            rec=chosen[int(np.argmax(y))]
            rows.append(dict(seed=seed,method=name,recommendation_regret=float(latent.max()-latent[rec]),simple_regret=float(latent.max()-latent[chosen].max()),selected_noise_bias=float(max(y)-latent[rec])))
            observations.extend(dict(seed=seed,method=name,step=i+1,pool_index=int(idx),latent=float(latent[idx]),observed=float(y[i])) for i,idx in enumerate(chosen))
        if (seed+1)%5==0:print(f'BO audit {seed+1}/{seeds}',flush=True)
    df=pd.DataFrame(rows);df.to_csv(OUT/'task_a_paired_runs.csv',index=False)
    pd.DataFrame(observations).to_csv(OUT/'task_a_observations.csv',index=False)
    stats=df.groupby('method').recommendation_regret.agg(['mean','std','median']).reset_index();stats.to_csv(OUT/'task_a_summary.csv',index=False)
    pivot=df.pivot(index='seed',columns='method',values='recommendation_regret')
    best=pool[np.isclose(latent,latent.max())].copy();best['noiseless_yield']=latent.max();best.to_csv(OUT/'task_a_tied_grid_optima.csv',index=False)
    sourceobs=pd.read_csv(ROOT/'results/repaired/task_a_observations.csv');winner=sourceobs.iloc[sourceobs.observed_yield.argmax()]
    point=latent_yield(pd.DataFrame([winner]))[0]
    return dict(evidence='paired synthetic comparison; no measured reaction yields',seeds=seeds,methods=stats.to_dict('records'),pool_size=len(pool),source_evaluations=len(sourceobs),source_noiseless_at_winner=point,source_observed_minus_noiseless=float(winner.observed_yield-point),noiseless_grid_max=float(latent.max()),number_tied_optima=len(best),continuous_formula_max=92.4865,paired={f'{m}_minus_random':bootstrap(pivot[m]-pivot.random) for m in pivot.columns if m!='random'}|{'physical_minus_onehot':bootstrap(pivot.scaled_physical_gp-pivot.scaled_onehot_gp)},kernel_warnings=warning_count,property_provenance='All solvent/electrolyte dictionary numbers are source-assigned and lack measurement conditions/citations; not independently validated.')

def make_molecule():
    mol=Chem.AddHs(Chem.MolFromSmiles(SMILES));AllChem.ComputeGasteigerCharges(mol);return mol
def molecular_ensemble():
    mol=make_molecule();params=AllChem.ETKDGv3();params.randomSeed=20260927;params.numThreads=1
    ids=list(AllChem.EmbedMultipleConfs(mol,numConfs=32,params=params));opt=AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=1,maxIters=1000)
    records=[];conf_records=[];atom_records=[];radii=[Chem.GetPeriodicTable().GetRvdw(a.GetAtomicNum()) for a in mol.GetAtoms()]
    rawradii=rdFreeSASA.classifyAtoms(mol)
    opts=rdFreeSASA.SASAOpts();opts.algorithm=rdFreeSASA.LeeRichards;opts.probeRadius=1.4
    for cid,(status,energy) in zip(ids,opt):
        if status!=0:
            conf_records.append(dict(conformer_id=cid,mmff_status=status,mmff_energy_kcal_mol=energy));continue
        total=rdFreeSASA.CalcSASA(mol,radii,confIdx=cid,opts=opts)
        atom_sum=sum(a.GetDoubleProp('SASA') for a in mol.GetAtoms())
        conf_records.append(dict(conformer_id=cid,mmff_status=status,mmff_energy_kcal_mol=energy,total_sasa_A2=total,atomic_sasa_sum=atom_sum))
        for atom in mol.GetAtoms():
            idx=atom.GetIdx();hs=[n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum()==1]
            atom_records.append(dict(conformer_id=cid,atom_index=idx,element=atom.GetSymbol(),sasa_A2=atom.GetDoubleProp('SASA'),gasteiger_charge=float(atom.GetProp('_GasteigerCharge'))))
            if atom.GetAtomicNum()!=6 or not hs:continue
            h_sasa=sum(mol.GetAtomWithIdx(h).GetDoubleProp('SASA') for h in hs)
            q=float(atom.GetProp('_GasteigerCharge'));aromatic=atom.GetIsAromatic()
            # Per-H average makes multi-H methyl bookkeeping explicit; ranking below is aromatic-only.
            score=-.6*q+.4*(h_sasa/len(hs))/15
            records.append(dict(conformer_id=cid,atom_index=idx,site_label=CORE_LABELS.get(idx,'methoxy_methyl' if idx==0 else f'phenyl_atom_{idx}'),aromatic=aromatic,attached_H=len(hs),gasteiger_charge=q,h_sasa_sum_A2=h_sasa,h_sasa_mean_A2=h_sasa/len(hs),heuristic_score=score))
    conf=pd.DataFrame(conf_records);conf.to_csv(OUT/'task_b_conformers.csv',index=False)
    sites=pd.DataFrame(records);sites.to_csv(OUT/'task_b_conformer_sites.csv',index=False)
    pd.DataFrame(atom_records).to_csv(OUT/'task_b_conformer_atoms.csv',index=False)
    good=conf[conf.mmff_status==0];best=int(good.loc[good.mmff_energy_kcal_mol.idxmin()].conformer_id)
    (OUT/'target_lowest_mmff.xyz').write_text(Chem.MolToXYZBlock(mol,confId=best),encoding='utf-8')
    writer=Chem.SDWriter(str(OUT/'target_conformers.sdf'))
    for cid in good.conformer_id:writer.write(mol,confId=int(cid))
    writer.close()
    identity=dict(name='5-methoxy-2-phenyl-1H-indole',input_smiles=SMILES,canonical_smiles=Chem.MolToSmiles(Chem.RemoveHs(mol)),formula=rdMolDescriptors.CalcMolFormula(mol),heavy_atoms=mol.GetNumHeavyAtoms(),atoms_with_H=mol.GetNumAtoms(),core_index_0based_to_label=CORE_LABELS,atom_order=[dict(index_0based=a.GetIdx(),element=a.GetSymbol(),neighbors=[x.GetIdx() for x in a.GetNeighbors()]) for a in mol.GetAtoms()],formal_charge=0,neutral_multiplicity=1,lowest_mmff_conformer=best,radii_A=radii,probe_A=1.4)
    (OUT/'target_identity.json').write_text(json.dumps(identity,indent=2),encoding='utf-8')
    arom=sites[sites.aromatic];winners=arom.loc[arom.groupby('conformer_id').heuristic_score.idxmax()].atom_index.value_counts()
    grouped=arom.groupby(['atom_index','site_label']).agg(score_mean=('heuristic_score','mean'),score_sd=('heuristic_score','std'),sasa_mean=('h_sasa_mean_A2','mean'),sasa_sd=('h_sasa_mean_A2','std')).reset_index()
    grouped['top_count']=grouped.atom_index.map(winners).fillna(0).astype(int);grouped.to_csv(OUT/'task_b_site_summary.csv',index=False)
    sensitivity=[]
    for probe in [1.2,1.4,1.8]:
        for algorithm in [rdFreeSASA.LeeRichards,rdFreeSASA.ShrakeRupley]:
            opts=rdFreeSASA.SASAOpts();opts.probeRadius=probe;opts.algorithm=algorithm
            total=rdFreeSASA.CalcSASA(mol,radii,confIdx=best,opts=opts)
            score=[]
            for a in mol.GetAtoms():
                hs=[n for n in a.GetNeighbors() if n.GetAtomicNum()==1]
                if a.GetAtomicNum()==6 and a.GetIsAromatic() and hs:score.append((-.6*float(a.GetProp('_GasteigerCharge'))+.4*sum(h.GetDoubleProp('SASA') for h in hs)/len(hs)/15,a.GetIdx()))
            v,idx=max(score);sensitivity.append(dict(probe_A=probe,algorithm=str(algorithm),total_sasa_A2=total,top_atom_index=idx,top_score=v))
    pd.DataFrame(sensitivity).to_csv(OUT/'task_b_sasa_sensitivity.csv',index=False)
    depiction=Chem.RemoveHs(mol);AllChem.Compute2DCoords(depiction)
    for a in depiction.GetAtoms():a.SetProp('atomNote',f'{a.GetIdx()}:'+CORE_LABELS.get(a.GetIdx(),''))
    drawer=Draw.rdMolDraw2D.MolDraw2DSVG(1100,600);drawer.drawOptions().addAtomIndices=False;drawer.DrawMolecule(depiction,highlightAtoms=[2,7,8,15]);drawer.FinishDrawing()
    (ROOT/'reports/figures/target_atom_map.svg').write_text(drawer.GetDrawingText(),encoding='utf-8')
    png=Draw.MolToImage(depiction,size=(1100,600),highlightAtoms=[2,7,8,15]);png.save(ROOT/'reports/figures/target_atom_map.png')
    return dict(evidence='real molecular force-field geometry and geometric SASA; heuristic ranking unvalidated',requested=32,embedded=len(ids),converged=len(good),formula=identity['formula'],MMFF_energy_range=good.mmff_energy_kcal_mol.agg(['min','max']).to_dict(),SASA_range_A2=good.total_sasa_A2.agg(['min','max']).to_dict(),classified_radii_unique=sorted(set(rawradii)),replacement_radii_by_element={a.GetSymbol():radii[a.GetIdx()] for a in mol.GetAtoms()},core_CH_sites=[label for idx,label in CORE_LABELS.items() if mol.GetAtomWithIdx(idx).GetAtomicNum()==6 and any(n.GetAtomicNum()==1 for n in mol.GetAtomWithIdx(idx).GetNeighbors())],C2_H_count=sum(n.GetAtomicNum()==1 for n in mol.GetAtomWithIdx(7).GetNeighbors()),C5_H_count=sum(n.GetAtomicNum()==1 for n in mol.GetAtomWithIdx(2).GetNeighbors()),aromatic_CH_count=arom.atom_index.nunique(),all_C_H_carbons=sites.atom_index.nunique(),top_counts={str(int(k)):int(v) for k,v in winners.items()},probe_algorithm_sensitivity=sensitivity,limitations=['Gasteiger charges are graph-based neutral partial charges, not radical-cation spin or frontier orbitals.','Conformer counts include repeated wells; rank frequencies are not Boltzmann populations or product probabilities.','Explicit-H SASA partitions are radius/probe dependent; 1.4 A is a convention, not a calibrated MeCN solvent model.'])

def audit_sac():
    data=pd.read_csv(ROOT/'results/repaired/task_c_all_sites.csv');cols=['d_electrons','metal_EN','coord_N','coord_O','d_band_center_proxy_eV']
    rows=[];predictions=[]
    for group in ['metal','coordination_type']:
        for left in data[group].unique():
            train=data[group]!=left;test=~train
            for name,model in [('extra_trees',ExtraTreesRegressor(n_estimators=100,random_state=42)),('mean',DummyRegressor()),('ridge',make_pipeline(StandardScaler(),Ridge(alpha=1.)))]:
                model.fit(data.loc[train,cols],data.loc[train,'activation_barrier_kcal_mol']);p=model.predict(data.loc[test,cols]);y=data.loc[test,'activation_barrier_kcal_mol'].to_numpy()
                rows.append(dict(split=group,held_out=left,model=name,n_train=int(train.sum()),n_test=int(test.sum()),test_mae=mean_absolute_error(y,p),test_rmse=np.sqrt(mean_squared_error(y,p))))
                for idx,pred in zip(data.index[test],p):predictions.append(dict(split=group,held_out=left,model=name,row_id=int(idx),target=float(data.loc[idx,'activation_barrier_kcal_mol']),prediction=float(pred)))
    metrics=pd.DataFrame(rows);metrics.to_csv(OUT/'task_c_group_holdouts.csv',index=False);pd.DataFrame(predictions).to_csv(OUT/'task_c_predictions.csv',index=False)
    source_model=ExtraTreesRegressor(n_estimators=100,random_state=42).fit(data[cols],data.activation_barrier_kcal_mol)
    all_rank=data.sort_values('activation_barrier_kcal_mol');all_rank.to_csv(OUT/'task_c_ranking.csv',index=False)
    stats=metrics.groupby(['split','model'])[['test_mae','test_rmse']].mean().reset_index()
    stats.to_csv(OUT/'task_c_summary.csv',index=False)
    return dict(evidence='synthetic scores; no experimentally synthesized dataset supplied',rows=len(data),in_sample_r2=source_model.score(data[cols],data.activation_barrier_kcal_mol),group_holdout_macro_means=stats.to_dict('records'),rng_scope='Task C has no local seed; sequential master execution inherits NumPy state from Task A. Standalone C is not independently seeded.',source_targets='metal-specific base_barrier + 4.5*(proxy+1.85)^2 + Gaussian noise SD=0.35; two-decimal rounding')

def flow_solution(q_uL,eta=.45,j0=.05,diffusion=1.2e-9):
    F=96485.33;R=8.314;T=298.15;C0=50.;width=.01;height=.0005;length=.1
    q=q_uL*1e-9/60;u=q/(width*height);area_density=1/height;gamma=6*u/height
    km=.67*diffusion*(gamma/(diffusion*length))**(1/3)
    ea=np.exp(.5*F*eta/(R*T));ec=np.exp(-.5*F*eta/(R*T));a=F*km;b=j0/C0*ea;r=j0*ec
    ceq=r/b;keff=km*b/(a+b);decay=keff*area_density/u
    z=np.linspace(0,length,201);cb=ceq+(C0-ceq)*np.exp(-decay*z);cs=(a*cb+r)/(a+b)
    flux=km*(cb-cs);j=F*flux
    out=float(cb[-1]);consumed=(C0-out)*q;volume_L=width*height*length*1000
    ode=solve_ivp(lambda z,y:-area_density/u*km*(y-(a*y+r)/(a+b)),(0,length),[C0],method='DOP853',rtol=1e-11,atol=1e-12)
    integral=np.trapezoid(j,z)*width
    row=dict(q_uL_min=q_uL,eta_V=eta,km_m_s=km,ks_m_s=b/F,keff_m_s=keff,mass_transfer_resistance_fraction=b/(a+b),Da=decay*length,Ceq_mol_m3=ceq,Cout_mol_m3=out,conversion=(C0-out)/C0,residence_s=length/u,STY_mmol_L_h=consumed*1000*3600/volume_L,consumption_mol_s=consumed,current_from_molar_balance_A=F*consumed,current_integrated_A=integral,current_relative_error=abs(integral-F*consumed)/max(abs(F*consumed),1e-30),analytic_ode_error=float(abs(ode.y[0,-1]-out)),delta_effective_over_height=(diffusion/km)/height)
    profile=pd.DataFrame(dict(z_m=z,Cb_mol_m3=cb,Cs_mol_m3=cs,j_A_m2=j,flux_mol_m2_s=flux))
    return row,profile
def audit_flow():
    rows=[]
    for q in [100,300,600,1200]:
        row,profile=flow_solution(q);rows.append(row);profile.to_csv(OUT/f'task_d_profile_{q}.csv',index=False)
    pd.DataFrame(rows).to_csv(OUT/'task_d_analytic_audit.csv',index=False)
    sweep=[]
    for q in [100,300,600,1200]:
        for eta in np.linspace(.05,.60,12):sweep.append(flow_solution(q,float(eta))[0])
    pd.DataFrame(sweep).to_csv(OUT/'task_d_eta_sweep.csv',index=False)
    sensitivity=[]
    for label,j0,diff in [('base',.05,1.2e-9),('half_j0',.025,1.2e-9),('double_j0',.1,1.2e-9),('half_D',.05,.6e-9),('double_D',.05,2.4e-9)]:
        for q in [100,300,600,1200]:sensitivity.append(dict(scenario=label,**flow_solution(q,j0=j0,diffusion=diff)[0]))
    pd.DataFrame(sensitivity).to_csv(OUT/'task_d_parameter_sensitivity.csv',index=False)
    return dict(evidence='reduced transport/kinetic ODE with prescribed uncalibrated parameters, not a resolved PDE',rows=rows,max_analytic_ode_error=max(r['analytic_ode_error'] for r in rows),max_current_integration_relative_error=max(r['current_relative_error'] for r in rows),volume_mL=.5,electrode_area_cm2=10,assumed_electrons=1,flow_sweep_rows=len(sweep),interpretation='km increases with flow; relative mass-transfer resistance decreases while residence time falls. STY is reactant disappearance under assumed unit selectivity, not measured isolated product yield.',limitations=['The code uses inlet concentration in the BV normalization, despite a bulk-concentration docstring.','The constant reverse term assumes fixed product activity; no product transport or counterelectrode balance is solved.','The implemented shear-rate correlation uses 0.67; the stated 1.85 Sherwood formula is not the executed expression.','Effective diffusion thickness is comparable to channel height; Leveque applicability needs validation.','No potential-field PDE, ohmic drop, fouling, Faradaic efficiency, calibrated kinetics, or cell-energy calculation.'])

def target_xtb():
    sys.path.insert(0,str(ROOT.parent/'scripts'));from molecular_descriptors import run_xtb,HARTREE_EV
    import shutil
    executable=os.environ.get('XTB_EXE') or shutil.which('xtb')
    if not executable:return {'status':'unavailable','reason':'xTB not on PATH'}
    folder=ROOT/'results/target_xtb';neutral=folder/'neutral';neutral.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(OUT/'target_lowest_mmff.xyz',neutral/'input.xyz')
    states={};states['neutral']=run_xtb(executable,neutral,['input.xyz','--gfn','2','--alpb','acetonitrile','--chrg','0','--uhf','0','--opt','tight','--json'],180)
    if not(states['neutral']['completed'] and states['neutral'].get('optimization_converged')):return {'status':'failed','states':states}
    for name,chrg in [('cation',1),('anion',-1)]:
        sub=folder/name;sub.mkdir(exist_ok=True);shutil.copyfile(neutral/'xtbopt.xyz',sub/'input.xyz')
        states[name]=run_xtb(executable,sub,['input.xyz','--gfn','2','--alpb','acetonitrile','--chrg',str(chrg),'--uhf','1','--scc','--json'],180)
    if not all(r['completed'] for r in states.values()):return {'status':'failed','states':states}
    mol=make_molecule();charges={n:np.loadtxt(folder/n/'charges') for n in states}
    atoms=[]
    for a in mol.GetAtoms():
        idx=a.GetIdx();atoms.append(dict(atom_index=idx,element=a.GetSymbol(),site_label=CORE_LABELS.get(idx,''),q0=charges['neutral'][idx],q_plus=charges['cation'][idx],q_minus=charges['anion'][idx],removal_response=charges['cation'][idx]-charges['neutral'][idx],addition_response=charges['neutral'][idx]-charges['anion'][idx]))
    pd.DataFrame(atoms).to_csv(folder/'charge_responses.csv',index=False)
    return dict(status='completed',evidence='executed semiempirical molecular calculation, not DFT or measured Eox',jobs=3,method='GFN2-xTB/ALPB(acetonitrile), neutral optimized, +/-1 doublets at fixed neutral nuclei',energies_Eh={n:r['energy_Eh'] for n,r in states.items()},removal_difference_eV=(states['cation']['energy_Eh']-states['neutral']['energy_Eh'])*HARTREE_EV,addition_difference_eV=(states['neutral']['energy_Eh']-states['anion']['energy_Eh'])*HARTREE_EV,charge_sums={n:float(v.sum()) for n,v in charges.items()},C3_removal_response=float(atoms[8]['removal_response']),limitations='No frequency, thermal, reference-electrode calibration or competing transition states; separate equilibrium ALPB response for each charge state.')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--seeds',type=int,default=30);args=parser.parse_args();start=time.perf_counter()
    s={'metadata':{'rdkit':rdBase.rdkitVersion,'threads':1,'seeds':args.seeds}}
    for key,fn in [('A',lambda:audit_bo(args.seeds)),('B',molecular_ensemble),('C',audit_sac),('D',audit_flow),('target_xTB',target_xtb)]:
        print('Starting audit '+key,flush=True);s[key]=fn();(OUT/'audit_summary.json').write_text(json.dumps(s,indent=2),encoding='utf-8')
    s['metadata']['elapsed_seconds']=time.perf_counter()-start
    (OUT/'audit_summary.json').write_text(json.dumps(s,indent=2),encoding='utf-8');print('Audit finished.',flush=True)
if __name__=='__main__':main()
