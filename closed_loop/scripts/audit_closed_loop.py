"""Reproducible arithmetic and mock-feedback audit, with explicit evidence roles.

No wet experiment, DFT calculation, or hardware control is performed.
"""
from pathlib import Path
import argparse,ast,csv,hashlib,importlib.metadata,json,math,sys,warnings
import numpy as np
from scipy.special import ndtr
from scipy.spatial import Delaunay
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel,Matern,WhiteKernel
from sklearn.preprocessing import StandardScaler
from rdkit import Chem
from rdkit.Chem import Descriptors,rdMolDescriptors

ROOT=Path(__file__).resolve().parents[1]
F=96485.33
SUBSTRATE='COc1ccc2[nH]c(cc2c1)c3ccccc3'
PARTNER='CSc1ccccc1'
FEATURES=['current_density_mA_cm2','electrolyte_concentration_M','temperature_K']

def save_json(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

def save_csv(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)

def source_arrays():
    tree=ast.parse((ROOT/'source/run_closed_loop_platform.py').read_text(encoding='utf-8'))
    f=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='assimilate_lab_feedback_and_update')
    out={}
    for n in f.body:
        if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Name) and n.targets[0].id in ['X_prior','y_prior_pred','y_wet_lab_actual']:
            out[n.targets[0].id]=np.asarray(ast.literal_eval(n.value.args[0]),dtype=float)
    return out['X_prior'],out['y_prior_pred'],out['y_wet_lab_actual']

def stoichiometry(scale_mmol=.2,volume_mL=6,j_mA_cm2=12.5,area_cm2=1.5,charge_factor=2.2,electrolyte_M=.1,meCN_fraction=.8):
    values=[scale_mmol,volume_mL,j_mA_cm2,area_cm2,charge_factor,electrolyte_M]
    if not all(math.isfinite(v) and v>0 for v in values) or not 0<meCN_fraction<1:raise ValueError('Invalid positive finite input or solvent fraction')
    mol=Chem.MolFromSmiles(SUBSTRATE);mw=Descriptors.MolWt(mol)
    I=j_mA_cm2*area_cm2;Q=scale_mmol/1000*F*charge_factor
    return dict(scale_mmol=scale_mmol,volume_mL=volume_mL,substrate_mg=scale_mmol*mw,substrate_M=scale_mmol/volume_mL,substrate_MW=mw,
        electrolyte_M=electrolyte_M,electrolyte_mmol=electrolyte_M*volume_mL,electrolyte_mg=electrolyte_M*volume_mL*387.43,
        MeCN_mL=volume_mL*meCN_fraction,HFIP_mL=volume_mL*(1-meCN_fraction),current_density_mA_cm2=j_mA_cm2,area_cm2=area_cm2,
        current_mA=I,charge_F_per_mol=charge_factor,charge_C=Q,duration_min=Q/(I/1000)/60,
        max_product_mmol_at_2e=Q/(2*F)*1000,full_conversion_FE_ceiling_pct=100*2/charge_factor)

def _evidence_file(root,name,digest):
    p=(root/name).resolve()
    if not name or not p.is_relative_to(root.resolve()) or not p.is_file():raise ValueError('Missing or nonlocal evidence file')
    if len(digest)!=64 or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:raise ValueError('Evidence hash mismatch')

def load_feedback(path,role='mock'):
    """Check role, units, identity and evidence linkage before fitting; not assay certification."""
    with Path(path).open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
    if len(rows)<4:raise ValueError('At least four observations required for this demonstration')
    seen=set();identities=set()
    for r in rows:
        if role not in {'mock','experimental'} or r.get('evidence_role')!=role:raise ValueError('Mixed or incorrect evidence roles')
        if not r.get('run_id') or r['run_id'] in seen:raise ValueError('Missing or duplicate run_id')
        seen.add(r['run_id'])
        for k in FEATURES+['yield_pct','prior_prediction_pct']:
            v=float(r[k])
            if not math.isfinite(v):raise ValueError('Nonfinite numeric field')
            r[k]=v
        if not 0<=r['yield_pct']<=100 or not 0<=r['prior_prediction_pct']<=100:raise ValueError('Yield outside percent scale')
        if any(r[k]<=0 for k in FEATURES) or r['temperature_K']<200:raise ValueError('Invalid condition or non-Kelvin temperature')
        if not r.get('reaction_id'):raise ValueError('Missing reaction identity')
        identity=[r['reaction_id']]
        if role=='experimental':
            for k in ['substrate_smiles','partner_smiles','product_smiles']:
                mol=Chem.MolFromSmiles(r.get(k,''))
                if mol is None or mol.GetNumAtoms()==0:raise ValueError('Missing valid reaction structure')
                identity.append(Chem.MolToSmiles(mol))
            if r.get('assay')!='HPLC' or not r.get('operator') or not r.get('recorded_at'):raise ValueError('Missing assay metadata')
            for kind in ['assay','calibration','electrolysis']:_evidence_file(Path(path).parent,r.get(kind+'_file',''),r.get(kind+'_sha256',''))
        identities.add(tuple(identity))
    if len(identities)!=1:raise ValueError('Multiple reaction identities require a different model')
    return rows

def make_mock_feedback(path):
    X,p,y=source_arrays();rows=[]
    for i,(x,yp,yy) in enumerate(zip(X,p,y)):
        row=dict(run_id=f'mock-{i+1:02}',evidence_role='mock',reaction_id='unspecified_demo',**dict(zip(FEATURES,map(float,x))),yield_pct=float(yy),prior_prediction_pct=float(yp),substrate_smiles='',partner_smiles='',product_smiles='',assay='HPLC_label_only',operator='',recorded_at='')
        for kind in ['assay','calibration','electrolysis']:row[kind+'_file']='';row[kind+'_sha256']=''
        rows.append(row)
    save_csv(path,rows)

def fit_gp(X,y,method='scaled_mle',noise_sd=2.):
    scaler=None
    if method=='original':kernel=Matern(nu=2.5)+WhiteKernel(noise_level=.5);normalize=False;restarts=5;optimizer='fmin_l_bfgs_b'
    else:
        scaler=StandardScaler().fit(X);X=scaler.transform(X);normalize=True
        if method=='scaled_fixed':
            noise_var=(noise_sd/max(float(np.std(y)),1e-12))**2
            kernel=ConstantKernel(1,'fixed')*Matern(1,'fixed',nu=2.5)+WhiteKernel(noise_var,'fixed');restarts=0;optimizer=None
        else:
            kernel=ConstantKernel(1,(.01,100))*Matern(1,(.05,20),nu=2.5)+WhiteKernel(.02,(1e-5,1));restarts=2;optimizer='fmin_l_bfgs_b'
    with warnings.catch_warnings(record=True) as messages:
        warnings.simplefilter('always')
        model=GaussianProcessRegressor(kernel=kernel,normalize_y=normalize,n_restarts_optimizer=restarts,random_state=42,optimizer=optimizer).fit(X,y)
    return model,scaler,[str(m.message) for m in messages]

def predictions(model,scaler,X):
    xx=X if scaler is None else scaler.transform(X)
    mean,obs_std=model.predict(xx,return_std=True)
    yscale=np.asarray(model._y_train_std).item()
    noise_var=float(model.kernel_.k2.noise_level)*yscale**2
    latent_std=np.sqrt(np.maximum(obs_std**2-noise_var,0))
    return mean,latent_std,obs_std

def expected_improvement(mu,sigma,best,xi=0.):
    mu,sigma=np.broadcast_arrays(np.asarray(mu,float),np.asarray(sigma,float))
    if np.any(sigma<0):raise ValueError('Negative uncertainty')
    improvement=mu-best-xi;positive=sigma>1e-14;safe=np.where(positive,sigma,1);z=improvement/safe
    ei=improvement*ndtr(z)+safe*np.exp(-z*z/2)/np.sqrt(2*np.pi)
    return np.where(positive,np.maximum(ei,0),np.maximum(improvement,0))

def candidate_pool():return np.array([[j,c,t] for j in np.linspace(8,22,15) for c in [.08,.10,.12] for t in [298.15,303.15]])

def unseen_mask(pool,X):return ~np.any(np.all(np.isclose(pool[:,None,:],X[None,:,:],atol=1e-10,rtol=0),axis=2),axis=1)

def active_learning(path,out,role='mock'):
    rows=load_feedback(path,role);X=np.array([[r[k] for k in FEATURES] for r in rows]);y=np.array([r['yield_pct'] for r in rows]);prior=np.array([r['prior_prediction_pct'] for r in rows])
    folds=[];metrics={};fits={};pool=candidate_pool();valid=unseen_mask(pool,X);candidate_rows=[]
    for method in ['original','scaled_mle','scaled_fixed','mean_baseline']:
        errors=[]
        for i in range(len(y)):
            keep=np.arange(len(y))!=i
            if method=='mean_baseline':pred=float(np.mean(y[keep]));sigma=None;w=[]
            else:
                gp,sc,w=fit_gp(X[keep],y[keep],method);preds,_,obs=predictions(gp,sc,X[i:i+1]);pred=float(preds[0]);sigma=float(obs[0])
            errors.append(pred-y[i]);folds.append(dict(method=method,held_out_run=rows[i]['run_id'],observed_yield_pct=float(y[i]),predicted_yield_pct=pred,error_pp=float(pred-y[i]),observation_sigma_pp=sigma,warnings=' | '.join(w)))
        metrics[method]=dict(MAE_pp=float(np.mean(np.abs(errors))),RMSE_pp=float(np.sqrt(np.mean(np.square(errors)))))
        if method=='mean_baseline':continue
        gp,sc,w=fit_gp(X,y,method);mu,latent,obs=predictions(gp,sc,pool);ei=expected_improvement(mu,latent,float(y.max()));ucb=mu+2.5*latent
        fits[method]=dict(kernel=str(gp.kernel_),log_marginal_likelihood=float(gp.log_marginal_likelihood_value_),warnings=w,noise_variance_pp2=float(gp.kernel_.k2.noise_level)*np.asarray(gp._y_train_std).item()**2)
        for i,x in enumerate(pool):candidate_rows.append(dict(method=method,candidate=i,**dict(zip(FEATURES,map(float,x))),temperature_C=float(x[2]-273.15),unseen=bool(valid[i]),mu_pct=float(mu[i]),latent_sigma_pp=float(latent[i]),observation_sigma_pp=float(obs[i]),EI_pp=float(ei[i]),UCB_latent_pct=float(ucb[i]),UCB_observation_pct=float(mu[i]+2.5*obs[i])))
    save_csv(out/'loo_predictions.csv',folds);save_csv(out/'candidate_predictions.csv',candidate_rows)
    selections={}
    for method in ['original','scaled_mle','scaled_fixed']:
        candidates=[r for r in candidate_rows if r['method']==method and r['unseen']]
        selections[method]={policy:sorted(candidates,key=lambda r:(-r[policy],r['candidate']))[:5] for policy in ['EI_pp','UCB_latent_pct','UCB_observation_pct']}
    scenarios=[]
    for noise in [.5,2.,5.]:
        gp,sc,w=fit_gp(X,y,'scaled_fixed',noise);mu,sigma,obs=predictions(gp,sc,pool)
        for xi in [0.,1.,3.]:
            ei=expected_improvement(mu,sigma,float(y.max()),xi);idx=int(np.argmax(np.where(valid,ei,-np.inf)))
            scenarios.append(dict(noise_sd_pp=noise,xi_pp=xi,candidate=idx,**dict(zip(FEATURES,map(float,pool[idx]))),mu_pct=float(mu[idx]),latent_sigma_pp=float(sigma[idx]),EI_pp=float(ei[idx])))
    save_csv(out/'acquisition_sensitivity.csv',scenarios)
    # Linear barycentric hull membership is a domain diagnostic, not chemical support.
    scaled=StandardScaler().fit(X)
    try:in_hull=Delaunay(scaled.transform(X)).find_simplex(scaled.transform(pool))>=0;hull_count=int(in_hull.sum())
    except Exception:hull_count=None
    best=selections['scaled_mle']['EI_pp'][0]
    # Independent Monte Carlo sanity check at the reported EI maximizer.
    samples=np.random.default_rng(4321).normal(best['mu_pct'],best['latent_sigma_pp'],200000)
    gains=np.maximum(samples-y.max(),0)
    result=dict(evidence_role=role,feedback_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),n=len(y),given_prediction_MAE_pp=float(np.mean(np.abs(prior-y))),LOO=metrics,kernels=fits,pool_count=len(pool),unseen_count=int(valid.sum()),pool_inside_training_hull=hull_count,best_observed_label_pct=float(y.max()),selections=selections,noise_scenarios=len(scenarios),EI_monte_carlo=dict(analytic_pp=best['EI_pp'],estimated_pp=float(gains.mean()),standard_error_pp=float(gains.std(ddof=1)/np.sqrt(len(gains))),draws=len(gains)),acquisition='Plug-in analytic one-point EI on latent variance; not noise-integrated EI or joint-batch EI',validation_scope='Role/identity/file-hash linkage and numerical tests; no certification of assay validity or scientific generalization')
    save_json(out/'active_learning.json',result);return result

def arithmetic_audit(out):
    base=stoichiometry();sweep=[]
    for scale in [.1,.2,.5]:
        for volume in [3.,6.,12.]:
            for area in [.75,1.5,3.]:
                for j in [8.,12.5,20.]:sweep.append(stoichiometry(scale,volume,j,area))
    save_csv(out/'stoichiometry_sweep.csv',sweep)
    equivalents=[]
    for name,smiles in [('thioanisole',PARTNER),('thiophenol','Sc1ccccc1')]:
        mol=Chem.MolFromSmiles(smiles);mw=Descriptors.MolWt(mol)
        for eq in [.5,1.,1.5,2.]:equivalents.append(dict(identity=name,smiles=smiles,formula=rdMolDescriptors.CalcMolFormula(mol),MW_g_mol=mw,assumed_equivalents=eq,mmol=.2*eq,mass_mg=.2*eq*mw,status='Arithmetic scenario only; partner and loading are not selected'))
    save_csv(out/'partner_scenarios.csv',equivalents)
    fe=[]
    for yy in [25,50,75,100]:
        fe.append(dict(assumed_yield_pct=yy,assumed_electrons_per_product=2,product_mmol=.2*yy/100,FE_pct=100*2*F*(.2e-3*yy/100)/base['charge_C']))
    save_csv(out/'charge_balance_scenarios.csv',fe)
    # Independent bounded quadrature of the stated constant-current profile.
    duration=base['duration_min']*60;time=np.linspace(0,duration,1001);charge=float(np.trapezoid(np.full_like(time,base['current_mA']/1000),time))
    outdata=dict(base=base,sweep_rows=len(sweep),partner_scenarios=len(equivalents),charge_quadrature_C=charge,charge_quadrature_absolute_error_C=abs(charge-base['charge_C']),ready_for_wet_lab=False,missing=['Defined product and balanced transformation','Selected coupling partner and equivalents','Validated assay and response factor','Quench compatibility and preparation','Electrode geometry/area convention, gap, stirring and temperature','Solubility and isolation method'],findings=['Source coupling_partner_smiles is unused','Thioanisole is not thiophenol','Source solvent volumes do not respond to changed total volume','Source current density is an assumption, not a linked earlier optimum','Source ORCA input uses CPCM, not SMD','Source MAE is not cross-validation; labels are mock','Source acquisition is UCB, not EI'])
    save_json(out/'stoichiometry.json',outdata);return outdata

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'results/audit');p.add_argument('--feedback',type=Path);p.add_argument('--role',choices=['mock','experimental'],default='mock');a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    if a.feedback is None:
        if a.role!='mock':raise ValueError('Experimental mode requires an explicit evidence-linked CSV')
        feedback=ROOT/'data/mock_feedback.csv';make_mock_feedback(feedback)
    else:feedback=a.feedback
    sto=arithmetic_audit(a.output);al=active_learning(feedback,a.output,a.role)
    versions={name:importlib.metadata.version(name) for name in ['numpy','scipy','pandas','scikit-learn','rdkit','matplotlib']}
    summary=dict(date='2026-09-28',evidence='Executed arithmetic/model audit; source feedback is mock; no wet experiment or DFT execution',stoichiometry=sto['base'],GP_LOO=al['LOO'],source_comparison_MAE_pp=al['given_prediction_MAE_pp'],selected_EI_condition=al['selections']['scaled_mle']['EI_pp'][0],packages=versions)
    save_json(a.output/'audit_summary.json',summary);print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
