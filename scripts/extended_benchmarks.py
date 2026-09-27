"""Synthetic robustness experiments. Results are NOT experimental chemistry.

python scripts/extended_benchmarks.py --output results/extended
All random streams, per-run measurements, and assumptions are saved.
"""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ.setdefault(key,'1')
import argparse, json, time, warnings
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import scipy, sklearn
from scipy.integrate import solve_ivp
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, ConstantKernel
from sklearn.ensemble import GradientBoostingRegressor, RandomForestClassifier
from sklearn.linear_model import LinearRegression
from sklearn.dummy import DummyRegressor, DummyClassifier
from sklearn.metrics import mean_squared_error,mean_absolute_error,r2_score,accuracy_score,balanced_accuracy_score,precision_score,recall_score,average_precision_score,brier_score_loss,log_loss
from sklearn.exceptions import ConvergenceWarning

BOUNDS=np.array([[5.,30.],[.05,.30],[20.,60.]])
F=96485.33212
def surface(x):
    x=np.asarray(x)
    return 92-.25*(x[...,0]-16.5)**2-400*(x[...,1]-.15)**2-.08*(x[...,2]-37.5)**2

def summary(values):
    v=np.asarray(values,dtype=float)
    return dict(n=len(v),mean=float(v.mean()),sd=float(v.std(ddof=1)),median=float(np.median(v)),q025=float(np.quantile(v,.025)),q975=float(np.quantile(v,.975)))

def paired_bootstrap(values,seed=999):
    v=np.asarray(values);rng=np.random.default_rng(seed)
    means=v[rng.integers(0,len(v),(10000,len(v)))].mean(axis=1)
    return dict(mean_difference=float(v.mean()),bootstrap_95pct_ci=np.quantile(means,[.025,.975]).tolist(),n_pairs=len(v),bootstrap_replicates=10000)

def run_bo(out,n_seeds):
    rows=[];trajectories=[];observations=[];warning_count=0
    for seed in range(n_seeds):
        # Shared initial design, future candidate pools, and stepwise noise.
        design=np.random.default_rng(10000+seed)
        initial=design.uniform(BOUNDS[:,0],BOUNDS[:,1],(8,3))
        pools=design.uniform(BOUNDS[:,0],BOUNDS[:,1],(12,500,3))
        noise=np.random.default_rng(20000+seed).normal(0,1.5,20)
        for method in ['random','original_gp','scaled_gp']:
            x=initial.copy();y=np.clip(surface(x)+noise[:8],0,99)
            t0=time.perf_counter()
            for step in range(12):
                candidates=pools[step]
                if method=='random':idx=0
                else:
                    if method=='original_gp':
                        gp=GaussianProcessRegressor(kernel=Matern(nu=2.5),alpha=.01,n_restarts_optimizer=5,random_state=30000+seed*12+step)
                        xx,cc,yy=x,candidates,y
                    else:
                        # Fixed observed training mean/scale, no test leakage.
                        scale=max(float(y.std()),1.)
                        xx=(x-BOUNDS[:,0])/(BOUNDS[:,1]-BOUNDS[:,0]);cc=(candidates-BOUNDS[:,0])/(BOUNDS[:,1]-BOUNDS[:,0]);yy=(y-y.mean())/scale
                        kernel=ConstantKernel(1.,(.01,100))*Matern(length_scale=[.3,.3,.3],length_scale_bounds=(.03,10),nu=2.5)
                        gp=GaussianProcessRegressor(kernel=kernel,alpha=2.25/scale**2,n_restarts_optimizer=1,random_state=30000+seed*12+step)
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter('always',ConvergenceWarning)
                        gp.fit(xx,yy)
                    warning_count+=sum(isinstance(w.message,ConvergenceWarning) for w in caught)
                    mu,sigma=gp.predict(cc,return_std=True)
                    idx=int(np.argmax(mu+1.96*sigma))
                chosen=candidates[idx]
                x=np.vstack([x,chosen]);y=np.append(y,np.clip(surface(chosen)+noise[8+step],0,99))
                latent=surface(x)
                trajectories.append(dict(seed=seed,method=method,evaluations=len(x),oracle_best_queried=float(latent.max()),observed_best=float(y.max()),oracle_simple_regret=float(92-latent.max())))
            latent=surface(x);selected=int(np.argmax(y))
            rows.append(dict(seed=seed,method=method,oracle_simple_regret=float(92-latent.max()),recommendation_regret=float(92-latent[selected]),selected_noise_bias=float(y[selected]-latent[selected]),observed_best=float(y.max()),elapsed_s=time.perf_counter()-t0))
            for i,(point,obs,true) in enumerate(zip(x,y,latent)):
                observations.append(dict(seed=seed,method=method,evaluation=i+1,j=float(point[0]),concentration=float(point[1]),epsilon=float(point[2]),observed=float(obs),latent=float(true)))
        if (seed+1)%5==0:print(f'BO seeds completed: {seed+1}/{n_seeds}',flush=True)
    df=pd.DataFrame(rows);df.to_csv(out/'bo_runs.csv',index=False)
    pd.DataFrame(trajectories).to_csv(out/'bo_trajectories.csv',index=False)
    pd.DataFrame(observations).to_csv(out/'bo_observations.csv',index=False)
    metrics={m:{k:summary(g[k]) for k in ['oracle_simple_regret','recommendation_regret','selected_noise_bias']} for m,g in df.groupby('method')}
    pivot=df.pivot(index='seed',columns='method',values='recommendation_regret')
    pairs={f'{m}_minus_random':paired_bootstrap(pivot[m]-pivot['random']) for m in ['original_gp','scaled_gp']}
    pairs['scaled_gp_minus_original_gp']=paired_bootstrap(pivot['scaled_gp']-pivot['original_gp'])
    return dict(evidence='synthetic optimizer comparison',seeds=n_seeds,evaluations_per_campaign=20,methods=metrics,paired_recommendation_regret=pairs,kernel_boundary_warning_count=warning_count,scope='Original GP hyperparameters reproduced, but this is a newly seeded paired design, not a replay of the original RNG stream. Bootstrap intervals concern synthetic campaign means, not wet-lab benefit.')

def redox_data(n,rng):
    x=rng.normal(size=(n,4));latent=-.85*x[:,0]+.3*x[:,2]+1.65
    return x,latent+rng.normal(0,.08,n)

def run_regression(out,n_seeds):
    records=[]
    xtest,ytest=redox_data(5000,np.random.default_rng(555))
    pd.DataFrame(np.column_stack([xtest,ytest]),columns=['x0','x1','x2','x3','synthetic_target']).to_csv(out/'regression_test.csv',index=False)
    for seed in range(n_seeds):
        for n in [20,60,120,240]:
            x,y=redox_data(n,np.random.default_rng(40000+seed*1000+n))
            for name,model in [('mean',DummyRegressor()),('linear',LinearRegression()),('gradient_boosting',GradientBoostingRegressor(n_estimators=50,random_state=42))]:
                model.fit(x,y);pred=model.predict(xtest)
                records.append(dict(seed=seed,n_train=n,model=name,train_r2=float(model.score(x,y)),test_r2=r2_score(ytest,pred),test_rmse=float(np.sqrt(mean_squared_error(ytest,pred))),test_mae=mean_absolute_error(ytest,pred)))
    df=pd.DataFrame(records);df.to_csv(out/'regression_runs.csv',index=False)
    # OOD test uses the synthetic target formula only, not voltage semantics.
    x,y=redox_data(60,np.random.default_rng(40404))
    oo=np.array([[-5.8,1.42,.88,22.4],[-5.8,.95,.54,18.1],[-5.8,.32,.12,35.8]])
    expected=-.85*oo[:,0]+.3*oo[:,2]+1.65
    ood={}
    for name,model in [('linear',LinearRegression()),('gradient_boosting',GradientBoostingRegressor(n_estimators=50,random_state=42))]:
        model.fit(x,y);pred=model.predict(oo)
        ood[name]=dict(predictions=pred.tolist(),absolute_errors_against_synthetic_formula=np.abs(pred-expected).tolist())
    return dict(evidence='synthetic regression validation; targets have no calibrated voltage units',independent_test_rows=5000,training_sizes=[20,60,120,240],training_seeds=n_seeds,at_n60={m:{k:summary(g[k]) for k in ['train_r2','test_r2','test_rmse','test_mae']} for m,g in df[df.n_train==60].groupby('model')},ood=ood,ood_reference=expected.tolist(),scope='All methods share one independent test set; variation across seeds is training-sample sensitivity, not independent-test replication.')

def run_sac(out):
    metals=['Cu','Co','Ni','Pd'];coords=['N4','N3C1','N2O2']
    labels=[m+'-'+c for m in metals for c in coords]
    means=np.array([18.5-.8*({'Cu':9,'Co':7,'Ni':8,'Pd':8}[m]-8)**2+.3*({'N4':12.1,'N3C1':10.8,'N2O2':11.5}[c]-11)**2 for m in metals for c in coords])
    records=[]
    for sd in [.1,.4,.8]:
        draws=means+np.random.default_rng(8800+int(sd*10)).normal(0,sd,(100000,12))
        ranks=np.argsort(np.argsort(draws,axis=1),axis=1)+1
        for i,label in enumerate(labels):
            p=float(np.mean(ranks[:,i]==1))
            records.append(dict(site=label,noise_sd=sd,noiseless_mean=float(means[i]),probability_rank1=p,mean_rank=float(ranks[:,i].mean()),monte_carlo_se=float(np.sqrt(p*(1-p)/len(draws)))))
    df=pd.DataFrame(records);df.to_csv(out/'sac_rank_sensitivity.csv',index=False)
    return dict(evidence='synthetic generator resampling, not a posterior about real catalysts',draws_per_noise_level=100000,noise_sds=[.1,.4,.8],at_original_sd=df[df.noise_sd==.4].to_dict('records'),original_winner='Co-N4',original_winner_probability=float(df[(df.noise_sd==.4)&(df.site=='Co-N4')].probability_rank1.iloc[0]),pair_CoN4_beats_CuN4_analytic=.5,scope='Probabilities use unconditioned generator means and independent Gaussian noise; they are not confidence in the experimentally best material.')

def classify_label(x):return ((x[:,0]<1.2)&(x[:,1]<1.8)&(x[:,2]>.8)).astype(int)
def class_metrics(y,p):
    pred=(p>=.5).astype(int)
    return dict(accuracy=accuracy_score(y,pred),balanced_accuracy=balanced_accuracy_score(y,pred),precision=precision_score(y,pred,zero_division=0),recall=recall_score(y,pred,zero_division=0),average_precision=average_precision_score(y,p),brier=brier_score_loss(y,p),log_loss=log_loss(y,p,labels=[0,1]))

def run_classification(out,n_seeds):
    rng=np.random.default_rng(777);test=rng.uniform(.1,2.5,(5000,4));ytest=classify_label(test)
    # Separate hard cases: points close to any decision boundary.
    distances=np.min(np.abs(test[:,:3]-[1.2,1.8,.8]),axis=1);boundary=distances<.08
    records=[];calibration=[]
    for seed in range(n_seeds):
        x=np.random.default_rng(70000+seed).uniform(.1,2.5,(80,4));y=classify_label(x)
        rf=RandomForestClassifier(n_estimators=30,random_state=42);rf.fit(x,y)
        p=rf.predict_proba(test)[:,1]
        for subset,mask in [('all',np.ones(len(test),dtype=bool)),('near_boundary',boundary)]:
            records.append(dict(seed=seed,subset=subset,model='random_forest',train_accuracy=float(rf.score(x,y)),**class_metrics(ytest[mask],p[mask])))
        baseline=np.full(len(test),y.mean())
        records.append(dict(seed=seed,subset='all',model='prevalence',train_accuracy=float(np.mean(y==int(y.mean()>=.5))),**class_metrics(ytest,baseline)))
        if seed==0:
            # Assign each row once; floating arange edges can omit p=0.3, 0.6, ...
            bins=np.minimum((p*10).astype(int),9)
            for bin_id in range(10):
                low=bin_id/10
                mask=bins==bin_id
                if mask.any():calibration.append(dict(lower=float(low),n=int(mask.sum()),mean_probability=float(p[mask].mean()),observed_frequency=float(ytest[mask].mean())))
            pd.DataFrame(dict(y=ytest,probability=p,near_boundary=boundary)).to_csv(out/'classification_test_predictions_seed0.csv',index=False)
    df=pd.DataFrame(records);df.to_csv(out/'classification_runs.csv',index=False)
    pd.DataFrame(calibration).to_csv(out/'classification_reliability_seed0.csv',index=False)
    return dict(evidence='held-out evaluation of a synthetic threshold rule',training_rows=80,test_rows=5000,test_positive_count=int(ytest.sum()),boundary_test_rows=int(boundary.sum()),training_seeds=n_seeds,summary={f'{m}/{s}':{k:summary(g[k]) for k in ['train_accuracy','accuracy','balanced_accuracy','precision','recall','average_precision','brier','log_loss']} for (m,s),g in df.groupby(['model','subset'])},hardcoded_0892_used=False,scope='No actual reaction labels or calibrated chemical success probabilities are available. Reliability bins are descriptive for seed 0, not a post-hoc calibration model.')

def sequential_solution(t,k1,k2):
    a=np.exp(-k1*t)
    if np.isclose(k1,k2,rtol=1e-10,atol=0):p=k1*t*a
    else:p=k1/(k2-k1)*(np.exp(-k1*t)-np.exp(-k2*t))
    return np.array([a,p,1-a-p])

def flow_metrics(q,k1=.015,k2=.002,volume_ml=1.,c0=.1,mw=200.,current=.2):
    tau=60*volume_ml/q;a,p,d=sequential_solution(tau,k1,k2);x=1-a
    feed_mol_s=q/1000/60*c0
    # Both A->P and P->D assumed to require two electrons per molecule.
    i_required=2*F*feed_mol_s*(x+d)
    rate=q/1000*60*c0*p*mw
    return dict(q_ml_min=float(q),residence_s=float(tau),conversion=float(x),product_yield=float(p),degraded_fraction=float(d),selectivity=float(p/x),product_g_h=float(rate),sty_g_L_h=float(rate/(volume_ml/1000)),minimum_current_A=float(i_required),available_current_A=current,charge_feasible=bool(i_required<=current),candidate_FE_at_available_current=float(2*F*feed_mol_s*p/current),energy_kWh_kg_at_assumed_3V=float(3*current/(rate/1000)/1000),mass_balance_residual=float(a+p+d-1))

def run_flow(out):
    grid=[flow_metrics(q) for q in np.linspace(.1,2,191)]
    df=pd.DataFrame(grid);df.to_csv(out/'flow_grid.csv',index=False)
    comparisons=[]
    for q in [.2,.5,1.,2.]:
        tau=60/q
        sol=solve_ivp(lambda t,y:[-.015*y[0],.015*y[0]-.002*y[1],.002*y[1]],(0,tau),[1.,0.,0.],rtol=1e-10,atol=1e-12)
        analytic=sequential_solution(tau,.015,.002)
        comparisons.append(dict(**flow_metrics(q),max_species_analytic_error=float(np.max(np.abs(sol.y[:,-1]-analytic)))))
    pd.DataFrame(comparisons).to_csv(out/'flow_selected.csv',index=False)
    rng=np.random.default_rng(5555);k1s=.015*np.exp(rng.normal(0,.2,10000));k2s=.002*np.exp(rng.normal(0,.2,10000))
    uncertain=[]
    for q in [.2,.5,1.,2.]:
        vals=[flow_metrics(q,k1,k2) for k1,k2 in zip(k1s,k2s)]
        uncertain.append(dict(q_ml_min=q,sty_q025=float(np.quantile([v['sty_g_L_h'] for v in vals],.025)),sty_q975=float(np.quantile([v['sty_g_L_h'] for v in vals],.975)),conversion_q025=float(np.quantile([v['conversion'] for v in vals],.025)),conversion_q975=float(np.quantile([v['conversion'] for v in vals],.975)),charge_feasible_fraction=float(np.mean([v['charge_feasible'] for v in vals]))))
    pd.DataFrame(uncertain).to_csv(out/'flow_parameter_sensitivity.csv',index=False)
    feasible=df[df.charge_feasible&(df.conversion>=.8)&(df.selectivity>=.8)]
    best=None if feasible.empty else feasible.loc[feasible.sty_g_L_h.idxmax()].to_dict()
    return dict(evidence='dimensionally consistent uncalibrated A->P->D design scenario',assumptions=dict(volume_ml=1.,feed_M=.1,product_MW_g_mol=200.,k1_s_inverse=.015,k2_s_inverse=.002,current_A=.2,cell_voltage_V=3.,electrons_per_each_step=2,unit_stoichiometry=True,product_initial_concentration=0.,plug_flow=True),selected=comparisons,best_grid_point_under_constraints=best,constraints=dict(min_conversion=.8,min_selectivity=.8,current_required_at_most_A=.2),parameter_sensitivity=uncertain,uncertainty_assumptions='10,000 independent lognormal k1 and k2 draws, log-SD 0.2; these are scenario intervals, not experimental confidence intervals.',scope='Charge-infeasible rows are counterfactual kinetic outputs and must not be reported as realizable production. No Butler-Volmer, measured transport, calibration, process energy balance, or industrial optimum is established.')

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=Path('results/extended'));p.add_argument('--seeds',type=int,default=30);args=p.parse_args()
    out=args.output;out.mkdir(parents=True,exist_ok=True);t0=time.perf_counter()
    result={'metadata':dict(started_utc=datetime.now(timezone.utc).isoformat(),numpy=np.__version__,scipy=scipy.__version__,sklearn=sklearn.__version__,n_seeds=args.seeds,threads=1)}
    for name,fn in [('bo',lambda:run_bo(out,args.seeds)),('regression',lambda:run_regression(out,args.seeds)),('sac',lambda:run_sac(out)),('classification',lambda:run_classification(out,args.seeds)),('flow',lambda:run_flow(out))]:
        print('Starting '+name,flush=True);result[name]=fn();(out/'summary.json').write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
    result['metadata'].update(elapsed_seconds=time.perf_counter()-t0,completed_utc=datetime.now(timezone.utc).isoformat())
    (out/'summary.json').write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
    print('All extended synthetic calculations completed.',flush=True)

if __name__=='__main__':main()
