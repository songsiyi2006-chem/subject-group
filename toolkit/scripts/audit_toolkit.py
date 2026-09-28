"""Analytical bookkeeping, synthetic recovery and figure-claim audit.

No raw instrument data or electronic-structure calculations are represented as measured.
"""
from pathlib import Path
import argparse,ast,csv,hashlib,importlib.metadata,json,math,warnings
import numpy as np
from scipy.stats import pearsonr
ROOT=Path(__file__).resolve().parents[1]

def jwrite(path,data):
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def cwrite(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def source_literals(function):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',SyntaxWarning);tree=ast.parse((ROOT/'source/run_deployment_and_figures.py').read_text(encoding='utf-8'))
    fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==function);out={}
    for n in fn.body:
        if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Name):
            try:out[n.targets[0].id]=ast.literal_eval(n.value)
            except (ValueError,TypeError):pass
    return out

def finite_positive(**values):
    if not all(math.isfinite(v) and v>0 for v in values.values()):raise ValueError('Inputs must be finite and positive')
def hplc(area,slope,intercept,dilution,volume_mL,substrate_mmol,calibration_range_mM=None):
    finite_positive(slope=slope,dilution=dilution,volume=volume_mL,substrate=substrate_mmol)
    if not math.isfinite(area) or not math.isfinite(intercept):raise ValueError('Nonfinite signal')
    measured=(area-intercept)/slope;actual=measured*dilution;amount=actual*volume_mL/1000;yield_pct=100*amount/substrate_mmol
    flags=[]
    if measured<0:flags.append('below_intercept_not_a_negative_product_amount')
    if not 0<=yield_pct<=100:flags.append('outside_stoichiometric_yield_range_for_1_to_1_product')
    if calibration_range_mM is None:flags.append('calibration_range_unknown')
    elif not calibration_range_mM[0]<=measured<=calibration_range_mM[1]:flags.append('outside_calibration_range')
    return dict(measured_mM=measured,reaction_mM=actual,product_mmol=amount,yield_pct=yield_pct,flags=flags,numeric_range_pass=not flags,measurement_validated=False)

def qnmr(product_integral,standard_integral,product_H,standard_H,standard_mass_mg,standard_MW,standard_purity,substrate_mmol,sample_fraction=1.):
    finite_positive(standard_integral=standard_integral,product_H=product_H,standard_H=standard_H,mass=standard_mass_mg,MW=standard_MW,substrate=substrate_mmol)
    if not math.isfinite(product_integral) or product_integral<0:raise ValueError('Product integral must be finite and nonnegative')
    if not 0<standard_purity<=1 or not 0<sample_fraction<=1:raise ValueError('Purity and sampled fraction must be in (0,1]')
    standard_mmol=standard_mass_mg/standard_MW*standard_purity
    amount=(product_integral/product_H)/(standard_integral/standard_H)*standard_mmol/sample_fraction
    yield_pct=100*amount/substrate_mmol
    return dict(standard_mmol=standard_mmol,product_mmol=amount,yield_pct=yield_pct,flags=[] if 0<=yield_pct<=100 else ['outside_stoichiometric_yield_range_for_1_to_1_product'],measurement_validated=False)

def parse_nmr_csv(path,role):
    """Parse assigned, pre-integrated peaks; no FID processing or peak assignment."""
    with Path(path).open(encoding='utf-8-sig',newline='') as f:records=list(csv.DictReader(f))
    out=[];seen=set()
    for r in records:
        if role not in {'source_example','synthetic','experimental'} or r.get('evidence_role')!=role:raise ValueError('Evidence role mismatch')
        if not r.get('sample_id') or r['sample_id'] in seen:raise ValueError('Missing or duplicate sample ID')
        seen.add(r['sample_id'])
        if r.get('mass_unit')!='mg' or r.get('amount_unit')!='mmol':raise ValueError('Explicit mg and mmol units required')
        if role=='experimental':
            if any(not r.get(k) for k in ['product_peak_assignment','product_identity','standard_identity','relaxation_record','analyst','recorded_at']):raise ValueError('Missing quantitative NMR metadata')
            for kind in ['raw','standard_certificate']:
                name=r.get(kind+'_file','');p=(Path(path).parent/name).resolve()
                if not name or not p.is_relative_to(Path(path).parent.resolve()) or not p.is_file():raise ValueError('Missing local NMR evidence')
                if hashlib.sha256(p.read_bytes()).hexdigest()!=r.get(kind+'_sha256'):raise ValueError('NMR evidence hash mismatch')
        values=[float(r[k]) for k in ['product_integral','standard_integral','product_H','standard_H','standard_mass_mg','standard_MW','standard_purity','substrate_mmol','sample_fraction']]
        out.append(dict(sample_id=r['sample_id'],evidence_role=role,**qnmr(*values)))
    if not out:raise ValueError('No NMR records')
    return out

def parse_integrated_csv(path,role):
    """Read pre-integrated numerical data; not a vendor chromatogram or FID parser."""
    with Path(path).open(encoding='utf-8-sig',newline='') as f:records=list(csv.DictReader(f))
    out=[];seen=set()
    for r in records:
        if role not in {'source_example','synthetic','experimental'} or r.get('evidence_role')!=role:raise ValueError('Evidence role mismatch')
        if not r.get('sample_id') or r['sample_id'] in seen:raise ValueError('Missing or duplicate sample ID')
        seen.add(r['sample_id'])
        if r.get('area_mode')!='absolute' or r.get('concentration_unit')!='mM' or r.get('volume_unit')!='mL':raise ValueError('Only absolute area, mM and mL supported')
        lo=r.get('calibration_min_mM','');hi=r.get('calibration_max_mM','')
        if bool(lo)!=bool(hi):raise ValueError('Both calibration range limits must be supplied together')
        bounds=(float(lo),float(hi)) if lo and hi else None
        if bounds is not None and (not all(math.isfinite(v) for v in bounds) or bounds[0]<0 or bounds[1]<=bounds[0]):raise ValueError('Invalid calibration range')
        if role=='experimental':
            if bounds is None or not r.get('product_identity') or not r.get('analyst') or not r.get('recorded_at'):raise ValueError('Missing experimental metadata')
            for kind in ['raw','calibration']:
                name=r.get(kind+'_file','');p=(Path(path).parent/name).resolve()
                if not name or not p.is_relative_to(Path(path).parent.resolve()) or not p.is_file():raise ValueError('Missing local evidence')
                if hashlib.sha256(p.read_bytes()).hexdigest()!=r.get(kind+'_sha256'):raise ValueError('Evidence hash mismatch')
        v=hplc(*[float(r[k]) for k in ['area','slope_area_per_mM','intercept','dilution_factor','volume_mL','substrate_mmol']],calibration_range_mM=bounds)
        out.append(dict(sample_id=r['sample_id'],evidence_role=role,**v))
    if not out:raise ValueError('No records')
    return out

def analytical(out):
    s=source_literals('process_analytical_data');h=hplc(s['measured_peak_area'],s['calib_slope'],s['calib_intercept'],s['aliquot_dilution_factor'],s['reaction_total_volume_mL'],s['initial_substrate_mmol'])
    original_nmr=5.12*.0357/.2*100;recomputed=qnmr(5.12,3,1,3,6,168.19,1,.2)
    scenarios=[]
    for purity in [1.,.995,.99,.98]:
        for fraction in [1.,.5,.25]:
            r=qnmr(5.12,3,1,3,6,168.19,purity,.2,fraction);scenarios.append(dict(purity=purity,sample_fraction=fraction,**r))
    cwrite(out/'qnmr_scenarios.csv',scenarios)
    relaxation=[]
    for delay in [1,3,5,10,15,30]:
        ratio=(1-np.exp(-delay/3))/(1-np.exp(-delay/1));relaxation.append(dict(repetition_time_s=delay,assumed_product_T1_s=3,assumed_standard_T1_s=1,signal_ratio_bias_factor=float(ratio),apparent_yield_pct=float(recomputed['yield_pct']*ratio),relative_bias_pct=float(100*(ratio-1))))
    cwrite(out/'relaxation_scenarios.csv',relaxation)
    n=.2*original_nmr/100;counter=[dict(changed_field='area',value=(n/.006/25)*14250+120,interpretation='One-field value needed if NMR label were correct; not a correction'),dict(changed_field='dilution_factor',value=n/(h['measured_mM']*.006),interpretation='One-field value needed if NMR label were correct; not a correction'),dict(changed_field='reaction_volume_mL',value=n*1000/h['reaction_mM'],interpretation='One-field value needed if NMR label were correct; not a correction')]
    cwrite(out/'concordance_counterfactuals.csv',counter)
    row=dict(sample_id='source-example-01',evidence_role='source_example',area_mode='absolute',concentration_unit='mM',volume_unit='mL',area=245600,slope_area_per_mM=14250,intercept=120,dilution_factor=25,volume_mL=6,substrate_mmol=.2,calibration_min_mM='',calibration_max_mM='',product_identity='',analyst='',recorded_at='',raw_file='',raw_sha256='',calibration_file='',calibration_sha256='')
    cwrite(ROOT/'data/source_integrated_area.csv',[row]);ingestion=parse_integrated_csv(ROOT/'data/source_integrated_area.csv','source_example');jwrite(out/'ingestion_example.json',ingestion)
    nmrrow=dict(sample_id='source-nmr-01',evidence_role='source_example',mass_unit='mg',amount_unit='mmol',product_integral=5.12,standard_integral=3,product_H=1,standard_H=3,standard_mass_mg=6,standard_MW=168.19,standard_purity=1,substrate_mmol=.2,sample_fraction=1,product_peak_assignment='',product_identity='',standard_identity='1,3,5-trimethoxybenzene',relaxation_record='',analyst='',recorded_at='',raw_file='',raw_sha256='',standard_certificate_file='',standard_certificate_sha256='')
    cwrite(ROOT/'data/source_nmr_integrals.csv',[nmrrow]);jwrite(out/'nmr_ingestion_example.json',parse_nmr_csv(ROOT/'data/source_nmr_integrals.csv','source_example'))
    # Independent-assumption propagation, not an uncertainty claim for unavailable measurements.
    rng=np.random.default_rng(1842);N=100000
    area=rng.normal(245600,2456,N);slope=rng.normal(14250,285,N);intercept=rng.normal(120,20,N);df=rng.normal(25,.25,N);vol=rng.normal(6,.03,N);n0=rng.normal(.2,.001,N)
    yh=100*(area-intercept)/slope*df*vol/1000/n0
    yn=100*(rng.normal(5.12,.0512,N)/rng.normal(3,.015,N)*3)*(rng.normal(6,.03,N)/168.19)/n0
    uncertainty=dict(assumptions='Independent normal SD: area 1%, slope 2%, intercept 20 area units, dilution 1%, volume 0.5%, n0 0.5%; NMR product integral 1%, IS integral 0.5%, IS mass 0.5%; common n0 draws; IS purity fixed at 1.',draws=N,HPLC_percentiles_2p5_50_97p5=np.quantile(yh,[.025,.5,.975]).tolist(),NMR_percentiles_2p5_50_97p5=np.quantile(yn,[.025,.5,.975]).tolist(),HPLC_fraction_above_100=float(np.mean(yh>100)),difference_percentiles_pp=np.quantile(yh-yn,[.025,.5,.975]).tolist(),meaning='Assumption sensitivity intervals, not validated measurement uncertainty or confidence intervals')
    result=dict(HPLC=h,NMR_source_rounded_amount_yield_pct=original_nmr,NMR_mass_recomputed=recomputed,disagreement_pp=abs(h['yield_pct']-original_nmr),maximum_area_for_100pct_given_other_inputs=120+14250*(.2/6*1000/25),NMR_peak_assignment='Unverified. A C3-substituted indole does not retain the replaced C3-H; choose an independently assigned nonoverlapping product resonance.',concordance_validated=False,uncertainty_scenario=uncertainty)
    jwrite(out/'analytical_audit.json',result);return result

def gaussian(t,mu,sigma):return np.exp(-.5*((t-mu)/sigma)**2)/(sigma*np.sqrt(2*np.pi))
def deconvolve_known_shapes(t,signal,mu1,mu2,sigma1,sigma2):
    design=np.column_stack([np.ones_like(t),t-6.5,gaussian(t,mu1,sigma1),gaussian(t,mu2,sigma2)])
    beta=np.linalg.lstsq(design,signal,rcond=None)[0];resid=signal-design@beta
    variance=float(resid@resid/(len(t)-4));cov=variance*np.linalg.inv(design.T@design)
    return dict(area1=float(beta[2]),area2=float(beta[3]),area1_OLS_SE=float(np.sqrt(cov[2,2])),condition_number=float(np.linalg.cond(design)),fit=design@beta,residual=resid,baseline=beta[0]+beta[1]*(t-6.5))

def synthetic_recovery(out):
    rng=np.random.default_rng(2709);levels=np.array([0,.25,.5,.75,1,1.5,2]);x=np.repeat(levels,3);sd=np.full(len(x),100.);y=14250*x+120+rng.normal(0,sd)
    design=np.column_stack([x,np.ones_like(x)]);beta=np.linalg.lstsq(design,y,rcond=None)[0];resid=y-design@beta;cov=float(resid@resid/(len(x)-2))*np.linalg.inv(design.T@design)
    cwrite(out/'synthetic_calibration.csv',[dict(evidence_role='synthetic',concentration_mM=float(xx),area=float(yy),assumed_area_SD=100.) for xx,yy in zip(x,y)])
    t=np.linspace(5.8,7.2,1401);A1=17500.;A2=7000.;recover=[]
    for separation in [.03,.06,.13,.25]:
        for mismatch in [False,True]:
            for seed in range(10):
                width=.08 if mismatch else .07
                signal=50+8*(t-6.5)+A1*gaussian(t,6.42,width)+A2*gaussian(t,6.42+separation,.09)+np.random.default_rng(seed).normal(0,300,len(t))
                r=deconvolve_known_shapes(t,signal,6.42,6.42+separation,.07,.09)
                window=(t>=6.27)&(t<=6.57);naive=float(np.trapezoid(signal[window]-(50+8*(t[window]-6.5)),t[window]))
                recover.append(dict(evidence_role='synthetic',separation_min=separation,width_mismatch=mismatch,seed=seed,true_area1=A1,fitted_area1=r['area1'],relative_error_pct=100*(r['area1']/A1-1),area1_OLS_SE=r['area1_OLS_SE'],condition_number=r['condition_number'],window_area=naive,window_error_pct=100*(naive/A1-1)))
                if separation==.13 and not mismatch and seed==0:
                    cwrite(out/'synthetic_chromatogram.csv',[dict(evidence_role='synthetic',time_min=float(tt),signal_au=float(yy),fitted_au=float(ff),residual_au=float(rr),baseline_au=float(bb)) for tt,yy,ff,rr,bb in zip(t,signal,r['fit'],r['residual'],r['baseline'])])
    cwrite(out/'deconvolution_recovery.csv',recover)
    summaries=[]
    for mismatch in [False,True]:
        v=[r['relative_error_pct'] for r in recover if r['width_mismatch']==mismatch];summaries.append(dict(width_mismatch=mismatch,mean_absolute_area_error_pct=float(np.mean(np.abs(v))),max_absolute_area_error_pct=float(np.max(np.abs(v)))))
    result=dict(calibration_slope=float(beta[0]),calibration_intercept=float(beta[1]),calibration_covariance_slope_intercept=cov.tolist(),calibration_rows=len(x),calibration_range_mM=[0,2],deconvolution_cases=len(recover),summaries=summaries,model='Two known Gaussian peak shapes plus linear baseline; linear least squares; widths and centers assumed known; no chemical peak identification or vendor raw/FID parsing')
    jwrite(out/'synthetic_recovery.json',result);return result

def pareto_mask(values):
    """All objectives minimized; exact duplicates are both nondominated."""
    values=np.asarray(values,float)
    if values.ndim!=2 or not np.isfinite(values).all():raise ValueError('Finite objective matrix required')
    dominance=np.all(values[:,None,:]<=values[None,:,:],axis=2)&np.any(values[:,None,:]<values[None,:,:],axis=2)
    return ~np.any(dominance,axis=0)

def figure_data(out):
    rng=np.random.RandomState(42);n=80;yy=rng.uniform(40,95,n);fe=np.clip(115-.7*yy+rng.normal(0,4,n),30,95);const=2*96485.33*3.2/(3.6e6*.223);sec=const/(fe/100)+rng.normal(0,.05,n)
    threshold=(yy>80)&(fe>65)&(sec<2.5);true=pareto_mask(np.column_stack([-yy,-fe,sec]));noiseless=pareto_mask(np.column_stack([-yy,-fe,const/(fe/100)]))
    rng2=np.random.default_rng(92);count=np.zeros(n);draws=2000
    for _ in range(draws):
        perturbed=np.column_stack([-(yy+rng2.normal(0,2,n)),-(fe+rng2.normal(0,2,n)),sec+rng2.normal(0,.03,n)])
        count+=pareto_mask(perturbed)
    cwrite(out/'pareto_conditions.csv',[dict(evidence_role='synthetic',condition=i,yield_pct=float(yy[i]),FE_pct=float(fe[i]),SEC_kWh_kg=float(sec[i]),SEC_noiseless_kWh_kg=float(const/(fe[i]/100)),source_threshold=bool(threshold[i]),nondominated=bool(true[i]),nondominated_noiseless=bool(noiseless[i]),assumed_noise_front_fraction=float(count[i]/draws)) for i in range(n)])
    literals=source_literals('generate_publication_figures');pred=np.array(literals['predicted_yields'],float);mock=np.array(literals['experimental_yields'],float);names=literals['substrates'];err=pred-mock
    cwrite(out/'scope_mock.csv',[dict(evidence_role='mock_not_experimental',substrate=name,given_prediction_pct=float(p),given_mock_label_pct=float(y),error_pp=float(p-y)) for name,p,y in zip(names,pred,mock)])
    loo=[]
    for i in range(6):loo.append(dict(omitted=names[i],pearson_r=float(pearsonr(np.delete(pred,i),np.delete(mock,i)).statistic),MAE_pp=float(np.mean(np.abs(np.delete(err,i))))))
    cwrite(out/'scope_leave_one_out_descriptive.csv',loo)
    levels=literals['delta_G'];cwrite(out/'hypothetical_energy.csv',[dict(evidence_role='illustrative_not_DFT',state_index=i,source_state=literals['steps'][i].replace('\n',' '),relative_energy_kcal_mol=g) for i,g in enumerate(levels)])
    result=dict(Pareto=dict(n=n,source_threshold_count=int(threshold.sum()),nondominated_count=int(true.sum()),threshold_true_positives=int((threshold&true).sum()),threshold_false_positives=int((threshold&~true).sum()),threshold_false_negatives=int((~threshold&true).sum()),noiseless_nondominated_count=int(noiseless.sum()),assumed_noise_draws=draws,noise_model='Independent unbounded Gaussian SD: yield 2 pp, FE 2 pp, SEC 0.03 kWh/kg; sensitivity only',energy_numerator_kWh_kg=const,product_MW_assumed_kg_mol=.223,electrons_assumed=2,voltage_assumed_V=3.2),scope=dict(n=6,MAE_pp=float(np.abs(err).mean()),RMSE_pp=float(np.sqrt(np.mean(err**2))),bias_pp=float(err.mean()),pearson_r=float(pearsonr(pred,mock).statistic),R2_relative_to_mock_mean=float(1-sum(err**2)/sum((mock-mock.mean())**2)),LOO_r_range=[min(r['pearson_r'] for r in loo),max(r['pearson_r'] for r in loo)],chemical_prediction_validated=False),energy=dict(source_levels_kcal_mol=levels,adjacent_differences_kcal_mol=np.diff(levels).tolist(),assumed_addition_barrier_from_cation_kcal_mol=levels[2]-levels[1],DFT_executed=False,mechanism_validated=False,notes='Relative states lack balanced species, electron/proton reservoirs, reference potential, TS/IRC evidence and thermal provenance. Do not derive a rate or mechanism.'))
    jwrite(out/'figure_audit.json',result);return result

def budget(out):
    items=[('Reagents and solvents',1,2500),('Electrodes and cell consumables',1,1000),('Analytical standards and columns',1,2000),('HPLC injections',120,15),('NMR samples',20,100),('CPU pilot allocation',1,200),('Archiving and presentation materials',1,500)]
    rows=[dict(item=n,quantity=q,assumed_unit_CNY=p,total_CNY=q*p,evidence_role='planning_estimate_not_quote') for n,q,p in items];cwrite(out/'grant_budget.csv',rows)
    result=dict(total_CNY=sum(r['total_CNY'] for r in rows),confirmed_funding_CNY=None,proposed_duration_months=12,source_calendar_month_difference=20,source_calendar_months_inclusive=21,source_claim='2026-09 to 2028-05 called two years; dates span 20 months or 21 named calendar months',next_call_deadline='unverified',approval_or_grant_submission=False);jwrite(out/'grant_planning.json',result);return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'results/audit');p.add_argument('--input-csv',type=Path);p.add_argument('--assay',choices=['hplc','nmr'],default='hplc');p.add_argument('--role',choices=['source_example','synthetic','experimental'],default='source_example');a=p.parse_args();out=a.output;out.mkdir(parents=True,exist_ok=True)
    if a.input_csv:
        parser=parse_integrated_csv if a.assay=='hplc' else parse_nmr_csv
        jwrite(out/'ingested_quantification.json',dict(source_sha256=hashlib.sha256(a.input_csv.read_bytes()).hexdigest(),assay=a.assay,records=parser(a.input_csv,a.role)));print('Saved pre-integrated quantification; no measurement certification.');return
    an=analytical(out);syn=synthetic_recovery(out);fig=figure_data(out);plan=budget(out)
    summary=dict(date='2026-09-28',evidence='Source arithmetic, synthetic recovery and planning estimates; no measured chromatograms, NMR spectra or DFT outputs',analytical=an,synthetic_recovery=syn,figures=fig,grant=plan,versions={n:importlib.metadata.version(n) for n in ['numpy','scipy','pandas','matplotlib','rdkit']})
    jwrite(out/'audit_summary.json',summary);print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
