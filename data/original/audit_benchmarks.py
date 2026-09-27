"""Independent arithmetic and source-code audit; does not edit original results.

Run beside benchmark_results.json and simulate_all_topics.py.
Outputs benchmark_audit.json and benchmark_overview.png/svg.
No electronic-structure calculations or experiments are performed.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import odeint
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
raw = (ROOT/'benchmark_results.json').read_bytes()
d = json.loads(raw)
t1, t2, t3, t4, t5 = d.values()
cond = t1['optimal_conditions']
j, c, eps = cond.values()
noiseless = 92 - .25*(j-16.5)**2 - 400*(c-.15)**2 - .08*(eps-37.5)**2
np.random.seed(42)
bounds = np.array([[5.,30.],[.05,.30],[20.,60.]])
x0 = np.random.uniform(bounds[:,0], bounds[:,1], size=(8,3))
y0 = np.clip(92 - .25*(x0[:,0]-16.5)**2 - 400*(x0[:,1]-.15)**2 - .08*(x0[:,2]-37.5)**2 + np.random.normal(0,1.5,8),0,99)
np.random.seed(101)
train = np.random.randn(60,4)
sites = np.array([[-5.8,1.42,.88,22.4],[-5.8,.95,.54,18.1],[-5.8,.32,.12,35.8]])
np.random.seed(7)
class_x = np.random.uniform(.1,2.5,(80,4))
class_y = ((class_x[:,0]<1.2)&(class_x[:,1]<1.8)&(class_x[:,2]>.8)).astype(int)
rank = sorted(t3['screened_configurations'], key=lambda x:x['activation_barrier_kcal_mol'])
gap = rank[1]['activation_barrier_kcal_mol'] - rank[0]['activation_barrier_kcal_mol']
flow = []
for q in [.2,.5,1.,2.]:
    analytic = 100*(1-np.exp(-3.5/q))
    numeric = 100*(1-odeint(lambda C,z: -3.5/q*C,[1.],np.linspace(0,1,50))[-1,0])
    entry = t5['flow_rate_sweep'][f'flow_{q}_mL_min']
    assert abs(numeric - analytic)<1e-4
    assert entry['conversion_pct'] == round(float(numeric),2)
    assert entry['sty_relative_val'] == round(float(numeric*q*.45),2)
    flow.append({'q_numeric':q, 'analytic_conversion_pct':float(analytic), 'odeint_conversion_pct':float(numeric), 'absolute_error_percentage_points':float(abs(numeric-analytic)), **entry})
audit = {
 'evidence_scope':'synthetic benchmarks and code/arithmetic audit only',
 'original_json_sha256':hashlib.sha256(raw).hexdigest(),
 'topic1':{'n_initial':8,'n_sequential':12,'n_total':20,'initial_best':float(y0.max()),'noiseless_yield_at_reported_best':float(noiseless),'noise_residual_at_reported_best':float(t1['final_best_yield']-noiseless),'noiseless_global_maximum':92.0,'simulated_noise_standard_deviation':1.5,'gp_alpha_variance':.01,'normalized_input':False,'normalized_target':False,'independent_random_search_baseline_executed':False},
 'topic2':{'training_size':60,'training_feature_min':train.min(axis=0).tolist(),'training_feature_max':train.max(axis=0).tolist(),'sites_outside_training_range':((sites<train.min(axis=0))|(sites>train.max(axis=0))).tolist(),'score_scope':'in-sample R2','site_C3_label_hardcoded':True,'real_DFT_or_CV_executed':False},
 'topic3':{'lowest':rank[0],'second_lowest':rank[1],'gap_kcal_mol':float(gap),'assigned_noise_sd_kcal_mol':.4,'hypothetical_eyring_ratio_at_298K':float(np.exp(gap/(.00198720425864083*298))),'real_transition_states':0},
 'topic4':{'training_size':80,'synthetic_positive_count':int(class_y.sum()),'synthetic_negative_count':int((class_y==0).sum()),'score_scope':'training accuracy','feasibility_score_hardcoded':True,'risk_statement_hardcoded':True,'prospective_validation_executed':False},
 'topic5':{'flow_checks':flow,'max_proxy_flow_in_tested_grid':2.0,'units_in_original_ODE_consistent':False,'physical_reactor_volume_or_product_molecular_weight_provided':False,'Butler_Volmer_implemented':False,'independent_mass_transfer_model_implemented':False,'original_summary_inconsistent':True},
 'original_json_unchanged': hashlib.sha256((ROOT/'benchmark_results.json').read_bytes()).hexdigest()==hashlib.sha256(raw).hexdigest()
}
(ROOT/'benchmark_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'#ffffff','axes.titleweight':'bold','axes.labelcolor':'#263548','text.color':'#263548','axes.edgecolor':'#94a3b8','svg.fonttype':'none'})
fig, ax = plt.subplots(2,2,figsize=(12,8.2),layout='constrained')
fig.suptitle('Executed synthetic benchmarks | No experimental or DFT validation',fontsize=15,fontweight='bold')
rounds=[0]+[r['round'] for r in t1['bo_trajectory']]
bests=[float(y0.max())]+[r['best_yield'] for r in t1['bo_trajectory']]
ax[0,0].plot(rounds,bests,marker='o',color='#087f8c',lw=2)
ax[0,0].axhline(92,color='#bb7824',ls='--',label='Noiseless surface maximum = 92')
ax[0,0].set(title='A. BO: best noisy observation',xlabel='Sequential round (0 = initial 8 points)',ylabel='Synthetic yield (%)',ylim=(80,98))
ax[0,0].legend(fontsize=8,loc='lower right')
names=[r['metal']+'-'+r['coordination'] for r in rank]
vals=[r['activation_barrier_kcal_mol'] for r in rank]
ax[0,1].barh(names,vals,color=['#087f8c']+['#c2d8df']*11)
ax[0,1].invert_yaxis()
ax[0,1].set(title='B. SAC: empirical mock barriers',xlabel='Assigned barrier (kcal/mol)',xlim=(0,20.2))
for i,v in enumerate(vals):ax[0,1].text(v+.1,i,f'{v:.2f}',va='center',fontsize=8)
qs=np.array([.2,.5,1.,2.]); conv=np.array([x['conversion_pct'] for x in flow]); proxy=np.array([x['sty_relative_val'] for x in flow])
ax[1,0].plot(qs,conv,'o-',color='#087f8c',lw=2)
ax[1,0].set(title='C. ODE: conversion decreases with q',xlabel='q (input labelled mL/min; units unresolved)',ylabel='Conversion (%)',ylim=(0,110))
ax[1,1].plot(qs,proxy,'s-',color='#bb7824',lw=2)
ax[1,1].set(title='D. Throughput proxy increases with q',xlabel='q (input labelled mL/min; units unresolved)',ylabel='Relative STY proxy (arbitrary units)',ylim=(0,90))
for a,ys in [(ax[1,0],conv),(ax[1,1],proxy)]:
    for x,y in zip(qs,ys):a.annotate(f'{y:.2f}',(x,y),xytext=(0,7),textcoords='offset points',ha='center',fontsize=9)
for a in ax.ravel():a.grid(axis='y' if a is not ax[0,1] else 'x',alpha=.15);a.set_axisbelow(True)
fig.savefig(ROOT/'benchmark_overview.png',dpi=180)
fig.savefig(ROOT/'benchmark_overview.svg')
print(json.dumps(audit,ensure_ascii=False,indent=2))
