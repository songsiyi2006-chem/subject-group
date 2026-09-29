"""Post hoc arithmetic and bilingual plots of archived, executed model studies.

No quantum calculation, training, PDE solve, molecular dynamics or experiment
is launched. Each input is identified by its file hash and frozen commit.
"""
from pathlib import Path
import csv, hashlib, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'manuscript/expanded'
RESULTS=OUT/'results'; FIG=OUT/'figures'
INPUTS={}
COMMIT='073bb3872d9112ec1ca46b3b9e80476d00ffbf97'

def read(name):
    p=ROOT/name;INPUTS[name]=hashlib.sha256(p.read_bytes()).hexdigest()
    if p.suffix=='.json':return json.loads(p.read_text('utf-8'))
    with p.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))

def savecsv(name,rows):
    with (RESULTS/name).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main(only_figure=None):
    RESULTS.mkdir(parents=True,exist_ok=True);FIG.mkdir(exist_ok=True)
    metrics=read('electrograph/results/learning/regression_metrics.csv')
    net=read('electratwin/results/reaction_network/study_summary.json')['baseline']
    budget=read('electratwin/results/benchmark_extension/budget_prefixes.csv')
    holdout=read('electratwin/results/benchmark_extension/holdout_grouped_metrics.csv')
    kinetic=read('electrograph/results/kinetics/potential_summary.csv')
    thermostat=read('synthapore/results/dynamics/thermostat_summary.csv')
    nve=read('synthapore/results/dynamics/nve_timestep_summary.csv')
    analytical=read('toolkit/results/audit/deconvolution_recovery.csv')
    hydrate=[r for r in metrics if r['split']=='test']
    savecsv('hydration_test_metrics.csv',hydrate)
    methods=['gp_mc_ehvi','random','maximin_spacefill']
    budget_keys=[(r['method'],int(r['seed']),int(r['evaluation'])) for r in budget]
    assert len(budget_keys)==len(set(budget_keys)), 'Duplicate policy/seed/budget record'
    final={m:{int(r['seed']):float(r['full_pool_HV_fraction']) for r in budget if r['method']==m and int(r['evaluation'])==25} for m in methods}
    assert all(len(v)==64 for v in final.values())
    hvs=[]
    for m in methods:
        a=np.array([final[m][s] for s in sorted(final[m])])
        hvs.append({'method':m,'seeds':len(a),'mean_HV_fraction':float(a.mean()),'sd_across_seeds':float(a.std(ddof=1)),'min_HV_fraction':float(a.min()),'max_HV_fraction':float(a.max())})
    savecsv('optimization_final.csv',hvs)
    paired=[]
    for other in methods[1:]:
        d=np.array([final[methods[0]][s]-final[other][s] for s in range(64)])
        paired.append({'left':methods[0],'right':other,'n':64,'mean_fraction_difference':float(d.mean()),'wins':int((d>1e-12).sum()),'losses':int((d< -1e-12).sum()),'ties':int((np.abs(d)<=1e-12).sum())})
    savecsv('optimization_paired.csv',paired)
    recover=[]
    for mismatch in ['False','True']:
        for separation in sorted({float(r['separation_min']) for r in analytical}):
            g=[r for r in analytical if r['width_mismatch']==mismatch and float(r['separation_min'])==separation]
            recover.append({'width_mismatch':mismatch,'separation_min':separation,'n':len(g),'mean_fitted_error_pct':float(np.mean([float(r['relative_error_pct']) for r in g])),'mean_absolute_fitted_error_pct':float(np.mean([abs(float(r['relative_error_pct'])) for r in g])),'mean_absolute_window_error_pct':float(np.mean([abs(float(r['window_error_pct'])) for r in g]))})
    savecsv('chromatography_grouped.csv',recover)
    savecsv('thermostat_comparison.csv',thermostat)
    savecsv('kinetic_potential_comparison.csv',kinetic)
    networkrows=[{'observable':k,'value':net[k],'unit':'percent'} for k in ['conversion_A_pct','net_P_yield_pct','net_P_selectivity_pct','net_P_faradaic_efficiency_pct','gross_AP_charge_fraction_pct']]
    savecsv('network_observables.csv',networkrows)
    available={x.name for x in font_manager.fontManager.ttflist}
    chinese=next((n for n in ['Microsoft YaHei','SimHei','Noto Sans CJK SC','SimSun'] if n in available),None)
    if not chinese:raise RuntimeError('A CJK font is required')
    assets=[]
    for lang in ['english','chinese']:
        zh=lang=='chinese'
        plt.rcParams.update({'font.family':['DejaVu Sans',chinese] if zh else ['DejaVu Sans'],'font.size':9,'axes.unicode_minus':False,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
        def finish(fig,name):
            if only_figure is not None and name != only_figure:
                for ext in ['png','svg']:
                    p=FIG/f'{name}_{lang}.{ext}'
                    assets.append({'path':p.relative_to(ROOT).as_posix(),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
                plt.close(fig)
                return
            fig.tight_layout(pad=1.3)
            for ext in ['png','svg']:
                p=FIG/f'{name}_{lang}.{ext}';fig.savefig(p,dpi=300,bbox_inches='tight');assets.append({'path':p.relative_to(ROOT).as_posix(),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
            plt.close(fig)
        fig,ax=plt.subplots(figsize=(9,4.1));labels=[r['model'].replace('mpnn_seed','MPNN ').replace('shuffled_train_','Shuffled ') for r in hydrate]
        ax.barh(labels,[float(r['RMSE_kcal_mol']) for r in hydrate],color=['#24799a']*3+['#bc7044']+['#565d70']*(len(hydrate)-4))
        ax.set_xlabel('测试集 RMSE / kcal mol⁻¹' if zh else 'Test RMSE / kcal mol⁻¹');ax.set_title('FreeSolv 实测标签：三个种子、打乱标签对照与简单基线' if zh else 'Measured FreeSolv labels: three seeds, shuffled-label control and baseline')
        finish(fig,'M1_Hydration')
        fig,axes=plt.subplots(1,2,figsize=(9,4.1))
        axes[0].bar(list(net['outlet_mM']),list(net['outlet_mM'].values()),color=['#555d70','#24799a','#ce9850','#bb594e']);axes[0].set_ylabel('出口浓度 / mM' if zh else 'Outlet concentration / mM');axes[0].set_title('a  四物种模型' if zh else 'a  Four-species model')
        cats=['A 转化率','P 净产率','P 选择性','P 净 FE','A→P 毛电荷占比'] if zh else ['A conversion','P net yield','P selectivity','P net FE','Gross A→P charge']
        axes[1].barh(cats,[r['value'] for r in networkrows],color='#24799a');axes[1].set_xlabel('%');axes[1].set_title('b  不同分母，不可互换' if zh else 'b  Different denominators, different observables')
        finish(fig,'M2_Network')
        fig,axes=plt.subplots(1,2,figsize=(9,4.2))
        for i,m in enumerate(methods):
            a=np.array([final[m][s] for s in range(64)]);axes[0].scatter(np.full(64,i)+np.linspace(-.12,.12,64),a,s=12,alpha=.65);axes[0].plot([i-.18,i+.18],[a.mean()]*2,color='black',lw=2)
        axes[0].set_xticks(range(3),['GP-EHVI','Random','Maximin']);axes[0].set_ylabel('最终 HV / 候选池 HV' if zh else 'Final HV / pool HV');axes[0].set_title('a  64 个配对种子，预算 25' if zh else 'a  64 paired seeds; budget 25')
        # The archived names are used directly, never silently substituted.
        selected=[r for r in holdout if r['target']=='STY_div_50000' and r['split_kind'].startswith('block')]
        short_models={'fixed_matern_gp':'GP','training_mean':'均值' if zh else 'Mean','quadratic_regression':'二次' if zh else 'Quad.'}
        labs=[short_models[r['model']] for r in selected]
        positions=[0,1,2,4,5,6]
        axes[1].bar(positions,[float(r['RMSE']) for r in selected],color='#bc7044');axes[1].set_xticks(positions,labs,fontsize=9.5);axes[1].set_ylabel('STY / 50000 的 RMSE' if zh else 'RMSE of STY / 50000');axes[1].set_title('b  分块留出误差' if zh else 'b  Blocked-input holdouts')
        axes[1].text(1,-.17,'η 分块' if zh else 'η block',ha='center',transform=axes[1].get_xaxis_transform(),fontsize=10)
        axes[1].text(5,-.17,'流量分块' if zh else 'Flow block',ha='center',transform=axes[1].get_xaxis_transform(),fontsize=10)
        finish(fig,'M3_Optimization')
        fig,axes=plt.subplots(1,2,figsize=(9,4.1));eta=np.array([float(r['eta_V']) for r in kinetic]);means=np.array([float(r['mean_TOF_s_1']) for r in kinetic]);reference=np.array([float(r['finite_window_CTMC_TOF_s_1']) for r in kinetic])
        axes[0].errorbar(eta,means,yerr=[float(r['sd_TOF_s_1']) for r in kinetic],fmt='o',label='SSA ± SD');axes[0].plot(eta,reference,label='Finite-window CTMC');axes[0].set_xlabel('过电位 / V' if zh else 'Overpotential / V');axes[0].set_ylabel('TOF / s⁻¹');axes[0].legend(fontsize=8,loc='lower right')
        states=['empty','substrate','radical','intermediate','product'];axes[1].stackplot(eta,*[[float(r['steady_'+s]) for r in kinetic] for s in states],labels=states);axes[1].set_xlabel('过电位 / V' if zh else 'Overpotential / V');axes[1].set_ylabel('稳态占据概率' if zh else 'Steady occupation probability');axes[1].legend(loc='center left',bbox_to_anchor=(1,.5),fontsize=7)
        fig.suptitle('执行的抽象动力学模型；非校准反应速率' if zh else 'Executed abstract kinetic model; rates are not chemically calibrated',fontsize=10);finish(fig,'M4_Kinetics')
        fig,axes=plt.subplots(1,2,figsize=(9,4.1));names=['BAOAB','Berendsen','Source DOF']
        axes[0].bar(names,[float(r['replicate_mean_temperature_K']) for r in thermostat],color='#24799a');axes[0].axhline(300,color='black',ls='--');axes[0].set_ylabel('平均温度 / K' if zh else 'Mean temperature / K');axes[0].tick_params(axis='x',rotation=15)
        axes[1].bar(names,[float(r['variance_ratio_to_canonical']) for r in thermostat],color='#bc7044');axes[1].set_yscale('log');axes[1].axhline(1,color='black',ls='--');axes[1].set_ylabel('温度方差 / 正则参考方差' if zh else 'Temperature variance / canonical reference');axes[1].tick_params(axis='x',rotation=15)
        fig.suptitle('解析谐振簇：平均值接近不代表分布正确' if zh else 'Analytic harmonic cluster: a correct mean does not establish the distribution',fontsize=10);finish(fig,'M5_Thermostats')
        fig,axes=plt.subplots(1,2,figsize=(9,4.1))
        for mismatch,ax in zip(['False','True'],axes):
            rows=[r for r in recover if r['width_mismatch']==mismatch];x=[r['separation_min'] for r in rows]
            ax.plot(x,[r['mean_absolute_fitted_error_pct'] for r in rows],'o-',label='Fitted peak area');ax.plot(x,[r['mean_absolute_window_error_pct'] for r in rows],'s--',label='Window integral');ax.set_xlabel('峰间距 / min' if zh else 'Peak separation / min');ax.set_ylabel('平均绝对面积误差 / %' if zh else 'Mean absolute area error / %');ax.set_title(('峰宽正确' if mismatch=='False' else '峰宽误设') if zh else ('Correct width' if mismatch=='False' else 'Misspecified width'));ax.legend(fontsize=8)
        fig.suptitle('80 个合成色谱案例；不是实际仪器数据' if zh else '80 synthetic chromatographic cases; not measured instrument data',fontsize=10);finish(fig,'M6_Analytical')
    result={'scope':'Retrospective arithmetic and plots of executed saved studies; synthetic/model results remain labeled.','source_commit':COMMIT,'new_quantum_jobs':0,'new_training_runs':0,'new_PDE_solves':0,'new_experiments':0,'input_sha256':INPUTS,'counts':{'hydration_test_models':len(hydrate),'optimization_seed_policy_pairs':192,'optimization_paired_comparisons':len(paired),'chromatography_cases':len(analytical),'chromatography_groups':len(recover),'thermostat_summary_rows':len(thermostat),'NVE_trajectories':len(nve),'kinetic_potential_rows':len(kinetic)},'optimization':hvs,'paired':paired,'analytical':recover,'figures':assets,'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (RESULTS/'multiscale_analysis.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'figures':len(assets),'counts':result['counts'],'paired':paired},ensure_ascii=False))

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--only-figure',choices=['M1_Hydration','M2_Network','M3_Optimization','M4_Kinetics','M5_Thermostats','M6_Analytical'])
    args=parser.parse_args()
    main(args.only_figure)
