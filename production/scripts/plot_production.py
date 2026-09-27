"""Publication figures from recorded results; does not rerun models."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
ROOT=Path(__file__).resolve().parents[1];A=ROOT/'results/audit';OUT=ROOT/'reports/figures'
def main():
    p=argparse.ArgumentParser();p.add_argument('--zh-font');args=p.parse_args()
    family='DejaVu Sans'
    if args.zh_font:font_manager.fontManager.addfont(args.zh_font);family=font_manager.FontProperties(fname=args.zh_font).get_name()
    bo=pd.read_csv(A/'task_a_paired_runs.csv');sites=pd.read_csv(A/'task_b_site_summary.csv');cv=pd.read_csv(A/'task_c_summary.csv');flow=pd.read_csv(A/'task_d_eta_sweep.csv')
    colors=['#286C8D','#A978AD','#DC8A36','#318573']
    for zh in [False,True]:
        suf='_zh' if zh else '';plt.rcParams.update({'font.family':family if zh else 'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'path' if zh else 'none'})
        fig,axs=plt.subplots(2,2,figsize=(14,9.5),layout='constrained')
        methods=['random','source_gp','scaled_physical_gp','scaled_onehot_gp'];labels=['随机','原始 GP','缩放物性 GP','缩放独热 GP'] if zh else ['Random','Source GP','Scaled physical','Scaled one-hot']
        axs[0,0].boxplot([bo[bo.method==x].recommendation_regret for x in methods],tick_labels=labels)
        axs[0,0].set(title='A  配对 BO 检验：30 种子 × 14 次评价' if zh else 'A  Paired BO: 30 seeds × 14 evaluations',ylabel='合成推荐遗憾值；越低越好' if zh else 'Synthetic recommendation regret; lower is better')
        x=np.arange(len(sites));labels=[r.site_label if r.site_label.startswith('C') else f'Ph:{r.atom_index}' for r in sites.itertuples()]
        axs[0,1].bar(x,sites.score_mean,yerr=sites.score_sd,color=['#DC8A36' if n=='C3' else '#318573' for n in labels],capsize=3)
        axs[0,1].set(xticks=x,xticklabels=labels,title='B  32 个构象的启发式评分；均值 ± SD' if zh else 'B  Heuristic score over 32 conformers; mean ± SD',ylabel='未校准评分；非反应概率' if zh else 'Uncalibrated score; not reaction probability')
        ax=axs[1,0]
        for i,split in enumerate(['metal','coordination_type']):
            sub=cv[cv.split==split].set_index('model').loc[['mean','ridge','extra_trees']]
            ax.bar(np.arange(3)+(i-.5)*.32,sub.test_mae,width=.32,color=colors[i],label=('留金属' if i==0 else '留配位类型') if zh else ('Leave metal out' if i==0 else 'Leave coordination out'))
        ax.set(xticks=range(3),xticklabels=['均值','岭回归','极端随机树'] if zh else ['Mean','Ridge','Extra Trees'],title='C  SAC 分组留出误差；人工标签' if zh else 'C  SAC group holdouts; synthetic labels',ylabel='宏平均 MAE（源评分单位）' if zh else 'Macro MAE (source score units)');ax.legend(fontsize=9)
        ax=axs[1,1]
        for i,(q,g) in enumerate(flow.groupby('q_uL_min')):ax.plot(g.eta_V,g.conversion*100,'o-',color=colors[i],label=f'{q} µL/min')
        ax.set(title='D  过电位情景扫描；固定假设参数' if zh else 'D  Overpotential scenarios; prescribed parameters',xlabel='过电位 (V)' if zh else 'Overpotential (V)',ylabel='模型反应物转化率 (%)' if zh else 'Modeled reactant conversion (%)');ax.legend(fontsize=8)
        fig.suptitle('四任务定量核验：保留负结果与模型边界' if zh else 'Four-task quantitative audit: negative results and model boundaries',fontsize=17)
        for ext in ['png','svg']:fig.savefig(OUT/f'production_audit{suf}.{ext}',dpi=180)
        plt.close(fig)
        fig,axs=plt.subplots(1,2,figsize=(13,4.7),layout='constrained')
        rows=pd.read_csv(A/'task_d_analytic_audit.csv')
        axs[0].plot(rows.q_uL_min,rows.conversion*100,'o-',color=colors[0]);axs[0].set(xlabel='流量 (µL/min)' if zh else 'Flow (µL/min)',ylabel='转化率 (%)' if zh else 'Conversion (%)',title='A  停留时间下降，转化率下降' if zh else 'A  Less residence time, lower conversion')
        axs[1].plot(rows.q_uL_min,rows.STY_mmol_L_h,'s-',color=colors[3]);axs[1].set(xlabel='流量 (µL/min)' if zh else 'Flow (µL/min)',ylabel='反应物消耗 STY (mmol/L/h)' if zh else 'Reactant-disappearance STY (mmol/L/h)',title='B  消耗通量增加；产物选择性未建模' if zh else 'B  Higher consumption flux; no product selectivity model')
        fig.suptitle('传质／动力学约化 ODE；η = 0.45 V' if zh else 'Reduced transport/kinetic ODE; η = 0.45 V',fontsize=15)
        for ext in ['png','svg']:fig.savefig(OUT/f'flow_tradeoff{suf}.{ext}',dpi=180)
        plt.close(fig)
        tasks=[('基础与记录','Foundations and records',1,12),('CV 与分析校准','CV and analytical calibration',4,15),('实测数据与对照 BO','Measured data and controlled BO',10,24),('结构描述符及核验','Descriptors and validation',7,30),('催化剂／机理联合课题','Catalyst/mechanism project',19,36),('条件满足后开展流动','Flow after evidence gates',28,42),('竞赛与大创材料','Competition/innovation dossier',16,36),('论文、答辩与复现归档','Manuscript, thesis and reproducibility',31,48)]
        fig,ax=plt.subplots(figsize=(14,5.1),layout='constrained')
        for i,(cn,en,start,end) in enumerate(tasks):ax.barh(i,end-start+1,left=start-.5,color=colors[i%4],height=.6)
        ax.set(yticks=range(len(tasks)),yticklabels=[t[0] if zh else t[1] for t in tasks],xticks=[1,6,12,18,24,30,36,42,48],xlim=(.5,48.5),xlabel='自项目启动起的月份；进阶取决于证据门槛' if zh else 'Months from project start; progression requires evidence gates')
        ax.invert_yaxis();ax.grid(axis='x',alpha=.2);ax.set_axisbelow(True)
        for month in [12.5,24.5,36.5]:ax.axvline(month,color='gray',lw=.8,ls='--')
        fig.suptitle('本科四年实施甘特图：拟议计划，不保证竞赛或发表结果' if zh else 'Four-year undergraduate Gantt: proposed plan, not guaranteed awards or publication',fontsize=15)
        for ext in ['png','svg']:fig.savefig(OUT/f'undergraduate_gantt{suf}.{ext}',dpi=180)
        plt.close(fig)
if __name__=='__main__':main()
