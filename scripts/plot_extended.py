"""Render figures from saved results; does not rerun any calculation."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--zh-font',default=None);args=parser.parse_args()
    if args.zh_font:
        from matplotlib import font_manager
        font_manager.fontManager.addfont(args.zh_font)
        zh_family=font_manager.FontProperties(fname=args.zh_font).get_name()
    else:zh_family='DejaVu Sans'
    out=ROOT/'reports/figures';out.mkdir(exist_ok=True)
    ext=ROOT/'results/extended'
    bo=pd.read_csv(ext/'bo_runs.csv');reg=pd.read_csv(ext/'regression_runs.csv')
    sac=pd.read_csv(ext/'sac_rank_sensitivity.csv');cls=pd.read_csv(ext/'classification_runs.csv')
    flow=pd.read_csv(ext/'flow_grid.csv');cal=pd.read_csv(ext/'classification_reliability_seed0.csv')
    colors=['#496B83','#EB9B38','#288B80']
    for zh in [False,True]:
        plt.rcParams.update({'font.family':zh_family if zh else 'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'path' if zh else 'none'})
        fig,axs=plt.subplots(2,3,figsize=(16,9.4),layout='constrained')
        ax=axs[0,0];methods=['random','original_gp','scaled_gp'];labels=['随机','原始 GP','缩放 GP'] if zh else ['Random','Original GP','Scaled GP']
        ax.boxplot([bo[bo.method==m].recommendation_regret for m in methods],tick_labels=labels,showfliers=True)
        ax.set(title='A  30 个配对种子的推荐遗憾值' if zh else 'A  Recommendation regret, 30 paired seeds',ylabel='合成产率百分点；越低越好' if zh else 'Synthetic yield points; lower is better')
        ax=axs[0,1]
        for m,c,label in zip(['mean','gradient_boosting','linear'],colors,['均值基线','梯度提升','线性回归'] if zh else ['Mean baseline','Gradient boosting','Linear regression']):
            g=reg[reg.model==m].groupby('n_train').test_rmse
            avg=g.mean();sd=g.std();ax.plot(avg.index,avg.values,'o-',color=c,label=label);ax.fill_between(avg.index,avg-sd,avg+sd,color=c,alpha=.15)
        ax.set(title='B  独立测试集学习曲线' if zh else 'B  Independent-test learning curves',xlabel='训练样本数' if zh else 'Training rows',ylabel='合成目标 RMSE；均值 ± SD' if zh else 'Synthetic target RMSE; mean ± SD',yscale='log');ax.legend(fontsize=9)
        ax=axs[0,2];sub=sac[sac.noise_sd==.4].sort_values('probability_rank1')
        ax.barh(sub.site,sub.probability_rank1*100,color='#288B80');ax.set(title='C  SAC 生成器重复抽样' if zh else 'C  SAC generator resampling',xlabel='排名第一的概率 (%)；非实验排序' if zh else 'Rank-1 probability (%); not experimental')
        ax=axs[1,0];groups=[('prevalence','all'),('random_forest','all'),('random_forest','near_boundary')]
        means=[];sds=[]
        for m,s in groups:
            v=cls[(cls.model==m)&(cls.subset==s)].balanced_accuracy;means.append(v.mean());sds.append(v.std())
        ax.bar(range(3),means,yerr=sds,color=colors,capsize=4)
        ax.set(xticks=range(3),xticklabels=['基线','森林：全部','森林：边界附近'] if zh else ['Baseline','Forest: all','Forest: near boundary'],ylim=(0,1.05),title='D  分类器独立测试；均值 ± SD' if zh else 'D  Held-out classification; mean ± SD',ylabel='平衡准确率' if zh else 'Balanced accuracy')
        ax=axs[1,1];ax.plot([0,1],[0,1],'--',color='gray',label='理想参照' if zh else 'Ideal reference')
        ax.plot(cal.mean_probability,cal.observed_frequency,'o-',color='#288B80',label='种子 0；未校准' if zh else 'Seed 0; uncalibrated')
        ax.set(title='E  可靠性诊断；仅作描述' if zh else 'E  Reliability diagnostic; descriptive only',xlabel='分箱平均预测概率' if zh else 'Mean predicted probability in bin',ylabel='实测合成标签频率' if zh else 'Observed synthetic-label frequency',xlim=(0,1),ylim=(0,1.04));ax.legend(fontsize=9)
        ax=axs[1,2];good=flow.charge_feasible.astype(bool)
        ax.plot(flow.q_ml_min,flow.sty_g_L_h,color='#A4ADB5',linestyle='--',label='动力学假设值' if zh else 'Kinetic scenario')
        ax.plot(flow[good].q_ml_min,flow[good].sty_g_L_h,color='#288B80',label='仅电荷约束可行' if zh else 'Charge constraint satisfied')
        eligible=good&(flow.conversion>=.8)&(flow.selectivity>=.8)
        ax.plot(flow[eligible].q_ml_min,flow[eligible].sty_g_L_h,lw=5,color='#EB9B38',label='同时满足 X、S ≥ 0.8' if zh else 'Also X, S ≥ 0.8')
        ax.set(title='F  流动模型；假设参数' if zh else 'F  Flow model; assumed parameters',xlabel='流量 (mL/min)' if zh else 'Flow (mL/min)',ylabel='假设 STY (g/L/h)' if zh else 'Scenario STY (g/L/h)');ax.legend(fontsize=8)
        fig.suptitle('扩展数值验证：合成数据与未校准场景' if zh else 'Extended numerical validation: synthetic data and uncalibrated scenarios',fontsize=17)
        suffix='_zh' if zh else ''
        for extn in ['png','svg']:fig.savefig(out/f'extended_validation{suffix}.{extn}',dpi=190)
        plt.close(fig)
        mol=pd.read_csv(ROOT/'results/molecular/molecular_summary.csv')
        fig,axs=plt.subplots(1,2,figsize=(13,5.2),layout='constrained')
        names=['苯','吡啶','苯甲醚','吲哚','N-甲基吲哚','苯并呋喃'] if zh else ['Benzene','Pyridine','Anisole','Indole','N-methylindole','Benzofuran']
        for col,label,c in [('fixed_nuclei_removal_gap_eV','乙腈 ALPB' if zh else 'Acetonitrile ALPB',colors[2]),('gas_fixed_nuclei_removal_gap_eV','气相；相同核坐标' if zh else 'Gas; same nuclei',colors[1])]:
            axs[0].plot(mol[col],names,'o-',color=c,label=label)
        axs[0].set(xlabel='模型电荷移除能量差 (eV)' if zh else 'Model charge-removal energy difference (eV)',title='A  固定核坐标的溶剂敏感性' if zh else 'A  Solvent sensitivity at fixed nuclei');axs[0].legend(fontsize=9)
        axs[1].barh(names,mol.TPSA_A2,color=colors[0]);axs[1].set(xlabel='RDKit TPSA (Å²)',title='B  拓扑极性表面积' if zh else 'B  Topological polar surface area')
        fig.suptitle('真实分子的 GFN2-xTB / RDKit 计算；非实验氧化电位' if zh else 'Executed GFN2-xTB / RDKit descriptors; not experimental oxidation potentials',fontsize=14)
        for extn in ['png','svg']:fig.savefig(out/f'molecular_descriptors{suffix}.{extn}',dpi=190)
        plt.close(fig)

if __name__=='__main__':main()
