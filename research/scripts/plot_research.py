"""Create inspectable bilingual scientific figures from committed audit data."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'results/audit';FIG=ROOT/'reports/figures'
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--zh-font',required=True);args=parser.parse_args()
    font_manager.fontManager.addfont(args.zh_font);zhfont=font_manager.FontProperties(fname=args.zh_font).get_name()
    grid=pd.read_csv(DATA/'conditions_source_rounded.csv');front=grid[grid.pareto].sort_values('j_mA_cm2');best=front.loc[front.Green_Score.idxmax()]
    transport=pd.read_csv(DATA/'transport_audit.csv')
    for zh in [False,True]:
        plt.rcParams.update({'font.family':zhfont if zh else 'DejaVu Sans','font.size':11,'axes.unicode_minus':False,'svg.fonttype':'path'})
        suffix='_zh' if zh else ''
        def save(fig,name):
            fig.savefig(FIG/(name+suffix+'.png'),dpi=180);fig.savefig(FIG/(name+suffix+'.svg'));plt.close(fig)
        fig,ax=plt.subplots(1,2,figsize=(13,4.8),layout='constrained')
        p=ax[0].scatter(grid.Yield_pct,grid.FE_pct,c=grid.SEC_kWh_kg,cmap='viridis_r',s=24,alpha=.55)
        ax[0].plot(front.Yield_pct,front.FE_pct,'o-',color='#db8126',label='12 条前沿记录 / 10 个目标组合' if zh else '12 frontier rows / 10 objective vectors')
        ax[0].scatter([best.Yield_pct],[best.FE_pct],marker='*',s=190,color='#b52138',zorder=5)
        ax[0].set(xlabel='合成产率 (%)' if zh else 'Synthetic yield (%)',ylabel='给定 FE (%)' if zh else 'Assigned FE (%)');ax[0].legend(fontsize=9)
        fig.colorbar(p,ax=ax[0],label='SEC (kWh/kg)')
        ax[1].plot(front.j_mA_cm2,front.Green_Score,'o-',color='#226c8c');ax[1].scatter([best.j_mA_cm2],[best.Green_Score],marker='*',s=190,color='#b52138',zorder=5)
        ax[1].set(xlabel='j (mA/cm²)',ylabel='Y × FE / SEC（指定权重）' if zh else 'Y × FE / SEC (assigned weighting)')
        fig.suptitle('人工目标曲面：帕累托前沿与评分折中' if zh else 'Assigned objective surface: Pareto frontier and score-based compromise');save(fig,'pareto_tradeoffs')
        fig,ax=plt.subplots(1,2,figsize=(13,4.8),layout='constrained')
        for regime,color,label in [('Microporous_POP','#246e8e','微孔情景' if zh else 'Microporous scenario'),('Hierarchical_POP','#e1852f','分级孔情景' if zh else 'Hierarchical scenario')]:
            sub=transport[transport.regime==regime];r=np.linspace(10,200,400);phi=r*1e-6*np.sqrt(sub.kv_s.iloc[0]/sub.D_eff.iloc[0]);eta=3*(phi/np.tanh(phi)-1)/phi**2
            ax[0].plot(r,eta,label=label,color=color);ax[0].scatter(sub.radius_um,sub.eta,color=color,s=28)
            p=pd.read_csv(DATA/f'profile_{regime}_200.csv');ax[1].plot(p.r_over_R,p.c_analytic,label=label,color=color);ax[1].scatter(p.r_over_R[::20],p.c_fvm[::20],s=18,facecolors='none',edgecolors=color)
        ax[0].axhline(.4,ls='--',color='gray',lw=1);ax[0].axvline(129.3793,ls=':',color='gray',lw=1)
        ax[0].set(xlabel='颗粒半径 (µm)' if zh else 'Particle radius (µm)',ylabel='内部有效性因子 η' if zh else 'Internal effectiveness factor η',ylim=(0,1.06));ax[0].legend()
        ax[1].set(xlabel='r/R',ylabel='C(r)/C(surface)',title='R = 200 µm；实线解析，圆圈有限体积' if zh else 'R = 200 µm; line: analytical, circles: FVM');ax[1].legend()
        fig.suptitle('恒定参数球形颗粒：尚未经实验标定' if zh else 'Spherical particles with fixed parameters: no experimental calibration');save(fig,'pore_transport')
        labels_en=['Baseline reaction and metrology','Current/voltage/mass balances','Matched prospective optimization','Defined substrate validation','Particle-size and pore study','Independent batches and transfer','Grant dossier and milestones','Reproduction and data handover']
        labels_zh=['基准反应与测量校准','电流／电压／物料平衡','配对前瞻优化','明确反应的底物验证','粒径与孔结构研究','独立批次与条件迁移','大创材料与阶段验收','复现及数据交接']
        starts=[1,2,4,5,6,8,2,10];ends=[3,5,8,9,10,11,11,12]
        fig,ax=plt.subplots(figsize=(13,4.7),layout='constrained');ax.barh(range(8),np.array(ends)-starts+1,left=np.array(starts)-.5,color=['#276d89','#ad80ad','#e18c30','#318875']*2)
        ax.set_yticks(range(8),labels_zh if zh else labels_en);ax.invert_yaxis();ax.set_xticks(range(1,13));ax.set_xlim(.5,12.5)
        ax.set(xlabel='相对项目启动月份；按证据门槛推进' if zh else 'Month from project start; progression requires evidence gates',title='一年执行计划：不保证立项、获奖或发表' if zh else 'One-year execution plan: no guarantee of funding, awards or publication');ax.grid(axis='x',alpha=.2);ax.set_axisbelow(True);save(fig,'one_year_roadmap')
    print('Six bilingual plots generated; structure grid is saved by the audit.')
if __name__=='__main__':main()
