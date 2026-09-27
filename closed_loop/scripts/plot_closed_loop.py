"""Plot saved numerical results; no model fitting or synthetic relabeling."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]

def main():
    data=json.loads((ROOT/'results/audit/active_learning.json').read_text());pred=pd.read_csv(ROOT/'results/audit/candidate_predictions.csv');out=ROOT/'reports/figures';out.mkdir(parents=True,exist_ok=True)
    chinese=Path('C:/Windows/Fonts/msyh.ttc')
    for zh in [False,True]:
        lang='chinese' if zh else 'english'
        font=FontProperties(fname=str(chinese)).get_name() if zh and chinese.exists() else 'DejaVu Sans'
        plt.rcParams.update({'font.family':font,'font.size':10,'axes.unicode_minus':False,'svg.fonttype':'none'})
        fig,ax=plt.subplots(1,2,figsize=(11.5,4.5),layout='constrained')
        names=['original','scaled_mle','scaled_fixed','mean_baseline'];labels=['原始 GP','标准化 MLE GP','固定核 GP','均值基线'] if zh else ['Original GP','Scaled MLE GP','Fixed-kernel GP','Mean baseline']
        values=[data['LOO'][n]['MAE_pp'] for n in names]
        bars=ax[0].barh(labels,values,color=['#b45309','#2563eb','#0891b2','#64748b']);ax[0].invert_yaxis();ax[0].set_xlim(0,85)
        ax[0].bar_label(bars,fmt='%.2f',padding=4);ax[0].set_xlabel('留一 MAE / 产率百分点' if zh else 'Leave-one-out MAE / yield percentage points');ax[0].set_title('六条模拟反馈；非实验验证' if zh else 'Six mock observations; no experimental validation')
        for name,label,color in [('scaled_mle',labels[1],'#2563eb'),('scaled_fixed',labels[2],'#0891b2')]:
            d=pred[(pred.method==name)&np.isclose(pred.electrolyte_concentration_M,.1)&np.isclose(pred.temperature_C,30)]
            ax[1].plot(d.current_density_mA_cm2,d.EI_pp,'o-',ms=3,label=label,color=color)
        ax[1].axvline(12,color='#64748b',ls=':',label='已采样，排除' if zh else 'Previously sampled; excluded');ax[1].set_xlabel('电流密度 / mA cm$^{-2}$' if zh else 'Current density / mA cm$^{-2}$');ax[1].set_ylabel('EI / 产率百分点' if zh else 'EI / yield percentage points');ax[1].set_title('0.10 M，30 °C；潜在函数不确定度' if zh else '0.10 M, 30 °C; latent uncertainty');ax[1].legend(fontsize=8)
        fig.savefig(out/f'feedback_audit_{lang}.png',dpi=180);plt.close(fig)
        fig,ax=plt.subplots(figsize=(8.6,4.6),layout='constrained');sweep=pd.read_csv(ROOT/'results/audit/stoichiometry_sweep.csv')
        for area in [.75,1.5,3.]:
            d=sweep[np.isclose(sweep.scale_mmol,.2)&np.isclose(sweep.volume_mL,6)&np.isclose(sweep.area_cm2,area)]
            ax.plot(d.current_density_mA_cm2,d.duration_min,'o-',label=f'A = {area:.2f} cm²')
        ax.set(xlabel='电流密度 / mA cm$^{-2}$' if zh else 'Current density / mA cm$^{-2}$',ylabel='理想恒流时间 / min' if zh else 'Ideal constant-current duration / min',title='0.20 mmol，2.20 F/mol；面积约定改变时间' if zh else '0.20 mmol, 2.20 F/mol; area convention changes duration');ax.legend();ax.grid(alpha=.2)
        fig.savefig(out/f'charge_sensitivity_{lang}.png',dpi=180);plt.close(fig)
        fig,ax=plt.subplots(figsize=(10.6,4.7),layout='constrained')
        tasks=['导师审阅与反应定义','仪器培训与分析校准','中性态/阳离子 DFT 试算','小规模基线与重复测量','数据录入与模型检验','下一轮条件与独立复测','本科项目材料与复现归档'] if zh else ['Mentor review and reaction definition','Instrument training and assay calibration','Neutral/cation DFT pilot','Small-scale baseline and repeats','Data ingestion and model checking','Next-round conditions and independent repeats','Undergraduate proposal and reproducible archive']
        starts=[1,2,2,4,5,6,7];lengths=[1,2,3,2,1,2,2]
        for i,(start,length) in enumerate(zip(starts,lengths)):ax.barh(i,length,left=start-.5,height=.58,color='#2563eb' if i%2==0 else '#0891b2')
        ax.set_yticks(range(len(tasks)),tasks);ax.invert_yaxis();ax.set_xticks(range(1,9));ax.set_xlim(.5,8.5);ax.grid(axis='x',alpha=.2)
        ax.set_xlabel('拟议周次；以导师、仪器与数据条件为前提' if zh else 'Proposed week; conditional on supervision, access and evidence');ax.set_title('八周落地计划（尚未启动湿实验或 DFT）' if zh else 'Eight-week integration plan (wet work and DFT not started)')
        fig.savefig(out/f'integration_gantt_{lang}.png',dpi=180);plt.close(fig)
    print('Created six source-backed PNG figures.')
if __name__=='__main__':main()
