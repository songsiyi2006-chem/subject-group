"""Source-backed static PNG/SVG figures, explicitly labeling synthetic evidence."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.colors import Normalize
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]

def main():
    out=ROOT/'reports/figures';out.mkdir(parents=True,exist_ok=True);data=ROOT/'results/audit';summary=json.loads((data/'audit_summary.json').read_text());p=pd.read_csv(data/'pareto_conditions.csv');scope=pd.read_csv(data/'scope_mock.csv');meta=[]
    for zh in [False,True]:
        lang='chinese' if zh else 'english';pick=lambda en,cn:cn if zh else en
        font=FontProperties(fname='C:/Windows/Fonts/msyh.ttc').get_name() if zh and Path('C:/Windows/Fonts/msyh.ttc').exists() else 'Arial'
        plt.rcParams.update({'font.family':font,'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False,'axes.unicode_minus':False,'svg.fonttype':'none','svg.hashsalt':'ai4s-toolkit','axes.linewidth':.7,'lines.linewidth':1.2,'figure.facecolor':'white'})
        def save(fig,name,sourcefiles):
            for ext in ['png','svg']:
                path=out/f'{name}_{lang}.{ext}';kwargs={'dpi':300} if ext=='png' else {'metadata':{'Date':None}}
                fig.savefig(path,**kwargs)
                record=dict(file=path.name,evidence='synthetic/mock/illustrative as labeled; no experimental or DFT validation',sources={s:hashlib.sha256((data/s).read_bytes()).hexdigest() for s in sourcefiles},sha256=hashlib.sha256(path.read_bytes()).hexdigest(),format=ext)
                if ext=='png':
                    with Image.open(path) as im:record.update(pixels=list(im.size),dpi=list(im.info.get('dpi',[])))
                else:record.update(editable_text='<text' in path.read_text(encoding='utf-8'),embedded_raster='<image' in path.read_text(encoding='utf-8'))
                meta.append(record)
            plt.close(fig)
        fig=plt.figure(figsize=(7.2,4.0));ax=fig.add_subplot(121,projection='3d');bx=fig.add_subplot(122);front=p.nondominated
        ax.scatter(p.loc[~front,'yield_pct'],p.loc[~front,'FE_pct'],p.loc[~front,'SEC_kWh_kg'],c='#a8adb4',s=11,alpha=.55,label=pick('Dominated','被支配点'))
        ax.scatter(p.loc[front,'yield_pct'],p.loc[front,'FE_pct'],p.loc[front,'SEC_kWh_kg'],c='#d97706',marker='^',s=23,edgecolors='#222222',linewidths=.4,label=pick('Nondominated (23)','非支配点（23）'))
        ax.set_xlabel(pick('Yield / %','产率 / %'),labelpad=0);ax.set_ylabel('FE / %',labelpad=0);ax.set_zlabel('SEC / kWh/kg',labelpad=0);ax.tick_params(pad=0);ax.view_init(25,130);ax.legend(loc='upper left',fontsize=6.5,bbox_to_anchor=(0,1.02));ax.set_title(pick('a  Exact nondominance','a  严格非支配筛选'),pad=20)
        Y,F=np.meshgrid(np.linspace(40,95,60),np.linspace(40,95,60));Z=summary['figures']['Pareto']['energy_numerator_kWh_kg']/(F/100)
        cs=bx.contour(Y,F,Z,levels=[.85,1,1.2,1.5,1.8],colors='#4b5563',linewidths=.75);bx.clabel(cs,fmt='%.2f',fontsize=7,manual=[(55,90.5),(87,76.9),(47,64.1),(47,51.3),(47,42.7)])
        bx.scatter(p.yield_pct,p.FE_pct,s=10,facecolors='none',edgecolors='#9ca3af',linewidths=.5);bx.scatter(p.loc[front,'yield_pct'],p.loc[front,'FE_pct'],s=18,c='#d97706',marker='^',edgecolors='#222222',linewidths=.3)
        bx.set(xlabel=pick('Yield / %','产率 / %'),ylabel='FE / %',title=pick('b  Assumed SEC contours / kWh/kg','b  假设 SEC 等值线 / kWh/kg'),xlim=(40,95),ylim=(40,95));fig.subplots_adjust(left=.025,right=.985,bottom=.18,top=.79,wspace=.34)
        fig.suptitle(pick('SYNTHETIC electrosynthesis landscape — 80 generated conditions','模拟电合成情景——80 个生成条件'),fontsize=10,y=.98)
        fig.text(.5,.025,pick('Contours: 2 electrons, 3.2 V, assumed product MW 0.223 kg/mol; no SEC noise.','等值线假设：2 电子、3.2 V、产物分子量 0.223 kg/mol；不含 SEC 噪声。'),ha='center',fontsize=7)
        save(fig,'Fig1_Pareto_Electrosynthesis',['pareto_conditions.csv'])

        fig,ax=plt.subplots(1,2,figsize=(7.2,3.9),gridspec_kw={'width_ratios':[1.35,1]},layout='constrained');arr=scope[['given_prediction_pct','given_mock_label_pct']].to_numpy();cmap=plt.get_cmap('cividis')
        im=ax[0].pcolormesh(np.arange(3)-.5,np.arange(7)-.5,arr,vmin=0,vmax=100,cmap=cmap,rasterized=False);ax[0].set_ylim(5.5,-.5);ax[0].set_yticks(range(6),scope.substrate);ax[0].set_xticks([0,1],[pick('Given prediction','给定预测'),pick('Mock label','模拟标签')]);ax[0].set_title(pick('a  Illustrative values / %','a  示例数值 / %'))
        for (i,k),v in np.ndenumerate(arr):ax[0].text(k,i,f'{v:.1f}',ha='center',va='center',color='black' if v>55 else 'white',fontsize=8)
        cb=fig.colorbar(im,ax=ax[0],fraction=.05,pad=.03,label=pick('Yield label / %','产率标签 / %'));cb.solids.set_rasterized(False);cb.solids.set_edgecolor('face')
        corr=np.corrcoef(arr.T);im2=ax[1].pcolormesh(np.arange(3)-.5,np.arange(3)-.5,corr,vmin=0,vmax=1,cmap='Blues',rasterized=False);ax[1].set_ylim(1.5,-.5);ax[1].set_aspect('equal');ax[1].set_xticks([0,1],['P','M']);ax[1].set_yticks([0,1],['P','M']);ax[1].set_title(pick('b  Pearson matrix; n = 6','b  Pearson 相关矩阵；n = 6'))
        for (i,k),v in np.ndenumerate(corr):ax[1].text(k,i,f'{v:.3f}',ha='center',va='center',color='white',fontsize=10)
        ax[1].set_xlabel(pick('P = prediction; M = mock label\nMAE = 3.417 pp; no measured validation','P = 给定预测；M = 模拟标签\nMAE = 3.417 个百分点；未实测验证'),fontsize=7)
        fig.suptitle(pick('MOCK substrate scope — neither series is experimental evidence','模拟底物范围——两列均不是实验证据'),fontsize=10)
        save(fig,'Fig2_Substrate_Scope_Heatmap',['scope_mock.csv'])

        fig,ax=plt.subplots(figsize=(7.2,3.8));values=summary['figures']['energy']['source_levels_kcal_mol'];labels=pick(['Reactants*','Cation*','Proposed TS*','Intermediate*','Product*'],['反应物*','阳离子*','假设过渡态*','中间体*','产物*'])
        for i,value in enumerate(values):
            ax.plot([i-.24,i+.24],[value,value],color='#075985',lw=2);ax.text(i,value+1.8,f'{value:+.1f}',ha='center',fontsize=8)
            if i<4:ax.plot([i+.24,i+.76],[value,values[i+1]],'--',color='#9ca3af',lw=.8)
        ax.annotate('',xy=(1.6,22.3),xytext=(1.6,14.8),arrowprops=dict(arrowstyle='<->',lw=.7,color='#222222'));ax.text(1.66,17.5,pick('7.5 level gap','7.5 能级差'),fontsize=7)
        ax.axhline(0,c='#6b7280',lw=.7,ls=':');ax.set_xticks(range(5),labels);ax.set_ylabel(pick('Illustrative relative level / kcal/mol','示意相对能级 / kcal/mol'));ax.set_ylim(-25,32);ax.set_xlim(-.5,4.5);ax.set_title(pick('HYPOTHETICAL energy profile — no DFT or TS validation','假设能级图——没有 DFT 或过渡态验证'),fontsize=10)
        fig.subplots_adjust(left=.12,right=.98,bottom=.22,top=.88);fig.text(.5,.04,pick('* Source values only. Electron/proton references, species balance and stationary points are undefined.','* 仅使用原文示例值；电子/质子参照、物种配平和驻点性质均未确定。'),ha='center',fontsize=7)
        save(fig,'Fig3_Reaction_Energy_Profile',['hypothetical_energy.csv'])

        fig,ax=plt.subplots(2,2,figsize=(7.2,5.4),layout='constrained');cal=pd.read_csv(data/'synthetic_calibration.csv');trace=pd.read_csv(data/'synthetic_chromatogram.csv');rec=pd.read_csv(data/'deconvolution_recovery.csv');rel=pd.read_csv(data/'relaxation_scenarios.csv');coef=summary['synthetic_recovery']
        ax[0,0].scatter(cal.concentration_mM,cal.area,s=10,c='#075985');xx=np.linspace(0,2,50);ax[0,0].plot(xx,coef['calibration_slope']*xx+coef['calibration_intercept'],c='#333333',lw=.8);ax[0,0].set(xlabel=pick('Concentration / mM','浓度 / mM'),ylabel=pick('Area / arbitrary units','面积 / 任意单位'),title=pick('a  Synthetic calibration, n = 21','a  合成校准数据，n = 21'))
        ax[0,1].plot(trace.time_min,trace.signal_au,c='#9ca3af',lw=.6,label=pick('Simulated signal','模拟信号'));ax[0,1].plot(trace.time_min,trace.fitted_au,c='#075985',lw=1,label=pick('Known-shape fit','已知峰形拟合'));ax[0,1].set(xlabel=pick('Time / min','时间 / min'),ylabel=pick('Signal / arbitrary units','信号 / 任意单位'),title=pick('b  Two synthetic peaks','b  两个合成峰'));ax[0,1].legend(fontsize=6.5)
        for mismatch,color,marker,label in [(False,'#075985','o',pick('Correct widths','峰宽正确')),(True,'#d97706','^',pick('Wrong target width','目标峰宽错误'))]:
            subset=rec[rec.width_mismatch==mismatch];by=subset.groupby('separation_min').relative_error_pct.apply(lambda v:np.mean(np.abs(v)))
            ax[1,0].plot(by.index,by.values,marker=marker,color=color,label=label);ax[1,0].set(xlabel=pick('Peak separation / min','峰间距 / min'),ylabel=pick('Mean absolute area error / %','面积平均绝对误差 / %'),title=pick('c  Assumption sensitivity; 10 seeds','c  假设敏感性；10 个种子'));ax[1,0].legend(fontsize=6.5)
        ax[1,1].plot(rel.repetition_time_s,rel.relative_bias_pct,'o-',color='#075985');ax[1,1].axhline(0,color='#6b7280',lw=.7,ls=':');ax[1,1].set(xlabel=pick('Repetition time / s','重复间隔 / s'),ylabel=pick('Relative integral-ratio bias / %','积分比相对偏差 / %'),title=pick('d  Assumed T1: product 3 s, IS 1 s','d  假定 T1：产物 3 s，内标 1 s'))
        fig.suptitle(pick('SYNTHETIC analytical diagnostics — not instrument validation','合成分析诊断——不是仪器或测量验证'),fontsize=10)
        save(fig,'Fig4_Analytical_Quality_Control',['synthetic_calibration.csv','synthetic_chromatogram.csv','deconvolution_recovery.csv','relaxation_scenarios.csv'])
    (ROOT/'results/figure_manifest.json').write_text(json.dumps(dict(figures=meta,figure_count=4,languages=['english','chinese'],formats=['png','svg'],scope='PNG previews at 300 DPI and editable vector SVG; not a journal acceptance or measured-data claim'),indent=2)+'\n',encoding='utf-8',newline='\n')
    print('Created 8 PNG and 8 editable SVG figures.')
if __name__=='__main__':main()
