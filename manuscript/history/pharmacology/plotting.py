"""Figures from frozen public records; no new scientific simulation."""
from pathlib import Path
import csv, json, hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

HERE=Path(__file__).resolve().parent
SRC=HERE/'sources'/'ai4pharm-lead-developability-suite'
OUT=HERE/'figures'
OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':['DejaVu Sans','Microsoft YaHei'],'font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black','figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white','svg.fonttype':'none','svg.hashsalt':'pharmacology-history-v1','axes.spines.top':False,'axes.spines.right':False})
BLUE='#28628f';ORANGE='#ba6624';GREEN='#38735a';GRAY='#b6bdc3'
def data(path):return list(csv.DictReader(path.read_text(encoding='utf-8-sig').splitlines()))
def f(x,k):return float(x[k])
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
records=[]
FINISH_ONLY=None
def finish(fig,number,lang,sources,note):
    if FINISH_ONLY is not None and (number,lang)!=FINISH_ONLY:
        plt.close(fig)
        return
    fig.subplots_adjust(left=.085,right=.98,top=.85,bottom=.24,wspace=.42)
    fig.text(.5,.04,note,ha='center',fontsize=8,color='black')
    for ext in ['png','svg']:
        p=OUT/f'Figure_P{number}_{lang}.{ext}';fig.savefig(p,dpi=300)
        rec=dict(file=p.relative_to(HERE).as_posix(),sha256=sha(p),format=ext,sources={str(x.relative_to(HERE)).replace('\\','/'):sha(x) for x in sources})
        if ext=='png':
            with Image.open(p) as im:rec.update(pixels=list(im.size),dpi=list(im.info.get('dpi',[])))
            assert rec['pixels']==[2700,1260]
        else:
            text=p.read_text(encoding='utf-8');rec.update(editable_text='<text' in text,embedded_raster='<image' in text)
            assert rec['editable_text'] and not rec['embedded_raster']
        records.append(rec)
    plt.close(fig)

def main(only=None):
    global FINISH_ONLY
    FINISH_ONLY=only
    records.clear()
    p1=HERE/'sources/derived/mpo_plot.csv';m=data(p1)
    p2=HERE/'sources/derived/tpd_plot.csv';t=data(p2)
    pl=SRC/'projects/task02_tpd/data_task2/linker_conformers.csv';link=data(pl)
    p3=HERE/'sources/derived/covalent_plot.csv';c=data(p3)
    p4=HERE/'sources/derived/optimization_plot.csv';g=data(p4)
    pm=SRC/'projects/task19_metadynamics/outputs/atomistic_cv_trajectory.csv';md=data(pm)
    pd=SRC/'projects/task20_denovo/outputs/docking/candidate_vina_scores.csv';d=data(pd)
    for lang in ['english','chinese']:
        zh=lang=='chinese';labels=['口服参照','毒性比较物','bRo5'] if zh else ['Oral references','Toxicity comparators','bRo5']
        fig,ax=plt.subplots(1,2,figsize=(9,4.2));groups=['Oral Drugs','Toxic Dropouts','bRo5 Modalities'];tot=[];conv=[]
        for j,group in enumerate(groups):
            xs=[x for x in m if x['archetype']==group];xpos=j+np.linspace(-.22,.22,len(xs));scores=np.array([f(x,'cns_mpo_approx') for x in xs]);low=[f(x,'mpo_pka_sensitivity_min') for x in xs];high=[f(x,'mpo_pka_sensitivity_max') for x in xs]
            ax[0].vlines(xpos,low,high,color=GRAY,lw=1.3);ax[0].scatter(xpos,scores,s=19,color=BLUE,zorder=3);ax[0].hlines(np.median(scores),j-.3,j+.3,color='black',lw=1.3)
            tot.append(sum(int(x['conformers_embedded']) for x in xs));conv.append(sum(int(x['conformers_converged']) for x in xs))
        ax[0].set(xticks=range(3),xticklabels=labels,ylim=(-.05,6.3),ylabel='近似 CNS-MPO 分数' if zh else 'Approximate CNS-MPO score',title='a  评分与 pKa 假设范围' if zh else 'a  Scores and pKa assumption ranges')
        xx=np.arange(3);ax[1].bar(xx,conv,color=BLUE,label='收敛' if zh else 'Converged');ax[1].bar(xx,np.array(tot)-conv,bottom=conv,color=ORANGE,label='未收敛' if zh else 'Not converged')
        for j in range(3):ax[1].text(j,tot[j]+1,f'{conv[j]}/{tot[j]}',ha='center',fontsize=9)
        ax[1].set(xticks=xx,xticklabels=labels,ylim=(0,70),ylabel='抽样构象数' if zh else 'Sampled conformers',title='b  构象收敛' if zh else 'b  Conformer convergence');ax[1].legend(frameon=False,fontsize=8,loc='upper left')
        finish(fig,1,lang,[p1],'pKa 范围为假设情景；不是置信区间。评分未作药理校准。' if zh else 'pKa ranges are assumption scenarios, not confidence intervals. Scores are not pharmacologically calibrated.')

        fig,ax=plt.subplots(1,2,figsize=(9,4.2))
        for alpha,color in zip([.01,1,10,100],[GRAY,BLUE,ORANGE,GREEN]):
            xs=[x for x in t if f(x,'alpha')==alpha];ax[0].plot([f(x,'P_total_nM') for x in xs],[f(x,'EPT_nM') for x in xs],label=f'α = {alpha:g}',color=color,lw=1.6)
        ax[0].axvline(131.6227766016838,color='black',ls=':',lw=1);ax[0].set(xscale='log',xlabel='降解剂总浓度 / nM' if zh else 'Total degrader / nM',ylabel='三元复合体 / nM' if zh else 'Ternary complex / nM',title='a  质量作用钩状曲线' if zh else 'a  Mass-action hook curves');ax[0].legend(frameon=False,fontsize=7.5,loc='upper left');ax[0].set_ylim(0,100)
        flex=[f(x,'distance_A') for x in link if x['linker']=='Flexible PEG'];rigid=[f(x,'distance_A') for x in link if x['linker']=='Rigid alkynyl']
        ax[1].hist(flex,bins=np.arange(6,13.1,.35),color=BLUE,alpha=.8,label='柔性 PEG（100 次）' if zh else 'Flexible PEG (100 draws)');ax[1].axvline(np.mean(rigid),color=ORANGE,lw=2,label='刚性炔基均值（100 次）' if zh else 'Rigid mean (100 draws)');ax[1].set(xlabel='端点距离 / Å' if zh else 'Endpoint distance / Å',ylabel='柔性片段抽样频数' if zh else 'Flexible-fragment draw count',title='b  连接子代理几何' if zh else 'b  Linker proxy geometry',xlim=(6,13));ax[1].legend(frameon=False,fontsize=7.5)
        finish(fig,2,lang,[p2,pl],'合成平衡情景与未加权片段抽样；不能推断降解效力或结合熵。' if zh else 'Synthetic equilibrium scenarios and unweighted fragment sampling; neither degradation potency nor binding entropy.')

        fig,ax=plt.subplots(1,2,figsize=(9,4.2));xx=np.arange(len(c));ids=[x['id'] for x in c]
        ax[0].scatter(xx-.08,[f(x,'rapid_equilibrium_efficiency') for x in c],color=ORANGE,marker='o',label='快速平衡' if zh else 'Rapid equilibrium');ax[0].scatter(xx+.08,[f(x,'qssa_efficiency') for x in c],color=BLUE,marker='s',label='准稳态' if zh else 'Quasi-steady state')
        ax[0].set(yscale='log',xticks=xx,xticklabels=ids,ylabel='效率 / M⁻¹ s⁻¹' if zh else 'Efficiency / M⁻¹ s⁻¹',title='a  动力学近似的差异' if zh else 'a  Kinetic approximation changes');ax[0].legend(frameon=False,fontsize=8,loc='lower right' if zh else 'lower left')
        vals=[f(x,'center_fplus_proxy') for x in c];ax[1].vlines(xx,1e-9,vals,color=GRAY);ax[1].scatter(xx,vals,color=BLUE,s=27);ax[1].set(yscale='log',ylim=(1e-9,1),xticks=xx,xticklabels=ids,ylabel='反应中心 LUMO 布居代理' if zh else 'Reaction-centre LUMO population proxy',title='b  EHT 定域的限制' if zh else 'b  Limits of EHT localization')
        finish(fig,3,lang,[p3],'实际 EHT 描述符 + 假设速率常数；没有实验动力学或安全性校准。' if zh else 'Executed EHT descriptors plus assumed rate constants; no experimental kinetic or safety calibration.')

        fig,ax=plt.subplots(1,3,figsize=(9,4.2));ax[0].plot([f(x,'time_ps') for x in md],[f(x,'distance_nm') for x in md],color=BLUE,marker='.',ms=4);ax[0].ticklabel_format(axis='y',useOffset=False);ax[0].set(xlabel='总时间 / ps' if zh else 'Total time / ps',ylabel='距离变量 / nm' if zh else 'Distance CV / nm',title='a  仅 1 ps 偏置动力学' if zh else 'a  Only 1 ps biased dynamics',ylim=(3.94,3.98))
        best=[f(g[0],'best_geometry'),f(g[-1],'best_geometry'),f(g[-1],'random_best_geometry')];median=[f(g[0],'median_geometry'),f(g[-1],'median_geometry'),f(g[-1],'random_selected_median')];xx=np.arange(3)
        ax[1].bar(xx-.17,best,.34,color=BLUE,label='最佳' if zh else 'Best');ax[1].bar(xx+.17,median,.34,color=GRAY,label='中位数' if zh else 'Median');ax[1].set(xticks=xx,xticklabels=['初始','进化','随机'] if zh else ['Initial','Evolved','Random'],ylim=(0,25),ylabel='任意量纲几何分数' if zh else 'Arbitrary geometry score',title='b  等预算参照' if zh else 'b  Equal-budget reference');ax[1].legend(frameon=False,fontsize=7,ncol=2,loc='upper left');ax[1].text(.5,.62,'20.670 vs 20.652',transform=ax[1].transAxes,ha='center',fontsize=8)
        ax[2].scatter([f(x,'geometry_score') for x in d],[f(x,'vina_score_kcal_mol') for x in d],color=ORANGE,s=22);ax[2].set(xlabel='几何分数' if zh else 'Geometry score',ylabel='Vina / kcal mol⁻¹',title='c  12 个对接候选' if zh else 'c  12 docked candidates');ax[2].tick_params(labelsize=7.5)
        finish(fig,4,lang,[pm,p4,pd],'不同观测量采用不同轴范围。PMF 未收敛；几何/对接分数不是亲和力。' if zh else 'Different observables use different axis ranges. PMF is unconverged; geometry/docking scores are not affinities.')
    if only is not None:
        old=json.loads((OUT/'manifest.json').read_text(encoding='utf-8'))['figures']
        replacement={r['file']:r for r in records}
        records[:]=[replacement.get(r['file'],r) for r in old]
    (OUT/'manifest.json').write_text(json.dumps({'figures':records,'source':'frozen repository records and arithmetic-only derived CSVs','dimensions_inches':[9,4.2],'dpi':300},ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'files':len(records),'png':sum(x['format']=='png' for x in records),'svg':sum(x['format']=='svg' for x in records)}))

if __name__=='__main__':
    import sys
    main((3,'english') if '--figure3-english-only' in sys.argv else None)
