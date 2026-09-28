"""Bilingual, data-backed graphics for correlation, learning and nuclear models."""
from collections import defaultdict
import csv,hashlib,json
from pathlib import Path
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

BASE=Path(__file__).resolve().parents[1]
DATA=BASE/'results/extensions'
OUT=BASE/'reports/figures/extensions'
COLORS=['#0B5D87','#CF7900','#7C4D79','#757B43']
ARTIFACTS=[]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(p):
    with (DATA/p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def read(p):return json.loads((DATA/p).read_text(encoding='utf-8'))
def num(r,k):return float(r[k])
def settings(cn):
    plt.rcParams.update({'font.family':['Microsoft YaHei','DejaVu Sans'] if cn else ['Arial','DejaVu Sans'],'font.size':9,
        'axes.titlesize':10,'axes.labelsize':9,'legend.fontsize':7,'axes.spines.top':False,
        'axes.spines.right':False,'svg.fonttype':'none','svg.hashsalt':'quantumequi-extensions',
        'axes.unicode_minus':False,'figure.facecolor':'white','savefig.facecolor':'white'})
def canvas(title):
    f,a=plt.subplots(1,2,figsize=(9,4.8),layout='constrained');f.suptitle(title,fontsize=12);return f,a
def save(fig,name,lang,sources):
    for ext in ('png','svg'):
        path=OUT/f'{name}_{lang}.{ext}'
        fig.savefig(path,dpi=300,metadata={'Date':None} if ext=='svg' else {})
        r={'file':path.name,'format':ext,'sha256':sha(path),'sources':{s:sha(DATA/s) for s in sources}}
        if ext=='png':
            with Image.open(path) as im:r.update(pixels=list(im.size),dpi=list(im.info.get('dpi',[])))
        else:
            tree=ET.parse(path);ns={'s':'http://www.w3.org/2000/svg'}
            r.update(editable_text=bool(tree.findall('.//s:text',ns)),embedded_raster=bool(tree.findall('.//s:image',ns)))
        ARTIFACTS.append(r)
    plt.close(fig)


def dissociation(lang):
    cn=lang=='chinese';source='correlation/curve.csv';data=rows(source)
    f,axes=canvas('实际 H2 计算：方法与基组对解离曲线的影响' if cn else 'EXECUTED H2 CURVES: METHOD AND BASIS DEPENDENCE')
    for ax,basis,letter in zip(axes,['sto-3g','cc-pvdz'],'ab'):
        for method,color in zip(['RHF','UHF','FCI'],COLORS):
            p=sorted([r for r in data if r['basis'].lower()==basis and r['method']==method],key=lambda r:num(r,'R_A'))
            ax.plot([num(r,'R_A') for r in p],[num(r,'energy_Hartree') for r in p],'.-',label=method,color=color,ms=4)
        ax.set(title=f'{letter}  {basis.upper() if basis=="sto-3g" else "cc-pVDZ"}',xlabel='核间距 / Å' if cn else 'Internuclear distance / Angstrom',ylabel='总能量 / Hartree' if cn else 'Total energy / Hartree');ax.legend()
    save(f,'Ext1_Dissociation',lang,[source])


def correlation(lang):
    cn=lang=='chinese';source='correlation/curve.csv';data=rows(source);f,(a,b)=canvas('电子相关与 UHF 自旋污染' if cn else 'ELECTRONIC CORRELATION AND UHF SPIN CONTAMINATION')
    for i,basis in enumerate(['sto-3g','cc-pvdz']):
        basisname='STO-3G' if i==0 else 'cc-pVDZ'
        q=[r for r in data if r['basis'].lower()==basis];reference={round(num(r,'R_A'),8):num(r,'energy_Hartree') for r in q if r['method']=='FCI'}
        for j,method in enumerate(['RHF','UHF']):
            p=sorted([r for r in q if r['method']==method],key=lambda r:num(r,'R_A'))
            a.plot([num(r,'R_A') for r in p],[(num(r,'energy_Hartree')-reference[round(num(r,'R_A'),8)])*1000 for r in p],'.'+'-' if j==0 else '.--',color=COLORS[i],label=basisname+' '+method,ms=3)
        p=sorted([r for r in q if r['method']=='UHF'],key=lambda r:num(r,'R_A'))
        b.plot([num(r,'R_A') for r in p],[num(r,'S2') for r in p],'o-',color=COLORS[i],ms=3,label=basisname)
    a.set(title='a  同一基组内的近似误差' if cn else 'a  Within-basis approximation error',ylabel='E − E(FCI) / mHartree');b.set(title='b  UHF 自旋诊断' if cn else 'b  UHF spin diagnostic',ylabel='⟨S²⟩');b.axhline(0,color='#777',ls=':',lw=1)
    for ax in (a,b):ax.set_xlabel('核间距 / Å' if cn else 'Internuclear distance / Angstrom');ax.legend()
    save(f,'Ext2_Correlation_Spin',lang,[source])


def learning(lang):
    cn=lang=='chinese';curve='learning/learning_curves.csv';prediction='learning/ensemble_predictions.csv';c=rows(curve);p=rows(prediction)
    f,(a,b)=canvas('真实 RHF 标签上的监督学习：三个配对种子' if cn else 'SUPERVISED LEARNING ON RHF LABELS: THREE PAIRED SEEDS')
    for i,obj in enumerate(['energy_only','energy_gradient']):
        label=(['仅能量','能量与梯度'] if cn else ['Energy only','Energy + gradient'])[i]
        for j,seed in enumerate(sorted({r['seed'] for r in c})):
            q=[r for r in c if r['objective']==obj and r['seed']==seed]
            a.semilogy([num(r,'epoch') for r in q],[num(r,'selection_score') for r in q],color=COLORS[i],alpha=.55,lw=1,label=label if j==0 else None)
        q=sorted([r for r in p if r['objective']==obj],key=lambda r:num(r,'R_A'))
        x=[num(r,'R_A') for r in q];mean=np.array([num(r,'ensemble_gradient_mean_Hartree_A') for r in q]);std=np.array([num(r,'ensemble_gradient_std_Hartree_A') for r in q])
        b.plot(x,mean,color=COLORS[i],label=label);b.fill_between(x,mean-std,mean+std,color=COLORS[i],alpha=.14)
    reference=rows('learning/predictions.csv');q=sorted([r for r in reference if r['objective']=='energy_only' and r['seed']==sorted({r['seed'] for r in reference})[0]],key=lambda r:num(r,'R_A'))
    b.plot([num(r,'R_A') for r in q],[num(r,'reference_gradient_Hartree_A') for r in q],':',color='#222',label='RHF 参考' if cn else 'RHF reference')
    b.axvspan(1.9,2.7,color='#858b95',alpha=.1,label='拉伸 OOD' if cn else 'Stretched OOD')
    a.set(title='a  验证集联合误差' if cn else 'a  Joint validation score',xlabel='训练轮数' if cn else 'Epoch',ylabel='归一化 RMSE 之和' if cn else 'Sum of normalized RMSEs')
    b.set(title='b  三种子均值 ± 标准差' if cn else 'b  Three-seed mean ± sample SD',xlabel='核间距 / Å' if cn else 'Internuclear distance / Angstrom',ylabel='dE/dR / Hartree Å⁻¹')
    a.legend();b.legend();save(f,'Ext3_Learning',lang,[curve,prediction,'learning/predictions.csv'])


def benchmark(lang):
    cn=lang=='chinese';s='learning/metrics.csv';data=rows(s);baselinefile=BASE/'results/electronic/surrogate_metrics.csv'
    with baselinefile.open(encoding='utf-8-sig',newline='') as h:base=list(csv.DictReader(h))
    f,axes=canvas('泛化误差：所有种子与冻结基线' if cn else 'GENERALIZATION ERRORS: ALL SEEDS AND FROZEN BASELINES')
    for ax,col,title in zip(axes,['energy_RMSE_Hartree','gradient_RMSE_Hartree_A'],['a  能量误差' if cn else 'a  Energy error','b  梯度误差' if cn else 'b  Gradient error']):
        for i,split in enumerate(['test','ood_stretch']):
            for j,obj in enumerate(['energy_only','energy_gradient']):
                q=[num(r,col) for r in data if r['split']==split and r['objective']==obj];x=i+(-.18 if j==0 else .18)
                ax.scatter(np.full(len(q),x)+np.linspace(-.035,.035,len(q)),q,color=COLORS[j],s=28,label=(['神经：仅能量','神经：能量与梯度'] if cn else ['NN: energy only','NN: energy + gradient'])[j] if i==0 else None)
            for j,model in enumerate(['rbf_energy_gradient','linear']):
                # Exact names are checked against the frozen baseline ledger.
                names=sorted({r['model'] for r in base})
                chosen=next((n for n in names if (('rbf' in n.lower() and ('force' in n.lower() or 'gradient' in n.lower())) if j==0 else 'linear' in n.lower())),None)
                if chosen is None:raise ValueError('Missing required frozen baseline')
                q=[r for r in base if r['split']==split and r['model']==chosen];v=num(q[0],col)
                ax.plot([i-.38,i+.38],[v,v],color=COLORS[j+2],ls='--' if j==0 else ':',lw=1.6,label=(['RBF 联合拟合','线性基线'] if cn else ['Joint RBF','Linear baseline'])[j] if i==0 else None)
        ax.set_yscale('log');ax.set_xticks([0,1],['测试集' if cn else 'Test','拉伸 OOD' if cn else 'Stretched OOD']);ax.set_title(title);ax.set_ylabel('RMSE / '+('Hartree' if col.startswith('energy') else 'Hartree Å⁻¹'));ax.legend(loc='best')
    save(f,'Ext4_Generalization',lang,[s,'../electronic/surrogate_metrics.csv'])


def levels(lang):
    cn=lang=='chinese';models=read('vibration/model_parameters.json');model=next(m for m in models if m['basis'].lower()=='cc-pvdz');source='vibration/isotope_levels.csv';data=rows(source)
    f,axes=canvas('FCI 参数化 Morse 模型：非谐振动与同位素' if cn else 'FCI-PARAMETERIZED MORSE MODEL: ANHARMONIC ISOTOPES')
    x=np.linspace(.3,2.6,600);v=model['D_e_Hartree']*(1-np.exp(-model['a_A_inverse']*(x-model['R_e_A'])))**2
    for ax,iso,letter in zip(axes,['H2','D2'],'ab'):
        ax.plot(x,v,color='#333',label='Morse 势能' if cn else 'Morse potential')
        q=[r for r in data if r['model_id']==model['model_id'] and r['isotope']==iso and int(r['v'])<5]
        for i,r in enumerate(q):
            e=num(r,'E_numerical_Hartree');mask=v<e
            ax.hlines(e,x[mask].min(),x[mask].max(),color=COLORS[0],lw=1.8,label='数值束缚能级' if cn and i==0 else ('Numerical bound levels' if i==0 else None))
            harmonic=num(r,'E_harmonic_Hartree');ax.hlines(harmonic,1.72,2.15,color=COLORS[1],ls='--',label='谐振能级' if cn and i==0 else ('Harmonic levels' if i==0 else None))
            ax.text(2.2,harmonic,f'v={i}',fontsize=7,va='center')
        ax.set(xlim=(.3,2.55),ylim=(0,max(num(r,'E_harmonic_Hartree') for r in q)*1.16),title=f'{letter}  {iso} · cc-pVDZ',xlabel='核间距 / Å' if cn else 'Internuclear distance / Angstrom',ylabel='相对模型势能最低值 / Hartree' if cn else 'Energy above model minimum / Hartree');ax.legend(loc='upper left')
    save(f,'Ext5_Anharmonic_Levels',lang,['vibration/model_parameters.json',source])


def numerical_limits(lang):
    cn=lang=='chinese';grid='vibration/grid_and_domain_diagnostics.csv';thermal='vibration/partition_bound_only.csv';g=rows(grid);t=rows(thermal)
    f,(a,b)=canvas('核振动数值收敛与谐振近似偏差' if cn else 'NUCLEAR GRID CONVERGENCE AND HARMONIC APPROXIMATION')
    for i,iso in enumerate(['H2','D2']):
        q=sorted([r for r in g if r['model_id']=='analytic_control' and r['isotope']==iso and r['negative_boundary_control']=='False' and num(r,'left_A')==0 and num(r,'right_A')==8],key=lambda r:num(r,'intervals'))
        a.loglog([num(r,'intervals') for r in q],[num(r,'max_first3_error_cm_1') for r in q],'o-',color=COLORS[i],label=iso)
        model=next(r['model_id'] for r in t if 'cc-pvdz' in r['model_id'].lower())
        q=sorted([r for r in t if r['model_id']==model and r['isotope']==iso],key=lambda r:num(r,'T_K'))
        b.plot([num(r,'T_K') for r in q],[(num(r,'F_vib_bound_only_Hartree')-num(r,'F_vib_harmonic_Hartree'))*1000 for r in q],'o-',color=COLORS[i],label=iso)
    a.set(title='a  解析 Morse 对照，0–8 Å' if cn else 'a  Analytic Morse control, 0–8 Angstrom',xlabel='网格区间数' if cn else 'Grid intervals',ylabel='前三能级最大误差 / cm⁻¹' if cn else 'Maximum first-three-level error / cm⁻¹')
    b.set(title='b  cc-pVDZ 参数化模型，仅束缚态' if cn else 'b  cc-pVDZ model, bound states only',xlabel='温度 / K' if cn else 'Temperature / K',ylabel='F(vib, bound) − F(harmonic) / mHartree')
    a.legend();b.legend();save(f,'Ext6_Numerical_Limits',lang,[grid,thermal])


def error_budget(lang):
    cn=lang=='chinese';source='error_budget/matched_error_components.csv';data=rows(source);f,axes=canvas('学习误差与参考方法偏差：同一几何的分解' if cn else 'LEARNING ERROR VERSUS REFERENCE BIAS AT MATCHED GEOMETRIES')
    for ax,obj,letter in zip(axes,['energy_only','energy_gradient'],'ab'):
        q=[r for r in data if r['objective']==obj and r['split'] in ('test','ood_stretch')];groups=defaultdict(list)
        for r in q:groups[num(r,'R_A')].append(r)
        x=sorted(groups)
        for key,label,color,ls in [('learning_error_Hartree','学习误差' if cn else 'Learning error',COLORS[0],'-'),('RHF_minus_FCI_Hartree','RHF 参考偏差' if cn else 'RHF reference bias',COLORS[1],'--'),('total_error_vs_FCI_Hartree','相对 FCI 总误差' if cn else 'Total versus FCI',COLORS[2],'-')]:
            ax.plot(x,[np.mean([num(r,key) for r in groups[v]])*1000 for v in x],'.'+ls,color=color,label=label)
        ax.axhline(0,color='#777',lw=.8);ax.axvspan(1.9,2.7,color='#777',alpha=.08)
        ax.set(title=letter+'  '+(('仅拟合能量' if obj=='energy_only' else '能量与梯度') if cn else ('Energy-only network' if obj=='energy_only' else 'Joint network')),xlabel='核间距 / Å' if cn else 'Internuclear distance / Angstrom',ylabel='三种子平均有符号误差 / mHartree' if cn else 'Three-seed mean signed error / mHartree');ax.legend()
    save(f,'Ext7_Error_Budget',lang,[source])


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    for lang in ('english','chinese'):
        settings(lang=='chinese')
        for function in (dissociation,correlation,learning,benchmark,levels,numerical_limits,error_budget):function(lang)
    record={'figure_count':7,'languages':['english','chinese'],'formats':['png','svg'],
            'generator_sha256':sha(Path(__file__)),'source_hash_base':'quantumequi/results/extensions','figures':ARTIFACTS}
    (DATA/'figure_manifest.json').write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'files':len(ARTIFACTS),'panels_per_language':14}))


if __name__=='__main__':main()
