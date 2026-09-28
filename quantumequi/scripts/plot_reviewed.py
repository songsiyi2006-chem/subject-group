"""Source-backed figures separating quantum references from numerical controls."""
import csv
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

BASE = Path(__file__).resolve().parents[1]
OUT = BASE/'reports/figures'
COLORS = ['#0B5D87','#CF7900','#757B43','#7C4D79']
FILES = []


def rows(name):
    with (BASE/name).open(encoding='utf-8-sig',newline='') as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def settings(cn):
    plt.rcParams.update({'font.family':'Microsoft YaHei' if cn else 'Arial','font.size':9,
        'axes.titlesize':10,'axes.labelsize':9,'legend.fontsize':7,'axes.spines.top':False,
        'axes.spines.right':False,'svg.fonttype':'none','svg.hashsalt':'quantumequi-reviewed',
        'axes.unicode_minus':False,'figure.facecolor':'white','savefig.facecolor':'white'})


def canvas(title):
    fig, axes = plt.subplots(1,2,figsize=(8,4.6),layout='constrained')
    fig.suptitle(title,fontsize=11)
    return fig,axes


def save(fig,stem,lang,sources):
    for ext in ('png','svg'):
        path = OUT/f'{stem}_{lang}.{ext}'
        fig.savefig(path,dpi=300,metadata={'Date':None} if ext=='svg' else {})
        record={'file':path.name,'format':ext,'sha256':sha(path),
            'sources':{name:sha(BASE/name) for name in sources},
            'evidence':'Actual H2 RHF references, untrained source audits and designed analytic controls explicitly distinguished'}
        if ext=='png':
            with Image.open(path) as im:
                record.update(pixels=list(im.size),dpi=list(im.info['dpi']))
        else:
            tree=ET.parse(path);ns={'s':'http://www.w3.org/2000/svg'}
            record.update(editable_text=bool(tree.findall('.//s:text',ns)),embedded_raster=bool(tree.findall('.//s:image',ns)))
        FILES.append(record)
    plt.close(fig)


def quantum(lang):
    cn=lang=='chinese'
    fig,(a,b)=canvas('源矩阵缺陷与实际 H2 量子参考' if cn else 'SOURCE MATRIX DEFECT AND EXECUTED H2 QUANTUM REFERENCE')
    clipping='results/electronic/source_clipping_sensitivity.csv'
    data=sorted([r for r in rows(clipping) if r['case']=='reactant'],key=lambda r:float(r['clipping_floor']))
    a.loglog([float(r['clipping_floor']) for r in data],[float(r['occupied_eigenvalue_sum_times_two_eV']) for r in data],
             'o-',color=COLORS[0],markersize=4)
    a.axvline(1e-5,color='#666666',ls='--',lw=1,label='源码阈值' if cn else 'Source threshold')
    a.set(title='a  源 EHT 对数值截断的依赖' if cn else 'a  Source EHT clipping sensitivity',
          xlabel='重叠本征值下限' if cn else 'Overlap eigenvalue floor',
          ylabel='占据本征值和 ×2 / 名义 eV' if cn else 'Occupied eigenvalue sum x2 / nominal eV')
    a.legend(loc='upper right')
    reference='results/electronic/h2_reference.csv'
    data=sorted(rows(reference),key=lambda r:float(r['R_A']))
    b.plot([float(r['R_A']) for r in data],[float(r['energy_total_Hartree']) for r in data],color='#AAAAAA',lw=1)
    for split,label,color,marker in zip(['train','validation','test','ood_stretch'],
        ['训练','验证','测试','拉伸 OOD'] if cn else ['Train','Validation','Test','Stretched OOD'],COLORS,['o','s','^','D']):
        points=[r for r in data if r['split']==split]
        b.scatter([float(r['R_A']) for r in points],[float(r['energy_total_Hartree']) for r in points],s=14,label=label,color=color,marker=marker,zorder=3)
    b.set(title='b  已执行 RHF/STO-3G 扫描' if cn else 'b  Executed RHF/STO-3G scan',
          xlabel='H–H 距离 / Å' if cn else 'H-H distance / Angstrom',
          ylabel='含核排斥的总能 / Hartree' if cn else 'Total energy including repulsion / Hartree')
    b.legend(loc='lower right')
    save(fig,'Fig1_Quantum_Reference',lang,[clipping,reference])


def force(lang):
    cn=lang=='chinese'
    fig,axes=canvas('未训练 EGNN：微分一致性检验' if cn else 'UNTRAINED EGNN: DERIVATIVE CONSISTENCY')
    source='results/potential/force_finite_differences.csv'
    data=rows(source)
    for axis,dtype,panel in zip(axes,['float32','float64'],['a','b']):
        for geometry,label,color in zip(['reactant','product','source_candidate'],
            ['反应物坐标','产物坐标','源码候选点'] if cn else ['Reactant coordinates','Product coordinates','Source candidate'],COLORS):
            selected=sorted([r for r in data if r['dtype']==dtype and r['geometry']==geometry],key=lambda r:float(r['displacement_A']))
            axis.loglog([float(r['displacement_A']) for r in selected],[float(r['force_max_abs_error']) for r in selected],
                        'o-',color=color,label=label,markersize=4)
        axis.set(title=panel+'  '+dtype,xlabel='差分步长 / Å' if cn else 'Difference step / Angstrom',
                 ylabel='最大力差 / 名义单位' if cn else 'Maximum force discrepancy / nominal units')
        axis.legend(loc='best')
    save(fig,'Fig2_Force_Audit',lang,[source])


def neb(lang):
    cn=lang=='chinese'
    fig,axes=canvas('解析势 CI-NEB：数值收敛，不是化学势垒' if cn else 'ANALYTIC CI-NEB: NUMERICAL, NOT CHEMICAL BARRIERS')
    table='results/path/analytic_final_images.csv'
    refs='results/path/analytic_references.json'
    data=rows(table)
    references=json.loads((BASE/refs).read_text())
    for axis,surface,panel in zip(axes,['tilted_sine_2d','periodic_curve_3d'],['a','b']):
        for n,color in zip([7,11,17],COLORS):
            case=f'{surface}_n{n}_tol1e-05_seed11'
            selected=sorted([r for r in data if r['case']==case],key=lambda r:int(r['image']))
            columns=['x','y'] if surface=='tilted_sine_2d' else ['x','y','z']
            xyz=np.array([[float(r[k]) for k in columns] for r in selected])
            arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(xyz,axis=0),axis=1))];arc/=arc[-1]
            energy=np.array([float(r['energy']) for r in selected]);energy-=energy[0]
            axis.plot(arc,energy,'o-',color=color,markersize=3,label=f'{n} '+('图像' if cn else 'images'))
        exact=references[surface]['forward_barrier']
        axis.axhline(exact,color='#777777',ls='--',lw=1)
        axis.set(title=panel+'  '+(('非对称二维势' if cn else 'Asymmetric 2D potential') if panel=='a' else ('周期三维势' if cn else 'Periodic 3D potential')),
            xlabel='归一化路径弧长' if cn else 'Normalized path arc length',
            ylabel='相对势能 / 无量纲' if cn else 'Relative potential energy / dimensionless')
        axis.set_ylim(top=exact*1.30)
        axis.legend(loc='upper right')
    save(fig,'Fig3_NEB_Validation',lang,[table,refs])


def thermo(lang):
    cn=lang=='chinese'
    fig,(a,b)=canvas('频率与热化学审计：源候选点非驻点' if cn else 'FREQUENCY AND THERMOCHEMISTRY AUDIT: SOURCE IS NONSTATIONARY')
    spectrum='results/thermochemistry/source_spectra.csv'
    data=sorted([r for r in rows(spectrum) if r['dtype']=='torch.float64_autograd' and r['space']=='projected'],key=lambda r:int(r['mode']))
    frequencies=np.array([float(r['signed_frequency_nominal_cm']) for r in data])
    a.bar(np.arange(1,len(data)+1),frequencies,color=[COLORS[1] if v<0 else COLORS[0] for v in frequencies])
    a.axhline(0,color='#777777',lw=.8)
    a.set(title='a  刚体投影后的源候选频谱' if cn else 'a  Projected source-candidate spectrum',
          xlabel='内部模式序号' if cn else 'Internal mode index',
          ylabel='有符号频率 / 名义 '+r'cm$^{-1}$' if cn else r'Signed frequency / nominal cm$^{-1}$',xticks=[1,6,12,18])
    table='results/thermochemistry/ideal_gas_temperature_pressure.csv'
    data=rows(table)
    for case,label,color in zip(['linear_minimum','linear_saddle','nonlinear_minimum','nonlinear_saddle'],
        ['线性极小值','线性鞍点','非线性极小值','非线性鞍点'] if cn else ['Linear minimum','Linear saddle','Nonlinear minimum','Nonlinear saddle'],COLORS):
        selected=sorted([r for r in data if r['case']==case and float(r['pressure_Pa'])==100000],key=lambda r:float(r['temperature_K']))
        b.plot([float(r['temperature_K']) for r in selected],[float(r['G_correction_kcal_mol']) for r in selected],
               'o--' if 'saddle' in case else 'o-',label=label,color=color,markersize=3)
    b.set(title='b  人工控制体系：理想气体，1 bar' if cn else 'b  Designed controls: ideal gas, 1 bar',
          xlabel='温度 / K' if cn else 'Temperature / K',
          ylabel='G − E 校正 / '+r'kcal mol$^{-1}$' if cn else r'G - E correction / kcal mol$^{-1}$')
    b.legend(loc='best')
    save(fig,'Fig4_Thermochemistry',lang,[spectrum,table])


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    for language in ['english','chinese']:
        settings(language=='chinese')
        for plot in [quantum,force,neb,thermo]:
            plot(language)
    (BASE/'results/figure_manifest.json').write_text(json.dumps({'figure_count':4,'languages':['english','chinese'],
        'formats':['png','svg'],'figures':FILES},indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n')
    print('Generated eight 300-DPI PNGs and eight editable SVGs from saved evidence.')


if __name__=='__main__':
    main()
