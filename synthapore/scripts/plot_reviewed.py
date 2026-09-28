"""Scientific figures from saved numerical benchmarks, not chemical observations."""
import csv
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np
from PIL import Image

BASE=Path(__file__).resolve().parents[1]
OUT=BASE/'reports/figures'
BLUE,GOLD,OLIVE,GRAY='#0B5D87','#CF7900','#757B43','#696969'
FILES=[]

def rows(name):
    with (BASE/name).open(encoding='utf-8-sig',newline='') as stream:
        return list(csv.DictReader(stream))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def settings(cn):
    plt.rcParams.update({'font.family':'Microsoft YaHei' if cn else 'Arial','font.size':9,
       'axes.titlesize':10,'axes.labelsize':9,'legend.fontsize':7.3,'axes.spines.top':False,
       'axes.spines.right':False,'svg.fonttype':'none','svg.hashsalt':'synthapore-reviewed',
       'axes.unicode_minus':False,'figure.facecolor':'white','savefig.facecolor':'white'})

def canvas(title):
    fig,axes=plt.subplots(1,2,figsize=(8,4.6),layout='constrained')
    fig.suptitle(title,fontsize=11.5)
    return fig,axes

def save(fig,stem,lang,sources):
    for ext in ['png','svg']:
        path=OUT/f'{stem}_{lang}.{ext}'
        fig.savefig(path,dpi=300,metadata={'Date':None} if ext=='svg' else {})
        record={'file':path.name,'format':ext,'sha256':sha(path),
          'sources':{name:sha(BASE/name) for name in sources},
          'evidence':'Synthetic shapes or explicitly assumed analytical models; no molecular or material validation'}
        if ext=='png':
            with Image.open(path) as im:
                record.update(pixels=list(im.size),dpi=list(im.info['dpi']))
        else:
            tree=ET.parse(path);ns={'s':'http://www.w3.org/2000/svg'}
            record.update(editable_text=bool(tree.findall('.//s:text',ns)),embedded_raster=bool(tree.findall('.//s:image',ns)))
        FILES.append(record)
    plt.close(fig)

def denoising(lang):
    cn=lang=='chinese'
    fig,(a,b)=canvas('合成几何去噪：不是分子扩散生成' if cn else 'SYNTHETIC GEOMETRY DENOISING: not molecular generation')
    data=rows('results/equivariant/denoising_metrics.csv')
    models=list(dict.fromkeys(r['model'] for r in data))
    models=[m for m in models if m.startswith('EGNN') or m.startswith('egnn')]+[m for m in models if not (m.startswith('EGNN') or m.startswith('egnn'))]
    labels=[]
    for m in models:
        if '444' in m: labels.append('EGNN '+m[-4:])
        elif m=='identity': labels.append('恒等' if cn else 'Identity')
        elif 'PCA' in m: labels.append('PCA 平面' if cn else 'PCA plane')
        else: labels.append('训练拟合平滑' if cn else 'Fitted smoother')
    y=np.arange(len(models))
    for offset,split,color,label in [(-.18,'test',BLUE,'测试集' if cn else 'Test'),(.18,'ood_warped',GOLD,'扭曲 OOD' if cn else 'Warped OOD')]:
        values=[float(next(r for r in data if r['model']==m and r['split']==split)['coordinate_RMSE']) for m in models]
        a.barh(y+offset,values,height=.32,color=color,label=label)
    a.set_yticks(y,labels);a.invert_yaxis()
    a.set(xlim=(0,.24),ylim=(6.2,-.5),xlabel='坐标 RMSE / 合成长度单位' if cn else 'Coordinate RMSE / synthetic length unit',
          title='a  独立测试与分布变化' if cn else 'a  Held-out and shifted-shape errors')
    a.legend(frameon=False,loc='lower center',ncol=2)
    losses=rows('results/equivariant/training_losses.csv')
    for seed,color in zip(['4441','4442','4443'],[BLUE,GOLD,OLIVE]):
        selected=[r for r in losses if r['seed']==seed]
        b.plot([int(r['epoch']) for r in selected],[float(r['validation_coordinate_MSE']) for r in selected],color=color,label=seed)
    b.set(xlabel='训练轮次' if cn else 'Epoch',ylabel='验证坐标 MSE' if cn else 'Validation coordinate MSE',
          title='b  验证损失；不按测试集选轮次' if cn else 'b  Validation loss selects checkpoints',yscale='log')
    b.legend(frameon=False)
    save(fig,'Fig1_Equivariance_Denoising',lang,['results/equivariant/denoising_metrics.csv','results/equivariant/training_losses.csv'])

def md(lang):
    cn=lang=='chinese'
    fig,(a,b)=canvas('解析谐振子：积分误差与恒温涨落' if cn else 'ANALYTIC OSCILLATOR: integration error and thermal fluctuations')
    data=rows('results/dynamics/nve_timestep_summary.csv')
    for seed,color in zip(['101','102','103','104'],[BLUE,GOLD,OLIVE,GRAY]):
        selected=[r for r in data if r['seed']==seed]
        a.loglog([float(r['timestep_fs']) for r in selected],[float(r['max_relative_total_energy_change']) for r in selected],
                 'o-',color=color,markersize=3,label=seed)
    a.set(xlabel='时间步 / fs' if cn else 'Timestep / fs',ylabel='最大相对总能量变化' if cn else 'Maximum relative total-energy change',
          title='a  固定时长 2000 fs；四个种子' if cn else 'a  Fixed 2000 fs duration; four seeds')
    a.set_xticks([.25,.5,1,2,4],['0.25','0.5','1','2','4'])
    a.xaxis.set_minor_formatter(NullFormatter())
    a.legend(frameon=False,ncol=2)
    stats=rows('results/dynamics/thermostat_summary.csv')
    values=[float(r['variance_ratio_to_canonical']) for r in stats]
    b.scatter(range(3),values,s=45,color=[BLUE,GOLD,OLIVE],zorder=3)
    b.axhline(1,color=GRAY,ls='--',lw=1)
    for i,r in enumerate(stats):
        b.annotate(f"{float(r['replicate_mean_temperature_K']):.2f} K",(i,values[i]),xytext=(0,10),textcoords='offset points',ha='center',fontsize=8)
    b.set_xticks(range(3),['BAOAB','Berendsen\n21 DOF','Berendsen\n24 DOF'])
    b.set(xlim=(-.5,2.5),ylim=(1e-8,10),yscale='log',
          ylabel='温度方差 / 正则参考方差' if cn else 'Temperature variance / canonical reference',
          title='b  八个副本；标注为温度均值' if cn else 'b  Eight replicas; labels show mean T')
    save(fig,'Fig2_MD_Validation',lang,['results/dynamics/nve_timestep_summary.csv','results/dynamics/thermostat_summary.csv'])

def neb(lang):
    cn=lang=='chinese'
    fig,(a,b)=canvas('解析双势阱 CI-NEB：不对应目标化学反应' if cn else 'ANALYTIC DOUBLE-WELL CI-NEB: not a chemical reaction')
    data=rows('results/dynamics/neb_final_bands.csv')
    summary=rows('results/dynamics/neb_summary.csv')
    chosen=[r for r in summary if r['seed']=='1' and float(r['tolerance_eV_A'])==1e-5]
    x=np.linspace(-1.2,1.2,160);y=np.linspace(-.15,1.05,120)
    xx,yy=np.meshgrid(x,y);zz=(xx**2-1)**2+4*(yy-.65*(1-xx**2))**2
    a.contour(xx,yy,zz,levels=[.1,.3,.6,1.,1.5,2.5,4.],colors='#C0C0C0',linewidths=.7)
    for case,color in zip(chosen,[BLUE,GOLD,OLIVE]):
        band=sorted([r for r in data if r['case']==case['case']],key=lambda r:int(r['image_index']))
        positions=np.array([[float(r['x_A']),float(r['y_A'])] for r in band])
        arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(positions,axis=0),axis=1))];arc/=arc[-1]
        label=case['n_images']+(' 个图像' if cn else ' images')
        a.plot(*positions.T,'o-',color=color,ms=3,lw=1,label=label)
        b.plot(arc,[float(r['energy_eV']) for r in band],'o-',color=color,ms=3,lw=1,label=label)
    a.scatter([0],[.65],marker='*',s=85,color='#C0392B',zorder=5)
    a.set(xlabel='x / Å',ylabel='y / Å',title='a  收敛路径与已知鞍点' if cn else 'a  Converged paths and known saddle',ylim=(-.15,1.05))
    b.axhline(1,color=GRAY,ls='--',lw=1)
    b.set(xlabel='归一化路径弧长' if cn else 'Normalized path arc length',ylabel='解析势能 / eV' if cn else 'Analytic potential energy / eV',
          title='b  力阈值 $10^{-5}$ eV/Å；种子 1' if cn else 'b  Force tolerance $10^{-5}$ eV/Å; seed 1',ylim=(-.03,1.18))
    b.legend(frameon=False,loc='upper right')
    save(fig,'Fig3_NEB_Validation',lang,['results/dynamics/neb_final_bands.csv','results/dynamics/neb_summary.csv'])

def pore(lang):
    cn=lang=='chinese'
    fig,(a,b)=canvas('几何与合成吸附基准：没有原子级 POP 结构' if cn else 'GEOMETRY AND SYNTHETIC ADSORPTION: no atomistic POP')
    data=rows('results/pore/geometry_grid_convergence.csv')
    geometries=list(dict.fromkeys(r['geometry'] for r in data))
    for geom,label,color in zip(geometries,['圆孔' if cn else 'Circle','六边形孔' if cn else 'Hexagon','周期障碍物' if cn else 'Periodic obstacles'],[BLUE,GOLD,OLIVE]):
        selected=[r for r in data if r['geometry']==geom and r['resolution']=='640']
        a.plot([float(r['probe_A']) for r in selected],[100*float(r['grid_accessible_fraction']) for r in selected],'o-',ms=3,color=color,label=label)
    a.set(xlabel='探针半径 / Å' if cn else 'Probe radius / Å',ylabel='探针中心可达面积分数 / %' if cn else 'Probe-center accessible area / %',
          title='a  640 × 640 中点网格' if cn else 'a  640 × 640 midpoint grid',ylim=(-2,102))
    a.legend(frameon=False,loc='upper right')
    points=rows('results/pore/source_isotherm_points.csv')
    b.plot([float(r['relative_pressure']) for r in points],[float(r['source_rounded_volume_cm3_g']) for r in points],'o-',ms=3,color=GOLD,label='源码合成点' if cn else 'Synthetic source points')
    fit=next(r for r in rows('results/pore/BET_window_fits.csv') if r['data']=='source_rounded' and float(r['lower_requested'])==.01 and float(r['upper_requested'])==.3)
    vm,c=float(fit['fitted_Vm_cm3_g']),float(fit['fitted_C'])
    pressure=np.linspace(.01,.6,200)
    b.plot(pressure,vm*c*pressure/((1-pressure)*(1+(c-1)*pressure)),ls='--',color=BLUE,label='低压拟合及外推' if cn else 'Low-p fit and extrapolation')
    b.axvline(.35,color=GRAY,ls=':',lw=1)
    b.set(xlabel='相对压力 $p/p_0$' if cn else 'Relative pressure $p/p_0$',ylabel='合成吸附量 / cm$^3$ g$^{-1}$' if cn else 'Synthetic uptake / cm$^3$ g$^{-1}$',
          title='b  $p/p_0$ = 0.35 人工拼接' if cn else 'b  Artificial branch switch at $p/p_0$ = 0.35',ylim=(0,600))
    b.legend(frameon=False,loc='lower right')
    save(fig,'Fig4_Pore_Adsorption',lang,['results/pore/geometry_grid_convergence.csv','results/pore/source_isotherm_points.csv','results/pore/BET_window_fits.csv'])

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    for lang in ['english','chinese']:
        settings(lang=='chinese')
        for function in [denoising,md,neb,pore]:
            function(lang)
    (BASE/'results/figure_manifest.json').write_text(json.dumps(dict(figure_count=4,languages=['english','chinese'],formats=['png','svg'],figures=FILES),indent=2)+'\n',encoding='utf-8',newline='\n')
    print('Saved eight PNG and eight editable SVG figures.')

if __name__=='__main__':
    main()
