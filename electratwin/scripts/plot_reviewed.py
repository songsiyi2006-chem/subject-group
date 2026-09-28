"""Static publication graphics for unvalidated ElectraTwin numerical models."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]

def main():
    out=ROOT/'reports/figures';out.mkdir(parents=True,exist_ok=True)
    def j(p):return json.loads((ROOT/p).read_text(encoding='utf-8'))
    base=j('results/transport/baseline_solution.json');source=j('results/compatibility/captured_model_arrays.json')['pde_result'];audit=j('results/source_audit.json')
    grid=pd.read_csv(ROOT/'results/transport/grid_convergence.csv');pool=pd.read_csv(ROOT/'results/control/candidate_pool.csv');paired=pd.read_csv(ROOT/'results/control/paired_budget_comparison.csv');manifest=[]
    for zh in [False,True]:
        lang='chinese' if zh else 'english';pick=lambda en,cn:cn if zh else en
        plt.rcParams.update({'font.family':'Microsoft YaHei' if zh else 'Arial','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False,'axes.unicode_minus':False,'svg.fonttype':'none','svg.hashsalt':'electratwin-reviewed','figure.facecolor':'white','axes.linewidth':.7})
        def save(fig,name,sources):
            for ext in ['png','svg']:
                p=out/f'{name}_{lang}.{ext}';fig.savefig(p,dpi=300,**({'metadata':{'Date':None}} if ext=='svg' else {}))
                r=dict(file=p.name,format=ext,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),evidence='Uncalibrated numerical model; no measured or industrial validation',sources={s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources})
                if ext=='png':
                    with Image.open(p) as im:r.update(pixels=list(im.size),dpi=list(im.info['dpi']))
                else:
                    t=p.read_text(encoding='utf-8');r.update(editable_text='<text' in t,embedded_raster='<image' in t)
                manifest.append(r)
            plt.close(fig)
        fig,axes=plt.subplots(1,2,figsize=(7.2,3.8),layout='constrained');f=base['field'];c=np.array(f['A_mol_m3'])
        im=axes[0].pcolormesh(np.linspace(0,60,c.shape[0]+1),np.linspace(0,300,c.shape[1]+1),c.T,cmap='cividis',vmin=0,vmax=50,rasterized=False)
        im.set_edgecolor('face')
        cb=fig.colorbar(im,ax=axes[0],orientation='horizontal',pad=.20,label=pick('A concentration / mM','A 浓度 / mM'));cb.solids.set_rasterized(False);cb.solids.set_edgecolor('face')
        axes[0].set(xlabel=pick('Channel position / mm','流道位置 / mm'),ylabel=pick('Gap coordinate / micrometre','间隙坐标 / 微米'),title=pick('a  Conservative 2D model','a  守恒二维模型'))
        axes[1].semilogy(source['channel_coords_x_mm'],source['anode_current_profile_x'],color='#9ca3af',lw=1.2,label=pick('Source model (n = 1 flux)','原模型（通量 n = 1）'))
        axes[1].semilogy(f['x_mm'],f['wall_current_A_m2'],color='#075985',lw=1.5,label=pick('Reviewed model (n = 2)','修订模型（n = 2）'))
        axes[1].set(xlabel=pick('Channel position / mm','流道位置 / mm'),ylabel=pick('Wall current density / A/m2','壁面电流密度 / A/m2'),title=pick('b  Different model definitions','b  不同模型定义'));axes[1].legend(fontsize=6.5,loc='upper right')
        fig.suptitle(pick('UNVALIDATED MODEL: 450 uL/min, imposed potential index 0.48 V','未验证模型：450 微升/分，给定电位指标 0.48 V'),fontsize=10)
        save(fig,'Fig1_Transport_Field',['results/transport/baseline_solution.json','results/compatibility/captured_model_arrays.json'])

        fig,axes=plt.subplots(1,2,figsize=(7.2,3.8),layout='constrained');pairs={(50,12),(100,24),(200,48),(400,96)};selected=grid[[tuple(v) in pairs for v in grid[['nx','ny']].to_numpy()]].sort_values('nx')
        axes[0].plot(selected.nx,selected.conversion_pct,'o-',color='#075985');axes[0].set_xscale('log',base=2);axes[0].set_xticks([50,100,200,400],['50 x 12','100 x 24','200 x 48','400 x 96']);axes[0].tick_params(axis='x',rotation=20);axes[0].set(xlabel=pick('Joint mesh refinement (nx x ny)','联合网格加密（nx x ny）'),ylabel=pick('Model conversion / %','模型转化率 / %'),title=pick('a  Discretization convergence','a  离散化收敛'))
        rows=audit['campaign']['rows'];xx=np.arange(1,len(rows)+1);raw=[r['reconstructed_unclipped_FE_pct'] for r in rows];clipped=[r['faradaic_efficiency'] for r in rows]
        axes[1].plot(xx,raw,'o-',color='#d97706',ms=4,label=pick('Recomputed, no clipping','重算原值，未裁剪'));axes[1].plot(xx,clipped,'s--',color='#075985',ms=3,label=pick('Source reported','原文报告值'));axes[1].axhline(100,c='#555555',lw=.7,ls=':');axes[1].set(xlabel=pick('Source grid condition index','原网格条件序号'),ylabel='FE / %',title=pick('b  Source FE inconsistency','b  原 FE 不一致'));axes[1].legend(fontsize=6.5);axes[1].set_ylim(0,205)
        fig.suptitle(pick('NUMERICAL AUDIT: convergence does not validate chemistry','数值审计：离散收敛不等于化学验证'),fontsize=10)
        save(fig,'Fig2_Numerical_Audit',['results/transport/grid_convergence.csv','results/source_audit.json'])

        fig,axes=plt.subplots(1,2,figsize=(7.2,3.8),layout='constrained');front=pool.pareto_STY_negative_SEC
        axes[0].scatter(pool.STY_assumed_product_kg_m3_day/1000,pool.SEC_assumed_product_kwh_kg,s=13,c='#adb2bb',label=pick('Dominated','被支配点'))
        axes[0].scatter(pool.loc[front,'STY_assumed_product_kg_m3_day']/1000,pool.loc[front,'SEC_assumed_product_kwh_kg'],s=21,c='#d97706',marker='^',label=pick('Nondominated','非支配点'))
        axes[0].set(xlabel='STY / (1000 kg/m3/day)',ylabel='SEC / kWh/kg',title=pick('a  81-point model pool','a  81 点模型候选池'));axes[0].legend(fontsize=6.5)
        for _,r in paired.iterrows():
            win=r.EHVI_full_pool_HV_fraction>r.random_full_pool_HV_fraction
            axes[1].plot([0,1],[r.random_full_pool_HV_fraction,r.EHVI_full_pool_HV_fraction],'o-',color='#075985' if win else '#d97706',lw=.8,ms=3,alpha=.75)
        axes[1].set_xticks([0,1],['Random','GP MC-EHVI']);axes[1].set_xlim(-.2,1.2);axes[1].set_ylim(.88,1.005);axes[1].set_ylabel(pick('Fraction of full-pool hypervolume','占全池超体积的比例'));axes[1].set_title(pick('b  Equal budget: 15 evaluations','b  相同预算：15 次评价'))
        axes[1].text(.5,.89,pick('8 paired seeds; EHVI wins 5, loses 3','8 个配对种子；EHVI 5 胜 3 负'),ha='center',fontsize=7)
        fig.suptitle(pick('SIMULATION ONLY: STY versus SEC; FE fixed by model identity','仅模拟：STY 与 SEC 折衷；FE 由模型恒等式固定'),fontsize=10)
        save(fig,'Fig3_Sequential_Optimization',['results/control/candidate_pool.csv','results/control/paired_budget_comparison.csv'])
    (ROOT/'results/figure_manifest.json').write_text(json.dumps(dict(figures=manifest,figure_count=3,languages=['english','chinese'],formats=['png','svg']),indent=2)+'\n',encoding='utf-8',newline='\n');print('Created 6 PNG and 6 editable vector SVG figures.')

if __name__=='__main__':main()
