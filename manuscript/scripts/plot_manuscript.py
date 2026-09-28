"""Publication figures from archived data and declared post hoc arithmetic."""
from pathlib import Path
import csv, hashlib, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

ROOT=Path(__file__).resolve().parents[2]
HERE=ROOT/'manuscript'; OUT=HERE/'figures'
def rows(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
BLUE,ORANGE,PURPLE='#176b8e','#cf7a11','#785584'

def main():
    OUT.mkdir(parents=True,exist_ok=True);records=[]
    gain=rows(HERE/'results/paired_gain_transfer.csv');mse=rows(HERE/'results/mse_attribution.csv')
    tail=rows(HERE/'results/residual_observable_comparison.csv');rich=rows(HERE/'results/tail_richardson.csv')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','axes.unicode_minus':True})
    for lang in ('english','chinese'):
        cn=lang=='chinese';plt.rcParams['font.family']=['Microsoft YaHei','DejaVu Sans'] if cn else ['Arial','DejaVu Sans']
        def save(fig,stem,sources):
            for ext in ('png','svg'):
                p=OUT/f'{stem}_{lang}.{ext}';fig.savefig(p,dpi=320,facecolor='white')
                records.append({'file':p.relative_to(HERE).as_posix(),'sha256':sha(p),'sources':{s:sha(ROOT/s) for s in sources}})
            plt.close(fig)
        fig,ax=plt.subplots(figsize=(9,4.6));ax.set_xlim(0,10);ax.set_ylim(0,5);ax.axis('off')
        boxes=[(0.2,3.15,2.65,1.3,'实际量子参考\nRHF / UHF / FCI\n同一几何与基组' if cn else 'Executed quantum references\nRHF / UHF / FCI\nMatched geometries and bases'),
               (3.65,3.15,2.65,1.3,'监督学习\n仅 RHF 能量与梯度\n三个配对种子' if cn else 'Supervised learning\nRHF energy and gradient labels\nThree paired seeds'),
               (7.1,3.15,2.65,1.3,'有符号误差归因\n学习误差与参考偏差\n保留交叉项' if cn else 'Error attribution\nLearning and reference bias\nRetain the cross term'),
               (0.2,0.75,2.65,1.3,'FCI 极小值与曲率\nMorse 参数化\n电子模型误差仍存在' if cn else 'FCI minimum and curvature\nMorse parameterization\nElectronic model error remains'),
               (3.65,0.75,2.65,1.3,'核振动离散化\nH2 / D2，J = 0\n网格与边界分别变化' if cn else 'Nuclear discretization\nH2 / D2, J = 0\nVary grid and boundary separately'),
               (7.1,0.75,2.65,1.3,'按观测量验证\n浅束缚态与配分函数\n不以矩阵残差代替' if cn else 'Observable-specific validation\nShallow state and partition sum\nBeyond the matrix residual')]
        for i,(x,y,w,h,t) in enumerate(boxes):
            ax.add_patch(Rectangle((x,y),w,h,facecolor='#f1f6f8' if i<3 else '#faf4e9',edgecolor='#b0bec5',lw=1));ax.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=10,linespacing=1.6)
        for start,end in [((2.9,3.8),(3.55,3.8)),((6.35,3.8),(7,3.8)),((2.9,1.4),(3.55,1.4)),((6.35,1.4),(7,1.4)),((1.52,3.05),(1.52,2.17))]:
            ax.add_patch(FancyArrowPatch(start,end,arrowstyle='-|>',mutation_scale=13,color='#536575'))
        ax.text(5,4.8,'两条受控计算链，共同检验误差是否传递到目标观测量' if cn else 'Two controlled calculation chains test errors at the target observable',ha='center',fontsize=11)
        ax.text(5,.22,'支持案例：路径与热化学、分子数据、传输与动力学、分析定量；证据类型分别保留' if cn else 'Supporting cases: paths, thermochemistry, molecular data, transport, kinetics and analytical quantification',ha='center',fontsize=8.8)
        fig.subplots_adjust(left=.01,right=.99,bottom=.03,top=.98)
        save(fig,'Fig1_Evidence_Design',['manuscript/results/analysis_summary.json'])

        fig,axs=plt.subplots(1,2,figsize=(9,4.6),layout='constrained')
        for ax,split,letter in zip(axs,['test','ood'],'ab'):
            for i,r in enumerate(v for v in gain if v['split']==split):
                ax.plot([0,1],[float(r['learning_RMSE_reduction_percent']),float(r['total_FCI_RMSE_reduction_percent'])],marker=['o','s','^'][i],color=[BLUE,ORANGE,PURPLE][i],label=r['seed'],lw=1.7)
            ax.set_xticks([0,1],['对 RHF 标签' if cn else 'Against RHF labels','对 FCI 参考' if cn else 'Against FCI']);ax.set_xlim(-.25,1.25);ax.set_ylim(-3,105)
            ax.set_ylabel('能量 RMSE 降低 / %' if cn else 'Reduction in energy RMSE / %');ax.set_title(letter+'  '+(('测试集，8 个匹配几何' if split=='test' else '域外集，5 个匹配几何') if cn else ('Test, 8 matched geometries' if split=='test' else 'OOD, 5 matched geometries')))
            ax.legend(title='种子' if cn else 'Seed',frameon=False);ax.axhline(0,color='#999999',lw=.7)
        save(fig,'Fig3_Gain_Transfer',['manuscript/results/paired_gain_transfer.csv'])

        fig,axs=plt.subplots(1,2,figsize=(9,4.6),layout='constrained')
        for ax,split,letter in zip(axs,['test','ood'],'ab'):
            labels=[]
            for i,r in enumerate(v for v in mse if v['split']==split and v['objective']=='energy_gradient'):
                vals=[float(r[k])*1e6 for k in ['learning_MSE_Hartree2','reference_bias_MSE_Hartree2','cross_term_Hartree2','total_MSE_Hartree2']]
                for j,val in enumerate(vals):ax.scatter(j+(i-1)*.08,val,color=[BLUE,ORANGE,PURPLE][i],marker=['o','s','^'][i],s=35,label=r['seed'] if j==0 else None,zorder=3)
            ax.set_yscale('symlog',linthresh=.03);ax.set_xticks(range(4),['学习项','参考项','交叉项','总误差'] if cn else ['Learning','Reference','Cross term','Total']);ax.set_ylabel('均方误差分量 / mHartree²' if cn else 'Mean-square component / mHartree²')
            ax.set_title(letter+'  '+(('测试集' if split=='test' else '匹配域外集') if cn else ('Test' if split=='test' else 'Matched OOD')));ax.axhline(0,color='#999999',lw=.7);ax.legend(frameon=False)
        save(fig,'Fig4_MSE_Attribution',['manuscript/results/mse_attribution.csv'])

        fig,axs=plt.subplots(1,2,figsize=(9,4.6),layout='constrained')
        t=[r for r in tail if int(r['states'])==17];x=np.arange(len(t));labels=[r['right_A'].split('.')[0]+' / '+str(round(float(r['intervals'])/1000,1)) for r in t]
        axs[0].semilogy(x,[float(r['eigen_residual_Hartree']) for r in t],marker='o',color=BLUE,label='矩阵残差' if cn else 'Matrix residual')
        axs[0].semilogy(x,[float(r['same_state_binding_absolute_error_Hartree']) for r in t],marker='s',color=ORANGE,label='最浅态束缚能误差' if cn else 'Shallow binding error')
        axs[0].set_xticks(x,labels,rotation=32);axs[0].set_xlabel('右边界 Å / 区间数 ×10³' if cn else 'Right boundary Å / intervals ×10³');axs[0].set_ylabel('Hartree');axs[0].set_title('a  '+('不同层次的数值误差' if cn else 'Errors at different numerical levels'));axs[0].legend(frameon=False)
        for i,r in enumerate(rich):
            label=f"{r['right_A'].split('.')[0]} Å, {int(r['fine_intervals'])}"
            axs[1].plot([0,1],[float(r['fine_signed_relative_error_percent']),float(r['extrapolated_signed_relative_error_percent'])],marker=['o','s','^'][i],color=[BLUE,ORANGE,PURPLE][i],label=label)
        axs[1].axhline(0,lw=.8,color='#999999');axs[1].set_xticks([0,1],['细网格','Richardson 估计'] if cn else ['Fine grid','Richardson estimate']);axs[1].set_xlim(-.2,1.2);axs[1].set_ylabel('束缚能有符号偏差 / %' if cn else 'Signed binding-energy error / %');axs[1].set_title('b  '+('固定边界的事后外推' if cn else 'Post hoc extrapolation at fixed boundary'));axs[1].legend(frameon=False,fontsize=8)
        save(fig,'Fig6_Observable_Convergence',['manuscript/results/residual_observable_comparison.csv','manuscript/results/tail_richardson.csv'])
    record={'generator_sha256':sha(Path(__file__)),'figure_groups':4,'files':records,'scope':'Four new bilingual figure groups; no regenerated quantum or neural data.'}
    (HERE/'results/figure_manifest.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'figures':len(records)}))

if __name__=='__main__':main()
