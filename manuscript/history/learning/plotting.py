"""Five paired publication figures from frozen historical CSVs, no scientific jobs."""
from pathlib import Path
import csv, hashlib, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'sources'; OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
font_manager.fontManager.addfont('C:/Windows/Fonts/msyh.ttc')
plt.rcParams.update({'font.family':['DejaVu Sans','Microsoft YaHei'],'font.size':10,
 'axes.titlesize':11,'axes.labelsize':10,'axes.edgecolor':'black','axes.labelcolor':'black',
 'text.color':'black','xtick.color':'black','ytick.color':'black','figure.facecolor':'white',
 'axes.facecolor':'white','savefig.facecolor':'white','svg.fonttype':'none','svg.hashsalt':'history-learning-v1',
 'axes.spines.top':False,'axes.spines.right':False,'axes.unicode_minus':False})
BLUE='#4575a5'; GOLD='#cb9d49'; PALE='#eaf0f6'
CMAP=LinearSegmentedColormap.from_list('light_blue',['white','#d9e5f1','#8fb1d0'])
manifest=[]
def rows(name):
 with (DATA/name).open(encoding='utf8') as h:return list(csv.DictReader(h))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def finish(fig,stem,lang,inputs,note):
 fig.text(.015,.025,note,fontsize=8.5,va='bottom')
 for ext in ['png','svg']:
  f=OUT/f'{stem}_{lang}.{ext}'
  meta={'Software':'Matplotlib'} if ext=='png' else {'Date':None,'Creator':'Matplotlib'}
  fig.savefig(f,dpi=300,metadata=meta)
  manifest.append({'path':str(f.relative_to(ROOT)).replace('\\','/'),'sha256':digest(f),
    'language':lang,'width_px':2700 if ext=='png' else None,'height_px':1260 if ext=='png' else None,
    'dpi':300 if ext=='png' else None,'width_inches':9,'height_inches':4.2,
    'inputs_sha256':{'sources/'+x:digest(DATA/x) for x in inputs},'caption_note':note})
 plt.close(fig)

for lang in ['english','chinese']:
 zh=lang=='chinese'
 # Figure L1: one bar per molecular object, range starts at zero.
 rs=rows('derived_conformers.csv');fig,ax=plt.subplots(figsize=(9,4.2))
 fig.subplots_adjust(left=.09,right=.98,bottom=.21,top=.85)
 bars=ax.bar(np.arange(len(rs)),[float(r['range_kcal_mol']) for r in rs],color=BLUE,width=.65,edgecolor='black',linewidth=.6)
 for i,(r,b) in enumerate(zip(rs,bars)):
  if r['force_field']=='UFF':b.set_facecolor('#dadada');b.set_hatch('///')
  ax.text(i,b.get_height()+.65,f"{float(r['range_kcal_mol']):.2f}",ha='center',fontsize=9)
 ax.set_xticks(np.arange(len(rs)),[r['id'] for r in rs]);ax.set_ylim(0,49)
 ax.set_ylabel('构象能量范围 / kcal mol⁻¹' if zh else 'Conformer energy range / kcal mol⁻¹')
 ax.set_xlabel('结构标识' if zh else 'Structure identifier')
 ax.set_title('历史构象筛选：每个有效结构有 50 条验收记录' if zh else 'Historical conformer screen: 50 accepted records per valid structure',pad=12)
 ax.text(.99,.95,'M07: UFF; 其余: MMFF94' if zh else 'M07: UFF; all others: MMFF94',transform=ax.transAxes,ha='right',va='top',fontsize=9)
 finish(fig,'Figure_L1_Conformer',lang,['derived_conformers.csv'],
 '500 条优化记录；不是 500 个独立极小点。M09 解析失败；M09R 为修复参考。' if zh else '500 optimization records, not 500 independent minima. M09 failed parsing; M09R is a repaired reference.')
 # Figure L2: SE has no invented confidence coverage or reference uncertainty.
 rs=rows('derived_vmc.csv');fig,ax=plt.subplots(figsize=(9,4.2))
 fig.subplots_adjust(left=.09,right=.98,bottom=.24,top=.84)
 x=np.arange(5);y=[float(r['signed_error_mEh']) for r in rs];se=[float(r['block_SE_mEh']) for r in rs]
 ax.errorbar(x,y,yerr=se,fmt='o',color=BLUE,ecolor='black',capsize=5,markersize=7,lw=1.2)
 ax.plot(x[-1],y[-1],'s',ms=8,mfc='white',mec='black')
 ax.axhline(1.6,color='black',ls='--',lw=1);ax.axhline(0,color='#999',lw=.6)
 ax.text(3.9,1.64,'1.6 mEh',ha='right',fontsize=9)
 ax.set_xticks(x,['H₂\n1.4011 bohr','H₂\n2.5 bohr','H₂\n4.0 bohr','H₂\n6.0 bohr','He'])
 ax.set_ylim(-.12,2.5);ax.set_xlim(-.45,4.45)
 ax.set_ylabel('VMC − 计算参考 / mEh' if zh else 'VMC − computed reference / mEh')
 ax.set_title('相同运行的内部检查通过，但 He 点精度目标未达' if zh else 'Internal checks pass in the same run; He misses the point-accuracy target',pad=12)
 finish(fig,'Figure_L2_VMC',lang,['derived_vmc.csv'],
 '误差棒：20 个保存块均值的 ±1 SE；参考为计算的 FCI CBS 外推；未计参考不确定性。' if zh else 'Bars: ±1 SE from 20 saved block means. Reference: computed FCI CBS extrapolation; reference uncertainty excluded.')
 # Figure L3: separate source-unit current scale, shared mV scales for voltage panels.
 rs=rows('phase28_window_robustness.csv');fig,axes=plt.subplots(1,3,figsize=(9,4.2))
 fig.subplots_adjust(left=.075,right=.985,bottom=.23,top=.75,wspace=.45)
 trace_ids=['nature2023_constant_potential','nature2023_constant_current','natsynth2026_constant_current']
 title_en=['2023 constant potential\nSource current ordinate','2023 constant current\nPotential change / mV','2026 constant current\nPotential change / mV']
 title_zh=['2023 恒电位\n源电流纵坐标变化','2023 恒电流\n电位变化 / mV','2026 恒电流\n电位变化 / mV']
 for k,(ax,tid) in enumerate(zip(axes,trace_ids)):
  starts=[0,1,5,10,20,30];wins=[1,3,6,12]
  grid=np.array([[float(next(r['change'] for r in rs if r['trace_id']==tid and float(r['startup_exclusion_h'])==s and float(r['window_h'])==w)) for w in wins] for s in starts])
  if k: grid*=1000
  im=ax.imshow(grid,cmap=CMAP,aspect='auto',vmin=-.16 if k==0 else -5,vmax=.05 if k==0 else 45)
  for i in range(6):
   for j in range(4):ax.text(j,i,f'{grid[i,j]:.3f}' if k==0 else f'{grid[i,j]:.1f}',ha='center',va='center',fontsize=8.2,color='black')
  ax.set_xticks(range(4),wins);ax.set_yticks(range(6),starts)
  ax.set_xlabel('窗口宽度 / h' if zh else 'Window width / h')
  if k==0:ax.set_ylabel('启动排除时长 / h' if zh else 'Startup exclusion / h')
  ax.set_title((title_zh if zh else title_en)[k],fontsize=10,pad=9)
  ax.axhline(2.5,color='black',lw=1.2)
 fig.suptitle('相同已发表轨迹：窗口定义改变表观漂移' if zh else 'The same published traces: window definitions change apparent drift',y=.97,fontsize=12)
 finish(fig,'Figure_L3_Trace_Windows',lang,['phase28_window_robustness.csv'],
 '每格为时间戳平衡的端点中位数差；不是独立实验或置信区间。电位两图共用色标范围。' if zh else 'Cells: timestamp-balanced endpoint median differences, not independent experiments or intervals. Voltage panels share a scale.')
 # Figure L4: raw offsets and constant-reference-subtracted contrasts.
 rs=rows('derived_reference_alignment.csv');fig,axes=plt.subplots(1,2,figsize=(9,4.2))
 fig.subplots_adjust(left=.075,right=.98,bottom=.23,top=.77,wspace=.3)
 for g,col,marker in [(1,BLUE,'o'),(2,GOLD,'s')]:
  rr=[r for r in rs if int(r['gfn'])==g]
  for ax,key in zip(axes,['raw_GFNmPBE0_eV','double_difference_eV']):
   ys=[float(r[key]) for r in rr];xx=np.arange(3)+(g-1.5)*.12
   ax.plot(xx,ys,marker=marker,linestyle='none',color=col,ms=7,label='GFN'+str(g))
 axes[0].set_ylim(-6.1,-4.5);axes[1].axhline(0,color='black',lw=.75);axes[1].set_ylim(-.06,.31)
 for ax in axes:
  ax.set_xticks(range(3),['Q01/P01','Q02/P04','Q03/P09']);ax.set_xlim(-.3,2.3)
  ax.set_ylabel('能量差 / eV' if zh else 'Energy contrast / eV');ax.legend(frameon=False,loc='upper left',ncol=2,fontsize=9)
 axes[0].set_title('原始 GFN − PBE0 偏移' if zh else 'Raw GFN − PBE0 offset')
 axes[1].set_title('以 Q01 为中心后的双重差值' if zh else 'Double difference after Q01 centering')
 fig.suptitle('同几何方法比较需要明确电荷态能量参考' if zh else 'Same-geometry comparison requires explicit charge-state energy references',y=.96,fontsize=12)
 finish(fig,'Figure_L4_Reference_Alignment',lang,['derived_reference_alignment.csv'],
 '固定核坐标、气相、PBE0/def2-SVP。PBE0 不是实验真值；原始偏移不是物理预测误差。' if zh else 'Fixed nuclei; gas phase; PBE0/def2-SVP. PBE0 is not experimental truth; raw offsets are not physical prediction errors.')
 # Figure L5: complete pairwise method/environment sensitivity, denominator explicit.
 rs=rows('derived_rank_reversals.csv');settings=['GFN1_gas','GFN1_acetonitrile','GFN1_dmf','GFN2_gas','GFN2_acetonitrile','GFN2_dmf']
 grid=np.zeros((6,6),int)
 for r in rs:
  i=settings.index(r['setting_a']);j=settings.index(r['setting_b']);grid[i,j]=grid[j,i]=int(r['rank_reversals'])
 fig,ax=plt.subplots(figsize=(9,4.2));fig.subplots_adjust(left=.19,right=.83,bottom=.26,top=.80)
 ax.imshow(grid,cmap=CMAP,aspect='auto',vmin=0,vmax=12)
 labels=['GFN1\ngas','GFN1\nMeCN','GFN1\nDMF','GFN2\ngas','GFN2\nMeCN','GFN2\nDMF']
 ax.set_xticks(range(6),labels);ax.set_yticks(range(6),[s.replace('\n',' / ') for s in labels])
 for i in range(6):
  for j in range(6):ax.text(j,i,str(grid[i,j]),ha='center',va='center',fontsize=10)
 ax.set_title('排序反转数：每种比较包含 12 个分子、66 个分子对' if zh else 'Rank reversals: 12 molecules and 66 molecular pairs per comparison',pad=14)
 finish(fig,'Figure_L5_Rank_Reversals',lang,['derived_rank_reversals.csv'],
 '15 组方法/环境比较；对角线为自身比较。MeCN = 乙腈；没有用此矩阵评估外部准确率。' if zh else '15 method/environment comparisons; diagonal is self-comparison. MeCN = acetonitrile. No external accuracy is evaluated.')

(OUT/'figure_manifest.json').write_text(json.dumps({'schema_version':'1.0','font_families':['DejaVu Sans','Microsoft YaHei'],
 'script_sha256':digest(Path(__file__)),'figures':manifest,'visual_qa':'pending'},ensure_ascii=False,indent=2)+'\n',encoding='utf8')
print(json.dumps({'files':len(manifest),'figures':5,'languages':2,'png_shape':[2700,1260]}))
