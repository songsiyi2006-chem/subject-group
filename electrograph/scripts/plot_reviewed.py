"""Render source-backed research figures from saved outputs; perform no new fits."""
from collections import Counter
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
from rdkit import Chem

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / 'reports/figures'
BLUE, GOLD, OLIVE, GRAY, PURPLE = '#0B5D87', '#CF7900', '#757B43', '#696969', '#8466A6'
FILES = []


def rows(name):
    with (BASE / name).open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def settings(chinese):
    plt.rcParams.update({'font.family': 'Microsoft YaHei' if chinese else 'Arial',
        'font.size': 9, 'axes.titlesize': 10.5, 'axes.labelsize': 9, 'legend.fontsize': 7.5,
        'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .8,
        'svg.fonttype': 'none', 'svg.hashsalt': 'electrograph-reviewed-20260928',
        'axes.unicode_minus': False, 'figure.facecolor': 'white', 'savefig.facecolor': 'white'})


def canvas(title):
    fig, axes = plt.subplots(1, 2, figsize=(8, 4.6), layout='constrained')
    fig.suptitle(title, fontsize=11.5)
    return fig, axes


def save(fig, stem, lang, sources):
    source_hashes = {name: sha(BASE / name) for name in sources}
    for suffix in ['png', 'svg']:
        path = OUT / f'{stem}_{lang}.{suffix}'
        fig.savefig(path, dpi=300, metadata={'Date': None} if suffix == 'svg' else {})
        item = dict(file=path.name, format=suffix, sha256=sha(path), sources=source_hashes,
                    evidence='Saved source audit, benchmark or assumed model; no electrochemical validation')
        if suffix == 'png':
            with Image.open(path) as im:
                item.update(pixels=list(im.size), dpi=list(im.info['dpi']))
        else:
            tree = ET.parse(path); ns = {'s': 'http://www.w3.org/2000/svg'}
            item.update(editable_text=bool(tree.findall('.//s:text', ns)),
                        embedded_raster=bool(tree.findall('.//s:image', ns)))
        FILES.append(item)
    plt.close(fig)


def source_audit(lang):
    cn = lang == 'chinese'
    fig, (a, b) = canvas('源码审计：硬编码曲线与未训练网络' if cn else 'SOURCE AUDIT: hardcoded curve and untrained network')
    data = rows('results/hardcoded_curve_audit.csv')
    x = [float(r['eta_V']) for r in data]
    a.plot(x, [float(r['source_hardcoded_TOF_s_1']) for r in data], 's-', color=GOLD,
           label='源码硬编码列表' if cn else 'Hardcoded source list')
    a.plot(x, [float(r['analytical_independent_cycle_TOF_s_1']) for r in data], 'o--', color=BLUE,
           label='给定速率的解析期望' if cn else 'Analytic rate-model expectation')
    a.set(xlabel='给定电势指标 / V' if cn else 'Imposed potential index / V', ylabel='TOF / s$^{-1}$',
          title='a  硬编码曲线不符合给定速率' if cn else 'a  Hardcoded curve versus specified rates', ylim=(0, 280))
    a.legend(loc='upper left', frameon=False)
    names_en = ['Caffeine', 'Melatonin', 'Tryptophol', 'Phenylquinoline', 'Benzofuran ester', 'Indoline']
    names_cn = ['咖啡因', '褪黑素', '色醇', '苯基喹啉', '苯并呋喃酯', '吲哚啉']
    smiles = ['Cn1cnc2c1c(=O)n(C)c(=O)n2C', 'CC(=O)NCCC1=CNc2c1cc(OC)cc2',
              'OCCc1c[nH]c2ccccc12', 'c1ccc(cc1)-c2ccc3ccccc3n2', 'CCOC(=O)c1cc2ccccc2o1', 'c1ccc2c(c1)CCN2']
    frequencies = Counter(r['smiles'] for r in rows('results/random_priority_probes.csv') if int(r['rank']) == 1)
    bars = b.barh(np.arange(6), [frequencies[s] for s in smiles], color=BLUE)
    b.bar_label(bars, padding=3, fontsize=8)
    b.set_yticks(np.arange(6), names_cn if cn else names_en); b.invert_yaxis()
    b.set(xlim=(0, 10), xlabel='32 次随机优先级中的首选次数' if cn else 'Top-choice count over 32 priority seeds',
          title='b  网络固定，仅随机探索项改变' if cn else 'b  Fixed network; only uniform term varies')
    save(fig, 'Fig1_Source_Audit', lang, ['results/hardcoded_curve_audit.csv', 'results/random_priority_probes.csv'])


def graph_benchmark(lang):
    cn = lang == 'chinese'
    fig, (a, b) = canvas('FreeSolv 水合自由能基准：不是氧化电位' if cn else 'FreeSolv HYDRATION BENCHMARK: not oxidation potential')
    data = {r['model']: r for r in rows('results/learning/regression_metrics.csv') if r['split'] == 'test'}
    models = ['mpnn_seed20260928', 'mpnn_seed20260929', 'mpnn_seed20260930', 'descriptor_ridge', 'training_mean', 'shuffled_train_mpnn_seed20261001']
    labels = ['MPNN 28', 'MPNN 29', 'MPNN 30', '岭回归' if cn else 'Ridge', '训练均值' if cn else 'Mean', '置乱训练' if cn else 'Shuffled']
    values = [float(data[m]['RMSE_kcal_mol']) for m in models]
    bars = a.barh(range(6), values, color=[BLUE]*3+[GOLD, GRAY, OLIVE])
    a.bar_label(bars, labels=[f'{v:.3f}' for v in values], padding=3, fontsize=8)
    a.set_yticks(range(6), labels); a.invert_yaxis()
    a.set(xlim=(0, 5.7), xlabel='测试 RMSE / kcal mol$^{-1}$' if cn else 'Test RMSE / kcal mol$^{-1}$',
          title='a  固定 51 分子测试集' if cn else 'a  Fixed 51-molecule test set')
    curves = rows('results/learning/campaign_learning_curves.csv')
    for method, label, color, marker, style in [
        ('latent_gp_ucb', '学习特征 GP' if cn else 'Learned-feature GP', BLUE, 'o', '-'),
        ('descriptor_gp_ucb', '描述符 GP' if cn else 'Descriptor GP', GOLD, 's', '--'),
        ('random', '随机' if cn else 'Random', OLIVE, '^', ':')]:
        run = sorted([r for r in curves if r['method'] == method], key=lambda r: int(r['step']))
        if not run:
            raise ValueError('Missing method ' + method)
        b.plot([int(r['step']) for r in run], [float(r['mean_simple_regret_kcal_mol']) for r in run],
               color=color, marker=marker, linestyle=style, markersize=3, label=label)
    b.set(xlabel='缓存查询次数（含 4 个初始点）' if cn else 'Cached queries (including 4 initial points)',
          ylabel='均值 simple regret / kcal mol$^{-1}$' if cn else 'Mean simple regret / kcal mol$^{-1}$',
          title='b  每方法 8 个随机种子' if cn else 'b  Eight search seeds per method', ylim=(-.5, 26), xlim=(1,16))
    b.set_xticks([1,4,8,12,16]); b.legend(frameon=False, loc='upper right')
    save(fig, 'Fig2_Graph_Benchmark', lang, ['results/learning/regression_metrics.csv', 'results/learning/campaign_learning_curves.csv'])


def conformers(lang):
    cn = lang == 'chinese'
    fig = plt.figure(figsize=(8, 4.6), layout='constrained')
    fig.suptitle('MMFF94 构象几何：权重不是溶液自由能布居' if cn else 'MMFF94 GEOMETRIES: minima weights are not solution populations', fontsize=11.5)
    a = fig.add_subplot(1,2,1,projection='3d'); b = fig.add_subplot(1,2,2)
    sdf = BASE / 'results/structure/target_optimized.sdf'
    mols = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=False) if m is not None]
    mol = min(mols, key=lambda m: float(m.GetProp('energy_kcal_mol')))
    positions = mol.GetConformer().GetPositions(); positions -= positions.mean(axis=0)
    for bond in mol.GetBonds():
        p = positions[[bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]]
        a.plot(p[:,0],p[:,1],p[:,2],color='#808080',linewidth=1,zorder=1)
    colors = {'C':'#414141','H':'#F5F5F5','N':BLUE,'O':GOLD}
    for element in ['C','H','N','O']:
        selected = np.array([positions[i] for i,atom in enumerate(mol.GetAtoms()) if atom.GetSymbol()==element])
        if len(selected):
            a.scatter(*selected.T,s=20 if element=='H' else 55,c=colors[element],edgecolors='#444444',
                      linewidths=.4,depthshade=False,label=element)
    a.set_proj_type('ortho'); a.set_box_aspect([1,1,.55]); a.view_init(elev=22, azim=35)
    extent = float(np.max(np.ptp(positions,axis=0))) / 2 + .4
    a.set(xlim=(-extent,extent),ylim=(-extent,extent),zlim=(-extent*.55,extent*.55))
    a.set_axis_off(); a.legend(loc='lower center', bbox_to_anchor=(.5,.10), ncol=4, frameon=False,handletextpad=.2,columnspacing=.8)
    a.set_title('a  目标最低采样能量几何' if cn else 'a  Lowest sampled target geometry',pad=8)
    a.text2D(.5,.02,'连接线省略键级；H 为显式原子' if cn else 'Connectivity shown; bond orders omitted',transform=a.transAxes,ha='center',fontsize=7.5)
    stats = [r for r in rows('results/structure/ensemble_statistics.csv') if r['scope']=='pooled' and float(r['temperature_K'])==298.15]
    order = ['target','melatonin','caffeine','phenylquinoline','tryptophol','indoline','benzofuran_ester']
    lookup = {r['molecule_id']:r for r in stats}
    labels_cn = ['目标','褪黑素','咖啡因','苯基喹啉','色醇','吲哚啉','苯并呋喃酯']
    labels_en = ['Target','Melatonin','Caffeine','Phenylquinoline','Tryptophol','Indoline','Benzofuran ester']
    values = [float(lookup[m]['weighted_sasa_A2']) for m in order]
    bars = b.barh(range(7),values,color=BLUE)
    b.bar_label(bars,labels=[f'{v:.1f}' for v in values],padding=3,fontsize=7.5)
    b.set_yticks(range(7),labels_cn if cn else labels_en); b.invert_yaxis()
    b.set(xlim=(0,550),xlabel='298.15 K 加权 SASA / Å$^2$' if cn else '298.15 K weighted SASA / Å$^2$',
          title='b  正半径与 24 个方向平均' if cn else 'b  Positive radii; 24-orientation average')
    save(fig,'Fig3_Conformer_Ensemble',lang,['results/structure/target_optimized.sdf','results/structure/ensemble_statistics.csv'])


def kinetics(lang):
    cn = lang == 'chinese'
    fig,(a,b)=canvas('假定位点动力学：数值核验，不是机理验证' if cn else 'ASSUMED SITE KINETICS: numerical, not mechanistic validation')
    scan=rows('results/kinetics/potential_summary.csv')
    x=np.array([float(r['eta_V']) for r in scan]);y=np.array([float(r['mean_TOF_s_1']) for r in scan])
    lo=np.array([float(r['mean_MC_lower_95pct_t']) for r in scan]);hi=np.array([float(r['mean_MC_upper_95pct_t']) for r in scan])
    a.errorbar(x,y,yerr=np.array([y-lo,hi-y]),fmt='o',capsize=3,color=BLUE,label='SSA 均值（32 种子）' if cn else 'SSA mean (32 seeds)')
    a.plot(x,[float(r['finite_window_CTMC_TOF_s_1']) for r in scan],'s--',color=GOLD,markersize=3,label='有限时域 CTMC' if cn else 'Finite-window CTMC')
    a.set(xlabel='给定电势指标 / V' if cn else 'Imposed potential index / V',ylabel='TOF / s$^{-1}$',
          title='a  均值的 95% Monte Carlo 区间' if cn else 'a  95% Monte Carlo intervals for mean',ylim=(26.25,26.85))
    a.legend(loc='upper right',frameon=False)
    trace=rows('results/kinetics/baseline_trajectory.csv')
    startup=[r for r in trace if float(r['time_s']) <= .2]
    sample=startup[::7]+[startup[-1]]
    analytic=rows('results/kinetics/analytic_transient.csv')
    states=['empty','substrate','radical','intermediate','product']
    labels=['空位','底物','自由基','中间体','产物'] if cn else ['Empty','Substrate','Radical','Intermediate','Product']
    for state,label,color in zip(states,labels,[BLUE,GRAY,GOLD,PURPLE,OLIVE]):
        b.plot([float(r['time_s']) for r in sample],[int(r['n_'+state])/256 for r in sample],color=color,lw=.65,alpha=.65)
        b.plot([float(r['time_s']) for r in analytic],[float(r[state]) for r in analytic],color=color,ls='--',lw=1.4,label=label)
    b.axvline(.2,color=GRAY,ls=':',lw=.8)
    b.set(xlabel='时间 / s' if cn else 'Time / s',ylabel='位点占据分数' if cn else 'Site fraction',ylim=(-.02,1.04),xlim=(0,.2),
          title='b  启动期：SSA 实线；CTMC 虚线' if cn else 'b  Startup: SSA solid; CTMC dashed')
    b.legend(frameon=False,loc='upper right',ncol=2)
    save(fig,'Fig4_Kinetic_Validation',lang,['results/kinetics/potential_summary.csv','results/kinetics/baseline_trajectory.csv','results/kinetics/analytic_transient.csv'])


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    for lang in ['english','chinese']:
        settings(lang=='chinese')
        for plot in [source_audit,graph_benchmark,conformers,kinetics]:
            plot(lang)
    record=dict(figure_count=4,languages=['english','chinese'],formats=['png','svg'],figures=FILES)
    (BASE/'results/figure_manifest.json').write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n')
    print('Created eight PNG and eight editable vector SVG figures from saved outputs.')


if __name__=='__main__':
    main()
