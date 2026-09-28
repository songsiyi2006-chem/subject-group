"""Publication figures for hypothetical network, input scenarios and cached ML tests."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
COLORS = ['#075985', '#d97706', '#737a46']


def main():
    out = ROOT / 'reports/figures'
    out.mkdir(parents=True, exist_ok=True)
    def frame(name):
        return pd.read_csv(ROOT / ('results/' + name))
    pool = frame('reaction_network/parameter_pool.csv')
    sensitivity = frame('uncertainty/sobol_prefix_estimates.csv')
    learning = frame('benchmark_extension/learning_curve_summary.csv')
    paired = frame('benchmark_extension/paired_HV_differences.csv')
    holdout = frame('benchmark_extension/holdout_grouped_metrics.csv')
    manifest = []
    for zh in [False, True]:
        lang = 'chinese' if zh else 'english'
        pick = lambda en, cn: cn if zh else en
        plt.rcParams.update({'font.family': 'Microsoft YaHei' if zh else 'Arial', 'font.size': 8,
                             'axes.titlesize': 9, 'axes.labelsize': 8, 'xtick.labelsize': 7,
                             'ytick.labelsize': 7, 'axes.spines.top': False, 'axes.spines.right': False,
                             'axes.unicode_minus': False, 'svg.fonttype': 'none',
                             'svg.hashsalt': 'electratwin-extensions', 'figure.facecolor': 'white',
                             'axes.linewidth': .7})
        def save(fig, name, sources):
            for ext in ['png', 'svg']:
                p = out / f'{name}_{lang}.{ext}'
                fig.savefig(p, dpi=300, **({'metadata': {'Date': None}} if ext == 'svg' else {}))
                row = dict(file=p.name, format=ext, sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                           evidence='Hypothetical uncalibrated numerical models; no experimental validation',
                           sources={s: hashlib.sha256((ROOT/'results'/s).read_bytes()).hexdigest() for s in sources})
                if ext == 'png':
                    with Image.open(p) as im:
                        row.update(pixels=list(im.size), dpi=list(im.info['dpi']))
                else:
                    text = p.read_text(encoding='utf-8')
                    row.update(editable_text='<text' in text, embedded_raster='<image' in text)
                manifest.append(row)
            plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.8), layout='constrained')
        matrix = pool.pivot(index='flow_rate_uL_min', columns='potential_index_V', values='net_P_faradaic_efficiency_pct')
        im = axes[0].pcolormesh(np.linspace(.2-(.75-.2)/12, .75+(.75-.2)/12, 8),
                                np.linspace(100-(1500-100)/12, 1500+(1500-100)/12, 8),
                                matrix.to_numpy(), cmap='cividis', vmin=0, vmax=100, rasterized=False)
        im.set_edgecolor('face')
        cb = fig.colorbar(im, ax=axes[0], orientation='horizontal', pad=.2, label=pick('Net P FE / %', '净 P FE / %'))
        cb.solids.set_rasterized(False); cb.solids.set_edgecolor('face')
        axes[0].set(xlabel=pick('Imposed potential index / V', '给定电位指标 / V'),
                    ylabel=pick('Flow / uL/min', '流量 / 微升每分'), xlim=(.2,.75), ylim=(100,1500),
                    title=pick('a  49-point network pool', 'a  49 点反应网络池'))
        middle = pool[np.isclose(pool.flow_rate_uL_min, 800)].sort_values('potential_index_V')
        for key, label, style, color in [('conversion_A_pct', pick('A conversion', 'A 转化率'), 'o-', COLORS[0]),
                                          ('net_P_yield_pct', pick('Net P yield', '净 P 收率'), 's--', COLORS[1]),
                                          ('net_P_faradaic_efficiency_pct', 'Net P FE' if not zh else '净 P FE', '^:', COLORS[2])]:
            axes[1].plot(middle.potential_index_V, middle[key], style, color=color, lw=1.3, ms=3, label=label)
        axes[1].set(xlabel=pick('Imposed potential index / V', '给定电位指标 / V'), ylabel='%', ylim=(0, 105),
                    title=pick('b  Slice: 800 uL/min', 'b  截面：800 微升每分'))
        axes[1].legend(fontsize=7, loc='center left')
        fig.suptitle(pick('ASSUMED NETWORK: A to P / A to B / P to D', '假设反应网络：A 到 P / A 到 B / P 到 D'), fontsize=10)
        save(fig, 'Fig4_Reaction_Network', ['reaction_network/parameter_pool.csv'])

        fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.8), layout='constrained')
        parameters = ['flow_rate_uL_min', 'overpotential_V', 'diffusivity_m2_s', 'current_scale_A_m2', 'height_m']
        names = ['Q', 'eta', 'D', 'j-scale', 'H']
        for ax, target, title in zip(axes, ['conversion_pct', 'STY_assumed_product_kg_m3_day'],
                                    [pick('a  Conversion sensitivity', 'a  转化率敏感性'), pick('b  STY sensitivity', 'b  STY 敏感性')]):
            records = sensitivity[(sensitivity.base_N == 256) & (sensitivity.output == target)].set_index('parameter').loc[parameters]
            xx = np.arange(5)
            ax.bar(xx-.18, records.first_order_S, width=.36, color=COLORS[0], label=pick('First order S', '一阶 S'))
            ax.bar(xx+.18, records.total_order_ST, width=.36, color=COLORS[1], label=pick('Total effect ST', '总效应 ST'))
            ax.set_xticks(xx, names); ax.set(ylabel=pick('Estimated variance fraction', '估计的方差占比'), title=title, ylim=(-.02, .9))
            ax.axhline(0, c='#555555', lw=.5); ax.legend(fontsize=7)
        fig.suptitle(pick('ASSUMED INPUTS: N = 256; finite-sample S may exceed ST', '假设输入分布：N = 256；有限样本估计可能 S > ST'), fontsize=10)
        save(fig, 'Fig5_Parameter_Sensitivity', ['uncertainty/sobol_prefix_estimates.csv'])

        fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.8), layout='constrained')
        methods = ['gp_mc_ehvi', 'random', 'maximin_spacefill']
        labels = ['GP MC-EHVI', pick('Random', '随机'), 'Maximin']
        for method, label, color, marker in zip(methods, labels, COLORS, ['o-', 's--', '^:']):
            r = learning[learning.method == method].sort_values('budget')
            axes[0].plot(r.budget, r.mean_full_pool_HV_fraction, marker, color=color, lw=1.3, ms=3, label=label)
        axes[0].set(xlabel=pick('Budget including 5 initial points', '含 5 个初始点的预算'),
                    ylabel=pick('Mean full-pool HV fraction', '平均全池超体积覆盖比例'), ylim=(.85, 1.005),
                    title=pick('a  64 seeds, cached evaluations', 'a  64 种子，缓存评价'))
        axes[0].set_xticks([5, 10, 15, 20, 25]); axes[0].legend(fontsize=7, loc='lower right')
        for comparator, label, color, shift, marker in [('random', pick('GP minus random', 'GP 减随机'), COLORS[0], -.3, 'o'),
                                                       ('maximin_spacefill', pick('GP minus maximin', 'GP 减 Maximin'), COLORS[1], .3, 's')]:
            r = paired[(paired.left_method == 'gp_mc_ehvi') & (paired.right_method == comparator) & (paired.budget > 5)]
            axes[1].errorbar(r.budget+shift, r['mean'], yerr=[r['mean']-r.lower_95pct_t, r.upper_95pct_t-r['mean']],
                             fmt=marker+'-', color=color, ms=3, lw=1, capsize=2, label=label)
        axes[1].axhline(0, color='#555555', lw=.8, ls=':'); axes[1].set_xticks([10,15,20,25])
        axes[1].set(xlabel=pick('Budget', '预算'), ylabel=pick('Paired mean raw HV difference', '配对原始超体积均值差'),
                    title=pick('b  Descriptive 95% t intervals', 'b  描述性 95% t 区间'))
        axes[1].legend(fontsize=6.5, loc='upper right')
        fig.suptitle(pick('SIMULATION ONLY: 81-point pool; seed effects, not physical uncertainty', '仅模拟：81 点候选池；种子差异并非物理不确定性'), fontsize=9.4)
        save(fig, 'Fig6_Optimization_Extension', ['benchmark_extension/learning_curve_summary.csv', 'benchmark_extension/paired_HV_differences.csv'])

        fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.8), layout='constrained')
        groups = ['random_54_27', 'block_flow', 'block_eta']
        xlabels = [pick('Random', '随机'), pick('Flow block', '流量分块'), pick('Potential block', '电位分块')]
        xx = np.arange(3)
        for i, (model, label, color) in enumerate(zip(['fixed_matern_gp', 'training_mean', 'quadratic_regression'],
                                                    ['GP', pick('Mean', '均值'), pick('Quadratic', '二次回归')], COLORS)):
            r = holdout[(holdout.model == model) & (holdout.target == 'STY_div_50000')].set_index('split_kind').loc[groups]
            axes[0].bar(xx+(i-1)*.24, r.RMSE, width=.24, color=color, label=label)
        axes[0].set_xticks(xx, xlabels); axes[0].set(ylabel='RMSE (STY / 50000)', title=pick('a  Heldout predictive error', 'a  留出预测误差')); axes[0].legend(fontsize=7)
        for i, (target, label, color) in enumerate(zip(['STY_div_50000', 'negative_SEC_kwh_kg'], ['STY / 50000', '-SEC'], COLORS[:2])):
            r = holdout[(holdout.model == 'fixed_matern_gp') & (holdout.target == target)].set_index('split_kind').loc[groups]
            axes[1].bar(xx+(i-.5)*.3, 100*r.latent_interval_coverage_95pct, width=.3, color=color, label=label)
        axes[1].axhline(95, ls=':', c='#555555', lw=.8); axes[1].set_xticks(xx, xlabels)
        axes[1].set(ylabel=pick('GP latent 95% interval coverage / %', 'GP 潜函数 95% 区间覆盖率 / %'), ylim=(0, 110),
                    title=pick('b  Nominal coverage: 95%', 'b  名义覆盖率：95%')); axes[1].legend(fontsize=7, loc='upper right')
        fig.suptitle(pick('MODEL HOLDOUT: 26 splits; pooled uses are not independent experiments', '模型留出：26 个拆分；合并预测并非独立实验'), fontsize=9.4)
        save(fig, 'Fig7_Heldout_Prediction', ['benchmark_extension/holdout_grouped_metrics.csv'])
    record = dict(figure_count=4, languages=['english','chinese'], formats=['png','svg'], figures=manifest)
    (ROOT/'results/extension_figure_manifest.json').write_text(json.dumps(record, indent=2)+'\n', encoding='utf-8', newline='\n')
    print('Created 8 PNG and 8 editable vector SVG extension figures.')


if __name__ == '__main__':
    main()
