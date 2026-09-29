"""Populate bilingual tables directly from frozen-data post hoc CSV files."""
import csv
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
def read(name):
    with (BASE/'results'/name).open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))
def table(caption, headers, rows):
    return '**'+caption+'**\n\n'+'| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+''.join('| '+' | '.join(map(str,r))+' |\n' for r in rows)
def fmt(x,n=6):return f'{float(x):.{n}f}'
for lang in ('english','chinese'):
    zh=lang=='chinese';tables={}
    rows=read('hydration_test_metrics.csv')
    tables[1]=table('Table M1. '+('相同测试分子的水合自由能预测。单位：kcal mol⁻¹。' if zh else 'Hydration free-energy prediction on identical test molecules. Units: kcal mol⁻¹.'),['模型' if zh else 'Model','RMSE','MAE'],[[r['model']+((' / '+r['seed']) if r.get('seed') else ''),fmt(r['RMSE_kcal_mol']),fmt(r['MAE_kcal_mol'])] for r in rows])
    rows=read('network_observables.csv')
    names=['A 转化率','P 净产率','P 净选择性','P 净法拉第效率','A→P 毛电荷占比'] if zh else ['A conversion','Net P yield','Net P selectivity','Net P Faradaic efficiency','Gross A→P charge fraction']
    tables[2]=table('Table M2. '+('四物种模型的不同观测量。' if zh else 'Distinct observables of the four-species model.'),['观测量' if zh else 'Observable','%'],[[n,fmt(r['value'])] for n,r in zip(names,rows)])
    rows=read('optimization_final.csv')
    tables[3]=table('Table M3. '+('64 个种子的最终归一化超体积（预算 25）。SD 为种子间标准差。' if zh else 'Final normalized hypervolume over 64 seeds at budget 25. SD describes variation across seeds.'),['策略' if zh else 'Policy','均值' if zh else 'Mean','SD','最小值' if zh else 'Minimum','最大值' if zh else 'Maximum'],[[r['method'],*[fmt(r[k]) for k in ['mean_HV_fraction','sd_across_seeds','min_HV_fraction','max_HV_fraction']]] for r in rows])
    rows=read('kinetic_potential_comparison.csv')
    tables[4]=table('Table M4. '+('各设定过电位下 32 条轨迹；TOF 单位为 s⁻¹。' if zh else '32 trajectories at each assigned overpotential; TOF in s⁻¹.'),['η / V','SSA 均值' if zh else 'SSA mean','SD','有限窗 CTMC' if zh else 'Finite-window CTMC'],[[fmt(r['eta_V'],1),*[fmt(r[k]) for k in ['mean_TOF_s_1','sd_TOF_s_1','finite_window_CTMC_TOF_s_1']]] for r in rows])
    rows=read('thermostat_comparison.csv')
    tables[5]=table('Table M5. '+('每种恒温器 8 个重复的温度矩。' if zh else 'Temperature moments for eight replicas per thermostat.'),['恒温器' if zh else 'Thermostat','均值 / K' if zh else 'Mean / K','方差 / K²' if zh else 'Variance / K²','方差比' if zh else 'Variance ratio'],[[r['thermostat'],fmt(r['replicate_mean_temperature_K']),fmt(r['mean_within_trajectory_temperature_variance_K2']),f"{float(r['variance_ratio_to_canonical']):.6g}"] for r in rows])
    rows=read('chromatography_grouped.csv')
    tables[6]=table('Table M6. '+('合成色谱的平均绝对峰面积误差；每行 10 个案例。' if zh else 'Mean absolute peak-area error in synthetic chromatography; ten cases per row.'),['峰宽误设' if zh else 'Wrong width','间距 / min' if zh else 'Separation / min','拟合误差 / %' if zh else 'Fitted error / %','积分误差 / %' if zh else 'Window error / %'],[[r['width_mismatch'],fmt(r['separation_min'],2),fmt(r['mean_absolute_fitted_error_pct']),fmt(r['mean_absolute_window_error_pct'])] for r in rows])
    path=BASE/f'multiscale_{lang}.md'
    text=path.read_text('utf-8')
    for number,value in tables.items():text=text.replace('{{TABLE_M'+str(number)+'}}',value.rstrip())
    path.write_bytes(text.replace('\r\n','\n').encode('utf-8'))
print('Six CSV-backed tables populated in each language.')
