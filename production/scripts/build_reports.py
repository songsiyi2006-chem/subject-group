"""Build separate and section-aligned bilingual reports from saved numerical records."""
import json,re
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'results/repaired';A=ROOT/'results/audit';REPORTS=ROOT/'reports'
def tab(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(str(x) for x in row)+' |' for row in rows])
def tokens(zh):
    d=json.loads((R/'production_benchmark_results.json').read_text());s=json.loads((A/'audit_summary.json').read_text());values={}
    def put(key,en,cn,rows):values[key]=tab(cn if zh else en,rows)
    pool=pd.read_csv(R/'task_a_candidates.csv')
    solvents=pool.drop_duplicates('solvent');elecs=pool.drop_duplicates('electrolyte')
    put('SOLVENTS',['Solvent','Assigned ε','Viscosity (mPa·s)','Donor number'],['溶剂','给定 ε','黏度 (mPa·s)','给体数'],[[r.solvent,r.feat_eps,r.feat_visc,r.feat_DN] for r in solvents.itertuples()])
    put('ELECTROLYTES',['Electrolyte','Assigned radius (pm)','Assigned oxidation limit (V vs SCE)'],['电解质','给定半径 (pm)','给定氧化阈值 (V vs SCE)'],[[r.electrolyte,r.feat_radius,r.feat_E_ox] for r in elecs.itertuples()])
    obs=pd.read_csv(R/'task_a_observations.csv')
    put('TRAJECTORY',['BO step','Solvent','Electrolyte','j (mA/cm²)','T (°C)','Observed yield (%)','Best so far (%)'],['BO 步数','溶剂','电解质','j (mA/cm²)','T (°C)','含噪产率 (%)','累计最好值 (%)'],[[r['iteration'],r['chosen_solvent'],r['chosen_electrolyte'],r['current_density'],f"{obs.iloc[5+r['iteration']].temp_K-273.15:.1f}",r['acquired_yield'],r['cumulative_best_yield']] for r in d['Task_A_Physical_BO']['optimization_trajectory']])
    put('BO_COMPARE',['Method','Mean regret','SD','Median'],['方法','遗憾值均值','SD','中位数'],[[v['method'],f"{v['mean']:.4f}",f"{v['std']:.4f}",f"{v['median']:.4f}"] for v in s['A']['methods']])
    put('BO_CI',['Paired contrast','Mean difference','95% bootstrap interval'],['配对差值','平均差','95% 自助区间'],[[k,f"{v['mean']:.4f}",f"[{v['ci95'][0]:.4f}, {v['ci95'][1]:.4f}]"] for k,v in s['A']['paired'].items()])
    labels={0:'methoxy-CH3',3:'C6',4:'C7',8:'C3',10:'C4',12:'phenyl:12',13:'phenyl:13',14:'phenyl:14',15:'phenyl:15',16:'phenyl:16'}
    sites=pd.read_csv(R/'task_b_all_ch_sites.csv')
    put('SITES',['Atom index (0-based)','Mapped label','Hybridization','Gasteiger q (e)','One H SASA (Å²)','Heuristic index'],['零基原子编号','映射标签','杂化','Gasteiger q (e)','一个 H 的 SASA (Å²)','启发式评分'],[[r.atom_index,labels[r.atom_index],r.hybridization,f'{r.gasteiger_charge:.4f}',f'{r.h_sasa_angstrom2:.3f}',f'{r.anodic_activation_index:.4f}'] for r in sites.itertuples()])
    site_summary=pd.read_csv(A/'task_b_site_summary.csv')
    put('CONFORMERS',['Mapped label','Mean score ± SD','Mean H SASA ± SD (Å²)','Highest-score count / 32'],['映射标签','评分均值 ± SD','H SASA 均值 ± SD (Å²)','评分第一次数 / 32'],[[r.site_label,f'{r.score_mean:.4f} ± {r.score_sd:.4f}',f'{r.sasa_mean:.3f} ± {r.sasa_sd:.3f}',r.top_count] for r in site_summary.itertuples()])
    put('SASA_SENS',['Probe (Å)','Algorithm','Total SASA (Å²)','Top atom index','Top score'],['探针 (Å)','算法','总 SASA (Å²)','最高分原子编号','最高分'],[[v['probe_A'],v['algorithm'],f"{v['total_sasa_A2']:.3f}",v['top_atom_index'],f"{v['top_score']:.4f}"] for v in s['B']['probe_algorithm_sensitivity']])
    put('XTB_ENERGY',['State','Charge / multiplicity','Energy (Eh)'],['状态','电荷／多重度','能量 (Eh)'],[[name,cm,f"{s['target_xTB']['energies_Eh'][name]:.12f}"] for name,cm in [('neutral','0 / 1'),('cation','+1 / 2'),('anion','−1 / 2')]])
    ranking=pd.read_csv(A/'task_c_ranking.csv')
    put('SAC_RANK',['Rank','Metal','Coordination label','d count','Assigned proxy (eV)','Synthetic score (kcal/mol)'],['排名','金属','配位标签','d 电子数','给定代理值 (eV)','合成评分 (kcal/mol)'],[[i+1,r.metal,r.coordination_type,r.d_electrons,f'{r.d_band_center_proxy_eV:.3f}',f'{r.activation_barrier_kcal_mol:.2f}'] for i,r in enumerate(ranking.itertuples())])
    put('IMPORTANCE',['Descriptor','Extra Trees importance'],['描述符','极端随机树重要性'],[[k,f'{v:.4f}'] for k,v in d['Task_C_SAC_Informatics']['descriptor_importances'].items()])
    put('HOLDOUT',['Held-out group','Model','Macro MAE','Macro RMSE'],['留出分组','模型','宏平均 MAE','宏平均 RMSE'],[[v['split'],v['model'],f"{v['test_mae']:.4f}",f"{v['test_rmse']:.4f}"] for v in s['C']['group_holdout_macro_means']])
    put('FLOW',['Flow (µL/min)','Velocity (mm/s)','Residence (s)','k_m (m/s)','Conversion (%)','STY (mmol/L/h)'],['流量 (µL/min)','速度 (mm/s)','停留时间 (s)','k_m (m/s)','转化率 (%)','STY (mmol/L/h)'],[[k.replace('_uL_min',''),v['flow_velocity_mm_s'],v['residence_time_sec'],v['mass_transfer_coeff_m_s'],v['conversion_pct'],v['space_time_yield_mmol_L_h']] for k,v in d['Task_D_Flow_Electrochemistry']['flow_sweep_kinetics'].items()])
    put('FLOW_AUDIT',['Flow (µL/min)','Damköhler number','Mass-transfer resistance fraction','Implied current (mA)','(D/k_m)/h'],['流量 (µL/min)','Damköhler 数','传质阻力占比','隐含电流 (mA)','(D/k_m)/h'],[[v['q_uL_min'],f"{v['Da']:.4f}",f"{v['mass_transfer_resistance_fraction']:.4f}",f"{v['current_from_molar_balance_A']*1000:.4f}",f"{v['delta_effective_over_height']:.4f}"] for v in s['D']['rows']])
    record=json.loads((ROOT/'source/source_record.json').read_text());values['SOURCE_SHA']=record['extracted_script_sha256']
    values['SUFFIX']='_zh' if zh else ''
    values['XTB_REMOVAL']=f"{s['target_xTB']['removal_difference_eV']:.4f}";values['XTB_ADDITION']=f"{s['target_xTB']['addition_difference_eV']:.4f}"
    return values

def main():
    texts={}
    for lang,zh in [('english',False),('chinese',True)]:
        text=(REPORTS/'templates'/f'{lang}.md').read_text(encoding='utf-8')
        for key,value in tokens(zh).items():text=text.replace('@@'+key+'@@',value)
        if re.search(r'@@[A-Z_]+@@',text):raise ValueError('Unresolved report placeholder')
        texts[lang]=text.rstrip()+'\n';(REPORTS/f'production_technical_report_{lang}.md').write_text(texts[lang],encoding='utf-8',newline='\n')
    # Match corresponding main sections and subsections, keeping only four numbered main headings.
    def main_parts(text):return re.split(r'^## ',text,flags=re.M)
    en,zh=main_parts(texts['english']),main_parts(texts['chinese']);assert len(en)==len(zh)
    combined='# Quantitative Organic Electrochemistry & Catalysis / 定量有机电合成与配位催化\n\nDate / 日期：2026-09-27。Audited, section-aligned bilingual edition / 经核验、按章节对齐的双语合并版。\n\n'
    # Keep each evidence statement but omit metadata describing a separate edition.
    combined+='\n\n'.join(en[0].split('\n\n')[2:])+'\n\n'+'\n\n'.join(zh[0].split('\n\n')[2:])
    for ep,zp in zip(en[1:],zh[1:]):
        eh,eb=ep.split('\n',1);zhh,zb=zp.split('\n',1)
        combined+='\n## '+eh+' / '+zhh+'\n\n'
        ea=re.split(r'^### ',eb,flags=re.M);za=re.split(r'^### ',zb,flags=re.M)
        assert len(ea)==len(za),(eh,len(ea),len(za))
        combined+=ea[0]+za[0]
        for epart,zpart in zip(ea[1:],za[1:]):
            ehead,econtent=epart.split('\n',1);zhead,zcontent=zpart.split('\n',1)
            combined+='\n### '+ehead+' / '+zhead+'\n\n**English**\n'+econtent+'\n**中文**\n'+zcontent
    (REPORTS/'production_technical_report_bilingual.md').write_text(combined.rstrip()+'\n',encoding='utf-8',newline='\n')
    print('Built English, Chinese, and section-aligned bilingual editions.')
if __name__=='__main__':main()
