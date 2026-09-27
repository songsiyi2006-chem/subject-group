"""Build separate and section-aligned reports from a single set of saved tables."""
import json,re
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];A=ROOT/'results/audit';REPORT=ROOT/'reports'
def table(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(str(v) for v in row)+' |' for row in rows])
def tokens(zh):
    data=json.loads((ROOT/'results/original/research_grade_results.json').read_text(encoding='utf-8'));summary=json.loads((A/'audit_summary.json').read_text(encoding='utf-8'))
    d={};names=['褪黑素','咖啡因','2-苯基喹啉','色醇','2-(噻吩-2-基)乙醇','苯并呋喃-2-羧酸乙酯','二氢吲哚','咔唑'] if zh else ['Melatonin','Caffeine','2-Phenylquinoline','Tryptophol','2-(Thiophen-2-yl)ethanol','Ethyl benzofuran-2-carboxylate','Indoline','Carbazole']
    def put(key,en,cn,rows):d[key]=table(cn if zh else en,rows)
    p=pd.read_csv(A/'pareto_complete.csv')
    put('PARETO',['ID','j (mA/cm²)','Catalyst (mol%)','Stirring (rpm)','U (V)','Y (%)','FE (%)','SEC (kWh/kg)','Green score'],['ID','j (mA/cm²)','催化剂 (mol%)','搅拌 (rpm)','U (V)','Y (%)','FE (%)','SEC (kWh/kg)','绿色评分'],[[r.condition_id,f'{r.j_mA_cm2:.2f}',f'{r.catalyst_mol_pct:.1f}',r.stirring_rpm,f'{r.U_cell_V:.3f}',f'{r.Yield_pct:.2f}',f'{r.FE_pct:.2f}',f'{r.SEC_kWh_kg:.3f}',f'{r.Green_Score:.6f}'] for r in p.itertuples()])
    w=pd.read_csv(A/'scalarization_choices.csv').sort_values('j_mA_cm2')
    put('WEIGHTS',['ID','j (mA/cm²)','Selections / 231'],['ID','j (mA/cm²)','被选择次数 / 231'],[[r.condition_id,f'{r.j_mA_cm2:.2f}',r.count] for r in w.itertuples()])
    c=summary['multiobjective']['constraint_scenarios']
    put('CONSTRAINTS',['Minimum Y (%)','Minimum FE (%)','Maximum SEC (kWh/kg)','Feasible rows','Best score ID'],['最低 Y (%)','最低 FE (%)','最高 SEC (kWh/kg)','可行记录','最高评分 ID'],[[r['min_yield'],r['min_FE'],r['max_SEC'],r['feasible_count'],r['best_condition_id'] if r['best_condition_id'] is not None else '—'] for r in c])
    put('AREA',['Area (cm²)','Selected j (mA/cm²)','U (V)','SEC (kWh/kg)'],['面积 (cm²)','所选 j (mA/cm²)','U (V)','SEC (kWh/kg)'],[[r['area_cm2'],f"{r['j_mA_cm2']:.4f}",f"{r['U_cell_V']:.4f}",f"{r['SEC_kWh_kg']:.4f}"] for r in summary['multiobjective']['area_sensitivity']])
    put('GP',['Model','Objective','Macro MAE','Macro RMSE'],['模型','目标','宏平均 MAE','宏平均 RMSE'],[[r['model'],r['objective'],f"{r['MAE']:.6f}",f"{r['RMSE']:.6f}"] for r in summary['multiobjective']['gp_macro_metrics']])
    desc=pd.read_csv(A/'scope_descriptors.csv');source=data['Module2_Substrate_Scope']['scope_summary_table']
    put('SCOPE',['Substrate','Formula','MW (g/mol)','logP','Assigned Eox proxy*','Atom index','Random yield label (%)**','Source category'],['底物','分子式','MW (g/mol)','logP','指定 Eox 代理量*','原子编号','随机产率标签 (%)**','源分类'],[[names[i],desc.iloc[i].formula,f"{r['molecular_weight']:.2f}",f"{r['calculated_logP']:.2f}",f"{r['estimated_E_ox_V_vs_SCE']:.3f}",r['most_reactive_carbon_idx'],f"{r['predicted_electrochemical_yield_pct']:.1f}",('高' if zh else 'High') if r['estimated_E_ox_V_vs_SCE']<1.45 else ('中' if zh else 'Moderate')] for i,r in enumerate(source)])
    put('LOCAL',['Substrate','TPSA (Å²)','Selected C charge (e)','Hybridization','Source MMFF status'],['底物','TPSA (Å²)','所选 C 电荷 (e)','杂化','源 MMFF 状态'],[[names[i],f'{r.TPSA:.2f}',f'{r.primary_q:.6f}',r.primary_hybridization,r.source_mmff_status] for i,r in enumerate(desc.itertuples())])
    xtb=pd.read_csv(A/'scope_xtb.csv')
    put('XTB',['Substrate','Neutral energy (Eh)','Cation energy (Eh)','Removal difference (eV)','Source-site charge response (e)'],['底物','中性能量 (Eh)','阳离子能量 (Eh)','移电子能量差 (eV)','源位点电荷响应 (e)'],[[names[i],f'{r.neutral_Eh:.12f}',f'{r.cation_Eh:.12f}',f'{r.removal_difference_eV:.6f}',f'{r.source_site_removal_response:.8f}'] for i,r in enumerate(xtb.itertuples())])
    rand=pd.read_csv(A/'scope_random_yield_summary.csv').set_index('molecule')
    put('RANDOM',['Substrate','Mean of 1,000 random labels (%)','SD (percentage points)'],['底物','1,000 次随机标签均值 (%)','SD（百分点）'],[[names[i],f"{rand.loc[r.molecule,'mean']:.4f}",f"{rand.loc[r.molecule,'std']:.4f}"] for i,r in enumerate(desc.itertuples())])
    m=data['Module3_POP_Pore_Transport']['diffusion_kinetic_evaluation'];rows=[]
    for regime,v in m.items():
        for r in v['transport_sweep']:rows.append([regime,r['pellet_radius_um'],f"{r['thiele_modulus_phi']:.3f}",f"{r['internal_effectiveness_factor_eta']:.4f}",f"{r['pore_diffusion_utilization_pct']:.2f}"])
    put('TRANSPORT',['Scenario','Radius (µm)','φ','η','Internal utilization (%)'],['情景','半径 (µm)','φ','η','内部利用率 (%)'],rows)
    put('THRESHOLDS',['Scenario','Target η','φ at threshold','Radius (µm)'],['情景','目标 η','阈值 φ','半径 (µm)'],[[r['regime'],r['target_eta'],f"{r['phi']:.6f}",f"{r['radius_um']:.6f}"] for r in summary['transport']['thresholds']])
    e=pd.read_csv(A/'external_film_scenarios.csv');e=e[(e.radius_um==200)&(e.Bi.isin([.1,1,10,100]))]
    put('EXTERNAL',['Scenario (R=200 µm)','Bi','Surface / bulk concentration','Overall effectiveness'],['情景（R=200 µm）','Bi','表面／主体浓度','总有效性因子'],[[r.regime,r.Bi,f'{r.Cs_over_Cbulk:.6f}',f'{r.eta_overall:.6f}'] for r in e.itertuples()])
    d['SUFFIX']='_zh' if zh else '';d['SHA']=json.loads((ROOT/'source/source_record.json').read_text())['script_sha256']
    return d
def main():
    texts={}
    for lang,zh in [('english',False),('chinese',True)]:
        text=(REPORT/'templates'/f'{lang}.md').read_text(encoding='utf-8')
        for key,val in tokens(zh).items():text=text.replace('@@'+key+'@@',val)
        assert not re.search(r'@@[A-Z_]+@@',text)
        texts[lang]=text.rstrip()+'\n';(REPORT/f'research_grade_technical_report_{lang}.md').write_text(texts[lang],encoding='utf-8',newline='\n')
    en=re.split(r'^## ',texts['english'],flags=re.M);zh=re.split(r'^## ',texts['chinese'],flags=re.M);assert len(en)==len(zh)==6
    combined='# Comprehensive Technical Report: Multiscale Machine Learning & Chemical Engineering for Green Organic Electrocatalysis\n# 综合技术报告：面向绿色有机电催化的多尺度机器学习与化学工程集成研究\n\n2026-09-27 · Section-aligned bilingual edition / 按章节对齐的双语版\n\n'
    for ep,zp in zip(en[1:],zh[1:]):
        eh,eb=ep.split('\n',1);zhh,zb=zp.split('\n',1);combined+='\n## '+eh+' / '+zhh+'\n\n'
        es=re.split(r'^### ',eb,flags=re.M);zs=re.split(r'^### ',zb,flags=re.M);assert len(es)==len(zs)
        combined+='**English**\n'+es[0]+'\n**中文**\n'+zs[0]
        for e,z in zip(es[1:],zs[1:]):
            eh,eb=e.split('\n',1);zhh,zb=z.split('\n',1);combined+='\n### '+eh+' / '+zhh+'\n\n**English**\n'+eb+'\n**中文**\n'+zb
    (REPORT/'research_grade_technical_report_bilingual.md').write_text(combined.rstrip()+'\n',encoding='utf-8',newline='\n');print('Built three research reports.')
if __name__=='__main__':main()
