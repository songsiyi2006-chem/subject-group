"""Build separate editions and a full section-aligned bilingual report from saved data."""
from pathlib import Path
import csv,json,re
ROOT=Path(__file__).resolve().parents[1]
def j(p):return json.loads((ROOT/p).read_text(encoding='utf-8'))
def rows(p):
    with (ROOT/p).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
def table(headers,data):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,r))+' |' for r in data])
def tokens(zh):
    choose=lambda en,cn:cn if zh else en
    al=j('results/audit/active_learning.json');sto=j('results/audit/stoichiometry.json')['base'];orig=j('results/original/closed_loop_results.json')['Active_Learning_Feedback']['round_2_recommended_condition'];dft=j('inputs/reviewed/input_manifest.json')
    t={'SOURCE_HASH':j('source/source_record.json')['script_sha256']}
    findings=[('Wet-lab feedback is experimental','Six hardcoded mock labels; no raw assay files','湿实验反馈为实验数据','实际为六条硬编码模拟标签，没有原始分析文件'),('MAE is cross-validation','2.883333 pp is a direct label comparison; LOO is computed separately','MAE 即交叉验证','2.883333 个百分点仅为给定标签比较；另做留一验证'),('Acquisition is EI','Source uses UCB; extension implements analytic EI','采集函数为 EI','原式为 UCB；扩展实现解析 EI'),('Shared DFT target and method','Indoline differs from target indole; source ORCA uses CPCM','DFT 目标和方法一致','吲哚啉不同于目标吲哚；原 ORCA 使用 CPCM'),('Ready-to-run SOP','Partner equivalents, product, quench and assay are unspecified','SOP 可直接执行','偶联试剂当量、产物、淬灭及分析方法未确定'),('Computed DFT energies','Input files only; no DFT executed','已得到 DFT 能量','仅生成输入，未执行 DFT')]
    t['FINDINGS']=table([choose('Source claim','原文主张'),choose('Audited status','核验状态')],[[r[2],r[3]] if zh else r[:2] for r in findings])
    t['REAGENTS']=table([choose('Material','物料'),'mmol','mg','mL',choose('Definition','定义')],[
        [choose('Target indole','目标吲哚'),'0.200',f'{sto["substrate_mg"]:.3f}',choose('N/A','不适用'),choose('1.00 equivalent; 0.033333 M nominal','1.00 当量；名义浓度 0.033333 M')],
        ['nBu4NPF6',f'{sto["electrolyte_mmol"]:.3f}',f'{sto["electrolyte_mg"]:.3f}',choose('N/A','不适用'),'0.100 M'],
        ['MeCN',choose('unknown','未知'),choose('unknown','未知'),'4.800',choose('80.0% v/v of nominal solvent','名义溶剂体积的 80.0%')],
        ['HFIP',choose('unknown','未知'),choose('unknown','未知'),'1.200',choose('20.0% v/v of nominal solvent','名义溶剂体积的 20.0%')],
        [choose('Coupling partner','偶联试剂'),choose('unknown','未知'),choose('unknown','未知'),choose('unknown','未知'),choose('Identity and equivalents unresolved','实际选用身份及当量未确定')]])
    partners=rows('results/audit/partner_scenarios.csv')
    t['PARTNERS']=table([choose('Scenario reagent','情景试剂'),choose('Assumed equivalents','假定当量'),'mmol','mg'],[[choose(r['identity'],{'thioanisole':'苯甲硫醚','thiophenol':'苯硫酚'}[r['identity']]),f'{float(r["assumed_equivalents"]):.2f}',f'{float(r["mmol"]):.3f}',f'{float(r["mass_mg"]):.3f}'] for r in partners])
    t['FE']=table([choose('Assumed yield / %','假定产率 / %'),choose('Product / mmol','产物 / mmol'),'FE / %'],[[f'{float(r["assumed_yield_pct"]):.1f}',f'{float(r["product_mmol"]):.3f}',f'{float(r["FE_pct"]):.4f}'] for r in rows('results/audit/charge_balance_scenarios.csv')])
    state_rows=[]
    for name,info in dft['molecules'].items():
        for state,q,m in [('neutral',0,1),('cation',1,2)]:state_rows.append([choose(name,{'target_indole':'目标吲哚','indoline_control':'吲哚啉对照'}[name]),info['formula'],choose(state,{'neutral':'中性态','cation':'阳离子'}[state]),q,m,info['neutral_electrons']-q,info['atoms']])
    t['STATES']=table([choose('Molecule','分子'),choose('Formula','分子式'),choose('State','状态'),choose('Charge','电荷'),choose('Multiplicity','多重度'),choose('Electrons','电子数'),choose('Atoms','原子数')],state_rows)
    t['FEEDBACK']=table([choose('Run','记录'),'j / mA cm⁻²','c / M','T / °C',choose('Given prediction / %','给定预测 / %'),choose('Mock label / %','模拟标签 / %')],[[r['run_id'],f'{float(r["current_density_mA_cm2"]):.1f}',f'{float(r["electrolyte_concentration_M"]):.3f}',f'{float(r["temperature_K"])-273.15:.1f}',f'{float(r["prior_prediction_pct"]):.1f}',f'{float(r["yield_pct"]):.1f}'] for r in rows('data/mock_feedback.csv')])
    names={'original':choose('Original GP','原始 GP'),'scaled_mle':choose('Scaled MLE GP','标准化 MLE GP'),'scaled_fixed':choose('Fixed-kernel GP','固定核 GP'),'mean_baseline':choose('Training-mean baseline','训练均值基线')}
    t['LOO']=table([choose('Method','方法'),'MAE / pp','RMSE / pp'],[[names[k],f'{v["MAE_pp"]:.4f}',f'{v["RMSE_pp"]:.4f}'] for k,v in al['LOO'].items()])
    rec=[[choose('Source UCB (observation SD)','原 UCB（观测标准差）'),f'{orig["current_density_mA_cm2"]:.2f}',f'{orig["electrolyte_concentration_M"]:.3f}',f'{orig["temperature_C"]:.1f}',f'{orig["model_predicted_yield_pct"]:.4f}',f'{orig["prediction_uncertainty_sigma"]:.4f}']]
    for method,policy,label in [('scaled_mle','EI_pp',choose('Scaled MLE EI','标准化 MLE EI')),('scaled_fixed','EI_pp',choose('Fixed-kernel EI','固定核 EI')),('scaled_fixed','UCB_latent_pct',choose('Fixed-kernel UCB','固定核 UCB'))]:
        r=al['selections'][method][policy][0];rec.append([label,f'{r["current_density_mA_cm2"]:.2f}',f'{r["electrolyte_concentration_M"]:.3f}',f'{r["temperature_C"]:.1f}',f'{r["mu_pct"]:.4f}',f'{r["latent_sigma_pp"]:.4f}'])
    t['RECOMMENDATIONS']=table([choose('Policy','策略'),'j / mA cm⁻²','c / M','T / °C','μ / %',choose('σ / pp (latent unless stated)','σ / pp（默认潜在函数）')],rec)
    gantt=[('1','Supervisor and student','Define product, partner and an explicit chemical question','导师与学生','定义产物、偶联试剂和明确化学问题'),('2–3','Trained mentor and student','Instrument training; assay calibration; confirm reagent and cell compatibility','经过训练的指导人员与学生','仪器培训、分析校准、核对物料及池体相容性'),('2–4','Computational mentor','Run one paired DFT pilot; inspect convergence/frequencies/spin','计算指导人员','执行一组配对 DFT 试算，检查收敛、频率与自旋'),('4–5','Supervised bench team','Baseline and independent repeats; archive full raw data','受监督实验人员','基线及独立重复，完整归档原始数据'),('5','Student with mentor review','Validate linked feedback, LOO and baseline comparison','学生及指导人员','核验反馈关联、留一误差与基线比较'),('6–7','Bench and analysis team','Predeclare next candidates; repeat comparator; evaluate held-out outcomes','实验与分析人员','预先声明候选，重复对照，评估留出结果'),('7–8','Student and supervisor','Evidence-limited undergraduate proposal and versioned release','学生与导师','形成证据边界清晰的大创申请和版本化材料')]
    t['GANTT']=table([choose('Proposed week','拟议周次'),choose('Responsible role','责任角色'),choose('Output required before advancing','进入下一步前所需产出')],[[r[0],r[3],r[4]] if zh else r[:3] for r in gantt])
    return t

def parse(text):
    result=[]
    for block in re.split(r'(?=^## \d+\.)',text,flags=re.M)[1:]:
        title,body=block.split('\n',1);parts=[]
        for sub in re.split(r'(?=^### \d+\.)',body,flags=re.M):
            if sub.strip():parts.append(sub.strip())
        result.append((title,parts))
    return result

def main():
    reports={}
    for lang,zh in [('english',False),('chinese',True)]:
        text=(ROOT/f'reports/{lang}.template.md').read_text(encoding='utf-8')
        for key,value in tokens(zh).items():text=text.replace('@@'+key+'@@',value)
        if '@@' in text:raise ValueError('Unresolved report token')
        (ROOT/f'reports/closed_loop_experimental_report_{lang}.md').write_text(text,encoding='utf-8',newline='\n');reports[lang]=text
    both=['# Closed-Loop AI4S Report: Integrating Data-Driven Optimization with Wet-Lab Execution & Quantum Mechanical Modeling','# 闭环 AI4S 报告：数据驱动优化与湿实验操作、量子化学计算的全流程集成','', '2026-09-28 · Complete section-aligned bilingual edition / 完整章节对齐双语版','']
    for (et,ep),(ct,cp) in zip(parse(reports['english']),parse(reports['chinese'])):
        both.append(et+' / '+re.sub(r'^## \d+\. ','',ct));both.append('')
        if len(ep)!=len(cp):raise ValueError('Unmatched sections')
        for e,c in zip(ep,cp):
            eh,eb=e.split('\n',1);ch,cb=c.split('\n',1)
            both.extend([eh+' / '+re.sub(r'^### \d+\.\d+ ','',ch),'',eb.strip(),'',cb.strip(),''])
    (ROOT/'reports/closed_loop_experimental_report_bilingual.md').write_text('\n'.join(both),encoding='utf-8',newline='\n')
    print('Built English, Chinese and complete section-aligned bilingual reports.')
if __name__=='__main__':main()
