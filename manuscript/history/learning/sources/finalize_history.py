"""Deterministic publication metadata and bilingual-table checks; no science jobs."""
from pathlib import Path
import hashlib,json,re,struct
ROOT=Path(__file__).resolve().parent.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(f):return json.loads((ROOT/f).read_text(encoding='utf8'))
def write(f,x):(ROOT/f).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf8')

caps={
1:('Conformer','Historical conformer energy ranges. Every valid structure contributes 50 accepted optimization records; the hatched M07 bar uses UFF and all others use MMFF94. The ranges do not measure entropy or establish unique minima.','历史构象能量范围。每个有效结构有 50 条验收优化记录；阴影线 M07 使用 UFF，其余使用 MMFF94。能量范围不是熵，也不能建立不同极小点的数量。'),
2:('VMC','VMC errors against the saved computed FCI CBS reference. Error bars are ±1 block standard error from 20 saved production blocks per system; they exclude reference uncertainty. The 1.6 mEh line is the unchanged point-estimate gate. Distances are in bohr.','VMC 相对于保存的 FCI CBS 计算参考的误差。误差棒为每体系 20 个保存生产块计算的 ±1 块标准误，未包含参考误差。1.6 mEh 虚线是未改变的点估计验收门槛。键长单位为 bohr。'),
3:('Trace_Windows','Endpoint-change sensitivity over all 72 saved window specifications, using the same three published traces. Voltage panels share the −5 to 45 mV color scale; the source-current panel has an independent −0.16 to 0.05 scale with unresolved current normalization. The black divider marks startup exclusion of at least 10 h. Values are descriptive contrasts, not independent-replicate confidence intervals [@yang2023; @yang2026].','同样三条已发表轨迹在全部 72 组保存窗口设置下的端点差值。电位图共用 −5 至 45 mV 色标；源电流图采用独立的 −0.16 至 0.05 色标，其电流归一化未解决。黑色分界线标记启动排除至少 10 h。数值是描述性差值，不是独立重复实验的置信区间 [@yang2023; @yang2026]。'),
4:('Reference_Alignment','Fixed-geometry gas-phase comparison with PBE0/def2-SVP. Left: raw charge-state offsets. Right: double differences after subtracting each method’s Q01 value. The different vertical scales are intentional. PBE0 is a method comparator, and the raw offsets are not calibrated physical prediction errors.','固定几何、气相条件下与 PBE0/def2-SVP 的比较。左图为原始电荷态偏移；右图为各方法减去自身 Q01 数值之后的双重差值。两图有意采用不同纵轴范围。PBE0 是方法比较对象，原始偏移不是经过校准的物理预测误差。'),
5:('Rank_Reversals','All 15 method/environment comparisons from the 72 state differences. Each off-diagonal cell counts reversed orderings among the same 66 pairs of 12 molecules, not new calculations or errors against experiment. Diagonal zeros are identity comparisons.','由 72 个态差值形成的全部 15 组方法/环境比较。每个非对角格统计相同 12 个分子的 66 个配对中的排序反转数，不是新计算数或相对实验的错误数。对角线零值来自自身比较。')}

for lang in ['english','chinese']:
 p=ROOT/f'section_{lang}.md';s=p.read_text(encoding='utf8');zh=lang=='chinese'
 s=re.sub(r'\[Sorkun, Khetan and Er, Scientific Data \(2019\)\]\(https://www.nature.com/articles/s41597-019-0151-1\)', '[@aqsoldb]',s)
 s=re.sub(r'\[Sorkun、Khetan 与 Er，Scientific Data（2019）\]\(https://www.nature.com/articles/s41597-019-0151-1\)', '[@aqsoldb]',s)
 s=s.replace('[Yang et al., Nature (2023)](https://www.nature.com/articles/s41586-023-05886-z); [Yang et al., Nature Synthesis (2026)](https://www.nature.com/articles/s44160-026-01039-y).','[@yang2023; @yang2026].')
 s=s.replace('[Yang 等，Nature（2023）](https://www.nature.com/articles/s41586-023-05886-z)；[Yang 等，Nature Synthesis（2026）](https://www.nature.com/articles/s44160-026-01039-y)。','[@yang2023; @yang2026]。')
 if '<!-- historical figures -->' not in s:
  block='\n\n<!-- historical figures -->\n'
  for n,(stem,en,cn) in caps.items():
   name=f'Figure_L{n}_{stem}_{lang}'
   label=f'图 L{n}' if zh else f'Figure L{n}'
   block+=f'\n![{label}](figures/{name}.png)\n\n{label}. {cn if zh else en} [SVG](figures/{name}.svg).\n'
  s+=block
 p.write_text(s,encoding='utf8')

phases=[
 ('conformer screening','executed_forcefield','Scoped checks pass; one original parse failure and one UFF substitution; untrained graph smoke is not predictive accuracy.'),
 ('torsion and dynamics','executed_forcefield','Scoped checks pass with force-field/charge substitutions.'),
 ('KRAS complex','incomplete','Complete MD and MM-GBSA acceptance pending; incomplete historical frames not production evidence.'),
 ('reaction path','partial_quantum','Independent saddle refinement passes stationarity and one imaginary-mode checks; original NEB not converged and IRC connectivity not verified.'),
 ('reaction network','failed_acceptance','Multiple imaginary modes invalidate TS; near-zero yield retained.'),
 ('explicit-solvent enhanced sampling','incomplete_sampling','Trajectory completed but product undersampled; no accepted free-energy barrier.'),
 ('strong correlation','incomplete','Unified-basis complete-point review pending; fallback results separated.'),
 ('photochemistry','incomplete','Full QC, dynamics and NTO acceptance pending.'),
 ('self-driving protocol','simulation','Protocol simulation, not physical robot execution.'),
 ('flow twin','model','Control-model execution, not validated hardware.'),
 ('neural wavefunction','executed_vmc','Five-system fresh run complete; He point gate fails; actual trained-network cusp unverified; weights absent.'),
 ('law discovery','incomplete','Checkpoint repair does not establish completion of full training.'),
 ('metalloenzyme PCET','incomplete','Full scientific computation pending acceptance.'),
 ('condensates','historical_snapshot','Historical running state is not a live process or complete scientific acceptance.'),
 ('spin/allostery','incomplete','Full production and window statistics pending acceptance.'),
 ('NPC transport','historical_snapshot','Historical running snapshot; module acceptance pending.'),
 ('relativistic quantum','bounded_model','Specified atomic/finite-basis checks; not a full actinide complex validation.'),
 ('cell metabolism','numerical_model','Specified numerical checks, not experimental whole-cell validation.'),
 ('enzyme geometry','bounded_model','Geometry tests and repairs, not complete enzyme evolution validation.'),
 ('CISS transport','effective_model','Effective transport with assumed OER and analysis responses.'),
 ('cavity QED','effective_model','Harmonic Pauli-Fierz model; optical splitting does not prove catalysis.'),
 ('molecular spin qubits','effective_model','Specified spin model and synthetic bath, not measurement.'),
 ('molecular reservoir','failed_target','Effective circuit model; original performance target not reached.'),
 ('HTE learning','synthetic_response','312 simulated measurements, zero new chemical experiments or hardware execution.'),
 ('asymmetric catalysis','partial_quantum','Preliminary minima/isolated-imine DFT; selectivity undetermined, no accepted DFT TS/IRC.'),
 ('constant-potential interface','design','48 conditions and 3 dry slabs, zero complete solvated interfaces.'),
 ('Ni photodynamics','design','Literature coordinates and protocol; zero new nonadiabatic trajectories.'),
 ('organoelectrocatalysis','published_data_and_scenarios','Real published-data reanalysis kept separate from synthetic lifetime/mechanism scenarios.'),
 ('pyridine site control','executed_quantum_and_scenarios','Actual molecular xTB/DFT; warnings and reference limits retained; no demonstrated same-pair switch.'),
 ('CO2RR twin','uncalibrated_scenario','Numerical conservation checks pass; no DFT/AIMD, 3D CFD, operando or industrial validation.'),
 ('macrocycle pipeline','bounded_forcefield','18 MMFF convergences; QM/MM and FEP blocked; zero accepted candidates, zero actual provider tokens.')]
ledger=[]
for n,(topic,typ,boundary) in enumerate(phases,1):
 source=('complex_PHASE1_19_CLOSEOUT.md' if n<=19 else 'complex_PHASE_INDEX.md')
 if n==11:source='phase11_acceptance.json'
 elif n in [25,26,27]:source='phase25_27_acceptance.json'
 elif n==24:source='phase24_metrics.json'
 elif n in [29,30]:source=f'phase{n}_acceptance.json'
 elif n==31:source='phase31_status.json'
 text=(ROOT/'sources'/source).read_text(encoding='utf8')
 needles=[f'Phase{n}',f'Phase {n}',f'phase{n:02d}',f'phase{n}']
 lines=[i for i,l in enumerate(text.splitlines(),1) if any(x in l for x in needles)]
 ledger.append({'phase':n,'topic':topic,'evidence_type':typ,'boundary':boundary,
   'source_artifact_id':source,'matching_text_lines':lines[:6],
   'inspection_scope':'Latest scoped status and selected outputs; not exhaustive all historical file validation'})
write('sources/phase_coverage.json',ledger)

audit=read('audit.json');summary=read('sources/posthoc_summary.json')
audit['coverage_ledger']='sources/phase_coverage.json'
audit['principal_claims']=[
 {'id':'shared_ancestry','evidence':'exact_git_object_and_byte_hash','source':'sources/deduplication.json','pointer':'/common_root_commit','value':summary['deduplication']['shared_ancestral_commit'],'boundary':'Two repositories are not independent replications.'},
 {'id':'reported_solubility_RMSE','evidence':'report_only_not_recomputed','source':'sources/aqueous_README.md','text_selector':'0.635','value':0.635,'boundary':'No frozen predictions, labels, split IDs or weights in the committed tree.'},
 {'id':'random_not_scaffold','evidence':'source_protocol','source':'sources/aqueous_aqsol_model.py','text_selector':'KFold(n_splits=5, shuffle=True, random_state=42)','boundary':'No historical scaffold-split result; canonical identity exclusion does not establish analogue/source independence.'},
 {'id':'conformer_count','evidence':'executed_forcefield_recomputed','source':'sources/posthoc_summary.json','pointer':'/phase01/accepted_optimization_records','value':500,'boundary':'Records, not unique minima or QM calls; same 10 structures in original and audited output.'},
 {'id':'VMC_He_accuracy','evidence':'executed_vmc_recomputed','source':'sources/phase11_results.json','pointer':'/systems/He/d_fci','value':1.6608865319498456,'unit':'mEh','boundary':'Fails fixed 1.6 mEh point gate; computed CBS reference has uncertainty.'},
 {'id':'VMC_internal_check','evidence':'saved_selftest_same_run','source':'sources/phase11_results.json','pointer':'/self_tests/laplacian_fd_relerr','value':3.515304116655698e-7,'boundary':'Not variational accuracy or a trained-wavefunction cusp certificate.'},
 {'id':'trace_record_count','evidence':'published_aggregate_recomputed','source':'sources/posthoc_summary.json','pointer':'/phase28/published_trace_rows','value':196496,'boundary':'Only 39726 distinct within-trace timestamps; electrode count unknown.'},
 {'id':'trace_window_dependence','evidence':'same_trace_estimand_sensitivity','source':'sources/derived_trace_summary.csv','selector':{'trace_id':'nature2023_constant_current'},'columns':['timestamp_balanced_1h_change','startup_ge10h_min','startup_ge10h_max'],'values':[0.0442499999999999,-0.0040000000000000036,0.00924999999999998],'unit':'V','boundary':'Specification sensitivity, not between-electrode uncertainty.'},
 {'id':'xtb_count_dedup','evidence':'job_ledger_recomputed','source':'sources/posthoc_summary.json','pointer':'/phase29','boundary':'156 executions=12 initial+144 extension; 144 unique settings; all extension calls retain warnings.'},
 {'id':'DFT_calls','evidence':'individual_result_and_ledger','source':'sources/phase29_dft_manifest.json','pointer':'/jobs','count':8,'boundary':'Four charge-state pairs, not eight independent molecules.'},
 {'id':'reference_centered_Q03_GFN1','evidence':'matched_geometry_recomputed','source':'sources/derived_reference_alignment.csv','selector':{'molecule_id':'Q03','gfn':'1'},'column':'double_difference_eV','value':0.25673639503660706,'unit':'eV','boundary':'PBE0 comparator; no physical redox calibration.'},
 {'id':'conservation_without_validation','evidence':'numerical_scenario','source':'sources/phase30_acceptance.json','pointer':'/checks/carbon_balance_mol_s','value':-4.1165801236558996e-19,'boundary':'Uncalibrated scenario despite conservation.'}]
audit['paired_counterexamples']=[
 {'identity':'same frozen Phase11 seed11 run','passed':'antisymmetry and Laplacian check','failed_or_missing':'He 1.6mEh gate; trained-wavefunction cusp validation','pairing_limit':'Run-level checks, not every trained-weight state revalidated; weights absent.'},
 {'identity':'same Phase29 geometry_sha256 within molecular method pairs','passed':'SCF convergence and stored arithmetic','failed_or_missing':'Physical electron-reference alignment and external chemistry accuracy','pairing_limit':'PBE0 is not truth; three molecules and one-parent basis test only.'},
 {'identity':'same three published traces','passed':'preserved data identity and numerical aggregation','failed_or_missing':'Window-independent degradation estimate and independent lifetime observations','pairing_limit':'No independently replicated electrodes inferred from rows.'}]
audit['proposed_figures']=[{'number':f'L{i}','files':[f'figures/Figure_L{i}_{v[0]}_{lang}.png' for lang in ['english','chinese']],
 'purpose':v[1]} for i,v in caps.items()]

# Every table's numeric tokens must be equal across the two languages, including units' exponents.
texts=[(ROOT/f'section_{l}.md').read_text(encoding='utf8') for l in ['english','chinese']]
def tables(s):
 out=[];cur=[]
 for line in s.splitlines()+['']:
  if line.startswith('|'):cur.append(line)
  elif cur:out.append(cur);cur=[]
 return out
en,cn=map(tables,texts);assert len(en)==len(cn)==6
num=lambda s:re.findall(r'[-−]?\d+(?:,\d{3})*(?:\.\d+)?|[⁰¹²³⁴⁵⁶⁷⁸⁹]+',s)
for i,(e,c) in enumerate(zip(en,cn),1):
 assert len(e)==len(c),(i,'row count')
 for j,(a,b) in enumerate(zip(e[2:],c[2:]),1):assert num(a)==num(b),(i,j,num(a),num(b))
for s in texts:
 assert '\\(' not in s and '\\)' not in s
 for link in re.findall(r'\]\(([^)]+)\)',s):
  if not link.startswith('http'):assert (ROOT/link).is_file(),link
fm=read('figures/figure_manifest.json')
for f in fm['figures']:
 p=ROOT/f['path'];assert sha(p)==f['sha256']
 if p.suffix=='.png':assert struct.unpack('>II',p.read_bytes()[16:24])==(2700,1260)
fm['visual_qa']={'status':'passed_sampled_render_review','reviewed_pngs':[
 'Figure_L1_Conformer_chinese.png','Figure_L2_VMC_english.png','Figure_L3_Trace_Windows_chinese.png',
 'Figure_L4_Reference_Alignment_english.png','Figure_L5_Rank_Reversals_chinese.png'],
 'scope':'One rendered language per figure inspected for units, clipping, labels and scientific meaning. All20 files hash/dimension checked; bilingual numeric tables share inputs. No claim of full PDF review.'}
write('figures/figure_manifest.json',fm)
audit['publication_validation']={'status':'passed','tables_each_language':6,'table_numeric_parity':True,
 'english_word_count':len(re.findall(r"\b[A-Za-z][A-Za-z'’-]*\b",texts[0])),
 'chinese_character_count':len(re.findall('[\u4e00-\u9fff]',texts[1])),
 'figures_each_language':5,'png_files':10,'svg_files':10,'new_science_jobs':0,
 'DFT_difference_max_error_eV':summary['phase29']['DFT_energy_arithmetic_max_error_eV'],
 'xTB_difference_max_error_eV':summary['phase29']['matrix_energy_arithmetic_max_error_eV'],
 'history_categories':'All54 commit subjects inventoried; latest scoped status covers all31 phases.',
 'visual_qa':'figures/figure_manifest.json'}
audit['integrity_hashes']={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in sorted(ROOT.rglob('*'))
 if p.is_file() and p.name!='audit.json' and '__pycache__' not in str(p)}
write('audit.json',audit)
print(json.dumps(audit['publication_validation']))
