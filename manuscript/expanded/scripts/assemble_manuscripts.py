"""Assemble two main manuscripts, preserving source fragments and v1 publication."""
import hashlib,json,re,os,subprocess
from pathlib import Path
BASE=Path(__file__).resolve().parents[1];OLD=BASE.parent;ROOT=OLD.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rebase(text,origin):
    def sub(m):
        url=m.group(2)
        if re.match(r'\w+://',url) or url.startswith('#'):return m.group(0)
        target=(origin/url.strip('<>')).resolve()
        return m.group(1)+'('+Path(os.path.relpath(target,BASE)).as_posix()+')'
    return re.sub(r'(!?\[[^\]\n]*\])\(([^)\n]+)\)',sub,text)
def fragment(path):
    return rebase(path.read_text('utf-8-sig').strip(),path.parent)
def numbered_chapter(text,number,title):
    lines=text.splitlines()
    if lines and lines[0].startswith('#'):
        lines=lines[1:]
    body='\n'.join(lines).strip()
    idx=0
    def heading(m):
        nonlocal idx
        idx+=1
        h=re.sub(r'^(?:[A-Z]\.)?\d+(?:\.\d+)*\.?\s+','',m.group(1))
        return f'### {number}.{idx} {h}'
    body=re.sub(r'^#{2,4}\s+(.+)$',heading,body,flags=re.M)
    return f'## {number} {title}\n\n'+body
abstracts={
'english':'''Computational chemistry projects combine molecular representations, electronic references, statistical sampling and downstream observables, but success at one stage does not establish accuracy at the next. We integrate selected numerical evidence from six public repository histories into a reproducible retrospective study, deduplicating inherited studies and distinguishing executed calculations from synthetic controls and proposed work. Two controlled H₂ chains provide the central quantitative analysis. Across paired neural fits, gradient supervision reduces RHF-label test RMSE by 60.24–73.65%, but error against matched finite-basis FCI by only 0.00838–0.02168%. Signed attribution identifies reference bias and retains error cancellation. In a separate FCI-parameterized Morse model, near-machine-precision algebraic residuals coexist with a missing shallow bound state; domain extension and refinement recover it, and a post hoc fixed-domain extrapolation reaches a 0.02688% discrepancy from the analytic model. Historical conformer, variational, spin-state, molecular-dynamics and docking records extend the analysis to reference alignment, state quality and incomplete sampling. Transport and analytical controls show why conservation, a correct mean or a small fitting residual can miss the requested product, distribution or peak area. The expanded evidence includes all recorded seed-level comparisons and unsuccessful outcomes, rather than only favorable results. New figures and tables are recalculated from archived records, with fixed commits and source hashes; no new experiments are claimed. The contribution is an auditable, observable-specific interpretation of coupled computation, not a new general-purpose chemistry algorithm or a validated catalyst mechanism.''',
'chinese':'''计算化学项目将分子表示、电子参考、统计采样及下游观测量相结合，但某一阶段的成功不能证明下一阶段的准确性。本文整合六个公开仓库历史中选取的数值证据，去除继承研究的重复计数，并区分已执行计算、合成对照与拟议工作，形成可复现的回顾性研究。两条受控 H₂ 计算链提供核心定量分析。在配对神经拟合中，梯度监督使相对 RHF 标签的测试 RMSE 降低 60.24–73.65%，但相对匹配有限基组 FCI 的误差只降低 0.00838–0.02168%；有符号归因识别参考偏差并保留误差抵消。在独立的 FCI 参数化 Morse 模型中，接近机器精度的代数残差与浅束缚态遗漏并存；扩展区间及加密网格恢复该态，事后固定区间外推相对解析模型的偏差达到 0.02688%。历史构象、变分、自旋态、分子动力学及对接记录将分析拓展到参考对齐、态质量与不充分采样。传输和分析对照说明，守恒、正确均值或较小拟合残差，仍可能遗漏所需产物、分布或峰面积。扩展证据保留全部已记录种子比较和未成功结果，不只选择有利结果。新增图表由归档记录重算，并附固定提交与文件摘要；本文不宣称新增实验。贡献在于对耦合计算给出可审计、针对观测量的解释，而非提出新的通用化学算法或经验证的催化机理。'''}
commits={
'aqueous-solubility-ml-benchmark':'5939547a936ed4baf7fa7ee9fac6fa6e90bfa8fe',
'ai4chem-complex-scaffolds-benchmark':'5e19c519dbd11430926828fddf5375b919af67f0',
'Paper-Analysis-ChemMLLM':'c1d309a0c2f876ccc8b4d2ce753ece971d58fe97',
'ai4pharm-lead-developability-suite':'a7a2a4df501d5dad030dc488691b02c73fdb423f',
'pincer-catmech-ai':'2edffb123791bd61acbfdeff763f603ee50c2287',
'subject-group':'073bb3872d9112ec1ca46b3b9e80476d00ffbf97'}
roles_en=['Inherited solubility code; no row-level prediction archive','Conformers, variational blocks and matched reference calculations','Literature-analysis provenance, not new molecular calculations','Executed sampling and docking controls','State-quality, baseline and pathway audits','Electronic, nuclear, transport, kinetic and analytical controls']
roles_cn=['继承溶解度代码；无逐行预测档案','构象、变分区块与匹配参考计算','文献分析来源，不计新增分子计算','已执行采样与对接对照','态质量、基线与路径审计','电子、核、传输、动力学和分析对照']
inputs={};documents={}
for lang in ['english','chinese']:
    zh=lang=='chinese';p=OLD/f'manuscript_{lang}.md';inputs[str(p.relative_to(ROOT))]=sha(p)
    old=fragment(p)
    intro_pos=re.search(r'^## 1 ',old,re.M).start()
    abstract_heading='## 摘要' if zh else '## Abstract'
    title='分子学习与多尺度计算化学中的参考和观测量依赖性' if zh else 'Reference and observable dependence in molecular learning and multiscale computational chemistry'
    front=f'# {title}\n\n'+('宋思毅\n\n广西师范大学，中国' if zh else 'Siyi Song\n\nGuangxi Normal University, China')+'\n\n'+('扩展研究论文 · 2026 年 9 月 29 日' if zh else 'Expanded research manuscript · 29 September 2026')+'\n\n'+abstract_heading+'\n\n'+abstracts[lang]+'\n\n'+('关键词：计算化学；参考偏差；分子学习；采样；数值收敛；可复现性' if zh else 'Keywords: computational chemistry; reference bias; molecular learning; sampling; numerical convergence; reproducibility')+'\n\n'
    core=old[intro_pos:re.search(r'^### 3\.4 ',old,re.M).start()].strip()
    method_end=re.search(r'^## 3 ',core,re.M).start()
    scopepath=BASE/f'scope_{lang}.md';scope=fragment(scopepath);inputs[str(scopepath.relative_to(ROOT))]=sha(scopepath)
    roles=roles_cn if zh else roles_en
    ht=('**表 H1. 历史取样框及冻结提交。**\n\n| 仓库 | 提交（前 8 位） | 选入证据的作用 |\n|---|---|---|\n' if zh else '**Table H1. Historical sampling frame and frozen commits.**\n\n| Repository | Commit prefix | Role of selected evidence |\n|---|---|---|\n')
    ht+='\n'.join(f'| [{name}](https://github.com/songsiyi2006-chem/{name}/tree/{commit}) | {commit[:8]} | {role} |' for (name,commit),role in zip(commits.items(),roles))
    scope=scope.replace('{{HISTORY_TABLE}}',ht)
    core=core[:method_end]+scope+'\n\n'+core[method_end:]
    # The original source remains frozen; the additional sampling frame is explicit.
    core=core.replace('The source collection is frozen at repository commit','The original source collection is frozen at repository commit')
    parts=[front+core]
    titles={4:('Molecular ensembles and reference alignment','分子系综与参考对齐'),5:('Catalytic state quality and limits of pathway inference','催化态质量与路径推断边界'),6:('Sampling and docking in molecular developability','分子可开发性中的采样与对接'),7:('Multiscale observables and controlled model failures','多尺度观测量与受控模型失败')}
    for number,group in [(4,'learning'),(5,'catalysis'),(6,'pharmacology'),(7,None)]:
        path=OLD/'history'/group/f'section_{lang}.md' if group else BASE/f'multiscale_{lang}.md'
        inputs[str(path.relative_to(ROOT))]=sha(path)
        parts.append(numbered_chapter(fragment(path),number,titles[number][int(zh)]))
    discussion=old[re.search(r'^## 4 ',old,re.M).start():re.search(r'^## 5 ',old,re.M).start()]
    discussion=re.sub(r'^(##(?:#)?) 4(?=[ .])',r'\g<1> 8',discussion,flags=re.M)
    syn=BASE/f'synthesis_{lang}.md';inputs[str(syn.relative_to(ROOT))]=sha(syn)
    parts.append(discussion.strip()+'\n\n'+fragment(syn))
    end=old[re.search(r'^## 5 ',old,re.M).start():]
    end=re.sub(r'^## 5 ', '## 9 ',end,flags=re.M)
    addition=('扩展稿的历史提取、派生数值、制图脚本及固定提交登记位于 manuscript/history 与 manuscript/expanded。所有新增图表均有可追溯的归档输入；公开实测数据的实验归属仍为原始生产者。' if zh else 'Historical extracts, derived values, plotting scripts and fixed-commit records for this expansion are supplied under manuscript/history and manuscript/expanded. All new figures and tables have traceable archived inputs; empirical datasets remain attributable to their original experimental producers.')
    end=end.replace('## 作者与稿件准备声明',addition+'\n\n## 作者与稿件准备声明') if zh else end.replace('## Author and preparation statements',addition+'\n\n## Author and preparation statements')
    parts.append(end)
    text='\n\n'.join(parts).strip()+'\n'
    # Editorial integration of independently retained fragments; source bytes stay intact.
    text=text.replace('A useful manuscript figure should display the whole window grid and mark the startup policy explicitly.','Figure L3 displays the whole window grid and marks the startup policy explicitly.')
    text=text.replace('因此，论文图应展示完整窗口网格，并明确标出启动阶段排除规则，而不是只展示最有利的一项统计。','图 L3 展示完整窗口网格，并明确标出启动阶段排除规则，保留全部统计口径。')
    text=text.replace("Paths are relative to this section's directory.", 'Run these two commands from the repository subdirectory `manuscript/history/catalysis`.')
    text=text.replace('路径均以本节目录为相对起点。','上述两条命令均从仓库子目录 `manuscript/history/catalysis` 执行。')
    text=text.replace(" The accompanying editorial assessment evaluates those gaps separately from the manuscript's demonstrated results.",'')
    text=text.replace('随附投稿评估将这些缺口与本文已展示的结果分别评价。','')
    for doi,key in [('10.1016/j.xcrp.2026.103310','chem_mllm'),('10.1016/j.chempr.2025.102922','cu_data_driven'),('10.1016/j.chempr.2025.102838','single_atom_review')]:
        text=re.sub(r'\[[^\]]+\]\(https://doi.org/'+re.escape(doi)+r'\)', '[@'+key+']',text)
    # Give figures and tables one global sequence, retaining mapping to fragment IDs.
    mapping={}
    for kind,cn in [('Figure','图'),('Table','表')]:
        pattern=re.compile(r'^(?:\*\*)?(?:'+kind+'|'+cn+r')\s*([A-Z]?\d+)\s*[.．。]',re.M)
        ids=[m.group(1) for m in pattern.finditer(text)]
        if len(set(ids))!=len(ids):raise ValueError((lang,kind,'duplicate caption',ids))
        mapping[kind]={key:str(i+1) for i,key in enumerate(ids)}
        if kind=='Figure':
            text=re.sub(r'Figures\s+([A-Z]\d+)[–—-]([A-Z]\d+)',lambda m:'Figures '+mapping[kind][m.group(1)]+'–'+mapping[kind][m.group(2)],text)
            text=re.sub(r'图\s*([A-Z]\d+)[–—-]([A-Z]\d+)',lambda m:'图 '+mapping[kind][m.group(1)]+'—'+mapping[kind][m.group(2)],text)
        text=re.sub(r'\b'+kind+r'\s+([A-Z]?\d+)\b',lambda m:kind+' '+mapping[kind].get(m.group(1),m.group(1)),text)
        text=re.sub(cn+r'\s*([A-Z]?\d+)(?!\d)',lambda m:cn+' '+mapping[kind].get(m.group(1),m.group(1)),text)
    if '{{' in text:raise ValueError('Unresolved table marker')
    target=BASE/f'manuscript_{lang}.md';target.write_bytes(text.encode())
    documents[target.stem]={'source_sha256':sha(target),'figures':len(mapping['Figure']),'tables':len(mapping['Table']),'caption_mapping':mapping,'characters':len(text),'whitespace_words':len(text.split())}
inventory=json.loads((OLD/'history/repository_inventory.json').read_text('utf-8-sig'))
for row in inventory['repositories']:row['frozen_commit']=commits[row['name']]
(OLD/'history/repository_inventory.json').write_bytes((json.dumps(inventory,ensure_ascii=False,indent=2)+'\n').encode())
(BASE/'results/assembly.json').write_bytes((json.dumps({'documents':documents,'input_sha256':inputs,'script_sha256':sha(Path(__file__)),'no_new_scientific_solver_calls':True},ensure_ascii=False,indent=2)+'\n').encode())
print(json.dumps(documents,ensure_ascii=False,indent=2))
