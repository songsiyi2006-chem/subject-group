"""Release checks for documents, source provenance and figures, not assay certification."""
from pathlib import Path
import csv,hashlib,json,re,sys,xml.etree.ElementTree as ET
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]

def main():
    checks=[]
    def check(name,ok):checks.append(dict(check=name,passed=bool(ok)))
    def j(name):return json.loads((ROOT/name).read_text(encoding='utf-8-sig'))
    def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
    def rows(name):
        with (ROOT/name).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
    source=j('source/source_record.json');run=j('results/original/execution.json');audit=j('results/audit/audit_summary.json');fig=j('results/figure_manifest.json');qa=j('results/figure_qa.json')
    check('exact archived specification',digest(ROOT/'source/specification.md')==source['specification_sha256'])
    check('exact complete source extraction and executed script',digest(ROOT/'source/run_deployment_and_figures.py')==source['script_sha256']==run['script_sha256'])
    check('failed naive regex retained',source['naive_regex_extraction']['parses'] is False)
    check('original executed but no measurement or DFT or submission',run['returncode']==0 and run['result_written'] and not any(run[k] for k in ['raw_instrument_data_supplied','DFT_executed','grant_submitted']))
    an=audit['analytical'];check('1292 percent failure retained',an['HPLC']['yield_pct']==1292 and not an['HPLC']['numeric_range_pass'] and not an['concordance_validated'])
    check('mass-based NMR value retained',abs(an['NMR_mass_recomputed']['yield_pct']-91.32528687793567)<1e-10)
    check('MC assumption status',an['uncertainty_scenario']['draws']==100000 and an['uncertainty_scenario']['HPLC_fraction_above_100']==1)
    for name,n in [('qnmr_scenarios.csv',12),('relaxation_scenarios.csv',6),('concordance_counterfactuals.csv',3),('synthetic_calibration.csv',21),('synthetic_chromatogram.csv',1401),('deconvolution_recovery.csv',80),('pareto_conditions.csv',80),('scope_mock.csv',6),('scope_leave_one_out_descriptive.csv',6),('hypothetical_energy.csv',5),('grant_budget.csv',7)]:check('row count '+name,len(rows('results/audit/'+name))==n)
    for name in ['synthetic_calibration.csv','synthetic_chromatogram.csv','deconvolution_recovery.csv','pareto_conditions.csv']:check('synthetic roles '+name,all(r['evidence_role']=='synthetic' for r in rows('results/audit/'+name)))
    check('scope roles not experiments',all(r['evidence_role']=='mock_not_experimental' for r in rows('results/audit/scope_mock.csv')))
    check('energy roles not DFT',all(r['evidence_role']=='illustrative_not_DFT' for r in rows('results/audit/hypothetical_energy.csv')))
    check('Pareto correction retained',audit['figures']['Pareto']['nondominated_count']==23 and audit['figures']['Pareto']['source_threshold_count']==0 and audit['figures']['Pareto']['assumed_noise_draws']==2000)
    check('width mismatch negative result retained',audit['synthetic_recovery']['summaries'][1]['mean_absolute_area_error_pct']>14)
    check('unapproved planning budget',audit['grant']['total_CNY']==10000 and not audit['grant']['approval_or_grant_submission'] and audit['grant']['confirmed_funding_CNY'] is None)
    check('16 figures in two formats/languages',len(fig['figures'])==16 and len({r['file'] for r in fig['figures']})==16)
    for r in fig['figures']:
        path=ROOT/'reports/figures'/r['file'];check('figure hash '+r['file'],digest(path)==r['sha256'])
        for name,sha in r['sources'].items():check(r['file']+' data '+name,digest(ROOT/'results/audit'/name)==sha)
        if r['format']=='png':
            with Image.open(path) as im:
                check('PNG size and resolution '+r['file'],list(im.size)==r['pixels'] and im.width>=2100 and all(abs(v-300)<.1 for v in im.info.get('dpi',[0,0])))
        else:
            tree=ET.parse(path);ns={'s':'http://www.w3.org/2000/svg'};texts=tree.findall('.//s:text',ns);text=' '.join(''.join(t.itertext()) for t in texts)
            check('editable pure vector '+r['file'],len(texts)>10 and not tree.findall('.//s:image',ns) and r['editable_text'] and not r['embedded_raster'])
            check('visible evidence label '+r['file'],any(word in text for word in ['SYNTHETIC','MOCK','HYPOTHETICAL','模拟','假设','合成分析']))
    check('visual QA covers current file bytes',qa['reviewed_figure_hashes']=={r['file']:r['sha256'] for r in fig['figures']} and qa['png_files_visually_reviewed']==8 and qa['svg_files_rendered_and_visually_reviewed']==8)
    contents={}
    for lang in ['english','chinese']:
        path=ROOT/'reports'/f'deployment_toolkit_report_{lang}.md';text=path.read_text(encoding='utf-8');contents[lang]=text
        check(lang+' five sections',re.findall(r'^## (\d+)\.',text,re.M)==['1','2','3','4','5'])
        for sec,n in [(1,3),(2,5),(3,4),(4,3),(5,3)]:check(lang+f' full subsections {sec}',all(re.search(rf'^### {sec}\.{i} ',text,re.M) for i in range(1,n+1)))
        check(lang+' numeric anchors',all(v in text for v in ['1292.00','91.392','91.3252869','1200.608','0.027363','14.614337','0.982099','0.7691904733']))
        check(lang+' evidence labels',all('['+x+']' in text for x in ['L','S','A','M','P','U']))
        check(lang+' source hash',source['script_sha256'] in text)
        check(lang+' complete figure captions',all('**'+('Figure' if lang=='english' else '图')+' '+str(i) in text for i in range(1,5)))
    def nums(t):return [re.findall(r'[-−]?\d+\.\d+',line) for line in t.splitlines() if line.startswith('|')]
    check('English Chinese report table numeric parity',nums(contents['english'])==nums(contents['chinese']))
    check('English report no Chinese prose',not re.search(r'[\u4e00-\u9fff]',contents['english']))
    proposals={}
    for lang,suffix in [('chinese',''),('english','_English')]:
        p=ROOT/'proposals'/f'National_Undergraduate_Grant_Proposal_PanTang_Lab{suffix}.md';t=p.read_text(encoding='utf-8');proposals[lang]=t
        check(lang+' complete twelve-section proposal',len(re.findall(r'^## ',t,re.M))==12)
        check(lang+' proposal planned counts',all(s in t for s in ['40','24','120','10000','1292.00','14.61']))
    check('English proposal no Chinese prose',not re.search(r'[\u4e00-\u9fff]',proposals['english']))
    check('English Chinese proposal numeric table parity',nums(proposals['english'])==nums(proposals['chinese']))
    reviewed=[ROOT/'README.md',ROOT/'data/README.md',ROOT/'reports/figures/README.md',*list((ROOT/'reports').glob('*.md')),*list((ROOT/'proposals').glob('*.md'))]
    for p in reviewed:
        t=p.read_text(encoding='utf-8');check(p.name+' balanced code fences',t.count('```')%2==0);check(p.name+' no unexpanded template','@@' not in t)
        for link in re.findall(r'\]\(([^)]+)\)',t):
            if not link.startswith(('https://','http://','#')):check(p.name+' link '+link,(p.parent/link).exists() or (p.parent/link).resolve()==(ROOT/'results/publication_validation.json').resolve())
    for p in ROOT.rglob('*.json'):
        if p.name=='publication_validation.json':continue
        try:json.loads(p.read_text(encoding='utf-8-sig'));ok=True
        except (ValueError,OSError):ok=False
        check('JSON parses '+p.relative_to(ROOT).as_posix(),ok)
    result=dict(passed=all(c['passed'] for c in checks),checks=len(checks),failed=[c for c in checks if not c['passed']],scope='Record, source, publication and numerical integrity; not analytical, chemical or grant validation',items=checks)
    (ROOT/'results/publication_validation.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:result[k] for k in ['passed','checks','failed']},indent=2));return 0 if result['passed'] else 1

if __name__=='__main__':sys.exit(main())
