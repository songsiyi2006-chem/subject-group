"""Deterministic publication validation; does not claim chemical validity."""
import csv,hashlib,json,re,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
    checks=[]
    def check(name,ok):checks.append(dict(check=name,passed=bool(ok)))
    def j(path):return json.loads((ROOT/path).read_text(encoding='utf-8-sig'))
    def rows(path):
        with (ROOT/path).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
    source=j('source/source_record.json');run=j('results/original/execution.json');data=j('results/original/research_grade_results.json');summary=j('results/audit/audit_summary.json')
    digest=hashlib.sha256((ROOT/'source/run_research_engine_original.py').read_bytes()).hexdigest()
    check('source hash and execution hash agree',digest==source['script_sha256']==run['script_sha256'])
    check('unchanged original completed all three modules',run['returncode']==0 and run['result_written'] and len(data)==3 and all(v['status']=='completed' for v in data.values()))
    for filename,n in [('conditions_source_rounded.csv',135),('conditions_unrounded.csv',135),('pareto_complete.csv',12),('scalarization_weights.csv',231),('gp_predictions.csv',270),('scope_descriptors.csv',8),('scope_conformers.csv',64),('scope_xtb.csv',8),('scope_random_yield_seeds.csv',8000),('transport_audit.csv',10),('transport_grid_convergence.csv',40),('external_film_scenarios.csv',50),('transport_parameter_sensitivity.csv',60)]:check('rows '+filename,len(rows('results/audit/'+filename))==n)
    for r in summary['scope']['rows']:
        for state in ['neutral','cation']:check('completed xTB '+r['molecule']+'/'+state,j('results/molecules/'+r['molecule']+'/'+state+'/run.json')['completed'])
    contents={}
    for lang in ['english','chinese','bilingual']:
        p=ROOT/'reports'/f'research_grade_technical_report_{lang}.md';text=p.read_text(encoding='utf-8');contents[lang]=text
        check(lang+' exactly five sections',re.findall(r'^## (\d+)\.',text,re.M)==['1','2','3','4','5'])
        for section,count in [(1,2),(2,4),(3,4),(4,4),(5,5)]:check(lang+f' section {section} subsections',all(re.search(rf'^### {section}\.{i} ',text,re.M) for i in range(1,count+1)))
        check(lang+' source and audit numerical anchors',all(v in text for v in ['43.47','74.49','0.779','4156.7141206675215','129.379256','260.666929','0.733333','0.916667','0.00029997','−0.02830746']))
        check(lang+' evidence labels',all('['+v+']' in text for v in ['L','S','M','Q','A','P','U']))
        check(lang+' source identity',digest in text)
        check(lang+' balanced fences',text.count('```')%2==0);check(lang+' balanced display math',text.count('$$')%2==0)
        check(lang+' no unresolved template',not re.search(r'@@[A-Z_]+@@',text))
        for link in re.findall(r'\]\(([^)]+)\)',text):
            if not link.startswith(('http://','https://','#')):check(lang+' link '+link,(p.parent/link).is_file())
    def decimals(text):return [re.findall(r'[-−]?\d+\.\d+',line) for line in text.splitlines() if line.startswith('|')]
    check('English/Chinese table decimal parity',decimals(contents['english'])==decimals(contents['chinese']))
    check('English no Chinese prose',not re.search(r'[\u4e00-\u9fff]',contents['english']))
    for lang in ['english','chinese']:
        bodies=re.split(r'^### [^\n]+\n',contents[lang],flags=re.M)[1:]
        bodies=[re.split(r'^## ',b,flags=re.M)[0].strip() for b in bodies]
        check(lang+' complete subsection bodies in combined edition',all(b in contents['bilingual'] for b in bodies))
        intro=re.split(r'^### ',re.split(r'^## 1\. [^\n]+\n',contents[lang],flags=re.M)[1],flags=re.M)[0].strip()
        check(lang+' abstract preserved in combined edition',intro in contents['bilingual'])
    for p in sorted((ROOT/'results').rglob('*.json')):
        if p.name=='publication_validation.json':continue
        try:json.loads(p.read_text(encoding='utf-8-sig'));ok=True
        except (ValueError,OSError):ok=False
        check('JSON parses '+p.relative_to(ROOT).as_posix(),ok)
    for p in [ROOT/'README.md',ROOT.parent/'README.md']:
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
            if not link.startswith(('http://','https://','#')) and not link.endswith('publication_validation.json'):check('README link '+link,(p.parent/link).exists())
    result=dict(passed=all(v['passed'] for v in checks),checks=len(checks),scope='Saved record and publication integrity, not chemical prediction accuracy or experimental validation.',failed=[v for v in checks if not v['passed']],items=checks)
    (ROOT/'results/publication_validation.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:result[k] for k in ['passed','checks','failed']},indent=2));return 0 if result['passed'] else 1
if __name__=='__main__':sys.exit(main())
