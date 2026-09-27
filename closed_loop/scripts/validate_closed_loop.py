"""Validate saved publication and record integrity; never certify chemistry."""
from pathlib import Path
import csv,hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[1]
def main():
    checks=[]
    def check(name,ok):checks.append(dict(check=name,passed=bool(ok)))
    def j(p):return json.loads((ROOT/p).read_text(encoding='utf-8-sig'))
    def rows(p):
        with (ROOT/p).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
    source=j('source/source_record.json');run=j('results/original/execution.json');audit=j('results/audit/active_learning.json');manifest=j('inputs/reviewed/input_manifest.json')
    check('original extraction and execution hashes',source['script_sha256']==run['script_sha256']==hashlib.sha256((ROOT/'source/run_closed_loop_platform.py').read_bytes()).hexdigest())
    check('source specification bytes preserved',source['specification_sha256']==hashlib.sha256((ROOT/'source/specification.md').read_bytes()).hexdigest())
    check('original completed without DFT or wet-lab execution',run['returncode']==0 and run['result_written'] and not run['dft_executed'] and not run['wet_lab_executed'])
    check('all feedback remains mock',audit['evidence_role']=='mock' and all(r['evidence_role']=='mock' for r in rows('data/mock_feedback.csv')))
    check('feedback provenance hash',audit['feedback_sha256']==hashlib.sha256((ROOT/'data/mock_feedback.csv').read_bytes()).hexdigest())
    for p,n in [('stoichiometry_sweep.csv',81),('partner_scenarios.csv',8),('charge_balance_scenarios.csv',4),('loo_predictions.csv',24),('candidate_predictions.csv',270),('acquisition_sensitivity.csv',9)]:check('rows '+p,len(rows('results/audit/'+p))==n)
    check('corrected input file count',len(manifest['files'])==12)
    for r in manifest['files']:check('input SHA256 '+r['path'],hashlib.sha256((ROOT/'inputs/reviewed'/r['path']).read_bytes()).hexdigest()==r['sha256'])
    contents={}
    for lang in ['english','chinese','bilingual']:
        p=ROOT/'reports'/f'closed_loop_experimental_report_{lang}.md';text=p.read_text(encoding='utf-8');contents[lang]=text
        check(lang+' exactly five main sections',re.findall(r'^## (\d+)\.',text,re.M)==['1','2','3','4','5'])
        for sec,count in [(1,3),(2,5),(3,5),(4,6),(5,4)]:check(lang+f' subsections {sec}',all(re.search(rf'^### {sec}\.{i} ',text,re.M) for i in range(1,count+1)))
        check(lang+' numeric anchors',all(v in text for v in ['42.4535452','37.7364846','44.655','232.458','2.883333','15.2560','13.9800','78.3720','11.6784','1.344386']))
        check(lang+' evidence labels',all('['+x+']' in text for x in ['L','S','A','M','Q','P','U']))
        check(lang+' source hash',source['script_sha256'] in text)
        check(lang+' no unexpanded token','@@' not in text);check(lang+' balanced fences',text.count('```')%2==0);check(lang+' balanced math',text.count('$$')%2==0)
        for link in re.findall(r'\]\(([^)]+)\)',text):
            if not link.startswith(('http://','https://','#')):
                target=(p.parent/link).resolve()
                check(lang+' link '+link,target.is_file() or target==(ROOT/'results/publication_validation.json').resolve())
    def decimals(text):return [re.findall(r'[-−]?\d+\.\d+',line) for line in text.splitlines() if line.startswith('|')]
    check('English Chinese table numeric parity',decimals(contents['english'])==decimals(contents['chinese']))
    check('English edition no Chinese prose',not re.search(r'[\u4e00-\u9fff]',contents['english']))
    for lang in ['english','chinese']:
        bodies=[re.split(r'^## ',b,flags=re.M)[0].strip() for b in re.split(r'^### [^\n]+\n',contents[lang],flags=re.M)[1:]]
        check(lang+' full subsection bodies retained in bilingual',all(b in contents['bilingual'] for b in bodies))
    for p in [ROOT/'README.md',ROOT/'inputs/reviewed/README.md',ROOT/'data/README.md']:
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
            if not link.startswith(('http://','https://','#')) and not link.endswith('publication_validation.json'):check('README link '+link,(p.parent/link).exists())
    for p in ROOT.rglob('*.json'):
        if p.name=='publication_validation.json':continue
        try:json.loads(p.read_text(encoding='utf-8-sig'));ok=True
        except (ValueError,OSError):ok=False
        check('JSON parses '+p.relative_to(ROOT).as_posix(),ok)
    result=dict(passed=all(c['passed'] for c in checks),checks=len(checks),failed=[c for c in checks if not c['passed']],scope='Record, structure, bilingual completeness and numeric integrity; not chemical or assay validation',items=checks)
    (ROOT/'results/publication_validation.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:result[k] for k in ['passed','checks','failed']},indent=2));return 0 if result['passed'] else 1
if __name__=='__main__':sys.exit(main())
