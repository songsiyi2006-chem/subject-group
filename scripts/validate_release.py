"""Validate the saved publication, without interpreting software success as chemistry."""
from pathlib import Path
import hashlib,json,re,unittest,sys,io,subprocess
ROOT=Path(__file__).resolve().parents[1]
def main():
    checks=[]
    def check(name,condition):checks.append({'check':name,'passed':bool(condition)})
    reports=[ROOT/'reports'/f'technical_report_{lang}.md' for lang in ['english','chinese']]
    contents=[p.read_text(encoding='utf-8') for p in reports]
    for path,txt in zip(reports,contents):
        tag=path.stem
        check(tag+' main architecture',all(re.search(rf'^## {i}\. ',txt,re.M) for i in [1,2,3]))
        check(tag+' 20 topic subsections',all(re.search(rf'^#### 2\.{i}\.{j} ',txt,re.M) for i in range(1,6) for j in range(1,5)))
        check(tag+' four undergraduate stages',all(re.search(rf'^### 3\.{i} ',txt,re.M) for i in range(1,5)))
        check(tag+' added methods',all(re.search(rf'^### {section}\.{i} ',txt,re.M) for section,limit in [('D',7),('E',4)] for i in range(1,limit+1)))
        check(tag+' original numerical anchors',all(x in txt for x in ['95.04474861433413','0.9984612448984918','18.09','0.892','96.98','74.36']))
        check(tag+' extension numerical anchors',all(x in txt for x in ['3.5557','1.3570','1.3208','0.9287','0.9921','5.841','0.9129','0.7276','0.55','10.7446','12.2279']))
        check(tag+' scientific evidence labels',all(x in txt for x in ['[L]','[C]','[Q]','[A]','[P]','[U]']))
        check(tag+' no unresolved template',not re.search(r'@@[A-Z]+@@',txt))
        check(tag+' balanced fenced code',txt.count('```')%2==0)
        check(tag+' balanced display math',txt.count('$$')%2==0)
        for link in re.findall(r'\]\(([^)]+)\)',txt):
            if not link.startswith(('http://','https://','#')):check(tag+' local link '+link,(path.parent/link).exists())
    # Decimal measurements must match row by row. Integer topic names can be
    # written as Chinese words (e.g. 2D NMR vs its Chinese translation).
    table_numbers=[]
    for txt in contents:
        nums=[]
        for line in txt.splitlines():
            if line.startswith('|'):nums.append(re.findall(r'[-−]?\d+\.\d+',line))
        table_numbers.append(nums)
    check('English/Chinese table decimal-value parity by row',table_numbers[0]==table_numbers[1])
    check('English edition has no Chinese prose',not re.search(r'[\u4e00-\u9fff]',contents[0]))
    production=subprocess.run([sys.executable,str(ROOT/'production/scripts/validate_production.py')],capture_output=True,text=True,encoding='utf-8')
    check('four-task production publication validation',production.returncode==0)
    if production.returncode:print(production.stdout,production.stderr)
    research=subprocess.run([sys.executable,str(ROOT/'research/scripts/validate_research.py')],capture_output=True,text=True,encoding='utf-8')
    check('three-module research publication validation',research.returncode==0)
    if research.returncode:print(research.stdout,research.stderr)
    closed_loop=subprocess.run([sys.executable,str(ROOT/'closed_loop/scripts/validate_closed_loop.py')],capture_output=True,text=True,encoding='utf-8')
    check('closed-loop publication validation',closed_loop.returncode==0)
    if closed_loop.returncode:print(closed_loop.stdout,closed_loop.stderr)
    toolkit=subprocess.run([sys.executable,str(ROOT/'toolkit/scripts/validate_toolkit.py')],capture_output=True,text=True,encoding='utf-8')
    check('analytical toolkit publication validation',toolkit.returncode==0)
    if toolkit.returncode:print(toolkit.stdout,toolkit.stderr)
    electratwin=subprocess.run([sys.executable,str(ROOT/'electratwin/scripts/validate_electratwin.py')],capture_output=True,text=True,encoding='utf-8')
    check('ElectraTwin publication validation',electratwin.returncode==0)
    if electratwin.returncode:print(electratwin.stdout,electratwin.stderr)
    electrograph=subprocess.run([sys.executable,str(ROOT/'electrograph/scripts/validate_electrograph.py')],capture_output=True,text=True,encoding='utf-8')
    check('ElectroGraph publication validation',electrograph.returncode==0)
    if electrograph.returncode:print(electrograph.stdout,electrograph.stderr)
    synthapore=subprocess.run([sys.executable,str(ROOT/'synthapore/scripts/validate_synthapore.py'),'--require-cross-review'],capture_output=True,text=True,encoding='utf-8')
    check('SynthaPore publication validation',synthapore.returncode==0)
    if synthapore.returncode:print(synthapore.stdout,synthapore.stderr)
    quantumequi=subprocess.run([sys.executable,str(ROOT/'quantumequi/scripts/validate_quantumequi.py'),'--require-cross-review'],capture_output=True,text=True,encoding='utf-8')
    check('QuantumEqui publication validation',quantumequi.returncode==0)
    if quantumequi.returncode:print(quantumequi.stdout,quantumequi.stderr)
    extensions=subprocess.run([sys.executable,str(ROOT/'quantumequi/scripts/validate_extensions.py')],capture_output=True,text=True,encoding='utf-8')
    check('QuantumEqui extended calculations validation',extensions.returncode==0)
    if extensions.returncode:print(extensions.stdout,extensions.stderr)
    manuscript=subprocess.run([sys.executable,str(ROOT/'manuscript/scripts/validate_manuscript.py'),'--require-documents'],capture_output=True,text=True,encoding='utf-8')
    check('integrated manuscript numerical and publication validation',manuscript.returncode==0)
    expansion=subprocess.run([sys.executable,str(ROOT/'manuscript/expanded/scripts/validate_expansion.py')],capture_output=True,text=True,encoding='utf-8')
    check('expanded historical manuscript validation',expansion.returncode==0)
    if expansion.returncode:print(expansion.stdout,expansion.stderr)
    if manuscript.returncode:print(manuscript.stdout,manuscript.stderr)
    files=[p for d in ['results','provenance','data/original','production/results','production/source','research/results','research/source','closed_loop/results','closed_loop/source','closed_loop/inputs','toolkit/results','toolkit/source','electratwin/results','electratwin/source','electrograph/results','electrograph/source','electrograph/data','synthapore/results','synthapore/source','quantumequi/results','quantumequi/source','manuscript'] for p in (ROOT/d).rglob('*.json') if p.name!='release_validation.json']
    for f in files:
        try:json.loads(f.read_text(encoding='utf-8-sig'));ok=True
        except Exception:ok=False
        check('JSON parses '+f.relative_to(ROOT).as_posix(),ok)
    # Check only publishable files, excluding git metadata, ignored caches and private scratch.
    sensitive=[]
    for p in ROOT.rglob('*'):
        if not p.is_file() or any(x in {'.git','__pycache__','work'} for x in p.relative_to(ROOT).parts):continue
        if p.suffix in {'.md','.json','.log','.txt','.csv','.py','.svg'} and p.name!='validate_release.py':
            txt=p.read_text(encoding='utf-8-sig',errors='replace')
            if re.search(r'[A-Z]:[\\/]Users[\\/]|[A-Z]:[\\/]Codex[\\/]|github_pat_|ghp_[A-Za-z0-9]{20}',txt):sensitive.append(p.relative_to(ROOT).as_posix())
    check('public files omit machine user paths and credential patterns',not sensitive)
    manifest=ROOT/'provenance/SHA256SUMS.txt'
    if manifest.exists():
        for line in manifest.read_text().splitlines():
            digest,name=line.split('  ',1);p=ROOT/name
            check('SHA256 '+name,p.exists() and hashlib.sha256(p.read_bytes()).hexdigest()==digest)
    # Running this file puts scripts/ on sys.path; namespace-package tests need
    # the repository root just as they do under python -m unittest.
    sys.path.insert(0,str(ROOT))
    suite=unittest.defaultTestLoader.discover(str(ROOT/'tests'))
    stream=io.StringIO();result=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    check('scientific unittest suite',result.wasSuccessful())
    # Preserve concise portable successful test output; on failure, surface diagnostics.
    if result.wasSuccessful():(ROOT/'provenance/test_results.txt').write_text(stream.getvalue(),encoding='utf-8')
    else:print(stream.getvalue())
    record={'passed':all(c['passed'] for c in checks),'checks':len(checks),'scientific_tests_run':result.testsRun,'scientific_failures':len(result.failures),'scientific_errors':len(result.errors),'failed':[c for c in checks if not c['passed']],'path_scan_findings':sensitive,'scope':'Record, architecture and numerical integrity; not chemical or experimental validation.','items':checks}
    (ROOT/'provenance/release_validation.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print(json.dumps({k:record[k] for k in ['passed','checks','scientific_tests_run','failed','path_scan_findings']},indent=2))
    return 0 if record['passed'] else 1
if __name__=='__main__':sys.exit(main())
