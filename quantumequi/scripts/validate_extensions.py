"""Read-only scientific arithmetic and artifact checks for the new extension.

No SCF, neural training, or nuclear eigensolver study is rerun by this validator.
This validates the frozen historical release, including its recovery ledger;
new scientific runs need their own accounting and renewed QA records.
"""
from collections import defaultdict
import csv,hashlib,json,re,sys
from pathlib import Path
from urllib.parse import unquote
import xml.etree.ElementTree as ET
import numpy as np
from scipy.special import logsumexp
from PIL import Image

BASE=Path(__file__).resolve().parents[1]
DATA=BASE/'results/extensions'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(s):return json.loads((DATA/s).read_text(encoding='utf-8-sig'))
def rows(s):
    with (DATA/s).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def n(r,k):return float(r[k])
def close(a,b):return bool(np.isclose(float(a),float(b),rtol=1e-8,atol=1e-11))


def main():
    checks=[]
    def check(name,ok):checks.append({'check':name,'passed':bool(ok)})
    def hashes(mapping,base,label):
        for name,value in mapping.items():
            path=base/name;check(label+' '+name,path.is_file() and sha(path)==value)
    def guard(name,fn):
        try:fn()
        except Exception as exc:checks.append({'check':name,'passed':False,'detail':str(exc).replace(str(BASE),'quantumequi')})

    def correlation():
        for part in ('correlation','correlation/pilot','correlation/pilot_fci_recovery'):
            s=read(part+'/summary.json');snapshot=DATA/part/'executed_code.py.txt'
            hashes(s['outputs_sha256'],DATA/part,part+' output')
            check(part+' exact executed source',sha(snapshot)==s['inputs_sha256']['scripts/correlation_extension.py'])
            check(part+' current finalizer hash',sha(BASE/'scripts/correlation_extension.py')==s['finalizer_source_sha256'])
        curve=rows('correlation/curve.csv')
        check('150 converged method-basis-geometry points',len(curve)==150 and all(r['status']=='converged' for r in curve)
              and len({(r['basis'],r['method'],r['R_A']) for r in curve})==150)
        by=defaultdict(dict)
        for r in curve:by[(r['basis'],r['R_A'])][r['method']]=r
        for key,values in by.items():
            e={m:n(r,'energy_Hartree') for m,r in values.items()}
            check('within-basis variational order '+str(key),set(e)=={'RHF','UHF','FCI'} and e['FCI']<=min(e['RHF'],e['UHF'])+1e-8 and e['UHF']<=e['RHF']+1e-8)
        check('UHF stretched spin breaking is recorded',any(n(r,'S2')>.9 for r in curve if r['method']=='UHF' and n(r,'R_A')>=3))

    def learning():
        for part in ('learning','learning/pilot'):
            s=read(part+'/summary.json');hashes(s['output_sha256'],DATA/part,part+' output')
            check(part+' source snapshot',(sha(DATA/part/'executed_code.py.txt')==s['code_sha256']))
        s=read('learning/summary.json');check('current learning source frozen',sha(BASE/'scripts/learning_extension.py')==s['code_sha256'])
        check('six matched neural runs and no new quantum calls',s['counts']['training_runs']==6 and s['counts']['optimizer_steps']==4800 and s['counts']['new_quantum_jobs']==0)
        for seed in (7301,7302,7303):
            runs=[r for r in s['runs'] if r['seed']==seed]
            check('paired initialization '+str(seed),len(runs)==2 and len({r['initial_state_sha256'] for r in runs})==1)
        curves=rows('learning/learning_curves.csv');pred=rows('learning/predictions.csv');metrics=rows('learning/metrics.csv')
        check('learning table sizes',len(curves)==4800 and len(pred)==252 and len(metrics)==24)
        for run in s['runs']:
            c=[r for r in curves if r['objective']==run['objective'] and int(r['seed'])==run['seed']]
            selected=min(c,key=lambda r:n(r,'selection_score'))
            check('validation-only checkpoint '+run['objective']+str(run['seed']),int(selected['epoch'])==run['selected_epoch'] and close(n(selected,'selection_score'),run['selected_validation_score']))
        for row in metrics:
            q=[r for r in pred if all(r[k]==row[k] for k in ('objective','seed','split'))]
            for label in ('energy','gradient'):
                unit='Hartree' if label=='energy' else 'Hartree_A'
                e=np.array([n(r,'predicted_'+label+'_'+unit)-n(r,'reference_'+label+'_'+unit) for r in q])
                check('independent '+label+' metrics '+str((row['objective'],row['seed'],row['split'])),len(e)==int(row['points']) and close(np.sqrt(np.mean(e**2)),row[label+'_RMSE_'+unit]) and close(np.mean(abs(e)),row[label+'_MAE_'+unit]))
        for row in rows('learning/ensemble_predictions.csv'):
            q=[r for r in pred if r['objective']==row['objective'] and r['point_id']==row['point_id']]
            check('three-seed ensemble '+row['objective']+row['point_id'],len(q)==3 and all(close(np.mean([n(r,'predicted_'+label+'_'+unit) for r in q]),row['ensemble_'+label+'_mean_'+unit]) and close(np.std([n(r,'predicted_'+label+'_'+unit) for r in q],ddof=1),row['ensemble_'+label+'_std_'+unit]) for label,unit in [('energy','Hartree'),('gradient','Hartree_A')]))
        pilot=rows('learning/pilot/predictions.csv')
        check('pilot excluded held-out predictions',{r['split'] for r in pilot}<={'train','validation'})

    def vibration():
        for part in ('vibration','vibration/pilot','vibration/tail_followup'):
            s=read(part+'/summary.json');hashes(s['outputs_sha256'],DATA/part,part+' output')
            hashes(s['inputs_sha256'],BASE.parent,part+' input')
            snapshot=s['executed_source_snapshot'];check(part+' executed source',sha(DATA/part/snapshot['file'])==snapshot['sha256'])
            if part=='vibration/tail_followup':check('current vibration code frozen',sha(BASE/'scripts/vibration_extension.py')==snapshot['sha256'])
        s=read('vibration/summary.json');table=rows('vibration/isotope_levels.csv');spectra=rows('vibration/isotope_summary.csv');parts=rows('vibration/partition_bound_only.csv')
        check('six isotope-model summaries',len(spectra)==6)
        for r in spectra:
            q=sorted([v for v in table if v['model_id']==r['model_id'] and v['isotope']==r['isotope']],key=lambda v:int(v['v']))
            e=np.array([n(v,'E_numerical_Hartree') for v in q]);check('spectrum identities '+r['model_id']+r['isotope'],len(q)==int(r['bound_states_numerical']) and np.all(np.diff(e)>0) and close(e[0],r['ZPE_numerical_Hartree']) and close((e[1]-e[0])*s['units']['Hartree_to_cm_inverse'],r['gap01_numerical_cm_1']))
        for r in parts:
            q=[v for v in table if v['model_id']==r['model_id'] and v['isotope']==r['isotope']]
            e=np.array([n(v,'E_numerical_Hartree') for v in q]);kt=1.380649e-23*n(r,'T_K')/4.359744722206e-18
            lq=logsumexp(-e/kt);prob=np.exp(-e/kt-lq);u=np.dot(prob,e);f=-kt*lq
            entropy=(u-f)/n(r,'T_K')*4.359744722206e-18*6.02214076e23
            check('independent bound-state partition '+r['model_id']+r['isotope']+r['T_K'],close(lq,r['log_q_relative_to_PES_minimum']) and close(u,r['U_vib_Hartree']) and close(f,r['F_vib_bound_only_Hartree']) and close(entropy,r['S_vib_J_mol_K']) and r['continuum_excluded']=='True')

    def decomposition():
        s=read('error_budget/summary.json');hashes(s['inputs_sha256'],BASE,'error budget input');hashes(s['outputs_sha256'],DATA/'error_budget','error budget output')
        table=rows('error_budget/matched_error_components.csv')
        predicted={(r['objective'],r['seed'],r['point_id']):r for r in rows('learning/predictions.csv')}
        with (BASE/'results/electronic/h2_reference.csv').open(encoding='utf-8-sig',newline='') as f:
            reference={r['point_id']:r for r in csv.DictReader(f)}
        fresh={(r['basis'],r['method'],round(n(r,'R_A'),8)):r for r in rows('correlation/curve.csv')}
        check('matched geometry arithmetic rows',len(table)==s['rows'] and len({r['R_A'] for r in table})==s['unique_geometries'])
        for row in table:
            check('energy-error identity '+row['objective']+row['seed']+row['point_id'],abs(n(row,'total_error_vs_FCI_Hartree')-n(row,'learning_error_Hartree')-n(row,'reference_replay_drift_Hartree')-n(row,'RHF_minus_FCI_Hartree'))<1e-12)
            p=predicted[(row['objective'],row['seed'],row['point_id'])];old=reference[row['point_id']]
            radius=round(n(row,'R_A'),8);rhf=fresh[('sto-3g','RHF',radius)];fci=fresh[('sto-3g','FCI',radius)]
            ep,eo,er,ef=n(p,'predicted_energy_Hartree'),n(old,'energy_total_Hartree'),n(rhf,'energy_Hartree'),n(fci,'energy_Hartree')
            check('independent four-energy reconstruction '+row['objective']+row['seed']+row['point_id'],
                  close(n(p,'R_A'),radius) and close(n(old,'R_A'),radius) and row['split']==old['split']==p['split'] and
                  all(close(actual,row[key]) for actual,key in [(ep-eo,'learning_error_Hartree'),(eo-er,'reference_replay_drift_Hartree'),(er-ef,'RHF_minus_FCI_Hartree'),(ep-ef,'total_error_vs_FCI_Hartree')]))

    def publication():
        manifest=read('figure_manifest.json');figures=manifest['figures'];check('28 static figures',len(figures)==28)
        check('current plotting code',sha(BASE/'scripts/plot_extensions.py')==manifest['generator_sha256'])
        for row in figures:
            file=BASE/'reports/figures/extensions'/row['file'];check('figure '+row['file'],sha(file)==row['sha256']);hashes(row['sources'],DATA,row['file']+' data')
            if row['format']=='png':
                with Image.open(file) as im:check('PNG layout '+row['file'],list(im.size)==[2700,1440] and all(abs(x-300)<.1 for x in im.info['dpi']))
            else:
                tree=ET.parse(file);ns={'s':'http://www.w3.org/2000/svg'};check('editable SVG '+row['file'],len(tree.findall('.//s:text',ns))>10 and not tree.findall('.//s:image',ns))
        qa=read('figure_qa.json');check('all static figures visually reviewed',qa['png_files_visually_reviewed']==14 and qa['svg_files_rendered_and_visually_reviewed']==14 and qa['reviewed_figure_hashes']=={r['file']:r['sha256'] for r in figures})
        explorer=read('explorer_manifest.json');hashes(explorer['inputs_sha256'],BASE,'explorer input');check('explorer output',sha(BASE/'reports/quantum_explorer.html')==explorer['output_sha256'] and explorer['rows']==150)
        html=(BASE/'reports/quantum_explorer.html').read_text(encoding='utf-8')
        embedded=json.loads(re.search(r'<script id="dataset" type="application/json">(.*?)</script>',html,re.S).group(1))
        saved={(r['basis'],r['method'],n(r,'R_A')):r for r in rows('correlation/curve.csv')}
        check('embedded quantum point identities',len(embedded)==len(saved)==150 and {(r['basis'],r['method'],r['R_A']) for r in embedded}==set(saved))
        for point in embedded:
            original=saved[(point['basis'],point['method'],point['R_A'])]
            check('embedded energy and spin '+str((point['basis'],point['method'],point['R_A'])),
                  close(point['energy_Hartree'],original['energy_Hartree']) and
                  ((point['S2'] is None and original['S2']=='') or (point['S2'] is not None and close(point['S2'],original['S2']))))
        browser=read('explorer_qa.json');check('interactive page reviewed',browser['html_sha256']==explorer['output_sha256'] and browser['passed'])
        editions=[]
        for lang in ('english','chinese'):
            p=BASE/'reports'/('extension_report_'+lang+'.md');text=p.read_text(encoding='utf-8');editions.append(text)
            check(lang+' six sections',re.findall(r'^## (\d+)\.',text,re.M)==list('123456'))
            check(lang+' seven PNG and SVG links',len(re.findall(r'\]\(figures/extensions/[^)]+\.png\)',text))==7 and len(re.findall(r'\]\(figures/extensions/[^)]+\.svg\)',text))==7)
            check(lang+' balanced code',text.count('```')%2==0)
            for link in re.findall(r'\]\(([^)]+)\)',text):
                if not link.startswith(('https://','http://','#')):check(lang+' link '+link,(p.parent/unquote(link.split('#')[0])).exists())
        numbers=lambda t:[re.findall(r'[-−]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?',l) for l in t.splitlines() if l.startswith('|')]
        check('bilingual table numerical parity',numbers(editions[0])==numbers(editions[1]));check('English report has no Chinese prose',not re.search(r'[\u4e00-\u9fff]',editions[0]))
        review=read('cross_review.json');check('cross-review no unresolved blockers',bool(review['reviews']) and all(r['blocking_findings']==[] for r in review['reviews']));hashes(review['reviewed_artifact_sha256'],BASE,'cross-review')

    for name,fn in [('correlation',correlation),('learning',learning),('vibration',vibration),('decomposition',decomposition),('publication',publication)]:guard(name,fn)
    record={'passed':all(r['passed'] for r in checks),'checks':len(checks),'failed':[r for r in checks if not r['passed']],
            'scope':'Saved numerical evidence, publication and QA; no chemical or experimental validation.',
            'validator_sha256':sha(Path(__file__)),'items':checks}
    (DATA/'validation.json').write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({k:record[k] for k in ('passed','checks','failed')},indent=2,ensure_ascii=False));return 0 if record['passed'] else 1


if __name__=='__main__':sys.exit(main())
