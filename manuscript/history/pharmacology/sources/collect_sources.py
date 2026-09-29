"""Freeze selected public git blobs; no scientific jobs or original checkout edits."""
import hashlib
import json
from pathlib import Path
import subprocess
from concurrent.futures import ThreadPoolExecutor

HERE = Path(__file__).resolve().parent
CLONES = HERE.parents[4] / 'history_sources'
REPOSITORIES = {
    'ai4pharm-lead-developability-suite': 'a7a2a4df501d5dad030dc488691b02c73fdb423f',
    'Paper-Analysis-ChemMLLM': 'c1d309a0c2f876ccc8b4d2ce753ece971d58fe97',
}

def git(repo, *args):
    for attempt in range(3):
        q = subprocess.run(['git', '-c', 'http.sslBackend=openssl', *args], cwd=CLONES/repo, capture_output=True)
        if q.returncode == 0:
            return q.stdout
    raise RuntimeError(q.stderr.decode('utf-8', 'replace'))

PHARM_PATHS = [
    'LICENSE', 'README.md',
    'projects/task01_lead_developability/MODEL_CARD.md',
    'projects/task01_lead_developability/data/reference_panel.json',
    'projects/task01_lead_developability/results_task1/developability_results.csv',
    'projects/task01_lead_developability/results_task1/developability_results.json',
    'projects/task01_lead_developability/results_task1/group_summary.csv',
    'projects/task01_lead_developability/results_task1/run_manifest.json',
    'projects/task01_lead_developability/run_task1_mpo_admet_developability.py',
    'projects/task02_tpd/run_task2_tpd_ternary_cooperativity.py',
    'projects/task02_tpd/manifest_task2.json',
    'projects/task02_tpd/data_task2/equilibrium_metrics.csv',
    'projects/task02_tpd/data_task2/equilibrium_species.csv',
    'projects/task02_tpd/data_task2/linker_summary.csv',
    'projects/task02_tpd/data_task2/linker_conformers.csv',
    'projects/task03_covalent_kinetics/run_task3_covalent_kinetics_residence_time.py',
    'projects/task03_covalent_kinetics/data_task3/warhead_results.csv',
    'projects/task03_covalent_kinetics/data_task3/kobs_concentration.csv',
    'projects/task03_covalent_kinetics/data_task3/run_manifest.json',
    'projects/task03_covalent_kinetics/data_task3/lfer_regression.json',
    'projects/task19_metadynamics/outputs/executed_driver.py',
    'projects/task19_metadynamics/outputs/summary.json',
    'projects/task19_metadynamics/outputs/run_metadata.json',
    'projects/task19_metadynamics/outputs/config.json',
    'projects/task19_metadynamics/outputs/atomistic_cv_trajectory.csv',
    'projects/task19_metadynamics/outputs/hill_history.csv',
    'projects/task19_metadynamics/outputs/verification.json',
    'projects/task19_metadynamics/inputs/5T35.pdb',
    'projects/task19_metadynamics/EXECUTION_AUDIT.json',
    'projects/task20_denovo/driver.py',
    'projects/task20_denovo/outputs/summary.json',
    'projects/task20_denovo/outputs/generation_summary.csv',
    'projects/task20_denovo/outputs/evaluated_candidates.csv',
    'projects/task20_denovo/outputs/random_search_control.csv',
    'projects/task20_denovo/outputs/random_search_summary.json',
    'projects/task20_denovo/outputs/docking/candidate_vina_scores.csv',
    'projects/task20_denovo/outputs/verification.json',
    'data_omnibus/recompute/projects/task01_lead_developability/results_task1/developability_results.csv',
    'data_omnibus/recompute/projects/task02_tpd/data_task2/equilibrium_metrics.csv',
    'data_omnibus/recompute/projects/task03_covalent_kinetics/data_task3/warhead_results.csv',
]

def collect(item):
    repo, path = item
    commit = REPOSITORIES[repo]
    data = git(repo, 'show', commit+':'+path)
    target = HERE/repo/path
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(data)
    return dict(repository=repo, commit=commit, path=path, copied_path=target.relative_to(HERE).as_posix(),
                sha256=hashlib.sha256(data).hexdigest(), bytes=len(data),
                url=f'https://github.com/songsiyi2006-chem/{repo}/blob/{commit}/{path}')

def main():
    repos=[]
    for repo, commit in REPOSITORIES.items():
        paths=git(repo,'ls-tree','-rz','--name-only',commit).decode().split('\0')[:-1]
        log=git(repo,'log','--format=%H%x09%aI%x09%s',commit).decode().splitlines()
        repos.append(dict(repository=repo,commit=commit,tracked_paths=paths,
                          tracked_files=len(paths),commit_count=len(log),
                          history=[dict(zip(['commit','author_time','subject'],x.split('\t',2))) for x in log],
                          license='MIT' if repo.startswith('ai4pharm') else 'No license file found in tracked tree'))
    requests=[('ai4pharm-lead-developability-suite',p) for p in PHARM_PATHS]
    requests.append(('Paper-Analysis-ChemMLLM','README.md'))
    files=list(ThreadPoolExecutor(2).map(collect,requests))
    result=dict(repositories=repos,files=files,note='Exact public blobs. Published papers/PDFs/figures were not copied. Paper README retained solely as minimal user-owned repository provenance; no open license inferred.')
    (HERE/'source_manifest.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'files':len(files),'bytes':sum(x['bytes'] for x in files),'repositories':[{k:r[k] for k in ['repository','commit','commit_count','tracked_files']} for r in repos]}))

if __name__ == '__main__':
    main()
