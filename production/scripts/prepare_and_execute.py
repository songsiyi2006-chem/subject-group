"""Preserve the supplied failure and run a documented, minimally repaired copy."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import difflib,hashlib,json,subprocess,sys,time
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def clean(text):
    for value in sorted({str(ROOT.resolve()),sys.executable,sys.prefix,sys.base_prefix},key=len,reverse=True):
        text=text.replace(value,'<local-runtime-or-project>')
    return text
def execute(script,folder):
    folder.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    proc=subprocess.run([sys.executable,str(script)],cwd=folder,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=300)
    (folder/'stdout.log').write_text(clean(proc.stdout),encoding='utf-8')
    (folder/'stderr.log').write_text(clean(proc.stderr),encoding='utf-8')
    record=dict(script=script.relative_to(ROOT).as_posix(),sha256=hashlib.sha256(script.read_bytes()).hexdigest(),returncode=proc.returncode,elapsed_s=time.perf_counter()-start,timestamp_utc=datetime.now(timezone.utc).isoformat(),result_written=(folder/'production_benchmark_results.json').exists(),threads=1,timeout_s=300)
    (folder/'execution.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    return record
def main():
    original=ROOT/'source/run_production_pipeline_original.py'
    original_record=execute(original,ROOT/'results/original_attempt')
    print('Original exit:',original_record['returncode'],flush=True)
    source=original.read_text(encoding='utf-8');repaired=source
    substitutions={
        'AllChem.EmbedMolecule(mol, params)':'embed_status = AllChem.EmbedMolecule(mol, params)\n    if embed_status != 0: raise RuntimeError("ETKDG embedding failed")',
        'AllChem.MMFFOptimizeMolecule(mol, maxIters=500)':'mmff_status = AllChem.MMFFOptimizeMolecule(mol, maxIters=500)\n    if mmff_status != 0: raise RuntimeError("MMFF geometry did not converge")',
        'radii = rdFreeSASA.classifyAtoms(mol)':'radii = [Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()) for atom in mol.GetAtoms()]',
        'sasa_per_atom = rdFreeSASA.calcSASA(mol, radii)':'total_sasa = rdFreeSASA.CalcSASA(mol, radii)',
        'h_sasa = rdFreeSASA.calcSASA(mol, radii, whichAtoms=[h_atom.GetIdx()])':'h_sasa = h_atom.GetDoubleProp("SASA")',
        'best_idx_overall = np.argmax(y_train)':'pd.DataFrame(df_space).to_csv("task_a_candidates.csv", index=False)\n    sampled = df_space.iloc[init_indices].copy()\n    sampled["observed_yield"] = y_train\n    sampled.to_csv("task_a_observations.csv", index=False)\n    best_idx_overall = np.argmax(y_train)',
        'top_site = df_carbons.iloc[0].to_dict()':'df_carbons.to_csv("task_b_all_ch_sites.csv", index=False)\n    Chem.MolToMolFile(mol, "task_b_mmff_geometry.mol")\n    pd.DataFrame([dict(atom_index=a.GetIdx(), element=a.GetSymbol(), radius_A=radii[a.GetIdx()], sasa_A2=a.GetDoubleProp("SASA"), gasteiger_charge=float(a.GetProp("_GasteigerCharge"))) for a in mol.GetAtoms()]).to_csv("task_b_all_atoms.csv", index=False)\n    top_site = df_carbons.iloc[0].to_dict()',
        '"heavy_atom_count": mol.GetNumHeavyAtoms(),':'"heavy_atom_count": mol.GetNumHeavyAtoms(),\n        "total_sasa_A2": total_sasa,\n        "embed_status": embed_status,\n        "mmff_status": mmff_status,\n        "sasa_method": "Lee-Richards; RDKit periodic-table vdW radii; explicit H; probe 1.4 A",\n        "score_evidence": "Uncalibrated heuristic; not a quantum or experimental site prediction",',
        'df_sac = pd.DataFrame(records)':'df_sac = pd.DataFrame(records)\n    df_sac.to_csv("task_c_all_sites.csv", index=False)',
        'Production calculations completed without mock primitives.':'Compatibility-repaired demonstrations completed; A/C are synthetic and B scoring is heuristic.',
        'High-fidelity scientific metrics exported to:':'Mixed-evidence numerical outputs exported to:'
    }
    for old,new in substitutions.items():
        if repaired.count(old)!=1:raise ValueError('Patch anchor mismatch: '+old)
        repaired=repaired.replace(old,new)
    header='# AUDITED COPY: A/C targets are synthetic, B ranking is heuristic, D is a reduced ODE.\n# Original source and explicit patch are preserved; no experimental validation is claimed.\n'
    repaired=header+repaired
    fixed=ROOT/'run_production_pipeline.py';fixed.write_text(repaired,encoding='utf-8',newline='\n')
    diff=''.join(difflib.unified_diff(source.splitlines(True),repaired.splitlines(True),fromfile='source/run_production_pipeline_original.py',tofile='run_production_pipeline.py'))
    (ROOT/'source/repair.patch').write_text(diff,encoding='utf-8',newline='\n')
    record=execute(fixed,ROOT/'results/repaired');print('Repaired exit:',record['returncode'],flush=True)
    if record['returncode']:raise SystemExit(record['returncode'])
if __name__=='__main__':main()
