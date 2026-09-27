"""Small real-molecule RDKit/GFN2-xTB pilot. No experimental Eox or selectivity.

Set XTB_EXE if xtb is not on PATH. Conda Windows users must include its
Library/bin directory on PATH. Each xTB job is serial and timeout bounded.
"""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ.setdefault(key,'1')
import argparse,hashlib,json,platform,re,shutil,subprocess,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem,rdBase
from rdkit.Chem import AllChem,Descriptors,rdMolDescriptors

MOLECULES=[('benzene','c1ccccc1'),('pyridine','c1ccncc1'),('anisole','COc1ccccc1'),('indole','c1ccc2[nH]ccc2c1'),('N_methylindole','Cn1ccc2ccccc21'),('benzofuran','c1ccc2occc2c1')]
HARTREE_EV=27.211386245988

def get_energy(text):
    vals=re.findall(r'TOTAL ENERGY\s+([-+\d.]+)\s+Eh',text)
    if not vals:raise ValueError('No total-energy record in xTB output')
    return float(vals[-1])

def run_xtb(executable,folder,args,timeout):
    folder.mkdir(parents=True,exist_ok=True)
    t0=time.perf_counter()
    command=[executable,*args,'--parallel','1','--norestart']
    try:
        proc=subprocess.run(command,cwd=folder,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,encoding='utf-8',errors='replace',timeout=timeout)
        stdout,stderr=proc.stdout,proc.stderr
        # Keep portable raw logs: remove only this exact scratch path if printed.
        stdout=stdout.replace(executable,'xtb');stderr=stderr.replace(executable,'xtb')
        for path in [str(folder.resolve()),str(folder.resolve()).replace('\\','/')]:
            stdout=stdout.replace(path,'<job-directory>');stderr=stderr.replace(path,'<job-directory>')
        (folder/'stdout.log').write_text(stdout,encoding='utf-8');(folder/'stderr.log').write_text(stderr,encoding='utf-8')
        # xTB 6.7.1 on Windows writes unescaped backslashes in program-call JSON.
        # Replace only the machine-specific executable path; numeric output is unchanged.
        json_path=folder/'xtbout.json'
        if json_path.exists():
            json_text=json_path.read_text(encoding='utf-8').replace(executable,'xtb')
            json.loads(json_text)
            json_path.write_text(json_text,encoding='utf-8')
        ok=proc.returncode==0 and 'normal termination of xtb' in (stdout+stderr) and 'abnormal termination' not in (stdout+stderr).lower()
        record=dict(args=['xtb',*args,'--parallel','1','--norestart'],returncode=proc.returncode,elapsed_s=time.perf_counter()-t0,completed=ok)
        record['public_log_normalization']='Executable and job-directory paths replaced by portable names; scientific values unchanged.'
        if '--opt' in args:record['optimization_converged']='GEOMETRY OPTIMIZATION CONVERGED' in stdout
        if ok:record['energy_Eh']=get_energy(stdout)
    except subprocess.TimeoutExpired as exc:
        stdout=(exc.stdout or b'');stderr=(exc.stderr or b'')
        (folder/'stdout.log').write_bytes(stdout if isinstance(stdout,bytes) else stdout.encode())
        (folder/'stderr.log').write_bytes(stderr if isinstance(stderr,bytes) else stderr.encode())
        record=dict(args=['xtb',*args],completed=False,error='timeout',timeout_s=timeout,elapsed_s=time.perf_counter()-t0)
    (folder/'run.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    return record

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=Path('results/molecular'));p.add_argument('--timeout',type=int,default=180);p.add_argument('--limit',type=int,default=6);args=p.parse_args()
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    executable=os.environ.get('XTB_EXE') or shutil.which('xtb')
    if not executable:raise SystemExit('xTB not found; set XTB_EXE. Existing RDKit inputs are not quantum results.')
    version=subprocess.run([executable,'--version'],capture_output=True,text=True,encoding='utf-8',errors='replace').stdout
    version_line=next((x.strip() for x in version.splitlines() if 'xtb version' in x),'unknown')
    rows=[];atom_rows=[];failures=[]
    for name,smiles in MOLECULES[:args.limit]:
        print('Molecular calculation: '+name,flush=True)
        folder=out/name;folder.mkdir(exist_ok=True)
        base=Chem.MolFromSmiles(smiles);mol=Chem.AddHs(base)
        params=AllChem.ETKDGv3();params.randomSeed=31415;params.numThreads=1
        ids=AllChem.EmbedMultipleConfs(mol,numConfs=4,params=params)
        if not ids:failures.append(dict(molecule=name,error='embedding failed'));continue
        opt=AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=1,maxIters=500)
        good=[(energy,int(cid)) for cid,(status,energy) in zip(ids,opt) if status==0]
        if not good:failures.append(dict(molecule=name,error='MMFF conformers did not converge'));continue
        mmff,confid=min(good)
        xyz=Chem.MolToXYZBlock(mol,confId=confid)
        neutral=folder/'neutral';neutral.mkdir(exist_ok=True);(neutral/'input.xyz').write_text(xyz,encoding='utf-8')
        identity=dict(name=name,input_smiles=smiles,canonical_isomeric_smiles=Chem.MolToSmiles(base,canonical=True,isomericSmiles=True),formula=rdMolDescriptors.CalcMolFormula(base),formal_charge=0,neutral_multiplicity=1,atom_count=mol.GetNumAtoms(),heavy_atoms=base.GetNumAtoms(),rdkit_seed=31415,requested_conformers=4,returned_conformers=len(ids),selected_mmff_energy_kcal_mol=float(mmff),selected_conformer_id=confid,atom_order=[dict(index_1based=a.GetIdx()+1,element=a.GetSymbol(),neighbors_1based=[n.GetIdx()+1 for n in a.GetNeighbors()]) for a in mol.GetAtoms()])
        (folder/'identity.json').write_text(json.dumps(identity,indent=2),encoding='utf-8')
        nrun=run_xtb(executable,neutral,['input.xyz','--gfn','2','--alpb','acetonitrile','--chrg','0','--uhf','0','--opt','tight','--json'],args.timeout)
        if not(nrun['completed'] and nrun.get('optimization_converged') and (neutral/'xtbopt.xyz').exists()):
            failures.append(dict(molecule=name,stage='neutral',details=nrun));continue
        states={0:nrun}
        for charge,subname in [(1,'cation'),(-1,'anion')]:
            sub=folder/subname;sub.mkdir(exist_ok=True);shutil.copyfile(neutral/'xtbopt.xyz',sub/'input.xyz')
            states[charge]=run_xtb(executable,sub,['input.xyz','--gfn','2','--alpb','acetonitrile','--chrg',str(charge),'--uhf','1','--scc','--json'],args.timeout)
        if not all(s['completed'] for s in states.values()):
            failures.append(dict(molecule=name,stage='charged_singlepoints',details=states));continue
        q0=np.loadtxt(neutral/'charges');qp=np.loadtxt(folder/'cation'/'charges');qm=np.loadtxt(folder/'anion'/'charges')
        removal=qp-q0;addition=q0-qm
        for i,atom in enumerate(mol.GetAtoms()):
            atom_rows.append(dict(molecule=name,atom_index_1based=i+1,element=atom.GetSymbol(),q_neutral=float(q0[i]),q_cation=float(qp[i]),q_anion=float(qm[i]),electron_removal_response=float(removal[i]),electron_addition_response=float(addition[i])))
        neutral_text=(neutral/'stdout.log').read_text(encoding='utf-8')
        homo_lines=[x for x in neutral_text.splitlines() if '(HOMO)' in x]
        lumo_lines=[x for x in neutral_text.splitlines() if '(LUMO)' in x]
        def orbital(lines):
            if not lines:return None
            # Orbital table prints Hartree then eV; last numeric token is eV.
            values=re.findall(r'[-+]?\d+\.\d+',lines[-1]);return float(values[-1]) if values else None
        row=dict(molecule=name,smiles=identity['canonical_isomeric_smiles'],formula=identity['formula'],mol_weight_g_mol=Descriptors.MolWt(base),TPSA_A2=rdMolDescriptors.CalcTPSA(base),Crippen_logP=Descriptors.MolLogP(base),neutral_energy_Eh=states[0]['energy_Eh'],cation_energy_Eh=states[1]['energy_Eh'],anion_energy_Eh=states[-1]['energy_Eh'],fixed_nuclei_removal_gap_eV=(states[1]['energy_Eh']-states[0]['energy_Eh'])*HARTREE_EV,fixed_nuclei_addition_affinity_eV=(states[0]['energy_Eh']-states[-1]['energy_Eh'])*HARTREE_EV,HOMO_eV=orbital(homo_lines),LUMO_eV=orbital(lumo_lines),neutral_charge_sum=float(q0.sum()),cation_charge_sum=float(qp.sum()),anion_charge_sum=float(qm.sum()),removal_response_sum=float(removal.sum()),addition_response_sum=float(addition.sum()),geometry_sha256=hashlib.sha256((neutral/'xtbopt.xyz').read_bytes()).hexdigest())
        # Solvent sensitivity at exactly the same nuclei; no gas-phase reoptimization.
        gas={}
        for charge,subname in [(0,'gas_neutral'),(1,'gas_cation'),(-1,'gas_anion')]:
            sub=folder/subname;sub.mkdir(exist_ok=True);shutil.copyfile(neutral/'xtbopt.xyz',sub/'input.xyz')
            gas[charge]=run_xtb(executable,sub,['input.xyz','--gfn','2','--chrg',str(charge),'--uhf',str(abs(charge)),'--scc','--json'],args.timeout)
        if all(s['completed'] for s in gas.values()):
            row.update(gas_neutral_energy_Eh=gas[0]['energy_Eh'],gas_cation_energy_Eh=gas[1]['energy_Eh'],gas_anion_energy_Eh=gas[-1]['energy_Eh'],gas_fixed_nuclei_removal_gap_eV=(gas[1]['energy_Eh']-gas[0]['energy_Eh'])*HARTREE_EV,gas_fixed_nuclei_addition_affinity_eV=(gas[0]['energy_Eh']-gas[-1]['energy_Eh'])*HARTREE_EV)
            row['alpb_minus_gas_removal_gap_eV']=row['fixed_nuclei_removal_gap_eV']-row['gas_fixed_nuclei_removal_gap_eV']
        else:failures.append(dict(molecule=name,stage='gas_sensitivity',details=gas))
        rows.append(row)
        pd.DataFrame(rows).to_csv(out/'molecular_summary.csv',index=False)
        pd.DataFrame(atom_rows).to_csv(out/'atom_charge_responses.csv',index=False)
    result=dict(evidence='Executed real-molecule semiempirical GFN2-xTB and RDKit descriptors',timestamp_utc=datetime.now(timezone.utc).isoformat(),xtb_version=version_line,rdkit_version=rdBase.rdkitVersion,method='GFN2-xTB/ALPB(acetonitrile); neutral tight optimization; neutral geometry held fixed for cation/anion SCC',charge_multiplicity=dict(neutral=[0,1],cation=[1,2],anion=[-1,2]),threads=1,timeout_per_job_s=args.timeout,requested=len(MOLECULES[:args.limit]),completed=len(rows),rows=rows,failures=failures,limitations=['Semiempirical model, not DFT.','No Hessian or thermal correction; optimized structures not certified minima.','Each charge state uses its own equilibrium ALPB response; fixed-nuclei gaps are not rigorously nonequilibrium vertical ionization energies.','No calibrated reference electrode, experimental Eox, reaction barrier, or regioselectivity prediction.','Atom indices refer to saved identity.json and coordinates, not unverified IUPAC site labels.','Mulliken finite differences are model-dependent charge responses, not site-specific oxidation potentials.'])
    (out/'summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(dict(completed=len(rows),failures=len(failures)),indent=2),flush=True)
    return 0 if len(rows)==len(MOLECULES[:args.limit]) and not failures else 1

if __name__=='__main__':raise SystemExit(main())
