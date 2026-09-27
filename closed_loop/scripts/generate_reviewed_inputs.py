"""Write paired-state DFT input templates from existing, hash-verified xTB geometries.

Templates target Gaussian 16 and ORCA 6.1. No quantum engine is called here.
"""
from pathlib import Path
import argparse,hashlib,json,shutil
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from audit_closed_loop import ROOT,SUBSTRATE,save_json

REPO=ROOT.parent
MOLECULES={
    'target_indole':(SUBSTRATE,'production/results/target_xtb/neutral','production/results/audit/target_identity.json'),
    'indoline_control':('c1ccc2c(c1)CCN2','research/results/molecules/indoline/neutral','research/results/molecules/indoline/identity.json')}

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def read_verified_geometry(smiles,folder,identity_path):
    manifest=dict((name,h) for h,name in (line.split('  ',1) for line in (REPO/'provenance/SHA256SUMS.txt').read_text().splitlines()))
    paths=[folder+'/xtbopt.xyz',folder+'/run.json',identity_path]
    for p in paths:
        if digest(REPO/p)!=manifest[p]:raise ValueError('Prior calculation changed: '+p)
    run=json.loads((REPO/paths[1]).read_text());identity=json.loads((REPO/identity_path).read_text())
    if not run['completed']:raise ValueError('Prior geometry optimization not completed')
    mol=Chem.AddHs(Chem.MolFromSmiles(smiles))
    if Chem.MolToSmiles(Chem.RemoveHs(mol))!=identity['canonical_smiles']:raise ValueError('Molecular identity mismatch')
    lines=(REPO/paths[0]).read_text().splitlines();n=int(lines[0]);atoms=[l.split() for l in lines[2:] if l.strip()]
    if n!=len(atoms) or n!=mol.GetNumAtoms():raise ValueError('XYZ atom count mismatch')
    symbols=[a.GetSymbol() for a in mol.GetAtoms()]
    if symbols!=[a[0] for a in atoms] or symbols!=[a['element'] for a in identity['atom_order']]:raise ValueError('Atom order mismatch')
    xyz=np.asarray([[float(v) for v in a[1:4]] for a in atoms]);dist=np.linalg.norm(xyz[:,None]-xyz[None,:],axis=2);np.fill_diagonal(dist,np.inf)
    if not np.isfinite(xyz).all() or float(dist.min())<.6:raise ValueError('Invalid geometry or colliding nuclei')
    coords='\n'.join(f'{a[0]:<2} {float(a[1]): .10f} {float(a[2]): .10f} {float(a[3]): .10f}' for a in atoms)
    bonds=[float(dist[b.GetBeginAtomIdx(),b.GetEndAtomIdx()]) for b in mol.GetBonds()]
    info=dict(smiles=smiles,canonical_smiles=Chem.MolToSmiles(Chem.RemoveHs(mol)),formula=rdMolDescriptors.CalcMolFormula(mol),atoms=n,neutral_electrons=sum(a.GetAtomicNum() for a in mol.GetAtoms()),minimum_pair_distance_A=float(dist.min()),bond_distance_range_A=[min(bonds),max(bonds)],geometry_method='Previously executed neutral GFN2-xTB/ALPB(acetonitrile) optimization; reused, not newly run',source_files={p:manifest[p] for p in paths})
    return coords,info

def gaussian(base,coords,charge,multiplicity,nproc=2,memory_GB=3):
    prefix='u' if multiplicity==2 else ''
    return f'''%chk={base}.chk
%nprocshared={nproc}
%mem={memory_GB}GB
#p {prefix}B3LYP/def2SVP empiricaldispersion=gd3bj opt=(tight,maxcycles=200) freq scrf=(smd,solvent=acetonitrile) int=ultrafine scf=(tight,xqc,maxcycle=512)

{base}: proposed solution optimization and harmonic frequencies

{charge} {multiplicity}
{coords}

--Link1--
%chk={base}.chk
%nprocshared={nproc}
%mem={memory_GB}GB
#p {prefix}M062X/def2TZVP scrf=(smd,solvent=acetonitrile) geom=check guess=read int=ultrafine scf=(tight,xqc,maxcycle=512)

{base}: proposed single point on optimized geometry

{charge} {multiplicity}

'''

def orca(base,coords,charge,multiplicity,stage,nproc=2,maxcore_MB=1000):
    spin='UKS' if multiplicity==2 else 'RKS'
    method='B3LYP/G D3BJ def2-SVP Opt Freq' if stage=='opt' else 'M062X def2-TZVP'
    geometry=f'* xyz {charge} {multiplicity}\n{coords}\n*' if stage=='opt' else f'* xyzfile {charge} {multiplicity} {base}_opt.xyz'
    comment='# Requires successfully converged opt geometry and frequency inspection.\n' if stage=='sp' else ''
    return f'''# ORCA 6.1 input template; engine parsing and DFT results remain unverified.
{comment}! {spin} {method} RIJCOSX def2/J TightSCF DEFGRID3 CPCM
%pal nprocs {nproc} end
%maxcore {maxcore_MB}
%cpcm
  smd true
  SMDsolvent "acetonitrile"
end
%scf MaxIter 512 end
{geometry}
'''

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'inputs/reviewed');a=p.parse_args();out=a.output;out.mkdir(parents=True,exist_ok=True)
    manifest=dict(status='input_templates_generated_not_engine_validated',engines={'Gaussian':'16','ORCA':'6.1'},engine_executables_on_checked_PATH={n:bool(shutil.which(n)) for n in ['g16','orca']},dft_executed=False,profile=dict(nproc=2,Gaussian_memory_GB=3,ORCA_maxcore_MB_per_process=1000,ORCA_nominal_total_MB=2000,suggested_scheduler_memory_GB=4,walltime='unknown; measure pilot',CPU_only=True),molecules={},files=[],job_dependencies=[])
    for name,(smiles,folder,identity) in MOLECULES.items():
        coords,info=read_verified_geometry(smiles,folder,identity);manifest['molecules'][name]=info
        (out/f'{name}_starting.xyz').write_text(f'{info["atoms"]}\nReused neutral GFN2-xTB/ALPB(MeCN) geometry; not DFT\n{coords}\n',encoding='utf-8',newline='\n')
        for state,charge,mult in [('neutral',0,1),('cation',1,2)]:
            electrons=info['neutral_electrons']-charge
            if (electrons-(mult-1))%2:raise ValueError('Electron/spin parity mismatch')
            base=f'{name}_{state}'
            content={base+'.gjf':gaussian(base,coords,charge,mult),base+'_opt.inp':orca(base,coords,charge,mult,'opt'),base+'_sp.inp':orca(base,coords,charge,mult,'sp')}
            for filename,text in content.items():
                path=out/filename;path.write_text(text,encoding='ascii',newline='\n')
                manifest['files'].append(dict(path=filename,sha256=digest(path),molecule=name,state=state,charge=charge,multiplicity=mult,electrons=electrons))
            manifest['job_dependencies'].append(dict(sp_input=base+'_sp.inp',requires_output_from=base+'_opt.inp',required_xyz=base+'_opt.xyz',required_checks=['SCF convergence','Geometry convergence','No significant imaginary frequency','Spin and wavefunction assessment'],geometry_exists=False))
    save_json(out/'input_manifest.json',manifest)
    print(json.dumps(dict(input_files=len(manifest['files']),molecules=len(manifest['molecules']),dft_executed=False,PATH_detection=manifest['engine_executables_on_checked_PATH'])))
if __name__=='__main__':main()
