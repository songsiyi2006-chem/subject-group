"""5T35 experimental-structure audit and bounded explicit-solvent WTMetaD pilot.

The default run is an atomistic ENGINE PILOT, never a converged PMF or alpha.
All heavy simulation dependencies are imported only when preparing/running MD.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform
import shutil
import time
from datetime import datetime, timezone

import numpy as np

HERE = Path(__file__).resolve().parent
EXECUTED_SOURCE = Path(__file__).read_bytes()
CONFIG = dict(seed=190921, temperature_K=300., timestep_ps=.001,
              equilibration_steps=200, metadynamics_steps=1000,
              hill_frequency_steps=50, initial_hill_kJ_mol=.25,
              bias_factor=10., distance_sigma_nm=.05, angle_sigma_rad=.15,
              solvent_padding_nm=1., ionic_strength_M=.15,
              cutoff_nm=.9, cpu_threads=1,
              protein_forcefield='amber14-all.xml', water_forcefield='amber14/tip3p.xml',
              ligand_forcefield='openff-2.3.0.offxml',
              ligand_charge_method='openff-gnn-am1bcc-1.0.0.pt',
              distance_min_nm=1.5, distance_max_nm=8.,
              distance_grid=96, angle_grid=96,
              minimization_max_iterations=50)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def table(path, rows):
    with Path(path).open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def structure_atoms(path):
    atoms = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(('ATOM  ', 'HETATM')) and line[16] in (' ', 'A'):
            atoms.append(dict(record=line[:6].strip(), chain=line[21],
                              resname=line[17:20].strip(), resid=int(line[22:26]),
                              atom=line[12:16].strip(), element=line[76:78].strip(),
                              xyz_A=[float(line[30:38]), float(line[38:46]), float(line[46:54])]))
    return atoms


def orientation_angle(points):
    """Signed four-centroid dihedral in radians; undefined axes are rejected."""
    p = np.asarray(points, float)
    b0, b1, b2 = -(p[1]-p[0]), p[2]-p[1], p[3]-p[2]
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0-np.dot(b0,b1)*b1, b2-np.dot(b2,b1)*b1
    if min(np.linalg.norm(v), np.linalg.norm(w)) < 1e-10:
        raise ValueError('Undefined orientation: collinear centroid axes')
    return float(np.arctan2(np.dot(np.cross(b1,v), w), np.dot(v,w)))


def structure_audit(out):
    atoms=structure_atoms(HERE/'inputs/5T35.pdb')
    v=[a for a in atoms if a['record']=='ATOM' and a['chain']=='D']
    t=[a for a in atoms if a['record']=='ATOM' and a['chain']=='A']
    ligand=[a for a in atoms if a['resname']=='759' and a['chain']=='D']
    xyzv,xyzt=np.array([a['xyz_A'] for a in v]),np.array([a['xyz_A'] for a in t])
    d=np.linalg.norm(xyzv[:,None,:]-xyzt[None,:,:],axis=2)
    pairs=[]
    residues={}
    for i,j in zip(*np.where(d<=4.5)):
        a,b=v[i],t[j]
        pairs.append(dict(VHL_residue=f"{a['resname']}{a['resid']}",VHL_atom=a['atom'],
                          BRD4_residue=f"{b['resname']}{b['resid']}",BRD4_atom=b['atom'],distance_A=float(d[i,j])))
        key=(a['resid'],a['resname'],b['resid'],b['resname'])
        residues[key]=min(residues.get(key,1e10),float(d[i,j]))
    table(out/'interface_atom_contacts.csv',sorted(pairs,key=lambda r:r['distance_A']))
    table(out/'interface_residue_contacts.csv',[
        dict(VHL_resid=k[0],VHL_resname=k[1],BRD4_resid=k[2],BRD4_resname=k[3],minimum_distance_A=val,
             mutation_free_energy_status='not_computed') for k,val in sorted(residues.items())])
    summary=dict(structure='5T35',experimental_method='X-ray diffraction',resolution_A=2.7,
                 assembly='author chains A/B/C/D and ligand 759 author chain D',
                 chain_identity={'A':'BRD4 BD2','B':'Elongin B','C':'Elongin C','D':'VHL'},
                 ligand='MZ1; CCD 759',ligand_heavy_atoms=len(ligand),
                 VHL_heavy_atoms=len(v),BRD4_heavy_atoms=len(t),
                 PPI_atom_pairs_within_4_5_A=len(pairs),PPI_residue_pairs_within_4_5_A=len(residues),
                 closest_PPI_distance_A=float(d.min()),
                 source_sha256={p.name:sha(p) for p in (HERE/'inputs').glob('*') if p.is_file()})
    dump(out/'structure_audit.json',summary)
    np.savez_compressed(out/'crystal_interface.npz',VHL_A=xyzv,BRD4_A=xyzt,
                        MZ1_A=np.array([a['xyz_A'] for a in ligand]))
    return summary


def prepare_ligand(out,c):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from openff.toolkit import Molecule, ForceField
    from openmm import unit
    template=Chem.MolFromMolFile(str(HERE/'inputs/759_ideal.sdf'),removeHs=True)
    lines=[s for s in (HERE/'inputs/5T35.pdb').read_text().splitlines()
           if s.startswith('HETATM') and s[17:20]=='759' and s[21]=='D']
    raw=Chem.MolFromPDBBlock('\n'.join(lines)+'\nEND\n',removeHs=True,sanitize=False)
    mol=AllChem.AssignBondOrdersFromTemplate(template,raw)
    Chem.AssignStereochemistryFrom3D(mol)
    a,b=Chem.MolToSmiles(mol),Chem.MolToSmiles(template)
    if a!=b:
        raise ValueError(f'Crystal vs CCD stereochemistry mismatch: {a} != {b}')
    mol=Chem.AddHs(mol,addCoords=True)
    writer=Chem.SDWriter(str(out/'MZ1_prepared.sdf'));writer.write(mol);writer.close()
    offmol=Molecule.from_rdkit(mol,allow_undefined_stereo=False)
    offmol.assign_partial_charges(c['ligand_charge_method'])
    ff=ForceField(c['ligand_forcefield'])
    system=ff.create_openmm_system(offmol.to_topology(),charge_from_molecules=[offmol])
    top=offmol.to_topology().to_openmm()
    for chain in top.chains(): chain.id='L'
    for residue in top.residues(): residue.name='MZ1';residue.id='1'
    positions=np.asarray(mol.GetConformer().GetPositions())*.1*unit.nanometer
    charges=offmol.partial_charges.m_as('elementary_charge')
    if abs(charges.sum())>1e-6:
        raise ValueError('Current neutralization protocol only supports neutral MZ1')
    dump(out/'ligand_parameterization.json',dict(atoms=offmol.n_atoms,canonical_isomeric_smiles=a,
             crystal_stereochemistry_matches_CCD=True,charge_sum_e=float(charges.sum()),
             forcefield=c['ligand_forcefield'],charges=c['ligand_charge_method'],
             charge_e=charges.tolist(),charge_status='NAGL learned AM1-BCC surrogate; not fresh quantum AM1-BCC'))
    return system,top,positions


def merge_systems(protein,ligand):
    """Merge disjoint standard AMBER/OpenFF systems with full cross LJ/Coulomb terms."""
    import openmm as mm
    offset=protein.getNumParticles()
    for i in range(ligand.getNumParticles()): protein.addParticle(ligand.getParticleMass(i))
    for i in range(ligand.getNumConstraints()):
        a,b,d=ligand.getConstraintParameters(i);protein.addConstraint(a+offset,b+offset,d)
    pnb=next(f for f in protein.getForces() if isinstance(f,mm.NonbondedForce))
    supported=(mm.HarmonicBondForce,mm.HarmonicAngleForce,mm.PeriodicTorsionForce,mm.NonbondedForce,mm.CMMotionRemover)
    if any(ligand.isVirtualSite(i) for i in range(ligand.getNumParticles())):
        raise ValueError('Virtual sites require an explicit merge implementation')
    for force in ligand.getForces():
        if not isinstance(force,supported): raise ValueError(f'Unsupported force merge: {type(force)}')
        if isinstance(force,mm.NonbondedForce):
            if force.getNumParticleParameterOffsets() or force.getNumExceptionParameterOffsets():
                raise ValueError('Offsets need explicit merge implementation')
            for i in range(force.getNumParticles()): pnb.addParticle(*force.getParticleParameters(i))
            for i in range(force.getNumExceptions()):
                a,b,q,s,e=force.getExceptionParameters(i);pnb.addException(a+offset,b+offset,q,s,e)
        elif isinstance(force,mm.HarmonicBondForce):
            new=mm.HarmonicBondForce()
            for i in range(force.getNumBonds()):
                a,b,r,k=force.getBondParameters(i);new.addBond(a+offset,b+offset,r,k)
            protein.addForce(new)
        elif isinstance(force,mm.HarmonicAngleForce):
            new=mm.HarmonicAngleForce()
            for i in range(force.getNumAngles()):
                a,b,c,t,k=force.getAngleParameters(i);new.addAngle(a+offset,b+offset,c+offset,t,k)
            protein.addForce(new)
        elif isinstance(force,mm.PeriodicTorsionForce):
            new=mm.PeriodicTorsionForce()
            for i in range(force.getNumTorsions()):
                a,b,c,d,n,t,k=force.getTorsionParameters(i);new.addTorsion(a+offset,b+offset,c+offset,d+offset,n,t,k)
            protein.addForce(new)
    return protein


def prepare_system(out,c,ignore_cache=False):
    import openmm as mm
    from openmm import app,unit
    from scipy.spatial import cKDTree
    import random
    random.seed(c['seed']);np.random.seed(c['seed'])
    cached=HERE/'inputs/prepared'
    if (cached/'manifest.json').exists() and not ignore_cache:
        for name,digest in json.loads((cached/'manifest.json').read_text())['files_sha256'].items():
            if sha(cached/name)!=digest: raise ValueError(f'Prepared input hash mismatch: {name}')
        system=mm.XmlSerializer.deserialize(gzip.decompress((cached/'unbiased_system.xml.gz').read_bytes()).decode())
        pdb=app.PDBFile(io.StringIO(gzip.decompress((cached/'prepared_explicit.pdb.gz').read_bytes()).decode()))
        mod=app.Modeller(pdb.topology,pdb.positions)
        for name in ('ligand_parameterization.json','preparation.json','MZ1_prepared.sdf'):
            (out/name).write_bytes((cached/name).read_bytes())
        dump(out/'preparation_reuse.json',dict(mode='verified precomputed atomistic parameterization',
            source='inputs/prepared',manifest_sha256=sha(cached/'manifest.json'),
            coordinates_precision='PDB restart, 0.001 angstrom',
            configuration_compatibility='Cached preparation is fixed to the archived force fields, box and ion parameters.'))
        return system,mod
    from pdbfixer import PDBFixer
    ligand,ltop,lpos=prepare_ligand(out,c)
    selected=[s for s in (HERE/'inputs/5T35.pdb').read_text().splitlines()
              if (s.startswith(('ATOM  ','TER   ')) and s[21] in 'ABCD') or (s.startswith('SEQRES') and s[11] in 'ABCD')]
    fixer=PDBFixer(pdbfile=io.StringIO('\n'.join(selected)+'\nEND\n'))
    fixer.findMissingResidues()
    missing=dict(fixer.missingResidues)
    chains=list(fixer.topology.chains())
    # Unresolved terminal tails remain absent; internal gaps are explicitly modeled.
    for key in list(fixer.missingResidues):
        if key[1] in (0,len(list(chains[key[0]].residues()))): del fixer.missingResidues[key]
    fixer.findMissingAtoms();fixer.addMissingAtoms(seed=c['seed'])
    fixer.addMissingHydrogens(pH=7.)
    ff=app.ForceField(c['protein_forcefield'],c['water_forcefield'])
    mod=app.Modeller(fixer.topology,fixer.positions)
    mod.addSolvent(ff,model='tip3p',padding=c['solvent_padding_nm']*unit.nanometer,
                   boxShape='dodecahedron',ionicStrength=c['ionic_strength_M']*unit.molar)
    # Solvate the protein then remove any solvent/ions overlapping the real ligand.
    coords=np.array(mod.positions.value_in_unit(unit.nanometer))
    lcoords=np.array(lpos.value_in_unit(unit.nanometer))
    tree=cKDTree(lcoords)
    deletion=[]
    for res in mod.topology.residues():
        if res.name in ('HOH','NA','CL'):
            if min(tree.query(coords[[a.index for a in res.atoms()]])[0])<.30: deletion.append(res)
    if any(r.name!='HOH' for r in deletion):
        raise ValueError('A solvated ion overlaps MZ1; regenerate solvent rather than silently change net charge')
    mod.delete(deletion)
    psystem=ff.createSystem(mod.topology,nonbondedMethod=app.PME,nonbondedCutoff=c['cutoff_nm']*unit.nanometer,
                            constraints=app.HBonds,rigidWater=True,ewaldErrorTolerance=.0005)
    system=merge_systems(psystem,ligand)
    mod.add(ltop,lpos)
    app.PDBFile.writeFile(mod.topology,mod.positions,(out/'prepared_explicit.pdb').open('w'),keepIds=True)
    (out/'unbiased_system.xml').write_text(mm.XmlSerializer.serialize(system))
    dump(out/'preparation.json',dict(particles=system.getNumParticles(),constraints=system.getNumConstraints(),
         solvent_removed_overlapping_ligand=len(deletion),missing_residues_detected={str(k):v for k,v in missing.items()},
         internal_missing_residues_modeled={str(k):v for k,v in fixer.missingResidues.items()},
         hydrogens='PDBFixer pH 7 default variants; not a pKa calculation',
         protein_construct='Resolved A/B/C/D with missing atoms repaired; terminal tails omitted',
         periodic_box_nm=np.array(mod.topology.getPeriodicBoxVectors().value_in_unit(unit.nanometer)).tolist(),
         neutralization='Protein-based neutralization with 0.15 M NaCl; MZ1 net charge verified zero',
         missing_atom_reconstruction='PDBFixer modeled geometry, not experimental observation'))
    return system,mod


def cv_forces(system,topology):
    import openmm as mm
    atoms=list(topology.atoms())
    groups={ch:[a.index for a in atoms if a.residue.chain.id==ch and a.element.symbol!='H'] for ch in ('D','A')}
    d=mm.CustomCentroidBondForce(2,'distance(g1,g2)')
    for ch in ('D','A'): d.addGroup(groups[ch])
    d.addBond([0,1],[])
    # Nonperiodic centroids of whole solute; no wrapping is used in this short pilot.
    d.setUsesPeriodicBoundaryConditions(False)
    orient=mm.CustomCentroidBondForce(4,'dihedral(g1,g2,g3,g4)')
    anchor_groups=[]
    for ch,low,high in [('D',65,80),('D',120,135),('A',350,365),('A',430,445)]:
        indices=[a.index for a in atoms if a.residue.chain.id==ch and a.name=='CA' and low<=int(a.residue.id)<=high]
        if len(indices)<3: raise ValueError(f'Orientation anchor too small: {ch}:{low}-{high}')
        orient.addGroup(indices);anchor_groups.append(indices)
    orient.addBond([0,1,2,3],[]);orient.setUsesPeriodicBoundaryConditions(False)
    return d,orient,groups,anchor_groups


def make_bias_force(c,groups,anchors):
    """Analytic mixed-periodicity bias, with real atomic centroid derivatives."""
    import openmm as mm
    if c['bias_factor']<=1 or c['hill_frequency_steps']<=0:
        raise ValueError('Invalid WTMetaD bias factor or deposition frequency')
    if c['metadynamics_steps']%c['hill_frequency_steps']:
        raise ValueError('Biased steps must contain whole deposition intervals')
    nhills=c['metadynamics_steps']//c['hill_frequency_steps']
    if not 1<=nhills<=100:
        raise ValueError('Analytic local engine supports 1..100 hills; production needs a scalable bias representation')
    expressions=[f'h{i}*exp(-0.5*((d-r{i})/{c["distance_sigma_nm"]})^2-0.5*(atan2(sin(phi-p{i}),cos(phi-p{i}))/{c["angle_sigma_rad"]})^2)' for i in range(nhills)]
    bias=mm.CustomCentroidBondForce(6,'+'.join(expressions)+';d=distance(g1,g2);phi=dihedral(g3,g4,g5,g6)')
    for g in [groups['D'],groups['A']]+anchors: bias.addGroup(g)
    bias.addBond(list(range(6)),[]);bias.setUsesPeriodicBoundaryConditions(False)
    for i in range(nhills):
        for name,val in [(f'h{i}',0.),(f'r{i}',0.),(f'p{i}',0.)]: bias.addGlobalParameter(name,val)
    bias.setForceGroup(31)
    return bias


def sync_bias_archive(out):
    """Audit a serialized bias against measured hill history; no trajectory edits."""
    import openmm as mm
    path=Path(out)/'biased_system.xml.gz'
    before=sha(path)
    system=mm.XmlSerializer.deserialize(gzip.decompress(path.read_bytes()).decode())
    bias=next(f for f in system.getForces() if isinstance(f,mm.CustomCentroidBondForce) and f.getForceGroup()==31)
    lookup={bias.getGlobalParameterName(i):i for i in range(bias.getNumGlobalParameters())}
    values={}
    for row in csv.DictReader(io.StringIO((Path(out)/'hill_history.csv').read_text(encoding='utf-8'))):
        n=int(row['hill'])-1
        values.update({f'h{n}':float(row['height_kJ_mol']),f'r{n}':float(row['center_distance_nm']),f'p{n}':float(row['center_angle_rad'])})
    changed=0
    for name,value in values.items():
        i=lookup[name]
        if bias.getGlobalParameterDefaultValue(i)!=value:changed+=1
        bias.setGlobalParameterDefaultValue(i,value)
    if changed:path.write_bytes(gzip.compress(mm.XmlSerializer.serialize(system).encode(),mtime=0))
    result=dict(parameters_verified=len(values),default_values_synchronized=changed,
        prior_system_sha256=before,final_system_sha256=sha(path),hill_history_sha256=sha(Path(out)/'hill_history.csv'),
        operation='Archive scalar bias defaults are synchronized to the recorded live-Context hills; trajectories and energies are unchanged.',
        postprocessing_driver_sha256=sha(__file__))
    dump(Path(out)/'bias_archive_audit.json',result)
    return result


def save_checkpoint(out,sim,c,stage):
    """Publish a hash-verified restart only after numerical files are on disk."""
    path=Path(out)/'latest.checkpoint.gz'
    temporary=path.with_suffix('.tmp')
    temporary.write_bytes(gzip.compress(sim.context.createCheckpoint(),compresslevel=1,mtime=0))
    os.replace(temporary,path)
    candidates=['latest.checkpoint.gz','minimization.csv','minimization_energies.json','thermal_reference.npz',
                'atomistic_cv_trajectory.csv','hill_history.csv','solute_trajectory.npz','unconverged_bias_free_energy.npz']
    metadata=dict(stage=stage,step=sim.currentStep,
        config={k:v for k,v in c.items() if k!='restart_from'},
        prepared_manifest_sha256=sha(HERE/'inputs/prepared/manifest.json'),
        files_sha256={name:sha(Path(out)/name) for name in candidates if (Path(out)/name).exists()},
        updated_utc=datetime.now(timezone.utc).isoformat(),
        portability='OpenMM binary checkpoint: same OpenMM version, CPU platform and compatible hardware required')
    tempmeta=Path(out)/'checkpoint_metadata.tmp'
    dump(tempmeta,metadata);os.replace(tempmeta,Path(out)/'checkpoint_metadata.json')


def load_checkpoint(out,sim,c):
    source=Path(c['restart_from'])
    metadata=json.loads((source/'checkpoint_metadata.json').read_text(encoding='utf-8'))
    if metadata['config']!={k:v for k,v in c.items() if k!='restart_from'}:
        raise ValueError('Restart configuration differs from the recorded simulation')
    if metadata['prepared_manifest_sha256']!=sha(HERE/'inputs/prepared/manifest.json'):
        raise ValueError('Restart atomistic preparation differs')
    for name,digest in metadata['files_sha256'].items():
        if sha(source/name)!=digest:raise ValueError(f'Incomplete/interrupted checkpoint transaction: {name}')
    sim.context.loadCheckpoint(gzip.decompress((source/'latest.checkpoint.gz').read_bytes()))
    if sim.currentStep!=metadata['step']:raise ValueError('Checkpoint step does not match its metadata')
    for name in metadata['files_sha256']:shutil.copyfile(source/name,Path(out)/name)
    dump(Path(out)/'restart_provenance.json',dict(source_directory=str(source.resolve()),
        checkpoint_metadata_sha256=sha(source/'checkpoint_metadata.json'),resumed_step=sim.currentStep))
    return metadata


def atomistic_pilot(out,c):
    import openmm as mm
    from openmm import app,unit
    np.random.seed(c['seed'])
    import random
    random.seed(c['seed'])
    system,mod=prepare_system(out,c)
    distance,angle,groups,anchors=cv_forces(system,mod.topology)
    # OpenMM's tabulated MetaD class rejects mixed periodic/nonperiodic CVs.
    # Exact analytic Gaussian hills support nonperiodic d and wrapped periodic phi.
    nhills=c['metadynamics_steps']//c['hill_frequency_steps']
    # One centroid force avoids creating a second all-atom Context for CV evaluation.
    bias=make_bias_force(c,groups,anchors);system.addForce(bias)
    integrator=mm.LangevinMiddleIntegrator(c['temperature_K']*unit.kelvin,1./unit.picosecond,c['timestep_ps']*unit.picosecond)
    integrator.setRandomNumberSeed(c['seed'])
    sim=app.Simulation(mod.topology,system,integrator,mm.Platform.getPlatformByName('CPU'),{'Threads':str(c['cpu_threads'])})
    print('Task19: explicit all-atom CPU context created',flush=True)
    if c.get('restart_from'):
        load_checkpoint(out,sim,c)
        energies=json.loads((out/'minimization_energies.json').read_text())
        initial,minimized=energies['initial_energy_kJ_mol'],energies['minimized_energy_kJ_mol']
        print(f'Task19: resumed verified checkpoint at step {sim.currentStep}',flush=True)
    else:
        sim.context.setPositions(mod.positions)
        initial=float(sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
        minimization_rows=[]
        class Reporter(mm.MinimizationReporter):
            def report(self,iteration,x,grad,args):
                minimization_rows.append(dict(iteration=iteration,**{str(k).replace(' ','_'):float(v) for k,v in args.items()}))
                if len(minimization_rows)%25==0:
                    print(f'Task19: minimizer report {len(minimization_rows)}',flush=True)
                return False
        sim.minimizeEnergy(maxIterations=c['minimization_max_iterations'],reporter=Reporter())
        if minimization_rows:table(out/'minimization.csv',minimization_rows)
        print('Task19: energy minimization finished',flush=True)
        minimized=float(sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
        dump(out/'minimization_energies.json',dict(initial_energy_kJ_mol=initial,minimized_energy_kJ_mol=minimized,
             minimization_callbacks=len(minimization_rows),max_iterations_per_constraint_stage=c['minimization_max_iterations']))
        save_checkpoint(out,sim,c,'minimized')
        shutil.copyfile(out/'latest.checkpoint.gz',out/'minimized.checkpoint.gz')
    if sim.currentStep<c['equilibration_steps']:
        if sim.currentStep==0:sim.context.setVelocitiesToTemperature(c['temperature_K']*unit.kelvin,c['seed'])
        sim.step(c['equilibration_steps']-sim.currentStep)
        print('Task19: short thermalization finished',flush=True)
    rows=[];positions=[];bias_arrays=[];hills=[]
    dgrid=np.linspace(c['distance_min_nm'],c['distance_max_nm'],c['distance_grid'])
    pgrid=np.linspace(-math.pi,math.pi,c['angle_grid'])
    dd,pp=np.meshgrid(dgrid,pgrid);bias_grid=np.zeros(dd.shape)
    minbox=min(np.linalg.norm(np.array(mod.topology.getPeriodicBoxVectors().value_in_unit(unit.nanometer)),axis=1))
    allsolute=[a.index for a in mod.topology.atoms() if a.residue.name not in ('HOH','NA','CL')]
    if (out/'thermal_reference.npz').exists():reference=np.load(out/'thermal_reference.npz')['positions_nm']
    else:
        reference=np.array(sim.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer))[allsolute]
        np.savez_compressed(out/'thermal_reference.npz',positions_nm=reference)
    if sim.currentStep==c['equilibration_steps']:save_checkpoint(out,sim,c,'thermalized')
    if (out/'hill_history.csv').exists():
        hills=[{k:(int(v) if k=='hill' else float(v)) for k,v in r.items()} for r in csv.DictReader(io.StringIO((out/'hill_history.csv').read_text()))]
        rows=[{k:(int(v) if k in ('step','hills_deposited') else float(v)) for k,v in r.items()} for r in csv.DictReader(io.StringIO((out/'atomistic_cv_trajectory.csv').read_text()))]
        positions=list(np.load(out/'solute_trajectory.npz')['positions_nm'])
        bias_arrays=list(np.load(out/'unconverged_bias_free_energy.npz')['bias_derived_free_energy_kJ_mol'])
        bias_grid=-bias_arrays[-1]*(c['bias_factor']-1)/c['bias_factor']
    if sim.currentStep!=c['equilibration_steps']+len(hills)*c['hill_frequency_steps']:
        raise ValueError('Checkpoint time and completed hill history disagree')
    masses=np.array([system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(system.getNumParticles())])
    for n in range(len(hills),c['metadynamics_steps']//c['hill_frequency_steps']):
        sim.step(c['hill_frequency_steps'])
        state=sim.context.getState(getEnergy=True,getPositions=True,enforcePeriodicBox=False)
        xyz=np.array(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
        energy=float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
        centroids=[np.average(xyz[g],axis=0,weights=masses[g]) for g in [groups['D'],groups['A']]+anchors]
        cv=[float(np.linalg.norm(centroids[0]-centroids[1])),orientation_angle(centroids[2:])]
        oldbias=float(sim.context.getState(getEnergy=True,groups={31}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
        height=c['initial_hill_kJ_mol']*math.exp(-oldbias/(.008314462618*c['temperature_K']*(c['bias_factor']-1)))
        sim.context.setParameter(f'r{n}',cv[0]);sim.context.setParameter(f'p{n}',cv[1]);sim.context.setParameter(f'h{n}',height)
        # Context changes do not update System XML defaults. Preserve the actual
        # deposited bias in the serialized artifact as well as the live Context.
        for index,value in [(3*n,height),(3*n+1,cv[0]),(3*n+2,cv[1])]:
            bias.setGlobalParameterDefaultValue(index,value)
        hills.append(dict(hill=n+1,time_ps=sim.currentStep*c['timestep_ps'],center_distance_nm=float(cv[0]),center_angle_rad=float(cv[1]),height_kJ_mol=height,bias_before_deposition_kJ_mol=oldbias))
        angular=np.arctan2(np.sin(pp-cv[1]),np.cos(pp-cv[1]))
        bias_grid+=height*np.exp(-.5*((dd-cv[0])/c['distance_sigma_nm'])**2-.5*(angular/c['angle_sigma_rad'])**2)
        if not (np.isfinite(xyz).all() and np.isfinite(energy)): raise ValueError('Nonfinite atomistic state')
        if np.max(np.linalg.norm(xyz[allsolute]-reference,axis=1))>minbox/4:
            raise ValueError('Whole-solute Cartesian CV validity gate failed: motion approaches box scale')
        if not c['distance_min_nm']<cv[0]<c['distance_max_nm']: raise ValueError('Distance CV left grid')
        rows.append(dict(step=sim.currentStep,time_ps=sim.currentStep*c['timestep_ps'],distance_nm=float(cv[0]),
                         orientation_rad=float(cv[1]),potential_with_bias_kJ_mol=energy,hills_deposited=n+1))
        positions.append(xyz[allsolute].astype(np.float32))
        bias_arrays.append(-c['bias_factor']/(c['bias_factor']-1)*bias_grid.copy())
        table(out/'atomistic_cv_trajectory.csv',rows);table(out/'hill_history.csv',hills)
        np.savez_compressed(out/'solute_trajectory.npz',positions_nm=np.array(positions),indices=allsolute)
        np.savez_compressed(out/'unconverged_bias_free_energy.npz',bias_derived_free_energy_kJ_mol=np.array(bias_arrays),distance_nm=dgrid,orientation_rad=pgrid)
        save_checkpoint(out,sim,c,'biased')
        print(f'Task19: hill {n+1}/{nhills}; d={cv[0]:.6f} nm; phi={cv[1]:.6f} rad',flush=True)
    table(out/'atomistic_cv_trajectory.csv',rows)
    table(out/'hill_history.csv',hills)
    np.savez_compressed(out/'solute_trajectory.npz',positions_nm=np.array(positions),indices=allsolute)
    np.savez_compressed(out/'unconverged_bias_free_energy.npz',bias_derived_free_energy_kJ_mol=np.array(bias_arrays),
                        distance_nm=dgrid,orientation_rad=pgrid)
    state=sim.context.getState(getPositions=True)
    pdbtext=io.StringIO();app.PDBFile.writeFile(mod.topology,state.getPositions(),pdbtext,keepIds=True)
    (out/'final_explicit.pdb.gz').write_bytes(gzip.compress(pdbtext.getvalue().encode(),mtime=0))
    (out/'biased_system.xml.gz').write_bytes(gzip.compress(mm.XmlSerializer.serialize(system).encode(),mtime=0))
    sync_bias_archive(out)
    (out/'integrator.xml').write_text(mm.XmlSerializer.serialize(integrator))
    dump(out/'cv_definition.json',dict(distance='VHL D vs BRD4 A heavy-atom mass-weighted COM; nm',
         orientation='signed dihedral of D:65-80 CA / D:120-135 CA / A:350-365 CA / A:430-445 CA mass centroids; rad',
         distance_atom_indices=groups,orientation_anchor_indices=anchors,periodic_angle=True,
         spatial_PBC='PME for interactions; CV uses unwrapped continuous whole-solute Cartesian positions',
         validity='Short pilot only: fail if any solute atom moves > minimum box-vector length / 4 after equilibration. Production needs molecule-whole imaging and explicit dissociation/PBC treatment.'))
    return dict(status='completed_short_all_atom_explicit_solvent_engine_pilot',particles=system.getNumParticles(),
       equilibration_ps=c['equilibration_steps']*c['timestep_ps'],biased_time_ps=c['metadynamics_steps']*c['timestep_ps'],
       hills=len(rows),initial_energy_kJ_mol=initial,minimized_energy_kJ_mol=minimized,
       distance_range_nm=[min(r['distance_nm'] for r in rows),max(r['distance_nm'] for r in rows)],
       orientation_range_rad=[min(r['orientation_rad'] for r in rows),max(r['orientation_rad'] for r in rows)],
       finite_positions_and_energies=True,pmf_converged=False,cooperativity_alpha=None,cooperativity_delta_delta_G_kJ_mol=None,
       free_energy_status='unconverged bias diagnostic only; no basin populations or binding free energy',
       forcefield='Amber14 protein / TIP3P / Sage 2.3 MZ1 with NAGL AM1-BCC surrogate charges',openmm_version=mm.__version__)


def figure(out,s):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(2,2,figsize=(13,10),constrained_layout=True)
    xyz=np.load(out/'crystal_interface.npz')
    for key,col in [('VHL_A','#1a7399'),('BRD4_A','#cd684a'),('MZ1_A','#7363af')]:
        a=xyz[key];ax[0,0].scatter(a[:,0],a[:,1],s=3 if key!='MZ1_A' else 13,c=col,alpha=.5,label=key.replace('_A',''))
    ax[0,0].set(xlabel='Crystal x (angstrom)',ylabel='Crystal y (angstrom)',title='A  Experimental 5T35 assembly: real coordinates');ax[0,0].legend()
    pilot=s['atomistic']
    if (out/'atomistic_cv_trajectory.csv').exists():
        rows=list(csv.DictReader(io.StringIO((out/'atomistic_cv_trajectory.csv').read_text(encoding='utf-8'))))
        t=np.array([float(r['time_ps']) for r in rows]);d=np.array([float(r['distance_nm']) for r in rows]);p=np.array([float(r['orientation_rad']) for r in rows])
        ax[0,1].plot(t,d,color='#1a7399');tw=ax[0,1].twinx();tw.plot(t,p,color='#cd684a');tw.set_ylabel('Orientation (rad)',color='#cd684a')
        ax[0,1].set(xlabel='Total simulation time (ps)',ylabel='COM distance (nm)',title='B  Atomistic WTMetaD: short engine pilot')
        b=np.load(out/'unconverged_bias_free_energy.npz');z=b['bias_derived_free_energy_kJ_mol'][-1];z=z-z.min()
        im=ax[1,0].pcolormesh(b['distance_nm'],b['orientation_rad'],z,shading='auto',cmap='viridis')
        ax[1,0].scatter(d,p,c='white',s=10,edgecolor='black',linewidth=.3)
        ax[1,0].set(xlabel='COM distance (nm)',ylabel='Orientation (rad)',title='C  Unconverged bias diagnostic: NOT a PMF')
        fig.colorbar(im,ax=ax[1,0],label='-gamma/(gamma-1) V + constant (kJ/mol)')
    else:
        for a in (ax[0,1],ax[1,0]):a.axis('off');a.text(.1,.6,'All-atom MD not completed.\nSee atomistic_preflight_error.json.',transform=a.transAxes)
    ax[1,1].axis('off')
    text=('D  Evidence gates\n\n'+f"Real PPI contacts: {s['structure']['PPI_atom_pairs_within_4_5_A']} atom pairs <= 4.5 A\n"+
         f"All-atom status: {'short engine pilot complete' if pilot['status'].startswith('completed_') else 'not completed'}\n"+f"Explicit atoms: {pilot.get('particles','unavailable')}\nBiased duration: {pilot.get('biased_time_ps',0):g} ps\n\n"+
         'Converged PMF: NO\nCooperativity alpha: NOT COMPUTED\nProductive degradation basins: NOT IDENTIFIED\n\nBinary and ternary standard-state binding cycles,\nreplicates and convergence are still required.')
    ax[1,1].text(0,1,text,va='top',fontsize=10,linespacing=1.7,wrap=True)
    fig.suptitle('Task 19 | Explicit atomistic computation with a strict free-energy gate',fontsize=16)
    (out/'figures').mkdir(exist_ok=True)
    fig.savefig(out/'figures/fig19_ternary_metadynamics_pmf_landscape.png',dpi=300)
    plt.close(fig)


def run(out:Path,restart_from:Path|None=None)->dict:
    out=Path(out)
    if out.exists() and any(out.iterdir()): raise FileExistsError(f'Refusing nonempty output: {out}')
    out.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();stamp=datetime.now(timezone.utc).isoformat()
    (out/'executed_driver.py').write_bytes(EXECUTED_SOURCE)
    c=dict(CONFIG)
    if restart_from is not None:c['restart_from']=str(Path(restart_from).resolve())
    dump(out/'config.json',c)
    structural=structure_audit(out)
    try:
        atomistic=atomistic_pilot(out,c)
    except Exception as exc:
        import traceback
        dump(out/'atomistic_preflight_error.json',dict(error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()))
        atomistic=dict(status='not_completed_preflight_or_execution_failure',cooperativity_alpha=None,
                       cooperativity_delta_delta_G_kJ_mol=None,pmf_converged=False)
    summary=dict(task=19,evidence_class='experimental structure plus bounded all-atom simulation pilot',structure=structural,
                 atomistic=atomistic,publication_claims_supported=False)
    dump(out/'summary.json',summary)
    figure(out,summary)
    dump(out/'verification.json',dict(structure_identity_checked=True,ligand_present=structural['ligand_heavy_atoms']>0,
          no_cooperativity_overclaim=atomistic['cooperativity_alpha'] is None,
          scientific_acceptance='PMF/cooperativity/productive-geometry gates not passed',
          atomistic_engine_passed=atomistic['status'].startswith('completed_')))
    dump(out/'run_metadata.json',dict(started_utc=stamp,elapsed_s=time.perf_counter()-start,python=platform.python_version(),
          platform=platform.platform(),seed=CONFIG['seed'],driver_sha256=hashlib.sha256(EXECUTED_SOURCE).hexdigest(),cpu_threads=1))
    dump(out/'manifest.json',dict(files_sha256={str(p.relative_to(out)).replace('\\','/'):sha(p)
          for p in sorted(out.rglob('*')) if p.is_file() and p.name!='manifest.json'}))
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);p.add_argument('--restart-from',type=Path)
    args=p.parse_args();result=run(args.out,args.restart_from);print(json.dumps(result,indent=2))
    if not result['atomistic']['status'].startswith('completed_'):
        raise SystemExit(2)
