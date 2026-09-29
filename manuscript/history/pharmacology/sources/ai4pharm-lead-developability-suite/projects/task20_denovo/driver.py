"""Task 20: reproducible fragment-graph NSGA-II with actual 3D pocket scoring.

100 candidates x 20 generations by default. Scores are hypotheses, not Kd.
Optional Vina evaluates a separate subset and an X77 redocking control.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,os,platform,subprocess,sys,time
from pathlib import Path
for key in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'):
    os.environ[key]='1'
import numpy as np
from scipy.spatial.transform import Rotation
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rdkit import Chem, RDConfig, rdBase
from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski, QED, Draw, rdMolDescriptors, ChemicalFeatures
try:
    from rdkit.Contrib.SA_Score import sascorer
except ImportError:
    sys.path.insert(0,str(Path(RDConfig.RDContribDir)/'SA_Score'))
    import sascorer

HERE=Path(__file__).resolve().parent
FEATURES=ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir)/'BaseFeatures.fdef'))
DEFAULT={'seed':20260921,'population':100,'generations':20,'poses_per_candidate':40,
         'crossover_probability':0.8,'mutation_probability':0.4,'vina_subset':12,
         'vina_exhaustiveness':4,'vina_timeout_s':180,'box_size_A':24.0,
         'score_status':'uncalibrated dimensionless geometry; no nanomolar conversion',
         'reference':'6W63 X77 redocking, 6YB7 unliganded active site as generation target'}
CORES=['[1*]c1ccc([2*])cc1','[1*]c1cc([2*])ccc1','[1*]c1ccnc([2*])c1',
       '[1*]c1ncc([2*])cn1','[1*]c1cc([2*])cs1','[1*]c1cc([2*])co1',
       '[1*]C1CCN([2*])CC1','[1*]N1CCN([2*])CC1']
LEFT=['[1*]C(=O)NC','[1*]C(=O)NCC','[1*]C(=O)Nc1ccccc1','[1*]C(=O)NC1CC1',
      '[1*]C(=O)NCCO','[1*]NC(=O)C','[1*]NC(=O)c1ccncc1','[1*]NC(=O)C1CC1',
      '[1*]CO','[1*]OC','[1*]OCC','[1*]N1CCOCC1','[1*]N1CCCCC1',
      '[1*]c1ccncc1','[1*]c1ccccc1','[1*]S(=O)(=O)NC','[1*]CCNC(=O)C',
      '[1*]C(=O)N1CCOCC1','[1*]C(=O)NC1CCCC1','[1*]NC(=O)CO']
RIGHT=['[2*]F','[2*]Cl','[2*]C','[2*]CC','[2*]C#N','[2*]OC','[2*]OCC',
       '[2*]CO','[2*]C(=O)N','[2*]C(=O)NC','[2*]NC(=O)C','[2*]C1CC1',
       '[2*]c1ccccc1','[2*]c1ccncc1','[2*]N1CCOCC1','[2*]C(F)(F)F']

def dump(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

def table(path,records):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(records[0]));w.writeheader();w.writerows(records)

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def join_fragment(a,b,label):
    combined=Chem.CombineMols(a,b)
    points=[x for x in combined.GetAtoms() if x.GetAtomicNum()==0 and x.GetIsotope()==label]
    if len(points)!=2 or any(x.GetDegree()!=1 for x in points):raise ValueError('Exactly two terminal anchors required')
    edit=Chem.RWMol(combined)
    edit.AddBond(points[0].GetNeighbors()[0].GetIdx(),points[1].GetNeighbors()[0].GetIdx(),Chem.BondType.SINGLE)
    for idx in sorted([p.GetIdx() for p in points],reverse=True):edit.RemoveAtom(idx)
    result=edit.GetMol();Chem.SanitizeMol(result);return result

def assemble(genes):
    m=join_fragment(Chem.MolFromSmiles(CORES[genes[0]]),Chem.MolFromSmiles(LEFT[genes[1]]),1)
    return join_fragment(m,Chem.MolFromSmiles(RIGHT[genes[2]]),2)

def desirability_low(x,a,b):return float(np.clip((b-x)/(b-a),0,1))

def descriptors(mol):
    mw=Descriptors.MolWt(mol);lp=Crippen.MolLogP(mol);tpsa=rdMolDescriptors.CalcTPSA(mol);hbd=Lipinski.NumHDonors(mol)
    # Explicitly a structural heuristic; no experimental pKa calibration.
    base=mol.HasSubstructMatch(Chem.MolFromSmarts('[N;X3;!$(N-C=O);!$(N-S(=O)=O)]'))
    aromatic=mol.HasSubstructMatch(Chem.MolFromSmarts('[nH0;+0]'))
    pka=9.0 if base else 5.0 if aromatic else 0.0
    ld=lp-math.log10(1+10**(pka-7.4))
    window=min(1-desirability_low(tpsa,20,40),desirability_low(tpsa,90,120))
    mpo=sum([desirability_low(lp,3,5),desirability_low(ld,2,4),desirability_low(mw,360,500),window,desirability_low(hbd,.5,3.5),desirability_low(pka,8,10)])
    return dict(mw_g_mol=mw,clogp=lp,tpsa_A2=tpsa,hbd=hbd,pka_assumed=pka,logd74_proxy=ld,
                cns_mpo_proxy=mpo,qed=float(QED.qed(mol)),sa_score=float(sascorer.calculateScore(mol)),
                rotatable_bonds=Lipinski.NumRotatableBonds(mol),fsp3=rdMolDescriptors.CalcFractionCSP3(mol))

def pdb_atoms(path):
    records=[];seen=set()
    for line in path.read_text().splitlines():
        if not line.startswith('ATOM  ') or line[21]!='A' or line[16] not in (' ','A'):continue
        element=line[76:78].strip()
        key=(int(line[22:26]),line[12:16].strip())
        if element=='H' or key in seen:continue
        seen.add(key);records.append((key,line[17:20],element,np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]),line))
    if len(records)<1000:raise ValueError('Protein structure not loaded completely')
    return records

def pocket(path):
    atoms=pdb_atoms(path);index={r[0]:r[3] for r in atoms}
    landmarks=[index[(143,'N')],index[(163,'NE2')],index[(166,'O')],index[(145,'SG')],index[(41,'NE2')]]
    center=np.mean(landmarks,axis=0)
    nearby=[r for r in atoms if np.linalg.norm(r[3]-center)<13]
    radii={'C':1.70,'N':1.55,'O':1.52,'S':1.80}
    return dict(atoms=np.array([r[3] for r in nearby]),radii=np.array([radii.get(r[2],1.7) for r in nearby]),
                center=center,acceptor_targets=np.array(landmarks[:2]),donor_targets=np.array([landmarks[2]]),
                protein_atoms=len(atoms),pocket_atoms=len(nearby))

def pose_score(coords,mol,p):
    xyz=coords if coords.ndim==3 else coords[None,:,:]
    dist=np.linalg.norm(xyz[:,:,None,:]-p['atoms'][None,None,:,:],axis=-1)
    radii=np.array([Chem.GetPeriodicTable().GetRvdw(a.GetAtomicNum()) for a in mol.GetAtoms()])
    gap=dist-radii[None,:,None]-p['radii'][None,None,:]
    clash=(np.maximum(-gap-.30,0)**2).sum(axis=(1,2))/mol.GetNumAtoms()
    contact=np.exp(-((gap-.6)/1.25)**2).sum(axis=(1,2))/mol.GetNumAtoms()
    features=FEATURES.GetFeaturesForMol(mol)
    accept=sorted({i for f in features if f.GetFamily()=='Acceptor' for i in f.GetAtomIds()})
    donor=sorted({i for f in features if f.GetFamily()=='Donor' for i in f.GetAtomIds()})
    feature=np.zeros(len(xyz))
    for ids,targets in ((accept,p['acceptor_targets']),(donor,p['donor_targets'])):
        if ids:
            d=np.linalg.norm(xyz[:,ids,None,:]-targets[None,None,:,:],axis=-1)
            feature+=np.exp(-((d-2.9)/.65)**2).max(axis=1).sum(axis=1)
    # All dimensionless. No conversion to kcal/mol, Ki, Kd or clinical effect.
    return contact+2*feature-12*clash,clash,feature

def evaluate(genes,p,c,cache):
    mol=assemble(genes);smiles=Chem.MolToSmiles(mol,isomericSmiles=True)
    if smiles in cache:return cache[smiles]
    des=descriptors(mol)
    if not 150<=des['mw_g_mol']<=550:raise ValueError('Outside prespecified molecular weight domain')
    seed=int(hashlib.sha256((smiles+str(c['seed'])).encode()).hexdigest()[:8],16)%2147483647
    explicit=Chem.AddHs(mol);prm=AllChem.ETKDGv3();prm.randomSeed=seed;prm.numThreads=1
    if AllChem.EmbedMolecule(explicit,prm)!=0:raise ValueError('Embedding failed')
    convergence=AllChem.MMFFOptimizeMolecule(explicit,mmffVariant='MMFF94s',maxIters=250)
    heavy=Chem.RemoveHs(explicit);xyz=np.array(heavy.GetConformer().GetPositions());xyz-=xyz.mean(axis=0)
    rng=np.random.default_rng(seed)
    rotations=Rotation.random(c['poses_per_candidate'],random_state=rng).as_matrix()
    trials=np.einsum('ai,nji->naj',xyz,rotations)+p['center'][None,None,:]+rng.normal(0,1.7,(len(rotations),1,3))
    scores,clashes,features=pose_score(trials,heavy,p);best=int(np.argmax(scores));chosen=trials[best]
    # Local translation refinement after the random orientation search.
    for scale in (1.0,.5,.25):
        for refinement in range(3):
            shifts=np.zeros((7,3));shifts[1:4]=np.eye(3)*scale;shifts[4:]=-np.eye(3)*scale
            proposal=chosen[None,:,:]+shifts[:,None,:]
            s,cl,ft=pose_score(proposal,heavy,p);k=int(np.argmax(s));chosen=proposal[k]
    score,clash,feature=pose_score(chosen,heavy,p)
    for i,coord in enumerate(chosen):heavy.GetConformer().SetAtomPosition(i,coord)
    record={'id':'G'+hashlib.sha256(smiles.encode()).hexdigest()[:12],'smiles':smiles,'genes':list(genes),**des,
            'pocket_geometry_score':float(score[0]),'steric_overlap_per_atom_A2':float(clash[0]),
            'pharmacophore_contact_proxy':float(feature[0]),'mmff_status':int(convergence),
            'predicted_kd_nM':None,'confirmed_novel':None,'synthetic_route_validated':False,
            '_mol':heavy}
    cache[smiles]=record;return record

def objective(r):return np.array([r['pocket_geometry_score'],r['cns_mpo_proxy']/6,r['qed'],-r['sa_score']/10])

def pareto_fronts(records):
    values=np.array([objective(r) for r in records]);n=len(records)
    dominates=np.all(values[:,None,:]>=values[None,:,:],axis=2)&np.any(values[:,None,:]>values[None,:,:],axis=2)
    counts=dominates.sum(axis=0);remaining=set(range(n));fronts=[]
    while remaining:
        front=sorted(i for i in remaining if counts[i]==0)
        if not front:raise RuntimeError('Dominance cycle')
        fronts.append(front);remaining.difference_update(front)
        counts-=dominates[front].sum(axis=0)
    return fronts

def select_population(records,n):
    unique={r['smiles']:r for r in records};records=list(unique.values());selected=[]
    for front in pareto_fronts(records):
        if len(selected)+len(front)<=n:selected.extend(records[i] for i in front);continue
        values=np.array([objective(records[i]) for i in front]);crowd=np.zeros(len(front))
        for column in range(values.shape[1]):
            order=np.argsort(values[:,column],kind='stable');span=np.ptp(values[:,column])
            if span==0:continue
            crowd[order[[0,-1]]]=np.inf
            crowd[order[1:-1]]+=(values[order[2:],column]-values[order[:-2],column])/span
        ranked=np.argsort(-crowd,kind='stable')[:n-len(selected)]
        selected.extend(records[front[k]] for k in ranked);break
    return selected

def evolve(p,c):
    rng=np.random.default_rng(c['seed']);cache={};rejected=[];limits=[len(CORES),len(LEFT),len(RIGHT)]
    def random_genes():return [int(rng.integers(n)) for n in limits]
    population=[];attempts=0
    while len({r['smiles'] for r in population})<c['population']:
        genes=random_genes();attempts+=1
        try:population.append(evaluate(genes,p,c,cache))
        except (ValueError,RuntimeError) as exc:rejected.append({'genes':str(genes),'reason':str(exc)})
        if attempts>10000:raise RuntimeError('Insufficient valid initial molecules')
    population=list({r['smiles']:r for r in population}.values())
    history=[];summaries=[];initial=[dict(r) for r in population]
    for generation in range(c['generations']):
        if generation:
            offspring=[];tries=0
            while len({r['smiles'] for r in offspring})<c['population']:
                a,b=[population[int(rng.integers(len(population)))] for _ in range(2)]
                genes=list(a['genes'])
                if rng.random()<c['crossover_probability']:
                    cut=int(rng.integers(1,3));genes[cut:]=b['genes'][cut:]
                if rng.random()<c['mutation_probability']:
                    k=int(rng.integers(3));genes[k]=int(rng.integers(limits[k]))
                if tries%5==0:genes=random_genes()  # explicit random immigrants
                try:offspring.append(evaluate(genes,p,c,cache))
                except (ValueError,RuntimeError) as exc:rejected.append({'genes':str(genes),'reason':str(exc)})
                tries+=1
                if tries>10000:raise RuntimeError('Genetic proposal budget exhausted')
            population=select_population(population+offspring,c['population'])
        front=set(pareto_fronts(population)[0])
        for i,r in enumerate(population):history.append({'generation':generation,'pareto':i in front,**{k:v for k,v in r.items() if not k.startswith('_')}})
        summaries.append(dict(generation=generation,population=len(population),front_size=len(front),
                              evaluated_unique=len(cache),best_pocket_score=max(r['pocket_geometry_score'] for r in population),
                              median_pocket_score=float(np.median([r['pocket_geometry_score'] for r in population])),
                              mean_mpo=float(np.mean([r['cns_mpo_proxy'] for r in population])),
                              sa_below_3p5=sum(r['sa_score']<3.5 for r in population)))
        print(f"Task20 generation {generation+1}/{c['generations']}: {len(cache)} unique molecules evaluated",flush=True)
    return population,initial,history,summaries,cache,rejected

def random_control(p,c,cache,budget):
    """Equal unique-graph budget, independent random order, identical pose evaluator.

    Cache reuse avoids repeating identical physics evaluations, not information
    sharing in the random selection rule. This is one seed, not a superiority test.
    """
    import itertools
    rng=np.random.default_rng(c['seed']+4001)
    genes=list(itertools.product(range(len(CORES)),range(len(LEFT)),range(len(RIGHT))))
    order=rng.permutation(len(genes));picked={};rejects=0
    for idx in order:
        try:
            r=evaluate(genes[int(idx)],p,c,cache);picked[r['smiles']]=r
        except (ValueError,RuntimeError):rejects+=1
        if len(picked)==budget:break
    if len(picked)!=budget:raise RuntimeError('Equal-budget random control space exhausted')
    records=list(picked.values());selected=select_population(records,c['population'])
    return records,{'unique_graph_budget':budget,'seed':c['seed']+4001,'rejected_proposals':rejects,
        'best_geometry_score':max(r['pocket_geometry_score'] for r in selected),
        'selected_median_geometry_score':float(np.median([r['pocket_geometry_score'] for r in selected])),
        'selected_sa_below_3p5':sum(r['sa_score']<3.5 for r in selected),
        'interpretation':'one equal-budget random-order control; not replicated optimization superiority'}

def clean_protein(source,destination):
    # pdb_atoms has selected blank/A conformers; normalize that explicit choice.
    lines=[r[-1][:16]+' '+r[-1][17:] for r in pdb_atoms(source)]
    destination.write_text('\n'.join(lines)+'\nTER\nEND\n',encoding='ascii',newline='\n')

def run_command(cmd,log,timeout=180):
    with log.open('w',encoding='utf-8') as handle:
        result=subprocess.run([str(x) for x in cmd],stdout=handle,stderr=subprocess.STDOUT,timeout=timeout,check=False)
    if result.returncode:raise RuntimeError(f'External calculation failed; see {log}')

def docking(out,population,p,c,vina,prep_python):
    if not vina:return {'status':'not_executed; provide --vina and --preparation-python'},[]
    folder=out/'docking';folder.mkdir();receptors={}
    for pid in ('6YB7','6W63'):
        clean=folder/(pid+'_protein.pdb');clean_protein(HERE/'inputs'/f'{pid}.pdb',clean)
        base=folder/(pid+'_receptor')
        run_command([prep_python,'-m','meeko.cli.mk_prepare_receptor','--read_pdb',clean,'-o',base,'-p'],folder/(pid+'_preparation.txt'))
        receptors[pid]=base.with_suffix('.pdbqt')
        if not receptors[pid].exists():
            alternate=folder/(pid+'_receptor_rigid.pdbqt')
            if not alternate.exists():raise RuntimeError('Prepared receptor absent')
            receptors[pid]=alternate
    ranked=sorted(population,key=lambda r:r['pocket_geometry_score'],reverse=True)[:c['vina_subset']]
    records=[]
    def dock_one(label,mol,pid,center,seed,exhaustiveness):
        input_sdf=folder/(label+'.sdf');writer=Chem.SDWriter(str(input_sdf));writer.write(Chem.AddHs(mol,addCoords=True));writer.close()
        ligand=folder/(label+'.pdbqt')
        run_command([prep_python,'-m','meeko.cli.mk_prepare_ligand','-i',input_sdf,'-o',ligand],folder/(label+'_prepare.txt'))
        output=folder/(label+'_pose.pdbqt')
        cmd=[vina,'--receptor',receptors[pid],'--ligand',ligand,'--cpu','1','--seed',str(seed),'--exhaustiveness',str(exhaustiveness),'--num_modes','3','--out',output]
        for axis,x in zip('xyz',center):cmd += ['--center_'+axis,str(float(x)),'--size_'+axis,str(c['box_size_A'])]
        run_command(cmd,folder/(label+'_vina.txt'),c['vina_timeout_s'])
        score=float(next(line.split()[3] for line in output.read_text().splitlines() if line.startswith('REMARK VINA RESULT:')))
        return score,output
    for i,r in enumerate(ranked):
        score,pose=dock_one(r['id'],r['_mol'],'6YB7',p['center'],c['seed']+i,c['vina_exhaustiveness'])
        records.append({'id':r['id'],'smiles':r['smiles'],'vina_score_kcal_mol':score,'geometry_score':r['pocket_geometry_score'],'cns_mpo_proxy':r['cns_mpo_proxy'],'qed':r['qed'],'sa_score':r['sa_score'],'predicted_kd_nM':None,'pose_file':pose.relative_to(out).as_posix()})
    template=Chem.SDMolSupplier(str(HERE/'inputs/X77_ideal.sdf'),removeHs=True)[0]
    lines=[line for line in (HERE/'inputs/6W63.pdb').read_text().splitlines() if line.startswith('HETATM') and line[17:20]=='X77' and line[76:78].strip()!='H']
    crystal=Chem.MolFromPDBBlock('\n'.join(lines)+'\nEND\n',removeHs=True,sanitize=False)
    crystal=AllChem.AssignBondOrdersFromTemplate(template,crystal)
    center=np.array(crystal.GetConformer().GetPositions()).mean(axis=0)
    control=[]
    # Independently embedded conformer: crystal coordinates do not seed the search.
    ligand=Chem.AddHs(Chem.RemoveHs(crystal));prm=AllChem.ETKDGv3();prm.randomSeed=c['seed'];prm.numThreads=1
    if AllChem.EmbedMolecule(ligand,prm):raise RuntimeError('X77 conformer failed')
    AllChem.MMFFOptimizeMolecule(ligand,maxIters=500);ligand=Chem.RemoveHs(ligand)
    for seed in (c['seed']+101,c['seed']+102):
        label='X77_redock_'+str(seed);score,pose=dock_one(label,ligand,'6W63',center,seed,8)
        # Meeko preserves original graph and atom mapping in the output PDBQT.
        exported=folder/(label+'_poses.sdf')
        run_command([prep_python,'-m','meeko.cli.mk_export',pose,'-s',exported],folder/(label+'_export.txt'))
        docked=Chem.SDMolSupplier(str(exported),removeHs=True)[0]
        matches=docked.GetSubstructMatches(crystal,uniquify=False,useChirality=False,maxMatches=1000)
        if not matches:raise RuntimeError('Cannot map X77 heavy atoms for redocking RMSD')
        xref=np.array(crystal.GetConformer().GetPositions());xpose=np.array(docked.GetConformer().GetPositions())
        rmsd=min(float(np.sqrt(np.mean(np.sum((xpose[list(match)]-xref)**2,axis=1)))) for match in matches)
        control.append({'seed':seed,'vina_score_kcal_mol':score,'symmetry_aware_heavy_atom_rmsd_A_without_superposition':rmsd,'within_2A':rmsd<=2})
    table(folder/'candidate_vina_scores.csv',records)
    return {'status':'executed','candidate_docks':len(records),'redocking_controls':control,'vina_executable_sha256':sha(Path(vina)),
            'vina_version':subprocess.check_output([str(vina),'--version'],text=True).strip(),
            'receptor':'rigid chain A, default Meeko residue protonation; crystal waters removed; no microstate ensemble',
            'nanomolar_affinity_validated':False},records

def figures(out,history,summaries,population,dockrows):
    folder=out/'figures';folder.mkdir()
    fig=plt.figure(figsize=(15,11));gs=fig.add_gridspec(3,3,height_ratios=[1,1,.85],hspace=.5,wspace=.34)
    ax=fig.add_subplot(gs[0,0]);gens=np.array([r['generation'] for r in summaries]);ax.plot(gens,[r['best_pocket_score'] for r in summaries],label='Best');ax.plot(gens,[r['median_pocket_score'] for r in summaries],label='Median');ax.set(xlabel='Generation (0–19)',ylabel='Geometry score (dimensionless)',title='A  Pocket-conditioned evolution');ax.legend()
    ax=fig.add_subplot(gs[0,1]);first=[r for r in history if r['generation']==0];last=[r for r in history if r['generation']==gens[-1]]
    ax.scatter([r['cns_mpo_proxy'] for r in first],[r['pocket_geometry_score'] for r in first],color='#b5c3cb',s=17,label='Generation 0')
    sc=ax.scatter([r['cns_mpo_proxy'] for r in last],[r['pocket_geometry_score'] for r in last],c=[r['sa_score'] for r in last],cmap='viridis_r',s=25);fig.colorbar(sc,ax=ax,label='SA score');ax.set(xlabel='CNS-MPO with assumed pKa',ylabel='Geometry score',title='B  Multicriterion candidate set');ax.legend(fontsize=8)
    ax=fig.add_subplot(gs[0,2]);ax.plot(gens,[r['front_size'] for r in summaries],label='Pareto front');ax.plot(gens,[r['sa_below_3p5'] for r in summaries],label='SA < 3.5');ax.set(xlabel='Generation',ylabel='Candidates / 100',title='C  Diversity and accessibility proxy');ax.legend(fontsize=8)
    ax=fig.add_subplot(gs[1,0]);
    if dockrows:
        ax.scatter([r['geometry_score'] for r in dockrows],[r['vina_score_kcal_mol'] for r in dockrows],c='#18857e');ax.set(xlabel='Geometry score',ylabel='Vina score (kcal/mol)',title='D  Independent subset docking')
    else:ax.text(.1,.5,'Vina not executed',transform=ax.transAxes);ax.set_title('D  Docking status')
    top=sorted(population,key=lambda r:r['pocket_geometry_score'],reverse=True)[:5]
    for n,r in enumerate(top):
        position=gs[1,1+n] if n<2 else gs[2,n-2]
        ax=fig.add_subplot(position);mol=Chem.MolFromSmiles(r['smiles']);AllChem.Compute2DCoords(mol)
        ax.imshow(Draw.MolToImage(mol,size=(600,280)));ax.axis('off');ax.set_title(f"{r['id']}\nMPO {r['cns_mpo_proxy']:.2f} · SA {r['sa_score']:.2f}",fontsize=10)
    fig.suptitle('20  |  Fragment-graph generation in a real 3D pocket',x=.055,ha='left',fontsize=19,fontweight='bold',color='#18364d')
    fig.text(.055,.935,'6YB7 Mpro unliganded active site • NSGA-II, 100 candidates × 20 generations • held-out X77 redocking control',fontsize=10,color='#526778')
    fig.text(.055,.022,'Geometry and docking are uncalibrated rankings. No Kd conversion, novelty claim, synthetic-route proof or clinical lead validation.\nCNS-MPO uses assumed pKa; Mpro inhibition does not itself require CNS exposure. All displayed structures are algorithm outputs.',fontsize=9,color='#526778')
    fig.subplots_adjust(top=.88,bottom=.085,left=.07,right=.965)
    fig.savefig(folder/'fig20_denovo_pareto_lead_optimization.png',dpi=300);plt.close(fig)

def run(out,config=None,vina=None,preparation_python=None):
    out=Path(out);c=dict(DEFAULT,**(config or {}))
    if out.exists():raise FileExistsError('Use a new output directory')
    if c['population']<4 or c['generations']<1 or c['poses_per_candidate']<4:raise ValueError('Insufficient configured run dimensions')
    out.mkdir(parents=True);start=time.perf_counter();p=pocket(HERE/'inputs/6YB7.pdb')
    dump(out/'config.json',c)
    population,initial,history,summaries,cache,rejected=evolve(p,c)
    evolution_evaluated=len(cache)
    baseline,baseline_summary=random_control(p,c,cache,evolution_evaluated)
    table(out/'random_search_control.csv',[{k:v for k,v in r.items() if not k.startswith('_')} for r in baseline])
    dump(out/'random_search_summary.json',baseline_summary)
    table(out/'population_history.csv',history);table(out/'generation_summary.csv',summaries)
    table(out/'evaluated_candidates.csv',[{k:v for k,v in r.items() if not k.startswith('_')} for r in cache.values()])
    dump(out/'rejected_proposals.json',rejected)
    writer=Chem.SDWriter(str(out/'final_candidates.sdf'))
    for r in population:
        m=r['_mol'];m.SetProp('_Name',r['id']);m.SetProp('canonical_smiles',r['smiles']);m.SetProp('geometry_score',str(r['pocket_geometry_score']));m.SetProp('evidence','generated graph and model pose, not experimentally validated');writer.write(m)
    writer.close()
    dockinfo,dockrows=docking(out,population,p,c,vina,preparation_python or sys.executable)
    figures(out,history,summaries,population,dockrows)
    check={'all_populations_exact':all(r['population']==c['population'] for r in summaries),
           'all_final_graphs_valid':all(Chem.MolFromSmiles(r['smiles']) is not None for r in population),
           'all_no_dummy_atoms':all(not any(a.GetAtomicNum()==0 for a in r['_mol'].GetAtoms()) for r in population),
           'finite_3d_scores':all(math.isfinite(r['pocket_geometry_score']) for r in cache.values()),
           'mpo_bounds':all(0<=r['cns_mpo_proxy']<=6 for r in cache.values()),
           'no_duplicate_final_smiles':len({r['smiles'] for r in population})==c['population'],
           'best_score_not_decreased':summaries[-1]['best_pocket_score']>=summaries[0]['best_pocket_score']}
    summary={'task':20,'population':c['population'],'generations':c['generations'],'history_rows':len(history),
             'unique_graphs_evaluated':evolution_evaluated,'unique_graphs_including_random_control':len(cache),
             'random_search_control':baseline_summary,'rejected_proposals':len(rejected),'final_pareto_size':len(pareto_fronts(population)[0]),
             'final_sa_below_3p5':sum(r['sa_score']<3.5 for r in population),'best_geometry_score_initial':summaries[0]['best_pocket_score'],
             'best_geometry_score_final':summaries[-1]['best_pocket_score'],'docking':dockinfo,
             'affinity_nM':None,'novelty_verified':False,'synthesis_proven':False,'clinical_recommendation':None,
             'scientific_acceptance':'nanomolar binding and novel synthesizable lead claims not established',
             'protein_atoms':p['protein_atoms'],'pocket_atoms':p['pocket_atoms'],'seconds':time.perf_counter()-start,
             'runtime':{'python':platform.python_version(),'rdkit':rdBase.rdkitVersion,'numpy':np.__version__},
             'source_sha256':sha(Path(__file__)),'inputs_sha256':{x.name:sha(x) for x in sorted((HERE/'inputs').glob('*')) if x.is_file()}}
    dump(out/'summary.json',summary);dump(out/'verification.json',{'passed':all(check.values()),'checks':check})
    dump(out/'manifest.json',{'files_sha256':{x.relative_to(out).as_posix():sha(x) for x in sorted(out.rglob('*')) if x.is_file()}})
    if not all(check.values()):raise RuntimeError('Numerical or graph verification failed')
    return summary

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);parser.add_argument('--config',type=Path);parser.add_argument('--vina',type=Path);parser.add_argument('--preparation-python',type=Path);a=parser.parse_args()
    config=json.loads(a.config.read_text()) if a.config else None
    print(json.dumps(run(a.out,config,a.vina,a.preparation_python),indent=2))
