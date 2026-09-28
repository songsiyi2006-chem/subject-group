"""Replay the supplied NEB with captured random weights; audit its saved claims.

All energy/force units below are nominal source conventions, not calibration.
The replay records algorithm state without changing the supplied NEB method.
"""
import ast
from collections import Counter
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time
import numpy as np
import torch

BASE = Path(__file__).resolve().parents[1]
OUT = BASE/'results/source_audit'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path,data):
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')


def write_csv(path,rows):
    with path.open('w',encoding='utf-8',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def minimum_pair(coords):
    d = torch.cdist(coords.detach(),coords.detach())
    d.fill_diagonal_(float('inf'))
    index = int(d.argmin())
    return float(d.min()), index//len(coords), index%len(coords)


def main():
    start=time.perf_counter(); OUT.mkdir(parents=True,exist_ok=True)
    torch.set_num_threads(1)
    source=BASE/'source/synthapore_engine_compat.py'
    spec=importlib.util.spec_from_file_location('synthapore_supplied_audit',source)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    capture_path=BASE/'results/compatibility/captured_results.json'
    captured=json.loads(capture_path.read_text())
    weights=BASE/'results/compatibility/untrained_source_weights.pt'
    model=module.NeuralInteratomicPotential(num_species=30,hidden_dim=48,num_layers=3)
    model.load_state_dict(torch.load(weights,map_location='cpu',weights_only=True)); model.eval()
    z=torch.tensor(captured['atomic_numbers'],dtype=torch.long)
    r0=torch.tensor(captured['r_init'],dtype=torch.float32)
    r1=torch.tensor(captured['r_final'],dtype=torch.float32)
    edges=torch.tensor(captured['edge_index'],dtype=torch.long)
    evaluations=[]
    def traced_potential(numbers,coords,index):
        energy,force=model(numbers,coords,index)
        close,i,j=minimum_pair(coords)
        evaluations.append(dict(iteration=len(evaluations)//7,image=len(evaluations)%7,
             nominal_energy_eV=float(energy.detach()),max_atomic_force_nominal_eV_A=float(force.detach().norm(dim=1).max()),
             min_pair_distance_A=close,min_pair_i=i,min_pair_j=j))
        return energy,force
    frames={}
    def profiler(frame,event,arg):
        if event=='return' and frame.f_code is module.ClimbingImageNEB.optimize_pathway.__code__:
            for name in ['band_coords','current_energies','climbing_idx','energy_profile_history']:
                frames[name]=frame.f_locals[name]
    prior=sys.getprofile()
    try:
        sys.setprofile(profiler)
        result=module.ClimbingImageNEB(traced_potential,num_images=7,spring_k=4.5).optimize_pathway(z,r0,r1,edges,max_iterations=45)
    finally:
        sys.setprofile(prior)
    assert len(evaluations)==315 and result==captured['neb_results']
    bands=frames['band_coords']; energy_final=[]; forces=[]
    for coords in bands:
        e,f=model(z,coords.clone().detach().requires_grad_(True),edges)
        energy_final.append(float(e.detach())); forces.append(f.detach())
    climbing=1+int(np.argmax(energy_final[1:-1]))
    residuals=[]; final_rows=[]
    for i in range(7):
        residual=None
        if 0<i<6:
            tangent=bands[i+1]-bands[i-1]; tangent=tangent/(torch.norm(tangent)+1e-8)
            parallel=torch.sum(forces[i]*tangent)*tangent
            if i==climbing:
                neb_force=forces[i]-2*parallel
            else:
                neb_force=forces[i]-parallel+4.5*(torch.norm(bands[i+1]-bands[i])-torch.norm(bands[i]-bands[i-1]))*tangent
            residual=float(neb_force.norm(dim=1).max()); residuals.append(residual)
        distance,j,k=minimum_pair(bands[i])
        final_rows.append(dict(image=i,reported_preupdate_energy_nominal_eV=frames['current_energies'][i],
                  reevaluated_postupdate_energy_nominal_eV=energy_final[i],
                  energy_change_nominal_eV=energy_final[i]-frames['current_energies'][i],
                  true_max_atomic_force_nominal_eV_A=float(forces[i].norm(dim=1).max()),
                  projected_max_atomic_force_nominal_eV_A=residual,min_pair_distance_A=distance,min_pair_i=j,min_pair_j=k))
    write_csv(OUT/'NEB_replay_evaluations.csv',evaluations)
    write_csv(OUT/'NEB_final_images.csv',final_rows)
    dump(OUT/'NEB_final_coordinates.json',{'atomic_numbers':z.tolist(),'positions_A':[r.tolist() for r in bands],
          'evidence':'Source random-weight band after final update; not a validated transition state'})
    md=captured['md_results']; md_rows=[]
    for t,ep,ek,temp in zip(md['time_series_fs'],md['potential_energy_eV'],md['kinetic_energy_eV'],md['temperature_profile_K']):
        factor=max(1+.5/50*(298.15/max(temp,1e-5)-1),.1)
        md_rows.append(dict(source_time_label_fs=t,actual_integration_time_fs=t+.5,
                 nominal_potential_energy_eV=ep,nominal_kinetic_energy_before_scaling_eV=ek,
                 source_temperature_before_scaling_K=temp,inferred_temperature_after_scaling_K=temp*factor))
    write_csv(OUT/'MD_telemetry_audit.csv',md_rows)
    tree=ast.parse(source.read_text(encoding='utf-8'))
    classes=[n.name for n in tree.body if isinstance(n,ast.ClassDef)]
    text=source.read_text(encoding='utf-8')
    report={
       'scope':'Exact source NEB replay and saved MD/source-code accounting; no new physical evidence',
       'source_sha256':sha(source),'captured_sha256':sha(capture_path),'weights_sha256':sha(weights),
       'script_sha256':sha(Path(__file__)),
       'counts':{'source_NEB_replays':1,'NEB_iterations':45,'NEB_images':7,'NEB_replay_energy_force_calls':315,
                 'final_band_energy_force_calls':7,'saved_MD_frames_audited':30,'new_MD_trajectories':0},
       'source_claims':{'classes':classes,'trained_parameters':False,'diffusion_noise_schedule':False,
                        'reverse_diffusion_sampler':False,'atomic_framework_coordinates':False,'DFT_or_experiments':0},
       'input':{'atomic_number_counts':dict(Counter(captured['atomic_numbers'])),'formula_inventory':'C4H2CuN',
                'directed_edges':len(captured['edge_index'][0]),'charge_spin_specified':False,
                'note':'Eight-node manually assigned fragment; not a full 2-phenylindole molecule or Cu-N4 material'},
       'NEB':{'exact_source_return_matches_capture':result==captured['neb_results'],
              'source_return':result,'reported_energy_max_image':int(np.argmax(frames['current_energies'])),
              'postupdate_energy_max_image':int(np.argmax(energy_final)),
              'postupdate_climbing_image':climbing,
              'max_stale_energy_difference_nominal_eV':float(np.max(np.abs(np.array(energy_final)-frames['current_energies']))),
              'final_max_projected_atomic_force_nominal_eV_A':max(residuals),
              'endpoint_max_atomic_forces_nominal_eV_A':[final_rows[i]['true_max_atomic_force_nominal_eV_A'] for i in [0,6]],
              'minimum_final_band_pair_distance_A':min(r['min_pair_distance_A'] for r in final_rows),
              'source_has_force_tolerance':False,'source_has_hessian_check':False,
              'barrier_is_free_energy':False,'source_NEB_uses_MD_external_field':False},
       'MD':{'reported_mean_K':md['mean_simulated_temp_K'],'recomputed_mean_K':float(np.mean(md['temperature_profile_K'])),
             'target_K':md['target_temperature_K'],'first_actual_frame_fs':.5,'last_actual_frame_fs':145.5,
             'reported_total_fs':150.,'source_temperature_denominator_dof':24,
             'internal_dof_if_COM_removed_consistently':21,
             'temperature_caveat':'Initial COM is removed once, but field can recreate drift; internal temperature cannot be recovered without velocities',
             'stored_trajectory_and_velocities':False,'equilibrium_demonstrated':False,
             'berendsen_canonical_sampling':False},
       'limitations':['Random neural energy/force scales are not calibrated eV or eV/Angstrom',
                      'Energy profiles are potential energy, not Gibbs free energy',
                      'NEB graph is fixed and endpoint charge/spin and stationary minima are unspecified',
                      'All field/model and geometry observations are implementation diagnostics'],
       'elapsed_seconds':time.perf_counter()-start,
       'output_sha256':{p.relative_to(BASE).as_posix():sha(p) for p in sorted(OUT.iterdir()) if p.is_file()},
    }
    dump(BASE/'results/source_audit.json',report)
    print(json.dumps({k:report[k] for k in ['counts','NEB','MD']},indent=2))


if __name__=='__main__':
    main()
