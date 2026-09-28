"""Separate interpolation error from the imperfect RHF reference at matched geometries.

This is arithmetic on frozen predictions and quantum records, not new training
or quantum chemistry. No FCI interpolation or label substitution is performed.
"""
from pathlib import Path
import csv
import hashlib
import json
import math

BASE=Path(__file__).resolve().parents[1]
OUT=BASE/'results/extensions/error_budget'


def energy_components(predicted,training_reference,fresh_reference,correlated_reference):
    """Keep a small reproducibility drift term instead of hiding it in learning error."""
    values=[predicted,training_reference,fresh_reference,correlated_reference]
    if not all(math.isfinite(float(v)) for v in values):
        raise ValueError('All four energies must be finite and use one unit/reference.')
    learning=predicted-training_reference
    replay=training_reference-fresh_reference
    approximation=fresh_reference-correlated_reference
    total=predicted-correlated_reference
    return dict(learning_error_Hartree=learning,reference_replay_drift_Hartree=replay,
                RHF_minus_FCI_Hartree=approximation,total_error_vs_FCI_Hartree=total,
                reconstruction_residual_Hartree=total-learning-replay-approximation)


def rows(path):
    with path.open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(prediction_path=None):
    # The final table adapter is explicit, because inference record names must
    # not silently substitute physical reference methods.
    prediction_path=prediction_path or BASE/'results/extensions/learning/predictions.csv'
    quantum_path=BASE/'results/extensions/correlation/curve.csv'
    reference_path=BASE/'results/electronic/h2_reference.csv'
    quantum=rows(quantum_path)
    reference={r['point_id']:r for r in rows(reference_path)}
    keyed={}
    for row in quantum:
        if row['status'] not in ('success','converged','ok','passed'):
            continue
        key=(row['basis'].lower(),row['method'].upper(),round(float(row['R_A']),8))
        if key in keyed:
            raise ValueError('Duplicate quantum method/geometry in main curve.')
        keyed[key]=float(row['energy_Hartree'])
    output=[]
    for row in rows(prediction_path):
        old=reference[row['point_id']]
        radius=round(float(old['R_A']),8)
        if abs(float(row['R_A'])-float(old['R_A']))>1e-10:
            raise ValueError('Point identity and radius disagree.')
        keys=[('sto-3g',m,radius) for m in ('RHF','FCI')]
        if not all(k in keyed for k in keys):
            continue
        predicted=float(row['predicted_energy_Hartree'])
        old_energy=float(old['energy_total_Hartree'])
        components=energy_components(predicted,old_energy,keyed[keys[0]],keyed[keys[1]])
        output.append(dict(objective=row['objective'],seed=row['seed'],point_id=row['point_id'],
                           split=old['split'],R_A=float(old['R_A']),**components))
    if not output:
        raise ValueError('No exact shared geometries; do not interpolate FCI labels.')
    OUT.mkdir(parents=True,exist_ok=True)
    with (OUT/'matched_error_components.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(output[0]));writer.writeheader();writer.writerows(output)
    summary={'scope':'Signed energy-error decomposition at identical saved geometries, STO-3G only.',
             'equation':'E_NN - E_FCI = (E_NN - E_old_RHF) + (E_old_RHF - E_fresh_RHF) + (E_fresh_RHF - E_FCI)',
             'rows':len(output),'unique_geometries':len({r['R_A'] for r in output}),
             'counts':{'new_quantum_jobs':0,'new_training_runs':0},
             'max_identity_residual_Hartree':max(abs(r['reconstruction_residual_Hartree']) for r in output),
             'max_reference_replay_drift_Hartree':max(abs(r['reference_replay_drift_Hartree']) for r in output),
             'limits':['FCI is exact only within the chosen finite one-electron basis and electronic Hamiltonian.',
                       'Small learning error against RHF does not remove its electronic approximation error.',
                       'Signed cancellation can reduce total error accidentally; no FCI labels entered training.'],
             'inputs_sha256':{p.relative_to(BASE).as_posix():sha(p) for p in [Path(__file__),prediction_path,quantum_path,reference_path]},
             'outputs_sha256':{'matched_error_components.csv':sha(OUT/'matched_error_components.csv')}}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:summary[k] for k in ('rows','unique_geometries','max_identity_residual_Hartree','max_reference_replay_drift_Hartree')}))


if __name__=='__main__':
    run()
