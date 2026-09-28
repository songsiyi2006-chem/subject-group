"""Post hoc error attribution on frozen calculations; no new quantum jobs or training."""
from pathlib import Path
import csv
import hashlib
import json
import math
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'manuscript/results'
SOURCES = [
    'quantumequi/results/extensions/error_budget/matched_error_components.csv',
    'quantumequi/results/extensions/vibration/tail_followup/diagnostics.csv',
]

def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def decompose_mse(learning, bias):
    learning, bias = np.asarray(learning, dtype=float), np.asarray(bias, dtype=float)
    if learning.shape != bias.shape or learning.ndim != 1 or not learning.size:
        raise ValueError('Equal nonempty one-dimensional matched arrays are required.')
    if not np.isfinite(learning).all() or not np.isfinite(bias).all():
        raise ValueError('Energies must be finite.')
    total = learning + bias
    mse_l, mse_b = np.mean(learning**2), np.mean(bias**2)
    cross = 2*np.mean(learning*bias)
    mse_t = np.mean(total**2)
    rms_l, rms_b, rms_t = map(math.sqrt, (mse_l, mse_b, mse_t))
    return {
        'n_geometries':len(learning), 'learning_RMSE_Hartree':rms_l,
        'reference_bias_RMSE_Hartree':rms_b, 'total_FCI_RMSE_Hartree':rms_t,
        'learning_MSE_Hartree2':float(mse_l), 'reference_bias_MSE_Hartree2':float(mse_b),
        'cross_term_Hartree2':float(cross), 'total_MSE_Hartree2':float(mse_t),
        'MSE_identity_residual_Hartree2':float(mse_t-mse_l-mse_b-cross),
        'triangle_lower_bound_Hartree':abs(rms_b-rms_l),
        'triangle_upper_bound_Hartree':rms_b+rms_l,
        'signed_cancellation_points':int(np.sum(learning*bias<0)),
    }

def richardson_binding(coarse, fine):
    """Second-order extrapolation at one fixed boundary; not a new eigensolve."""
    hc, hf = float(coarse['grid_spacing_A']), float(fine['grid_spacing_A'])
    if float(coarse['right_A']) != float(fine['right_A']) or not np.isclose(hc/hf, 2):
        raise ValueError('Require the same domain and a factor-two grid refinement.')
    if int(coarse['bound_states_returned']) != 17 or int(fine['bound_states_returned']) != 17:
        raise ValueError('Both grids must contain the same targeted bound state.')
    bc, bf = float(coarse['last_state_binding_Hartree']), float(fine['last_state_binding_Hartree'])
    exact = float(fine['last_full_line_exact_binding_Hartree'])
    extrapolated = (4*bf-bc)/3
    return {'right_A':float(fine['right_A']), 'coarse_intervals':int(coarse['intervals']),
            'fine_intervals':int(fine['intervals']), 'exact_binding_Hartree':exact,
            'fine_binding_Hartree':bf, 'extrapolated_binding_Hartree':extrapolated,
            'fine_signed_relative_error_percent':100*(bf/exact-1),
            'extrapolated_signed_relative_error_percent':100*(extrapolated/exact-1)}

def write_csv(name, rows):
    with (OUT/name).open('w', encoding='utf-8', newline='') as stream:
        writer=csv.DictWriter(stream, fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    source = read_csv(ROOT/SOURCES[0])
    if max(abs(float(r['reference_replay_drift_Hartree'])) for r in source) != 0:
        raise ValueError('This two-term analysis requires the separately checked zero replay drift.')
    attribution=[]
    for split in ('test','ood'):
        # Stored OOD spelling is preserved below; never combine train/validation with held-out metrics.
        eligible=[r for r in source if r['split'].lower()==split or (split=='ood' and 'ood' in r['split'].lower())]
        for seed in ('7301','7302','7303'):
            for objective in ('energy_only','energy_gradient'):
                rows=sorted((r for r in eligible if r['seed']==seed and r['objective']==objective), key=lambda r:float(r['R_A']))
                if not rows:raise ValueError(f'Missing matched records: {split}/{seed}/{objective}')
                components=decompose_mse([float(r['learning_error_Hartree']) for r in rows], [float(r['RHF_minus_FCI_Hartree']) for r in rows])
                total=np.array([float(r['total_error_vs_FCI_Hartree']) for r in rows])
                if not np.isclose(np.mean(total**2), components['total_MSE_Hartree2'], rtol=1e-12, atol=1e-15):raise ValueError('Frozen total mismatch')
                attribution.append(dict(split=split, seed=int(seed), objective=objective, **components))
    paired=[]
    for split in ('test','ood'):
        for seed in (7301,7302,7303):
            group={r['objective']:r for r in attribution if r['split']==split and r['seed']==seed}
            energy,joint=group['energy_only'],group['energy_gradient']
            paired.append(dict(split=split,seed=seed,n_geometries=energy['n_geometries'],
                learning_RMSE_reduction_percent=100*(1-joint['learning_RMSE_Hartree']/energy['learning_RMSE_Hartree']),
                total_FCI_RMSE_reduction_percent=100*(1-joint['total_FCI_RMSE_Hartree']/energy['total_FCI_RMSE_Hartree']),
                energy_only_FCI_RMSE_Hartree=energy['total_FCI_RMSE_Hartree'],
                joint_FCI_RMSE_Hartree=joint['total_FCI_RMSE_Hartree']))
    tail=read_csv(ROOT/SOURCES[1]);by={(float(r['right_A']),int(r['intervals'])):r for r in tail}
    extrapolation=[richardson_binding(by[(48.,7200)],by[(48.,14400)]),
                   richardson_binding(by[(48.,14400)],by[(48.,28800)]),
                   richardson_binding(by[(96.,57600)],by[(96.,115200)])]
    tail_comparison=[]
    for r in tail:
        residual=float(r['eigen_residual_max_Hartree'])
        exact=float(r['last_full_line_exact_binding_Hartree'])
        same_state=int(r['bound_states_returned'])==17
        binding_error=abs(float(r['last_state_binding_Hartree'])-exact) if same_state else None
        tail_comparison.append(dict(right_A=float(r['right_A']),intervals=int(r['intervals']),
            states=int(r['bound_states_returned']),eigen_residual_Hartree=residual,
            same_state_binding_absolute_error_Hartree=binding_error,
            error_to_matrix_residual_ratio=binding_error/residual if same_state else None,
            F4000K_error_Hartree=float(r['F4000K_numerical_minus_full_line_Morse_Hartree'])))
    write_csv('mse_attribution.csv',attribution);write_csv('paired_gain_transfer.csv',paired)
    write_csv('tail_richardson.csv',extrapolation);write_csv('residual_observable_comparison.csv',tail_comparison)
    summary={'analysis_type':'Post hoc arithmetic on frozen data, not new electronic calculations, training or eigensolves.',
        'source_commit':'1ad05c243155a9f64b18e0d58c665424fde919a0',
        'counts':{'attribution_groups':12,'paired_comparisons':6,'richardson_estimates':3,'new_quantum_jobs':0,'new_training_runs':0,'new_eigensolves':0},
        'MSE_identity':'MSE(total)=MSE(learning)+MSE(reference bias)+2*mean(learning*reference bias)',
        'matched_geometries':{split:next(r['n_geometries'] for r in attribution if r['split']==split) for split in ('test','ood')},
        'max_MSE_identity_residual_Hartree2':max(abs(r['MSE_identity_residual_Hartree2']) for r in attribution),
        'limitations':['Post hoc analysis; not an independent replication or a preregistered benchmark.',
        'Only 8 matched test and 5 matched OOD geometries; OOD here differs from the original 9-point OOD set.',
        'Three paired seeds characterize these runs; no population-level significance or confidence interval is claimed.',
        'The MSE cross term may be negative, so component percentages are not variance fractions.',
        'Richardson estimates assume the leading second-order regime at fixed boundary. The full-line Morse value is an analytic model comparator, not experiment.',
        'Small matrix residual certifies the discrete operator, not the continuum boundary problem.'],
        'inputs_sha256':{name:sha(ROOT/name) for name in SOURCES},
        'script_sha256':sha(Path(__file__)),
        'outputs_sha256':{name:sha(OUT/name) for name in ['mse_attribution.csv','paired_gain_transfer.csv','tail_richardson.csv','residual_observable_comparison.csv']}}
    (OUT/'analysis_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'paired':paired,'richardson':extrapolation},indent=2))

if __name__=='__main__':main()
