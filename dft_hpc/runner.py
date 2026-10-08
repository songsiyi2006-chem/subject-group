"""Explicitly invoked Psi4 single-point runner; no automatic submission/retry."""
import argparse
import json
import math
import os
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path

try:
    from .workflow import digest, read_xyz, validate_state, write_json, validate_campaign
except ImportError:
    from workflow import digest, read_xyz, validate_state, write_json, validate_campaign


def run(campaign_path, index):
    campaign_path = Path(campaign_path).resolve()
    root = campaign_path.parent
    config = json.loads(campaign_path.read_text())
    validate_campaign(config, root)
    if not 0 <= index < len(config['jobs']):
        raise ValueError('Invalid array index')
    for name, expected in config['code_sha256'].items():
        if digest(root/name) != expected:
            raise ValueError('Prepared runner has changed; prepare a new campaign')
    job = config['jobs'][index]
    xyz = (root/job['geometry']).resolve()
    if not xyz.is_relative_to(root) or digest(xyz) != job['geometry_sha256']:
        raise ValueError('Geometry identity mismatch')
    atoms = read_xyz(xyz)
    if validate_state(atoms, job['charge'], job['multiplicity']) != job['electrons']:
        raise ValueError('Electron count mismatch')
    if job['task'] != 'single_point' or job['environment'] != 'gas' or job['energy_kind'] != 'electronic':
        raise ValueError('Runner supports gas electronic single points only')
    attempt = root/'runs'/job['id']/uuid.uuid4().hex
    attempt.mkdir(parents=True, exist_ok=False)
    result = dict(job=job, code_sha256=config['code_sha256'], status='started',
                  geometry_sha256=digest(xyz), started_utc=datetime.now(timezone.utc).isoformat(),
                  scf_converged=False, wavefunction_stability='not_checked',
                  optimization='not_performed', frequencies='not_performed')
    write_json(attempt/'result.json', result)
    original = Path.cwd()
    try:
        import psi4
        os.chdir(attempt)
        psi4.core.set_output_file(str(attempt/'engine.out'), False)
        psi4.core.IOManager.shared_object().set_default_path(str(attempt))
        psi4.set_num_threads(config['resources']['cores'])
        psi4.set_memory(f"{int(config['resources']['memory_mb'] * 0.8)} MB")
        geometry = '\n'.join(' '.join(map(str, a)) for a in atoms)
        mol = psi4.geometry(f"{job['charge']} {job['multiplicity']}\n{geometry}\nunits angstrom\nsymmetry c1\nno_reorient\nno_com")
        options = dict(job['options'], basis=job['basis'],
                       reference='rks' if job['multiplicity'] == 1 else 'uks')
        psi4.set_options(options)
        energy, wfn = psi4.energy(job['method'], molecule=mol, return_wfn=True)
        import numpy as np
        ca = wfn.Ca().np[:, :wfn.nalpha()]
        cb = wfn.Cb().np[:, :wfn.nbeta()]
        sz = (wfn.nalpha()-wfn.nbeta())/2
        spin_squared = sz*(sz+1)+wfn.nbeta()-float(np.sum((ca.T @ wfn.S().np @ cb)**2))
        if not math.isfinite(float(energy)) or not math.isfinite(spin_squared):
            raise ValueError('Nonfinite engine result')
        result.update(status='completed', scf_converged=True, energy_hartree=float(energy),
                      engine_version=psi4.__version__, nalpha=wfn.nalpha(), nbeta=wfn.nbeta(),
                      basis_functions=wfn.nso(), spin_squared=spin_squared,
                      spin_diagnostic='reported_not_used_as_universal_acceptance_threshold')
    except Exception as exc:
        result.update(status='failed', error_type=type(exc).__name__)
        print(f'Calculation failed: {type(exc).__name__}', file=sys.stderr)
    finally:
        if 'psi4' in locals():
            try:
                psi4.core.close_outfile()
            except Exception as exc:
                result.update(status='failed', log_close_error_type=type(exc).__name__)
        os.chdir(original)
        result['finished_utc'] = datetime.now(timezone.utc).isoformat()
        if (attempt/'engine.out').exists():
            result['engine_output_sha256'] = digest(attempt/'engine.out')
        write_json(attempt/'result.json', result)
    return 0 if result['status'] == 'completed' else 1


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--campaign', required=True)
    p.add_argument('--index', type=int, required=True)
    args = p.parse_args()
    sys.exit(run(args.campaign, args.index))
