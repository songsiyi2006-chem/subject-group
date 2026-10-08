"""Standard-library-only preparation and pairing; never submits jobs."""
import argparse
import hashlib
import json
import math
from pathlib import Path

BASE = Path(__file__).resolve().parent
HARTREE_TO_EV = 27.211386245988
Z = dict(H=1, C=6, N=7, O=8, F=9, Cl=17, Br=35)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n',
                          encoding='utf-8', newline='\n')


def read_xyz(path):
    lines = Path(path).read_text('utf-8').splitlines()
    if len(lines) < 3:
        raise ValueError('XYZ header or atom rows missing')
    n = int(lines[0])
    if n < 1 or len(lines) != n + 2:
        raise ValueError('XYZ atom count mismatch')
    atoms = []
    for line in lines[2:]:
        fields = line.split()
        if len(fields) != 4 or fields[0] not in Z:
            raise ValueError('Unsupported atom or invalid XYZ row')
        coords = [float(x) for x in fields[1:]]
        if not all(math.isfinite(x) for x in coords):
            raise ValueError('Nonfinite coordinate')
        atoms.append((fields[0], *coords))
    for i, a in enumerate(atoms):
        for b in atoms[:i]:
            if math.dist(a[1:], b[1:]) < 0.2:
                raise ValueError('Overlapping nuclei')
    return atoms


def validate_state(atoms, charge, multiplicity):
    if type(charge) is not int or type(multiplicity) is not int:
        raise ValueError('Charge and multiplicity must be integers')
    electrons = sum(Z[a[0]] for a in atoms) - charge
    unpaired = multiplicity - 1
    if electrons < 1 or unpaired < 0 or unpaired > electrons or (electrons-unpaired) % 2:
        raise ValueError('Electron count / multiplicity mismatch')
    return electrons


def plan():
    jobs = []
    for source in json.loads((BASE/'geometry_sources.json').read_text()):
        xyz = BASE/source['file']
        if digest(xyz) != source['sha256']:
            raise ValueError('Archived geometry hash mismatch')
        atoms = read_xyz(xyz)
        for state, charge, mult in [('cation', 1, 1), ('radical', 0, 2)]:
            jobs.append(dict(id=f"{source['molecule_id']}_tzvp_{state}",
                molecule_id=source['molecule_id'], state=state, charge=charge,
                multiplicity=mult, electrons=validate_state(atoms, charge, mult),
                geometry=source['file'], geometry_sha256=source['sha256'],
                symbols=[a[0] for a in atoms], method='pbe0', basis='def2-tzvp',
                environment='gas', task='single_point', energy_kind='electronic',
                options=dict(scf_type='df', e_convergence=1e-9, d_convergence=1e-7,
                             maxiter=150, dft_radial_points=75, dft_spherical_points=302),
                source=source))
    return dict(schema_version=1, status='prepared_not_executed', jobs=jobs,
                deferred=[dict(id='quinazolinone_anions', reason='Original complete state and geometry inputs must be verified before preparation')])


def prepare(destination, cores=8, memory_mb=16000, hours=4, concurrency=2):
    for x in [cores, memory_mb, hours, concurrency]:
        if type(x) is not int or x < 1:
            raise ValueError('Resource limits must be positive integers')
    if cores > 64 or hours > 168 or concurrency > 16 or not 2000 <= memory_mb <= 512000:
        raise ValueError('Outside bounded preparation resource limits')
    data = plan()  # Validate everything before creating a delivery directory.
    dest = Path(destination)
    dest.mkdir(parents=True, exist_ok=False)
    (dest/'geometries').mkdir()
    for source in json.loads((BASE/'geometry_sources.json').read_text()):
        (dest/source['file']).write_bytes((BASE/source['file']).read_bytes())
    for filename in ['runner.py', 'workflow.py', 'geometry_sources.json']:
        (dest/filename).write_bytes((BASE/filename).read_bytes())
    data['resources'] = dict(cores=cores, memory_mb=memory_mb, hours=hours, concurrency=concurrency)
    data['code_sha256'] = {n: digest(dest/n) for n in ['runner.py', 'workflow.py']}
    write_json(dest/'campaign.json', data)
    script = f'''#!/bin/bash
#SBATCH --job-name=dft-supplement
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cores}
#SBATCH --mem={memory_mb}M
#SBATCH --time={hours}:00:00
#SBATCH --array=0-{len(data['jobs'])-1}%{concurrency}
#SBATCH --output=slurm-%A_%a.out
set -euo pipefail
cd "${{SLURM_SUBMIT_DIR:?Submit from the prepared campaign directory}}"
# Activate the site's Psi4 Python environment before sbatch.
# Add site account/partition settings on the sbatch command line.
export OMP_NUM_THREADS={cores}
export MKL_NUM_THREADS={cores}
export OPENBLAS_NUM_THREADS={cores}
python runner.py --campaign campaign.json --index "${{SLURM_ARRAY_TASK_ID}}"
'''
    (dest/'submit.slurm').write_text(script, encoding='utf-8', newline='\n')
    return data


def validate_campaign(campaign, root):
    """This version executes only the four reviewed, frozen task definitions."""
    if campaign.get('jobs') != plan()['jobs']:
        raise ValueError('Campaign differs from the reviewed four-job plan')
    r = campaign['resources']
    if any(type(r.get(k)) is not int for k in ['cores', 'memory_mb', 'hours', 'concurrency']):
        raise ValueError('Invalid resource schema')
    if not (1 <= r['cores'] <= 64 and 2000 <= r['memory_mb'] <= 512000
            and 1 <= r['hours'] <= 168 and 1 <= r['concurrency'] <= 16):
        raise ValueError('Resource bounds exceeded')
    for name in ['runner.py', 'workflow.py']:
        if digest(Path(root)/name) != campaign['code_sha256'].get(name):
            raise ValueError('Prepared code differs from manifest')
    for job in campaign['jobs']:
        if digest(Path(root)/job['geometry']) != job['geometry_sha256']:
            raise ValueError('Prepared geometry differs from manifest')


def protocol(job):
    keys = ['molecule_id', 'geometry_sha256', 'symbols', 'method', 'basis',
            'environment', 'task', 'energy_kind', 'options']
    return json.dumps({k: job[k] for k in keys}, sort_keys=True)


def assess(record, job, code_hashes):
    reasons = []
    if record.get('job') != job:
        reasons.append('job_identity_mismatch')
    if record.get('code_sha256') != code_hashes:
        reasons.append('runner_identity_mismatch')
    if record.get('status') != 'completed' or record.get('scf_converged') is not True:
        reasons.append('not_scf_completed')
    e = record.get('energy_hartree')
    if type(e) not in (int, float) or not math.isfinite(e):
        reasons.append('missing_or_nonfinite_energy')
    if not record.get('engine_version'):
        reasons.append('missing_engine_version')
    if record.get('geometry_sha256') != job['geometry_sha256']:
        reasons.append('geometry_mismatch')
    return reasons


def collect(campaign_path):
    campaign_path = Path(campaign_path)
    root = campaign_path.parent
    campaign = json.loads(campaign_path.read_text())
    validate_campaign(campaign, root)
    groups, audits = {}, []
    for job in campaign['jobs']:
        paths = sorted((root/'runs'/job['id']).glob('*/result.json'))
        eligible = []
        for path in paths:
            try:
                r = json.loads(path.read_text())
                if not isinstance(r, dict):
                    raise ValueError('Result must be an object')
                reasons = assess(r, job, campaign['code_sha256'])
                raw = path.parent/'engine.out'
                if not raw.is_file() or digest(raw) != r.get('engine_output_sha256'):
                    reasons.append('raw_log_missing_or_changed')
                if not reasons:
                    eligible.append(r)
            except (ValueError, KeyError, TypeError):
                reasons = ['invalid_result_record']
            audits.append(dict(job_id=job['id'], record=path.relative_to(root).as_posix(), reasons=reasons))
        status = 'available' if len(eligible) == 1 else ('ambiguous' if len(eligible) > 1 else 'missing_or_failed')
        groups.setdefault(protocol(job), []).append((job, eligible[0] if status == 'available' else None, status))
    pairs = []
    for key, entries in groups.items():
        entries.sort(key=lambda x: x[0]['charge'], reverse=True)
        pair = dict(protocol=json.loads(key), status='incomplete', energy_difference_hartree=None, energy_difference_eV=None)
        pair['states'] = [{'id': j['id'], 'status': s} for j, r, s in entries]
        if len(entries) == 2 and all(r is not None for j, r, s in entries):
            high, low = entries
            if high[0]['charge'] - low[0]['charge'] == 1 and high[1]['engine_version'] == low[1]['engine_version']:
                delta = low[1]['energy_hartree'] - high[1]['energy_hartree']
                pair.update(status='diagnostic_pair', energy_difference_hartree=delta,
                            energy_difference_eV=delta*HARTREE_TO_EV)
        pair['interpretation'] = 'Fixed-geometry electronic difference only; wavefunction stability, minimum, free energy and electrode potential not validated.'
        pairs.append(pair)
    return dict(schema_version=1, records=audits, pairs=pairs,
                accepted_chemical_predictions=0)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest='command', required=True)
    p = sub.add_parser('prepare'); p.add_argument('destination')
    for key, default in [('cores', 8), ('memory-mb', 16000), ('hours', 4), ('concurrency', 2)]:
        p.add_argument('--'+key, type=int, default=default)
    p = sub.add_parser('collect'); p.add_argument('campaign'); p.add_argument('--output', required=True)
    args = ap.parse_args()
    if args.command == 'prepare':
        data = prepare(args.destination, args.cores, args.memory_mb, args.hours, args.concurrency)
        print(json.dumps({'status': data['status'], 'jobs': len(data['jobs'])}))
    else:
        output = Path(args.output)
        if output.exists():
            raise FileExistsError('Choose a new report path to preserve previous audit')
        write_json(output, collect(args.campaign))


if __name__ == '__main__':
    main()
