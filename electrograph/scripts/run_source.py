"""Run the immutable source and an explicitly recorded API-only repair.

Original scientific claims and numerical choices remain unchanged in the repair.
Public logs replace local filesystem prefixes; no package is installed here.
"""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import difflib
import hashlib
import importlib.metadata
import json
import os
import re
import runpy
import subprocess
import sys
import time

BASE = Path(__file__).resolve().parents[1]
REPO = BASE.parent
ORIGINAL = BASE / 'source/electrograph_kmc_core.py'
COMPAT = BASE / 'source/electrograph_kmc_core_compat.py'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n',
                    encoding='utf-8', newline='\n')


def sanitize(value):
    for path, name in [(str(REPO), '<repository>'),
                       (str(Path(sys.executable).parent.parent), '<python-runtime>')]:
        value = value.replace(path, name).replace(path.replace('\\', '/'), name)
    return re.sub(r'[A-Za-z]:[\\/]Users[\\/][^\r\n\"\']+', '<local-runtime-path>', value)


def prepare_compatibility():
    original = ORIGINAL.read_bytes()
    newline = '\r\n' if b'\r\n' in original else '\n'
    text = original.decode('utf-8')
    repairs = [
        ('Remove nonexistent unused AROMATIC hybridization enum',
         'Chem.rdchem.HybridizationType.SP3,' + newline + '                      Chem.rdchem.HybridizationType.AROMATIC]',
         'Chem.rdchem.HybridizationType.SP3]'),
        ('Use ETKDG EmbedParameters pruning field and supported overload',
         '        cids = AllChem.EmbedMultipleConfs(' + newline +
         '            mol, numConfs=self.num_confs,' + newline +
         '            params=AllChem.ETKDGv3(),' + newline +
         '            pruneRmsdThreshold=self.rmsd_thresh' + newline + '        )',
         '        params = AllChem.ETKDGv3()' + newline +
         '        params.pruneRmsThresh = self.rmsd_thresh' + newline +
         '        cids = AllChem.EmbedMultipleConfs(mol, self.num_confs, params)'),
        ('Correct FreeSASA function and conformer keyword names',
         'rdFreeSASA.calcSASA(mol, radii_table, confId=cid)',
         'rdFreeSASA.CalcSASA(mol, radii_table, confIdx=cid)'),
        ('Correct force-field minimizer keyword', 'ff.Minimize(maxIters=400)', 'ff.Minimize(maxIts=400)'),
    ]
    for label, before, after in repairs:
        if text.count(before) != 1:
            raise ValueError('Source repair anchor changed: ' + label)
        text = text.replace(before, after)
    COMPAT.write_bytes(text.encode('utf-8'))
    diff = ''.join(difflib.unified_diff(original.decode('utf-8').splitlines(True),
                    text.splitlines(True), fromfile='original', tofile='API-only compatibility'))
    (BASE / 'source/compatibility.patch').write_text(diff, encoding='utf-8', newline='\n')
    versions = {p: importlib.metadata.version(p) for p in ['numpy', 'scipy', 'pandas', 'matplotlib', 'torch', 'rdkit']}
    record = dict(specification_sha256=digest(BASE / 'source/specification.md'),
                  original_python_sha256=digest(ORIGINAL), compatibility_python_sha256=digest(COMPAT),
                  repairs=[r[0] for r in repairs], python=sys.version.split()[0], packages=versions,
                  compatibility_scope='Four API repairs only; random weights, zero classified radii, source kinetics, hardcoded curve and CGR claims preserved.',
                  log_policy='Local repository and runtime path prefixes sanitized in published logs.',
                  conformer_seed_policy='Original and compatibility scripts leave ETKDG randomSeed at its default; the reviewed structure module records explicit seeds.')
    attempt = BASE / 'source/electrograph_kmc_core_api_attempt_1.py'
    if attempt.exists():
        record['archived_failed_intermediate_source_sha256'] = digest(attempt)
    write_json(BASE / 'source/source_record.json', record)


def child(mode):
    source = ORIGINAL if mode == 'original' else COMPAT
    values = runpy.run_path(str(source), run_name='__main__')
    names = ['conf_results', 'kmc_results', 'cgr_results', 'ranked_subs', 'substrate_library', 'target_smiles', 'rxn_test']
    write_json(Path.cwd() / 'captured_results.json', {name: values[name] for name in names})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--child', choices=['original', 'compatibility'])
    parser.add_argument('--stage', choices=['both', 'original', 'compatibility'], default='both')
    parser.add_argument('--timeout-seconds', type=int, default=600)
    args = parser.parse_args()
    if args.child:
        child(args.child)
        return 0
    if args.timeout_seconds <= 0:
        parser.error('timeout must be positive')
    prepare_compatibility()
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               PYTHONIOENCODING='utf-8', MPLBACKEND='Agg')
    statuses = []
    for mode in (['original', 'compatibility'] if args.stage == 'both' else [args.stage]):
        folder = BASE / 'results' / mode
        folder.mkdir(parents=True, exist_ok=True)
        started = datetime.now(timezone.utc).isoformat()
        start = time.perf_counter()
        timed_out = False
        try:
            run = subprocess.run([sys.executable, str(Path(__file__).resolve()), '--child', mode],
                cwd=folder, env=env, capture_output=True, text=True, encoding='utf-8',
                errors='replace', timeout=args.timeout_seconds)
            code, stdout, stderr = run.returncode, run.stdout, run.stderr
        except subprocess.TimeoutExpired as exc:
            timed_out, code = True, 124
            def decoded(s):
                return s.decode('utf-8', errors='replace') if isinstance(s, bytes) else (s or '')
            stdout, stderr = decoded(exc.stdout), decoded(exc.stderr)
        (folder / 'stdout.log').write_text(sanitize(stdout), encoding='utf-8', newline='\n')
        (folder / 'stderr.log').write_text(sanitize(stderr), encoding='utf-8', newline='\n')
        record = dict(mode=mode, started_utc=started, elapsed_seconds=time.perf_counter()-start,
            exit_code=code, timed_out=timed_out, source_sha256=digest(ORIGINAL if mode == 'original' else COMPAT),
            command='python electrograph/scripts/run_source.py --child ' + mode,
            working_directory='electrograph/results/' + mode, physical_instrument_connected=False,
            role='Captured source execution; printed scientific claims are not endorsements')
        write_json(folder / 'execution.json', record)
        statuses.append(record)
        print(json.dumps({k: record[k] for k in ['mode', 'exit_code', 'timed_out', 'elapsed_seconds']}), flush=True)
    # Original failure is preserved and reported; success here requires a successful compatibility run.
    relevant = [r for r in statuses if r['mode'] == 'compatibility'] or statuses
    return 0 if all(r['exit_code'] == 0 for r in relevant) else 1


if __name__ == '__main__':
    sys.exit(main())
