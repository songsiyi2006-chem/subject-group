"""Archive the source failure and a single species-table repair; no installation.

The repaired random neural potential has no physical energy/charge calibration.
Changing the embedding table also changes subsequent random initialization.
"""
import argparse
from datetime import datetime, timezone
import difflib
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import runpy
import subprocess
import sys
import time

BASE = Path(__file__).resolve().parents[1]
ORIGINAL = BASE / 'source/synthapore_engine.py'
COMPAT = BASE / 'source/synthapore_engine_compat.py'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8', newline='\n')


def prepare():
    old = ORIGINAL.read_bytes().decode('utf-8')
    anchor = 'mlip = NeuralInteratomicPotential(num_species=10, hidden_dim=48, num_layers=3)'
    assert old.count(anchor) == 1
    new = old.replace(anchor, anchor.replace('num_species=10', 'num_species=30'))
    COMPAT.write_bytes(new.encode('utf-8'))
    patch = ''.join(difflib.unified_diff(old.splitlines(True), new.splitlines(True), fromfile='original', tofile='species-table repair'))
    (BASE/'source/compatibility.patch').write_text(patch, encoding='utf-8', newline='\n')
    dump(BASE/'source/source_record.json', {
        'specification_sha256':sha(BASE/'source/specification.md'),
        'original_python_sha256':sha(ORIGINAL), 'compatibility_python_sha256':sha(COMPAT),
        'repairs':['Main model species-table length 10 -> 30 to admit direct atomic-number index 29'],
        'scope':'Single execution repair; no training, physical calibration, mathematical or scientific correction',
        'random_initialization_caveat':'Embedding expansion consumes additional random draws and changes downstream random weights',
        'python':sys.version.split()[0],
        'packages':{name:importlib.metadata.version(name) for name in ['numpy','scipy','pandas','torch','matplotlib']},
        'physical_instrument_connected':False,
        'log_policy':'Local repository and runtime prefixes sanitized; errors and printed claims retained',
    })


def child(mode):
    import torch
    torch.set_num_threads(1)
    values = runpy.run_path(str(ORIGINAL if mode=='original' else COMPAT),run_name='__main__')
    capture = {name:values[name] for name in ['md_results','neb_results','pop_results']}
    for name in ['atomic_numbers','masses','r_init','r_final','edge_index','e_field']:
        capture[name] = values[name].detach().cpu().tolist()
    capture['model_parameter_count'] = sum(p.numel() for p in values['mlip'].parameters())
    capture['evidence'] = 'Supplied random-weight demonstration; source field names do not validate chemical claims'
    torch.save(values['mlip'].state_dict(), 'untrained_source_weights.pt')
    dump(Path('captured_results.json'),capture)


def sanitize(value):
    for path,label in [(str(BASE.parent),'<repository>'),(str(Path(sys.executable).parent.parent),'<python-runtime>')]:
        value = value.replace(path,label).replace(path.replace('\\','/'),label)
    return re.sub(r'[A-Za-z]:[\\/]Users[\\/][^\r\n\"\']+', '<local-runtime-path>', value)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--child', choices=['original','compatibility'])
    parser.add_argument('--stage', choices=['both','original','compatibility'], default='both')
    parser.add_argument('--timeout-seconds',type=int,default=600)
    args = parser.parse_args()
    if args.child:
        child(args.child)
        return 0
    if args.timeout_seconds <= 0:
        parser.error('timeout must be positive')
    prepare()
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONIOENCODING='utf-8',MPLBACKEND='Agg')
    results = []
    for mode in ['original','compatibility'] if args.stage=='both' else [args.stage]:
        folder = BASE/'results'/mode
        folder.mkdir(parents=True,exist_ok=True)
        started = datetime.now(timezone.utc).isoformat()
        tic = time.perf_counter()
        try:
            result = subprocess.run([sys.executable,str(Path(__file__).resolve()),'--child',mode],cwd=folder,env=env,
                     capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=args.timeout_seconds)
            code,stdout,stderr = result.returncode,result.stdout,result.stderr
        except subprocess.TimeoutExpired as exc:
            def decode(s):
                return s.decode('utf-8',errors='replace') if isinstance(s,bytes) else (s or '')
            code,stdout,stderr = 124,decode(exc.stdout),decode(exc.stderr)
        (folder/'stdout.log').write_text(sanitize(stdout),encoding='utf-8',newline='\n')
        (folder/'stderr.log').write_text(sanitize(stderr),encoding='utf-8',newline='\n')
        record = dict(mode=mode,exit_code=code,started_utc=started,elapsed_seconds=time.perf_counter()-tic,
                      source_sha256=sha(ORIGINAL if mode=='original' else COMPAT),
                      command='python synthapore/scripts/run_source.py --child '+mode,
                      physical_instrument_connected=False,
                      role='Retained source demonstration; printed claims are not endorsements',
                      output_sha256={p.name:sha(p) for p in sorted(folder.iterdir()) if p.is_file() and p.name!='execution.json'})
        dump(folder/'execution.json',record)
        print(json.dumps({k:record[k] for k in ['mode','exit_code','elapsed_seconds']}),flush=True)
        results.append(record)
    relevant = [r for r in results if r['mode']=='compatibility'] or results
    return 0 if all(r['exit_code']==0 for r in relevant) else 1


if __name__=='__main__':
    sys.exit(main())
