"""Execute the unchanged supplied script and capture its actual random weights.

Source scientific labels are retained for audit, not endorsed. Partial state is
captured on failure; no mathematical or compatibility repair is made here.
"""
import argparse
from datetime import datetime, timezone
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
SOURCE = BASE/'source/quantum_egnn_neb_engine.py'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8', newline='\n')


def capture(values, complete):
    import numpy as np
    import torch
    saved = {'execution_complete': complete, 'evidence': 'Unchanged supplied script; untrained neural potential and source claims retained for audit'}
    for key in ['elements', 'eht_results', 'neb_results', 'rrho_results', 'final_payload']:
        if key in values:
            saved[key] = values[key]
    for key in ['z_atomic_numbers', 'R_reactant', 'R_product', 'ts_coords', 'e_val', 'f_val']:
        if key in values:
            value = values[key]
            saved[key] = value.detach().cpu().tolist() if isinstance(value, torch.Tensor) else np.asarray(value).tolist()
    if 'egnn_model' in values:
        model = values['egnn_model']
        torch.save(model.state_dict(), 'untrained_source_weights.pt')
        saved['model_parameter_count'] = sum(p.numel() for p in model.parameters())
    dump(Path('captured_results.json'), saved)


def child():
    import torch
    torch.set_num_threads(1)
    try:
        values = runpy.run_path(str(SOURCE), run_name='__main__')
    except Exception as exc:
        frame = exc.__traceback__
        while frame is not None:
            if Path(frame.tb_frame.f_code.co_filename) == SOURCE:
                capture(frame.tb_frame.f_globals, False)
                break
            frame = frame.tb_next
        raise
    capture(values, True)


def sanitize(value):
    for path, label in [(str(BASE.parent), '<repository>'), (str(Path(sys.executable).parent.parent), '<python-runtime>')]:
        value = value.replace(path, label).replace(path.replace('\\', '/'), label)
    return re.sub(r'[A-Za-z]:[\\/]Users[\\/][^\r\n\"\']+', '<local-runtime-path>', value)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--child', action='store_true')
    parser.add_argument('--timeout-seconds', type=int, default=600)
    args = parser.parse_args()
    if args.child:
        child()
        return 0
    if args.timeout_seconds <= 0:
        parser.error('timeout must be positive')
    dump(BASE/'source/source_record.json', {
        'specification_sha256': sha(BASE/'source/specification.md'), 'original_python_sha256': sha(SOURCE),
        'repairs': [], 'python': sys.version.split()[0],
        'packages': {name: importlib.metadata.version(name) for name in ['numpy', 'scipy', 'pandas', 'torch', 'matplotlib', 'rdkit']},
        'scope': 'Unchanged source execution and capture; no instrument connection or environment installation',
        'log_policy': 'Local runtime and repository prefixes sanitized; errors and printed claims retained',
    })
    output = BASE/'results/original'
    output.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONIOENCODING='utf-8', MPLBACKEND='Agg')
    tic = time.perf_counter()
    started = datetime.now(timezone.utc).isoformat()
    try:
        result = subprocess.run([sys.executable, str(Path(__file__).resolve()), '--child'], cwd=output, env=env,
                                capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=args.timeout_seconds)
        code, stdout, stderr = result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired as exc:
        def decode(value):
            return value.decode('utf-8', errors='replace') if isinstance(value, bytes) else (value or '')
        code, stdout, stderr = 124, decode(exc.stdout), decode(exc.stderr)
    (output/'stdout.log').write_text(sanitize(stdout), encoding='utf-8', newline='\n')
    (output/'stderr.log').write_text(sanitize(stderr), encoding='utf-8', newline='\n')
    record = {'exit_code': code, 'started_utc': started, 'elapsed_seconds': time.perf_counter()-tic,
              'source_sha256': sha(SOURCE), 'command': 'python quantumequi/scripts/run_source.py --child',
              'physical_instrument_connected': False, 'source_claims_endorsed': False,
              'output_sha256': {p.relative_to(output).as_posix(): sha(p) for p in sorted(output.rglob('*'))
                                if p.is_file() and p.name != 'execution.json'}}
    dump(output/'execution.json', record)
    print(json.dumps({k: record[k] for k in ['exit_code', 'elapsed_seconds']}), flush=True)
    return code


if __name__ == '__main__':
    sys.exit(main())
