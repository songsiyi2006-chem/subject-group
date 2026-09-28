"""Portable sequential entry point for the three calculation extensions.

Default execution regenerates result files. Preserve the release archive first.
This launcher never connects instruments and never updates figure-QA records.
"""
from pathlib import Path
import argparse
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
MODULES = {
    'network': 'electratwin/scripts/reaction_network.py',
    'uncertainty': 'electratwin/scripts/uncertainty_analysis.py',
    'benchmark': 'electratwin/scripts/benchmark_extension.py',
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modules', nargs='+', choices=list(MODULES), default=list(MODULES))
    parser.add_argument('--dry-run', action='store_true', help='Show commands without changing outputs')
    parser.add_argument('--timeout-seconds', type=int, default=1800, help='Maximum time per calculation module')
    args = parser.parse_args()
    if args.timeout_seconds <= 0:
        parser.error('--timeout-seconds must be positive')
    names = list(dict.fromkeys(args.modules))
    for name in names:
        if not (ROOT / MODULES[name]).is_file():
            parser.error('Missing module: ' + MODULES[name])
    for name in names:
        print('python ' + MODULES[name], flush=True)
    if args.dry_run:
        return 0
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONIOENCODING='utf-8')
    for name in names:
        print('\nRunning ' + name + ' (numerical model only)', flush=True)
        try:
            result = subprocess.run([sys.executable, MODULES[name]], cwd=ROOT, env=env,
                                    timeout=args.timeout_seconds, check=False)
        except subprocess.TimeoutExpired:
            print('Stopped at timeout: ' + name + '; inspect partial outputs before reuse.', file=sys.stderr)
            return 124
        if result.returncode:
            print('Stopped after failed module: ' + name, file=sys.stderr)
            return result.returncode
    print('Calculation modules completed. Regenerated figures require a new visual review; archived QA is unchanged.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
