"""Run archived-source, audit and reviewed modules in a fixed dependency order.

Re-execution replaces generated numerical files. Publication inspection remains
a separate step: regenerated figure bytes require a fresh figure QA record.
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
STAGES = {
    'source': 'run_source.py',
    'audit': 'audit_source.py',
    'learning': 'graph_learning.py',
    'structure': 'structure_reviewed.py',
    'kinetics': 'kinetics_reviewed.py',
    'figures': 'plot_reviewed.py',
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modules', nargs='+', choices=list(STAGES), default=list(STAGES))
    parser.add_argument('--timeout-seconds', type=int, default=1800)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    if args.timeout_seconds <= 0:
        parser.error('--timeout-seconds must be positive')
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONIOENCODING='utf-8')
    for name, script in STAGES.items():
        if name not in args.modules:
            continue
        relative = 'electrograph/scripts/' + script
        print(('DRY RUN: ' if args.dry_run else 'RUN: ') + 'python ' + relative, flush=True)
        if args.dry_run:
            continue
        try:
            subprocess.run([sys.executable, str(ROOT / relative)], cwd=ROOT,
                           env=env, check=True, timeout=args.timeout_seconds)
        except subprocess.TimeoutExpired:
            print(name + ' exceeded the per-stage timeout', file=sys.stderr)
            return 124
        except subprocess.CalledProcessError as exc:
            print(name + ' failed; subsequent stages were not run', file=sys.stderr)
            return exc.returncode or 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
