"""Recompute saved post-hoc analysis; no physics or scheduler calls."""
import csv
import io
import json
from pathlib import Path
from .paired_budget import analyze


def main():
    root = Path(__file__).resolve().parent/'results'
    report = json.loads((root/'analysis.json').read_text('utf-8'))
    expected = analyze()
    checks = {'deterministic_full_recomputation': report == expected,
              '64_paired_seeds': len(report['seeds']) == 64,
              '5_recorded_budgets': report['budgets'] == [5, 10, 15, 20, 25],
              '12_comparisons': len(report['comparisons']) == 12,
              '192_curves': len(report['curve_metrics']) == 192,
              'no_new_solver_runs': report['new_solver_runs'] == 0}
    for name, rows in [('curve_metrics', expected['curve_metrics']), ('comparisons', expected['comparisons'])]:
        stream = io.StringIO(newline='')
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
        checks[name+'_csv_exact'] = (root/(name+'.csv')).read_bytes() == stream.getvalue().encode('utf-8')
    print(json.dumps({'passed': all(checks.values()), 'checks': checks}))
    return 0 if all(checks.values()) else 1


if __name__ == '__main__':
    raise SystemExit(main())
