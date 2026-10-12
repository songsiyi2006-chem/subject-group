"""Paired-seed summaries of cached optimization campaigns (standard library)."""
import argparse
import csv
import hashlib
import json
import math
import random
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'electratwin/results/benchmark_extension/budget_prefixes.csv'
METHODS = ('gp_mc_ehvi', 'random', 'maximin_spacefill')
ROLE = 'cached_deterministic_PDE_pool_evaluation_use'
TIE_TOLERANCE = 1e-12


def finite(values):
    values = list(values)
    if not values or any(type(v) not in (int, float) or not math.isfinite(v) for v in values):
        raise ValueError('Expected nonempty finite numeric values')
    return values


def percentile(values, probability):
    values = sorted(finite(values))
    if not 0 <= probability <= 1:
        raise ValueError('Invalid percentile probability')
    position = (len(values)-1)*probability
    lo, hi = math.floor(position), math.ceil(position)
    return values[lo]+(values[hi]-values[lo])*(position-lo)


def exact_sign_p(wins, losses):
    """Two-sided exact Binomial(n, .5), excluding predefined numerical ties."""
    if any(type(x) is not int or x < 0 for x in [wins, losses]):
        raise ValueError('Sign counts must be nonnegative integers')
    n = wins+losses
    if n == 0:
        return 1.0
    return min(1.0, 2*sum(math.comb(n, k) for k in range(min(wins, losses)+1))/2**n)


def holm_adjust(pvalues):
    pvalues = finite(pvalues)
    if any(p < 0 or p > 1 for p in pvalues):
        raise ValueError('P-values outside [0,1]')
    order = sorted(range(len(pvalues)), key=lambda i: pvalues[i])
    adjusted = [0.0]*len(pvalues)
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, min(1.0, (len(order)-rank)*pvalues[index]))
        adjusted[index] = running
    return adjusted


def paired_summary(left, right, bootstrap_samples=4000, random_seed=20261012):
    """Resample whole seed pairs, never budget rows as independent samples."""
    if set(left) != set(right) or len(left) < 2:
        raise ValueError('At least two exactly matched seeds required')
    if type(bootstrap_samples) is not int or not 100 <= bootstrap_samples <= 100000:
        raise ValueError('Bootstrap samples must be an integer in [100,100000]')
    if type(random_seed) is not int:
        raise ValueError('Random seed must be an integer')
    keys = sorted(left)
    a, b = finite(left[k] for k in keys), finite(right[k] for k in keys)
    differences = finite(x-y for x, y in zip(a, b))
    n = len(keys)
    mean = statistics.mean(differences)
    rng = random.Random(random_seed)
    resampled = [statistics.mean(differences[rng.randrange(n)] for _ in range(n))
                 for _ in range(bootstrap_samples)]
    wins = sum(d > TIE_TOLERANCE for d in differences)
    losses = sum(d < -TIE_TOLERANCE for d in differences)
    leave_one_out = [(sum(differences)-d)/(n-1) for d in differences]
    return dict(n_pairs=n, left_mean=statistics.mean(a), right_mean=statistics.mean(b),
                mean_difference=mean, median_difference=statistics.median(differences),
                sd_difference=statistics.stdev(differences), wins=wins, losses=losses,
                ties=n-wins-losses, sign_p_raw=exact_sign_p(wins, losses),
                bootstrap_mean_ci95=[percentile(resampled, .025), percentile(resampled, .975)],
                bootstrap_samples=bootstrap_samples, bootstrap_seed=random_seed,
                leave_one_out_mean_min=min(leave_one_out), leave_one_out_mean_max=max(leave_one_out),
                max_absolute_leave_one_out_shift=max(abs(x-mean) for x in leave_one_out))


def normalized_auc(budgets, values):
    budgets, values = finite(budgets), finite(values)
    if len(budgets) != len(values) or len(budgets) < 2 or any(b <= a for a, b in zip(budgets, budgets[1:])):
        raise ValueError('AUC requires equal lengths and strictly increasing budgets')
    return sum((b-a)*(x+y)/2 for a, b, x, y in zip(budgets, budgets[1:], values, values[1:]))/(budgets[-1]-budgets[0])


def load_panel(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f))
    panel = {m: {} for m in METHODS}
    for row in rows:
        method = row['method']
        seed, budget = int(row['seed']), int(row['evaluation'])
        value = float(row['full_pool_HV_fraction'])
        if method not in panel or row['evidence_role'] != ROLE:
            raise ValueError('Unexpected method or evidence role')
        if seed < 0 or budget <= 0 or not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError('Invalid seed, budget or fraction')
        curve = panel[method].setdefault(seed, {})
        if budget in curve:
            raise ValueError('Duplicate method/seed/budget')
        curve[budget] = value
    seeds = sorted(panel[METHODS[0]])
    if len(seeds) < 2 or any(set(panel[m]) != set(seeds) for m in METHODS):
        raise ValueError('Unequal or insufficient seed sets')
    budgets = sorted(panel[METHODS[0]][seeds[0]])
    if len(budgets) < 2:
        raise ValueError('At least two recorded budgets required')
    for method in METHODS:
        for seed in seeds:
            curve = panel[method][seed]
            if sorted(curve) != budgets:
                raise ValueError('Unbalanced budget panel')
            if any(curve[b] < curve[a]-TIE_TOLERANCE for a, b in zip(budgets, budgets[1:])):
                raise ValueError('Cumulative hypervolume decreases')
    return panel, seeds, budgets


def analyze(path=SOURCE, bootstrap_samples=4000, random_seed=20261012, threshold=.99):
    if type(threshold) not in (int, float) or not math.isfinite(threshold) or not 0 < threshold <= 1:
        raise ValueError('Threshold must be in (0,1]')
    panel, seeds, budgets = load_panel(path)
    curve_metrics, attainment = [], []
    aucs = {m: {} for m in METHODS}
    for method in METHODS:
        for seed in seeds:
            curve = panel[method][seed]
            auc = normalized_auc(budgets, [curve[b] for b in budgets])
            aucs[method][seed] = auc
            reached = next((b for b in budgets if curve[b] >= threshold), None)
            index = budgets.index(reached) if reached is not None else None
            curve_metrics.append(dict(method=method, seed=seed, normalized_auc=auc,
                final_fraction=curve[budgets[-1]], first_recorded_budget_reaching_threshold=reached,
                crossing_lower_bound=budgets[index-1] if index is not None and index > 0 else (budgets[-1] if reached is None else None),
                crossing_upper_bound=reached,
                crossing_status='right_censored' if reached is None else ('left_censored' if index == 0 else 'interval_censored')))
        attainment.append(dict(method=method, total_seeds=len(seeds), threshold=threshold,
            counts_by_recorded_budget={str(b): sum(panel[method][s][b] >= threshold for s in seeds) for b in budgets},
            not_reached_by_final_budget=sum(panel[method][s][budgets[-1]] < threshold for s in seeds)))
    comparisons = []
    for other in METHODS[1:]:
        for budget in budgets:
            comparisons.append(dict(left=METHODS[0], right=other, endpoint=f'budget_{budget}',
                **paired_summary({s: panel[METHODS[0]][s][budget] for s in seeds},
                                 {s: panel[other][s][budget] for s in seeds}, bootstrap_samples, random_seed)))
        comparisons.append(dict(left=METHODS[0], right=other, endpoint='normalized_auc',
            **paired_summary(aucs[METHODS[0]], aucs[other], bootstrap_samples, random_seed)))
    for row, p in zip(comparisons, holm_adjust([r['sign_p_raw'] for r in comparisons])):
        row['sign_p_holm_all_endpoints'] = p
    report = dict(schema_version=1, source_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
        analyzer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        evidence_role=ROLE, seeds=seeds, budgets=budgets, comparisons=comparisons,
        curve_metrics=curve_metrics, attainment=attainment,
        interpretation=dict(resampling_unit='Whole matched seed pairs, conditional on one fixed cached candidate pool.',
            bootstrap='Pointwise percentile intervals, not simultaneous intervals; exchangeable seed replicates assumed.',
            sign_test='Two-sided probability of a positive non-tied difference equals 0.5; not a test of the mean difference.',
            multiplicity=f'Holm adjustment across all {len(comparisons)} archived exploratory endpoint comparisons.',
            auc='Piecewise-linear trapezoidal average over the observed budget range only; intermediate values are interpolated.',
            crossing='Recorded-prefix bounds only; never assign final budget to an unreached target.',
            scope='Post-hoc algorithm-seed sensitivity; no independent experiments or physical-model validation.'),
        new_solver_runs=0)
    return report


def write_report(destination, report):
    dest = Path(destination)
    dest.mkdir(parents=True, exist_ok=False)
    (dest/'analysis.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    for name, rows in [('curve_metrics', report['curve_metrics']), ('comparisons', report['comparisons'])]:
        with (dest/(name+'.csv')).open('w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=SOURCE)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--bootstrap-samples', type=int, default=4000)
    parser.add_argument('--seed', type=int, default=20261012)
    args = parser.parse_args()
    write_report(args.output, analyze(args.input, args.bootstrap_samples, args.seed))
