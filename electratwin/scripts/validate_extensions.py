"""Validate saved extension evidence and publication files without scientific reruns."""
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from urllib.parse import unquote
import xml.etree.ElementTree as ET
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
RESULT = ROOT / 'results/extension_validation.json'


def main():
    checks = []
    def check(name, ok, detail=None):
        item = dict(check=name, passed=bool(ok))
        if detail is not None:
            item['detail'] = str(detail).replace(str(ROOT), 'electratwin').replace(str(REPO), 'repository')
        checks.append(item)
    def read(path):
        return json.loads((ROOT / path).read_text(encoding='utf-8-sig'))
    def rows(path):
        with (ROOT / path).open(encoding='utf-8-sig', newline='') as stream:
            return list(csv.DictReader(stream))
    def digest(path):
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    def hashes(mapping, base, label):
        for path, expected in mapping.items():
            target = base / path
            check(label + ' ' + path, target.is_file() and digest(target) == expected)
    def guarded(name, function):
        try:
            function()
        except Exception as exc:
            check(name + ' readable and valid schema', False, f'{type(exc).__name__}: {exc}')
    def count(path, expected):
        data = rows(path)
        check(path + ' row count', len(data) == expected)
        return data
    def near(a, b):
        return math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-12)

    def network():
        n = read('results/reaction_network/study_summary.json')
        hashes(n['source_sha256'], REPO, 'network source')
        p = n['PDE_counts']
        check('network study PDE count 62 plus comparator 1', p['network_solves'] == 62 and p['frozen_comparison_solves'] == 1 and p['total_PDE_solves'] == 63)
        registry = count('results/reaction_network/solve_registry.csv', 62)
        check('network consecutive solve IDs', [int(r['network_solve_id']) for r in registry] == list(range(1, 63)))
        check('network study count partition', [p[k] for k in ['baseline', 'pool', 'joint_grid', 'limits', 'analytic_limits']] == [1, 49, 3, 5, 4])
        for name, size in [('parameter_pool', 49), ('grid_convergence', 3), ('limiting_cases', 5), ('analytic_sequential_limit', 4)]:
            count('results/reaction_network/' + name + '.csv', size)
        check('network numerical checks and status', n['all_network_cases_converged'] and all(r['converged'] == 'True' and r['nonnegative'] == 'True' for r in registry))
        check('network all saved balances', all(max(float(r[k]) for k in ['max_species_balance_relative_error', 'total_molar_balance_relative_error', 'charge_balance_relative_error']) < 1e-8 for r in registry))
        invalid = read('results/reaction_network/input_rejections.json')
        check('eight network inputs rejected before solving', len(invalid) == p['invalid_inputs_rejected_without_solving'] == 8 and all(r['rejected_before_PDE_solve'] for r in invalid))
        limits = rows('results/reaction_network/limiting_cases.csv')
        check('negative product-feed net FE retained', any(r['case'] == 'P_feed_overoxidation' and near(r['net_P_faradaic_efficiency_pct'], -100) for r in limits))
        s = n['baseline']
        check('baseline current branch sum', near(s['current_A'], sum(s['channel_current_A'].values())))
        check('baseline net-product FE independently reproduced', near(s['net_P_faradaic_efficiency_pct'], 100 * 2 * 96485.33212 * s['net_outlet_gain_mol_s']['P'] / s['current_A']))

    def uncertainty():
        u = read('results/uncertainty/summary.json')
        hashes(u['source_sha256'], REPO, 'uncertainty source')
        hashes(u['output_sha256'], ROOT / 'results/uncertainty', 'uncertainty output')
        sizes = {'design_evaluations': 1792, 'grid_check_evaluations': 16, 'grid_check_differences': 32, 'local_difference_evaluations': 21,
                 'local_derivatives': 40, 'base_AB_quantiles': 4, 'sobol_prefix_estimates': 60, 'sobol_resampling_intervals': 20, 'sobol_bootstrap_replicates': 10000}
        tables = {name: count('results/uncertainty/' + name + '.csv', size) for name, size in sizes.items()}
        c = u['counts']
        check('uncertainty 1792+16+21=1829 PDE count', [c[k] for k in ['main_PDE_solves', 'grid_check_PDE_solves', 'local_difference_PDE_solves', 'total_PDE_solves']] == [1792, 16, 21, 1829])
        check('bootstrap and prefixes use no extra PDE', c['additional_PDE_solves_for_prefix_or_bootstrap'] == 0 and c['bootstrap_output_rows'] == 10000)
        design = tables['design_evaluations']
        check('Sobol A/B/hybrid blocks each contain 256 rows', Counter(r['block'] for r in design) == {k: 256 for k in ['A', 'B', 'AB_0', 'AB_1', 'AB_2', 'AB_3', 'AB_4']})
        idx = {(r['block'], int(r['paired_row_index'])): r for r in design}
        columns = ['u_' + p['name'] for p in u['parameters']]
        check('hybrids replace exactly one matching B coordinate', all(idx[(f'AB_{i}', j)][key] == idx[('B' if k == i else 'A', j)][key] for i in range(5) for j in range(256) for k, key in enumerate(columns)))
        check('512 base evaluations define scenario quantiles', all(int(r['base_AB_count']) == 512 and 'not_experimental' in r['role'] for r in tables['base_AB_quantiles']))
        check('uncertainty statuses and balances', u['numerical_checks']['all_converged'] and all(r['converged'] == 'True' and float(r['material_balance_relative_error']) < 1e-8 for r in design))
        check('unclipped sensitivity and negative stability bounds retained', not u['jansen']['negative_or_above_one_estimates_clipped'] and all(r['indices_clipped'] == 'False' for r in tables['sobol_prefix_estimates']) and any(float(r['S_bootstrap_p025']) < 0 for r in tables['sobol_resampling_intervals']))
        check('QMC bootstrap limitations retained', 'not iid' in u['bootstrap']['interpretation'] and 'not a rigorous' in u['bootstrap']['interpretation'])

    def benchmark():
        b = read('results/benchmark_extension/summary.json')
        paths = {'pool': 'results/control/candidate_pool.csv', 'old_control': 'scripts/metrics_control.py', 'old_campaign': 'results/control/sequential_campaigns.csv', 'extension': 'scripts/benchmark_extension.py'}
        hashes({paths[k]: v for k, v in b['source_sha256'].items()}, ROOT, 'benchmark source')
        hashes(b['output_sha256'], ROOT / 'results/benchmark_extension', 'benchmark output')
        check('benchmark zero new PDE experimental and hardware runs', all(b[k] == 0 for k in ['new_PDE_solves', 'new_experimental_runs', 'hardware_connections']))
        c = b['campaign']
        check('benchmark declared 64 by 3 by 25', c['seeds'] == list(range(64)) and c['total_runs'] == 192 and c['cached_oracle_evaluation_uses'] == 4800 and c['budget_including_initial'] == 25 and c['shared_initial_count'] == 5)
        campaign = count('results/benchmark_extension/campaign_evaluations.csv', 4800)
        for name, size in [('budget_prefixes', 960), ('learning_curve_summary', 15), ('paired_HV_differences', 15), ('holdout_predictions', 4212), ('holdout_split_metrics', 156), ('holdout_grouped_metrics', 18)]:
            count('results/benchmark_extension/' + name + '.csv', size)
        grouped = defaultdict(list)
        for row in campaign:
            grouped[(int(row['seed']), row['method'])].append(row)
        methods = ['gp_mc_ehvi', 'random', 'maximin_spacefill']
        check('192 complete campaign identities', set(grouped) == {(seed, method) for seed in range(64) for method in methods})
        for seed in range(64):
            starts = []
            for method in methods:
                run = grouped[(seed, method)]
                check(f'campaign {seed} {method} budget uniqueness HV', [int(r['evaluation']) for r in run] == list(range(1, 26)) and len({r['candidate_index'] for r in run}) == 25 and all(float(y['hypervolume']) >= float(x['hypervolume']) - 1e-12 for x, y in zip(run, run[1:])))
                starts.append([r['candidate_index'] for r in run[:5]])
            check(f'campaign {seed} shared five initial points', starts[0] == starts[1] == starts[2])
        budget = read('results/benchmark_extension/oracle_budget_checks.json')
        check('192 observed oracle counts', len(budget) == 192 and all(r['passed'] and r['oracle_calls'] == r['unique_calls'] == r['requested_budget'] == 25 for r in budget))
        lookup = {(r['seed'], r['method'], r['evaluation']): r for r in campaign}
        old = rows('results/control/sequential_campaigns.csv')
        check('all original 240 saved prefixes independently match', len(old) == 240 and all(lookup[(r['seed'], r['method'], r['evaluation'])]['candidate_index'] == r['candidate_index'] and near(lookup[(r['seed'], r['method'], r['evaluation'])]['hypervolume'], r['hypervolume']) for r in old))
        splits = read('results/benchmark_extension/holdout_split_metadata.json')
        check('26 split identities and kinds', len(splits) == len({r['split_id'] for r in splits}) == 26 and Counter(r['kind'] for r in splits) == {'random_54_27': 20, 'block_flow': 3, 'block_eta': 3})
        pool = rows('results/control/candidate_pool.csv')
        for s in splits:
            train, test = set(s['train_ids']), set(s['test_ids'])
            check(s['split_id'] + ' training test separation', len(train) == 54 and len(test) == 27 and not train & test and train | test == set(range(81)) and not s['heldout_targets_seen_by_fit'])
            check(s['split_id'] + ' training-only normalization', all(near(s['input_min_from_training'][i], min(float(pool[j][col]) for j in train)) and near(s['input_max_from_training'][i], max(float(pool[j][col]) for j in train)) for i, col in enumerate(['flow_rate_uL_min', 'overpotential_V'])))
        test_sets = {r['split_id']: set(r['test_ids']) for r in splits}
        predictions = rows('results/benchmark_extension/holdout_predictions.csv')
        check('4212 unique heldout predictions never train members', len({(r['split_id'], r['model'], r['target'], r['candidate_index']) for r in predictions}) == 4212 and all(int(r['candidate_index']) in test_sets[r['split_id']] for r in predictions))
        groups = {(r['split_kind'], r['model'], r['target']): r for r in rows('results/benchmark_extension/holdout_grouped_metrics.csv')}
        check('negative blocked GP comparisons retained', all(float(groups[(kind, 'fixed_matern_gp', target)]['RMSE']) > float(groups[(kind, 'quadratic_regression', target)]['RMSE']) for kind in ['block_flow', 'block_eta'] for target in b['holdout']['targets']))
        check('eta GP latent interval undercoverage retained', all(float(groups[('block_eta', 'fixed_matern_gp', target)]['latent_interval_coverage_95pct']) < .8 for target in b['holdout']['targets']))
        paired = rows('results/benchmark_extension/paired_HV_differences.csv')
        check('budget15 GP maximin negative comparisons retained', any(r['budget'] == '15' and r['left_method'] == 'gp_mc_ehvi' and r['right_method'] == 'maximin_spacefill' and r['losses'] == '39' and float(r['lower_95pct_t']) < 0 < float(r['upper_95pct_t']) for r in paired))

    def figures():
        figures = read('results/extension_figure_manifest.json')['figures']
        check('16 distinct extension figure files', len(figures) == len({r['file'] for r in figures}) == 16 and Counter(r['format'] for r in figures) == {'png': 8, 'svg': 8})
        for r in figures:
            path = ROOT / 'reports/figures' / r['file']
            check('figure current hash ' + r['file'], path.is_file() and digest(path) == r['sha256'])
            hashes(r['sources'], ROOT / 'results', r['file'] + ' source')
            if r['format'] == 'png':
                with Image.open(path) as im:
                    check('PNG dimensions resolution ' + r['file'], list(im.size) == r['pixels'] == [2160, 1140] and all(abs(v - 300) < .1 for v in im.info.get('dpi', [0, 0])))
            else:
                xml = ET.parse(path); ns = {'s': 'http://www.w3.org/2000/svg'}
                nodes = xml.findall('.//s:text', ns)
                words = ' '.join(''.join(n.itertext()) for n in nodes)
                check('SVG editable pure vector ' + r['file'], len(nodes) > 10 and not xml.findall('.//s:image', ns) and r['editable_text'] and not r['embedded_raster'])
                check('SVG visible evidence label ' + r['file'], bool(re.search(r'HYPOTHETICAL|UNCALIBRATED|SCENARIO|SIMULATION ONLY|NUMERICAL|ASSUMED|MODEL HOLDOUT|假设|模拟|情景|模型留出', words)))
        qa = read('results/extension_figure_qa.json')
        check('visual QA matches all current PNG SVG bytes', qa['reviewed_figure_hashes'] == {r['file']: r['sha256'] for r in figures} and qa['png_files_visually_reviewed'] == 8 and qa['svg_files_rendered_and_visually_reviewed'] == 8)

    def reports():
        texts = {}
        for lang in ['english', 'chinese']:
            path = ROOT / 'reports' / f'extension_report_{lang}.md'
            text = path.read_text(encoding='utf-8'); texts[lang] = text
            check(lang + ' five main sections', re.findall(r'^## (\d+)\.', text, re.M) == ['1', '2', '3', '4', '5'])
            expected = [f'{s}.{i}' for s, count in [(1, 2), (2, 3), (3, 3), (4, 4), (5, 3)] for i in range(1, count + 1)]
            check(lang + ' fifteen ordered subsections', re.findall(r'^### (\d+\.\d+) ', text, re.M) == expected)
            check(lang + ' saved numerical anchors', all(v in text for v in ['1892', '4800', '4212', '11.3870659398', '0.1598292455', '0.992024355', '0.982781182', '0.989217617', '77.7778', '75.3086']))
            for number, base in [(4, 'Reaction_Network'), (5, 'Parameter_Sensitivity'), (6, 'Optimization_Extension'), (7, 'Heldout_Prediction')]:
                check(lang + f' figure {number} caption and two formats', ('**' + ('Figure ' if lang == 'english' else '图 ') + str(number)) in text and all(f'figures/Fig{number}_{base}_{lang}.{ext}' in text for ext in ['png', 'svg']))
            for link in re.findall(r'\]\(([^)]+)\)', text):
                if not link.startswith(('http://', 'https://', '#')):
                    target = path.parent / unquote(link.split('#', 1)[0].strip('<>'))
                    check(lang + ' local link ' + link, target.exists() or target.resolve() == RESULT.resolve())
            check(lang + ' no unfinished tokens and balanced fences', not re.search(r'@@[A-Z_]+@@', text) and text.count('```') % 2 == 0 and text.count('\\[') == text.count('\\]'))
        # A hyphen in an English compound such as Budget-25 is not a minus sign.
        numbers = lambda t: [re.findall(r'(?<![A-Za-z0-9.])[-−]?\d+(?:\.\d+)?(?:e[-+]?\d+)?', line) for line in t.splitlines() if line.startswith('|')]
        check('English Chinese numeric table parity', numbers(texts['english']) == numbers(texts['chinese']))
        check('English report no Chinese prose', not re.search(r'[\u4e00-\u9fff]', texts['english']))

    for name, function in [('network', network), ('uncertainty', uncertainty), ('benchmark', benchmark), ('figures', figures), ('reports', reports)]:
        guarded(name, function)
    result = dict(passed=all(c['passed'] for c in checks), checks=len(checks), failed=[c for c in checks if not c['passed']],
                  scope='Saved provenance, numerical accounting, negative results, report parity and publication QA records; no numerical reruns or physical validation.', items=checks)
    RESULT.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: result[k] for k in ['passed', 'checks', 'failed']}, indent=2, ensure_ascii=False))
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
