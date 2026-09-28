"""Check saved numerical provenance and publication artifacts, not physical validity."""
from pathlib import Path
import csv
import hashlib
import json
import math
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]


def main():
    checks = []

    def check(name, ok):
        checks.append(dict(check=name, passed=bool(ok)))

    def read(name):
        return json.loads((ROOT / name).read_text(encoding='utf-8-sig'))

    def digest(name):
        return hashlib.sha256((ROOT / name).read_bytes()).hexdigest()

    def rows(name):
        with (ROOT / name).open(encoding='utf-8', newline='') as f:
            return list(csv.DictReader(f))

    extensions = subprocess.run(
        [sys.executable, str(ROOT / 'scripts/validate_extensions.py')],
        capture_output=True, text=True, encoding='utf-8', timeout=120,
    )
    check('saved extension provenance and reports', extensions.returncode == 0)
    if extensions.returncode:
        print(extensions.stdout + extensions.stderr)

    source = read('source/source_record.json')
    original = read('results/original/execution.json')
    compat = read('results/compatibility/execution.json')
    audit = read('results/source_audit.json')
    check('archived specification bytes', digest('source/specification.md') == source['specification_sha256'])
    check('original source and execution bytes', digest('source/electratwin_core.py') == source['original_python_sha256'] == original['source_sha256'])
    check('compatibility source and execution bytes', digest('source/electratwin_core_numpy_compat.py') == source['compatibility_python_sha256'] == compat['source_sha256'])
    first = (ROOT / 'source/electratwin_core.py').read_bytes()
    second = (ROOT / 'source/electratwin_core_numpy_compat.py').read_bytes()
    check('only one compatibility call replacement', first.count(b'np.trapz(') == 1 and second == first.replace(b'np.trapz(', b'np.trapezoid('))
    check('original failure and compatibility success retained', original['exit_code'] == 1 and compat['exit_code'] == 0 and not original['timed_out'] and not compat['timed_out'])
    check('neither source run connected physical hardware', not original['physical_instrument_connected'] and not compat['physical_instrument_connected'])
    check('audit tied to captured arrays', digest('results/compatibility/captured_model_arrays.json') == audit['source_capture_sha256'])
    check('source conservation failure retained', abs(audit['single_case']['relative_current_mass_flux_mismatch_n1'] - .15035134329700786) < 1e-12)
    c = audit['campaign']
    check('source clipping and fixed grid retained', len(c['rows']) == 12 and c['clipped_fe_rows'] == 8 and set(c['budget_parameter_probe_using_stub_solver'].values()) == {12})
    check('no source EHVI or surrogate falsely credited', not c['ehvi_computed'] and not c['gp_surrogate_fitted'])

    verification = read('results/transport/verification.json')
    baseline = read('results/transport/baseline_solution.json')
    check('transport verification accounting', verification['total_transport_solves'] == 39 and verification['all_finite_volume_solutions_passed_numerical_checks'])
    check('baseline converged and concentrations nonnegative', baseline['summary']['converged'] and baseline['summary']['nonnegative'])
    check('baseline conversion from same result', math.isclose(baseline['summary']['conversion_pct'], verification['baseline']['conversion_pct'], abs_tol=1e-12))
    check('global numerical balances', verification['maximum_material_balance_relative_error'] < 1e-8 and verification['maximum_charge_balance_relative_error'] < 1e-8)
    check('bounded refinement evidence', verification['last_joint_refinement_delta_pp'] < .04 and verification['screening_domain_maximum_abs_conversion_delta_pp'] < .15)
    for name, count in [('grid_convergence.csv', 7), ('parameter_sweep.csv', 12), ('limiting_cases.csv', 5), ('analytic_plug_limit.csv', 4), ('screening_domain_grid_check.csv', 5)]:
        check('transport rows ' + name, len(rows('results/transport/' + name)) == count)

    optimization = read('results/control/optimization_summary.json')
    pool = rows('results/control/candidate_pool.csv')
    campaigns = rows('results/control/sequential_campaigns.csv')
    paired = rows('results/control/paired_budget_comparison.csv')
    check('control script provenance', digest('scripts/metrics_control.py') == optimization['source_sha256'])
    check('transport script provenance', digest('scripts/transport_reviewed.py') == optimization['transport_source_sha256'])
    check('pool provenance', digest('results/control/candidate_pool.csv') == optimization['candidate_pool_sha256'])
    check('81 unique candidates', len(pool) == 81 and len({r['candidate_index'] for r in pool}) == 81)
    check('81 pool evaluations and 240 recorded uses kept separate', optimization['pool_design']['unique_PDE_evaluations'] == 81 and optimization['budget']['total_recorded_evaluation_uses'] == len(campaigns) == 240)
    check('pool numerical balances', all(float(r['material_balance_relative_error']) < 1e-8 and float(r['charge_balance_relative_error']) < 1e-8 for r in pool))
    check('pool evidence role', all(r['evidence_role'] == 'uncalibrated_transport_simulation_with_assumed_product_mass' for r in pool))
    check('fixed FE objective degeneracy retained', optimization['FE_is_fixed_by_single_reaction_charge_identity'] and optimization['STY_FE_model_identity_front_count'] == 1)
    check('no unseen-target or experimental-noise claim', not optimization['method']['unseen_objectives_available_to_GP'] and optimization['method']['experimental_noise_model'] is None)
    cross = read('results/cross_review.json')
    check('400 differential HVI comparisons passed', cross['passed'] and cross['total_candidate_comparisons'] == 400 and cross['maximum_absolute_difference'] < 1e-12)
    for name, sha in cross['source_sha256'].items():
        check('cross-review source ' + name, hashlib.sha256((ROOT.parent / name).read_bytes()).hexdigest() == sha)
    for seed in range(8):
        starts = []
        for method in ['gp_mc_ehvi', 'random']:
            run = [r for r in campaigns if int(r['seed']) == seed and r['method'] == method]
            check(f'budget seed {seed} {method}', len(run) == 15 and [int(r['evaluation']) for r in run] == list(range(1, 16)))
            check(f'no duplicate seed {seed} {method}', len({r['candidate_index'] for r in run}) == 15)
            check(f'nondecreasing hypervolume seed {seed} {method}', all(float(b['hypervolume']) >= float(a['hypervolume']) - 1e-12 for a, b in zip(run, run[1:])))
            starts.append([r['candidate_index'] for r in run[:5]])
        check(f'shared initial data seed {seed}', starts[0] == starts[1])
    comp = optimization['comparison']
    interval = comp['paired_seed_mean_difference_95pct_t_interval']
    check('negative method comparisons retained', len(paired) == 8 and comp['EHVI_wins'] == 5 and comp['random_wins'] == 3 and interval[0] < 0 < interval[1])
    check('paired mean reproduced', abs(sum(float(r['paired_difference']) for r in paired) / 8 - comp['mean_paired_HV_difference']) < 1e-12)
    metrics = read('results/control/baseline_metrics.json')
    check('extra baseline PDE accounted for', metrics['additional_PDE_evaluations_outside_81_point_pool'] == 1)
    check('baseline metrics tied to reviewed current', math.isclose(metrics['transport_summary']['current_A'], baseline['summary']['current_A'], abs_tol=1e-12))
    check('incomplete inventories not formal metrics', metrics['metrics']['full_process_PMI'] is None and metrics['metrics']['atom_economy_pct'] is None)
    scpi = read('results/control/scpi_simulator_trace.json')
    check('simulator has zero physical connections', scpi['simulation_only'] and scpi['physical_connections'] == 0 and not scpi['vendor_compatibility_verified'])
    check('simulator ends disconnected with output off', scpi['final_state'] == dict(connected=False, output_enabled=False, current_A=0.0))
    engineering = read('results/engineering/summary.json')
    check('engineering tied to current transport bytes', engineering['transport_source_sha256'] == digest('results/transport/baseline_solution.json'))
    check('engineering rows', len(rows('results/engineering/hydraulic_scenarios.csv')) == engineering['hydraulic_cases'] == 27 and len(rows('results/engineering/thermal_scenarios.csv')) == engineering['thermal_cases'] == 9)
    check('no CFD thermal or industrial validation credited', not any(engineering[k] for k in ['thermal_field_solved', 'hydraulic_CFD_solved', 'industrial_validation']))

    figures = read('results/figure_manifest.json')['figures']
    qa = read('results/figure_qa.json')
    check('12 unique figure files', len(figures) == len({r['file'] for r in figures}) == 12)
    for r in figures:
        rel = 'reports/figures/' + r['file']
        check('figure hash ' + r['file'], digest(rel) == r['sha256'])
        for name, sha in r['sources'].items():
            check(r['file'] + ' input ' + name, digest(name) == sha)
        if r['format'] == 'png':
            with Image.open(ROOT / rel) as im:
                check('PNG size and resolution ' + r['file'], list(im.size) == r['pixels'] == [2160, 1140] and all(abs(v - 300) < .1 for v in im.info.get('dpi', [0, 0])))
        else:
            tree = ET.parse(ROOT / rel)
            ns = {'s': 'http://www.w3.org/2000/svg'}
            text_nodes = tree.findall('.//s:text', ns)
            text = ' '.join(''.join(t.itertext()) for t in text_nodes)
            check('editable pure vector ' + r['file'], len(text_nodes) > 10 and not tree.findall('.//s:image', ns) and r['editable_text'] and not r['embedded_raster'])
            check('visible evidence label ' + r['file'], any(w in text for w in ['UNVALIDATED', 'NUMERICAL AUDIT', 'SIMULATION ONLY', '未验证', '数值审计', '仅模拟']))
    check('visual QA covers current bytes', qa['reviewed_figure_hashes'] == {r['file']: r['sha256'] for r in figures} and qa['png_files_visually_reviewed'] == 6 and qa['svg_files_rendered_and_visually_reviewed'] == 6)

    contents = {}
    for lang in ['english', 'chinese']:
        path = ROOT / 'reports' / f'electratwin_report_{lang}.md'
        text = path.read_text(encoding='utf-8')
        contents[lang] = text
        check(lang + ' five report sections', re.findall(r'^## (\d+)\.', text, re.M) == ['1', '2', '3', '4', '5'])
        for section, count in [(1, 3), (2, 3), (3, 3), (4, 4), (5, 3)]:
            check(lang + f' full subsections {section}', all(re.search(rf'^### {section}\.{i} ', text, re.M) for i in range(1, count + 1)))
        check(lang + ' numerical anchors', all(v in text for v in ['15.03513433', '58.3872858750', '42.2513750193', '10.1373667653', '0.979543496', '0.959559623', '0.0084974830', '0.0693115316']))
        check(lang + ' complete figure captions', all('**' + ('Figure' if lang == 'english' else '图') + ' ' + str(i) in text for i in range(1, 4)))
        check(lang + ' no unresolved template', not re.search(r'@@[A-Z_]+@@', text))
        check(lang + ' balanced code and display math', text.count('```') % 2 == 0 and text.count('\\[') == text.count('\\]'))
    def table_numbers(text):
        return [re.findall(r'[-−]?\d+(?:\.\d+)?(?:e[-+]?\d+)?', line) for line in text.splitlines() if line.startswith('|')]
    check('English Chinese numeric table parity', table_numbers(contents['english']) == table_numbers(contents['chinese']))
    check('English report has no Chinese prose', not re.search(r'[\u4e00-\u9fff]', contents['english']))
    # Exclude verbatim source specification: its embedded claims are the audit target.
    for path in [ROOT / 'README.md', *sorted((ROOT / 'reports').rglob('*.md'))]:
        text = path.read_text(encoding='utf-8')
        for link in re.findall(r'\]\(([^)]+)\)', text):
            if not link.startswith(('http://', 'https://', '#')):
                target = path.parent / link
                generated_record = target.resolve() == (ROOT / 'results/publication_validation.json').resolve()
                check(path.name + ' link ' + link, generated_record or target.exists())
    result = dict(passed=all(c['passed'] for c in checks), checks=len(checks), failed=[c for c in checks if not c['passed']], scope='Saved numerical provenance, report consistency and publication artifacts; not experimental or industrial validation.', items=checks)
    (ROOT / 'results/publication_validation.json').write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps({k: result[k] for k in ['passed', 'checks', 'failed']}, indent=2))
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
