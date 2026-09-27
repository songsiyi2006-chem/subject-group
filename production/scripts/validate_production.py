"""Validate publication/data integrity without asserting experimental validity."""
import csv
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    checks = []

    def check(name, condition):
        checks.append({'check': name, 'passed': bool(condition)})

    def read_json(relative):
        return json.loads((ROOT / relative).read_text(encoding='utf-8-sig'))

    def rows(relative):
        with (ROOT / relative).open(encoding='utf-8', newline='') as f:
            return list(csv.DictReader(f))

    record = read_json('source/source_record.json')
    extracted = ROOT / 'source/run_production_pipeline_original.py'
    check('original extracted script SHA256', hashlib.sha256(extracted.read_bytes()).hexdigest() == record['extracted_script_sha256'])
    for variant, expected in [('original_attempt', 1), ('repaired', 0)]:
        execution = read_json(f'results/{variant}/execution.json')
        check(variant + ' exit status', execution['returncode'] == expected)
        check(variant + ' result presence', execution['result_written'] == (expected == 0))
        script = extracted if variant == 'original_attempt' else ROOT / 'run_production_pipeline.py'
        check(variant + ' executed script hash', hashlib.sha256(script.read_bytes()).hexdigest() == execution['sha256'])
    check('original full JSON absent', not (ROOT / 'results/original_attempt/production_benchmark_results.json').exists())
    check('original API error preserved', 'calcSASA' in (ROOT / 'results/original_attempt/stderr.log').read_text())

    result = read_json('results/repaired/production_benchmark_results.json')
    check('four task results', len(result) == 4 and all(any(k.startswith(f'Task_{letter}_') for k in result) for letter in 'ABCD'))
    for name, size in [('repaired/task_a_candidates.csv', 288), ('repaired/task_a_observations.csv', 14),
                       ('repaired/task_b_all_atoms.csv', 30), ('repaired/task_c_all_sites.csv', 16),
                       ('audit/task_a_paired_runs.csv', 120), ('audit/task_a_observations.csv', 1680), ('audit/task_b_conformers.csv', 32),
                       ('audit/task_d_eta_sweep.csv', 48), ('audit/task_d_parameter_sensitivity.csv', 20)]:
        check('row count ' + name, len(rows('results/' + name)) == size)
    summary = read_json('results/audit/audit_summary.json')
    check('target xTB completed', summary['target_xTB']['status'] == 'completed' and summary['target_xTB']['jobs'] == 3)
    for state in ['neutral', 'cation', 'anion']:
        check('xTB record ' + state, read_json(f'results/target_xtb/{state}/run.json')['completed'])

    contents = {}
    for language in ['english', 'chinese', 'bilingual']:
        path = ROOT / f'reports/production_technical_report_{language}.md'
        txt = path.read_text(encoding='utf-8')
        contents[language] = txt
        check(language + ' exactly four main sections', re.findall(r'^## (\d+)\.', txt, re.M) == ['1', '2', '3', '4'])
        for section, count in [(1, 3), (2, 5), (3, 4), (4, 3)]:
            check(language + f' section {section} coverage', all(re.search(rf'^### {section}\.{i} ', txt, re.M) for i in range(1, count + 1)))
        check(language + ' numerical anchors', all(s in txt for s in ['92.03', '91.3865', '92.4865', '0.7483', '0.5172', '17.09', '77.55', '1713.41', '10.5141', '6.0021']))
        check(language + ' evidence labels', all('[' + k + ']' in txt for k in ['L', 'S', 'M', 'Q', 'A', 'P', 'U']))
        check(language + ' source identity', record['extracted_script_sha256'] in txt)
        check(language + ' no unresolved placeholders', not re.search(r'@@[A-Z_]+@@', txt))
        check(language + ' balanced code fences', txt.count('```') % 2 == 0)
        check(language + ' balanced display math', txt.count('$$') % 2 == 0)
        for link in re.findall(r'\]\(([^)]+)\)', txt):
            if not link.startswith(('http://', 'https://', '#')):
                check(language + ' link ' + link, (path.parent / link).is_file())

    def table_decimals(txt):
        return [re.findall(r'[-−]?\d+\.\d+', line) for line in txt.splitlines() if line.startswith('|')]

    check('separate editions table decimal parity', table_decimals(contents['english']) == table_decimals(contents['chinese']))
    check('English edition omits Chinese prose', not re.search(r'[\u4e00-\u9fff]', contents['english']))
    for language in ['english', 'chinese']:
        # Each complete subsection body must survive the section-aligned merge.
        parts = re.split(r'^### [^\n]+\n', contents[language], flags=re.M)[1:]
        bodies = [re.split(r'^## ', p, flags=re.M)[0].strip() for p in parts]
        check(language + ' all subsection bodies preserved in combined edition', all(p in contents['bilingual'] for p in bodies))
    for p in sorted((ROOT / 'results').rglob('*.json')):
        if p.name == 'publication_validation.json':
            continue
        try:
            json.loads(p.read_text(encoding='utf-8-sig'))
            ok = True
        except (ValueError, OSError):
            ok = False
        check('JSON parses ' + p.relative_to(ROOT).as_posix(), ok)
    output = {'passed': all(c['passed'] for c in checks), 'checks': len(checks),
              'scope': 'Publication and saved-record integrity; scientific invariants are separately checked by tests/test_production.py. No experimental validity is asserted.',
              'failed': [c for c in checks if not c['passed']], 'items': checks}
    # Deterministic content so successful validation does not invalidate release hashes.
    (ROOT / 'results/publication_validation.json').write_text(json.dumps(output, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps({k: output[k] for k in ['passed', 'checks', 'failed']}, indent=2))
    return 0 if output['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
