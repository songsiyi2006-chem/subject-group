"""Independent saved-evidence checks; never executes a quantum engine."""
import json
import math
import re
from .workflow import BASE, digest, plan, write_json


def main():
    folder = BASE/'results'
    load = lambda name: json.loads((folder/name).read_text('utf-8'))
    checks = []
    def check(name, condition):
        checks.append(dict(name=name, passed=bool(condition)))
    saved = load('prepared_plan.json')
    check('frozen four-job preparation matches source', saved == plan() and len(saved['jobs']) == 4)
    result = load('local_pilot_result.json')
    provenance = load('local_pilot_provenance.json')
    check('pilot task matches frozen first task', result['job'] == saved['jobs'][0])
    check('one local completed attempt, zero remote jobs', provenance['local_attempts'] == provenance['completed_single_points'] == 1 and provenance['remote_submissions'] == 0)
    check('native result bytes retained', digest(folder/'local_pilot_result.json') == provenance['native_result_sha256'])
    check('original log identity retained', result['engine_output_sha256'] == provenance['original_engine_output_sha256'])
    check('portable log hash', digest(folder/'local_pilot_engine_portable.out') == provenance['portable_engine_output_sha256'])
    log = (folder/'local_pilot_engine_portable.out').read_text('utf-8')
    energies = re.findall(r'@DF-RKS Final Energy:\s*([-+0-9.Ee]+)', log)
    check('independent native energy extraction', len(energies) == 1 and math.isclose(float(energies[0]), result['energy_hartree'], rel_tol=0, abs_tol=1e-10))
    check('SCF completed', result['status'] == 'completed' and result['scf_converged'] is True)
    check('limitations preserved', result['wavefunction_stability'] == 'not_checked' and result['optimization'] == result['frequencies'] == 'not_performed')
    for name, expected in result['code_sha256'].items():
        check('executed code matches shipped code '+name, digest(BASE/name) == expected)
    collected = load('local_pilot_collection.json')
    check('missing partners not invented', len(collected['pairs']) == 2 and all(p['status'] == 'incomplete' and p['energy_difference_hartree'] is None for p in collected['pairs']))
    check('zero accepted chemical predictions', collected['accepted_chemical_predictions'] == 0)
    check('portable log omits local paths', not re.search(r'[A-Z]:[\\/](?:Users|Codex)[\\/]', log))
    audit = dict(passed=all(c['passed'] for c in checks), checks=checks,
                 scope='Preparation and saved single-point record consistency, not independent physical validation.')
    write_json(folder/'validation.json', audit)
    print(json.dumps({'passed': audit['passed'], 'checks': len(checks), 'failed': [c for c in checks if not c['passed']]}))
    return 0 if audit['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
