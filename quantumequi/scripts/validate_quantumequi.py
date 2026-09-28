"""Validate saved evidence and publication; do not rerun quantum or neural jobs."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re
import sys
from urllib.parse import unquote
import xml.etree.ElementTree as ET
import numpy as np
from PIL import Image

BASE = Path(__file__).resolve().parents[1]
REPO = BASE.parent
RESULT = BASE/'results/publication_validation.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(name):
    return json.loads((BASE/name).read_text(encoding='utf-8-sig'))


def rows(name):
    with (BASE/name).open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def number(row, key):
    return float(row[key])


def close(a, b, rtol=1e-7, atol=1e-10):
    return bool(np.isclose(float(a), float(b), rtol=rtol, atol=atol))


def main(require_cross_review=False):
    checks = []
    def check(name, condition):
        checks.append({'check': name, 'passed': bool(condition)})
    def guarded(name, function):
        try:
            function()
        except Exception as exc:
            checks.append({'check': name, 'passed': False, 'detail': type(exc).__name__+': '+str(exc)})
    def hashes(mapping, root, label):
        for name, expected in mapping.items():
            path = root/name
            check(label+' '+name, path.is_file() and digest(path) == expected)

    def source():
        record = read('source/source_record.json')
        code = (BASE/'source/quantum_egnn_neb_engine.py').read_bytes()
        embedded = re.findall(rb'```python\r?\n(.*?)\r?\n```', (BASE/'source/specification.md').read_bytes(), re.S)
        check('exact embedded script preserved', len(embedded) == 1 and embedded[0]+b'\n' == code)
        check('source hash preserved', digest(BASE/'source/quantum_egnn_neb_engine.py') == record['original_python_sha256'])
        check('specification hash preserved', digest(BASE/'source/specification.md') == record['specification_sha256'])
        check('no original code repairs', record['repairs'] == [])
        execution = read('results/original/execution.json')
        check('original completed without instruments', execution['exit_code'] == 0 and not execution['physical_instrument_connected'])
        hashes(execution['output_sha256'], BASE/'results/original', 'original output')
        captured = read('results/original/captured_results.json')
        payload = read('results/original/quantum_neb_results.json')
        check('completed source capture matches exported metrics', captured['execution_complete']
            and payload['Extended_Huckel_QM'] == captured['eht_results']
            and payload['CI_NEB_Reaction_Path'] == captured['neb_results']
            and payload['RRHO_Statistical_Thermodynamics'] == captured['rrho_results'])
        check('actual source has eight atoms and 24385 parameters', len(captured['elements']) == 8
            and captured['model_parameter_count'] == 24385 and sorted(captured['elements']) == ['C','C','Cu','H','N','N','N','N'])
        check('source itself does not confirm a saddle', captured['rrho_results']['imaginary_frequency_count'] == 2
              and not captured['rrho_results']['saddle_point_confirmed'])
        check('source nominal negative barrier retained', captured['neb_results']['activation_energy_barrier_kcal_mol'] < 0)

    def potential():
        summary = read('results/potential/summary.json')
        hashes(summary['inputs_sha256'], BASE, 'potential input')
        hashes(summary['outputs_sha256'], BASE/'results/potential', 'potential output')
        transforms = rows('results/potential/transformations.csv')
        fd = rows('results/potential/force_finite_differences.csv')
        baseline = rows('results/potential/baselines.csv')
        forces = rows('results/potential/coordinates_forces.csv')
        c = summary['counts']
        check('potential count ledger independently sums', len(transforms) == c['transformation_probes'] == 48
              and len(fd) == c['full_Cartesian_finite_difference_checks'] == 24 and len(baseline) == 6
              and c['energy_force_model_evaluations'] == len(transforms)+len(fd)*48+len(baseline)+3 == 1209
              and c['training_epochs'] == 0)
        for dtype in ('float32', 'float64'):
            values = [r for r in transforms if r['dtype'] == dtype]
            check(dtype+' symmetry maxima recompute', close(max(number(r,'energy_abs_error') for r in values), summary[dtype+'_max_energy_covariance_error'])
                  and close(max(number(r,'force_max_abs_error') for r in values), summary[dtype+'_max_force_covariance_error']))
        check('double finite differences agree with autograd', max(number(r,'force_max_abs_error') for r in fd
              if r['dtype']=='float64' and number(r,'displacement_A') == 1e-4) < 2e-11)
        check('single precision cancellation is visible', max(number(r,'force_max_abs_error') for r in fd
              if r['dtype']=='float32' and number(r,'displacement_A') == 1e-5) > .01)
        for b in baseline:
            selected = [r for r in forces if r['dtype'] == b['dtype'] and r['geometry'] == b['geometry']]
            coords = np.array([[number(r,k) for k in ('x_A','y_A','z_A')] for r in selected])
            f = np.array([[number(r,k) for k in ('force_x_nominal','force_y_nominal','force_z_nominal')] for r in selected])
            check('saved force and torque arithmetic '+b['dtype']+' '+b['geometry'], len(selected) == 8
                and close(np.max(abs(f)), b['max_force_component']) and close(np.linalg.norm(f.sum(0)), b['net_force_norm'])
                and close(np.linalg.norm(np.cross(coords-coords.mean(0),f).sum(0)), b['torque_norm']))
        check('unused final coordinate parameters and ablation', sum(r['count'] for r in summary['unused_energy_parameters'])
            == summary['unused_energy_parameter_count'] == 1088 and summary['last_coordinate_head_ablation']['energy_difference'] == 0
            and summary['last_coordinate_head_ablation']['force_max_difference'] == 0)

    def path():
        for suffix in ('','/pilot'):
            folder='results/path'+suffix
            summary=read(folder+'/summary.json')
            hashes(summary['source_sha256'],REPO,folder+' source')
            hashes(summary['output_sha256'],BASE/folder,folder+' output')
        ledger=read('results/path/summary.json')['counts']
        audit=read('results/path/read_only_audit.json')
        check('path saved-data audit has no failed checks',audit['all_passed'] and audit['checks_total']==audit['checks_passed']==246
              and all(r['passed'] for r in audit['checks']))
        source = read('results/path/source_audit.json')
        hashes(source['source_sha256'], REPO, 'path source')
        hashes(source['output_sha256'], BASE/'results/path', 'path audit output')
        forward, reverse = source['cases']
        check('exact source replay equals captured result', forward['capture_match'] and forward['source_return'] == read('results/original/captured_results.json')['neb_results'])
        check('source forward and reverse workloads', source['counts']['source_energy_force_calls'] == sum(r['source_energy_force_calls'] for r in source['cases']) == 420
              and source['counts']['final_band_energy_force_calls'] == 14)
        check('source sequential update is exactly reconstructed', all(r['sequential_update_reconstruction_max_error_A'] == 0 for r in source['cases'])
              and source['source_is_inplace_sequential'] and source['source_tangent_uses_already_updated_left_neighbor'])
        check('source nonstationarity and stale energy retained', forward['source_candidate_true_max_atomic_force_nominal'] > .01
              and forward['max_stale_energy_difference_nominal'] > 1e-5 and not source['source_has_force_stop'])
        cases = rows('results/path/analytic_cases.csv')
        bands = rows('results/path/analytic_final_images.csv')
        check('two distinct analytic surfaces tested', {r['surface'] for r in cases} == {'tilted_sine_2d','periodic_curve_3d'})
        check('all reviewed cases satisfy forces and saddle order', all(r['converged'] == 'True'
              and number(r,'max_NEB_force') <= number(r,'force_tolerance') and number(r,'climbing_true_force') <= number(r,'force_tolerance')
              and int(r['negative_hessian_modes']) == 1 for r in cases))
        check('forward reverse barriers use same energy reference', all(close(number(r,'forward_barrier')-number(r,'reverse_barrier'),r['reaction_energy']) for r in cases))
        check('final band row count matches all images', len(bands) == sum(int(r['n_images']) for r in cases))
        reversal=read('results/path/analytic_reversal.json')['runs']
        check('analytic NEB workload arithmetic',len(cases)==ledger['sweep_cases']==36 and len(reversal)==ledger['reversal_cases']==4
              and sum(int(r['iterations']) for r in cases)==ledger['sweep_iterations']==7774
              and sum(int(r['iterations']) for r in reversal)==ledger['reversal_iterations']==748
              and sum(int(r['band_energy_force_point_evaluations']) for r in cases)==ledger['sweep_band_point_evaluations']==95904
              and sum(int(r['band_energy_force_point_evaluations']) for r in reversal)==ledger['reversal_band_point_evaluations']==8272
              and ledger['sweep_endpoint_point_evaluations']+ledger['reversal_endpoint_point_evaluations']==80)

    def electronic():
        independent=read('results/electronic/saved_energy_identity_check.json')
        hashes(independent['inputs_sha256'],BASE/'results/electronic','independent electronic input')
        check('independent electronic energy and log audit',independent['passed'] and sum(r['jobs'] for r in independent['records'])==58
              and all(r['passed'] for r in independent['records'])
              and digest(BASE/'results/electronic/verify_saved_energy.py')==independent['script_sha256'])
        for suffix in ('','/pilot'):
            folder='results/electronic'+suffix
            summary=read(folder+'/summary.json')
            for name,value in summary['inputs_sha256'].items():
                actual=(BASE/folder/'executed_script.py.txt') if suffix and name=='scripts/electronic_reviewed.py' else BASE/name
                check(folder+' input '+name,digest(actual)==value)
            hashes(summary['outputs_sha256'],BASE/folder,folder+' output')
            snapshot=summary['executed_source_snapshot']
            check(folder+' historical executed source',digest(BASE/folder/snapshot['file'])==snapshot['sha256'])
            jobs=rows(folder+'/quantum_job_accounting.csv')
            matrices=read(folder+'/ao_matrices.json')
            count=summary['sections']['quantum_reference']['counts']['total_SCF_jobs']
            check(folder+' SCF ledger and AO matrix count',len(jobs)==len(matrices)==count==(3 if suffix else 55))
            for matrix in matrices:
                s,f,c,e=[np.array(matrix[k]) for k in ('S','Fock_Hartree','C_alpha','epsilon_alpha_Hartree')]
                residual=np.max(abs(f@c-(s@c)*e))
                orthogonality=np.max(abs(c.T@s@c-np.eye(len(e))))
                check(folder+' AO residual '+matrix['job_id'],residual<1e-11 and orthogonality<1e-11)
            for job in jobs:
                check(folder+' energy decomposition '+job['job_id'],close(number(job,'energy_Hartree'),number(job,'electronic_energy_Hartree')+number(job,'nuclear_repulsion_Hartree')))
            for row in rows(folder+'/h2_gradient_finite_differences.csv'):
                fd=(number(row,'plus_energy_Hartree')-number(row,'minus_energy_Hartree'))/(2*number(row,'step_A'))
                check(folder+' central difference '+row['R_A']+' '+row['step_A'],close(fd,row['finite_difference_dE_dR_Hartree_A'])
                      and close(abs(fd-number(row,'analytic_dE_dR_Hartree_A')),row['absolute_error_Hartree_A']))
        audit=read('results/electronic/source_audit.json')
        matrices=read('results/electronic/source_matrices.json')
        for case in audit['cases']:
            values=matrices[case['case']]
            s=np.array(values['S_dimensionless']);h=np.array(values['H_eV'])
            c=np.array(values['canonical_coefficients']);e=np.array(values['canonical_energies_eV'])
            check('source overlap rank and capacity '+case['case'],np.linalg.matrix_rank(s)==case['positive_overlap_rank']==8
                  and case['discarded_null_directions']==23 and 2*case['positive_overlap_rank']<case['valence_electrons']
                  and bool(case['strict_solver_rejection']) and bool(case['closed_shell_capacity_rejection']))
            check('canonical original residual stays nonzero '+case['case'],np.max(abs(h@c-(s@c)*e))>4
                  and np.max(abs(c.T@s@c-np.eye(8)))<1e-11)
        reference=rows('results/electronic/h2_reference.csv')
        groups=defaultdict(set)
        for row in reference:
            groups[row['split']].add(row['point_id'])
        check('42 disjoint H2 distances and 17 8 8 9 split',len(reference)==42 and len({r['R_A'] for r in reference})==42
              and {k:len(v) for k,v in groups.items()}=={'train':17,'validation':8,'test':8,'ood_stretch':9}
              and sum(map(len,groups.values()))==len(set.union(*groups.values())))
        predictions=rows('results/electronic/surrogate_predictions.csv')
        metrics=rows('results/electronic/surrogate_metrics.csv')
        check('five models on all 42 H2 points',len(predictions)==210 and len(metrics)==20)
        for metric in metrics:
            points=[r for r in predictions if r['model']==metric['model'] and r['split']==metric['split']]
            de=np.array([number(r,'predicted_energy_Hartree')-number(r,'reference_energy_Hartree') for r in points])
            dg=np.array([number(r,'predicted_dE_dR_Hartree_A')-number(r,'reference_dE_dR_Hartree_A') for r in points])
            check('independent surrogate metrics '+metric['model']+' '+metric['split'],len(points)==int(metric['points'])
                  and close(np.sqrt(np.mean(de**2)),metric['energy_RMSE_Hartree'])
                  and close(np.sqrt(np.mean(dg**2)),metric['gradient_RMSE_Hartree_A'])
                  and close(np.mean(abs(de)),metric['energy_MAE_Hartree']) and close(np.mean(abs(dg)),metric['gradient_MAE_Hartree_A']))

    def thermochemistry():
        for suffix in ('','/pilot'):
            folder='results/thermochemistry'+suffix
            summary=read(folder+'/summary.json')
            code=BASE/folder/'executed_code.py.txt' if suffix else BASE/'scripts/thermochemistry_reviewed.py'
            check(folder+' executed code hash',digest(code)==summary['code_sha256'])
            hashes(summary['output_sha256'],BASE/folder,folder+' output')
        summary=read('results/thermochemistry/summary.json')
        constants=summary['constants']
        factor=np.sqrt(4184/constants['Avogadro_mol_inverse']/1e-20/constants['atomic_mass_kg_CODATA2022'])/(2*np.pi*constants['speed_light_cm_s'])
        check('independent frequency conversion from SI constants',close(factor,constants['kcal_mol_A2_amu_to_cm_factor'],rtol=1e-12))
        audit=read('results/thermochemistry/source_rrho_audit.json')
        a=audit['autograd_float64_analysis']
        check('source nonstationary candidate is rejected',not a['stationary'] and a['stationary_point_classification']=='nonstationary'
              and a['maximum_atomic_gradient']>.01 and not audit['valid_source_RRHO_released'])
        expected=np.sign(a['projected_eigenvalues'])*np.sqrt(abs(np.array(a['projected_eigenvalues'])))*factor
        check('source projected frequencies recompute without blind deletion',np.allclose(expected,a['projected_signed_frequencies_cm'],rtol=1e-12)
              and len(expected)==18 and sum(expected<-.01)==1 and a['rigid_rank']==6 and audit['source_deleted_negative_modes']==2)
        check('source empirical corrections retained as defects',audit['source_empirical_translation_rotation_entropy_cal_mol_K']==78.5
              and audit['source_missing_thermal_enthalpy'] and audit['source_plot_hardcoded_TS_product_offsets_kcal_mol']==[3.2,-1.8])
        table=rows('results/thermochemistry/ideal_gas_temperature_pressure.csv')
        check('RRHO parameter and Hessian workload ledger',len(table)==summary['counts_exclude_pilot_and_tests']['rrho_parameter_cases']==80
              and len(rows('results/thermochemistry/analytic_hessian_convergence.csv'))==summary['counts_exclude_pilot_and_tests']['analytic_fd_hessians']==40
              and len(rows('results/thermochemistry/source_hessian_convergence.csv'))==summary['source_audit_counts']['fd_hessians']==8
              and len(rows('results/thermochemistry/low_frequency_sensitivity.csv'))==summary['counts_exclude_pilot_and_tests']['positive_mode_sensitivity_cases']==72)
        grouped=defaultdict(list)
        for row in table:
            t=number(row,'temperature_K')
            expected=number(row,'zpve_kcal_mol')+number(row,'thermal_enthalpy_excluding_ZPVE_kcal_mol')-t*number(row,'total_entropy_cal_mol_K')/1000
            check('G correction arithmetic '+row['case']+' '+row['temperature_K']+' '+row['pressure_Pa'],close(expected,row['G_correction_kcal_mol']))
            grouped[(row['case'],t)].append(row)
        for key, group in grouped.items():
            ref=group[0]
            check('ideal-gas pressure law '+str(key),all(close(number(r,'G_correction_kcal_mol')-number(ref,'G_correction_kcal_mol'),
                constants['gas_constant_cal_mol_K']*key[1]/1000*np.log(number(r,'pressure_Pa')/number(ref,'pressure_Pa'))) for r in group))

    def figures():
        manifest = read('results/figure_manifest.json')['figures']
        check('eight PNG and eight SVG figures', len(manifest)==16 and sum(r['format']=='png' for r in manifest)==8 and sum(r['format']=='svg' for r in manifest)==8)
        for row in manifest:
            path = BASE/'reports/figures'/row['file']
            check('figure hash '+row['file'], path.is_file() and digest(path)==row['sha256'])
            hashes(row['sources'], BASE, row['file']+' table')
            if row['format']=='png':
                with Image.open(path) as im:
                    check('PNG dimensions and DPI '+row['file'], list(im.size)==row['pixels']==[2400,1380]
                          and all(abs(d-300)<.1 for d in im.info.get('dpi',(0,0))))
            else:
                tree=ET.parse(path)
                ns={'s':'http://www.w3.org/2000/svg'}
                check('editable SVG '+row['file'], len(tree.findall('.//s:text',ns))>10 and not tree.findall('.//s:image',ns))
        qa=read('results/figure_qa.json')
        check('visual review covers current figures', qa['reviewed_figure_hashes']=={r['file']:r['sha256'] for r in manifest}
              and qa['png_files_visually_reviewed']==8 and qa['svg_files_rendered_and_visually_reviewed']==8
              and qa['generator_sha256']==digest(BASE/'scripts/plot_reviewed.py'))

    def publication():
        editions=[]
        for language in ('english','chinese'):
            path=BASE/'reports'/('quantumequi_report_'+language+'.md')
            text=path.read_text(encoding='utf-8')
            editions.append(text)
            check(language+' six report sections', re.findall(r'^## (\d+)\.',text,re.M)==list('123456'))
            check(language+' balanced delimiters', text.count('```')%2==0 and text.count('\\[')==text.count('\\]') and text.count('\\(')==text.count('\\)'))
            check(language+' four PNG and SVG references', len(re.findall(r'\]\(figures/[^)]+\.png\)',text))==4 and len(re.findall(r'\]\(figures/[^)]+\.svg\)',text))==4)
        nums=lambda text:[re.findall(r'(?<![A-Za-z0-9.])[-−]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?',line) for line in text.splitlines() if line.startswith('|')]
        check('bilingual numerical table parity', nums(editions[0])==nums(editions[1]))
        check('English report has no Chinese prose', not re.search(r'[\u4e00-\u9fff]',editions[0]))
        check('module README exists', (BASE/'README.md').is_file())
        for path in list(BASE.glob('*.md'))+list((BASE/'reports').rglob('*.md')):
            for link in re.findall(r'\]\(([^)]+)\)',path.read_text(encoding='utf-8')):
                if link.startswith(('https://','http://','#','data:','mailto:')):
                    continue
                target=path.parent/unquote(link.split('#',1)[0].strip('<>'))
                check(path.relative_to(BASE).as_posix()+' link '+link,target.exists() or target.resolve()==RESULT.resolve())
        path=BASE/'results/cross_review.json'
        if require_cross_review or path.exists():
            review=read('results/cross_review.json')
            check('cross-review has no unresolved blockers', bool(review['reviews']) and all(r['blocking_findings']==[] for r in review['reviews']))
            for group,mapping in review['reviewed_artifact_sha256'].items():
                hashes(mapping,BASE,'cross-review '+group)

    for name, fn in [('source',source),('potential',potential),('path',path),('electronic',electronic),('thermochemistry',thermochemistry),('figures',figures),('publication',publication)]:
        guarded(name,fn)
    result={'passed':all(c['passed'] for c in checks),'checks':len(checks),'failed':[c for c in checks if not c['passed']],
            'scope':'Saved evidence, hashes, independent arithmetic and publication integrity; no chemical validation or study reruns',
            'cross_review_required':require_cross_review,'validator_sha256':digest(Path(__file__)),'items':checks}
    RESULT.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:result[k] for k in ('passed','checks','failed')},indent=2,ensure_ascii=False))
    return 0 if result['passed'] else 1


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--require-cross-review',action='store_true')
    sys.exit(main(parser.parse_args().require_cross_review))
