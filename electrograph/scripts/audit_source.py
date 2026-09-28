"""Inspect captured source outputs and probe untrained-network/ranking behavior.

All neural outputs below are random-initialization diagnostics, never Eox labels.
"""
from collections import Counter
import ast
import csv
import hashlib
import json
from pathlib import Path
import runpy

import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import rdFreeSASA

BASE = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_csv(name, rows):
    with (BASE / 'results' / name).open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    torch.set_num_threads(1)
    source = BASE / 'source/electrograph_kmc_core.py'
    compatibility = BASE / 'source/electrograph_kmc_core_compat.py'
    captured_path = BASE / 'results/compatibility/captured_results.json'
    captured = json.loads(captured_path.read_text(encoding='utf-8'))
    tree = ast.parse(source.read_text(encoding='utf-8'))
    assignments = {target.id: node.value for node in ast.walk(tree) if isinstance(node, ast.Assign)
                   for target in node.targets if isinstance(target, ast.Name)}
    hardcoded = ast.literal_eval(assignments['simulated_tofs'])
    model = runpy.run_path(str(compatibility), run_name='source_audit_only')
    Graph, Net = model['MolecularGraphFeaturizer'], model['ElectroGraphMPNN']
    smiles = captured['substrate_library']
    graphs = [Graph.featurize_molecule(s) for s in smiles]
    symbols = [[a.GetSymbol() for a in Chem.AddHs(Chem.MolFromSmiles(s)).GetAtoms()] for s in smiles]
    rows = []
    for seed in range(12):
        torch.manual_seed(seed)
        net = Net(14, 6, hidden_dim=64, num_steps=3).eval()
        with torch.no_grad():
            for i, graph in enumerate(graphs):
                value, scores, _ = net(graph)
                idx = int(scores.argmax())
                rows.append(dict(seed=seed, molecule_index=i, smiles=smiles[i],
                    raw_untrained_molecular_output=float(value.item()),
                    argmax_atom_index=idx, argmax_element=symbols[i][idx],
                    max_untrained_atom_score=float(scores[idx]),
                    evidence='untrained_random_network_not_oxidation_potential_or_Fukui'))
    save_csv('untrained_initialization_probes.csv', rows)
    torch.manual_seed(42)
    fixed = Net(14, 6, hidden_dim=64, num_steps=3).eval()
    learner = model['DeepKernelActiveLearner'](fixed)
    rank_rows = []
    for seed in range(32):
        np.random.seed(seed)
        rankings = learner.rank_substrates(smiles)
        for rank, row in enumerate(rankings, 1):
            rank_rows.append(dict(uniform_seed=seed, rank=rank, **row,
                                 evidence='same_untrained_weights_random_uniform_exploration'))
    save_csv('random_priority_probes.csv', rank_rows)
    original_atoms = []
    for row in captured['ranked_subs']:
        mol = Chem.AddHs(Chem.MolFromSmiles(row['smiles']))
        original_atoms.append(dict(smiles=row['smiles'], source_atom_index=row['primary_radical_atom_idx'],
            actual_element=mol.GetAtomWithIdx(row['primary_radical_atom_idx']).GetSymbol()))
    target = Chem.AddHs(Chem.MolFromSmiles(captured['target_smiles']))
    radii = list(rdFreeSASA.classifyAtoms(target))
    F, R, T = 96485.33, 8.314, 298.15
    non_electrochemical_ceiling = 1 / (1 / 45 + 1 / 350 + 1 / 80)
    curve = []
    for eta, claimed in zip(np.linspace(.2, .7, 6), hardcoded):
        factor = np.exp(.5 * F * eta / (R * T))
        rates = np.array([45., 120 * factor, 350., 200 * factor, 80.])
        curve.append(dict(eta_V=float(eta), source_hardcoded_TOF_s_1=claimed,
            analytical_independent_cycle_TOF_s_1=float(1 / np.sum(1 / rates)),
            steady_expectation_ceiling_s_1=non_electrochemical_ceiling,
            hardcoded_exceeds_steady_expectation_ceiling=bool(claimed > non_electrochemical_ceiling)))
    save_csv('hardcoded_curve_audit.csv', curve)
    model_spread = []
    for i, smi in enumerate(smiles):
        subset = [r for r in rows if r['molecule_index'] == i]
        values = [r['raw_untrained_molecular_output'] for r in subset]
        model_spread.append(dict(smiles=smi, initialization_count=12, minimum=min(values),
            maximum=max(values), sample_sd=float(np.std(values, ddof=1)),
            argmax_elements=dict(Counter(r['argmax_element'] for r in subset))))
    summary = dict(
        evidence_role='Source-code and numerical behavior audit; no chemical calibration',
        source_sha256=digest(source), compatibility_sha256=digest(compatibility),
        captured_sha256=digest(captured_path), audit_script_sha256=digest(Path(__file__)),
        counts=dict(sensitivity_model_initializations=12, fixed_weight_model_initializations=1,
                    total_model_initializations=13, initialization_forward_passes=72,
                    fixed_weight_priority_runs=32, priority_forward_passes=192,
                    new_kMC_trajectories=0, new_conformer_embeddings=0),
        source_neural_model=dict(trainable_parameters=sum(p.numel() for p in fixed.parameters()),
            training_labels_supplied=False, optimizer_or_backpropagation=False,
            trained_checkpoint_loaded=False, measured_reference_electrode_calibration=False,
            GP_model_fitted=False, posterior_uncertainty_computed=False,
            uniform_exploration_range=[.12, .35], original_selected_atom_elements=original_atoms,
            initialization_sensitivity=model_spread,
            fixed_weight_uniform_seed_top_choices=dict(Counter(r['smiles'] for r in rank_rows if r['rank'] == 1))),
        source_conformers=dict(requested_conformers=8, returned_conformers=captured['conf_results']['generated_conformers_count'],
            source_SASA_A2=captured['conf_results']['boltzmann_averaged_sasa'],
            source_unweighted_Rg_A=captured['conf_results']['boltzmann_averaged_rg'],
            atoms=len(radii), classified_zero_radii=sum(r == 0 for r in radii),
            minimizer_return_codes_recorded=False, postoptimization_deduplication=False,
            explicit_ETKDG_seed=False, buried_volume_computed=False,
            weights_role='MMFF energies of returned candidates, not conformational free energies'),
        source_kinetics=dict(captured=captured['kmc_results'], hardcoded_potential_curve=curve,
            rate_model='Identical independent five-state sites; no neighbor events or periodic operations',
            final_coverage_role='Pre-final-event populations, rounded; not a demonstrated steady-state average',
            time_log_role='Pre-event populations paired with post-event time every 50 steps',
            event_budget_role='Fixed 3500 events; source returns requested step count even after early stop',
            true_steady_expectation_ceiling_s_1=non_electrochemical_ceiling,
            ceiling_scope='Expectation for these assumed independent-site rates; finite stochastic estimates may exceed it'),
        source_reaction=dict(captured=captured['cgr_results'], atom_mapping_supplied=False,
            bond_changes_computed=False, reaction_class_is_hardcoded=True,
            partner='CSc1ccccc1 is thioanisole; its S has no S-H bond',
            heavy_atom_equality_is_not_mass_conservation=True),
        output_sha256={n: digest(BASE / 'results' / n) for n in ['untrained_initialization_probes.csv', 'random_priority_probes.csv', 'hardcoded_curve_audit.csv']},
        primary_sources=[
            dict(url='https://proceedings.mlr.press/v70/gilmer17a.html', role='Supervised message-passing framework, not validation of source random weights'),
            dict(url='https://proceedings.mlr.press/v51/wilson16.html', role='Deep kernel learning uses learned kernels and GP marginal likelihood; source has no GP'),
            dict(url='https://www.rdkit.org/docs/source/rdkit.Chem.rdFreeSASA.html', role='FreeSASA options; installed docstring additionally specifies supplied radii and confIdx'),
        ])
    (BASE / 'results/source_audit.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(source_audit='written', zero_radii=summary['source_conformers']['classified_zero_radii'],
                         top_choice_counts=summary['source_neural_model']['fixed_weight_uniform_seed_top_choices'],
                         steady_TOF_ceiling_s_1=non_electrochemical_ceiling), indent=2))


if __name__ == '__main__':
    main()
