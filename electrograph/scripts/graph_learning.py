"""Bounded FreeSolv learning and frozen-pool GP-UCB benchmark.

Targets are experimental hydration free energies, kcal/mol; never Eox or Fukui.
The molecular encoder is trained only on the training partition. GP posterior
uncertainties are conditional model uncertainties, not measurement error bars.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import time
from dataclasses import dataclass

import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
import scipy
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import t as student_t
import sklearn
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
import torch
from torch import nn

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'electrograph/data/learning'
OUT = ROOT / 'electrograph/results/learning'
ELEMENTS = ('C', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'I', 'P', 'B', 'Si', 'H')
CHIRAL = (None, 'R', 'S')  # Absolute CIP, not atom-order-dependent CW/CCW tags.
BONDS = (Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE,
         Chem.BondType.AROMATIC)
STEREO = (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOE,
          Chem.BondStereo.STEREOZ)
NODE_DIM, EDGE_DIM = len(ELEMENTS) + 1 + 7 + len(CHIRAL) + 1, 11
DESCRIPTOR_NAMES = ['MolWt', 'MolLogP', 'TPSA', 'HBD', 'HBA', 'RotatableBonds',
                    'HeavyAtoms', 'Heteroatoms', 'Rings', 'FractionCSP3']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False,
                                     allow_nan=False) + '\n', encoding='utf-8')


def write_csv(path, rows):
    if not rows:
        raise ValueError('Refusing an empty result table')
    with Path(path).open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def one_hot_unknown(value, choices):
    return [float(value == x) for x in choices] + [float(value not in choices)]


@dataclass
class Graph:
    x: torch.Tensor
    edge_index: torch.Tensor
    edge_attr: torch.Tensor
    smiles: str = ''


def molecular_graph(smiles: str) -> Graph:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError('A nonempty valid molecule is required')
    # Heavy-atom graph with implicit H counts; no invented electronegativities.
    atoms = [one_hot_unknown(a.GetSymbol(), ELEMENTS) +
             [a.GetFormalCharge() / 3, a.GetDegree() / 4,
              a.GetTotalNumHs() / 4, float(a.GetIsAromatic()),
              float(a.IsInRing()), a.GetMass() / 200, a.GetIsotope() / 250] +
             one_hot_unknown(a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else None, CHIRAL)
             for a in mol.GetAtoms()]
    edges, attributes = [], []
    for b in mol.GetBonds():
        feat = one_hot_unknown(b.GetBondType(), BONDS) + [float(b.GetIsConjugated()),
                float(b.IsInRing())] + one_hot_unknown(b.GetStereo(), STEREO)
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        edges.extend([[i, j], [j, i]])
        attributes.extend([feat, feat])
    return Graph(torch.tensor(atoms, dtype=torch.float32),
                 torch.tensor(edges, dtype=torch.long).reshape(-1, 2).T.contiguous(),
                 torch.tensor(attributes, dtype=torch.float32).reshape(-1, EDGE_DIM),
                 Chem.MolToSmiles(mol, isomericSmiles=True))


def batch_graphs(graphs):
    if not graphs:
        raise ValueError('Empty graph batch')
    xs, ei, ea, memberships, offset = [], [], [], [], 0
    for i, graph in enumerate(graphs):
        xs.append(graph.x)
        ei.append(graph.edge_index + offset)
        ea.append(graph.edge_attr)
        memberships.append(torch.full((len(graph.x),), i, dtype=torch.long))
        offset += len(graph.x)
    return (torch.cat(xs), torch.cat(ei, 1), torch.cat(ea),
            torch.cat(memberships), len(graphs))


class EdgeConditionedLayer(nn.Module):
    def __init__(self, hidden_dim=24):
        super().__init__()
        self.hidden_dim = hidden_dim
        self.edge_net = nn.Sequential(nn.Linear(EDGE_DIM, 16), nn.SiLU(),
                                      nn.Linear(16, hidden_dim * hidden_dim))
        self.update = nn.GRUCell(hidden_dim, hidden_dim)

    def forward(self, x, edge_index, edge_attr):
        aggregated = torch.zeros_like(x)
        if edge_index.shape[1]:
            weights = self.edge_net(edge_attr).reshape(-1, self.hidden_dim, self.hidden_dim)
            messages = torch.bmm(weights, x[edge_index[0]].unsqueeze(-1)).squeeze(-1)
            aggregated.index_add_(0, edge_index[1], messages / math.sqrt(self.hidden_dim))
        # Isolated atoms receive a defined zero-message GRU update.
        return self.update(aggregated, x)


class HydrationMPNN(nn.Module):
    def __init__(self, hidden_dim=24, steps=2):
        super().__init__()
        self.embedding = nn.Sequential(nn.Linear(NODE_DIM, hidden_dim), nn.SiLU())
        self.layers = nn.ModuleList([EdgeConditionedLayer(hidden_dim) for _ in range(steps)])
        self.readout = nn.Sequential(nn.Linear(2 * hidden_dim, 32), nn.SiLU(), nn.Linear(32, 1))

    def forward(self, batch):
        x, edge_index, edge_attr, membership, count = batch
        h = self.embedding(x)
        for layer in self.layers:
            h = layer(h, edge_index, edge_attr)
        total = torch.zeros((count, h.shape[1]), dtype=h.dtype, device=h.device)
        total.index_add_(0, membership, h)
        sizes = torch.bincount(membership, minlength=count).to(h.dtype).unsqueeze(1)
        latent = torch.cat([total, total / sizes], dim=1)
        return self.readout(latent).squeeze(1), latent


def descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError('Invalid molecule')
    return np.array([Descriptors.MolWt(mol), Descriptors.MolLogP(mol),
                     rdMolDescriptors.CalcTPSA(mol), rdMolDescriptors.CalcNumHBD(mol),
                     rdMolDescriptors.CalcNumHBA(mol), rdMolDescriptors.CalcNumRotatableBonds(mol),
                     mol.GetNumHeavyAtoms(), rdMolDescriptors.CalcNumHeteroatoms(mol),
                     rdMolDescriptors.CalcNumRings(mol), rdMolDescriptors.CalcFractionCSP3(mol)])


def structure_key(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or not mol.GetNumAtoms():
        raise ValueError('Invalid molecule')
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def group_key(smiles):
    mol = Chem.MolFromSmiles(smiles)
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    return ('ring:' + scaffold) if scaffold else ('acyclic:' + structure_key(smiles))


def assign_splits(records):
    groups = {}
    for row in records:
        groups.setdefault(row['group_key'], []).append(row)
    # Largest groups first, hash breaks size ties; no outcomes used.
    names = ['train', 'validation', 'test']
    targets = np.array([.60, .20, .20]) * len(records)
    sizes = np.zeros(3)
    for key in sorted(groups, key=lambda k: (-len(groups[k]), hashlib.sha256(k.encode()).hexdigest())):
        choice = int(np.argmax(targets - sizes))
        for row in groups[key]:
            row['split'] = names[choice]
        sizes[choice] += len(groups[key])
    if min(sizes) < 4:
        raise ValueError('Insufficient distinct groups for split')
    return records


def prepare_dataset(source_csv, count=256):
    if count < 32:
        raise ValueError('At least 32 structures required')
    official = {}
    for line in (DATA / 'official_database.txt').read_text(encoding='utf-8').splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        parts = [p.strip() for p in line.split(';')]
        key = structure_key(parts[1])
        if key in official:
            raise ValueError('Official duplicate structure requires manual reconciliation')
        official[key] = {'official_id': parts[0], 'experimental_uncertainty_kcal_mol': float(parts[4]),
                         'experimental_reference': parts[7], 'official_expt': float(parts[3])}
    records, seen, discrepancies = [], set(), []
    with Path(source_csv).open(encoding='utf-8-sig') as handle:
        local_rows = list(csv.DictReader(handle))
    for i, row in enumerate(local_rows):
        key = structure_key(row['smiles'])
        if key in seen:
            raise ValueError('Local canonical duplicate requires manual reconciliation')
        seen.add(key)
        upstream = official.get(key)
        if upstream is None or abs(float(row['expt']) - upstream['official_expt']) > 1e-8:
            discrepancies.append({'local_row': i, 'smiles': key, 'local_expt': float(row['expt']),
                                  'official_expt': None if upstream is None else upstream['official_expt']})
            continue
        records.append({'molecule_id': upstream['official_id'], 'local_row': i, 'smiles': key,
                        'expt_kcal_mol': float(row['expt']),
                        'experimental_uncertainty_kcal_mol': upstream['experimental_uncertainty_kcal_mol'],
                        'experimental_reference': upstream['experimental_reference'],
                        'group_key': group_key(key), 'selection_hash': hashlib.sha256(key.encode()).hexdigest()})
    records = assign_splits(sorted(records, key=lambda row: row['selection_hash'])[:count])
    if len(records) != count:
        raise ValueError('Not enough matched unique structures')
    write_csv(DATA / 'selected_freesolv.csv', records)
    provenance = {
        'dataset': 'FreeSolv experimental hydration free energies, local SAMPL snapshot',
        'target_unit': 'kcal/mol', 'source_rows': len(local_rows), 'official_rows': len(official),
        'source_filename': Path(source_csv).name, 'source_sha256': sha(source_csv),
        'official_database_sha256': sha(DATA / 'official_database.txt'),
        'license_sha256': sha(DATA / 'FreeSolv_LICENSE.txt'),
        'selected_sha256': sha(DATA / 'selected_freesolv.csv'), 'selected_count': len(records),
        'selected_label_agreement_tolerance': 1e-8, 'excluded_discrepancies': discrepancies,
        'selection': 'First N canonical-structure SHA256 values, label blind, matched to official values',
        'split': '60/20/20 target, whole Murcko ring scaffolds, largest-group greedy deficit allocation',
        'acyclic_policy': 'Full canonical structure as group; does NOT establish acyclic scaffold novelty',
        'license': 'CC-BY-4.0, with upstream caveat that original third-party restrictions may apply',
        'attribution': 'David L. Mobley and J. Peter Guthrie, FreeSolv, JCAMD (2014), doi:10.1007/s10822-014-9747-x',
        'retrieved_date': '2026-09-28',
        'sources': ['https://github.com/MobleyLab/FreeSolv',
                    'https://github.com/MobleyLab/FreeSolv/blob/master/LICENSE',
                    'https://github.com/MobleyLab/FreeSolv/blob/master/database.txt'],
        'modifications': 'Canonicalized SMILES, selected subset, added split/group/hash columns; experimental values unchanged'}
    dump(DATA / 'provenance.json', provenance)
    return records


def load_records():
    with (DATA / 'selected_freesolv.csv').open(encoding='utf-8') as handle:
        records = list(csv.DictReader(handle))
    for row in records:
        row['expt_kcal_mol'] = float(row['expt_kcal_mol'])
    return records


def regression_metrics(truth, pred):
    residual = np.asarray(pred) - np.asarray(truth)
    return {'RMSE_kcal_mol': float(np.sqrt(np.mean(residual ** 2))),
            'MAE_kcal_mol': float(np.mean(abs(residual)))}


def fit_mpnn(graphs, y, train_ids, validation_ids, seed, epochs, out, shuffled=False):
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    model = HydrationMPNN()
    mean, scale = float(np.mean(y[train_ids])), float(np.std(y[train_ids]))
    if scale <= 0:
        raise ValueError('Constant training targets')
    yfit = y.copy()
    if shuffled:
        yfit[train_ids] = rng.permutation(yfit[train_ids])
    target = torch.tensor((yfit - mean) / scale, dtype=torch.float32)
    optimizer = torch.optim.Adam(model.parameters(), lr=.003, weight_decay=1e-4)
    validation_batch = batch_graphs([graphs[i] for i in validation_ids])
    losses, best_loss, best_state, best_epoch = [], float('inf'), None, None
    for epoch in range(1, epochs + 1):
        model.train()
        order = rng.permutation(train_ids)
        loss_sum = 0.0
        for start in range(0, len(order), 32):
            ids = order[start:start + 32]
            optimizer.zero_grad()
            pred, _ = model(batch_graphs([graphs[i] for i in ids]))
            loss = torch.mean((pred - target[ids]) ** 2)
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            loss_sum += float(loss.detach()) * len(ids)
        model.eval()
        with torch.no_grad():
            pred, _ = model(validation_batch)
            validation_loss = float(torch.mean((pred - target[validation_ids]) ** 2))
        losses.append({'seed': seed, 'model': 'shuffled_train_mpnn' if shuffled else 'mpnn',
                       'epoch': epoch, 'training_standardized_MSE': loss_sum / len(order),
                       'validation_standardized_MSE': validation_loss})
        if validation_loss < best_loss:
            best_loss, best_epoch = validation_loss, epoch
            best_state = copy.deepcopy(model.state_dict())
    model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        pred, latent = model(batch_graphs(graphs))
    name = ('shuffled_train_mpnn' if shuffled else 'mpnn') + f'_seed{seed}'
    torch.save({'state_dict': best_state, 'target_mean': mean, 'target_scale': scale,
                'hidden_dim': 24, 'steps': 2, 'node_dim': NODE_DIM, 'edge_dim': EDGE_DIM,
                'best_epoch': best_epoch, 'seed': seed, 'target': 'hydration_kcal_mol'},
               out / (name + '.pt'))
    return pred.numpy() * scale + mean, latent.numpy(), losses, {
        'model': name, 'seed': seed, 'best_epoch': best_epoch, 'epochs_executed': epochs,
        'validation_standardized_MSE': best_loss, 'target_mean': mean, 'target_scale': scale,
        'parameter_count': sum(p.numel() for p in model.parameters()),
        'training_count': len(train_ids), 'validation_count': len(validation_ids)}


def matern52(x, z):
    # Dimension-normalized Euclidean distance; fixed length scale one.
    distance = np.sqrt(np.maximum(np.sum((x[:, None, :] - z[None, :, :]) ** 2, axis=2), 0) / x.shape[1])
    root5 = math.sqrt(5) * distance
    return (1 + root5 + 5 * distance ** 2 / 3) * np.exp(-root5)


def gp_posterior(features, observed, values, noise_variance=1e-5):
    x, observed, values = np.asarray(features, float), np.asarray(observed, int), np.asarray(values, float)
    if (x.ndim != 2 or not len(observed) or len(observed) != len(values)
            or len(set(observed.tolist())) != len(observed) or noise_variance <= 0
            or not np.all(np.isfinite(x)) or not np.all(np.isfinite(values))
            or min(observed) < 0 or max(observed) >= len(x)):
        raise ValueError('Invalid GP conditioning data')
    ymean = float(values.mean())
    yscale = float(values.std())
    if yscale < 1e-12:
        yscale = 1.0
    K = matern52(x[observed], x[observed]) + noise_variance * np.eye(len(observed))
    cross = matern52(x, x[observed])
    chol = cho_factor(K, lower=True)
    mean = ymean + yscale * (cross @ cho_solve(chol, (values - ymean) / yscale))
    variance = 1 - np.sum(cross * cho_solve(chol, cross.T).T, axis=1)
    if min(variance) < -1e-8:
        raise ArithmeticError('Negative posterior variance exceeds roundoff')
    return mean, yscale * np.sqrt(np.maximum(variance, 0))


def choose_ucb(features, observed, observed_values, beta=1.96):
    if len(observed) >= len(features):
        raise ValueError('Candidate pool exhausted')
    mean, sd = gp_posterior(features, observed, observed_values)
    scores = mean + beta * sd
    scores[np.asarray(observed, int)] = -np.inf
    return int(np.argmax(scores)), mean, sd


def campaign(features, target, seed, budget=16, initial_count=4, method='latent_gp_ucb'):
    target = np.asarray(target, float)
    if not (1 <= initial_count <= budget <= len(target)) or not np.all(np.isfinite(target)):
        raise ValueError('Invalid campaign budget or target')
    rng = np.random.default_rng(seed)
    initial = rng.choice(len(target), size=initial_count, replace=False).tolist()
    selected, values, rows = [], [], []
    for step in range(budget):
        mu = sd = acquisition = None
        if step < initial_count:
            choice = initial[step]
        elif method == 'random':
            choice = int(rng.choice(np.setdiff1d(np.arange(len(target)), selected)))
        elif method in ('latent_gp_ucb', 'descriptor_gp_ucb'):
            choice, mean, sigma = choose_ucb(features, selected, values)
            mu, sd = float(mean[choice]), float(sigma[choice])
            acquisition = mu + 1.96 * sd
        else:
            raise ValueError('Unknown campaign method')
        # The full target array is accessed only at the selected oracle index.
        value = float(target[choice])
        selected.append(choice)
        values.append(value)
        rows.append({'seed': seed, 'method': method, 'step': step + 1,
                     'candidate_index': choice, 'objective_negative_hydration_kcal_mol': value,
                     'best_observed_objective': max(values), 'observed_count_before_selection': step,
                     'posterior_mean_selected': mu, 'posterior_sd_selected': sd, 'UCB_selected': acquisition,
                     'initial_point': step < initial_count})
    return rows


def run_study(out=OUT, epochs=60, seeds=(20260928, 20260929, 20260930), pilot=False):
    started = time.perf_counter()
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(1)
    torch.use_deterministic_algorithms(True)
    records = load_records()
    graphs = [molecular_graph(row['smiles']) for row in records]
    X = np.stack([descriptors(row['smiles']) for row in records])
    y = np.array([row['expt_kcal_mol'] for row in records])
    partitions = {name: np.array([i for i, row in enumerate(records) if row['split'] == name])
                  for name in ('train', 'validation', 'test')}
    train, validation, test = (partitions[name] for name in ('train', 'validation', 'test'))
    scaler = StandardScaler().fit(X[train])
    Xstd = scaler.transform(X)
    predictions, losses, fits, metrics, pred_arrays, primary_latent = [], [], [], [], {}, None
    for seed in seeds:
        pred, latent, logs, info = fit_mpnn(graphs, y, train, validation, seed, epochs, out)
        pred_arrays[info['model']] = pred
        losses.extend(logs)
        fits.append(info)
        if primary_latent is None:
            primary_latent = latent
    shuffled_seed = 20261001
    pred, _, logs, info = fit_mpnn(graphs, y, train, validation, shuffled_seed, epochs, out, shuffled=True)
    pred_arrays[info['model']] = pred
    losses.extend(logs)
    fits.append(info)
    ridge_trials = []
    for alpha in (.1, 1., 10., 100.):
        ridge = Ridge(alpha=alpha).fit(Xstd[train], y[train])
        ridge_trials.append((regression_metrics(y[validation], ridge.predict(Xstd[validation]))['RMSE_kcal_mol'], alpha, ridge))
    _, best_alpha, ridge = min(ridge_trials, key=lambda row: (row[0], row[1]))
    pred_arrays['descriptor_ridge'] = ridge.predict(Xstd)
    pred_arrays['training_mean'] = np.full(len(y), y[train].mean())
    dump(out / 'descriptor_ridge.json', {'descriptor_names': DESCRIPTOR_NAMES,
        'training_scaler_mean': scaler.mean_.tolist(), 'training_scaler_scale': scaler.scale_.tolist(),
        'alpha_selected_on_validation': best_alpha, 'alpha_trials': [{'alpha': a, 'validation_RMSE': m} for m, a, _ in ridge_trials],
        'coef': ridge.coef_.tolist(), 'intercept': float(ridge.intercept_)})
    for model, pred in pred_arrays.items():
        for split, ids in partitions.items():
            metrics.append({'model': model, 'split': split, 'count': len(ids), **regression_metrics(y[ids], pred[ids])})
        predictions.extend({'model': model, 'molecule_id': row['molecule_id'], 'split': row['split'],
                            'observed_hydration_kcal_mol': y[i], 'predicted_hydration_kcal_mol': float(pred[i])}
                           for i, row in enumerate(records))
    latent_scaler = StandardScaler().fit(primary_latent[train])
    latent_std = latent_scaler.transform(primary_latent)
    dump(out / 'encoder_scaler.json', {'predeclared_seed': int(seeds[0]), 'mean': latent_scaler.mean_.tolist(),
                                      'scale': latent_scaler.scale_.tolist(), 'fit_ids': train.tolist()})
    embedding_rows = [{'molecule_id': row['molecule_id'], 'split': row['split'],
                       **{f'z{j}': float(value) for j, value in enumerate(primary_latent[i])}}
                      for i, row in enumerate(records)]
    budget = min(16, len(test))
    campaign_rows, campaign_checks = [], []
    al_seeds = range(2 if pilot else 8)
    for seed in al_seeds:
        initials = []
        for method, features in [('latent_gp_ucb', latent_std[test]), ('descriptor_gp_ucb', Xstd[test]),
                                 ('random', Xstd[test])]:
            rows = campaign(features, -y[test], seed, budget=budget, method=method)
            for row in rows:
                row['molecule_id'] = records[test[row['candidate_index']]]['molecule_id']
                row['simple_regret_kcal_mol'] = float(max(-y[test]) - row['best_observed_objective'])
            initials.append([row['candidate_index'] for row in rows[:4]])
            campaign_rows.extend(rows)
            campaign_checks.append({'seed': seed, 'method': method, 'budget': budget,
                                    'oracle_calls': len(rows), 'unique_calls': len(set(row['candidate_index'] for row in rows)),
                                    'test_pool_only': all(row['molecule_id'] in {records[i]['molecule_id'] for i in test} for row in rows)})
        assert initials[0] == initials[1] == initials[2]
    curves = []
    for method in ('latent_gp_ucb', 'descriptor_gp_ucb', 'random'):
        for step in range(1, budget + 1):
            selected = [row for row in campaign_rows if row['method'] == method and row['step'] == step]
            regret = np.array([row['simple_regret_kcal_mol'] for row in selected])
            curves.append({'method': method, 'step': step, 'seed_count': len(selected),
                           'mean_simple_regret_kcal_mol': float(regret.mean()),
                           'sd_simple_regret_kcal_mol': float(regret.std(ddof=1)),
                           'optimum_found_count': int(sum(regret < 1e-10))})
    paired = []
    for left, right in [('latent_gp_ucb', 'random'), ('latent_gp_ucb', 'descriptor_gp_ucb'), ('descriptor_gp_ucb', 'random')]:
        differences = []
        for seed in al_seeds:
            a = next(row for row in campaign_rows if row['seed'] == seed and row['method'] == left and row['step'] == budget)
            b = next(row for row in campaign_rows if row['seed'] == seed and row['method'] == right and row['step'] == budget)
            differences.append(b['simple_regret_kcal_mol'] - a['simple_regret_kcal_mol'])
        delta = np.array(differences)
        half = float(student_t.ppf(.975, len(delta) - 1) * delta.std(ddof=1) / math.sqrt(len(delta)))
        paired.append({'left': left, 'right': right, 'budget': budget, 'n': len(delta),
                       'mean_regret_reduction_left_vs_right': float(delta.mean()),
                       'lower_95pct_t': float(delta.mean() - half), 'upper_95pct_t': float(delta.mean() + half),
                       'left_wins': int(sum(delta > 1e-10)), 'ties': int(sum(abs(delta) <= 1e-10)),
                       'left_losses': int(sum(delta < -1e-10))})
    for filename, rows in [('predictions.csv', predictions), ('training_losses.csv', losses),
                           ('regression_metrics.csv', metrics), ('learned_embeddings.csv', embedding_rows),
                           ('campaign_evaluations.csv', campaign_rows), ('campaign_learning_curves.csv', curves),
                           ('campaign_paired_comparisons.csv', paired)]:
        write_csv(out / filename, rows)
    dump(out / 'campaign_checks.json', campaign_checks)
    dump(out / 'splits.json', {name: {'row_indices': ids.tolist(), 'molecule_ids': [records[i]['molecule_id'] for i in ids],
                                     'groups': sorted(set(records[i]['group_key'] for i in ids))} for name, ids in partitions.items()})
    summary = {
        'schema_version': 1, 'evidence': 'supervised hydration benchmark and retrospective frozen-pool search; not electrochemistry',
        'pilot': pilot, 'target': 'Experimental hydration free energy, kcal/mol',
        'dataset_count': len(records), 'split_counts': {name: len(ids) for name, ids in partitions.items()},
        'dataset_provenance_sha256': sha(DATA / 'provenance.json'),
        'selected_data_sha256': sha(DATA / 'selected_freesolv.csv'),
        'source_sha256': {str(Path(__file__).relative_to(ROOT)).replace('\\', '/'): sha(__file__)},
        'training': {'seeds': list(seeds), 'shuffled_seed': shuffled_seed, 'epochs_per_fit': epochs,
                     'optimizer': 'Adam(lr=0.003, weight_decay=1e-4), batch=32, grad_norm_clip=5',
                     'hidden_dim': 24, 'message_steps': 2, 'latent_dim': 48,
                     'selection': 'lowest validation MSE epoch; fixed seeds; test labels never select checkpoints',
                     'shuffled_control': 'training labels permuted only; validation labels kept true for checkpoint selection',
                     'fits': fits},
        'architecture': 'Directed reciprocal bond edges, edge-conditioned full matrices, GRU, sum+mean invariant pooling',
        'atom_output': 'No atom-reactivity/Fukui head: no atom-level labels available',
        'baseline': '10 RDKit descriptors, training-only StandardScaler, ridge alpha selected on validation',
        'counts': {'neural_training_runs': len(fits), 'training_epochs_executed': len(losses),
                   'prediction_rows': len(predictions), 'campaigns': len(campaign_checks),
                   'cached_oracle_calls': len(campaign_rows), 'new_experiments': 0, 'new_DFT': 0},
        'active_learning': {'pool': 'test partition only', 'pool_count': len(test),
                            'objective': 'maximize negative experimental hydration free energy',
                            'budget': budget, 'initial_shared_points': 4, 'seeds': list(al_seeds),
                            'latent_encoder_predeclared_seed': int(seeds[0]),
                            'encoder_and_scaler_label_access': 'train; validation only selects checkpoint; no test labels',
                            'GP': 'Fixed Matern-5/2, dimension-normalized distances, unit amplitude and lengthscale, noise variance=1e-5 standardized units',
                            'conditioning': 'Only already selected oracle values; outcome mean/std refitted on observed values',
                            'acquisition': 'UCB=posterior mean+1.96*posterior latent sd; deterministic lowest-index tie',
                            'uncertainty': 'Conditional GP latent sd; neither experimental uncertainty nor calibrated confidence',
                            'extra_pretraining_cost': 'Learned feature method consumes supervised train labels and validation checkpoint choices, not counted as pool calls; descriptor method has no such label pretraining'},
        'limitations': ['Single small split; acyclic groups are individual structures, not novel scaffold families',
                        'No electrochemical measurements or atom labels; no oxidation or reactivity validation',
                        'Three neural seeds are not a broad hyperparameter search or independent data replicates',
                        'Fixed-pool eight-seed t intervals are descriptive and have no multiplicity correction',
                        'Experimental uncertainty metadata is preserved but the unweighted point-label loss does not model it',
                        'Frozen learned features plus GP is not jointly trained end-to-end deep kernel learning'],
        'runtime_seconds': time.perf_counter() - started,
        'versions': {'python': platform.python_version(), 'numpy': np.__version__, 'torch': torch.__version__,
                     'rdkit': rdkit.__version__, 'scipy': scipy.__version__, 'sklearn': sklearn.__version__},
        'threads': {'torch': torch.get_num_threads(), **{k: os.environ.get(k) for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']}},
        'method_sources': ['https://proceedings.mlr.press/v70/gilmer17a.html',
                           'https://docs.pytorch.org/docs/stable/generated/torch.nn.GRUCell.html',
                           'https://scikit-learn.org/stable/common_pitfalls.html'],
        'output_sha256': {p.name: sha(p) for p in sorted(out.iterdir()) if p.suffix in ('.csv', '.pt', '.json') and p.name != 'summary.json'}}
    dump(out / 'summary.json', summary)
    print(json.dumps({'dataset': len(records), 'split_counts': summary['split_counts'],
                      'counts': summary['counts'], 'test_metrics': [m for m in metrics if m['split'] == 'test'],
                      'runtime_seconds': summary['runtime_seconds']}, indent=2))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare-from', type=Path, help='Prepare 256 rows from local SAMPL.csv, matching official snapshot')
    parser.add_argument('--pilot', action='store_true', help='Two-epoch smoke training; separate output subdirectory')
    args = parser.parse_args()
    if args.prepare_from:
        prepare_dataset(args.prepare_from)
    if args.pilot:
        run_study(out=OUT / 'pilot', epochs=2, seeds=(20260928,), pilot=True)
    else:
        run_study()


if __name__ == '__main__':
    main()
