"""Reproducible force-field minima descriptors and explicit-map reaction differences.

These are geometry/force-field calculations, not solution free energies, redox
potentials, reaction mechanisms, or experimental predictions.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdFreeSASA, rdMolAlign, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "electrograph/results/structure"
R_KCAL_MOL_K = 8.31446261815324 / 4184.0
TARGET = "COc1ccc2[nH]c(cc2c1)c3ccccc3"
SOURCE_REACTION = TARGET + ".CSc1ccccc1>>COc1ccc2[nH]c(c(Sc3ccccc3)c2c1)c4ccccc4"
MAPPED_EXAMPLE = "[CH3:1][CH:2]([H:4])[O:3][H:5]>>[CH3:1][CH:2]=[O:3].[H:4][H:5]"
MOLECULES = [
    ("target", "5-methoxy-2-phenyl-1H-indole", TARGET),
    ("melatonin", "Melatonin", "CC(=O)NCCC1=CNc2c1cc(OC)cc2"),
    ("caffeine", "Caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("phenylquinoline", "2-Phenylquinoline", "c1ccc(cc1)-c2ccc3ccccc3n2"),
    ("tryptophol", "Tryptophol", "OCCc1c[nH]c2ccccc12"),
    ("indoline", "Indoline", "c1ccc2c(c1)CCN2"),
    ("benzofuran_ester", "Ethyl benzofuran-2-carboxylate", "CCOC(=O)c1cc2ccccc2o1"),
]
SOURCES = [
    {"url": "https://www.rdkit.org/docs/RDKit_Book.html", "use": "ETKDG and pruning parameters"},
    {"url": "https://www.rdkit.org/docs/source/rdkit.Chem.rdForceFieldHelpers.html", "use": "MMFF optimization; installed API docstring specifies status 0 as converged"},
    {"url": "https://www.rdkit.org/docs/source/rdkit.Chem.rdFreeSASA.html", "use": "CalcSASA explicit radii, confIdx, algorithm and probe radius"},
    {"url": "https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html", "use": "PeriodicTable.GetRvdw and atom properties"},
    {"url": "https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html", "use": "Symmetry-aware molecular alignment"},
    {"url": "https://www.rdkit.org/docs/Overview.html", "use": "RDKit BSD software license; online documentation CC BY-SA 4.0"},
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError("Cannot export an empty table without a schema")
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_molecule(smiles: str, preserve_hydrogens: bool = False) -> Chem.Mol:
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("SMILES must be nonempty")
    params = Chem.SmilesParserParams()
    params.removeHs = not preserve_hydrogens
    mol = Chem.MolFromSmiles(smiles, params)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError("Invalid SMILES")
    return mol


def minima_weights(energy_kcal_mol, temperature_K: float) -> np.ndarray:
    """Unit-degeneracy Boltzmann-like weights of sampled MMFF minima only."""
    e = np.asarray(energy_kcal_mol, dtype=float)
    if e.ndim != 1 or e.size == 0 or not np.isfinite(e).all():
        raise ValueError("Energies must be a finite, nonempty vector")
    if not np.isfinite(temperature_K) or temperature_K <= 0:
        raise ValueError("Temperature must be positive and finite")
    factors = np.exp(-(e - np.min(e)) / (R_KCAL_MOL_K * temperature_K))
    return factors / factors.sum()


def radius_of_gyration(positions, masses) -> float:
    points = np.asarray(positions, dtype=float)
    mass = np.asarray(masses, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or mass.shape != (len(points),):
        raise ValueError("Positions and masses must have shapes (n,3), (n,)")
    if not len(points) or not np.isfinite(points).all() or not np.isfinite(mass).all() or np.any(mass <= 0):
        raise ValueError("Coordinates must be finite and masses positive")
    center = np.average(points, axis=0, weights=mass)
    return float(np.sqrt(np.average(np.sum((points - center) ** 2, axis=1), weights=mass)))


def vdw_radii(mol: Chem.Mol) -> list[float]:
    values = [float(Chem.GetPeriodicTable().GetRvdw(a.GetAtomicNum())) for a in mol.GetAtoms()]
    if not values or not np.isfinite(values).all() or min(values) <= 0:
        raise ValueError("Every atom must have a positive explicit vdW radius")
    return values


def heavy_atom_rmsd(mol: Chem.Mol, first: int, second: int) -> float:
    """Symmetry-aware aligned RMSD; work on copies to preserve saved coordinates."""
    heavy = Chem.RemoveHs(Chem.Mol(mol))
    return float(rdMolAlign.GetBestRMS(Chem.Mol(heavy), heavy, prbId=first, refId=second, maxMatches=10000))


def sasa_rotation_average(mol: Chem.Mol, conf_id: int, radii: list[float], rotations: int = 24) -> tuple[dict, list[float]]:
    """Reduce finite surface-quadrature orientation artifacts, retaining their range.

    Lee-Richards slice resolution is the installed FreeSASA default. Its numerical
    quadrature is not exactly rotation invariant. The fixed SO(3) sample is an
    integration diagnostic, not a statistical conformer/physical confidence band.
    """
    if rotations < 2:
        raise ValueError("At least two orientations are needed")
    copy = Chem.Mol(mol)
    xyz = np.array(mol.GetConformer(conf_id).GetPositions())
    xyz -= xyz.mean(axis=0)
    conf = copy.GetConformer(conf_id)
    options = rdFreeSASA.SASAOpts(rdFreeSASA.SASAAlgorithm.LeeRichards, rdFreeSASA.SASAClassifier.Protor, 1.4)
    rng = np.random.default_rng(17381)
    values = []
    for _ in range(rotations):
        rotation, upper = np.linalg.qr(rng.normal(size=(3, 3)))
        rotation = rotation @ np.diag(np.sign(np.diag(upper)))
        rotation[:, 0] *= np.linalg.det(rotation)
        for i, point in enumerate(xyz @ rotation):
            conf.SetAtomPosition(i, point)
        values.append(float(rdFreeSASA.CalcSASA(copy, radii, confIdx=conf_id, opts=options)))
    stats = {"sasa_A2": float(np.mean(values)),
             "sasa_orientation_sd_A2": float(np.std(values, ddof=1)),
             "sasa_orientation_range_A2": float(np.ptp(values)),
             "sasa_first_half_vs_all_abs_delta_A2": float(abs(np.mean(values[:rotations//2]) - np.mean(values)))}
    return stats, values


def select_minima(mol: Chem.Mol, records: list[dict], threshold_A: float = 0.35) -> tuple[list[int], list[dict]]:
    if not np.isfinite(threshold_A) or threshold_A < 0:
        raise ValueError("RMSD threshold must be nonnegative")
    eligible = [r for r in records if r["optimization_status"] == 0 and np.isfinite(r["energy_kcal_mol"])]
    eligible.sort(key=lambda r: (r["energy_kcal_mol"], r["conformer_id"]))
    selected, assignments = [], []
    for row in eligible:
        cid = row["conformer_id"]
        rmsds = [heavy_atom_rmsd(mol, cid, kept) for kept in selected]
        nearest = int(np.argmin(rmsds)) if rmsds else None
        duplicate = nearest is not None and rmsds[nearest] <= threshold_A
        representative = selected[nearest] if duplicate else cid
        assignments.append({"conformer_id": cid, "representative_id": representative,
                            "rmsd_to_representative_A": rmsds[nearest] if duplicate else 0.0,
                            "selected": not duplicate})
        if not duplicate:
            selected.append(cid)
    return selected, assignments


def aggregate_minima(records: list[dict], selected: list[int], temperature_K: float) -> tuple[dict, list[dict]]:
    selected_records = {r["conformer_id"]: r for r in records}
    if not selected:
        return {"status": "no_converged_minima", "temperature_K": temperature_K}, []
    values = [selected_records[cid] for cid in selected]
    weights = minima_weights([r["energy_kcal_mol"] for r in values], temperature_K)
    summary = {"status": "success", "temperature_K": temperature_K,
               "unique_minima": len(values), "effective_minima_count": float(1 / np.sum(weights ** 2)),
               "minimum_energy_kcal_mol": min(r["energy_kcal_mol"] for r in values),
               "energy_span_kcal_mol": max(r["energy_kcal_mol"] for r in values) - min(r["energy_kcal_mol"] for r in values)}
    for key in ("sasa_A2", "mass_weighted_Rg_A", "unweighted_Rg_A"):
        x = np.array([r[key] for r in values])
        mean = float(weights @ x)
        summary["weighted_" + key] = mean
        summary["within_minima_sd_" + key] = float(np.sqrt(weights @ ((x - mean) ** 2)))
    rows = [{"conformer_id": r["conformer_id"], "temperature_K": temperature_K,
             "relative_energy_kcal_mol": r["energy_kcal_mol"] - summary["minimum_energy_kcal_mol"],
             "forcefield_minima_weight": float(w)} for r, w in zip(values, weights)]
    return summary, rows


def save_sdf(path: Path, mol: Chem.Mol, records: list[dict], phase: str) -> None:
    writer = Chem.SDWriter(str(path))
    for row in records:
        copy = Chem.Mol(mol)
        copy.SetProp("_Name", f"{row['molecule_id']}_{row['seed']}_{row['conformer_id']}")
        copy.SetProp("coordinate_phase", phase)
        for k, v in row.items():
            copy.SetProp(k, str(v))
        writer.write(copy, confId=row["conformer_id"])
    writer.close()


def profile_molecule(smiles: str, molecule_id: str, name: str, output_dir: Path,
                     seeds=(20260928, 20260929, 20260930), attempts: int = 24,
                     temperatures=(250.0, 298.15, 350.0), max_iters: int = 1000) -> dict:
    if attempts < 1 or max_iters < 1 or not seeds or len(set(seeds)) != len(seeds):
        raise ValueError("Positive sampling and iteration budgets and distinct seeds required")
    for t in temperatures:
        minima_weights([0], t)
    output_dir.mkdir(parents=True, exist_ok=True)
    base = parse_molecule(smiles)
    mol = Chem.AddHs(base)
    result = {"molecule_id": molecule_id, "name": name, "source_smiles": smiles,
              "canonical_smiles": Chem.MolToSmiles(base), "formula": rdMolDescriptors.CalcMolFormula(base),
              "formal_charge": Chem.GetFormalCharge(base),
              "radical_electrons": sum(a.GetNumRadicalElectrons() for a in base.GetAtoms()),
              "requested_conformers": attempts * len(seeds), "attempts_per_seed": attempts,
              "seeds": list(seeds), "records": [], "replicas": [], "summaries": [], "weights": [], "surface_quadrature": []}
    if len(Chem.GetMolFrags(mol)) != 1 or not AllChem.MMFFHasAllMoleculeParams(mol):
        result.update(status="unsupported", reason="Disconnected structure or incomplete MMFF94 parameters")
        return result
    pooled = Chem.Mol(mol)
    pooled.RemoveAllConformers()
    pooled_records, initial_records = [], []
    initial_pool = Chem.Mol(pooled)
    radii = vdw_radii(mol)
    classified = rdFreeSASA.classifyAtoms(Chem.Mol(mol))
    result["radii_audit"] = {"explicit_vdw_A": radii,
        "classifier_zero_count": sum(float(v) == 0 for v in classified) if classified is not None else None,
        "classifier_returned_none": classified is None,
        "method": "RDKit periodic-table vdW radii, all atoms including explicit H",
        "probe_radius_A": 1.4, "algorithm": "LeeRichards", "rotation_average_count": 24, "rotation_seed": 17381}
    masses = [a.GetMass() for a in mol.GetAtoms()]
    for seed in seeds:
        replica = Chem.Mol(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = int(seed)
        params.numThreads = 1
        params.pruneRmsThresh = 0.15
        params.onlyHeavyAtomsForRMS = True
        params.useSymmetryForPruning = True
        params.trackFailures = True
        ids = list(AllChem.EmbedMultipleConfs(replica, numConfs=attempts, params=params))
        global_ids = []
        for cid in ids:
            gid = initial_pool.AddConformer(Chem.Conformer(replica.GetConformer(cid)), assignId=True)
            global_ids.append(gid)
            initial_records.append({"molecule_id": molecule_id, "seed": seed, "conformer_id": gid})
        optimized = list(AllChem.MMFFOptimizeMoleculeConfs(replica, numThreads=1, maxIters=max_iters,
                          mmffVariant="MMFF94", nonBondedThresh=100.0)) if ids else []
        replica_rows = []
        for cid, gid, (status, energy) in zip(ids, global_ids, optimized):
            assigned = pooled.AddConformer(Chem.Conformer(replica.GetConformer(cid)), assignId=True)
            if gid != assigned:
                raise RuntimeError("Initial and optimized coordinate IDs diverged")
            points = replica.GetConformer(cid).GetPositions()
            surface, rotation_values = sasa_rotation_average(replica, cid, radii)
            row = {"molecule_id": molecule_id, "seed": seed, "conformer_id": gid,
                   "optimization_status": int(status), "energy_kcal_mol": float(energy),
                   **surface,
                   "mass_weighted_Rg_A": radius_of_gyration(points, masses),
                   "unweighted_Rg_A": radius_of_gyration(points, np.ones(len(points)))}
            result["surface_quadrature"].extend({"molecule_id": molecule_id, "seed": seed, "conformer_id": gid,
                "orientation_index": i, "sasa_A2": value} for i, value in enumerate(rotation_values))
            if not all(np.isfinite(row[k]) for k in ("energy_kcal_mol", "sasa_A2", "mass_weighted_Rg_A", "unweighted_Rg_A")):
                raise RuntimeError("Nonfinite geometry result: retained run must be inspected")
            pooled_records.append(row)
            replica_rows.append(row)
        selected, assignments = select_minima(pooled, replica_rows)
        result["replicas"].append({"seed": seed, "requested": attempts,
            "returned_after_embedding_pruning": len(ids), "unreturned_pruned_or_failed": attempts-len(ids),
            "embedding_failure_event_counts": list(params.GetFailureCounts()),
            "converged": sum(r["optimization_status"] == 0 for r in replica_rows),
            "nonconverged": sum(r["optimization_status"] != 0 for r in replica_rows),
            "selected_minima": selected, "deduplication": assignments})
        for t in temperatures:
            summary, weights = aggregate_minima(pooled_records, selected, float(t))
            result["summaries"].append(dict(molecule_id=molecule_id, scope=f"seed_{seed}", **summary))
            result["weights"].extend(dict(molecule_id=molecule_id, scope=f"seed_{seed}", **w) for w in weights)
    selected, assignments = select_minima(pooled, pooled_records)
    result["pooled_deduplication"] = assignments
    result["pooled_selected_minima"] = selected
    for t in temperatures:
        summary, weights = aggregate_minima(pooled_records, selected, float(t))
        result["summaries"].append(dict(molecule_id=molecule_id, scope="pooled", **summary))
        result["weights"].extend(dict(molecule_id=molecule_id, scope="pooled", **w) for w in weights)
    result["seed_variability"] = []
    for t in temperatures:
        rows = [s for s in result["summaries"] if s["scope"].startswith("seed_") and s["temperature_K"] == t and s["status"] == "success"]
        if not rows:
            continue
        variation = {"molecule_id": molecule_id, "temperature_K": t, "successful_seeds": len(rows)}
        for key in ("minimum_energy_kcal_mol", "weighted_sasa_A2", "weighted_mass_weighted_Rg_A", "effective_minima_count"):
            x = [r[key] for r in rows]
            variation[key + "_mean"] = float(np.mean(x))
            variation[key + "_range"] = float(np.ptp(x))
            variation[key + "_sample_sd"] = float(np.std(x, ddof=1)) if len(x) > 1 else None
        result["seed_variability"].append(variation)
    result["records"] = pooled_records
    result["returned_after_embedding_pruning"] = len(pooled_records)
    result["converged"] = sum(r["optimization_status"] == 0 for r in pooled_records)
    result["nonconverged"] = len(pooled_records) - result["converged"]
    result["pooled_unique_minima"] = len(selected)
    result["status"] = "success" if selected else "no_converged_minima"
    save_sdf(output_dir / (molecule_id + "_embedded.sdf"), initial_pool, initial_records, "ETKDG_before_MMFF")
    save_sdf(output_dir / (molecule_id + "_optimized.sdf"), pooled, pooled_records, "MMFF94_after_optimization_status_retained")
    result["coordinate_sha256"] = {p.name: sha256(p) for p in output_dir.glob(molecule_id + "_*.sdf")}
    return result


def side_inventory(smiles: str) -> dict:
    mol = parse_molecule(smiles, preserve_hydrogens=True)
    with_h = Chem.AddHs(mol)
    elements = Counter(a.GetSymbol() for a in with_h.GetAtoms())
    return {"formula": rdMolDescriptors.CalcMolFormula(mol), "elements": dict(sorted(elements.items())),
            "formal_charge": Chem.GetFormalCharge(mol), "heavy_atoms": mol.GetNumHeavyAtoms()}


def reaction_inventory(reaction: str) -> dict:
    if reaction.count(">>") != 1:
        raise ValueError("Use exactly one reactants>>products separator")
    left, right = reaction.split(">>")
    reactants, products = side_inventory(left), side_inventory(right)
    elements = set(reactants["elements"]) | set(products["elements"])
    delta = {e: products["elements"].get(e, 0) - reactants["elements"].get(e, 0) for e in sorted(elements)}
    return {"reactants": reactants, "products": products, "element_delta_products_minus_reactants": delta,
            "charge_delta": products["formal_charge"] - reactants["formal_charge"],
            "elementally_balanced": all(v == 0 for v in delta.values()),
            "charge_balanced": products["formal_charge"] == reactants["formal_charge"],
            "heavy_atom_count_equal": products["heavy_atoms"] == reactants["heavy_atoms"]}


def _mapped_side(smiles: str) -> tuple[dict, dict]:
    mol = parse_molecule(smiles, preserve_hydrogens=True)
    atoms = {}
    for atom in mol.GetAtoms():
        num = atom.GetAtomMapNum()
        if num <= 0:
            raise ValueError("Every explicit atom must have a positive atom map")
        if num in atoms:
            raise ValueError("Duplicate atom map")
        atoms[num] = {"atomic_number": atom.GetAtomicNum(), "isotope": atom.GetIsotope(),
                      "formal_charge": atom.GetFormalCharge(), "implicit_or_bracket_H": atom.GetTotalNumHs(),
                      "radical_electrons": atom.GetNumRadicalElectrons(), "aromatic": atom.GetIsAromatic()}
    bonds = {}
    for bond in mol.GetBonds():
        pair = tuple(sorted((bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum())))
        bonds[pair] = {"order": bond.GetBondTypeAsDouble(), "aromatic": bond.GetIsAromatic()}
    return atoms, bonds


def mapped_bond_changes(reaction: str) -> dict:
    """Explicit-map topological CGR; no automapping, chemistry classification, or mechanism inference."""
    inventory = reaction_inventory(reaction)
    left, right = reaction.split(">>")
    ra, rb = _mapped_side(left)
    pa, pb = _mapped_side(right)
    if set(ra) != set(pa):
        raise ValueError("Atom-map sets differ between sides")
    for key in ra:
        if (ra[key]["atomic_number"], ra[key]["isotope"]) != (pa[key]["atomic_number"], pa[key]["isotope"]):
            raise ValueError("Mapped atom element or isotope identity mismatch")
    changes = []
    for pair in sorted(set(rb) | set(pb)):
        before, after = rb.get(pair), pb.get(pair)
        if before != after:
            changes.append({"atom_maps": list(pair), "before": before, "after": after,
                            "kind": "formed" if before is None else "broken" if after is None else "order_changed"})
    atom_changes = [{"atom_map": k, "before": ra[k], "after": pa[k]} for k in sorted(ra) if ra[k] != pa[k]]
    return {"status": "mapped_topology_computed", "reaction_smiles": reaction,
            "atom_count": len(ra), "bond_changes": changes, "atom_property_changes": atom_changes,
            "reaction_center_atom_maps": sorted({i for c in changes for i in c["atom_maps"]} | {a["atom_map"] for a in atom_changes}),
            "inventory": inventory, "mechanistic_or_reaction_class_claim": None}


def audit_target_identity() -> dict:
    mol = parse_molecule(TARGET)
    numbering = {"1": 6, "2": 7, "3": 8, "3a": 9, "4": 10, "5": 2, "6": 3, "7": 4, "7a": 5}
    edges = [("1", "2"), ("2", "3"), ("3", "3a"), ("3a", "7a"), ("7a", "1"),
             ("3a", "4"), ("4", "5"), ("5", "6"), ("6", "7"), ("7", "7a")]
    checks = {"indole_ring_edges_exist": all(mol.GetBondBetweenAtoms(numbering[a], numbering[b]) is not None for a, b in edges),
              "position_1_is_pyrrolic_NH": mol.GetAtomWithIdx(6).GetAtomicNum() == 7 and mol.GetAtomWithIdx(6).GetTotalNumHs() == 1,
              "methoxy_bonds_C5_O_CH3": mol.GetBondBetweenAtoms(2, 1) is not None and mol.GetBondBetweenAtoms(1, 0) is not None and mol.GetAtomWithIdx(1).GetAtomicNum() == 8 and mol.GetAtomWithIdx(0).GetTotalNumHs() == 3,
              "phenyl_bond_at_C2": mol.GetBondBetweenAtoms(7, 11) is not None and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in range(11, 17))}
    return {"name_verified_by_graph": "5-methoxy-2-phenyl-1H-indole", "checks": checks,
            "source_zero_based_atom_index_by_indole_position": numbering,
            "canonical_smiles": Chem.MolToSmiles(mol), "formula": rdMolDescriptors.CalcMolFormula(mol)}


def reaction_audits() -> dict:
    source = reaction_inventory(SOURCE_REACTION)
    return {"source_unmapped_reaction": {"reaction_smiles": SOURCE_REACTION,
            "status": "insufficient_atom_mapping", "inventory": source,
            "reaction_center": None, "mechanism": None,
            "sulfur_reagent_identity": "CSc1ccccc1 is thioanisole (methyl phenyl sulfide), not thiophenol (Sc1ccccc1)",
            "interpretation": "Unbalanced source input retained; no inferred mapping or replacement chemistry."},
            "educational_balanced_example": dict(mapped_bond_changes(MAPPED_EXAMPLE),
                 interpretation="Ethanol to acetaldehyde plus H2 is atom bookkeeping only, not an electrocatalytic mechanism or a prediction."),
            "equal_heavy_count_counterexample": dict(reaction_inventory("CC>>CO"), reaction_smiles="CC>>CO"),
            "target_identity": audit_target_identity()}


def run_study(output_dir: Path = DEFAULT_OUTPUT, pilot: bool = False) -> dict:
    started = time.perf_counter()
    output_dir.mkdir(parents=True, exist_ok=True)
    molecules = MOLECULES[:1] if pilot else MOLECULES
    seeds = (20260928,) if pilot else (20260928, 20260929, 20260930)
    attempts = 6 if pilot else 24
    profiles = []
    for mid, name, smiles in molecules:
        print(f"Profiling {mid}: {attempts} requested conformers x {len(seeds)} seed(s)", flush=True)
        profile = profile_molecule(smiles, mid, name, output_dir, seeds=seeds, attempts=attempts)
        profiles.append(profile)
        write_json(output_dir / (mid + ".json"), profile)
    tables = {"conformers": [r for p in profiles for r in p["records"]],
              "ensemble_statistics": [r for p in profiles for r in p["summaries"]],
              "minima_weights": [r for p in profiles for r in p["weights"]],
              "surface_quadrature": [r for p in profiles for r in p["surface_quadrature"]],
              "seed_variability": [r for p in profiles for r in p.get("seed_variability", [])]}
    for name, rows in tables.items():
        if rows:
            write_csv(output_dir / (name + ".csv"), rows)
    write_json(output_dir / "reaction_audit.json", reaction_audits())
    # Explicit unsupported and malformed examples do not enter the chemistry pool.
    negative_cases = []
    for label, smi in [("unsupported_metal", "[Fe]"), ("disconnected", "CC.O"), ("invalid", "C1(")]:
        try:
            value = profile_molecule(smi, label, label, output_dir, seeds=(1,), attempts=1)
            negative_cases.append({"case": label, "status": value["status"], "reason": value.get("reason")})
        except ValueError as exc:
            negative_cases.append({"case": label, "status": "rejected", "reason": str(exc)})
    write_json(output_dir / "input_rejections.json", negative_cases)
    summary = {"schema_version": 1, "pilot": pilot,
        "evidence": "Executed ETKDGv3/MMFF94 force-field calculations; not quantum chemistry or solution thermodynamics",
        "method": {"embedding": "ETKDGv3", "embedding_pruneRmsThresh_A": 0.15,
          "postoptimization_heavy_atom_RMSD_threshold_A": 0.35, "postoptimization_selection": "Ascending-energy greedy symmetry-aware RMSD, converged minima only",
          "optimization": "MMFF94", "max_iterations": 1000, "nonbonded_threshold_A": 100.0,
          "energy_unit": "kcal/mol", "R_kcal_mol_K": R_KCAL_MOL_K, "temperatures_K": [250.0, 298.15, 350.0],
          "weights": "exp(-(E-Emin)/(RT)), unit degeneracy per sampled distinct minimum; no vibrational, solvation or basin entropy",
          "probe_radius_A": 1.4, "SASA_algorithm": "LeeRichards averaged over 24 fixed SO(3) orientations",
          "SASA_rotation_seed": 17381, "SASA_resolution": "Installed FreeSASA default slice resolution; no claim of converged quadrature",
          "SASA_radii": "RDKit atomic vdW radii including explicit H",
          "Rg": "Mass weighted over all explicit atoms, with unweighted source-equivalent comparator",
          "seeds": list(seeds), "attempts_per_seed": attempts, "num_threads": 1,
          "embedding_count_boundary": "Requested is numConfs budget; unreturned includes pruning and possible failures; failure event counters are not counts of missing conformers"},
        "counts": {"molecules": len(profiles), "replica_runs": len(profiles)*len(seeds),
          "requested_conformers": sum(p["requested_conformers"] for p in profiles),
          "returned_after_embedding_pruning": sum(p.get("returned_after_embedding_pruning", 0) for p in profiles),
          "MMFF_minimizations": sum(len(p["records"]) for p in profiles),
          "converged": sum(p.get("converged", 0) for p in profiles), "nonconverged": sum(p.get("nonconverged", 0) for p in profiles),
          "pooled_unique_minima": sum(p.get("pooled_unique_minima", 0) for p in profiles),
          "ensemble_statistic_rows": len(tables["ensemble_statistics"]), "minima_weight_rows": len(tables["minima_weights"]),
          "SASA_orientation_evaluations": len(tables["surface_quadrature"]),
          "coordinate_files": 2*len(profiles), "negative_input_cases": len(negative_cases),
          "QM_calculations": 0, "experimental_measurements": 0},
        "molecule_summary": [{k: p[k] for k in ("molecule_id", "formula", "status", "requested_conformers", "returned_after_embedding_pruning", "converged", "nonconverged", "pooled_unique_minima")} for p in profiles],
        "limitations": ["Sparse conformer search with fixed protonation/tautomer and no explicit solvent/electrode.",
            "Minima counts and descriptor weights depend on RMSD threshold and seed; no basin-volume degeneracy correction.",
            "No Hessians: force-field convergence is an optimizer criterion, not proof of a true local minimum.",
            "Three seeds measure search variability only and are not calibrated uncertainty intervals.",
            "SASA depends on probe and radii; mass-weighted Rg is not buried volume.",
            "Finite SASA integration is orientation sensitive: 24-orientation averages reduce but do not eliminate it; raw values and half-versus-full means retained.",
            "CGR is explicit-map bond comparison, not automated atom mapping or reaction-mechanism inference."],
        "runtime_versions": {"python": platform.python_version(), "rdkit": rdBase.rdkitVersion, "numpy": np.__version__},
        "elapsed_seconds": time.perf_counter()-started, "sources": SOURCES,
        "license_note": "RDKit software BSD; documentation CC BY-SA 4.0. No external dataset/checkpoint downloaded or redistributed.",
        "source_sha256": {"electrograph/scripts/structure_reviewed.py": sha256(Path(__file__)),
                          "electrograph/source/electrograph_kmc_core.py": sha256(ROOT / "electrograph/source/electrograph_kmc_core.py")},
        "output_sha256": {p.name: sha256(p) for p in sorted(output_dir.iterdir()) if p.is_file() and p.name != "summary.json"}}
    write_json(output_dir / "summary.json", summary)
    print(json.dumps(summary["counts"], indent=2), flush=True)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    destination = args.output_dir or (DEFAULT_OUTPUT / "pilot" if args.pilot else DEFAULT_OUTPUT)
    run_study(destination, args.pilot)
