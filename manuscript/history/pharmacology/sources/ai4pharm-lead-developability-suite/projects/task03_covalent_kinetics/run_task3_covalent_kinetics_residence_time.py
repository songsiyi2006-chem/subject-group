#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Task 3: reproducible covalent-inhibitor kinetics and reactivity workflow.

Python >=3.10. Install: python -m pip install numpy scipy matplotlib rdkit pillow
Run: python run_task3_covalent_kinetics_residence_time.py
Optional: --quantum xtb --xtb /path/to/xtb --git-sync
All output defaults to this script's project directory; use --output-dir to change it.

Scientific contract
-------------------
The eight molecules are hypothetical, not marketed-drug structures. RDKit/YAeHMOP
extended Hueckel (EHT) orbital calculations are real low-level calculations, NOT
DFT. Optional GFN2-xTB replaces the electronic descriptors. Default kinetic
constants and the barrier model are explicitly illustrative, not fitted to
experimental observations. No transition-state search or clinical prediction
is claimed. The script preserves this provenance in reports, tables and plots.

Output: four required 300-dpi PNGs, bilingual Markdown reports, CSVs, SDF/XYZ,
parameter JSON, validation and provenance JSON, and an idempotent README index.
Git sync is opt-in, checks the repository and main branch, stages only generated
deliverables, commits, and performs an ordinary (never forced) push.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import logging
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Any

try:
    import numpy as np
    import scipy
    from scipy.integrate import solve_ivp
    from scipy.optimize import least_squares
    from scipy.stats import linregress
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from PIL import Image
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, Crippen, rdMolDescriptors, rdEHTTools
except ImportError as exc:
    raise SystemExit(
        f"Missing dependency: {exc}. Install with:\n"
        "python -m pip install numpy scipy matplotlib rdkit pillow\n"
        "RDKit must include the rdEHTTools / YAeHMOP extension."
    ) from exc

LOG = logging.getLogger("task3")
LN2 = math.log(2.0)
R = 8.31446261815324  # J mol^-1 K^-1
KB = 1.380649e-23     # J K^-1
PLANCK = 6.62607015e-34
CORE = "c1ccc(Nc2ncccn2)cc1"
COMMIT_MESSAGE = "feat(task3): covalent inhibitor kinetics, residence time & GSH reactivity calibration"
REFERENCES = [
    ("Flanagan et al. (2014), experimental GSH reactivity and computational methods",
     "https://doi.org/10.1021/jm501412a"),
    ("Bradshaw et al. (2015), reversible covalent inhibitors and tunable residence",
     "https://doi.org/10.1038/nchembio.1817"),
    ("RDKit: rdEHTTools / YAeHMOP API",
     "https://www.rdkit.org/docs/source/rdkit.Chem.rdEHTTools.html"),
    ("xTB: orbital properties and machine-readable JSON output",
     "https://xtb-docs.readthedocs.io/en/latest/properties.html"),
    ("Assay Guidance Manual: mechanism-of-action assays",
     "https://www.ncbi.nlm.nih.gov/books/NBK92001/"),
]


@dataclass(frozen=True)
class Warhead:
    id: str
    name: str
    smiles: str
    mechanism: str
    kon_M_inv_s: float
    koff_s: float
    kinact_s: float
    krev_s: float
    barrier_baseline_kJ_mol: float
    kgsh_reverse_s: float = 0.0
    kgsh_M_inv_s: float | None = None
    provenance: str = "Illustrative scenario; no experimental calibration"


def default_panel() -> list[Warhead]:
    """Atom map 901 is the reacting carbon; 902 is alpha for Michael acceptors.

    Cyanoacrylamide reversibility is a scenario assumption requiring measurement;
    it is not inferred solely from a substructure. The common arylaminopyrimidine
    core has no asserted potency against a specific protein.
    """
    return [
        Warhead("W01", "Acrylamide", "[CH2:901]=[CH:902]C(=O)N" + CORE,
                "Michael", 1e6, .10, .003, 0, 77),
        Warhead("W02", "Dimethylaminomethyl acrylamide", "CN(C)C/[CH:901]=[CH:902]/C(=O)N" + CORE,
                "Michael", 8e5, .08, .002, 0, 80),
        Warhead("W03", "Alpha-cyanoacrylamide", "[CH2:901]=[C:902](C#N)C(=O)N" + CORE,
                "reversible Michael", 6e5, .12, .012, .0002, 70, .0003),
        Warhead("W04", "Chloroacetamide", "Cl[CH2:901]C(=O)N" + CORE,
                "SN2", 5e5, .20, .018, 0, 69),
        Warhead("W05", "Vinyl sulfone", "[CH2:901]=[CH:902]S(=O)(=O)" + CORE,
                "Michael", 4e5, .30, .008, 0, 72),
        Warhead("W06", "Methacrylamide", "[CH2:901]=[C:902](C)C(=O)N" + CORE,
                "Michael", 2e5, .40, .0008, 0, 83),
        Warhead("W07", "Beta-methyl cyanoacrylamide", "C/[CH:901]=[C:902](C#N)/C(=O)N" + CORE,
                "reversible Michael", 9e5, .06, .006, .00003, 76, .0001),
        Warhead("W08", "Beta-isopropyl cyanoacrylamide", "CC(C)/[CH:901]=[C:902](C#N)/C(=O)N" + CORE,
                "reversible Michael", 1.2e6, .03, .003, .000005, 79, .00004),
    ]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def write_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(value, encoding="utf-8", newline="\n")
    temp.replace(path)


def write_json(path: Path, data: Any) -> None:
    write_text(path, json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def write_csv(path: Path, rows: list[dict]) -> None:
    require(bool(rows), f"Cannot write empty table: {path}")
    keys = list(dict.fromkeys(k for row in rows for k in row))
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    with temp.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)


def load_panel(path: Path | None) -> list[Warhead]:
    panel = default_panel()
    if path:
        config = json.loads(path.read_text(encoding="utf-8-sig"))
        require(isinstance(config, dict) and set(config) == {"warheads"},
                "Parameter JSON must contain exactly a 'warheads' object.")
        updates = config["warheads"]
        require(isinstance(updates, dict), "'warheads' must map IDs to parameter objects.")
        require(not (set(updates) - {w.id for w in panel}), "Unknown warhead ID.")
        allowed = {"kon_M_inv_s", "koff_s", "kinact_s", "krev_s", "barrier_baseline_kJ_mol",
                   "kgsh_reverse_s", "kgsh_M_inv_s", "provenance"}
        replaced = []
        for w in panel:
            u = updates.get(w.id, {})
            require(isinstance(u, dict) and not (set(u) - allowed), f"Invalid fields for {w.id}")
            if u:
                require(bool(str(u.get("provenance", "")).strip()),
                        f"{w.id}: parameter overrides need an explicit provenance string.")
            replaced.append(replace(w, **u))
        panel = replaced
    for w in panel:
        for field in ("kon_M_inv_s", "koff_s", "kinact_s", "barrier_baseline_kJ_mol"):
            val = getattr(w, field)
            require(isinstance(val, (int, float)) and math.isfinite(val) and val > 0,
                    f"{w.id}: {field} must be finite and positive.")
        for field in ("krev_s", "kgsh_reverse_s"):
            val = getattr(w, field)
            require(isinstance(val, (int, float)) and math.isfinite(val) and val >= 0,
                    f"{w.id}: {field} must be finite and nonnegative.")
        if w.kgsh_M_inv_s is not None:
            require(isinstance(w.kgsh_M_inv_s, (int, float)) and
                    math.isfinite(w.kgsh_M_inv_s) and w.kgsh_M_inv_s > 0,
                    f"{w.id}: kgsh_M_inv_s must be positive or null.")
    return panel


def prepare_molecule(w: Warhead, seed: int, conformers: int) -> tuple[Any, int, int | None]:
    base = Chem.MolFromSmiles(w.smiles)
    require(base is not None, f"Invalid SMILES for {w.id}")
    common = Chem.MolFromSmarts("c1ccc(Nc2ncccn2)cc1")
    require(base.HasSubstructMatch(common), f"Shared core missing in {w.id}")
    indices = {a.GetAtomMapNum(): a.GetIdx() for a in base.GetAtoms() if a.GetAtomMapNum()}
    require(901 in indices, f"No mapped reaction center in {w.id}")
    center, alpha = indices[901], indices.get(902)
    for a in base.GetAtoms():
        a.SetAtomMapNum(0)
    mol = Chem.AddHs(base)
    options = AllChem.ETKDGv3()
    options.randomSeed = seed
    options.numThreads = 1
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=conformers, params=options))
    require(bool(ids), f"3D embedding failed for {w.id}")
    require(AllChem.MMFFHasAllMoleculeParams(mol), f"MMFF parameters missing for {w.id}")
    results = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=1, maxIters=2000)
    valid = [(energy, cid) for cid, (status, energy) in zip(ids, results)
             if status == 0 and math.isfinite(energy)]
    require(bool(valid), f"No converged MMFF conformer for {w.id}")
    _, best = min(valid)
    conf = Chem.Conformer(mol.GetConformer(best))
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    mol.SetProp("_Name", w.id + " " + w.name)
    mol.SetProp("reaction_center_zero_based", str(center))
    return mol, center, alpha


def eht_descriptors(mol: Any) -> tuple[float, float, np.ndarray]:
    ok, result = rdEHTTools.RunMol(mol)
    require(ok, "YAeHMOP extended Hueckel calculation failed")
    require(result.numElectrons % 2 == 0, "EHT workflow requires a closed-shell molecule")
    nocc = result.numElectrons // 2
    energies = np.asarray(result.GetOrbitalEnergies(), dtype=float)
    populations = np.asarray(result.GetReducedChargeMatrix(), dtype=float)
    require(populations.shape == (mol.GetNumAtoms(), result.numOrbitals),
            "Unsupported EHT reduced charge matrix shape")
    local = populations[:, nocc]
    require(np.all(np.isfinite(local)) and abs(local.sum() - 1) < 1e-5,
            "Invalid EHT LUMO populations")
    return float(energies[nocc - 1]), float(energies[nocc]), local


def xtb_descriptors(mol: Any, executable: str, folder: Path,
                    timeout: int) -> tuple[float, float, np.ndarray]:
    """Neutral optimization + vertical electron attachment at that geometry.

    f+ = q(N) - q(N+1), with Mulliken charges, fixed nuclei and ALPB water.
    Raw outputs are retained. Failure never silently falls back to another method.
    """
    folder.mkdir(parents=True, exist_ok=True)
    xyz = folder / "input.xyz"
    Chem.MolToXYZFile(mol, str(xyz))
    control = "$write\n  json=true\n$end\n"
    env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

    def run(directory: Path, geometry: Path, charge: int, optimize: bool) -> dict:
        directory.mkdir(parents=True, exist_ok=True)
        shutil.copy2(geometry, directory / "molecule.xyz")
        write_text(directory / "xcontrol.inp", control)
        command = [executable, "molecule.xyz", "--gfn", "2", "--alpb", "water",
                   "--chrg", str(charge), "--uhf", "1" if charge == -1 else "0",
                   "--input", "xcontrol.inp"]
        if optimize:
            command += ["--opt", "tight"]
        # A fresh directory per invocation avoids accepting stale result files.
        run_result = subprocess.run(command, cwd=directory, env=env, capture_output=True,
                                    text=True, errors="replace", timeout=timeout, check=False)
        write_text(directory / "stdout.log", run_result.stdout)
        write_text(directory / "stderr.log", run_result.stderr)
        require(run_result.returncode == 0 and "normal termination of xtb" in
                (run_result.stdout + run_result.stderr).lower(),
                f"xTB failed; see {directory}")
        out = directory / "xtbout.json"
        require(out.is_file(), f"xTB did not generate JSON: {directory}")
        data = json.loads(out.read_text(encoding="utf-8"))
        return {re.sub(r"\s+", "", key).lower(): value for key, value in data.items()}

    with tempfile.TemporaryDirectory(prefix="run_", dir=folder) as temp:
        scratch = Path(temp)
        try:
            neutral_opt = scratch / "neutral_opt"
            run(neutral_opt, xyz, 0, True)
            optimized = neutral_opt / "xtbopt.xyz"
            require(optimized.is_file(), "xTB optimized geometry missing")
            neutral = run(scratch / "neutral_sp", optimized, 0, False)
            anion = run(scratch / "anion_sp", optimized, -1, False)
        finally:
            # Preserve diagnostics even when an external calculation fails.
            for name in ("neutral_opt", "neutral_sp", "anion_sp"):
                source = scratch / name
                if not source.exists():
                    continue
                dest = folder / name
                dest.mkdir(exist_ok=True)
                for f in source.iterdir():
                    if f.is_file() and f.suffix in {".json", ".log", ".xyz", ".inp"}:
                        shutil.copy2(f, dest / f.name)
    energies = np.asarray(neutral["orbitalenergies/ev"], dtype=float)
    occup = np.asarray(neutral["fractionaloccupation"], dtype=float)
    nocc = int(round(float(occup.sum()) / 2))
    require(0 < nocc < len(energies), "Invalid xTB occupation count")
    q0 = np.asarray(neutral["partialcharges"], dtype=float)
    q1 = np.asarray(anion["partialcharges"], dtype=float)
    require(q0.shape == q1.shape == (mol.GetNumAtoms(),), "xTB atom-count mismatch")
    local = q0 - q1
    require(np.all(np.isfinite(local)) and abs(local.sum() - 1) < .02,
            "xTB Fukui charge difference does not sum to one")
    return float(energies[nocc - 1]), float(energies[nocc]), local


def quantum_panel(panel: list[Warhead], args: argparse.Namespace, out: Path) -> tuple[list[dict], list[dict]]:
    records, atoms = [], []
    structures = out / "data_task3" / "structures"
    structures.mkdir(parents=True, exist_ok=True)
    with Chem.SDWriter(str(structures / "warhead_panel.sdf")) as writer:
        for index, w in enumerate(panel):
            LOG.info("Electronic descriptors: %s (%s)", w.id, args.quantum)
            mol, center, alpha = prepare_molecule(w, args.seed + index, args.conformers)
            writer.write(mol)
            Chem.MolToXYZFile(mol, str(structures / f"{w.id}.xyz"))
            if args.quantum == "xtb":
                homo, lumo, local = xtb_descriptors(mol, args.xtb,
                    out / "data_task3" / "quantum" / w.id, args.xtb_timeout)
                method = "GFN2-xTB ALPB water; vertical condensed Fukui f+"
            else:
                homo, lumo, local = eht_descriptors(mol)
                method = "EHT/YAeHMOP; frozen-orbital LUMO population proxy"
            require(math.isfinite(homo) and math.isfinite(lumo) and lumo > homo,
                    f"Invalid frontier orbital gap: {w.id}")
            # Convention: eta=(I-A)/2, S=1/(2 eta); omega=mu^2/(2 eta).
            eta, mu = (lumo - homo) / 2, (homo + lumo) / 2
            heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
            rank = 1 + sum(local[j] > local[center] for j in heavy if j != center)
            status = ("SN2 center; beta-carbon comparison not applicable" if alpha is None else
                      "beta greater than alpha in this descriptor" if local[center] > local[alpha] else
                      "beta localization NOT supported by this descriptor")
            if alpha is not None and local[center] < .01:
                status += "; weak center contribution (<0.01); localization not established"
            records.append(dict(id=w.id, name=w.name, smiles=w.smiles, mechanism=w.mechanism,
                electronic_method=method, HOMO_eV=homo, LUMO_eV=lumo, mu_eV=mu, eta_eV=eta,
                omega_eV=mu * mu / (2 * eta), softness_eV_inv=1 / (2 * eta),
                reaction_center_zero_based=center, alpha_zero_based=alpha,
                center_fplus_proxy=float(local[center]), center_softness_eV_inv=float(local[center] / (2 * eta)),
                alpha_fplus_proxy=None if alpha is None else float(local[alpha]),
                center_heavy_atom_rank=int(rank), localization_check=status,
                MW_g_mol=Descriptors.MolWt(mol), cLogP=Crippen.MolLogP(mol),
                TPSA_A2=rdMolDescriptors.CalcTPSA(mol), kinetic_provenance=w.provenance))
            for a in mol.GetAtoms():
                i = a.GetIdx()
                atoms.append(dict(id=w.id, atom_zero_based=i, element=a.GetSymbol(),
                    reaction_center=i == center, alpha_carbon=i == alpha,
                    fplus_proxy=float(local[i]), local_softness_eV_inv=float(local[i] / (2 * eta)),
                    electronic_method=method))
    return records, atoms


def slow_rate(w: Warhead, concentration: float) -> tuple[float, float]:
    """Stable exact eigenvalues for irreversible, constant-free-I binding."""
    a = w.kon_M_inv_s * concentration
    total = a + w.koff_s + w.kinact_s
    fast = .5 * (total + math.sqrt(max(0.0, total * total - 4 * a * w.kinact_s)))
    return a * w.kinact_s / fast, fast


def checked_solve(fun: Any, times: np.ndarray, y0: list | np.ndarray, **kwargs: Any) -> np.ndarray:
    solution = solve_ivp(fun, (float(times[0]), float(times[-1])), y0,
                         t_eval=times, method="BDF", rtol=2e-8, atol=2e-11, **kwargs)
    require(solution.success, f"ODE integration failed: {solution.message}")
    require(np.all(np.isfinite(solution.y)), "ODE returned non-finite states")
    require(solution.y.min() >= -2e-7, "ODE returned materially negative states")
    return solution.y


def concentration_kinetics(panel: list[Warhead], records: list[dict], count: int) -> tuple[list[dict], list[dict]]:
    """ODE-derived slow-phase kobs; hyperbolic regression is an approximation.

    Forward-only diagnostics (krev=0) permit like-for-like kinact comparison.
    Actual krev is used in dynamic exposure and washout, including cyano compounds.
    """
    samples, trajectories = [], []
    concentrations = np.geomspace(1e-9, 1e-4, count)
    for w, record in zip(panel, records):
        obs = []
        kd = w.koff_s / w.kon_M_inv_s
        ki = (w.koff_s + w.kinact_s) / w.kon_M_inv_s
        max_eigen_error = 0.0
        for dose in concentrations:
            slow, fast = slow_rate(w, float(dose))
            a = w.kon_M_inv_s * dose
            matrix = np.array([[-a, w.koff_s, 0],
                               [a, -w.koff_s - w.kinact_s, 0],
                               [0, w.kinact_s, 0]])
            times = np.unique(np.r_[0, np.geomspace(.001 / fast, 6 / slow, 160)])
            y = checked_solve(lambda t, state: matrix @ state, times, [1., 0., 0.], jac=matrix)
            require(np.max(abs(y.sum(axis=0) - 1)) < 1e-6, "Clamped-I target mass balance failed")
            survival = y[0] + y[1]  # unmodified target, not instantaneous enzyme activity
            mask = (times >= 12 / fast) & (survival > 1e-4)
            tfit, sfit = times[mask], survival[mask]
            require(len(tfit) >= 8, "Insufficient slow-phase samples")
            # Scale time and residuals; estimate positive amplitude and rate in log space.
            ts = (tfit - tfit[0]) * slow
            normalized = sfit / sfit[0]
            fit = least_squares(lambda p: np.exp(p[0] - np.exp(p[1]) * ts) - normalized,
                                [0., 0.], bounds=([-3., -5.], [3., 5.]),
                                xtol=1e-12, ftol=1e-12, gtol=1e-12)
            require(fit.success, "Nonlinear slow-phase regression failed")
            kobs = float(math.exp(fit.x[1]) * slow)
            error = abs(kobs / slow - 1)
            max_eigen_error = max(max_eigen_error, error)
            obs.append(kobs)
            samples.append(dict(id=w.id, concentration_M=float(dose), kobs_ODE_fit_s=kobs,
                kobs_exact_s=slow, kobs_rapid_equilibrium_s=w.kinact_s * dose / (kd + dose),
                kobs_QSSA_s=w.kinact_s * dose / (ki + dose),
                relative_fit_eigen_error=error, model="forward-only clamped free inhibitor"))
            for t, f, c, adduct in zip(times, *y):
                trajectories.append(dict(id=w.id, concentration_M=float(dose), time_s=float(t),
                    E_fraction=float(f), EI_fraction=float(c), covalent_fraction=float(adduct)))
        observed = np.asarray(obs)
        # Fit log(kobs): equal relative weighting across five orders of concentration.
        def residual(p: np.ndarray) -> np.ndarray:
            kmax, half_uM = np.exp(p)
            predicted = kmax * (concentrations * 1e6) / (half_uM + concentrations * 1e6)
            return np.log(predicted) - np.log(observed)
        fit = least_squares(residual, np.log([w.kinact_s, ki * 1e6]),
                            bounds=(np.log([1e-12, 1e-9]), np.log([1e3, 1e9])),
                            xtol=1e-12, ftol=1e-12, gtol=1e-12)
        require(fit.success, f"Hyperbola regression failed for {w.id}")
        kfit, kifit_uM = np.exp(fit.x)
        efficiency = w.kinact_s / kd
        category = ("weak/slow (<1e3)" if efficiency < 1e3 else
                    "intermediate (1e3-1e4)" if efficiency < 1e4 else
                    "requested screening band (1e4-1e6)" if efficiency <= 1e6 else
                    "above requested screening band (>1e6)")
        record.update(kon_M_inv_s=w.kon_M_inv_s, koff_s=w.koff_s, kinact_s=w.kinact_s,
            krev_s=w.krev_s, KD_M=kd, KI_QSSA_M=ki, kinact_over_KD_M_inv_s=efficiency,
            kinact_over_KI_QSSA_M_inv_s=w.kinact_s / ki,
            kinact_fit_s=float(kfit), KI_fit_M=float(kifit_uM * 1e-6),
            efficiency_fit_M_inv_s=float(kfit / (kifit_uM * 1e-6)),
            hyperbola_log_RMSE=float(np.sqrt(np.mean(fit.fun ** 2))),
            kinact_over_koff=w.kinact_s / w.koff_s,
            max_kobs_eigen_relative_error=max_eigen_error,
            efficiency_category=category, noncovalent_residence_s=1 / w.koff_s)
    return samples, trajectories


def reactivity(panel: list[Warhead], records: list[dict], args: argparse.Namespace) -> list[dict]:
    """Uncalibrated LFER barrier scenario + Eyring conversion, not a TS calculation.

    DeltaG = baseline + 2 kJ/mol/eV * (LUMO - panel mean LUMO), clipped to +/-4.
    k2_thiolate = (kBT/h)/Cstandard exp(-DeltaG/RT), with Cstandard=1 M.
    kGSH apparent at the stated pH = thiolate fraction * k2_thiolate.
    Overrides of kgsh bypass this model and retain their declared provenance.
    """
    mean_lumo = float(np.mean([d["LUMO_eV"] for d in records]))
    fraction = 1 / (1 + 10 ** (args.gsh_pka - args.ph))
    curves = []
    for w, d in zip(panel, records):
        barrier = w.barrier_baseline_kJ_mol + float(np.clip(2 * (d["LUMO_eV"] - mean_lumo), -4, 4))
        kthiolate = KB * args.temperature / PLANCK * math.exp(-1000 * barrier / (R * args.temperature))
        kgsh = w.kgsh_M_inv_s if w.kgsh_M_inv_s is not None else fraction * kthiolate
        forward = kgsh * args.gsh_mM * 1e-3
        relaxation = forward + w.kgsh_reverse_s
        plateau = forward / relaxation
        half_forward = LN2 / forward
        # Time to 50% total conversion may never occur for a reversible GSH adduct.
        half_actual = -math.log1p(-.5 / plateau) / relaxation if plateau > .5 else None
        d.update(barrier_surrogate_kJ_mol=barrier,
            barrier_provenance="Uncalibrated descriptor LFER scenario; not a transition-state barrier",
            thiolate_fraction=fraction, kgsh_M_inv_s=kgsh, kgsh_reverse_s=w.kgsh_reverse_s,
            kgsh_provenance=w.provenance if w.kgsh_M_inv_s is not None else "Eyring/LFER scenario; uncalibrated",
            gsh_forward_half_life_min=half_forward / 60,
            gsh_relaxation_half_life_min=LN2 / relaxation / 60,
            gsh_time_to_50pct_adduct_min=None if half_actual is None else half_actual / 60,
            gsh_equilibrium_adduct_fraction=plateau,
            safety_ratio_requested=d["kinact_over_KD_M_inv_s"] / kgsh,
            safety_ratio_QSSA=d["kinact_over_KI_QSSA_M_inv_s"] / kgsh,
            hyperreactive_screen_flag=half_forward < 1800,
            requested_screening_sweet_spot=(1e4 <= d["kinact_over_KD_M_inv_s"] <= 1e6 and half_forward >= 1800))
        times = np.unique(np.r_[np.linspace(0, 6 * 3600, 100), np.geomspace(1, 6 / relaxation, 100)])
        y = checked_solve(lambda t, z: [-forward * z[0] + w.kgsh_reverse_s * z[1],
                                       forward * z[0] - w.kgsh_reverse_s * z[1]], times, [1., 0.])
        exact = plateau * (-np.expm1(-relaxation * times))
        require(np.max(abs(y[1] - exact)) < 2e-6, "GSH ODE/analytic mismatch")
        require(np.max(abs(y.sum(axis=0) - 1)) < 1e-6, "GSH mass balance failure")
        for t, free, adduct in zip(times, *y):
            curves.append(dict(id=w.id, time_s=float(t), free_warhead_fraction=float(free),
                               GSH_adduct_fraction=float(adduct), analytic_adduct_fraction=float(
                                   plateau * (-np.expm1(-relaxation * t)))))
    efficiency = np.array([d["kinact_over_KD_M_inv_s"] for d in records])
    half = np.array([d["gsh_forward_half_life_min"] for d in records])
    for i, d in enumerate(records):
        dominated = np.any((efficiency >= efficiency[i]) & (half >= half[i]) &
                           ((efficiency > efficiency[i]) | (half > half[i])))
        d["pareto_nondominated"] = not bool(dominated)
    return curves


def dynamic_exposure(panel: list[Warhead], records: list[dict], args: argparse.Namespace) -> list[dict]:
    """Finite inhibitor bolus: target binding, GSH consumption and drug clearance.

    All concentrations are scaled by Etotal for numerical conditioning. GSH is
    buffered and constant. Target is conserved over this separate 6-hour assay.
    I+EI+C+SG+cleared is conserved, including reversible C -> EI and SG -> I.
    """
    rows = []
    e0 = args.target_nM * 1e-9
    clearance = LN2 / (args.drug_half_life_h * 3600)
    times = np.linspace(0, 6 * 3600, 241)
    for w, d in zip(panel, records):
        kg = d["kgsh_M_inv_s"] * args.gsh_mM * 1e-3
        for dose in np.geomspace(1e-9, 1e-4, 7):
            ratio = dose / e0
            def fun(t: float, y: np.ndarray) -> list[float]:
                f, c, a, drug, sg, lost = y
                bind = w.kon_M_inv_s * e0 * f * drug
                diss = w.koff_s * c
                cov, rev = w.kinact_s * c, w.krev_s * a
                gsh, back = kg * drug, w.kgsh_reverse_s * sg
                return [-bind + diss, bind - diss - cov + rev, cov - rev,
                        -bind + diss - gsh + back - clearance * drug,
                        gsh - back, clearance * drug]
            y = checked_solve(fun, times, [1, 0, 0, ratio, 0, 0])
            require(np.max(abs(y[:3].sum(axis=0) - 1)) < 2e-6, "Dynamic target mass balance failed")
            require(np.max(abs(y[1:].sum(axis=0) - ratio)) / max(1, ratio) < 2e-6,
                    "Dynamic inhibitor mass balance failed")
            for j, t in enumerate(times):
                rows.append(dict(id=w.id, initial_drug_M=float(dose), time_s=float(t),
                    E_free_M=float(y[0, j] * e0), EI_M=float(y[1, j] * e0),
                    covalent_M=float(y[2, j] * e0), free_drug_M=float(y[3, j] * e0),
                    GSH_adduct_M=float(y[4, j] * e0), cleared_drug_M=float(y[5, j] * e0)))
    return rows


def washout_model(w: Warhead, turnover_h: float, args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, float]:
    """2-hour clamped pulse; exact switch to a perfect free-drug sink at t=0.

    Normalized target total starts at 1. All states degrade at kdeg and new free
    target is synthesized at kdeg, so d(total)/dt=kdeg*(1-total).
    Bond reversal returns C to EI, allowing intracomplex re-covalentization.
    """
    kd = LN2 / (turnover_h * 3600)
    concentration = args.pulse_uM * 1e-6

    def rhs(dose: float):
        def fun(t: float, y: np.ndarray) -> list[float]:
            f, c, a = y
            bind = w.kon_M_inv_s * dose * f
            return [kd - kd * f - bind + w.koff_s * c,
                    bind - (w.koff_s + w.kinact_s + kd) * c + w.krev_s * a,
                    w.kinact_s * c - (w.krev_s + kd) * a]
        return fun
    pre = np.linspace(-7200, 0, 101)
    pulse = checked_solve(rhs(concentration), pre, [1, 0, 0])
    post = np.unique(np.r_[0, np.geomspace(.01, args.washout_hours * 3600, 360)])
    recovery = checked_solve(rhs(0), post, pulse[:, -1])
    y = np.column_stack([pulse[:, :-1], recovery])
    t = np.r_[pre[:-1], post]
    require(np.max(abs(y.sum(axis=0) - 1)) < 2e-6, "Washout target balance failed")
    # Mean first passage time out of the two bound states, including degradation.
    bound = np.array([[-w.koff_s - w.kinact_s - kd, w.krev_s],
                      [w.kinact_s, -w.krev_s - kd]])
    occupancy = float(pulse[1:, -1].sum())
    require(occupancy > 0, "Pulse produced zero occupancy")
    mean_life_s = float(np.ones(2) @ np.linalg.solve(-bound, pulse[1:, -1] / occupancy))
    return t, y, mean_life_s


def washout(panel: list[Warhead], args: argparse.Namespace) -> tuple[list[dict], list[dict], list[dict]]:
    control = Warhead("REV", "Reversible noncovalent control", "", "reversible noncovalent",
                     2e4, 2e-4, 0, 0, 0, provenance="Illustrative control; KD=10 nM")
    rows, summary, reduced = [], [], []
    for w in [control] + panel:
        for half in args.turnover_hours:
            t, y, mean = washout_model(w, half, args)
            post_mask = t >= 0
            post_t, inhibition = t[post_mask], (y[1] + y[2])[post_mask]
            initial = inhibition[0]
            hit = np.flatnonzero(inhibition <= .5 * initial)
            recovery_half = None
            if len(hit):
                j = int(hit[0])
                recovery_half = float(np.interp(.5 * initial,
                    inhibition[j-1:j+1][::-1], post_t[j-1:j+1][::-1]) / 3600) if j else 0.
            summary.append(dict(id=w.id, name=w.name, turnover_half_life_h=half,
                turnover_mean_lifetime_h=half / LN2, occupancy_at_washout=float(initial),
                effective_bound_mean_lifetime_h=mean / 3600,
                recovery_half_of_initial_inhibition_h=recovery_half,
                noncovalent_intrinsic_residence_h=1 / w.koff_s / 3600,
                intrinsic_covalent_bond_lifetime_h=None if w.krev_s == 0 else 1 / w.krev_s / 3600,
                bond_lifetime_status="no chemical reversal" if w.krev_s == 0 else "finite"))
            for j, time in enumerate(t):
                rows.append(dict(id=w.id, turnover_half_life_h=half, time_since_washout_h=float(time / 3600),
                    free_active_fraction=float(y[0, j]), EI_fraction=float(y[1, j]),
                    covalent_fraction=float(y[2, j]), inhibited_fraction=float(y[1, j] + y[2, j])))
            if w.id == "W01":
                # Explicit implementation of the requested reduced turnover equation.
                kd = LN2 / (half * 3600)
                kobs = w.kinact_s * args.pulse_uM * 1e-6 / (w.koff_s / w.kon_M_inv_s + args.pulse_uM * 1e-6)
                pre = np.linspace(-7200, 0, 101)
                free_pre = checked_solve(lambda time, z: [kd - (kd + kobs) * z[0]], pre, [1.])[0]
                post = np.linspace(0, args.washout_hours * 3600, 241)
                free_post = checked_solve(lambda time, z: [kd * (1-z[0])], post, [free_pre[-1]])[0]
                analytic = 1 - (1-free_pre[-1]) * np.exp(-kd * post)
                require(np.max(abs(free_post - analytic)) < 1e-6, "Reduced turnover analytic mismatch")
                for time, f in zip(np.r_[pre[:-1], post], np.r_[free_pre[:-1], free_post]):
                    reduced.append(dict(id=w.id, turnover_half_life_h=half,
                        time_since_washout_h=float(time/3600), free_active_fraction=float(f),
                        model="reduced rapid-equilibrium irreversible turnover"))
    return rows, summary, reduced


def save_figure(fig: Any, directory: Path, name: str) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    fig.savefig(directory / name, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    with Image.open(directory / name) as im:
        require(min(im.size) > 1000, f"Figure unexpectedly small: {name}")
        require(abs(im.info.get("dpi", (0, 0))[0] - 300) < 1, f"Figure DPI incorrect: {name}")


def make_figures(records: list[dict], samples: list[dict], recovery: list[dict],
                 args: argparse.Namespace, directory: Path) -> dict:
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
        "axes.spines.top": False, "axes.spines.right": False, "axes.titlesize": 11,
        "axes.labelsize": 10, "legend.fontsize": 8, "figure.dpi": 120,
        "axes.grid": True, "grid.alpha": .18, "savefig.dpi": 300})
    colors = plt.get_cmap("tab10")(np.arange(8))
    fig, axes = plt.subplots(2, 4, figsize=(15, 7.5), layout="constrained")
    for d, color, ax in zip(records, colors, axes.flat):
        data = [v for v in samples if v["id"] == d["id"]]
        dose = np.array([v["concentration_M"] for v in data])
        ax.semilogx(dose * 1e6, np.array([v["kobs_ODE_fit_s"] for v in data]) * 60,
                    "o", color=color, ms=3.5, label="ODE slow-phase fit")
        dense = np.geomspace(1e-9, 1e-4, 300)
        ax.semilogx(dense * 1e6, 60 * d["kinact_fit_s"] * dense / (d["KI_fit_M"] + dense),
                    color=color, label="Hyperbolic approximation")
        ax.set(title=d["id"] + " | " + d["name"].replace(" cyanoacrylamide", "\ncyanoacrylamide"),
               xlabel="Free inhibitor (µM)", ylabel=r"$k_{obs}$ (min$^{-1}$)")
        ax.text(.04, .95, f"Fitted $k_{{inact}}$ = {d['kinact_fit_s']:.3g} s⁻¹\n"
                f"Fitted $K_I$ = {d['KI_fit_M']*1e9:.1f} nM", transform=ax.transAxes,
                va="top", fontsize=8)
        ax.set_ylim(0, max(v["kobs_ODE_fit_s"] for v in data) * 75)
    axes.flat[0].legend(loc="lower right")
    fig.suptitle("Two-step covalent inactivation | illustrative kinetics\n"
                 "Forward-only diagnostic; fitted $K_I$ is not assumed equal to $K_D$", fontsize=14)
    save_figure(fig, directory, "fig1_kobs_concentration_hyperbola.png")

    fig, axes = plt.subplots(1, len(args.turnover_hours), figsize=(5 * len(args.turnover_hours), 5),
                             squeeze=False, layout="constrained")
    comparison = [("REV", "Reversible control", "#555555", "--"),
                  ("W01", "W01: irreversible acrylamide", colors[0], "-"),
                  ("W08", "W08: reversible cyanoacrylamide", colors[7], "-")]
    for half, ax in zip(args.turnover_hours, axes.flat):
        for wid, label, color, linestyle in comparison:
            data = [v for v in recovery if v["id"] == wid and v["turnover_half_life_h"] == half]
            ax.plot([v["time_since_washout_h"] for v in data],
                    [100*v["inhibited_fraction"] for v in data], label=label,
                    color=color, ls=linestyle, lw=2)
        ax.axvspan(-2, 0, color="#dbeafe", alpha=.8)
        ax.axvline(0, color="#333333", lw=.8, ls=":")
        ax.set(title=f"Protein turnover half-life: {half:g} h",
               xlabel="Time from washout (h)", ylabel="Target inhibition (%)",
               xlim=(-2, args.washout_hours), ylim=(-2, 103))
        ax.legend(loc="upper right")
    fig.suptitle(f"Cellular recovery | 2 h pulse at {args.pulse_uM:g} µM, then a perfect drug sink\n"
                 "Illustrative rates; occupancy = noncovalent complex + covalent adduct", fontsize=14)
    save_figure(fig, directory, "fig2_target_recovery_washout.png")

    fig, ax = plt.subplots(figsize=(9, 6.6), layout="constrained")
    x = np.array([d["kinact_over_KD_M_inv_s"] for d in records])
    y = np.array([d["gsh_forward_half_life_min"] for d in records])
    xmin, xmax = min(100, x.min() / 2), max(2e6, x.max()*2)
    ymin, ymax = min(1, y.min()/2), max(300, y.max()*2)
    ax.set(xscale="log", yscale="log", xlim=(xmin, xmax), ylim=(ymin, ymax))
    ax.add_patch(Rectangle((1e4, 30), 1e6-1e4, ymax-30, color="#ccebd8", alpha=.65, zorder=0))
    ax.axhline(30, color="#b91c1c", ls="--", lw=1, label="30 min heuristic screening threshold")
    frontier = sorted((d for d in records if d["pareto_nondominated"]), key=lambda d: d["kinact_over_KD_M_inv_s"])
    ax.plot([d["kinact_over_KD_M_inv_s"] for d in frontier],
            [d["gsh_forward_half_life_min"] for d in frontier], ":", color="#333333", label="Pareto frontier")
    for d, color in zip(records, colors):
        ax.scatter(d["kinact_over_KD_M_inv_s"], d["gsh_forward_half_life_min"], s=95,
                   color=color, marker="s" if d["krev_s"] else "o", edgecolors="#222222", lw=.7)
        offset = (-7, -13) if d["id"] == "W04" else (6, 7)
        ax.annotate(d["id"], (d["kinact_over_KD_M_inv_s"], d["gsh_forward_half_life_min"]),
                    xytext=offset, textcoords="offset points", fontsize=9,
                    ha="right" if d["id"] == "W04" else "left")
    ax.text(.98, .04, "Shading: requested screening 'sweet spot'\nSquares: reversible covalent scenarios\n"
            "Forward GSH half-life; longer is favored", ha="right", va="bottom", transform=ax.transAxes, fontsize=9)
    ax.set(xlabel=r"Requested on-target efficiency $k_{inact}/K_D$ (M$^{-1}$ s$^{-1}$)",
           ylabel="Forward GSH conjugation half-life (min)",
           title="On-target / GSH trade-off | illustrative screening scenarios\n"
                 "This index does not establish clinical viability or predict DILI")
    ax.legend(loc="lower left")
    save_figure(fig, directory, "fig3_gsh_reactivity_safety_radar.png")

    fig, ax = plt.subplots(figsize=(9, 6), layout="constrained")
    energies = np.array([d["LUMO_eV"] for d in records])
    logeff = np.log10(x)
    require(np.ptp(energies) > 1e-8, "LUMO values have no variance; cannot fit LFER")
    regression = linregress(energies, logeff)
    dense = np.linspace(energies.min()-.2, energies.max()+.2, 200)
    ax.plot(dense, regression.intercept + regression.slope*dense, color="#444444", lw=1.6,
            label=f"Exploratory OLS: slope {regression.slope:.2f} eV⁻¹; R²={regression.rvalue**2:.3f}")
    for d, color, xv, yv in zip(records, colors, energies, logeff):
        ax.scatter(xv, yv, s=90, color=color, marker="s" if d["krev_s"] else "o",
                   edgecolor="#222222", lw=.6)
        ax.annotate(d["id"], (xv, yv), xytext=(5, 6), textcoords="offset points")
    ax.set(xlabel=f"Calculated LUMO energy (eV; {args.quantum.upper()})",
           ylabel=r"$\log_{10}[(k_{inact}/K_D)/(1\ M^{-1}s^{-1})]$",
           title="Exploratory LFER | calculated orbitals vs illustrative kinetics\n"
                 "Mixed mechanisms and uncalibrated rates: no predictive inference")
    ax.legend(loc="best")
    save_figure(fig, directory, "fig4_warhead_lumo_reactivity_correlation.png")
    return dict(slope_per_eV=float(regression.slope), intercept=float(regression.intercept),
                R_squared=float(regression.rvalue**2), n=8,
                interpretation="Exploratory only: independent illustrative kinetic inputs, mixed mechanisms")


def markdown_table(records: list[dict]) -> str:
    lines = ["| ID | Warhead | LUMO (eV) | kinact/KD (M^-1 s^-1) | GSH forward t1/2 (min) | Ratio | Flag <30 min |",
             "|---|---|---:|---:|---:|---:|---|"]
    for d in records:
        lines.append(f"| {d['id']} | {d['name']} | {d['LUMO_eV']:.3f} | {d['kinact_over_KD_M_inv_s']:.3g} | "
                     f"{d['gsh_forward_half_life_min']:.2f} | {d['safety_ratio_requested']:.3g} | "
                     f"{'yes' if d['hyperreactive_screen_flag'] else 'no'} |")
    return "\n".join(lines)


def reports(records: list[dict], summary: list[dict], lfer: dict, args: argparse.Namespace, out: Path) -> None:
    table = markdown_table(records)
    refs = "\n".join(f"- [{title}]({url})" for title, url in REFERENCES)
    pareto = ", ".join(d["id"] for d in records if d["pareto_nondominated"])
    flags = ", ".join(d["id"] for d in records if d["hyperreactive_screen_flag"]) or "none / 无"
    unsupported = ", ".join(d["id"] for d in records if "NOT supported" in d["localization_check"]
                            or "not established" in d["localization_check"]) or "none / 无"
    method = records[0]["electronic_method"]
    maxerr = max(d["max_kobs_eigen_relative_error"] for d in records)
    provenance = "\n".join(f"- {d['id']}: kinetics = {d['kinetic_provenance']}; GSH = {d['kgsh_provenance']}"
                            for d in records)
    en = r"""# Task 3 — Covalent inhibitor kinetics, residence and GSH reactivity

## Scope and evidence status

This is a reproducible computational demonstration, not a calibrated drug discovery
prediction. Eight hypothetical warheads share an arylaminopyrimidine core; no target
protein, binding pose or measured kinetic dataset was supplied. They are not
sotorasib, osimertinib or ibrutinib. Orbital descriptors are calculated; default
binding and chemical rates are illustrative. The barrier model is an uncalibrated
descriptor surrogate, not a located transition state. No clinical efficacy,
therapeutic window, safety or DILI probability is established by these outputs.

## 3A — Electronic descriptors and thiolate barrier scenario

RDKit validates structures, tracks the reacting atom, samples ETKDGv3 conformers,
and chooses the lowest-energy converged MMFF conformer. EHT/YAeHMOP is a genuine
semiempirical orbital calculation, but its absolute energies are not DFT energies
and must not be mixed with xTB values in a calibrated correlation. The optional
xTB mode optimizes with GFN2-xTB/ALPB water and computes neutral/anion single points
at the same neutral geometry. Logs, geometries and JSON are retained.

Using the Koopmans-style approximation I≈−EHOMO and A≈−ELUMO:

$$\mu=(E_H+E_L)/2,\quad\eta=(E_L-E_H)/2,\quad
\omega=\mu^2/(2\eta),\quad S=1/(2\eta).$$

The hardness/softness convention is explicitly fixed above. EHT local f+ is the
atom-resolved LUMO population (a frozen-orbital proxy); xTB f+ is q(N)−q(N+1), a
vertical finite difference. Local softness is S f+. Neither is a direct validation
of a reaction pathway. The beta carbon is compared with alpha, and its rank among
all heavy atoms is exported. Negative population values are retained, not clipped
to fabricate localization. Chloroacetamide undergoes SN2 at the carbon bearing Cl;
calling that atom a Michael beta carbon would be incorrect. Neutral microstates,
one selected conformer and unverified anion binding limit these descriptors;
protonation of the dimethylamino group requires separate pH-dependent modeling.

The optional barrier scenario is fully specified:

$$\Delta G^\ddagger_{scenario}=B_{warhead}+
\mathrm{clip}[2\,(E_L-\overline{E_L}),-4,4]\quad\mathrm{kJ/mol}.$$

The coefficient is 2 kJ mol⁻¹ eV⁻¹. The eight baselines are declared in the parameter
JSON, not learned from data. For a 1 M standard concentration and transmission
coefficient one:

$$k_{thiolate}=(k_BT/h)/C^\circ\exp[-\Delta G^\ddagger/(RT)],\quad
f_{thiolate}=1/(1+10^{pK_a-pH}),\quad k_{GSH}=f_{thiolate}k_{thiolate}.$$

An explicit kGSH override bypasses this Eyring scenario. The retained barrier column
then remains a separate scenario and is not inferred from that measured rate.
Real barrier validation requires a solvent/microstate-consistent thiolate addition
or SN2 pathway, transition-state optimization, one appropriate imaginary frequency,
IRC confirmation, free-energy corrections and experimental calibration.

## 3B — Two-step inactivation and regression

For free target F, noncovalent complex B and covalent complex C:

$$E+I\rightleftharpoons B\rightleftharpoons C,$$
$$\dot F=-k_{on}IF+k_{off}B,$$
$$\dot B=k_{on}IF-(k_{off}+k_{inact})B+k_{rev}C,$$
$$\dot C=k_{inact}B-k_{rev}C.$$

For the irreversible limit krev=0, rapid pre-equilibrium gives
B/(F+B)=I/(KD+I), where KD=koff/kon. Thus for unmodified target U=F+B,
Udot=−kinact I U/(KD+I), yielding the requested saturating kobs equation.
This reduction needs constant free inhibitor, negligible depletion, and binding
equilibration faster than chemical inactivation. KD is not generally the fitted KI.
The quasi-steady-state kinetic constant is KI,QSSA=(koff+kinact)/kon. The exact slow
decay eigenvalue at clamped I is:

$$a=k_{on}I,\quad S_r=a+k_{off}+k_{inact},\quad
\lambda_{slow}=\frac{2ak_{inact}}{S_r+\sqrt{S_r^2-4ak_{inact}}}.$$

Its low-I slope is kinact/KI,QSSA. At high I it tends to kinact. An exact transient
is biexponential; a hyperbola is not generally its exact slow eigenvalue.

The code numerically solves full mass action at 1 nM–100 µM. After 12 fast-mode time
constants, nonlinear least squares fits A exp(−kobs t) to U. It then fits a positive
hyperbola in log-rate residuals, giving equal relative weighting across doses.
The CSV distinguishes KD, KI,QSSA, fitted KI, fitted kinact, efficiency definitions
and approximation error. These are noiseless synthetic trajectories, so statistical
confidence intervals would not represent experimental uncertainty and are omitted.
Acquisition times adapt to slow rates and can exceed realistic assay durations.
Experiments require fixed observation windows, replicates and identifiability checks.

For comparable forward reactivity, figure 1 sets krev=0 for all warheads; cyano
compounds' reversibility is restored in dynamic exposure and washout. The forward
efficiency does not quantify their equilibrium occupancy or residence. Instantaneous
inhibition is B+C, whereas U=F+B is the covalently unmodified fraction.

A separate finite-bolus simulation includes actual inhibitor depletion, clearance,
GSH trapping and reverse GSH release. With target normalized by Etotal,
I+ B+C+SG+cleared is conserved. This assay has no protein turnover; the washout
experiment below handles turnover separately. No substrate competition is included.
The requested efficiency categories (<10³, 10³–10⁴, 10⁴–10⁶, >10⁶ M⁻¹s⁻¹) are
screening labels, not universal thresholds for clinical optimality.

## 3C — Residence, target turnover and rapid washout

Affinity contributes to recognition and exposure-dependent occupancy; residence
alone does not dictate clinical efficacy. Pharmacodynamics also depend on unbound
exposure, target synthesis, target vulnerability, tissue distribution and selectivity.
Long-lived covalent occupancy can persist after free drug clearance, as demonstrated
experimentally for reversible covalent kinase inhibitors [Bradshaw et al., 2015].

All target states degrade at kdeg=ln2/t1/2,protein; synthesis supplies free target
at ksyn=kdeg Etotal,baseline. Accordingly:

$$\dot F=k_{syn}-k_{deg}F-k_{on}IF+k_{off}B,$$
$$\dot B=k_{on}IF-(k_{off}+k_{inact}+k_{deg})B+k_{rev}C,$$
$$\dot C=k_{inact}B-(k_{rev}+k_{deg})C.$$

The cells receive a 2-hour constant free-drug pulse, followed by an instantaneous
perfect sink for free inhibitor. Bound drug remains; released drug is immediately
removed. The experiment has no residual intracellular depot or external rebinding.
Reversible C→B can still re-form C within the bound complex.

The requested reduced irreversible model is also integrated:
Fdot=ksyn−(kdeg+kobs)F. After washout,
F(t)=Etotal−[Etotal−F(0)]exp(−kdeg t).
For a single reversible noncovalent state intrinsic residence is 1/koff, while the
bound cellular lifetime is 1/(koff+kdeg). An irreversible bond has infinite chemical
lifetime in this model; 1/kdeg is the mean protein/adduct lifetime, not a chemical
off-rate, and ln2/kdeg is the recovery half-life. For reversible covalent states,
1/krev is just the bond-opening lifetime and differs from complete target release.

The full bound-state matrix Q is used to compute the occupancy-weighted mean
first-passage lifetime after washout: tau=1ᵀ(−Q)⁻¹p_bound(0). This accounts for
bond reversal, re-covalentization, dissociation and turnover. The numerical recovery
half-time is measured relative to inhibition at washout; it is not the same as tau.
An unreached half-time is exported as null rather than guessed or extrapolated.

## 3D — GSH twin and screening trade-off

At buffered 5 mM GSH (configurable), the irreversible pseudo-first-order model is
I(t)=I0 exp(−kGSH[GSH]t), so t1/2=ln2/(kGSH[GSH]). The code numerically integrates
this twin and checks against its analytic solution. For reversible GSH adducts,
the adduct fraction is p[1−exp(−(a+b)t)], with a=kGSH[GSH], b=kreverse,GSH,
p=a/(a+b). It exports the forward half-life, relaxation half-life, equilibrium
fraction and actual time to 50% adduct. The latter may not exist if p≤0.5.

The requested safety ratio (kinact/KD)/kGSH is dimensionless; a QSSA ratio is also
provided. It compares two rate efficiencies, not concentrations, selectivity across
proteins, exposure or clinical safety. Figure 3 maximizes both on-target efficiency
and forward GSH half-life, marks the Pareto frontier, and shades the requested
10⁴–10⁶ M⁻¹s⁻¹ / ≥30 min screening region. Its historical filename contains 'radar',
but it is the requested Pareto scatter plot.

A forward GSH half-life below 30 min triggers a hyperreactivity screening flag.
This is a user-specified heuristic, not a validated idiosyncratic DILI classifier.
GSH conjugation can detoxify electrophiles; GSH depletion and protein modification
depend on exposure, regeneration, metabolism and compartment. A buffered-GSH
assay does not itself simulate cellular GSH depletion. GSH reactivity measurements
provide intrinsic-reactivity information [Flanagan et al., 2014], not clinical outcomes.

## Reproducibility and limits

No random noise is added to kinetic data. The seed controls conformer generation;
library versions and SHA-256 hashes are recorded. Numerical checks cover mass
balance, positivity, ODE/analytic agreement, limiting kinetics, turnover recovery
and PNG DPI. Determinism is expected within the same software build, not necessarily
bitwise across platforms. EHT and xTB remain approximate; no reaction barrier,
binding potency or rate is claimed as experimentally validated. Figure 4 is an
exploratory OLS relationship across mixed mechanisms, with independently specified
kinetic inputs; neither a strong nor weak correlation establishes causality.

Use `python run_task3_covalent_kinetics_residence_time.py --help` for all options.
Edit `data_task3/parameters_template.json` and pass `--parameters` to supply rates
with provenance. Partial overrides retain the other illustrative values; provenance
must describe this. GSH rate overrides are apparent rates at the configured assay
conditions and are not automatically pH-corrected. `--self-test-only` checks core
mathematics. `--git-sync` requires an existing main checkout and origin remote;
generated paths alone are staged and pushed, without force or branch switching.
The script can be rerun: its README section is replaced idempotently and artifacts
are overwritten only at their dedicated paths. Avoid concurrent runs to one folder.
"""
    en += f"\n## Results of this run\n\n- Electronic method: {method}.\n" \
          f"- Temperature: {args.temperature:g} K; pH {args.ph:g}; GSH pKa assumption {args.gsh_pka:g}; " \
          f"GSH {args.gsh_mM:g} mM.\n- Pareto candidates: {pareto}.\n- Hyperreactivity flags: {flags}.\n" \
          f"- Beta localization unsupported or weak: {unsupported}.\n- LFER R²: {lfer['R_squared']:.4f} (illustrative).\n" \
          f"- Maximum ODE kobs / eigenvalue relative error: {maxerr:.3g}.\n\n{table}\n\n" \
          f"### Input provenance by compound\n\n{provenance}\n\n## References\n\n{refs}\n"
    zh = r"""# 任务3：共价抑制剂动力学、靶点驻留与谷胱甘肽反应性

## 研究范围与证据等级

本工作是可复现的计算演示，并非经过实验校准的候选药物预测。八种假想分子共享
芳基氨基嘧啶核心；未提供明确靶蛋白、结合构象或实测动力学数据。这些结构不是
索托拉西布、奥希替尼或伊布替尼。轨道描述符来自实际低层级量子计算，默认结合
速率与化学反应速率则是演示参数。能垒使用未经校准的描述符替代模型，没有执行
过渡态搜索。本结果不确立临床疗效、治疗窗、安全性或药物性肝损伤概率。

## 3A：亲电性与硫负离子反应能垒情景

RDKit 检查结构，以原子映射定位反应碳，生成 ETKDGv3 构象，并从收敛的 MMFF
构象中选择最低能者。默认 EHT/YAeHMOP 扩展休克尔法是真实半经验轨道计算，
不是 DFT；绝对轨道能量不可与 xTB 数值混用。可选 GFN2-xTB/ALPB 水溶剂模式
先优化中性分子，再在相同几何上计算中性和阴离子单点，保存原始日志与 JSON。

采用 I≈−EHOMO、A≈−ELUMO 的轨道近似：

$$\mu=(E_H+E_L)/2,\quad\eta=(E_L-E_H)/2,\quad
\omega=\mu^2/(2\eta),\quad S=1/(2\eta).$$

这里明确使用上述硬度与软度约定。EHT 的局部 f+ 是原子分辨的 LUMO 布居，
属于冻结轨道近似；xTB 的 f+=q(N)−q(N+1) 来自垂直有限差分。局部软度为 S f+。
程序比较迈克尔受体 β 碳与 α 碳的局部值，并导出反应碳在全部重原子中的排名。
未支持 β 位局域化时如实报告；负布居不截断、不修饰。氯乙酰胺在连氯碳发生
SN2 取代，不适用迈克尔 β 碳描述。中性微观状态、单一选定构象及阴离子束缚
状态的不确定性限制了解释；二甲氨基的质子化需要另行建模。

能垒情景为：

$$\Delta G^\ddagger=B_{warhead}+\mathrm{clip}[2(E_L-\overline{E_L}),-4,4]
\quad\mathrm{kJ/mol}.$$

系数为 2 kJ mol⁻¹ eV⁻¹，各亲电反应基团的基线在参数 JSON 中公开，未由实验拟合。
采用 1 M 标准态、传递系数为 1 的 Eyring 换算：

$$k_{thiolate}=(k_BT/h)/C^\circ\exp[-\Delta G^\ddagger/(RT)],$$
$$f_{thiolate}=1/(1+10^{pK_a-pH}),\qquad k_{GSH}=f_{thiolate}k_{thiolate}.$$

显式输入 kGSH 可覆盖该模型，此时能垒列仍是独立演示情景，并非从实测速率推断。
真正的加成或 SN2 能垒验证需统一溶剂与微观状态、过渡态优化、单一恰当虚频、
IRC 路径核验、自由能修正及实验校准；本脚本不会将替代值冒充这些计算。

## 3B：两步共价失活与非线性回归

记游离靶点为 F、非共价复合物为 B、共价复合物为 C：

$$\dot F=-k_{on}IF+k_{off}B,$$
$$\dot B=k_{on}IF-(k_{off}+k_{inact})B+k_{rev}C,$$
$$\dot C=k_{inact}B-k_{rev}C.$$

不可逆极限 krev=0，且结合达到快速预平衡时，B/(F+B)=I/(KD+I)，其中
KD=koff/kon。于是未共价修饰靶点 U=F+B 满足 Udot=−kinact I U/(KD+I)，
得到题目中的饱和双曲线。该近似要求游离药物恒定、耗竭可忽略、结合平衡快于
化学失活。KD 一般不等于回归所得 KI；准稳态 KI,QSSA=(koff+kinact)/kon。

恒定游离浓度下的精确慢衰减特征值为：

$$a=k_{on}I,\quad S_r=a+k_{off}+k_{inact},\quad
\lambda_{slow}=\frac{2ak_{inact}}{S_r+\sqrt{S_r^2-4ak_{inact}}}.$$

低浓度斜率为 kinact/KI,QSSA，高浓度趋向 kinact。完整瞬态为双指数，慢特征值
通常不是严格双曲线。程序在 1 nM–100 µM 范围积分全质量作用方程，跳过 12 个
快模态时间常数，再以非线性最小二乘拟合 A exp(−kobs t)。第二阶段拟合正参数
双曲线，使用对数速率残差，使各浓度近似按相对误差等权。表格分别列出 KD、
KI,QSSA、拟合 KI、拟合 kinact、各效率定义及双曲线近似误差。

这些轨迹没有加入实验噪声，因而不将回归协方差冒充实验置信区间。模拟采样时间
随慢速率调整，可能远超实际实验时长；真实实验需要限定窗口、重复测量及参数
可辨识性检验。图1统一使用 krev=0 比较正向效率；细胞洗脱和动态浓度模拟恢复
各化合物的实际设定逆反应。可逆共价分子的正向效率不能替代平衡占有率或驻留。
瞬时抑制为 B+C，未共价修饰比例则为 F+B，两者不是同一实验读数。

另有有限初始剂量模型，同时包括药物耗竭、清除、GSH 捕获和逆向释放，并核验
I+B+C+SG+已清除药物的质量守恒。该六小时模型不含蛋白周转；周转由独立洗脱
模块处理。本任务未包含底物竞争。题设 <10³、10³–10⁴、10⁴–10⁶、>10⁶ M⁻¹s⁻¹
分级仅作筛选标签，不能将其中某一范围称为普遍适用的临床最优标准。

## 3C：靶点驻留、蛋白周转与快速洗脱

亲和力参与识别与浓度依赖的占有过程，驻留时间也不能独自决定临床疗效；游离
暴露、靶点合成、通路对占有率的敏感性、组织分布和选择性同样重要。共价占有
可在游离药物清除后持续，Bradshaw 等的可逆共价激酶抑制剂实验提供了实例。

全部靶点状态以 kdeg=ln2/t1/2,protein 降解，新蛋白以
ksyn=kdeg Etotal,baseline 合成到游离池：

$$\dot F=k_{syn}-k_{deg}F-k_{on}IF+k_{off}B,$$
$$\dot B=k_{on}IF-(k_{off}+k_{inact}+k_{deg})B+k_{rev}C,$$
$$\dot C=k_{inact}B-(k_{rev}+k_{deg})C.$$

先恒定游离药物处理两小时，随后瞬时将游离药物降为零。已结合药物保留；其后
释放的药物立即移除，不考虑细胞药物库与外源再结合。不过 C→B 后仍可在同一
复合物内部再次形成共价键。脚本也单独积分题设简化模型
Fdot=ksyn−(kdeg+kobs)F，洗脱后解析解为
F(t)=Etotal−[Etotal−F(0)]exp(−kdeg t)，并核验数值解。

单态非共价抑制剂本征驻留时间为 1/koff，含周转时的细胞结合寿命为
1/(koff+kdeg)。理想不可逆键的化学寿命无限；1/kdeg 是蛋白/加合物的平均寿命，
ln2/kdeg 才是周转控制的恢复半衰期。可逆共价分子的 1/krev 仅为键打开的平均
等待时间，并非完整解离驻留时间。程序通过结合态矩阵 Q 计算洗脱时初始占有
加权的平均首次离开时间 tau=1ᵀ(−Q)⁻¹p_bound(0)，包含解离、断键、再次成键和
蛋白降解。恢复半时指抑制降至洗脱初值的一半，未在窗口内达到则记为 null。

## 3D：GSH 共轭孪生模型与反应性权衡

缓冲恒定 GSH 默认 5 mM。不可逆拟一级动力学为
I(t)=I0 exp(−kGSH[GSH]t)，t1/2=ln2/(kGSH[GSH])。数值积分与解析结果互相核验。
可逆 GSH 加合物的比例为 p[1−exp(−(a+b)t)]，其中
a=kGSH[GSH]、b=kreverse,GSH、p=a/(a+b)。输出正向半衰期、弛豫半衰期、
平衡加合物比例及达到总量 50% 加合的时间；p≤0.5 时最后一个指标不存在。

题设安全比 (kinact/KD)/kGSH 无量纲，同时给出 QSSA 版本。该比值没有包含暴露、
蛋白组选择性或临床结局，因此并非真实治疗指数。图3以靶向效率和 GSH 正向
半衰期同时越大越好绘制 Pareto 前沿，阴影标记题设筛选区。虽然文件名保留 radar，
内容按要求采用散点图。正向半衰期低于 30 分钟触发高反应性筛选标记，而非经过
验证的特异质性 DILI 预测。GSH 共轭可能有解毒作用；耗竭与蛋白共价修饰取决于
暴露、再生、代谢及细胞区室。恒定 GSH 试验本身并未模拟细胞 GSH 耗竭。
Flanagan 等提供 GSH 本征反应性测量方法；这类测量并不直接给出临床安全结论。

## 可复现性、使用及限制

随机种子仅控制构象生成，动力学轨迹未加随机噪声。软件版本、文件 SHA-256
及验证结果均保存；检查包括浓度非负性、质量守恒、解析/数值一致性、速率极限、
蛋白恢复和图片 300 DPI。同一软件构建下应可复现，但不承诺跨平台逐字节一致。
图4是混合机理、独立演示速率与计算 LUMO 的探索性线性回归，无论相关性强弱均
不能证明因果或预测能力。xTB 与 EHT 均为近似方法，默认能垒和速率未作实验验证。

运行 `python run_task3_covalent_kinetics_residence_time.py --help` 查看参数。
可修改 `data_task3/parameters_template.json` 后通过 `--parameters` 导入带来源的
速率；仅覆盖部分字段时，其余仍为演示值，来源说明应写清这一点。kGSH 覆盖值
按当前实验条件的表观速率使用，不自动作 pH 修正。`--self-test-only` 执行核心
数学检查。`--git-sync` 要求已有 main 分支检出及 origin 远端，仅暂存生成文件，
正常提交推送，不强制推送或切换分支。重复运行幂等更新 README 区块并覆盖专用
输出路径；请勿同时向同一目录运行多个实例。
"""
    zh += f"\n## 本次结果\n\n- 电子结构方法：{method}。\n" \
          f"- 温度 {args.temperature:g} K；pH {args.ph:g}；GSH pKa 假设 {args.gsh_pka:g}；" \
          f"GSH {args.gsh_mM:g} mM。\n- Pareto 候选：{pareto}。\n- 高反应性标记：{flags}。\n" \
          f"- β 位局域化未获支持或贡献很弱的分子：{unsupported}。\n- 探索性 LFER R²：{lfer['R_squared']:.4f}。\n" \
          f"- ODE 拟合 kobs 相对特征值的最大误差：{maxerr:.3g}。\n\n{table}\n\n" \
          f"### 各化合物参数来源\n\n{provenance}\n\n## 参考资料\n\n{refs}\n"
    write_text(out / "COVALENT_DRUG_KINETICS_REPORT_EN.md", en)
    write_text(out / "COVALENT_DRUG_KINETICS_REPORT_ZH.md", zh)


def update_readme(out: Path) -> None:
    path = out / "README.md"
    existing = path.read_text(encoding="utf-8-sig") if path.exists() else "# Computational pharmacology\n"
    start, end = "<!-- TASK3:BEGIN -->", "<!-- TASK3:END -->"
    require(existing.count(start) == existing.count(end) and existing.count(start) <= 1,
            "README has malformed Task 3 markers; no overwrite performed")
    section = f"""{start}
## Task 3 — Covalent inhibitor kinetics and residence time

- [Standalone Python workflow](run_task3_covalent_kinetics_residence_time.py)
- [English technical report](COVALENT_DRUG_KINETICS_REPORT_EN.md)
- [中文技术报告](COVALENT_DRUG_KINETICS_REPORT_ZH.md)
- [Four publication-resolution figures](figures_task3/)
- [Results and parameter provenance](data_task3/warhead_results.csv)
- [Numerical validation](data_task3/validation.json)
- [Run manifest and file hashes](data_task3/run_manifest.json)

```sh
python -m pip install numpy scipy matplotlib rdkit pillow
python run_task3_covalent_kinetics_residence_time.py
python run_task3_covalent_kinetics_residence_time.py --self-test-only
# Optional external xTB installation:
python run_task3_covalent_kinetics_residence_time.py --quantum xtb --xtb xtb
# Optional commit/push from an existing main checkout with origin:
python run_task3_covalent_kinetics_residence_time.py --git-sync
```

Run these local commands from the directory containing this script.
Default outputs are written beside the script, including `figures_task3/`.
Use `--output-dir PATH` for a separate output root.
Orbital descriptors are calculated with EHT (or optional GFN2-xTB); kinetic inputs
and the barrier scenario are illustrative and uncalibrated. The 30-minute GSH flag
and shaded efficiency band are screening heuristics, not clinical/DILI predictions.
The reports distinguish KD, kinetic KI, fitted KI, chemical residence and turnover.
{end}"""
    if start in existing:
        first, rest = existing.split(start, 1)
        _, last = rest.split(end, 1)
        result = first + section + last
    else:
        result = existing.rstrip() + "\n\n" + section + "\n"
    write_text(path, result)


def self_tests() -> dict:
    w = default_panel()[0]
    slow, _ = slow_rate(w, 1e-15)
    expected = w.kinact_s * w.kon_M_inv_s / (w.koff_s + w.kinact_s)
    require(abs(slow/1e-15/expected - 1) < 1e-6, "Low-concentration kinetic limit failed")
    slow, _ = slow_rate(w, 1.)
    require(abs(slow/w.kinact_s - 1) < 1e-6, "Saturating kinetic limit failed")
    kgsh, gsh = .1, .005
    half = LN2 / (kgsh*gsh)
    require(abs(math.exp(-kgsh*gsh*half)-.5) < 1e-14, "GSH half-life formula failed")
    half_protein = 24 * 3600
    kdeg = LN2 / half_protein
    require(abs(math.exp(-kdeg*half_protein)-.5) < 1e-14, "Turnover half-life formula failed")
    require(abs(1/kdeg - half_protein/LN2) < 1e-8, "Mean lifetime / half-life distinction failed")
    # Full washout must match a single reversible state's known release kinetics.
    args = argparse.Namespace(pulse_uM=1., washout_hours=12.)
    ctrl = replace(w, kinact_s=0., krev_s=0., koff_s=.001)
    t, y, mean = washout_model(ctrl, 24., args)
    post = t >= 0
    initial = y[1, post][0]
    expected_curve = initial*np.exp(-(.001+kdeg)*t[post])
    require(np.max(abs(y[1, post]-expected_curve)) < 2e-6, "Reversible washout analytic limit failed")
    require(abs(mean*(.001+kdeg)-1) < 1e-8, "Bound-state first passage test failed")
    return {"status": "passed", "checks": ["low-dose QSSA slope", "saturating kinact",
        "GSH half-life", "protein half-life versus mean lifetime",
        "reversible washout analytic limit", "bound-state mean first passage"]}


def git_sync(out: Path, generated: list[Path]) -> dict:
    """Never stage unrelated pre-existing changes or force a remote update."""
    git_directory = out
    def git(*args: str, check: bool = True) -> subprocess.CompletedProcess:
        result = subprocess.run(["git", "-C", str(git_directory), *args], capture_output=True,
                                text=True, errors="replace", timeout=120, check=False)
        if check and result.returncode:
            raise RuntimeError(f"git {' '.join(args)} failed: {result.stderr.strip()}")
        return result
    require(shutil.which("git") is not None, "git executable not available")
    root_result = git("rev-parse", "--show-toplevel", check=False)
    require(root_result.returncode == 0,
            "Git sync unavailable: output directory is not inside an existing Git repository.")
    root = Path(root_result.stdout.strip()).resolve()
    # Stage repository-relative paths from the repository root, not the nested
    # project/output folder used only for repository discovery.
    git_directory = root
    require(git("branch", "--show-current").stdout.strip() == "main",
            "Git sync requires an existing main checkout; branch was not changed.")
    git("remote", "get-url", "origin")
    paths = [p.resolve().relative_to(root).as_posix() for p in generated]
    staged = git("diff", "--cached", "--name-only", "-z").stdout.split("\0")
    require(not (set(filter(None, staged)) - set(paths)),
            "Unrelated changes are already staged; commit left untouched.")
    # Stage in bounded chunks to respect Windows command-line length limits.
    for offset in range(0, len(paths), 30):
        git("add", "--", *paths[offset:offset+30])
    changes = git("diff", "--cached", "--quiet", check=False)
    require(changes.returncode in (0, 1), "Cannot inspect staged changes")
    committed = changes.returncode == 1
    if committed:
        git("commit", "-m", COMMIT_MESSAGE)
    git("push", "origin", "main")
    return {"status": "pushed", "created_commit": committed, "commit": git("rev-parse", "HEAD").stdout.strip()}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent,
                        help="Output directory (default: this project's directory)")
    parser.add_argument("--quantum", choices=("eht", "xtb"), default="eht")
    parser.add_argument("--xtb", default="xtb", help="xTB executable name or path")
    parser.add_argument("--xtb-timeout", type=int, default=600, help="Seconds per external xTB call")
    parser.add_argument("--conformers", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260910)
    parser.add_argument("--concentrations", type=int, default=25)
    parser.add_argument("--parameters", type=Path, help="Optional JSON parameter overrides with provenance")
    parser.add_argument("--temperature", type=float, default=310.15, help="Kelvin, barrier-to-rate scenario")
    parser.add_argument("--ph", type=float, default=7.4)
    parser.add_argument("--gsh-pka", type=float, default=8.7, help="Illustrative effective thiol pKa")
    parser.add_argument("--gsh-mM", dest="gsh_mM", type=float, default=5.)
    parser.add_argument("--target-nM", dest="target_nM", type=float, default=10.)
    parser.add_argument("--pulse-uM", dest="pulse_uM", type=float, default=1.)
    parser.add_argument("--drug-half-life-h", type=float, default=1.)
    parser.add_argument("--turnover-hours", type=float, nargs="+", default=[8., 24., 72.])
    parser.add_argument("--washout-hours", type=float, default=120.)
    parser.add_argument("--self-test-only", action="store_true")
    parser.add_argument("--git-sync", action="store_true", help="Commit generated files and push origin main")
    args = parser.parse_args(argv)
    for field in ("temperature", "gsh_mM", "target_nM", "pulse_uM", "drug_half_life_h", "washout_hours"):
        require(math.isfinite(getattr(args, field)) and getattr(args, field) > 0, f"{field} must be positive")
    require(200 <= args.temperature <= 400, "temperature must be 200–400 K")
    require(0 <= args.ph <= 14 and 0 <= args.gsh_pka <= 14, "pH and pKa must be in [0,14]")
    require(1 <= args.conformers <= 100, "conformers must be 1–100")
    require(8 <= args.concentrations <= 200, "concentrations must be 8–200")
    require(0 <= args.seed < 2**31-8, "seed outside RDKit supported integer range")
    require(args.xtb_timeout > 0, "xTB timeout must be positive")
    require(1 <= len(args.turnover_hours) <= 5 and len(set(args.turnover_hours)) == len(args.turnover_hours)
            and all(math.isfinite(v) and v > 0 for v in args.turnover_hours),
            "Provide 1–5 distinct positive turnover half-lives")
    if args.quantum == "xtb" and not args.self_test_only:
        binary = shutil.which(args.xtb)
        require(binary is not None, "xTB executable not found; install xTB or use --quantum eht")
        args.xtb = str(Path(binary).resolve())
    return args


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args(argv)
    validation = self_tests()
    if args.self_test_only:
        print(json.dumps(validation, indent=2))
        return 0
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    data = out / "data_task3"
    data.mkdir(exist_ok=True)
    panel = load_panel(args.parameters)
    # Copy the complete artifact into the selected output root for README portability.
    script = Path(__file__).resolve()
    destination = out / "run_task3_covalent_kinetics_residence_time.py"
    if script != destination:
        shutil.copy2(script, destination)
    LOG.info("Starting Task 3; default kinetics and barrier estimates are illustrative.")
    electronic, atoms = quantum_panel(panel, args, out)
    LOG.info("Integrating concentration series and fitting inactivation kinetics")
    samples, trajectories = concentration_kinetics(panel, electronic, args.concentrations)
    LOG.info("Integrating GSH conjugation, finite-drug exposure and washout")
    gsh_curves = reactivity(panel, electronic, args)
    dynamic = dynamic_exposure(panel, electronic, args)
    recovery, residence, reduced = washout(panel, args)
    tables = {"warhead_results.csv": electronic, "atomic_frontier_descriptors.csv": atoms,
        "kobs_concentration.csv": samples, "clamped_concentration_trajectories.csv": trajectories,
        "gsh_conjugation.csv": gsh_curves, "dynamic_drug_exposure.csv": dynamic,
        "washout_trajectories.csv": recovery, "residence_summary.csv": residence,
        "reduced_turnover_trajectories.csv": reduced}
    for name, rows in tables.items():
        write_csv(data / name, rows)
    lfer = make_figures(electronic, samples, recovery, args, out / "figures_task3")
    reports(electronic, residence, lfer, args, out)
    update_readme(out)
    template_fields = ("kon_M_inv_s", "koff_s", "kinact_s", "krev_s", "barrier_baseline_kJ_mol",
                       "kgsh_reverse_s", "kgsh_M_inv_s", "provenance")
    write_json(data / "parameters_template.json", {"warheads": {
        w.id: {key: getattr(w, key) for key in template_fields} for w in panel}})
    write_json(data / "lfer_regression.json", lfer)
    validation["checks"] += ["8 valid shared-core structures", "MMFF converged conformers",
        "finite frontier orbitals and atom populations", "clamped target mass conservation",
        "dynamic target and drug mass conservation", "GSH ODE / analytic agreement",
        "washout total target conservation", "reduced turnover analytic recovery", "four 300 DPI PNGs"]
    validation["max_kobs_eigen_relative_error"] = max(d["max_kobs_eigen_relative_error"] for d in electronic)
    require(validation["max_kobs_eigen_relative_error"] < 1e-4, "kobs/eigenvalue validation tolerance exceeded")
    validation["quantum_backend_exercised"] = args.quantum
    validation["xtb_backend_validated_this_run"] = args.quantum == "xtb"
    write_json(data / "validation.json", validation)
    # Enumerate only this workflow's owned artifacts, not unrelated files in output-dir.
    generated = [destination, out / "README.md", out / "COVALENT_DRUG_KINETICS_REPORT_EN.md",
                 out / "COVALENT_DRUG_KINETICS_REPORT_ZH.md"]
    generated += [data / name for name in tables]
    generated += [data / name for name in ("parameters_template.json", "lfer_regression.json", "validation.json")]
    generated += [data / "structures" / f"{w.id}.xyz" for w in panel]
    generated += [data / "structures" / "warhead_panel.sdf"]
    if args.quantum == "xtb":
        generated += [f for f in (data / "quantum").rglob("*") if f.is_file()]
    generated += [out / "figures_task3" / name for name in (
        "fig1_kobs_concentration_hyperbola.png", "fig2_target_recovery_washout.png",
        "fig3_gsh_reactivity_safety_radar.png", "fig4_warhead_lumo_reactivity_correlation.png")]
    manifest = {"created_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version, "versions": {"numpy": np.__version__, "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__, "rdkit": rdkit.__version__,
        "pillow": importlib.metadata.version("pillow")},
        "arguments": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        "panel_parameters": [asdict(w) for w in panel], "references": REFERENCES,
        "evidence_status": "Computed orbitals; illustrative default kinetics and uncalibrated barrier surrogate",
        "sha256": {f.relative_to(out).as_posix(): hashlib.sha256(f.read_bytes()).hexdigest()
                   for f in sorted(set(generated))}}
    write_json(data / "run_manifest.json", manifest)
    generated.append(data / "run_manifest.json")
    if args.git_sync:
        try:
            status = git_sync(out, generated)
        except (RuntimeError, OSError, subprocess.TimeoutExpired) as exc:
            write_json(data / "git_sync_status.json", {"status": "failed", "reason": str(exc),
                "note": "Scientific deliverables completed; a local commit may precede a failed push."})
            LOG.error("Deliverables are complete, but Git synchronization failed: %s", exc)
            return 3
        write_json(data / "git_sync_status.json", status)
    else:
        write_json(data / "git_sync_status.json", {"status": "not_requested",
            "note": "Use --git-sync from an existing main checkout with origin."})
    LOG.info("Completed: %s", out)
    LOG.info("Validation passed; 4 figures at 300 DPI, 2 reports, 9 CSV tables and provenance saved.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, ValueError, OSError, subprocess.TimeoutExpired, KeyError) as error:
        LOG.error("Task 3 failed: %s", error)
        raise SystemExit(1) from error
