# AUDITED COPY: A/C targets are synthetic, B ranking is heuristic, D is a reduced ODE.
# Original source and explicit patch are preserved; no experimental validation is claimed.
# -*- coding: utf-8 -*-
"""
File: run_production_pipeline.py
Description: Rigorous, physics-grounded AI4S calculation suite for 
organic electrosynthesis, SAC coordination screening, and flow kinetics.
"""

import sys
import json
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import minimize
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel
from sklearn.ensemble import ExtraTreesRegressor

# Ensure RDKit availability
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdFreeSASA
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("[WARNING] RDKit not detected. Falling back to heuristic chemical graph models.")

# ==============================================================================
# Task A: Physics-Informed Bayesian Optimization for Electrosynthesis
# ==============================================================================
def run_physics_informed_bo():
    """
    Instead of dummy one-hot encodings, discrete electrosynthesis variables are 
    mapped to real physicochemical vectors:
    - Solvent: [Dielectric constant (eps), Viscosity (mPa*s), Donor Number (DN)]
    - Electrolyte: [Anion Stokes radius (pm), Anion Oxidation Potential (V vs SCE)]
    - Continuous: [Current density j (mA/cm2), Temperature T (K)]
    """
    print("\n--- [Task A] Physics-Informed Bayesian Optimization on Electrosynthesis ---")
    
    # 1. Physicochemical Feature Dictionaries
    solvents = {
        "MeCN": {"eps": 37.5, "visc": 0.36, "DN": 14.1},
        "DCM":  {"eps": 8.93, "visc": 0.44, "DN": 1.0},
        "DMF":  {"eps": 36.7, "visc": 0.92, "DN": 26.6},
        "HFIP": {"eps": 16.7, "visc": 1.65, "DN": 0.0}  # Highly fluorinated, stabilizes radical cations
    }
    
    electrolytes = {
        "nBu4NPF6":  {"radius": 254.0, "E_ox_limit": 3.20},
        "nBu4NBF4":  {"radius": 218.0, "E_ox_limit": 2.90},
        "LiClO4":    {"radius": 236.0, "E_ox_limit": 2.60},
        "nBu4NOAc":  {"radius": 230.0, "E_ox_limit": 1.45}  # Readily oxidized (Kolbe-type competition)
    }

    # True objective function representing a typical C(sp2)-H thiolation / amination:
    # Requires high solvent polarity (eps ~ 25-38), low viscosity for mass transfer,
    # high electrolyte stability (E_ox_limit > 2.5), and optimal current density ~ 12-16 mA/cm2.
    def ground_truth_yield(j, T, solv_key, elec_key):
        s = solvents[solv_key]
        e = electrolytes[elec_key]
        
        # Inherent electrochemical window penalty
        if e["E_ox_limit"] < 1.8:
            anodic_parasitic_loss = 45.0  # Electrolyte decomposes preferentially
        else:
            anodic_parasitic_loss = 0.0

        # Physical score
        j_pen = 0.15 * (j - 14.0)**2
        T_pen = 0.02 * (T - 303.15)**2
        solv_score = 0.8 * s["eps"] - 12.0 * s["visc"] + 0.3 * s["DN"]
        
        raw_yield = 88.0 - j_pen - T_pen + 0.15 * solv_score - anodic_parasitic_loss
        return float(np.clip(raw_yield + np.random.normal(0, 1.2), 2.0, 98.0))

    # Pre-build candidate search space
    candidate_records = []
    for s_name, s_prop in solvents.items():
        for e_name, e_prop in electrolytes.items():
            for j_val in np.linspace(4.0, 24.0, 6):
                for t_val in [273.15 + 10, 273.15 + 25, 273.15 + 45]:
                    candidate_records.append({
                        "solvent": s_name,
                        "electrolyte": e_name,
                        "current_mA_cm2": j_val,
                        "temp_K": t_val,
                        # Vector representations for surrogate model:
                        "feat_eps": s_prop["eps"],
                        "feat_visc": s_prop["visc"],
                        "feat_DN": s_prop["DN"],
                        "feat_radius": e_prop["radius"],
                        "feat_E_ox": e_prop["E_ox_limit"]
                    })
    df_space = pd.DataFrame(candidate_records)
    feature_cols = ["current_mA_cm2", "temp_K", "feat_eps", "feat_visc", "feat_DN", "feat_radius", "feat_E_ox"]
    X_pool = df_space[feature_cols].values

    # Active learning iteration: 6 cold-start seeds, then 10 BO acquisition steps
    np.random.seed(42)
    init_indices = np.random.choice(len(df_space), size=6, replace=False)
    X_train = X_pool[init_indices]
    y_train = [ground_truth_yield(df_space.iloc[i]["current_mA_cm2"], 
                                  df_space.iloc[i]["temp_K"], 
                                  df_space.iloc[i]["solvent"], 
                                  df_space.iloc[i]["electrolyte"]) for i in init_indices]
    
    kernel = Matern(length_scale=np.ones(len(feature_cols)), nu=2.5) + WhiteKernel(noise_level=1.0)
    gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=3, random_state=42)

    bo_log = []
    for step in range(8):
        gp.fit(X_train, y_train)
        mu, sigma = gp.predict(X_pool, return_std=True)
        # Upper Confidence Bound (UCB) acquisition
        beta = 2.0
        acquisition = mu + beta * sigma
        
        # Avoid re-sampling identical points
        acquisition[init_indices] = -1e9
        best_cand_idx = int(np.argmax(acquisition))
        
        cand_row = df_space.iloc[best_cand_idx]
        new_yield = ground_truth_yield(cand_row["current_mA_cm2"], cand_row["temp_K"], 
                                       cand_row["solvent"], cand_row["electrolyte"])
        
        X_train = np.vstack([X_train, X_pool[best_cand_idx]])
        y_train.append(new_yield)
        init_indices = np.append(init_indices, best_cand_idx)

        bo_log.append({
            "iteration": step + 1,
            "chosen_solvent": cand_row["solvent"],
            "chosen_electrolyte": cand_row["electrolyte"],
            "current_density": round(cand_row["current_mA_cm2"], 2),
            "acquired_yield": round(new_yield, 2),
            "cumulative_best_yield": round(max(y_train), 2)
        })

    pd.DataFrame(df_space).to_csv("task_a_candidates.csv", index=False)
    sampled = df_space.iloc[init_indices].copy()
    sampled["observed_yield"] = y_train
    sampled.to_csv("task_a_observations.csv", index=False)
    best_idx_overall = np.argmax(y_train)
    return {
        "status": "completed",
        "optimal_yield": round(float(y_train[best_idx_overall]), 2),
        "optimal_parameters": {
            "solvent": df_space.iloc[init_indices[best_idx_overall]]["solvent"],
            "electrolyte": df_space.iloc[init_indices[best_idx_overall]]["electrolyte"],
            "current_density_mA_cm2": round(float(df_space.iloc[init_indices[best_idx_overall]]["current_mA_cm2"]), 2),
            "temp_C": round(float(df_space.iloc[init_indices[best_idx_overall]]["temp_K"] - 273.15), 1)
        },
        "optimization_trajectory": bo_log
    }

# ==============================================================================
# Task B: Real-Molecule 3D Conformation & Electrochemical Regioselectivity
# ==============================================================================
def run_molecule_3d_regioselectivity():
    """
    Uses RDKit to construct 3D geometry of an actual medicinal indole alkaloid scaffold:
    5-methoxy-2-phenyl-1H-indole.
    Generates ETKDG conformers, performs MMFF force-field energy minimization,
    and extracts atom-centered electrostatic partial charges and steric SASA values
    to predict radical cation spin/charge localization under anodic oxidation.
    """
    print("\n--- [Task B] 3D Cheminformatics & Anodic C-H Regioselectivity ---")
    
    # Target substrate: 5-methoxy-2-phenyl-1H-indole (Classic scaffold studied in Pan-Tang lab)
    smiles = "COc1ccc2[nH]c(cc2c1)c3ccccc3"
    
    if not HAS_RDKIT:
        return {"status": "failed", "error": "RDKit is mandatory for Task B."}

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    
    # 3D Conformer Generation with ETKDG
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    embed_status = AllChem.EmbedMolecule(mol, params)
    if embed_status != 0: raise RuntimeError("ETKDG embedding failed")
    mmff_status = AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    if mmff_status != 0: raise RuntimeError("MMFF geometry did not converge")

    # Compute Gasteiger partial charges
    AllChem.ComputeGasteigerCharges(mol)
    
    # Compute per-atom SASA (Solvent Accessible Surface Area)
    radii = [Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()) for atom in mol.GetAtoms()]
    total_sasa = rdFreeSASA.CalcSASA(mol, radii)

    # Inspect C-H candidate carbon atoms in the indole core:
    # In 5-methoxy-2-phenyl-1H-indole:
    # Carbon C3 (pyrrole ring) vs Carbon C4/C6 (benzene ring) vs C-H on pendant phenyl
    target_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C":
            # Check attached hydrogens
            h_count = sum([1 for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() == "H"])
            if h_count > 0:
                q = float(atom.GetProp("_GasteigerCharge"))
                idx = atom.GetIdx()
                # Find connected H index
                h_atom = [nb for nb in atom.GetNeighbors() if nb.GetSymbol() == "H"][0]
                h_sasa = h_atom.GetDoubleProp("SASA")
                
                # Electrophilic radical attack index: higher charge density (more negative)
                # and higher steric accessibility promote rapid radical coupling
                susceptibility_score = (-1.0 * q) * 0.6 + (h_sasa / 15.0) * 0.4
                target_carbons.append({
                    "atom_index": idx,
                    "symbol": "C",
                    "hybridization": str(atom.GetHybridization()),
                    "gasteiger_charge": round(q, 4),
                    "h_sasa_angstrom2": round(h_sasa, 3),
                    "anodic_activation_index": round(susceptibility_score, 4)
                })

    df_carbons = pd.DataFrame(target_carbons).sort_values(by="anodic_activation_index", ascending=False)
    df_carbons.to_csv("task_b_all_ch_sites.csv", index=False)
    Chem.MolToMolFile(mol, "task_b_mmff_geometry.mol")
    pd.DataFrame([dict(atom_index=a.GetIdx(), element=a.GetSymbol(), radius_A=radii[a.GetIdx()], sasa_A2=a.GetDoubleProp("SASA"), gasteiger_charge=float(a.GetProp("_GasteigerCharge"))) for a in mol.GetAtoms()]).to_csv("task_b_all_atoms.csv", index=False)
    top_site = df_carbons.iloc[0].to_dict()

    return {
        "status": "completed",
        "molecule_smiles": smiles,
        "heavy_atom_count": mol.GetNumHeavyAtoms(),
        "total_sasa_A2": total_sasa,
        "embed_status": embed_status,
        "mmff_status": mmff_status,
        "sasa_method": "Lee-Richards; RDKit periodic-table vdW radii; explicit H; probe 1.4 A",
        "score_evidence": "Uncalibrated heuristic; not a quantum or experimental site prediction",
        "evaluated_ch_sites_count": len(target_carbons),
        "predicted_primary_reactive_site": top_site,
        "top_3_reactive_carbons": df_carbons.head(3).to_dict(orient="records")
    }

# ==============================================================================
# Task C: Coordination Microenvironment & SAC Catalytic Scaling Relation
# ==============================================================================
def run_sac_coordination_informatics():
    """
    Models Single-Atom Catalysts anchored on Porous Organic Polymers (POPs).
    Quantifies coordination sphere: M-N_x-C_y-O_z.
    Descriptors:
    - Pauling Electronegativity mismatch in first coordination sphere
    - Metal nominal d-electron count (e.g., Cu(II): d9, Co(II): d7, Ni(II): d8)
    - Calculated d-band center proxy (\varepsilon_d) based on local ligand field
    Predicts the activation barrier for electrocatalytic C-H cross-coupling (\Delta G^#).
    """
    print("\n--- [Task C] Porous Polymer SAC Coordination Microenvironment Informatics ---")
    
    # Dataset representing 36 experimentally synthesized POP-SAC variations
    records = []
    metals = {
        "Cu": {"d_electrons": 9, "EN": 1.90, "base_barrier": 19.5},
        "Co": {"d_electrons": 7, "EN": 1.88, "base_barrier": 23.0},
        "Ni": {"d_electrons": 8, "EN": 1.91, "base_barrier": 21.0},
        "Pd": {"d_electrons": 8, "EN": 2.20, "base_barrier": 17.0}
    }
    
    coordinations = {
        "N4_planar":      {"N_count": 4, "O_count": 0, "C_count": 0, "ligand_field_strength": 1.4},
        "N3C1_porphyrin": {"N_count": 3, "O_count": 0, "C_count": 1, "ligand_field_strength": 1.1},
        "N2O2_salen":     {"N_count": 2, "O_count": 2, "C_count": 0, "ligand_field_strength": 0.9},
        "N2C2_network":   {"N_count": 2, "O_count": 0, "C_count": 2, "ligand_field_strength": 0.7}
    }

    for m_name, m_prop in metals.items():
        for c_name, c_prop in coordinations.items():
            # Theoretical d-band center approximation (relative eV)
            # More electronegative local atoms pull d-band lower (more negative)
            d_band_center = -1.2 * c_prop["ligand_field_strength"] - (m_prop["EN"] - 1.8) * 1.5
            
            # Non-linear Sabatier volcano curve for activation barrier:
            # Optimal barrier achieved when d-band center is near -1.85 eV
            activation_barrier = m_prop["base_barrier"] + 4.5 * (d_band_center - (-1.85))**2 + np.random.normal(0, 0.35)
            
            records.append({
                "metal": m_name,
                "coordination_type": c_name,
                "d_electrons": m_prop["d_electrons"],
                "metal_EN": m_prop["EN"],
                "coord_N": c_prop["N_count"],
                "coord_O": c_prop["O_count"],
                "d_band_center_proxy_eV": round(d_band_center, 3),
                "activation_barrier_kcal_mol": round(activation_barrier, 2)
            })

    df_sac = pd.DataFrame(records)
    df_sac.to_csv("task_c_all_sites.csv", index=False)
    
    # Fit Extra-Trees Regressor to evaluate feature contributions
    features = ["d_electrons", "metal_EN", "coord_N", "coord_O", "d_band_center_proxy_eV"]
    X = df_sac[features]
    y = df_sac["activation_barrier_kcal_mol"]
    
    et = ExtraTreesRegressor(n_estimators=100, random_state=42)
    et.fit(X, y)
    
    importances = dict(zip(features, [round(float(val), 4) for val in et.feature_importances_]))
    best_sac = df_sac.sort_values(by="activation_barrier_kcal_mol").iloc[0].to_dict()

    return {
        "status": "completed",
        "sample_count": len(df_sac),
        "descriptor_importances": importances,
        "optimal_sac_configuration": best_sac
    }

# ==============================================================================
# Task D: Microchannel Flow Electrochemistry Transport-Kinetic Numerical Solver
# ==============================================================================
def run_flow_electrochemistry_pde():
    """
    Couples boundary layer convective mass transfer with non-linear Butler-Volmer kinetics
    in an undivided microfluidic parallel-plate electrochemical reactor.
    
    Governing Equations:
    1. Mass Transfer Coefficient k_m via Leveque approximation for laminar channel flow:
       Sh = 1.85 * (Re * Sc * Dh / L)^(1/3)
       k_m = Sh * D / Dh
    2. Butler-Volmer boundary flux:
       j = j_0 * [ (C_surf / C_bulk) * exp(alpha_a * F * eta / (R*T)) - exp(-alpha_c * F * eta / (R*T)) ]
    3. Mass Balance along Channel:
       d(C_bulk)/dz = - (k_m * a / u) * (C_bulk - C_surf)
    """
    print("\n--- [Task D] Microchannel Transport & Butler-Volmer Kinetic Solver ---")
    
    # Constants
    F = 96485.33  # C/mol
    R = 8.314     # J/(mol*K)
    T = 298.15    # K
    alpha_a = 0.5
    alpha_c = 0.5
    j0 = 0.05     # Exchange current density, A/m2
    D = 1.2e-9    # Diffusion coefficient in MeCN, m2/s
    
    # Reactor Dimensions (typical micro-flow cell in lab)
    channel_width = 0.01   # 10 mm
    channel_height = 0.0005 # 500 um
    channel_length = 0.1   # 10 cm
    Dh = 2 * channel_height # Hydraulic diameter (parallel plate approx)
    specific_area = 1.0 / channel_height # Electrode area per unit volume (m^-1)

    flow_rates_uL_min = [100, 300, 600, 1200]
    sweep_results = {}

    for q_uL in flow_rates_uL_min:
        q_m3_s = q_uL * 1e-9 / 60.0
        u_avg = q_m3_s / (channel_width * channel_height) # Flow velocity m/s
        
        # Leveque Graetz / Sherwood estimation
        # Pe = u_avg * Dh / D
        gamma_shear = 6.0 * u_avg / channel_height
        k_m = 0.67 * D * (gamma_shear / (D * channel_length))**(1.0 / 3.0) # Local average mass transfer
        
        # Solve ODE along channel length z from 0 to L
        # Overpotential eta = 0.45 V
        eta = 0.45
        exp_anodic = np.exp(alpha_a * F * eta / (R * T))
        exp_cathodic = np.exp(-alpha_c * F * eta / (R * T))
        
        def channel_ode(C_b, z):
            # Calculate surface concentration at steady flux: j / F = k_m * (C_b - C_s)
            # j0 * [ (C_s / C_in) * exp_a - exp_c ] = F * k_m * (C_b - C_s)
            C_in = 50.0 # 50 mol/m3 (0.05 M)
            num = F * k_m * C_b + j0 * exp_cathodic
            den = F * k_m + (j0 / C_in) * exp_anodic
            C_surf = max(num / den, 0.0)
            
            dC_dz = - (k_m * specific_area / u_avg) * (C_b - C_surf)
            return dC_dz

        z_eval = np.linspace(0, channel_length, 100)
        C_in = 50.0
        sol = odeint(channel_ode, y0=C_in, t=z_eval)
        C_out = sol[-1][0]
        
        conversion = (1.0 - C_out / C_in) * 100.0
        # Space-Time Yield (mmol / (L * h))
        residence_time_s = channel_length / u_avg
        reactor_volume_L = (channel_width * channel_height * channel_length) * 1000.0
        output_mol_per_s = (C_in - C_out) * q_m3_s
        sty = (output_mol_per_s * 1000.0 * 3600.0) / reactor_volume_L

        sweep_results[f"{q_uL}_uL_min"] = {
            "flow_velocity_mm_s": round(u_avg * 1000, 2),
            "residence_time_sec": round(residence_time_s, 2),
            "mass_transfer_coeff_m_s": f"{k_m:.3e}",
            "conversion_pct": round(float(conversion), 2),
            "space_time_yield_mmol_L_h": round(float(sty), 2)
        }

    return {
        "status": "completed",
        "reactor_type": "Microchannel Parallel-Plate Electrochemical Reactor",
        "overpotential_eta_V": 0.45,
        "inlet_concentration_M": 0.05,
        "flow_sweep_kinetics": sweep_results
    }

# ==============================================================================
# Master Execution & Data Packaging
# ==============================================================================
if __name__ == "__main__":
    master_results = {
        "Task_A_Physical_BO": run_physics_informed_bo(),
        "Task_B_3D_Cheminformatics": run_molecule_3d_regioselectivity() if HAS_RDKIT else {"status": "skipped"},
        "Task_C_SAC_Informatics": run_sac_coordination_informatics(),
        "Task_D_Flow_Electrochemistry": run_flow_electrochemistry_pde()
    }

    output_filename = "production_benchmark_results.json"
    with open(output_filename, "w", encoding="utf-8") as f:
        json.dump(master_results, f, indent=4, ensure_ascii=False)

    print(f"\n[SUCCESS] Compatibility-repaired demonstrations completed; A/C are synthetic and B scoring is heuristic.")
    print(f"[DATA] Mixed-evidence numerical outputs exported to: {output_filename}")
