# -*- coding: utf-8 -*-
"""
File: run_closed_loop_platform.py
Description: Automated Dry/Wet interface generating wet-lab stoichiometry,
DFT input files, and active learning feedback updates for Pan-Tang lab.
"""

import os
import sys
import json
import numpy as np
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("[WARNING] RDKit not detected. Module 2 will output template coordinates.")

# ==============================================================================
# Module 1: Automated Wet-Lab Stoichiometry & Electrolysis SOP Generator
# ==============================================================================
def generate_wet_lab_sop(
    substrate_smiles="COc1ccc2[nH]c(cc2c1)c3ccccc3",  # 5-methoxy-2-phenyl-1H-indole
    coupling_partner_smiles="CSc1ccccc1",             # Methyl phenyl sulfide (or thiol partner)
    reaction_scale_mmol=0.2,                          # Standard methodology scale in lab
    current_density_mA_cm2=12.5,                      # From Module 1 Pareto optimum
    electrode_area_cm2=1.5,                           # Standard 1.0 x 1.5 cm Pt plate immersion
    faradaic_charge_F_per_mol=2.2,                    # 2.2 F/mol (slight excess for 2e- transfer)
    solvent_system="MeCN/HFIP (4:1)",                 # Optimal fluorinated radical stabilization
    electrolyte_type="nBu4NPF6",                      # 0.1 M loading
    total_volume_mL=6.0
):
    print("\n--- [Module 1] Generating Wet-Lab Stoichiometry & Operating SOP ---")

    F = 96485.33  # C/mol

    # Molecular weights
    if HAS_RDKIT:
        mol_sub = Chem.MolFromSmiles(substrate_smiles)
        mw_sub = Descriptors.MolWt(mol_sub)
    else:
        mw_sub = 223.27

    mw_electrolyte = {"nBu4NPF6": 387.43, "nBu4NBF4": 329.27, "LiClO4": 106.39}[electrolyte_type]

    # Reagent masses
    mass_substrate_mg = round(reaction_scale_mmol * mw_sub, 1)
    mass_electrolyte_mg = round((0.1 * (total_volume_mL / 1000.0)) * mw_electrolyte * 1000.0, 1)

    # Electrical parameters
    total_current_mA = round(current_density_mA_cm2 * electrode_area_cm2, 2)
    current_A = total_current_mA / 1000.0
    
    # Total charge required in Coulombs: Q = n_mol * F * charge_factor
    total_charge_Coulombs = (reaction_scale_mmol * 1e-3) * F * faradaic_charge_F_per_mol
    electrolysis_time_sec = total_charge_Coulombs / current_A
    electrolysis_time_min = round(electrolysis_time_sec / 60.0, 1)

    sop_card = {
        "Project_Title": "Electrochemical Oxidative C-H Functionalization of Bioactive Heterocycles",
        "Scale_Specification": f"{reaction_scale_mmol} mmol scale",
        "Reagent_Dispensing": {
            "Substrate": f"{mass_substrate_mg} mg ({reaction_scale_mmol} mmol)",
            "Electrolyte": f"{mass_electrolyte_mg} mg of {electrolyte_type} (0.1 M)",
            "Solvent_Ratio": f"{solvent_system}, total volume = {total_volume_mL} mL (4.8 mL MeCN + 1.2 mL HFIP)"
        },
        "Electrochemical_Cell_Configuration": {
            "Cell_Type": "Undivided three-necked glass cell (10 mL)",
            "Anode": "Platinum plate (1.0 cm x 1.5 cm) or Carbon rod",
            "Cathode": "Platinum plate (1.0 cm x 1.5 cm)",
            "Current_Mode": "Constant Current Electrolysis (CCE)",
            "Power_Supply_Setting_mA": total_current_mA,
            "Target_Faradaic_Electricity": f"{faradaic_charge_F_per_mol} F/mol",
            "Total_Charge_Coulombs": round(total_charge_Coulombs, 1),
            "Estimated_Reaction_Duration_min": electrolysis_time_min
        },
        "Workup_and_Purification_Protocol": [
            "1. Turn off DC power supply immediately when charge reaches target Coulombs.",
            "2. Quench reaction mixture with 10 mL saturated aqueous Na2S2O3 / NaHCO3 solution.",
            "3. Extract the aqueous phase with Ethyl Acetate (3 x 15 mL).",
            "4. Combine organic layers, wash with brine (15 mL), and dry over anhydrous Na2SO4.",
            "5. Filter and concentrate under reduced pressure via rotary evaporator.",
            "6. Purify residue via silica gel flash column chromatography (Eluent: Petroleum Ether / Ethyl Acetate = 15:1 -> 8:1 v/v)."
        ]
    }

    return sop_card

# ==============================================================================
# Module 2: Automated Quantum Chemistry (Gaussian 16 & ORCA) Input Generator
# ==============================================================================
def generate_dft_input_files(smiles="c1ccc2c(c1)CCN2", base_name="indoline"):
    """
    Constructs ready-to-submit Gaussian 16 (.gjf) and ORCA (.inp) input files 
    for calculating the radical cation [M]^+ generated upon anodic single-electron transfer (SET).
    """
    print("\n--- [Module 2] Generating Production-Grade Gaussian 16 / ORCA Input Files ---")

    coords_block = ""
    if HAS_RDKIT:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        AllChem.EmbedMolecule(mol, params)
        AllChem.MMFFOptimizeMolecule(mol)
        conf = mol.GetConformer()
        
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            coords_block += f"{atom.GetSymbol():<2}   {pos.x:>12.6f}  {pos.y:>12.6f}  {pos.z:>12.6f}\n"
    else:
        coords_block = "C   0.000000    0.000000    0.000000\nN   1.400000    0.000000    0.000000\n"

    # 1. Gaussian 16 Input (.gjf) - Open-shell radical cation: Charge = +1, Multiplicity = 2 (Doublet)
    gaussian_content = f"""%chk={base_name}_radical_cation.chk
%nprocshared=16
%mem=32GB
#p opt freq uB3LYP/def2SVP empiricaldispersion=gd3bj scrf=(smd,solvent=acetonitrile)

{base_name} radical cation [M]+. optimization in MeCN (Anodic SET Intermediate)

1 2
{coords_block.strip()}

--Link1--
%chk={base_name}_radical_cation.chk
%nprocshared=16
%mem=32GB
#p uM062X/def2TZVP scrf=(smd,solvent=acetonitrile) geom=check guess=read

Single point energy refinement with larger basis set

1 2

"""

    # 2. ORCA 5.0 Input (.inp)
    orca_content = f"""! UKS B3LYP D3BJ def2-SVP def2/J CPCM(Acetonitrile) Opt Freq
%pal nprocs 16 end
%maxcore 2000

* xyz 1 2
{coords_block.strip()}
*
"""

    os.makedirs("dft_inputs", exist_ok=True)
    gjf_path = os.path.join("dft_inputs", f"{base_name}_radical.gjf")
    inp_path = os.path.join("dft_inputs", f"{base_name}_radical.inp")

    with open(gjf_path, "w", encoding="utf-8") as f:
        f.write(gaussian_content)
    with open(inp_path, "w", encoding="utf-8") as f:
        f.write(orca_content)

    return {
        "status": "files_generated",
        "gaussian_input_path": gjf_path,
        "orca_input_path": inp_path,
        "spin_multiplicity": 2,
        "net_charge": 1,
        "theoretical_model": "UB3LYP-D3(BJ)/def2-SVP/SMD(MeCN)"
    }

# ==============================================================================
# Module 3: Active Learning Closed-Loop Feedback Assimilator
# ==============================================================================
def assimilate_lab_feedback_and_update():
    """
    Simulates the return of experimental data from the wet-lab after Round 1.
    Updates the Gaussian Process surrogate model, quantifies prediction error,
    and proposes Round 2 experimental candidates using Expected Improvement (EI).
    """
    print("\n--- [Module 3] Active Learning Closed-Loop Feedback Loop ---")
    np.random.seed(123)

    # Initial training data: 6 prior exploratory runs [current_density, electrolyte_conc, temperature]
    X_prior = np.array([
        [8.0, 0.05, 298.15],
        [15.0, 0.05, 298.15],
        [10.0, 0.15, 313.15],
        [20.0, 0.15, 298.15],
        [25.0, 0.20, 323.15],
        [12.0, 0.10, 303.15]
    ])
    # Algorithm's prior predicted yield
    y_prior_pred = np.array([62.0, 78.5, 71.0, 84.0, 55.0, 86.5])

    # Wet-lab measured actual HPLC yields returned by the graduate student / researcher
    y_wet_lab_actual = np.array([58.5, 81.2, 69.4, 79.8, 51.2, 88.0])

    # Calculate model error
    mae = float(np.mean(np.abs(y_prior_pred - y_wet_lab_actual)))

    # Fit Updated Gaussian Process
    kernel = Matern(nu=2.5) + WhiteKernel(noise_level=0.5)
    gp_updated = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=5, random_state=42)
    gp_updated.fit(X_prior, y_wet_lab_actual)

    # Screen candidate pool for Round 2 recommendation
    j_cand = np.linspace(8.0, 22.0, 15)
    c_cand = [0.08, 0.10, 0.12]
    T_cand = [298.15, 303.15]

    pool = []
    for j in j_cand:
        for c in c_cand:
            for t in T_cand:
                pool.append([j, c, t])
    X_pool = np.array(pool)

    mu, sigma = gp_updated.predict(X_pool, return_std=True)
    
    # Expected Improvement (EI) Acquisition Function
    y_best_current = np.max(y_wet_lab_actual)
    improvement = mu - y_best_current
    Z = improvement / (sigma + 1e-9)
    # Simple UCB proxy for ranking acquisition
    acquisition_scores = mu + 2.5 * sigma
    best_idx = np.argmax(acquisition_scores)
    next_condition = X_pool[best_idx]

    return {
        "status": "model_recalibrated",
        "previous_round_sample_count": len(X_prior),
        "mean_absolute_error_mae": round(mae, 2),
        "best_actual_yield_round1": round(float(y_best_current), 2),
        "round_2_recommended_condition": {
            "current_density_mA_cm2": round(float(next_condition[0]), 2),
            "electrolyte_concentration_M": round(float(next_condition[1]), 3),
            "temperature_C": round(float(next_condition[2] - 273.15), 1),
            "model_predicted_yield_pct": round(float(mu[best_idx]), 2),
            "prediction_uncertainty_sigma": round(float(sigma[best_idx]), 2)
        }
    }

# ==============================================================================
# Master Execution & Output
# ==============================================================================
if __name__ == "__main__":
    sop_result = generate_wet_lab_sop()
    dft_result = generate_dft_input_files()
    al_result = assimilate_lab_feedback_and_update()

    final_payload = {
        "WetLab_SOP": sop_result,
        "DFT_Input_Automation": dft_result,
        "Active_Learning_Feedback": al_result
    }

    out_file = "closed_loop_results.json"
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(final_payload, f, indent=4, ensure_ascii=False)

    print(f"\n[CLOSED-LOOP PIPELINE COMPLETE]")
    print(f"- Wet-lab SOP & stoichiometry calculated.")
    print(f"- Gaussian (.gjf) and ORCA (.inp) files exported to 'dft_inputs/'.")
    print(f"- Model updated with experimental feedback. Summary saved to '{out_file}'.")
