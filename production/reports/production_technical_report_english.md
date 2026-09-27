# Quantitative Organic Electrochemistry & Catalysis: Audited Technical Report

**Date:** 2026-09-27. **Audience:** Pan–Tang research group and undergraduate researchers, Guangxi Normal University. **Specification:** Four-section architecture in Section 3 of the newly supplied pipeline document. **Edition:** English; a separate Chinese edition and a section-aligned bilingual edition accompany it.

**Evidence statement.** The supplied script was attempted unchanged and failed in Task B before writing its final JSON. A separately archived repair completed all four tasks; exact reported source values refer to `results/repaired/production_benchmark_results.json`. A and C still use synthetic target equations. B computes actual molecular force-field geometry and geometric descriptors, but its site ranking is a heuristic. D is a reduced one-dimensional transport/kinetic ODE with assumed parameters. Additional studies comprise 30-seed BO comparisons, 32 conformers, grouped SAC holdouts, analytical/flux checks, and three actual GFN2-xTB jobs on the target molecule. No wet experiment, DFT transition state, calibrated oxidation potential, validated regioselectivity, or production readiness is established.

**Labels:** [L] literature; [S] executed synthetic calculation; [M] molecular mechanics/geometric calculation; [Q] executed semiempirical electronic structure; [A] code or independent numerical audit; [P] proposed experimental work; [U] unverified or missing evidence. A molecular input and physical units do not, by themselves, validate a predictive model.

## 1. Executive Scientific Context

### 1.1 Scientific alignment and current maturity

The proposed research connects three complementary questions: how to use electrical input to control organic bond formation, how medicinal heterocycle structure governs accessible redox pathways, and how a porous coordination environment changes catalytic function. Published work on porous-ligand Pd sites and single-atom Fe redox mediation supports this scientific intersection. These precedents motivate a laboratory-aligned program; they do not demonstrate that the present substrate is a standard group substrate or that every named instrument is available locally. [Porous-ligand Pd study](https://doi.org/10.1016/j.chempr.2020.06.020); [Fe-mediated electrosynthesis](https://onlinelibrary.wiley.com/doi/abs/10.1002/anie.202404295).

The revised workflow is more chemically explicit than the preceding generic benchmark: it distinguishes solvent/electrolyte identities, builds a defined molecular structure, introduces a coordination-dependent proxy, and balances electrode flux against convection. Its present maturity is nevertheless **an auditable computational prototype**. “Production-grade,” “without mock primitives,” and “high-fidelity” are claims in the supplied document, not demonstrated conclusions of the run.

| Task | Exact main output from the repaired execution | Supported interpretation |
|---|---|---|
| A | 92.03%; MeCN/LiClO4; 16.0 mA/cm²; 25.0 °C | Largest observed synthetic response among 14 evaluations |
| B | 17 heavy atoms; 10 C-bearing-H candidates; atom 15 score 0.7483 | Neutral geometric/charge heuristic favors a pendant-phenyl atom |
| C | 16 configurations; Pd/N2O2_salen score 17.09 kcal/mol | Lowest assigned synthetic score, not a computed activation free energy |
| D | 77.55→23.80% conversion; 465.30→1713.41 mmol/L/h STY | Reduced-model reactant-consumption trade-off as flow increases |

Three findings change the research interpretation. Physical-feature GP variants did not outperform random search in the added finite-budget comparison. C2 and C5 of the specified indole are substituted and have no C–H bonds. The flow model's mass-transfer coefficient increases with flow, while residence time falls; lower high-flow conversion cannot be attributed simply to worse mass transfer.

### 1.2 Execution, repairs, and provenance

[A] The existing environment was reused: Python 3.12.14, NumPy 2.4.6, SciPy 1.18.0, pandas 2.3.3, scikit-learn 1.9.0, RDKit 2026.03.5, Matplotlib 3.11.1 and xTB 6.7.1. Runs were CPU-only and limited to one thread. The untouched source exited with `AttributeError`: RDKit exposes `CalcSASA`, not `calcSASA`. The later `whichAtoms` keyword is also incompatible with the installed signature. Furthermore, `classifyAtoms` returned zero radii for all 30 atoms of this small molecule, so merely changing capitalization would not yield an appropriate van der Waals surface.

The repair uses `CalcSASA`, explicit RDKit periodic-table van der Waals radii, and atom-level `SASA` properties after the total calculation. It checks ETKDG/MMFF convergence, saves candidate/atom/configuration records, and replaces misleading terminal success text with an evidence-qualified message. A/C target equations, BO budget, source site-score coefficients, and D equations remain unchanged. The radius choice is a scientific convention change as well as a software repair; repaired B values must not be described as untouched-source results. The complete edit is recorded in `source/repair.patch`. [RDKit SASA API](https://www.rdkit.org/docs/source/rdkit.Chem.rdFreeSASA.html).

```text
Extracted original script SHA-256:
2b225a9ab8a66275ac276b11335b6d225660dfc8157d8f567f8ea7218a3349b6
```

The original failure, repaired logs, metadata, unrounded intermediate tables and additional audits are preserved separately. Source comments are retained as historical input and are explicitly evaluated here; they are not promoted into evidence. The original run produced no complete `production_benchmark_results.json`.

### 1.3 A falsifiable route to useful AI4S

The immediate research question is whether a traceable model improves a predeclared laboratory decision over a simple baseline on a chemically identified system. BO must be compared under equal experimental budgets; an oxidation predictor requires measured, reference-consistent labels; catalyst selection requires defined structures and comparable activity data; flow optimization requires measured selectivity and charge balance. Published experimental BO provides a precedent for this evaluation approach, not a guarantee that a particular embedding or GP will win. [Experimental Bayesian reaction optimization](https://www.nature.com/articles/s41586-021-03213-y).

Use four gates: G0, valid identity and executable records; G1, calibrated measurement and reproducibility; G2, frozen prospective predictions compared with baselines; G3, reproducible mechanistic or process value. This delivery improves G0 and performs computational checks. It does not pass G1–G3 for a new chemical reaction.

## 2. Rigorous Quantitative Breakdown

### 2.1 Task A — physical-property embeddings and mixed-variable BO

[S] The input space consists of four solvents, four electrolyte labels, six current densities (4, 8, 12, 16, 20, 24 mA/cm²), and three temperatures (10, 25, 45 °C): **288 candidates**. Although current and temperature are physically continuous, this implementation searches a finite grid. There are **6 initial evaluations plus 8 acquisition steps**, not the ten steps stated in a comment. Seven descriptor columns enter a Matérn-plus-WhiteKernel GP without feature scaling or target normalization; three optimizer restarts and UCB coefficient 2.0 are used.

The following values are **source-assigned descriptors**, retained exactly for reproduction. Their temperature dependence, measurement method, concentration, solvent/reference-electrode conventions, and supporting citations were not supplied. In particular, the electrolyte radii and oxidation limits are not verified universal material constants.

| Solvent | Assigned ε | Viscosity (mPa·s) | Donor number |
| --- | --- | --- | --- |
| MeCN | 37.5 | 0.36 | 14.1 |
| DCM | 8.93 | 0.44 | 1.0 |
| DMF | 36.7 | 0.92 | 26.6 |
| HFIP | 16.7 | 1.65 | 0.0 |

| Electrolyte | Assigned radius (pm) | Assigned oxidation limit (V vs SCE) |
| --- | --- | --- |
| nBu4NPF6 | 254.0 | 3.2 |
| nBu4NBF4 | 218.0 | 2.9 |
| LiClO4 | 236.0 | 2.6 |
| nBu4NOAc | 230.0 | 1.45 |

The response generator is

$$
y=\mathrm{clip}\{88-0.15(j-14)^2-0.02(T-303.15)^2+0.15[0.8\epsilon-12\mu+0.3DN]-45\mathbf{1}(E_{lim}<1.8)+\xi,2,98\},\quad \xi\sim N(0,1.2^2).
$$

It is an assigned scoring surface. No measured yield, overpotential, solvent oxidation current, electrolyte concentration, or calibrated transport law enters this objective. Radius is included in the GP input but absent from the target formula. PF6, BF4 and ClO4 labels receive exactly the same zero parasitic penalty; the generator cannot physically distinguish their performance. Assigned properties are held fixed while temperature varies.

| BO step | Solvent | Electrolyte | j (mA/cm²) | T (°C) | Observed yield (%) | Best so far (%) |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | MeCN | LiClO4 | 20.0 | 10.0 | 79.59 | 85.59 |
| 2 | MeCN | LiClO4 | 8.0 | 10.0 | 76.82 | 85.59 |
| 3 | MeCN | LiClO4 | 24.0 | 10.0 | 72.21 | 85.59 |
| 4 | MeCN | LiClO4 | 4.0 | 10.0 | 69.66 | 85.59 |
| 5 | MeCN | LiClO4 | 16.0 | 25.0 | 92.03 | 92.03 |
| 6 | MeCN | LiClO4 | 12.0 | 25.0 | 89.59 | 92.03 |
| 7 | MeCN | LiClO4 | 20.0 | 25.0 | 86.57 | 92.03 |
| 8 | MeCN | LiClO4 | 8.0 | 25.0 | 87.08 | 92.03 |

The reported optimum is **MeCN/LiClO4, 16.0 mA/cm², 25.0 °C, 92.03%**. The unrounded recorded observation is available in the CSV. Its latent value is **91.3865%**, with favorable sampled noise **0.6475310419 percentage points**. Six grid points tie at the latent maximum: MeCN, either 12 or 16 mA/cm², 25 °C, and any of the three non-acetate electrolyte labels. The continuous formula maximum, 92.4865%, occurs at 14 mA/cm² and 30 °C, neither included in the grid. Thus the selected label is not a unique chemical optimum. The model cannot substantiate the requested explanation that it “minimizes overpotential and competitive solvent oxidation.”

[A/S] An additional comparison uses 30 paired seeds, four methods and 14 evaluations per method: **120 campaigns, 1,680 observations**. Methods share six starting points and stepwise noise; sampling never repeats a complete candidate. The source GP is compared with random search and two normalized GP variants using either assigned physical descriptors or categorical one-hot features. The latter variants share an amplitude-scaled ARD Matérn kernel, one restart, target standardization from observations only, and noise variance 1.44 adjusted to that target scale. Because several changes accompany scaling, source-versus-normalized differences are not a single-factor ablation. The physical-versus-one-hot comparison also changes descriptor dimension and geometry.

Recommendation regret is the known grid maximum minus the latent response at the point selected by maximum noisy observation. Lower values are better; units are synthetic yield percentage points.

| Method | Mean regret | SD | Median |
| --- | --- | --- | --- |
| random | 2.6404 | 2.2697 | 4.0000 |
| scaled_onehot_gp | 4.8088 | 2.9055 | 4.5415 |
| scaled_physical_gp | 5.1545 | 2.8864 | 4.3517 |
| source_gp | 3.0202 | 3.5123 | 0.5415 |

| Paired contrast | Mean difference | 95% bootstrap interval |
| --- | --- | --- |
| scaled_onehot_gp_minus_random | 2.1684 | [0.9445, 3.5377] |
| scaled_physical_gp_minus_random | 2.5141 | [1.5242, 3.6860] |
| source_gp_minus_random | 0.3799 | [-0.9302, 1.8865] |
| physical_minus_onehot | 0.3457 | [-0.4807, 1.2498] |

Intervals are 10,000 paired-seed bootstrap intervals for the mean difference. The source-GP-minus-random interval crosses zero. Both normalized variants have greater mean regret than random search under this 14-evaluation protocol; physical versus one-hot is inconclusive. **The claimed advantage of physical embedding is not supported here.** There were 2,275 captured GP convergence warnings across the audit fits. The summary field named `kernel_warnings` counts all scikit-learn `ConvergenceWarning` instances; message-level categories were not retained. No broad hyperparameter search, new reaction distribution, experimentally informed property uncertainty, or GP calibration study was performed. Retaining these negative results is more informative than optimizing the analysis until one strategy appears favorable.

### 2.2 Task B — molecular identity, 3D geometry, SASA and site claims

[M/A] The SMILES `COc1ccc2[nH]c(cc2c1)c3ccccc3` defines **5-methoxy-2-phenyl-1H-indole, C15H13NO**, with 17 heavy atoms and 30 atoms after explicit hydrogens. ETKDGv3 seed 42 and MMFF94, maximum 500 iterations, both completed successfully. MMFF provides a force-field minimum, not an electronic-structure or vibrational certification of a chemical minimum.

The atom mapping is explicit: N1=6, C2=7, C3=8, C3a=9, C4=10, C5=2, C6=3, C7=4, C7a=5, using **zero-based RDKit indices**. C2 carries phenyl and C5 carries methoxy; both have zero attached hydrogens. The available indole-core C–H sites are C3, C4, C6 and C7. Five additional aromatic C–H sites belong to phenyl, and the methoxy methyl carbon adds a tenth C-bearing-H candidate. A C3-versus-C2/C5 C–H comparison is structurally inapplicable to this substrate.

![Mapped target structure](figures/target_atom_map.png)

Figure 1. Saved structure with atom-index/core-label annotations. Highlights include C5, C2, C3 and the source heuristic's top pendant-phenyl atom. Coordinates, connectivity and mapping are retained in `target_identity.json`; atom indices must not be mistaken for conventional ring locants.

The repaired SASA calculation uses explicit-H radii returned by this RDKit installation: C 1.70, O 1.55, N 1.60 and H 1.20 Å, Lee–Richards surface integration and probe radius 1.40 Å. Total SASA is **454.6229919137123 Å²**. This is a geometric convention, not a calibrated acetonitrile-solvation model. Atom areas sum to the molecular total. Gasteiger charges are computed from molecular connectivity and do not become radical-cation spin densities or frontier orbital populations merely because a conformer is present. [FreeSASA method](https://pmc.ncbi.nlm.nih.gov/articles/PMC4776673/); [RDKit charge implementation](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/PartialCharges/GasteigerCharges.cpp).

The retained source heuristic is $H_i=-0.6q_i+0.4A_{H_i}/15$. Its weighting and area scale are uncalibrated, and source code uses only the first attached H even for methyl. Exact repaired output is:

| Atom index (0-based) | Mapped label | Hybridization | Gasteiger q (e) | One H SASA (Å²) | Heuristic index |
| --- | --- | --- | --- | --- | --- |
| 15 | phenyl:15 | SP2 | -0.0616 | 26.676 | 0.7483 |
| 14 | phenyl:14 | SP2 | -0.0622 | 26.655 | 0.7481 |
| 13 | phenyl:13 | SP2 | -0.0616 | 26.647 | 0.7476 |
| 4 | C7 | SP2 | -0.0343 | 24.030 | 0.6614 |
| 10 | C4 | SP2 | -0.0106 | 21.407 | 0.5772 |
| 0 | methoxy-CH3 | SP3 | 0.0775 | 21.579 | 0.5289 |
| 8 | C3 | SP2 | -0.0285 | 18.755 | 0.5172 |
| 12 | phenyl:12 | SP2 | -0.0529 | 17.550 | 0.4997 |
| 16 | phenyl:16 | SP2 | -0.0529 | 17.406 | 0.4959 |
| 3 | C6 | SP2 | -0.0179 | 14.600 | 0.4001 |

The top atom, **15**, is on the pendant phenyl ring, not indole C3. The C3 score is **0.5172**, versus **0.7483** for atom 15. For the winner, approximately 95% of the score comes from the area term, illustrating how the assumed scale drives ranking. These outputs cannot explain a C3-selective SET/HAT pathway. SET is a molecular electron-transfer event, and later deprotonation, radical trapping, adsorption and competing barriers require separate mechanistic evidence.

[M] Additional ETKDG sampling (seed 20260927) generated 32 conformers; all converged within 1,000 MMFF94 iterations. MMFF energies span **32.214649–32.903951 kcal/mol** and total SASA **436.033643–457.348241 Å²**. These are force-field energies with their own reference, not reaction barriers. Repeated conformational wells are retained; the set is not an equilibrium ensemble. The audit averages over attached H atoms for methyl bookkeeping and restricts site ranking to the nine aromatic C–H carbons.

| Mapped label | Mean score ± SD | Mean H SASA ± SD (Å²) | Highest-score count / 32 |
| --- | --- | --- | --- |
| C6 | 0.4856 ± 0.1225 | 17.808 ± 4.593 | 0 |
| C7 | 0.6650 ± 0.0044 | 24.165 ± 0.164 | 0 |
| C3 | 0.5179 ± 0.0140 | 18.779 ± 0.524 | 0 |
| C4 | 0.4833 ± 0.1240 | 17.886 ± 4.649 | 0 |
| phenyl_atom_12 | 0.4990 ± 0.0133 | 17.521 ± 0.497 | 0 |
| phenyl_atom_13 | 0.7506 ± 0.0046 | 26.763 ± 0.173 | 2 |
| phenyl_atom_14 | 0.7558 ± 0.0051 | 26.944 ± 0.192 | 27 |
| phenyl_atom_15 | 0.7497 ± 0.0052 | 26.728 ± 0.196 | 3 |
| phenyl_atom_16 | 0.4996 ± 0.0155 | 17.544 ± 0.580 | 0 |

Pendant-phenyl atoms 14, 15 and 13 rank first in 27, 3 and 2 conformers respectively. C3 ranks first in none. These frequencies are sampling diagnostics, not regioisomer probabilities. Gasteiger charges remain invariant across conformers; SASA changes with geometry. On the lowest-MMFF conformer, altering probe radius and the SASA algorithm also changes the highest-scoring atom:

| Probe (Å) | Algorithm | Total SASA (Å²) | Top atom index | Top score |
| --- | --- | --- | --- | --- |
| 1.2 | LeeRichards | 421.495 | 14 | 0.6745 |
| 1.2 | ShrakeRupley | 416.629 | 15 | 0.6932 |
| 1.4 | LeeRichards | 454.086 | 14 | 0.7550 |
| 1.4 | ShrakeRupley | 449.548 | 15 | 0.7619 |
| 1.8 | LeeRichards | 521.291 | 13 | 0.9320 |
| 1.8 | ShrakeRupley | 515.741 | 15 | 0.9719 |

[Q] To add an actual electronic-structure calculation, the lowest-MMFF conformer underwent neutral **GFN2-xTB/ALPB(acetonitrile)** tight optimization, followed by +1 and −1 doublet single points at exactly the neutral nuclei. All three jobs terminated normally; neutral optimization converged. Input/output coordinates, partial charges, JSON and logs are retained.

| State | Charge / multiplicity | Energy (Eh) |
| --- | --- | --- |
| neutral | 0 / 1 | -45.759576727469 |
| cation | +1 / 2 | -45.373191206070 |
| anion | −1 / 2 | -45.980151497784 |

The fixed-nuclei model charge-removal difference is **10.5141 eV**, and the electron-addition difference is **6.0021 eV**, using 27.211386245988 eV/Eh. Total atomic charges match 0, +1 and −1 within 10⁻⁶ e. C3's Mulliken removal response $q_{+1}-q_0$ is **0.03769962 e**. This redistribution is not a site oxidation potential or a reaction barrier. Each charged state has its own equilibrium ALPB response; no nonequilibrium solvent, ionic relaxation, Hessian, thermal correction, reference electrode or measured Eox calibration is included. GFN2-xTB is semiempirical, not DFT. [GFN2-xTB primary method](https://pubs.acs.org/doi/10.1021/acs.jctc.8b01176).

### 2.3 Task C — coordination proxies, synthetic barriers and group transfer

[S/A] Four metals and four coordination labels produce **16 rows**, not the 36 experimentally synthesized variations claimed in a comment. No experimental dataset or catalyst structure is loaded. Nominal d counts, electronegativity values, ligand-field strengths and metal-specific baseline scores are assigned. The coordination labels are not fully defined POP structures and do not specify oxidation state, axial ligation, spin, protonation, support geometry or electrochemical potential.

$$
\varepsilon_d^{proxy}=-1.2L_f-1.5(\chi_M-1.8),\qquad
G_{score}=B_M+4.5(\varepsilon_d^{proxy}+1.85)^2+\zeta,\quad\zeta\sim N(0,0.35^2).
$$

The quadratic has its minimum at the assigned −1.85 eV proxy; in a barrier plot this is U-shaped. Calling it a Sabatier volcano is an analogy to activity maxima, not an adsorption-energy or rate derivation. For isolated sites, a single scalar called a d-band center does not replace a calculated local electronic spectrum. The metal-specific bases (Cu 19.5, Co 23.0, Ni 21.0, Pd 17.0) already favor Pd before learning.

| Rank | Metal | Coordination label | d count | Assigned proxy (eV) | Synthetic score (kcal/mol) |
| --- | --- | --- | --- | --- | --- |
| 1 | Pd | N2O2_salen | 8 | -1.680 | 17.09 |
| 2 | Pd | N3C1_porphyrin | 8 | -1.920 | 17.20 |
| 3 | Pd | N2C2_network | 8 | -1.440 | 17.88 |
| 4 | Pd | N4_planar | 8 | -2.280 | 17.95 |
| 5 | Cu | N4_planar | 9 | -1.830 | 19.45 |
| 6 | Cu | N3C1_porphyrin | 9 | -1.470 | 20.46 |
| 7 | Ni | N4_planar | 8 | -1.845 | 21.52 |
| 8 | Cu | N2O2_salen | 9 | -1.230 | 21.68 |
| 9 | Ni | N3C1_porphyrin | 8 | -1.485 | 21.74 |
| 10 | Ni | N2O2_salen | 8 | -1.245 | 22.02 |
| 11 | Cu | N2C2_network | 9 | -0.990 | 22.36 |
| 12 | Co | N4_planar | 7 | -1.800 | 23.31 |
| 13 | Co | N3C1_porphyrin | 7 | -1.440 | 23.89 |
| 14 | Ni | N2C2_network | 8 | -1.005 | 24.13 |
| 15 | Co | N2O2_salen | 7 | -1.200 | 24.90 |
| 16 | Co | N2C2_network | 7 | -0.960 | 26.01 |

The lowest observed source score is **Pd/N2O2_salen, 17.09 kcal/mol**, with proxy **−1.680 eV**. It is not a calculated $\Delta G^\ddagger$, and no Eyring rate or catalyst recommendation is inferred. The Extra Trees model has 100 trees and seed 42. Its reported feature importances are:

| Descriptor | Extra Trees importance |
| --- | --- |
| d_electrons | 0.1179 |
| metal_EN | 0.7578 |
| coord_N | 0.0445 |
| coord_O | 0.0170 |
| d_band_center_proxy_eV | 0.0628 |

Metal electronegativity importance, **0.7578**, exceeds d-count importance **0.1179** and proxy importance **0.0628**. These values reflect the fabricated target, correlated encodings and chosen fitted model. They are not causal measures of ligand-field control. In-sample R² is **1.0**; this is insufficient evidence of predictive transfer.

[A/S] Leave-one-metal-out and leave-one-coordination-type-out evaluations use four folds per split, 12 training rows and four test rows per fold. Extra Trees, a training-mean predictor and a standardized ridge baseline are fitted using training folds only. Results are unweighted means of the four fold metrics, in the source's synthetic score units:

| Held-out group | Model | Macro MAE | Macro RMSE |
| --- | --- | --- | --- |
| coordination_type | extra_trees | 0.8480 | 0.9303 |
| coordination_type | mean | 2.3249 | 2.7633 |
| coordination_type | ridge | 0.8663 | 1.0043 |
| metal | extra_trees | 2.7869 | 2.8272 |
| metal | mean | 2.9250 | 3.0952 |
| metal | ridge | 1.1440 | 1.2637 |

Extra Trees transfer is much poorer across unseen metals (MAE **2.7869**) than across unseen coordination labels (**0.8480**); the ridge baseline achieves lower metal-held-out MAE (**1.1440**) on this dataset. Four groups and artificial labels do not establish generalization to real catalysts. Task C also lacks its own NumPy seed: the full run is repeatable because Task A sets and consumes the shared random stream, but standalone C results depend on external RNG state. A future dataset should remove this coupling and preserve structural/provenance groups.

### 2.4 Task D — coupled boundary flux and a reduced channel ODE

[S/A] The numerical model assumes width 10 mm, height 0.5 mm and length 100 mm, giving **0.5 mL** volume and **10 cm²** of one active electrode. Inputs are $C_{in}=50$ mol/m³ (0.05 M), $D=1.2\times10^{-9}$ m²/s, $j_0=0.05$ A/m², $T=298.15$ K, $\alpha_a=\alpha_c=0.5$ and $\eta=0.45$ V. None is measured here. The code assumes a one-electron flux relation and treats the full reactant loss as useful output.

The executed transport correlation is $k_m=0.67D[\gamma/(DL)]^{1/3}$ with $\gamma=6u/h$. It differs from the header's 1.85 Sherwood expression. The implemented BV normalization uses inlet concentration, although the introductory equation names bulk concentration. Writing $a=Fk_m$, $b=(j_0/C_{in})e^{\alpha_aF\eta/RT}$ and $r=j_0e^{-\alpha_cF\eta/RT}$ gives

$$
C_s=\frac{aC_b+r}{a+b},\quad
\frac{dC_b}{dz}=-\frac{k_ma_s}{u}(C_b-C_s),\quad a_s=1/h.
$$

At fixed overpotential all coefficients are constant, so elimination of $C_s$ yields a linear ODE with analytical solution

$$
C_b(z)=C_{eq}+(C_{in}-C_{eq})e^{-k_{eff}a_sz/u},\quad
C_{eq}=r/b,\quad k_{eff}=\frac{k_mk_s}{k_m+k_s},\quad k_s=b/F.
$$

This is a reduced boundary-flux/axial-ODE coupling, not a solved velocity/concentration/potential-field PDE. BV is exponential in the prescribed overpotential; that fact does not make the fixed-parameter concentration equation nonlinear. The constant reverse term assumes fixed product activity, without product or counterelectrode transport.

| Flow (µL/min) | Velocity (mm/s) | Residence (s) | k_m (m/s) | Conversion (%) | STY (mmol/L/h) |
| --- | --- | --- | --- | --- | --- |
| 100 | 0.33 | 300.0 | 2.588e-06 | 77.55 | 465.3 |
| 300 | 1.0 | 100.0 | 3.732e-06 | 50.66 | 911.83 |
| 600 | 2.0 | 50.0 | 4.702e-06 | 35.52 | 1278.87 |
| 1200 | 4.0 | 25.0 | 5.924e-06 | 23.8 | 1713.41 |

The source STY calculation is dimensionally consistent as **reactant disappearance** per reactor volume and time. It becomes product STY only if unit selectivity and the stated electron/product stoichiometry are established. At increasing flow, $k_m$ rises by a factor of about 2.29 while residence time falls twelvefold. Conversion therefore decreases while disappearance throughput increases. This does not identify an industrial optimum or confirm useful product production.

| Flow (µL/min) | Damköhler number | Mass-transfer resistance fraction | Implied current (mA) | (D/k_m)/h |
| --- | --- | --- | --- | --- |
| 100 | 1.4939 | 0.9622 | 6.2354 | 0.9275 |
| 300 | 0.7064 | 0.9464 | 12.2192 | 0.6431 |
| 600 | 0.4389 | 0.9334 | 17.1378 | 0.5104 |
| 1200 | 0.2718 | 0.9175 | 22.9610 | 0.4051 |

The mass-transfer resistance fraction is $(1/k_m)/(1/k_m+1/k_s)=k_s/(k_m+k_s)$. It falls from **0.9622 to 0.9175** with increasing flow. All four conditions are predominantly mass-transfer-controlled under the assumptions, but high flow modestly reduces the relative transport resistance. The declining Damköhler number explains the conversion trend more accurately than claiming that high flow worsens mass transfer.

Independent DOP853 integration and the analytical outlet agree to **8.38×10⁻¹² mol/m³**. Integrating $j(z)$ over the electrode with a 201-point trapezoidal grid agrees with $Fq(C_{in}-C_{out})$ within relative error **4.65×10⁻⁶**; adaptive quadrature independently closes this current balance. Implied current rises from **6.2354 to 22.9610 mA**. This is internal consistency under a one-electron model, not measured Faradaic efficiency. Overpotential is not the total cell voltage, so the script cannot supply full electrical energy per product mass.

The effective thickness $D/k_m$ is 0.405–0.928 times channel height, which raises an applicability question for a thin developing boundary-layer approximation. The prefactor and geometry need experimental or higher-fidelity transport validation. Additional saved scenarios comprise 48 flow/overpotential combinations (four flows and 0.05–0.60 V in 0.05 V steps) and 20 cases varying $j_0$ or D by factors of one-half or two. These are parameter sensitivity calculations, not uncertainty intervals estimated from measurements.

![Flow conversion and disappearance throughput](figures/flow_tradeoff.png)

Figure 2. The same assumed flow model on separate axes: conversion decreases as reactant-disappearance STY rises. No product selectivity or process cost is included.

### 2.5 Cross-task interpretation and validation limits

![Four-task quantitative audit](figures/production_audit.png)

Figure 3. A: 30-seed regret distributions; boxes span quartiles, center lines are medians, whiskers extend to 1.5 IQR and outliers are shown. B: aromatic-site heuristic means with ±1 conformer SD; orange denotes C3 and Ph labels retain saved atom indices. C: macro mean errors on held-out synthetic groups. D: uncalibrated overpotential scenarios. These panels distinguish statistical, geometric and equation-solving evidence; none is an experimental validation.

The calculations support useful negative conclusions: the physical embedding is not yet an effective acquisition policy in this benchmark; the geometry/charge heuristic cannot justify the requested site assignment; the SAC regressor mostly learns assigned metal trends; and a dimensional transport model still requires calibration and product-specific balances. These are actionable design findings rather than evidence of a failed chemical reaction, since no reaction was conducted.

## 3. Wet-Lab Implementation & Instrumentation Protocol

### 3.1 Instrument roles and pre-experimental records

[P/U] The following is a proposed measurement program, not an inventory or a recovered group SOP. Confirm access, training and calibration with the laboratory before execution. Keep a single run registry linking canonical structures, atom maps, batch identities, reagent/solvent provenance, cell dimensions, electrode surface treatment, time/current/potential traces and raw analytical files. Actual preparative partner identity, reaction stoichiometry, catalyst preparation and verified baseline yield were not supplied; they remain prerequisites rather than invented recipes.

| Instrument or configuration | Specific role | Required evidence before interpretation |
|---|---|---|
| Potentiostat with three-electrode cell | CV and potential referencing | Blank stability, reference protocol, resistance measurement |
| Undivided preparative cell; glassy-carbon anode, Pt-plate cathode | Initial validation configuration | Measured immersed area/gap, current density, charge integral, temperature |
| Ag/Ag+ reference in MeCN; ferrocene internal reference check | Non-aqueous potential scale | Reference filling solution, junction, calibration before/after measurements |
| n-Bu4NPF6 electrolyte | Proposed baseline electrolyte for measurement | Solubility, conductivity and blank oxidation window in the actual medium |
| HPLC-MS | Time-resolved composition after calibrated sampling | Separation, response factors, internal standard and sampling/quench validation |
| 1H/13C NMR, HSQC/HMBC and selective NOESY | Product and regioisomer identity | Complete assignments and consistent connectivity constraints |
| N2 sorption at 77 K; ICP-MS; HAADF-STEM; complementary XAS/XPS | Surface area, loading, dispersion and coordination | Independent batches, representative fields and measurement uncertainty |

The proposed n-Bu4NPF6 baseline is a measurement choice, not a result that reverses or endorses the synthetic LiClO4 optimum. All electrolyte choices require measured compatibility with the defined reaction. Literature establishes that CV and cell configuration are central to electrosynthesis; exact implementation still depends on the chosen chemistry. [Electrosynthesis methodology](https://www.nature.com/articles/s41570-022-00372-y).

### 3.2 CV, prospective BO and product assignment

Begin with **1 mM target substrate and 0.10 M n-Bu4NPF6 in MeCN** as proposed analytical starting concentrations. Use a polished glassy-carbon disk, Pt counter electrode and a nonaqueous-compatible Ag/Ag+ reference. Record electrolyte-only, substrate-only where interpretable, each intended reaction partner, and reaction-mixture traces. Determine the usable potential window from blanks before scanning the analyte; no absolute preparative potential is inferred from the synthetic descriptor table. Acquire 50, 100 and 200 mV/s scans and independent solution preparations, retain repeated scans that reveal fouling, and measure reference behavior before and after the series. Report solvent, temperature, reference composition, uncompensated resistance and any correction method. Distinguish reversible midpoint, anodic peak and onset potentials.

For preparative work, first select a reaction with identified partners and a verified SOP. Reproduce its baseline in three independent preparations. Use an undivided glassy-carbon/Pt cell only while testing for crossover or cathodic product loss; introduce a divided control if those processes are plausible in the defined system. Record area, gap, stirring, concentration, temperature, atmosphere and $Q=\int I\,dt$. A constant-current density and a controlled anode potential are distinct operating modes; do not claim both are independently fixed with one controller. The source's 16 mA/cm² is a candidate to evaluate only after the measurement window and baseline are known, not a validated starting recipe.

Freeze a matched BO/random/space-filling comparison before outcome collection. Randomize feasible run order, use the same valid-evaluation budget and repeat reference conditions across days. Permit only actual solvent/electrolyte combinations that remain soluble and stable under the chosen operating conditions. Retain failed reactions and distinguish them from instrument failures. Report assay yield, isolated yield, conversion and selectivity separately; calculate FE only after electron stoichiometry and product amount are established. A claimed optimizer benefit must survive independent repeats and exceed analytical and between-run variability.

For kinetics, collect samples at specified fractions of the baseline charge endpoint, for example 0, 0.25, 0.5, 0.75 and 1.0 times the registered charge. Validate quenching and internal-standard recovery. HPLC-MS is time-resolved offline analysis unless an online sampling interface and latency are actually implemented. MS supports composition, while calibrated chromatography quantifies species. Confirm regioisomers with 1H/13C assignments and HSQC/HMBC connectivity; use NOESY as complementary spatial evidence rather than the sole site assignment. Evaluate C3, C4, C6, C7 and pendant-phenyl substitution as applicable to the defined reaction. C2 and C5 have no C–H in this substrate.

### 3.3 POP-SAC characterization and mechanism controls

Start from one verified support/metalation route rather than synthesizing every abstract label. Prepare independently replicated batches, report metal-loading basis by ICP-MS after validated digestion, and characterize pore structure by N2 sorption at 77 K with a material-compatible degassing SOP and documented BET fitting interval. Use sufficiently representative aberration-corrected HAADF-STEM fields to assess dispersion. Isolated bright features support spatial dispersion; they do not alone prove that all metal is atomically dispersed or determine the working coordination environment.

Add XANES/EXAFS where accessible for coordination hypotheses and XPS for surface-state trends. Do not infer a working M–N_x–C_y–O_z geometry solely from a precursor name, BET area, ICP loading or the absence of diffraction peaks. Compare the support, metal-free modified support, relevant soluble metal and a comparable nanoparticle control. Record leaching and post-reaction structure; filtration/poisoning can support but not independently settle active-species identity. Normalize performance by geometric area and catalyst mass; site-normalized rates require defensible active-site counts.

Mechanistic work should discriminate hypotheses with concentration/time profiles, competitive or isotopic measurements where interpretable, and product/intermediate identification. Targeted electronic-structure calculations require actual catalyst/substrate geometries, charge/spin choices and competing pathways. A trained model on the 16 synthetic labels is not a substitute for those structures or for measured catalytic rates.

### 3.4 Flow validation, acceptance gates and stopping conditions

Measure wetted volume and residence-time distribution before reaction fitting. If the proposed 10 mm × 0.5 mm × 100 mm one-face geometry is adopted, verify its actual 0.5 mL volume and electrode exposure rather than treating nominal dimensions as measured facts. Test the four source flows, 100, 300, 600 and 1200 µL/min, with the same defined inlet composition. Record pressure, temperature, current, electrode potential and cell voltage, then collect time-binned outlet samples until a prespecified stability criterion is met. Steady-state collection must be demonstrated; an arbitrary elapsed number of nominal residence times is not sufficient by itself.

Fit transport/kinetic parameters only after calibrating composition and obtaining product-specific mass and charge balances. Reserve at least one flow/overpotential condition for prediction, compare the reduced model against a residence-time baseline, and estimate parameter identifiability rather than reporting only a best fit. Determine whether a developing-layer correlation applies in the measured geometry and whether the reverse reaction/product activity assumption matters. Track fouling, current drift, pressure rise and product recovery during sustained operation.

Advance only when a defined product, repeatable quantification, defensible selectivity, and predictive validation are available. Stop or redesign when structure assignments are unresolved, apparent optimizer gains vanish on replication, current cannot account for product formation, or separation/electrode degradation eliminates a throughput benefit. Competition submission and journal ambition do not replace these gates.

## 4. Undergraduate Milestone Gantt

### 4.1 A realistic four-year schedule

[P] Months below are relative to project start, not promises of institutional competition deadlines. The schedule is designed for a freshman learning chemistry and computation together; progression depends on supervision, course workload, instrument access and the evidence gates above.

![Four-year undergraduate research Gantt](figures/undergraduate_gantt.png)

| Phase | Time | Main work | Deliverable and advancement criterion |
|---|---|---|---|
| Freshman | Months 1–12 | TLC, chromatography, cell operation, lab records, Python/RDKit, CV basics | Reproducible baseline and an independently rerunnable calculation; retain failed runs |
| Sophomore | Months 13–24 | Measured-data registry, analytical calibration, small prospective BO or descriptor project | A usable group tool with traceable labels and a matched baseline comparison |
| Junior | Months 25–36 | One focused wet/dry question, structural/mechanistic controls, optional catalyst/flow work | A discriminating dataset and draft manuscript; claims survive group/temporal holdout or independent experiments |
| Senior | Months 37–48 | Reproduction, uncertainty, thesis, sustained operation if justified, handover | Defended thesis, complete raw-data/code archive and clearly delimited conclusions |

During year one, aim for reliable records and interpretable measurements rather than an impressive model score. During year two, a small well-controlled study is preferable to combining all four themes. During year three, choose catalyst, mechanism or flow depth according to accumulated evidence, not the synthetic ranking. Year four should include independent reproduction and documentation time instead of assuming every experiment will succeed on schedule.

### 4.2 Competitions, innovation projects and paper preparation

Prepare a reusable dossier containing a falsifiable problem statement, measured baseline, provenance, comparison design, negative results, uncertainty, budget and division of labor. Adapt it to Challenge Cup, the innovation/entrepreneurship competition named “Internet+” in the brief, and an undergraduate innovation training proposal only after checking the current university notice. Their eligibility, official names and deadlines are not established by this calculation. Neither an award nor a high-impact acceptance is guaranteed.

A credible manuscript should center one experimentally supported advance: for example a prospective reduction in valid experiments, a transferable redox predictor with reference-consistent data, a defensible coordination/selectivity relationship, or sustained flow productivity with full balances. The current prototype can contribute software and an audit appendix. It cannot yet supply those scientific claims. If the evidence remains negative, report its useful methodological lesson accurately and adjust the project question.

### 4.3 Deliverables, reproduction and source boundaries

The package contains original source/failure logs, a repaired executable, full machine-readable results, audit scripts, grouped tests, atomic coordinates/charges, bilingual figures and three report editions. The English and Chinese editions contain the same numerical tables. The source-named combined file `production_technical_report_bilingual.md` aligns corresponding sections. Root repository tests cover both the earlier report and this extension.

Run from the repository root in the recorded environment:

```text
python production/scripts/prepare_and_execute.py
python production/scripts/audit_pipeline.py --seeds 30
python production/scripts/plot_production.py --zh-font PATH_TO_CJK_FONT
python production/scripts/build_reports.py
python -m unittest discover -s tests -v
python production/scripts/validate_production.py
```

These commands refresh generated files in the working copy; use a separate clone when preserving a release. Activate the environment first, particularly for numerical DLL resolution on Windows. Set `XTB_EXE` if xTB is not on PATH. Source extraction is preserved already; no repeated installation is required. Numerical, integrity and report-architecture tests are not claims of experimental validity. Raw machine-specific executable/directory paths are normalized in public logs, with scientific values retained.

Primary sources used for context and software interpretation are listed below. No external measured training set or supporting-information experiment was reproduced. The solvent/electrolyte property dictionary remains unverified and is explicitly treated as assigned input.

1. Huang and colleagues (2020), porous-ligand Pd catalysis. [Publisher DOI](https://doi.org/10.1016/j.chempr.2020.06.020). Context only; metadata/indexed summary available, full text/SI not reconstructed.
2. Wang and colleagues (2024), Fe single-atom redox mediation. [Publisher abstract](https://onlinelibrary.wiley.com/doi/abs/10.1002/anie.202404295). Scientific alignment; not a live equipment inventory.
3. Shields and colleagues (2021), experimental Bayesian reaction optimization. [Nature article](https://www.nature.com/articles/s41586-021-03213-y). Comparison-design precedent, not evidence for this embedding.
4. RDKit developers, [SASA API](https://www.rdkit.org/docs/source/rdkit.Chem.rdFreeSASA.html) and [charge implementation](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/PartialCharges/GasteigerCharges.cpp). API and graph-based-charge interpretation.
5. Mitternacht (2016), FreeSASA. [Primary software paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4776673/). Surface algorithms and configurable geometric conventions.
6. Bannwarth, Ehlert and Grimme (2019), GFN2-xTB. [Primary method paper](https://pubs.acs.org/doi/10.1021/acs.jctc.8b01176). Method provenance, not validation of this target's electron-removal energies.
7. xTB developers, [command documentation](https://xtb-docs.readthedocs.io/en/latest/commandline.html). Charge/spin, optimization and single-point settings.
8. Leech and Lam (2022), electrosynthesis practice. [Methodological review](https://www.nature.com/articles/s41570-022-00372-y). General analytical/cell context; exact proposed measurements above are newly designed, not copied group SOPs.

**Research position.** This is a complete execution, repair, numerical audit and proposed validation report under the supplied architecture. Production readiness remains unestablished. The next scientific step is a chemically defined, calibrated and prospective measurement program.
