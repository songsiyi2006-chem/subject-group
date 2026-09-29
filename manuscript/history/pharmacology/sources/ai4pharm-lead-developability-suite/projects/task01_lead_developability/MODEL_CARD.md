# Model card · Task 1 v1.0.0

## Intended use and validation status

An auditable medicinal-chemistry demonstration for hypothesis generation and software benchmarking. It has no labeled training set, fitted ADMET coefficients, calibration study, held-out accuracy estimate, or clinical qualification. Software tests establish implementation behavior, not biological predictive validity. Clinical archetype labels never enter scoring.

| Component | Method | Status |
|---|---|---|
| MW, cLogP, TPSA, HBD/HBA, RotB, aromatic rings, Fsp³ | RDKit descriptors | Computed descriptors; toolkit-dependent |
| CNS-MPO | Original six desirability transformations | Published method with approximate input substitution |
| Custom MPO | cLogP optimum 2–4, other terms as above | User-requested variant; not Pfizer CNS-MPO |
| ESOL | Delaney coefficients, aromatic atom fraction | Published empirical model; RDKit implementation not revalidated |
| pKa, logD at pH 7.4 | Chemical-class assumptions + neutral fraction | Unvalidated; optional sourced overrides |
| hERG, Caco-2, HIA | Explicit descriptor formulas | Unfitted, uncalibrated exploratory scenarios |
| Oral developability | Equal mean of six desirabilities | Custom ranking, not bioavailability |
| Internal H-bond vectors | ETKDGv3 + converged MMFF geometry | Finite gas-phase sample |
| Absolute chameleonic index | No valid identification from available data | Always `null` |

## Continuous score definitions

Let `clip(x)=min(1,max(0,x))`. Define a decreasing ramp and an increasing ramp:

$$L(x;a,b)=\operatorname{clip}\left(\frac{b-x}{b-a}\right),\qquad U(x;a,b)=1-L(x;a,b).$$

A plateau window is

$$W(x;a,b,c,d)=\min\{U(x;a,b),L(x;c,d)\},\quad a<b\le c<d.$$

Each contribution is continuous in its numerical input and lies in [0,1]. Integer structural counts remain discrete descriptors.

| CNS-MPO input | Desirability |
|---|---|
| cLogP | L(x;3,5) |
| cLogD₇.₄ | L(x;2,4) |
| MW, Da | L(x;360,500) |
| TPSA, Å² | W(x;20,40,90,120) |
| HBD | L(x;0.5,3.5) |
| Most basic pKa | L(x;8,10) |

`cns_mpo_approx = sum(six terms)`, range 0–6. A molecule with no recognized basic site receives a pKa contribution of one; its pKa value remains null rather than being invented as zero. This convention inherits the rule detector's incompleteness. The custom score substitutes `W(cLogP;0,2,4,6)` for the cLogP term. Original CNS-MPO does **not** impose an optimum of 2–4. Sources: [Wager et al., 2010](https://doi.org/10.1021/cn100008c) and [independent implementation comparison](https://pmc.ncbi.nlm.nih.gov/articles/PMC8260158/).

## Ionization assumptions and overrides

The dominant matched base uses the highest assumed pKa; the dominant acid uses the lowest. Approximate class values are amidine/guanidine 11.5, aliphatic amine 9.0, aniline-like 5.0, aromatic nitrogen 4.5; carboxylic acid 4.5, tetrazole 4.9, sulfonamide 7.0, imide 8.5, phenol 10.0, thiol 9.5. These numbers are explicitly **assumptions**, not compound-specific experimental values. The exact SMARTS and site indices are inspectable in code/results. Amides and many sulfonyl-adjacent nitrogens are excluded from generic amines.

At pH 7.4, let `b=10^(pKa_base−7.4)` and `a=10^(7.4−pKa_acid)`; absent matched sites contribute zero to the respective variable.

$$f_0=\frac{1}{(1+b)(1+a)},\quad \log D_{7.4}\approx\mathrm{cLogP}+\log_{10}f_0.$$

`f_cation=b/(1+b)` and `q=b/(1+b)+a/(1+a)`. Here `q` is expected absolute charge burden in a dominant-site approximation, not net molecular charge. It can approach two for ampholytes. No coupled microstates, secondary pKas, explicit ionic partitioning, or tautomer standardization are included. A missing matched site does not prove a molecule is nonionizable; enolic acids and substituted heterocycles are important weaknesses. Salts, disconnected mixtures, net-charged parents and permanent quaternary zwitterions are rejected.

Override file example (numbers below are an illustrative input format, **not evidence for any compound's pKa**):

```json
{
  "Diazepam": {
    "basic_pka": 3.0,
    "acidic_pka": null,
    "logd74": 2.5,
    "source": "REPLACE with verified assay record, method, temperature and citation"
  }
}
```

Omitted fields retain assumptions. Null pKa removes that dominant ionization site. Overrides do not alter structure or the detected topological nitrogen sites. `logd74` overrides only logD; it does not refit neutral fraction or charge burden. Values are validated for finite, plausible numerical bounds and unknown fields/names are rejected. The manifest stores the exact override file content. The pKa sensitivity minimum and maximum vary acid/base assumptions independently by ±1, retaining any direct logD override; these are scenario bounds, not confidence intervals.

## Solubility

$$\log_{10}S\,[\mathrm{mol/L}]=0.16-0.63\mathrm{cLogP}-0.0062MW+0.066RotB-0.74AP.$$

`AP` is aromatic heavy atoms divided by all heavy atoms, not aromatic ring count. RotB uses RDKit's strict definition. Conversion: `S[µg/mL] = 10^logS × MW × 1000`. No arbitrary extra ring penalty or pH correction is added. These are structure-derived water-solubility estimates without a controlled pH, polymorph, salt, or formulation specification. The original descriptor implementation differs from RDKit. Values for MW >500, HBA >10, or HBD >5 carry an extrapolation caution; this flag is a transparent coarse screen, not a validated applicability-domain classifier. Source: [Delaney, 2004](https://doi.org/10.1021/ci034243x).

## Explicit exploratory ADMET equations

**None of the coefficients in this section was fitted or validated.** They encode qualitative hypotheses so users can inspect and replace them. They must not be described as validated QSARs.

Volume `V` is RDKit grid volume in Å³ of the first converged conformer. If 3D fails or is disabled, `V=1.2×MW` is an explicitly recorded heuristic fallback. Define `Vq=V(1+0.5q)`.

$$\log_{10}\left(\frac{P_{app}}{10^{-6}\mathrm{cm/s}}\right)=1.8-0.012TPSA-0.35q-0.45\log_{10}(V_q/300).$$

The numeric Caco-2 output is `10^logPapp`, in units of 10⁻⁶ cm/s. Transporters, paracellular flux, donor/receiver pH, unstirred water layers, and assay protocols are absent. This is a permeability ranking scenario, not a validated Caco-2 assay estimate.

$$\Delta G_{proxy}=0.025TPSA+1.5q+0.15HBD-0.30\mathrm{cLogP},\quad HIA_{proxy}=100\sigma((3-\Delta G_{proxy})/1.2).$$

`σ(z)=1/(1+exp(−z))`. The energy-like score has a nominal kcal/mol scale; it is not a computed transfer free energy. The HIA mapping represents a desolvation-only scenario and ignores dissolution, transport, metabolism and dose. HIA is not oral bioavailability; the two can differ greatly.

For the hERG heuristic, count pairs between detected basic nitrogens and aromatic rings whose minimum shortest-path separation is 3–8 bonds. Let `T=min(pairs/2,1)`, `R=min(aromatic_rings,5)`, and `E=clip(max(ESP,0)/0.25)`. Neutral-form Gasteiger charges yield a point-charge potential at aromatic centroids, `Σqᵢ/max(rᵢ,1Å)`, with the largest centroid value used. This omits dielectric scaling, protonation microstates, quantum polarization and receptor interactions. If unavailable, E=0 and `herg_esp_available=false` explicitly records the omitted term.

$$z=-3+0.70(\mathrm{cLogP}-2)+1.5f_{cation}+0.5R+0.8T+0.25E,\qquad hERG_{proxy}=\sigma(z).$$

The score lies in [0,1] but is **not a probability**, IC₅₀, pIC₅₀, clinical arrhythmia risk, or exposure-adjusted safety margin. The architecture is motivated by lipophilic/basic ligand features and hERG binding experiments; the coefficients are original illustrative choices. [Mutagenesis experiments](https://pubmed.ncbi.nlm.nih.gov/18987434/) support the importance of particular pore residues, not this numeric formula.

## bRo5 flexibility and H-bond interpretation

RDKit strict RotB omits ring bonds; macrocycles can therefore have substantial conformational freedom despite modest RotB. Aromatic ring count and Fsp³ are independent descriptors, with reference preferences ≤3 and ≥0.42. They do not guarantee oral absorption.

The reproducible gas-phase ETKDGv3 ensemble is optimized with MMFF; only converged conformers contribute to geometric outputs. If initial embedding fails, a random-coordinate retry permits cis as well as trans amides; this matters for cyclic N-methyl peptides. Each embedding attempt receives 30% of the total worker time budget (rounded down to whole seconds, minimum one); the parent process still enforces the overall timeout. The retry method is recorded. This is not exhaustive cis/trans or solvent sampling. A donor–hydrogen vector counts as internally paired when donor–acceptor ≤3.5 Å, hydrogen–acceptor ≤2.6 Å, angle D–H···A ≥120°, and D–A graph distance ≥4 bonds. Each donor–hydrogen vector is counted once using its nearest qualifying acceptor. Acceptor sites paired to any vector are deduplicated. Donor/acceptor feature assignment uses RDKit `BaseFeatures.fdef`, not a claim about a particular aqueous microstate.

The outputs include mean/min/max internal vector counts, unpaired donor vectors, unpaired acceptor sites and `chi_geometry_range_proxy=max(internal vectors)−min(internal vectors)`. Conformers are equally weighted, not a thermodynamic ensemble. One conformer yields zero range, not evidence of rigidity. Unpaired sites are **potentially accessible**, not demonstrably solvent-exposed. No solvent-accessible surface calculation, matched polar/nonpolar ensemble, experimental exposed polar surface area, or free-energy weighting is available. Therefore `absolute_chameleonic_hb_index=null` for every compound. Source for implemented geometry tools: [RDKit Book](https://www.rdkit.org/docs/RDKit_Book.html).

## Oral ranking and radar anchors

| Quantity | Fixed desirability |
|---|---|
| ESOL logS, mol/L | U(logS;−7,−3) |
| log₁₀ numeric Papp in 10⁻⁶ cm/s | U(logPapp;−1,1) |
| hERG heuristic | 1−hERGproxy |
| TPSA, Å² | L(TPSA;90,200) |
| Fsp³ | U(Fsp³;0,0.42) |
| Aromatic rings | L(rings;3,6) |

`Oral developability score = 100 × mean(six desirabilities)`. These are custom trade-offs, not a clinically calibrated success threshold. They may penalize successful bRo5 drugs because formulation, metabolism and transport are absent. HIA is not added to avoid counting another closely related passive-transport surrogate. The radar replaces aromatic-ring desirability with `CNS-MPO/6`; it is a profile display, not a visualization of the oral score's exact component sum. No parameters depend on clinical archetype.

## Output contract and reproducibility

`results_task1/developability_results.json` is the authoritative, full-precision record set. CSV is a rounded interchange export; null JSON values become blank CSV cells. `_proxy`, `_approx`, and `_assumed` identify epistemic status. Dimensional suffixes specify units; dimensionless scores use documented ranges. Clinical metadata and structure source fields accompany every row. Individual CNS and oral desirability contributions allow score auditing.

`group_summary.csv` reports n, means/medians and Ro5 counts descriptively. `run_manifest.json` records dependency versions, platform, seed, time budget, script/panel hashes, exact overrides, geometry status counts and SHA-256 hashes of results and figures. A successful run writes the manifest last. If a run fails midway, existing files may contain artifacts from different runs: treat them as a complete set only after all manifest checksums match. Do not run two processes into the same output directory simultaneously.

Numerical validation rejects invalid or duplicate structures and incorrectly assigned panel counts, and requires finite primary scores. Geometry failures remain explicit and can be made fatal with `--strict-3d`. Plot exports use 300 DPI and fixed scales; the SVG hash salt and jitter seed are fixed. Wall-clock timeouts and toolkit/platform changes can change 3D availability or conformers, so the seed alone does not guarantee byte-identical images across environments.

For actual QSAR validation, the next stage needs standardized measured endpoints, assay/context metadata, scaffold-separated train/validation/test sets, leakage checks, calibration and uncertainty assessment, applicability-domain analysis, and prospective experiments. This 30-compound historical panel alone cannot supply those claims.
