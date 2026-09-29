# Historical molecular learning and quantum evidence: ancestry, reference choice and validation limits

## L.1. A history audit with explicit units of evidence

The historical repositories contribute a useful comparison between increasingly ambitious scientific workflows and the evidence needed to support their conclusions. The central result is not a count of completed software modules. It is a set of traceable cases in which an implementation check, a converged calculation or a favorable reported metric does not establish the next scientific claim. These cases complement the detailed H₂ error chain in the main study: they reveal similar distinctions in solubility learning, conformer screening, variational Monte Carlo, electrochemical time series and charged molecular calculations.

We fixed the public repositories at the remote HEADs verified on 29 September 2026: `aqueous-solubility-ml-benchmark` at `5939547a936ed4baf7fa7ee9fac6fa6e90bfa8fe`, and `ai4chem-complex-scaffolds-benchmark` at `5e19c519dbd11430926828fddf5375b919af67f0`. The former contains one commit; the latter contains 54. Crucially, the entire aqueous commit is the root commit of the complex-scaffolds history. Five inherited Python files and four inherited PNG files have identical byte hashes across the two pinned trees. They are a shared scientific record, not an independent replication. The later repository adds new stages, numerical corrections, fresh-run acceptance reviews, a folder migration, and Phases 28–31. Moving files into project folders does not itself constitute a new calculation. The complete commit-subject inventory and review categories are preserved in [the history CSV](sources/complex_history.csv); this inventory is not a claim that every historical diff or every large trajectory was exhaustively revalidated.

Table L1. Pinned evidence, counting units and reproducibility boundaries.

| Evidence object | Historical amount | Unit counted | Status in this audit |
|---|---:|---|---|
| Aqueous repository | 1 | commit | Shared ancestor; no stored prediction rows |
| Complex-scaffolds repository | 54 | commits | Stage and audit history; no additive sample count |
| Identical inherited code and figure files | 9 | files | Exact byte-level duplication |
| Phase 01 | 500 | accepted force-field optimization records | 10 valid structures; original and audited rounded energies agree |
| Phase 11 | 100 | saved VMC production block means | 5 systems; 4,050 training epochs; no saved wavefunction weights |
| Phase 28 | 196,496 | published trace rows | 39,726 distinct within-trace timestamps; electrode count unverified |
| Phase 29 | 156 | historical xTB executions | 144 unique settings; 12 executions repeated prior settings |
| Phase 29 DFT crosscheck | 8 | converged single-point jobs | 4 charge-state pairs; 8 matched GFN comparison rows |

The source package retains original MIT notices, pinned paths, original SHA256 hashes, JSON pointers and CSV selectors in [audit.json](audit.json). Selected publisher-derived numerical aggregates retain source DOIs and workbook hashes; their inclusion does not transfer third-party rights to the repository MIT license. Publisher metadata was checked, but the live rights text could not be retrieved through the publisher identity redirect. No publisher article text, original figures or full workbook is bundled here. The deterministic [post-hoc script](sources/recompute_history.py) uses only the Python standard library and does not import the original scientific programs. It checks source hashes and recomputes the displayed arithmetic without additional electronic-structure jobs, neural training, molecular optimization or experiments. Raw solver logs were not all reparsed in this historical review; the retained individual result files and execution ledgers support the narrower claims below.

## L.2. Solubility: a plausible protocol is not a reconstructed benchmark

The original aqueous workflow combines 11 RDKit descriptors with a 2,048-bit Morgan fingerprint of radius 2. Its total feature dimension is 2,059. The code canonicalizes SMILES, removes unparsable structures, averages labels sharing a canonical SMILES, and applies shuffled five-fold cross-validation with seed 42. A separate ESOL holdout uses a random 20% split with the same seed. Before training the two head-to-head models, the script removes AqSolDB molecules whose canonical SMILES match a molecule in that ESOL test set. This is an explicit exact-identity overlap control. It does not establish scaffold separation, independence of source laboratories, disjoint analog series or harmonized measurement conditions. No scaffold-split result or repeated-seed performance distribution is stored in this repository.

The README reports five-fold R² = 0.815 and RMSE = 1.019 log units, followed by a same-test-set comparison of RF RMSE = 0.803 and gradient-boosting RMSE = 0.635 for 226 molecules, with 222 overlapping AqSolDB molecules removed. The corresponding reported R² values are 0.864 and 0.915. These remain **historical report-level values**. The committed tree contains code and raster figures, but no molecular prediction table, raw label snapshot, frozen split identifiers, fitted model or execution log that allows those metrics to be reconstructed at the pinned commit. We did not recover predicted values from image pixels or rerun the training. Therefore the reported 20.9% improvement must not enter a pooled performance analysis alongside fully preserved benchmarks as if it were independently recomputed.

The code also resolves an ambiguity in the historical “old ESOL model” wording. In the head-to-head script, both the four-descriptor random forest and the descriptor-plus-fingerprint gradient booster are fitted to the filtered AqSolDB training rows. The standalone `esol_model.py` uses ESOL training data, but that is a different execution pathway. The head-to-head comparison changes features and model class together; it does not isolate the causal contribution of additional training data, fingerprint representation, or estimator architecture. A single held-out random split is useful evidence under that protocol, but it cannot establish industrial reliability or structural extrapolation to complex natural products.

The paclitaxel discussion further illustrates why chemical interpretation requires a separate evidence layer. The code uses a non-stereospecified input string and 2D features. Neither a low predicted solubility nor a fingerprint bit demonstrates a particular intramolecular hydrogen bond or the burial of polar surface in a conformer ensemble. Even the printed reference-unit conversion requires correction: using the reported molecular weight of 853.9 g mol⁻¹, the stated 0.3–1 μg mL⁻¹ range corresponds to log₁₀ S = −6.454286 to −5.931407 in mol L⁻¹, not the source program’s printed −6.1 to −5.6 interval. This is an arithmetic audit of a stated range, not independent validation of the experimental range itself. The primary AqSolDB descriptor article identifies the dataset, but it does not validate this repository’s particular fitted model or its molecular explanation. [@aqsoldb].

## L.3. Conformer screening measures sampling behavior, not a trained model’s accuracy

The first complex-scaffolds benchmark contains 11 input records. One original helicene input, M09, fails kekulization and is retained as a failed input; M09R is its separately identified repaired reference. Thus the completed panel consists of 10 valid structures. Each has 50 embedded and accepted, converged optimization records in the audited output. Nine structures use MMFF94, whereas the boron-containing M07 uses UFF because of force-field coverage. The total of 500 is a number of accepted optimization records, not a count of unique minima, independent experiments or quantum calculations. The old and audited JSON files have identical rounded minimum energies, maximum energies and energy ranges for all ten valid structures. A corrected acceptance policy matters for future inputs, but these preserved results do not show an energetic improvement caused by the correction.

Table L2. Audited conformer ranges, with force-field and fragment boundaries retained. Energies are in kcal mol⁻¹.

| Structure | Force field | Accepted records | Minimum energy | Maximum energy | Stored range | Geometric IMHB count |
|---|---|---:|---:|---:|---:|---:|
| M01 | MMFF94 | 50 | 67.37 | 110.42 | 43.06 | 2 |
| M02 | MMFF94 | 50 | 106.33 | 114.92 | 8.59 | 0 |
| M03 | MMFF94 | 50 | 29.16 | 41.33 | 12.17 | 1 |
| M04 | MMFF94 | 50 | 71.34 | 82.83 | 11.48 | 1 |
| M05 | MMFF94 | 50 | −0.15 | 7.72 | 7.87 | 0 |
| M06 | MMFF94 | 50 | 38.04 | 65.47 | 27.43 | 1 |
| M07 | UFF | 50 | 78.56 | 78.82 | 0.26 | 0 |
| M08 | MMFF94 | 50 | 8.43 | 15.27 | 6.84 | 0 |
| M09R | MMFF94 | 50 | 97.19 | 97.19 | 0.00 | 0 |
| M10 | MMFF94 | 50 | 123.99 | 139.09 | 15.11 | 0 |

The table uses the stored ranges rather than subtracting already rounded endpoints; independent rounding explains differences of 0.01. Absolute energies cannot be compared between different chemical compositions as a stability ranking, and the UFF and MMFF94 entries cannot be pooled as one energy scale. M03 is a disconnected two-fragment input: its graph and 3D optimization describe the largest fragment, whereas some bulk descriptors describe the full assembly. Calling that record a fully connected degrader conformer would combine incompatible molecular objects. M09R’s rounded zero range also does not prove that 50 independently distinct conformational minima were discovered.

The historical flag named `ENTROPY` is assigned from an energy-range threshold. An energy range is not a thermodynamic entropy: neither populations, degeneracies nor an equilibrium measure are established by the flag. Similarly, geometrically counted intramolecular hydrogen-bond contacts are structural diagnostics under a distance-and-angle rule, not measured bond free energies. The graph neural network checks confirm feature and tensor compatibility for the supplied graphs, but do not include trained property predictions. Consequently, this benchmark supports an input-applicability and sampling audit, not a claim that a 2D neural network has been quantitatively defeated or that 3D learning has won. The contrast with the fully supervised FreeSolv case in the main manuscript is deliberate: one checks representation plumbing, the other measures label prediction under a specified split.

## L.4. A neural-wavefunction counterexample with matched numerical and accuracy checks

Phase 11 provides the strongest historical match between passing internal checks and a retained accuracy failure within one completed run. The fresh run stores five systems, each with 2,048 walkers and 20 production block means. The training schedule comprises 1,200 epochs for equilibrium H₂, 550 each for three stretched H₂ cases, and 1,200 for He, for 4,050 epochs in total. The preserved result contains an antisymmetry maximum deviation of zero and an automatic-differentiation versus finite-difference Laplacian relative error of 3.515304 × 10⁻⁷. These diagnostics validate specified numerical operations. They do not bound the variational ansatz error or the stochastic uncertainty of a final energy.

For every system, we independently reconstruct Ē = Σᵢ Eᵢ/20 and SE = s(Eᵢ)/√20 from the saved block means. The comparison reference is a computed FCI aug-cc-pV(T,Q)Z complete-basis extrapolation preserved in the historical run, not an error-free continuum energy. Distances in the system labels are **bohr**, whereas the QuantumEqui bond-scan distances are in ångström. The two studies therefore cannot be joined using the numeric bond label alone, and their different basis/reference definitions preclude treating them as identical-reference replications.

Table L3. VMC point estimates and independently reconstructed block standard errors. Energies are Hartree; errors and SE are mEh.

| System | VMC energy | Computed reference | Signed error | Block SE | Absolute-error threshold | Point gate |
|---|---:|---:|---:|---:|---:|---|
| H₂, R = 1.4011 bohr | −1.173979102 | −1.174441940 | 0.462837 | 0.321004 | 1.600000 | Pass |
| H₂, R = 2.5 bohr | −1.093561268 | −1.093943014 | 0.381746 | 0.203228 | 1.600000 | Pass |
| H₂, R = 4.0 bohr | −1.015835784 | −1.016336946 | 0.501162 | 0.133787 | 1.600000 | Pass |
| H₂, R = 6.0 bohr | −0.999316513 | −1.000740584 | 1.424071 | 0.390502 | 1.600000 | Pass |
| He | −2.902038144 | −2.903699030 | 1.660887 | 0.564223 | 1.600000 | Fail |

Four point estimates pass the unchanged 1.6 mEh gate; He fails. Uncertainty does not convert the failing point estimate into a pass, and a passing point estimate does not establish a confidence-bound guarantee. The earlier acceptance audit, using the same 20 stored block means, reported an iid-block-bootstrap signed-error interval of [0.595197, 2.746738] mEh for He and [0.667166, 2.165752] mEh for stretched H₂ at 6 bohr. These historical intervals are retained as descriptive reanalyses, not new Monte Carlo observations. They assume independent blocks and omit reference extrapolation uncertainty. No claim of autocorrelation removal follows from having 20 block means.

The cusp evidence is weaker still. The saved electron–nucleus and electron–electron directional slopes are −2.430565 and −1.862306, but they were obtained from a finite-radius single-ray log-density fit, rather than spherical averaging followed by the coalescence limit. A later angular-average diagnostic was tested on an analytic function, not on the trained network. The wavefunction `state_dict` was not persisted, so this historical result cannot support a frozen-model cusp remeasurement or independent resampling without retraining. Corrected future source labels cannot retroactively validate an old figure. This is a concrete example of why saving loss histories and numerical energies alone is insufficient for a reproducible learned wavefunction claim.

## L.5. Published time series: repeated timestamps and estimand sensitivity

Phase 28 reanalyzes numerical source data accompanying two primary organoelectrocatalysis studies. The 2023 Nature record supplies constant-potential and constant-current traces; the 2026 Nature Synthesis record supplies a longer constant-current trace. These are previously published experimental observations, not experiments performed for this manuscript. Source workbook hashes and sheet/column coordinates are preserved, including Figure2g columns A:B and D:E for the earlier study and sheet 2h columns A:B for the later study. [@yang2023; @yang2026].

The three traces contain 196,496 numerical rows but only 39,726 unique timestamps when counted separately within each trace. Equal-time entries must not be called independent electrodes or independent degradation events. Their group sizes can change the weight of an endpoint statistic; the retained analysis therefore includes timestamp-balanced medians as well as raw-row medians. The 2023 constant-potential current ordinate has an unresolved normalization, so it is reported in source-ordinate units rather than silently converted to an absolute or area-normalized current. End of record is censoring of observation, not an observed failure time.

Table L4. Trace multiplicity and sensitivity of endpoint changes to analysis windows. For constant-current traces, change is in V; for the constant-potential trace, change uses its unresolved source current ordinate.

| Trace | Rows | Unique timestamps | Duplicate fraction | Balanced 1 h change | All-window change range | Startup ≥10 h change range |
|---|---:|---:|---:|---:|---|---|
| Nature 2023, constant potential | 43,112 | 13,246 | 0.692754 | −0.156250 | [−0.156250, 0.0428125] | [−0.023750, 0.021250] |
| Nature 2023, constant current | 18,750 | 7,707 | 0.588960 | 0.044250 | [−0.004000, 0.044250] | [−0.004000, 0.009250] |
| Nature Synthesis 2026, constant current | 134,634 | 18,773 | 0.860563 | 0.017000 | [0.003000, 0.017000] | [0.003000, 0.006000] |

The 72 saved window specifications expose a substantial estimand choice. For the 2023 constant-current trace, the default balanced one-hour difference is 44.25 mV, while windows excluding at least the first 10 h yield differences between −4.00 and 9.25 mV. For the 2026 trace, the analogous range contracts from 3–17 mV over all specifications to 3–6 mV after that exclusion. These changes use the same preserved traces. They are not between-electrode uncertainty intervals, and choosing the smallest value retrospectively would be an analysis choice, not evidence of improved durability. A useful manuscript figure should display the whole window grid and mark the startup policy explicitly.

Separate source sheets preserve 50 time-dependent mean/SD entries, with the published caption reporting four independent experiments and without individual replicate values. Their final chlorine-production mean is 0.338190 mol at 24 h, and the reported final Faradaic-efficiency mean is 99.792360%. Historical mean intervals rely on assumptions about replicate independence and distribution. Serial comparisons additionally depend on the unknown covariance between timepoints; the repository retains a correlation-sensitivity analysis rather than inventing replicate trajectories. High selectivity over an observed window therefore does not identify a degradation hazard, a lifetime distribution, a mechanism or an industrial replacement schedule. This distinction is particularly relevant when a digital-twin optimizer uses such quantities as objectives: the objective’s definition and sampling unit must precede its optimization.

## L.6. Charged molecular energies: a reference-aware method comparison

Phase 29 extends an initial 12-call molecular xTB pilot with a 144-call matrix: 12 molecular species, two GFN variants, three environments and two charge/spin states. The extension repeats the 12 initial settings, so there are 156 historical executions but only 144 unique settings across both rounds; 132 settings are new in the extension. All 144 extension calls have normal termination and converged SCC recorded, and all retain floating-point numerical warnings. The 72 state differences are reconstructed here from the saved electronic energies. Their arithmetic matches the preserved table, but matching subtraction does not prove physical accuracy in solution or at an electrode.

The independent crosscheck contains eight converged PBE0 jobs at fixed MMFF cation geometries: cation/radical pairs for Q01–Q03 in def2-SVP and an additional Q01 pair in def2-TZVP. The molecular correspondence is Q01/P01 pyridine, Q02/P04 3-methoxypyridine and Q03/P09 pyridine-3-carbonitrile. Geometry hashes agree across each paired method comparison. ΔE is defined as E(radical, charge 0) − E(cation, charge +1), with the same nuclei within a pair. These are model-state energy differences, not experimental redox potentials, interfacial barriers or product selectivities.

Raw GFN–PBE0 offsets are approximately −5.86 to −5.59 eV for GFN1 and −4.89 to −4.81 eV for GFN2. Because the methods use different electronic energy references for different charge states, those raw offsets must not be labeled prediction errors against physical ground truth. A more interpretable, limited diagnostic removes the common Q01 offset: δM(Q) = ΔEM(Q) − ΔEM(Q01). Its cross-method double difference is δGFN(Q) − δPBE0(Q). This subtraction cancels a method-wide constant charge-state offset while retaining molecule-dependent differences; it does not calibrate redox potentials.

Table L5. Matched def2-SVP reference contrasts and their double differences, in eV. PBE0 is a comparator, not experimental truth.

| Molecule | PBE0 ΔE | PBE0 Q01-centered | GFN1 Q01-centered | GFN2 Q01-centered | GFN1 double difference | GFN2 double difference |
|---|---:|---:|---:|---:|---:|---:|
| Q01/P01 | −4.957325 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 |
| Q02/P04 | −4.776458 | 0.180867 | 0.166602 | 0.162533 | −0.014265 | −0.018334 |
| Q03/P09 | −5.689793 | −0.732468 | −0.475732 | −0.678319 | 0.256736 | 0.054149 |

The two non-parent contrasts show that agreement is molecule dependent. The GFN1 double difference for Q03 is 0.256736 eV; GFN2 gives 0.054149 eV. Conversely, the Q02 double differences are −0.014265 and −0.018334 eV. Q01’s def2-SVP to def2-TZVP increment is −0.039814 eV. That single-molecule basis check cannot establish basis convergence for the other molecules. The radical S² values lie between 0.770825 and 0.776918, compared with 0.75 for a pure doublet, and wavefunction stability was not checked. No geometry relaxation, diffuse-basis study, thermal correction, counterion, electrode reference or explicit interface is included in this DFT comparison.

The larger matrix also preserves ranking instability: the gas-phase GFN1/GFN2 comparison reverses 6 of 66 molecular pairs. GFN1 gas versus acetonitrile reverses 7 of 66. These are within-panel sensitivity counts at fixed stored geometries, not out-of-sample accuracy statistics. They provide a useful companion to the main manuscript’s H₂ reference-bias decomposition: the H₂ chain has an explicit finite-basis FCI comparator, whereas the charged-molecule chain first requires a physically defensible reference alignment and still lacks external target validation.

## L.7. Historical acceptance gates and a restrained innovation claim

Table L6. Selected later-stage evidence boundaries; these units must not be added to a single “total calculation” count.

| Stage | Preserved numerical scope | Strongest defensible conclusion | Missing or failed layer |
|---|---|---|---|
| 24 | 312 simulated measurements; 312 QC passes | Analytical and acquisition workflow executes on synthetic responses | 0 new experiments; no hardware or chemical validation |
| 25 | 25 audited preliminary xTB minima across pilot jobs | Local molecular minima and partial imine calculations | Selectivity undetermined; accepted DFT TS/IRC absent |
| 26 | 48 condition entries; 3 dry slab seeds | Interface study design and starting objects | 0 complete solvated interfaces |
| 27 | 6-member design; 3 excitation windows | Literature-coordinate and analysis preparation | 0 new nonadiabatic trajectories |
| 29 | 8 converged PBE0 single points | Same-geometry method comparison | 0 validated TS/IRC; same-pair switching unproved |
| 30 | Carbon residual −4.11658 × 10⁻¹⁹ mol s⁻¹; heat residual −7.43849 × 10⁻¹⁵ W | Conservation and numerical checks of a reduced scenario | Uncalibrated; no DFT/AIMD or industrial validation |
| 31 | 18 converged MMFF conformers | Bounded local molecular engineering pilot | 0 FEP trajectories; 0 accepted lead candidates |

Earlier phases receive the same treatment in the machine-readable coverage ledger. Phase 04’s independent saddle refinement does not retroactively make the original NEB converged or supply a verified IRC; Phase 05 retains multiple imaginary modes and near-zero yield; Phase 06 finishes a trajectory but undersamples the product. Phases 03, 07, 08 and 12–16 have pending or incompletely accepted scientific scopes. Phases 09–10 are protocol/control models, Phases 17–19 have bounded atomic, numerical or geometric validations, and Phases 20–23 use effective transport, optical, spin or circuit models. Their presence in a Git history does not demonstrate laboratory implementation. Phase numbers identify modules, not maturity levels.

The historical material suggests a methodological contribution that can be tested and reused: a claim should be linked to a specific identity, reference, sampling unit, computational state and acceptance gate. This is stronger than a generic call for reproducibility because the preserved counterexamples are paired. The same VMC run passes internal checks but fails the He point-accuracy target; the same molecular geometry has converged methods but reference-sensitive energy comparisons; the same experimental trace changes its apparent drift under alternative startup windows. Other comparisons are explicitly unpaired: a solubility README without predictions cannot be equated with a frozen supervised benchmark, and a force-field conformer panel cannot be promoted into a GNN accuracy test. Distinguishing these cases prevents a persuasive narrative from outrunning its data.

Figures L1–L5 present the resulting historical comparisons: per-structure conformer ranges with UFF identified, VMC signed errors with saved-block standard errors and a fixed accuracy threshold, complete startup/window sensitivity for the three published traces, paired raw and Q01-centered molecular-energy contrasts, and the complete method/environment rank-reversal matrix. All five figures use the accompanying CSVs and retain negative outcomes. They describe newly organized and recomputed historical evidence; they do not add experiments, training runs or independent samples to the preserved record.


<!-- historical figures -->

![Figure L1](figures/Figure_L1_Conformer_english.png)

Figure L1. Historical conformer energy ranges. Every valid structure contributes 50 accepted optimization records; the hatched M07 bar uses UFF and all others use MMFF94. The ranges do not measure entropy or establish unique minima. [SVG](figures/Figure_L1_Conformer_english.svg).

![Figure L2](figures/Figure_L2_VMC_english.png)

Figure L2. VMC errors against the saved computed FCI CBS reference. Error bars are ±1 block standard error from 20 saved production blocks per system; they exclude reference uncertainty. The 1.6 mEh line is the unchanged point-estimate gate. Distances are in bohr. [SVG](figures/Figure_L2_VMC_english.svg).

![Figure L3](figures/Figure_L3_Trace_Windows_english.png)

Figure L3. Endpoint-change sensitivity over all 72 saved window specifications, using the same three published traces. Voltage panels share the −5 to 45 mV color scale; the source-current panel has an independent −0.16 to 0.05 scale with unresolved current normalization. The black divider marks startup exclusion of at least 10 h. Values are descriptive contrasts, not independent-replicate confidence intervals [@yang2023; @yang2026]. [SVG](figures/Figure_L3_Trace_Windows_english.svg).

![Figure L4](figures/Figure_L4_Reference_Alignment_english.png)

Figure L4. Fixed-geometry gas-phase comparison with PBE0/def2-SVP. Left: raw charge-state offsets. Right: double differences after subtracting each method’s Q01 value. The different vertical scales are intentional. PBE0 is a method comparator, and the raw offsets are not calibrated physical prediction errors. [SVG](figures/Figure_L4_Reference_Alignment_english.svg).

![Figure L5](figures/Figure_L5_Rank_Reversals_english.png)

Figure L5. All 15 method/environment comparisons from the 72 state differences. Each off-diagonal cell counts reversed orderings among the same 66 pairs of 12 molecules, not new calculations or errors against experiment. Diagonal zeros are identity comparisons. [SVG](figures/Figure_L5_Rank_Reversals_english.svg).
