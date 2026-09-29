## Pharmacological extensions: observable-specific evidence across scoring, kinetics and molecular simulation

### Historical scope and the meaning of a recorded result

The pharmacology history provides an independent setting in which computational execution and the intended scientific observable can diverge. We examined the complete nine-commit history of `ai4pharm-lead-developability-suite`, fixed at `a7a2a4df501d5dad030dc488691b02c73fdb423f`. The selected cases concern molecular descriptors, ternary binding, covalent kinetics, an atomistic ternary-complex pilot and structure-guided search. They were chosen for inspectable numerical records and their connection to the manuscript's validation argument, rather than as evidence of drug discovery. No new molecular dynamics, electronic-structure calculation, docking or model training was performed for this integration.

Historical identity was checked at the git-object level. The selected Task 1, Task 2 and Task 3 result tables are byte-identical to their original publication commits, despite subsequent directory reorganization. Their copies in the omnibus recomputation directory are also byte-identical. Consequently, these copies are not pooled as additional molecules, independent chemical replicates or new experimental evidence. The [audit](audit.json) links every reported numerical group to its commit, source path, SHA-256 and field names; [source metadata](sources/source_manifest.json) distinguish archived execution provenance from later presentation changes. The MIT license accompanies the copied pharmacology material, while third-party structures retain their attribution.

### Multiparameter scores preserve useful arithmetic but not clinical validity

The frozen developability panel contains thirty parent structures, each with a PubChem identifier and structural provenance. Three purposive strata contain ten compounds each. Those identities are public structure records, not thirty measured ADMET labels. The CNS multiparameter optimization score combines six desirability transformations; the historical implementation substitutes toolkit descriptors and chemical-class ionization assumptions. The distinction matters because exact reproduction of a scoring equation does not validate its inputs for a particular ionization state or assay. The original CNS-MPO method provides the methodological reference, whereas the repository's oral score is its own equal-weight ranking construction. [Wager et al., 2010](https://doi.org/10.1021/cn100008c)

**Table P1. Frozen developability strata. Scores are descriptive proxies; conformer counts are sampling outputs.**

| Stratum | Molecules | Median CNS-MPO | Median oral proxy | Embedded / converged conformers | Molecules with partial convergence |
|---|---:|---:|---:|---:|---:|
| Oral references | 10 | 4.989850 | 81.049384 | 48 / 48 | 0 |
| Toxicity-related comparators | 10 | 3.647260 | 75.028315 | 54 / 54 | 0 |
| bRo5 modalities | 10 | 1.617406 | 43.759364 | 60 / 54 | 4 |

All six-term CNS sums, the ESOL equation and the equal-weight oral score were independently recalculated from full-precision stored fields. Agreement was within 10⁻¹². This is an arithmetic check, not an ADMET accuracy result. In particular, the hERG, Caco-2 and absorption outputs use unfitted formulas; a bounded hERG score is not an event probability. Clinical stratum labels do not enter the scoring formula, but a hand-selected clinical panel still cannot estimate prospective discrimination. Differences between its medians are therefore descriptive, without significance or generalization claims.

Geometry introduces a separate limitation. Of the 162 embedded conformers, 156 converged; all six failures occur within four bRo5 molecules. The stored absolute chameleonic hydrogen-bond index is missing for all thirty structures. That is an appropriate non-identification result: a small, equally weighted gas-phase ensemble does not measure solvent-dependent exposure. The plotted pKa envelope varies assumptions and must be read as a scenario range, not a confidence interval. This case shows why provenance suffixes such as “proxy,” “approximate” and “assumed” need to remain attached to downstream charts.

![Developability scores and sampling limits](figures/Figure_P1_english.png)

**Figure P1.** Individual approximate CNS-MPO scores and their assumption-based pKa ranges; group medians and conformer convergence are shown separately. The selected strata and uncalibrated score do not establish predictive accuracy.

### Ternary equilibrium: a validated equation does not identify productive degradation

The ternary-complex model is especially useful because multiple observables can be separated within the same exact mass-action system. Its archived grid has 724 rows: four cooperativity scenarios and 181 total degrader concentrations per scenario. With E₀ = T₀ = 100 nM, Kᴅ,E = 10 nM and Kᴅ,T = 100 nM, the peak total concentration remains approximately 131.622777 nM while cooperativity changes the peak complex concentration and the width of the high-occupancy window. The alternative supplied peak expression gives 148.323970 nM, a 12.688680% excess for this parameter set. Its disagreement is a formula-level issue, not numerical solver noise.

**Table P2. Ternary-equilibrium sensitivity at fixed affinities and total protein concentrations.**

| Cooperativity α | Exact peak total degrader / nM | Peak ternary complex / nM | 80%-of-peak window / nM |
|---:|---:|---:|---:|
| 0.01 | 131.622777 | 0.570646 | 70.727187–236.327186 |
| 1 | 131.622777 | 29.053551 | 68.346071–275.136772 |
| 10 | 131.622777 | 66.147685 | 69.289276–439.603083 |
| 100 | 131.622777 | 87.675477 | 73.626624–1287.533850 |

Independent substitution of all stored species into the three mass balances and three equilibrium relations passed a 10⁻⁹ relative threshold. Yet neither a high peak nor a broad concentration window identifies ubiquitination-competent geometry, cellular permeability or degradation rate. The repository's HiBiT and western-blot-like readouts are explicitly synthetic. They cannot be relabeled as biological replication simply because their curves are compatible with the chosen kinetic equations.

The associated linker calculation further separates a geometric distribution from biological function. Each capped fragment retains one hundred ETKDG/MMFF draws, including duplicates; these are not one hundred independent conformational basins. The extremely small rigid-linker spread indicates repeated convergence near the same geometry. The assumed 8–12 Å compatibility interval contains no experimentally established exit-vector information, and its occupancy is not a binding-entropy estimate.

**Table P3. Unweighted linker-fragment sampling; the interval is an illustrative geometric criterion.**

| Capped linker | Draws / converged | Mean endpoint distance / Å | Sample SD / Å | Fraction within 8–12 Å |
|---|---:|---:|---:|---:|
| Flexible PEG | 100 / 100 | 9.671507 | 0.957658 | 0.96 |
| Rigid alkynyl | 100 / 100 | 9.604954 | 3.909640 × 10⁻⁷ | 1.00 |

![Ternary binding and linker distributions](figures/Figure_P2_english.png)

**Figure P2.** Mass-action hook curves and sampled linker endpoint distances. Cooperativity alters occupancy while these scenarios retain the same optimum total concentration. Linker fragments are geometric proxies, not simulated protein-bound degraders.

### Covalent kinetics: changing the approximation changes the claimed efficiency

Eight hypothetical warheads combine executed extended-Hückel orbital descriptors with assumed kinetic constants. These are two distinct evidence types. Their displayed kinetic provenance explicitly states that no experimental calibration was performed; the descriptor-derived barrier is not a located transition-state barrier. Optional xTB capability is not evidence that xTB was executed for this panel.

The choice of efficiency denominator is consequential even before biochemical validation. The dissociation constant is Kᴅ = kₒff/kₒn, whereas the forward quasi-steady-state denominator is Kᵢ = (kₒff + kᵢnact)/kₒn. Thus the rapid-equilibrium efficiency exceeds the quasi-steady-state efficiency by the relative factor kᵢnact/kₒff. These expressions are familiar approximations; the contribution here is preserving their distinct meanings in the archived comparison, rather than proposing new kinetics.

**Table P4. Assumed kinetic scenarios and executed EHT localization descriptors. Efficiencies are in M⁻¹ s⁻¹.**

| Warhead ID | kᵢnact/Kᴅ | kᵢnact/Kᵢ | Relative excess / % | Reaction-centre LUMO population proxy |
|---|---:|---:|---:|---:|
| W01 | 30000.000 | 29126.214 | 3.000 | 0.351241 |
| W02 | 20000.000 | 19512.195 | 2.500 | 0.348904 |
| W03 | 60000.000 | 54545.455 | 10.000 | 0.460621 |
| W04 | 45000.000 | 41284.404 | 9.000 | 8.903366 × 10⁻⁹ |
| W05 | 10666.667 | 10389.610 | 2.667 | 8.093850 × 10⁻⁷ |
| W06 | 400.000 | 399.202 | 0.200 | 2.981633 × 10⁻⁵ |
| W07 | 90000.000 | 81818.182 | 10.000 | 0.470255 |
| W08 | 120000.000 | 109090.909 | 10.000 | 0.471430 |

Across two hundred concentration records, the fitted forward-only slow rates agree with an independently evaluated matrix eigenvalue to within 7.82 × 10⁻⁸ relative error. This does not establish reversible residence times: the reverse reaction is deliberately suppressed in this diagnostic. Likewise, a computed global LUMO cannot be assumed to localize on the electrophilic atom. Several population proxies in Table P4 are very small. Numerical rate agreement and the chemical interpretation of orbital localization therefore require separate checks.

![Kinetic approximation and localization](figures/Figure_P3_english.png)

**Figure P3.** Rapid-equilibrium and quasi-steady-state efficiencies, alongside reaction-centre EHT population proxies. Rate constants are hypothetical inputs; neither panel measures target potency or GSH safety.

### Atomistic and docking execution still require observable-specific acceptance

The later atomistic case starts from the public 5T35 crystal structure. Independent heavy-atom distance calculations recover 112 atom pairs and fourteen residue pairs within 4.5 Å between the selected VHL and BRD4 chains, with a minimum distance of 2.798363 Å. This is geometric analysis of an experimental structure. The original structural study also used biophysical and cellular evidence; those experiments are external literature, not measurements reproduced by this repository. [Gadd et al., 2017](https://doi.org/10.1038/nchembio.2329)

The saved OpenMM pilot has 150,838 particles but only 0.2 ps of equilibration and 1.0 ps of biased dynamics. Twenty recorded Gaussian depositions satisfy the archived well-tempered height rule; the observed distance spans only 0.024018 nm. Finite coordinates, valid checkpoint metadata and a reproducible bias are useful execution checks. They provide no evidence of recurrent basin transitions or converged free energy. Accordingly, the archived PMF-convergence flag is false and cooperativity remains null. Particle count does not compensate for inadequate sampling of the desired observable. [@openmm; @wtmeta]

The structure-guided search has an equally important control. A population history of two thousand rows corresponds to 1,034 evaluated graphs, while the matched-budget random search also evaluates 1,034 graphs. Their union contains 1,597 identities, implying an overlap of 471; population rows are not new independent molecules. The final best geometric score improves by 7.543625% over the initial population but by only 0.084750% over the single random control. The random control also has a higher selected median score. These observations do not support replicated optimization superiority.

**Table P5. Executed molecular workflows and the observables they do not yet establish.**

| Recorded comparison | Numerical result | Remaining evidence boundary |
|---|---|---|
| Atomistic pilot | 150,838 particles; 0.2 + 1.0 ps; 20 hills | No converged PMF or cooperativity |
| Best geometry score: initial / evolved / random | 19.219952 / 20.669833 / 20.652330 | One random control; arbitrary score |
| Selected median: evolved / random | 14.515998 / 15.329122 | No replicated superiority estimate |
| Twelve candidate Vina scores | −6.368 to −4.405 kcal mol⁻¹ | Empirical docking scores, not measured affinity |
| Two X77 redocking heavy-atom RMSDs | 1.242364 / 1.054641 Å | Pose recovery does not calibrate affinity |

Twelve candidate docking records and two redocking controls confirm that docking was executed. However, recovering a known pose tests geometry under that receptor-preparation protocol, not binding thermodynamics across new ligands. The repository appropriately leaves predicted nanomolar affinity and confirmed novelty absent. The score cannot be converted into measured affinity by changing its label. [@vina]

![Molecular pilot and search controls](figures/Figure_P4_english.png)

**Figure P4.** The brief atomistic distance trace, equal-budget geometric-search comparison and the selected docking-score subset. Axes represent different observables. Executed dynamics and docking do not establish converged free energies, clinical effects or experimentally confirmed leads.

### Literature history contributes context, not additional calculations

`Paper-Analysis-ChemMLLM`, fixed at `c1d309a0c2f876ccc8b4d2ce753ece971d58fe97`, has six commits and twenty-seven tracked files. Its tree contains reading summaries and document artifacts, with no scientific driver, checkpoint or primary result dataset. The public article on ChemMLLM can support discussion of multimodal chemical interfaces; the copper-catalysis article and single-atom electrosynthesis review provide distinct experimental and review contexts. Their DOI records were checked against original publisher pages. None supplies an executed result of the present work. [@chem_mllm; @cu_data_driven; @single_atom_review]

The audit retains minimal repository provenance and original bibliographic citations, without copying journal PDFs or published figures. No open-license declaration was found in the reading repository, and no journal-quartile claim is inferred from reputation. Its multiple educational editions are not independent scientific replications. Together, these historical cases reinforce an observable-specific rule: structural identity, exact arithmetic, stable integration, pose recovery and biological effect occupy different evidential levels and must not be silently substituted for one another.
