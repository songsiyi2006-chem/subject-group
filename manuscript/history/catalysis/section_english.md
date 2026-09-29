# Catalysis records as a test of evidence propagation

## C.1. Scope, provenance, and the unit of evidence

The public `pincer-catmech-ai` archive provides a complementary test of the validation hierarchy developed in the main text. Whereas the H₂ study permits controlled comparisons along a well-defined coordinate, the catalysis archive contains several chemically and numerically different objects: molecular electronic-structure jobs, constructed solvent clusters, incomplete reaction paths, supervised conformer-energy labels, and synthetic kinetic fixtures. Their coexistence does not establish an integrated, experimentally validated catalytic mechanism. Our purpose is to determine which conclusions survive when the original calculation records are followed through successive acceptance criteria.

The audit fixes the public repository at commit **2edffb123791bd61acbfdeff763f603ee50c2287**, independently matched to the remote HEAD and main branch. This is later than the local historical checkout at 9aa9c08b15f2f64b3f108e4a90e0f8ec0a4352c0. The ten-commit history includes the original workflow, a native screening campaign, Phase 3 spin and microsolvation studies, Phase 4 precursor checks, and subsequent molecular preflights associated with a Cu–N–P research direction. The last description is a project context, not the elemental composition of every calculation. In particular, the latest molecular DFT model contains no Cu. Commit messages establish chronology; numerical claims below are anchored to saved results rather than inferred from those messages.

The [source catalog](sources/catalog.json) identifies each source by a C or N identifier and records its fixed commit, repository path, Git blob identifier, original-byte SHA-256, retained excerpt hash, and JSON-pointer or CSV-row selector. C01–C24 cover the curated data and attribution records; N01–N10 identify native logs used for independent energy checks. All numbers in the tables refer to these sources, and the six `posthoc_*.csv` files attach evidence identifiers and formulas to each derived row. Ten final electronic energies were re-extracted from native logs: the eight converged Phase 3 states and the two completed later molecular DFT states. All matched their corresponding machine-readable records at the archived precision, with a maximum absolute discrepancy of 0.0 Hartree. This establishes transcription consistency, not independent electronic-structure accuracy. No new quantum-chemistry job was run for the present audit.

Counts retain their original unit. A target slot is not a launched calculation; an SCF solution is not necessarily a minimum; a conformer pair is not an independent catalyst; and an input prepared for a remote computer is not a submitted job. The Phase 3 readiness table contains 24 designs and 40 temperature/base combinations per design, hence 960 rows, but records zero accepted physical predictions for the proposed 18-step kinetic mechanism [C01]. Its baseline of 383.15 K and total added tBuOK equivalent of 0.05 are scenario inputs. Total base does not determine the free-anion activity, and the alcohol activity was not measured. Multiplication of a design grid therefore measures bookkeeping scope rather than the number of experimentally constrained predictions.

## C.2. Spin screening: a converged state is not a crossing

The principal vertical spin matrix used PBE/STO-3G, a 35 × 110 integration grid, SOSCF, electronic and density convergence targets of 10⁻⁸ and 10⁻⁶, and nominal multiplicities of 1, 3, and 5 for 18 designs [C03]. Calculations were performed on supplied geometries; identity and connectivity checks do not by themselves validate the proposed active species. The archived record contains 54 target slots: 48 launched physical state attempts, 40 timeouts, six missing source geometries, and eight converged states. These counts describe the specified matrix rather than the sum of every exploratory protocol in the repository.

**Table C1. Acceptance counts in the archived vertical spin matrix.** Sources: C02–C06; exact selectors and arithmetic are in [posthoc_spin_gates.csv](sources/posthoc_spin_gates.csv). Rows are different acceptance stages, not disjoint samples.

| Record or acceptance stage | Count |
|---|---:|
| Target state slots | 54 |
| Physical state attempts | 48 |
| SCF-converged states | 8 |
| Spin- and identity-eligible states | 6 |
| Raw computable vertical gaps | 1 |
| Eligible diagnostic gaps | 0 |
| Accepted MECPs | 0 |

![Figure C1. Spin-screen acceptance counts.](figures/Figure_C1_Spin_Gates_english.png)

**Figure C1.** Counts at successive stages of the fixed-commit spin audit. The raw-gap row counts an energy difference, whereas preceding rows count states; the plot is an evidence progression and must not be interpreted as independent sampling or a statistical survival curve. [Vector figure](figures/Figure_C1_Spin_Gates_english.svg).

The distinction between convergence and a usable spin pair is visible in the expectation values of S². For Fe_bipyridine_pnnoh_iPr, the converged triplet had ⟨S²⟩ = 2.661445763547917, giving a contamination diagnostic of 0.661445763547917 relative to S(S+1) = 2. The Fe_macho_pnp_iPr quintet had ⟨S²⟩ = 6.166902016475518, or 0.166902016475518 above the ideal value of 6. Both were flagged by the archived screening rule. The other six converged states passed the recorded spin and identity checks; this is a necessary numerical filter, not proof that their electronic states represent the experimental catalyst [C03–C04].

![Figure C2. Spin contamination in converged states.](figures/Figure_C2_Spin_Quality_english.png)

**Figure C2.** Recorded ⟨S²⟩ − S(S+1) for the eight converged states, with archived spin-screen failures in orange. M denotes multiplicity. Values near zero include roundoff-level negative singlet values and are displayed as approximately zero. [Vector figure](figures/Figure_C2_Spin_Quality_english.svg).

Only one raw vertical difference can be formed from the available matrix: the Fe_bipyridine_pnnoh_iPr triplet-minus-singlet value is −0.0314591288624797 Hartree [C05]. Because its triplet fails the spin screen, this is not an accepted spin-ordering result or a crossing energy. The downstream MECP record explicitly reports `no_eligible_pair`, zero quantum pair evaluations, and zero quantum state attempts [C06]. Thus no MECP optimization was launched from this matrix. The zero result is a precondition rejection, not evidence that a crossing does not exist, and certainly not a converged minimum-energy crossing point. This distinction is central to error propagation: a downstream routine cannot supply missing state-pair information merely by accepting a numerical difference as input.

## C.3. Microsolvation: electronic association and thermodynamic association differ

The microsolvation data concern bounded constructed clusters rather than an experimentally measured solvent ensemble. Nine GFN2-xTB optimizations combined three cluster sizes with three starting seeds; all converged while retaining the recorded connectivity. Calculations used xTB 6.7.1, ALPB toluene, the extreme optimization setting, accuracy 0.5, and 300 K electronic smearing [C09–C12]. One sampled lowest-energy structure per cluster size was selected for a Hessian, yielding three selected-cluster Hessians rather than a comprehensive conformer partition function. The selected structures correspond to seeds 2, 1, and 0 for one, two, and three tBuOH molecules, respectively.

**Table C2. Selected microsolvation association quantities.** All entries are in kcal mol⁻¹. C09–C12 support the values; [posthoc_solvation_selected.csv](sources/posthoc_solvation_selected.csv) records the selected rows. The final column is the frozen-fragment contribution beyond the summed pair interactions.

| tBuOH count | ΔE association | ΔG, 298.15 K, 1 M | ΔG, 383.15 K, 1 M | Beyond-pair contribution |
|---|---:|---:|---:|---:|
| 1 | −10.994093 | 0.658254 | 3.797726 | 0.000000 |
| 2 | −21.166260 | 3.586112 | 10.220408 | 1.547073 |
| 3 | −31.716522 | 4.105509 | 13.685284 | 3.558360 |

Increasingly negative electronic association energies do not imply increasingly favorable standard association free energies. At 383.15 K, the corresponding 1 M free energies are positive and rise from 3.797726 to 13.685284 kcal mol⁻¹ across the sampled sizes. The archived thermochemistry applies ALPB solvation and one 1 atm-to-1 M correction, with every participating species defined at 1 M. Actual solution populations require concentrations or activities and adequate conformational sampling; they cannot be recovered from the table alone. Neither these clusters nor their association free energies include a catalyst, an explicit tBuOK ion pair, or an accepted transition state [C12].

The beyond-pair term also requires a narrow interpretation. For the two-alcohol cluster, the frozen interaction energy is −0.9529920930613116 eV, compared with a pairwise sum of −1.0200795115757728 eV. Their difference is positive, not extra stabilization. For the three-alcohol cluster, the corresponding quantities are −1.4095253966465862 and −1.5638304890707104 eV. Because fragments are evaluated with ALPB, this decomposition includes changes in continuum cavity and reaction-field contributions. It is not a pure isolated-molecule electronic many-body interaction [C11]. A later sensitivity calculation generated 405 arithmetic rows without new quantum jobs and exactly reproduced the 27 baseline thermochemistry rows. At 383.15 K, the reported low-frequency-cutoff ranges were 0.793212, 2.298506, and 3.465158 kcal mol⁻¹ for the three cluster sizes [C18]. This sensitivity is a model-choice diagnostic, not a statistical error bar or experimental validation.

## C.4. Proton-wire paths and kinetic fixtures: retaining the negative result

The proton-wire study contains two independently initialized paths and a continuation of the first path. It uses a neutral GFN2-xTB/ALPB-toluene model without explicit potassium or a free tBuO⁻ species [C13–C14]. The last force residual is therefore meaningful only for that modeled band, and total-base activity remains unresolved. None of the saved attempts meets the 0.07 eV Å⁻¹ force target. Continuing an unsuccessful path does not produce an additional independent chemical trial.

**Table C3. Unaccepted proton-wire path records.** The energy column is the sampled band maximum above its reactant reference, not an activation barrier. C13–C14 and [posthoc_proton_wire.csv](sources/posthoc_proton_wire.csv) provide exact values and force-ratio formulas. No row is an accepted transition state.

| Path record | Last step | Force / eV Å⁻¹ | Force/target | Band maximum / eV |
|---|---:|---:|---:|---:|
| Start 1, cluster seed 02 | 100 | 0.375882 | 5.369743 | 1.795702 |
| Start 2, cluster seed 01 | 100 | 1.347282 | 19.246886 | 4.568802 |
| Start 1 continuation | 200 | 1.206630 | 17.237571 | 2.391468 |

![Figure C4. Unaccepted proton-wire paths.](figures/Figure_C4_Unaccepted_Paths_english.png)

**Figure C4.** (a) Last recorded force residual divided by the convergence target; the dashed line is one. (b) Sampled band maxima. Start 1 and Start 2 denote independent initializations, not the numerical cluster-seed identifiers. The continuation reuses Start 1. Unconverged bands cannot establish activation free energies. [Vector figure](figures/Figure_C4_Unaccepted_Paths_english.svg).

The product-graph check passed only for the second initialization. It failed for the first and its continuation, while even the graph-passing attempt retained a force residual more than nineteen times the target. The archived elapsed time was 1182.0914974212646 s, with zero accepted transition states and zero activation free energies [C14]. A visible peak on a discretized path is insufficient: endpoint identity, band convergence, stationary-point characterization, and connection to the intended reaction remain separate questions. Here the saved evidence stops before those requirements are collectively satisfied.

The associated kinetic code has a different evidential role. Its synthetic fixture evaluates 40 grids, generated from eight temperatures and five initial free-base equivalents, at 500 time points. Rates are arbitrary test inputs. The maximum analytic-versus-finite-difference Jacobian error was 9.630873876176338 × 10⁻¹¹; the largest metal and base conservation errors were 3.642919299551295 × 10⁻¹⁷ M and 4.996003610813204 × 10⁻¹⁶ M. The metal-sensitivity conservation residual was 5.72390984436566 × 10⁻¹⁷ M [C15]. These are useful differential-equation implementation checks. They do not measure catalytic turnover, establish an 18-step network, or replace absent activation free energies. The ledger explicitly retains zero physical catalyst predictions.

## C.5. A trained EGNN with a negative baseline comparison

The archived auxiliary EGNN is trained, unlike the untrained neural potential examined elsewhere in this manuscript. Its evidence is nevertheless restricted to relative conformer energies. The dataset construction accepted 610 native-converged structures and rejected 202 records, producing 540 non-self pairs in 70 composition groups spanning six families [C07–C08]. The reference within each composition group is the lowest original conformer index, selected deterministically without an energy-based choice. Native optimization convergence does not establish that every structure is a Hessian-certified minimum.

The split holds out complete backbone/substituent families across metal, state, and conformer records. Training contains 374 pairs from four families, validation contains 80 pairs from `macho_pnp|iPr`, and testing contains 86 pairs from `macho_pnp|Ph`. This is a stronger separation than randomly splitting closely related atomistic records, but it is still a small, task-specific archive and does not demonstrate unrestricted catalyst transfer. The 12,675-parameter network was trained using PyTorch 2.10.0 on a CPU with two threads. The ledger records nine completed or partial epochs, a best-validation epoch of 0, and 85.032 s elapsed time [C07].

**Table C4. Auxiliary relative-energy prediction and simple baselines.** MAEs are in eV. The equal-reference baseline predicts zero relative energy; the second baseline uses the training target mean. C07 and [posthoc_egnn_metrics.csv](sources/posthoc_egnn_metrics.csv) contain full precision and checkpoint provenance.

| Split | Pairs | EGNN MAE | Equal-reference MAE | Training-mean MAE |
|---|---:|---:|---:|---:|
| Training | 374 | 0.286664 | 0.286438 | 0.297798 |
| Validation | 80 | 0.261617 | 0.259393 | 0.337521 |
| Test | 86 | 0.195200 | 0.193999 | 0.260416 |

![Figure C3. Auxiliary EGNN and baseline errors.](figures/Figure_C3_EGNN_Baselines_english.png)

**Figure C3.** Split-specific MAEs for the actual trained auxiliary network and two baselines. The targets are relative conformer electronic energies; no activation-barrier or MECP-gap labels were available for training those outputs. Bars do not represent repeated-training uncertainty. [Vector figure](figures/Figure_C3_EGNN_Baselines_english.svg).

The network improves on the training-mean predictor but not on the equal-reference baseline in any split. Its test MAE of 0.19519976927396127 eV exceeds the equal-reference value of 0.19399893976147659 eV, corresponding to an improvement of −0.0012008295124846802 eV [C07]. The absence of a positive result should not be concealed by reporting only the weaker baseline. The retained test RMSE is 0.3845663866708711 eV. Activation-barrier and MECP-gap heads have zero observed labels and remain untrained; the 540-pair auxiliary task cannot confer validity on those heads.

Checkpoint-level symmetry checks establish a different property. Under float64 rotations, improper rotations, permutations, and batching, the recorded maximum scalar and coordinate errors were both 1.7763568394002505 × 10⁻¹⁵, in eV and Å respectively [C07]. This is strong evidence for the tested transformations of that implementation, but invariant or equivariant errors can coexist with poor chemical prediction. The result illustrates why architecture verification, label validity, baseline comparison, and external validation must remain separate. No retrospective claim of a fully test-blind development history is made from the saved split labels alone.

## C.6. Later native calculations narrow the claim rather than completing the mechanism

Phase 4 added nine local quantum calls: three Psi4 jobs and six xTB calls. The two converged Psi4 jobs did not produce a usable same-method spin pair [C16–C17]. Their gradients also preclude interpreting SCF convergence as geometry stationarity.

**Table C5. Same-geometry Phase 4 spin pilot.** The molecular target is Fe_bipyridine_pnnoh_iPr; all calculations use PBE. Gradient maxima are in Hartree bohr⁻¹. A dash means unavailable, not zero. C17 and [posthoc_phase4_jobs.csv](sources/posthoc_phase4_jobs.csv) preserve energies, timings, and screening decisions.

| Basis; multiplicity | Outcome | ⟨S²⟩ | Maximum gradient | Spin screen |
|---|---|---:|---:|---|
| STO-3G; 5 | SCF converged | 6.193760 | 0.201560 | Fail |
| def2-SVP; 1 | SCF converged | approximately 0 | 0.035808 | Pass |
| def2-SVP; 5 | Timeout | — | — | Unavailable |

The quintet recovery energy with STO-3G was −2533.7232403662692 Hartree, whereas the def2-SVP singlet energy was −2562.3192161270376 Hartree. Their difference is not a spin gap because the basis sets differ. The def2-SVP quintet timed out after 450.0779999999795 s. None of these jobs establishes chemical validity of the proposed catalytic species [C17].

The six xTB calls had a distinct purpose: optimization, fresh gradient, and full Hessian for each of two 62-atom crystallographic precursor molecules. Both had 180 internal modes and no imaginary frequencies, with lowest frequencies of 33.46886841017446 and 33.81102119885175 cm⁻¹. Fresh maximum forces were 0.0007802942092995298 and 0.0007509869895041861 eV Å⁻¹. Their mapped structural RMSD of 0.00354291 Å does not establish distinct basins [C19]. These are useful model-level minimum checks for precursor structures, not validation of an active intermediate or spin ground state. Six remote-computing inputs were prepared and zero were submitted [C16]. The crystallographic provenance and its upstream CC BY-NC 4.0 notice are retained as C24; the audit copies the notice, not the complete external structural dataset.

The later molecular preflight concerns independently generated coordinates for 3-methylquinazolin-4(3H)-one, C₉H₈N₂O: 20 atoms, 84 neutral electrons, and no Cu [C22]. Twenty completed native xTB calls include a useful stationary-point counterexample. An initially optimized neutral structure had a significant imaginary frequency of −50.49 cm⁻¹. A rotor repair followed by very tight optimization produced a lowest positive frequency of 78.64 cm⁻¹ among 54 vibrational modes [C20]. The accepted result is an xTB-model minimum; it is not a DFT minimum or an experimental structure.

![Figure C5. Molecular preflight quality checks.](figures/Figure_C5_Molecular_Preflight_english.png)

**Figure C5.** (a) Lowest xTB vibrational frequency before and after the rotor repair. (b) Later fixed-nuclei DFT job counts, with preflight rejections excluded from native launches. The molecular model contains no Cu, and zero complete neutral/anion pairs are available. [Vector figure](figures/Figure_C5_Molecular_Preflight_english.svg).

At fixed nuclear geometry, the GFN1/GFN2 differences in anion-minus-neutral charging energies were 0.7816478999991645 eV in the gas phase and 0.8319797664252286 eV with equilibrium ALPB DMF, exceeding the archived 0.2 eV sensitivity threshold [C20]. These quantities cannot be relabeled as nonequilibrium vertical electron affinities, electrode potentials, or solution reaction free energies. The ensuing native DFT pilot used a diffuse basis but did not complete any charge pair.

**Table C6. Later fixed-nuclei molecular DFT launches.** All jobs use def2-SVPD. Sources C21–C23 and [posthoc_dft_jobs.csv](sources/posthoc_dft_jobs.csv) distinguish four native launches from two additional preflight rejections. A dash denotes an unavailable converged energy.

| Functional; charge | Outcome | Energy / Hartree | Wall time / s |
|---|---|---:|---:|
| PBE0; 0 | SCF completed | −531.566397414 | 108.659154 |
| PBE0; −1 | Timeout | — | 300.061065 |
| B3LYP; 0 | SCF completed | −532.173578168 | 117.232500 |
| B3LYP; −1 | Timeout | — | 300.134975 |

The neutral calculations use 307 basis functions. Neither the complete method spread for charging nor a converged electron-addition energy can be obtained from the two neutral energies; specifically, subtracting neutral energies from different functionals does not yield an electron affinity. SCF stability was not verified; the protocol has no solvent, DFT geometry optimization, or DFT Hessian [C21]. In this archive, increased computational sophistication therefore refines the location of missing evidence rather than completing a Cu-catalytic mechanism.

## C.7. Reproduction and implications for the main argument

The audit is reproduced from an existing Git object store, without checking out the complete large repository or rerunning chemistry. With a preconfigured Python interpreter, `python sources/reproduce_audit.py --repository <local-git-object-store>` regenerates the retained extracts and arithmetic tables at the pinned commit. `python plotting.py` regenerates the figures from those extracts; NumPy, Matplotlib, and the documented fonts must already be available. Paths are relative to this section's directory. The [audit record](audit.json), [native-energy comparison](sources/native_energy_checks.json), [commit history](sources/commit_history.json), and [figure manifest](figures/manifest.json) separate numerical provenance from rendering. PNGs receive AI visual inspection; SVGs receive structural inspection, which is not an independent rendered-image review.

The results supply several concrete stopping rules for evidence propagation. A converged but spin-contaminated state cannot automatically label a crossing problem. Attractive cluster electronic energies cannot replace concentration-dependent association thermodynamics. An unconverged band cannot label a reaction barrier. A symmetry-correct neural model can fail a simple baseline, and training one auxiliary head does not train another physical observable. Finally, a molecular preflight under a catalysis project title does not become a metal-site calculation without the metal and its environment. These are model-specific findings supported by archived native data, not assertions that the intended catalytic chemistry is impossible. Together they explain why the main paper treats acceptance criteria and observable definitions as part of the scientific result rather than as administrative annotations.
