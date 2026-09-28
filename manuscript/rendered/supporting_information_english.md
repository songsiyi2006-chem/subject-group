# Supporting information: layered validation and error propagation in computational chemistry workflows

This document consolidates the archived calculations accompanying the manuscript. Directory names identify successive software case studies, not validated levels of research or industrial readiness. No new scientific calculation was run for this consolidation. Numerical entries below are taken from retained machine-readable results and their methodological records; source-code assertions are distinguished from reviewed results.

## S1. Evidence classes, provenance and calculation ledger

Five evidence classes are kept separate: **public experimental labels**, used only in the FreeSolv hydration benchmark; **executed electronic-structure calculations**, comprising GFN2-xTB and Psi4 H/H₂ work; **force-field or analytic-model calculations**, including conformers, transport, stochastic cycles and Morse nuclei; **synthetic or supplied values**, including reaction yields, analytical examples and adsorption curves; and **candidate inputs**, including unexecuted Gaussian/ORCA files and proposed laboratory procedures. Public experimental labels do not mean experiments were performed in this study. No oxidation-potential, catalyst-activity, preparative-electrolysis or spectroscopy measurement was acquired.

The reproducibility chain links supplied source, preserved execution or failure, explicit compatibility changes, reviewed implementation, raw arrays, summaries and publication checks. API repairs are not assumed to preserve every downstream numerical value: enlarging an embedding table can change random-number consumption. Conversely, successful execution, small matrix residuals or invariant scalar energies and covariant forces do not supply missing chemical calibration. Peer-agent cross-checks are computational assistance, not independent human peer review.

**Table S1. Selected calculation units.** Main counts exclude pilots and tests unless the stage is explicitly named. Rows are not interchangeable units and should not be added into one performance score.

| Stage and unit | Count |
|---|---:|
| Initial benchmark: paired synthetic BO observations | 1800 |
| Initial molecular study: GFN2-xTB jobs | 36 |
| Production audit: paired BO observations | 1680 |
| Production target: GFN2-xTB jobs | 3 |
| Research molecular audit: GFN2-xTB jobs | 16 |
| Toolkit: synthetic peak-recovery cases | 80 |
| ElectraTwin original reviewed workflows: PDE solves | 121 |
| ElectraTwin extensions: PDE solves | 1892 |
| ElectraTwin main extension: cached selection uses | 4800 |
| ElectroGraph main: neural fits / SSA trajectories | 4 / 273 |
| ElectroGraph main: SSA events | 9950504 |
| SynthaPore main: MD trajectories / integration steps | 44 / 350000 |
| QuantumEqui original main: SCF jobs / analytic gradients | 55 / 42 |
| QuantumEqui original pilot: SCF jobs / analytic gradients | 3 / 1 |
| QuantumEqui correlation main: energy drivers | 183 |
| QuantumEqui correlation pilot / recovery: energy drivers | 14 / 6 |
| QuantumEqui correlation: pre-driver option rejections | 6 |
| QuantumEqui learning extension: main / pilot fits | 6 / 2 |
| QuantumEqui nuclear extension: main / pilot / tail eigenproblems | 38 / 6 / 7 |

There are 55 GFN2-xTB jobs across three distinct stages, with repeated chemical identities permitted; this is not 55 unique molecules. The original 58 Psi4 jobs and 203 subsequent energy-driver calls are separately retained. FCI calls include their underlying SCF phase, which is not counted again. Cache queries, nominal neural energies and force-field minimizations are not new quantum calculations.

## S2. Initial five-topic benchmark and molecular extension

The initial study tests whether attractive optimization, classification and catalyst-ranking outputs survive elementary controls. Its five source modules use generated targets: a noisy quadratic electrosynthesis yield surface, fabricated redox regression labels, empirical catalyst scores, threshold-based annulation labels and an ODE flow proxy. The source BO uses eight initial and twelve sequential observations with a Matérn GP/UCB acquisition. Its maximum observed yield is 95.0447486%, whereas the noiseless response there is 91.7937080% and the generator maximum is 92%. The positive discrepancy is sampled noise, not chemical improvement. Redox training R²=0.998461 and classification training accuracy=1 do not establish held-out performance; the reported annulation feasibility 0.892 is a constant.

The paired extension uses 30 seeds, three strategies and 20 observations per campaign. Mean recommendation regrets are 3.555742, 1.356960 and 1.320765 synthetic percentage points for random, original GP and scaled GP, respectively. A separate 5000-row regression test across four training sizes includes 360 fits: at training size 60, boosting test R²=0.928671 versus linear regression 0.992062. Repeated catalyst-score draws rank the original Co–N₄ winner first only 5.841% of the time at the assigned noise scale. Classification balanced accuracy falls from 0.912908 overall to 0.727614 near the generator boundary. These are robustness results for the generators.

The dimensionally explicit flow extension assumes A→P→D, first-order rates 0.015 and 0.002 s⁻¹, 1 mL reactor, 0.10 M feed, 200 g mol⁻¹ product and 0.200 A available current. Analytical and numerical species fractions agree within 6.41×10⁻¹². Required current includes both reaction extents, so some high-flow kinetic outputs are charge-infeasible. Among 191 grid points, the conditional qualifying maximum is 0.55 mL min⁻¹ and 464.00 g L⁻¹ h⁻¹. Its assumed-rate scenario intervals are not measured process uncertainty.

Actual molecular work covers benzene, pyridine, anisole, indole, N-methylindole and benzofuran. Each has a neutral GFN2-xTB/ALPB(acetonitrile) optimization, fixed-neutral-nuclei cation/anion calculations and three gas-phase single points. Neutral singlets and ±1 doublets give 36 jobs. These are semiempirical calculations, not DFT. ALPB-minus-gas sensitivity and charge-response sums are retained; separate equilibrium solvent responses prevent interpreting differences as rigorously nonequilibrium vertical ionizations or reference-electrode potentials. Frequencies, charged relaxation and chemical site validation are absent. Sources: [original outputs and audit](../../data/original/benchmark_audit.json), [extended summary](../../results/extended/summary.json), [molecular energies and identities](../../results/molecular/summary.json).

## S3. Production case: representation, geometry and reduced flow

This case investigates whether physically named features improve a synthetic optimizer, whether a geometry/charge score identifies reactive C–H sites, and whether fitted catalyst scores transfer. The original run fails at a SASA API; the preserved repair additionally specifies positive atomic radii, an explicit model-convention change. The optimizer searches 288 discrete conditions, using six initial and eight acquired values. Thirty paired seeds compare four strategies; mean regrets are 2.640387 for random, 3.020240 for the source GP, 4.808777 for scaled one-hot GP and 5.154517 for scaled physical-feature GP. Normalization and physically named descriptors do not confer superiority on this generator.

For 5-methoxy-2-phenylindole, 32 ETKDG/MMFF conformers and six probe/algorithm settings test heuristic stability. C2 and C5 are substituted and lack C–H bonds. The highest score belongs to pendant phenyl sites, not C3; Gasteiger charge is not radical-cation spin density. Three target GFN2-xTB/ALPB jobs give fixed-nuclei removal/addition differences of 10.514086/6.002145 eV, without calibrated electrochemical meaning. Sixteen synthetic catalyst labels produce an in-sample R² of one, but metal-group holdout RMSE is 2.827171 for Extra Trees versus 1.263652 for ridge, in the assigned score units.

The flow calculation is a linear axial ODE at fixed effective overpotential, not a resolved electrohydrodynamic PDE. Analytical/DOP853 comparisons and integrated current check its equations. Conversion falls from 77.55% to 23.80% while reactant-disappearance STY rises from 465.30 to 1713.41 mmol L⁻¹ h⁻¹ across the source flow range. Higher flow increases the assigned mass-transfer coefficient while shortening residence time; it does not simply worsen mass transfer. Product selectivity and calibrated kinetics remain missing. Data and patch: [audit summary](../../production/results/audit/audit_summary.json), [repair](../../production/source/repair.patch), [target jobs](../../production/results/target_xtb).

## S4. Research case: Pareto preferences, substrate labels and particle diffusion

The source enumerates 135 synthetic conditions with 12 nondominated records; it performs no GP fit. A “green score” selects a condition with yield 43.47%, FE 74.49% and SEC 0.779 kWh kg⁻¹, all generator outputs. Auditing 231 weight combinations yields ten distinct selections. The added GP is retrospective, using five contiguous current-band folds of 108 training/27 held-out rows, not sequential BO. Macro yield RMSE is 1.506982 percentage points versus 12.663336 for the training-mean baseline. Selection depends on preferences, area and cost assumptions.

Eight specified substrates receive 64 converged MMFF conformers, eight neutral optimizations and eight fixed-nuclei cation GFN2-xTB/ALPB calculations. The source yield labels remain random; its charge-based score does not depend on 3D coordinates, and one nominal “aromatic” selected site belongs to caffeine's nonaromatic atom. Structure calculations do not convert those labels into experimental scope data.

The prescribed first-order spherical pore model uses φ=R√(kv/Deff) and η=3(φ cothφ−1)/φ². Forty mesh, 50 external-film and 60 parameter cases test finite-volume and analytic limits. At 100 μm, the microporous model gives η=0.488961, contradicting the hardcoded statement that radii above 50 μm yield η<0.40. The actual model threshold is 129.379256 μm. Maximum finite-volume/analytic discrepancy is 2.93519×10⁻⁵. Assigned diffusivities and kinetics have no material-specific calibration. Source tables: [complete audit](../../research/results/audit/audit_summary.json), [molecular records](../../research/results/molecules), [original results](../../research/results/original/research_grade_results.json).

## S5. Closed-loop case: accounting, candidate inputs and feedback

This case checks an experimental plan without treating it as an executed experiment. At 0.2 mmol in 6 mL, substrate concentration is 0.0333333 M. A 12.5 mA cm⁻² setting on 1.5 cm² gives 18.75 mA; 2.2 F mol⁻¹ requires 42.4535452 C and 37.7364846 min. For a stipulated two-electron product, complete substrate conversion would cap FE at 90.9091%. The plan still lacks a fully defined product, coupling stoichiometry and verified quench/analysis procedure.

Six mock feedback labels support only a diagnostic study. The source prediction-label MAE of 2.883333 percentage points is not cross-validation. Twenty-four leave-one-out predictions compare four methods: scaled maximum-likelihood GP MAE is 15.255953 versus 13.980000 for the mean baseline. Ninety candidate conditions give 270 GP predictions; true EI is distinguished from the source's UCB expression. Nine acquisition scenarios and 200000 normal draws inspect formula and assumption sensitivity, without validating a prospective recommendation.

Four Gaussian and eight ORCA files cover two molecules and neutral/cation states, reusing two previously calculated geometries. Reviewed inputs correct an indoline/target identity mismatch and CPCM/SMD labeling. No engine or scheduler executed these inputs; input generation is not DFT evidence. The feedback schema separates mock, proposed and experimentally supported records. Sources: [arithmetic/model audit](../../closed_loop/results/audit/audit_summary.json), [feedback analysis](../../closed_loop/results/audit/active_learning.json), [input manifest and execution boundary](../../closed_loop/inputs/reviewed/README.md).

## S6. Analytical toolkit: measurement arithmetic and model misspecification

Pre-integrated HPLC/NMR examples test analytical accounting. The supplied HPLC inputs imply 1292.00% yield; mass-based qNMR gives 91.325287%, whereas the original rounded standard amount gives 91.392%. No raw chromatogram, FID, calibrated assignment or independently measured standard purity supports these as measurements. Every one of 100000 assumed-error HPLC draws remains above 100%. The parsers reject missing provenance, duplicate identities, unsupported units and invalid denominators, but even a correctly hashed file does not establish assay validity; `measurement_validated` remains false.

Twenty-one synthetic calibration rows span seven concentration levels. Eighty known-shape peak-recovery cases use two Gaussian peaks and a linear baseline. Mean absolute target-area error rises from 0.027363% with correct widths to 14.614337% with a misspecified width. Small conditional regression residuals cannot exclude a wrong peak model. Six ideal pulse-relaxation scenarios and twelve purity/sample-fraction cases likewise explore assumptions rather than real instrument performance.

Exact dominance yields 23 Pareto points among 80 synthetic conditions, whereas the source threshold rule selects zero. The six substrate prediction/mock-label pairs and five constant energy levels provide neither chemical generalization nor DFT pathways. The accompanying grant text, timeline and CNY 10000 budget are planning documents, not funding, supervisory consent or laboratory access. Source: [analytical, Pareto and planning audit](../../toolkit/results/audit/audit_summary.json), [CSV data dictionary](../../toolkit/data/README.md), [all recovery cases](../../toolkit/results/audit/deconvolution_recovery.csv).

## S7. ElectraTwin: conservation, reaction networks and sequential decisions

The reviewed transport model solves cell-centred finite-volume convection/diffusion for abstract A/P species, with axial upwinding, two-direction diffusion and a half-cell Robin wall flux. Original-source current and material-flow accounting differ by 15.0351%; eight of twelve source FE values exceed 100% before clipping. The compatibility change only replaces removed NumPy `trapz` with `trapezoid`. Reviewed baseline conversion/current are 58.3872859%/42.2513750 mA. Thirty-nine verification solves and 82 control-workflow solves give 121 PDE solves. Single-reaction desired-product FE is structurally 100%, so the useful optimization objectives are STY and electrical SEC, conditional on assumed molecular mass and voltage.

The initial eight-seed comparison gives five MC-EHVI wins and three random wins; its mean paired hypervolume difference is 0.030407 with descriptive interval [−0.008497,0.069312]. Twenty-seven hydraulic and nine electrical-heating scenarios use selected fluid properties, not CFD or a thermal PDE. The SCPI state machine operates only in memory, has zero physical connections and terminates with output off.

The extension adds explicit A→P, A→B and P→D wall kinetics. The four-species flux matrix has zero column sums, conserving abstract molecular amount; electron counters distinguish gross formation from intact net P. Sixty-two network solves plus one comparison give baseline conversion 58.4087%, net P yield 11.5438% and net P FE 11.3871%. The uncertainty and optimization modules retain the older single-reaction model: these side-reaction losses cannot be appended to their predictions.

Five chosen independent parameter distributions generate 1792 Sobol-design solves, 16 mesh checks and 21 derivative solves. Jansen estimates use N=256 with prefix checks; some first-order indices exceed total effects and first-order sums exceed one. These estimator defects remain untruncated. Five hundred paired-row resamples are exploratory diagnostics, not validated QMC confidence intervals. The 64-seed, three-policy, 25-selection extension uses 4800 cached values, not new PDE solves. Mean final hypervolume fractions are 0.992024, 0.982781 and 0.989218 for GP, random and maximin. Across 26 holdout splits, quadratic regression outperforms the fixed GP on both blocked-input groups. Sources: [transport verification](../../electratwin/results/transport/verification.json), [network](../../electratwin/results/reaction_network/study_summary.json), [uncertainty](../../electratwin/results/uncertainty/summary.json), [benchmark](../../electratwin/results/benchmark_extension/summary.json).

## S8. ElectroGraph: experimental hydration labels, conformers and stochastic cycles

The hydration-data source and direct stochastic-simulation method are attributed to FreeSolv and Gillespie, respectively [13, 17]. These references identify the data and algorithm, not the validity of the assumed electrochemical mechanism.

The supplied network is untrained, its exploration score is random, and source SASA radii are zero. Four explicit installed-API repairs permit execution without validating those assumptions. The reviewed work below changes both implementation and evidence source.

**Supervised data.** The 642-row FreeSolv/SAMPL snapshot matches official v0.52 structures and experimental values within 10⁻⁸ kcal mol⁻¹. A label-blind hash subset contains 256 molecules, partitioned 154/51/51 into training/validation/test. Ring structures share nonchiral Murcko groups; acyclic structures use complete canonical identities, so close acyclic analogues may cross partitions. A two-layer edge-conditioned GRU MPNN has 29369 parameters. Three seeds run 60 epochs each, plus a shuffled-training-label control; scaling and checkpoint selection exclude test labels. MPNN test RMSEs are 1.897239, 1.746786 and 1.249601 kcal mol⁻¹, versus descriptor-ridge 1.653425. Two seeds lose to the simpler baseline. These public hydration labels do not support oxidation potentials or atom-reactivity heads. Experimental uncertainty fields are retained but not modeled in the fitting loss. [Data, attribution and CC-BY 4.0 terms](../../electrograph/data/learning/README.md); [learning summary](../../electrograph/results/learning/summary.json).

Development isolation is limited: the historical pilot outputs 1024 predictions over all 256 molecules, including test structures. Main fitting and checkpoint selection remain partitioned, but the archive does not establish end-to-end blindness to test performance during development. The reported split therefore supports a bounded retrospective comparison, not an untouched prospective test.

**Molecular structure.** Seven molecules receive 24 ETKDGv3 requests for each of three seeds: 504 requests yield 120 pruned geometries, all MMFF94-converged, and 36 pooled distinct minima after 0.35 Å heavy-atom RMSD deduplication. Weights exp[−(E−Emin)/RT] at 250/298.15/350 K assume unit degeneracy and are force-field minima weights, not solution free energies. Positive vdW radii, a 1.4 Å probe and 24 orientations give 2880 SASA evaluations. For the confirmed C5-methoxy/C2-phenyl target, mass-weighted Rg=3.608485 Å versus 4.094013 Å under the unweighted definition. The source unmapped reaction loses C₁H₄; thioanisole and thiophenol are distinct. Only an explicitly mapped, balanced ethanol→acetaldehyde+H₂ educational reaction is assigned bond edits. [Structure and reaction records](../../electrograph/results/structure/summary.json).

**Kinetics.** Five irreversible, independent-site states are sampled by direct SSA and compared with CTMC transients and renewal identities. Main work comprises 273 trajectories and 9950504 events, plus 105 analytic rate perturbations; the separate pilot adds five trajectories and 169216 events. At fixed rates, TOF∞=(Σi1/ki)⁻¹. The assigned adsorption/coupling/desorption rates bound its expectation below 26.61034847 s⁻¹, contradicting a hardcoded source value of 260.2 s⁻¹. Maximum scan-mean disagreement with finite-window CTMC is 0.2043%. This verifies the specified stochastic model, without a spatial lattice, reversible thermodynamics, transport or fitted chemical rates. [Kinetic summary](../../electrograph/results/kinetics/summary.json).

## S9. SynthaPore: symmetry, sampling, paths and geometric pores

The weak bath-coupling method used in the thermostat comparison is attributed to Berendsen and colleagues [18]; the variance ratios below are results of the repository's own controlled model study.

The source requires changing a ten-entry species embedding to 30 for Cu; the extra random draws alter initialization. Its untrained 55394-parameter energy model acts on an eight-atom C₄H₂CuN fragment. Exact captured-weight auditing therefore accompanies, rather than replaces, the compatibility result.

Supervised denoising uses 368 artificial eight-node shapes split 192/48/64/64, with independent train/validation/test/warped-OOD generators. Three 13176-parameter networks train for 60 epochs each. Held-out performance improves over the fitted smoother, but seed 4441 is worse on warped shapes. The pilot already evaluates the same 368 shapes, including test/OOD, producing 1472 per-shape metric rows; validation-only main checkpoint selection does not establish end-to-end test blindness. There is no time-dependent score, reverse process or molecular diffusion generation. Original charge updates depend on the origin, and untrained forces lack chemical calibration. [Equivariance/denoising records](../../synthapore/results/equivariant/summary.json).

MD uses an analytic eight-particle harmonic cluster with 21 internal Cartesian degrees of freedom. Twenty NVE trajectories and 24 thermostat replicas total 350000 integration steps. Energy-error convergence orders span 1.96457–2.04409. At 300 K, BAOAB gives variance/canonical-reference ratio 0.982727; correctly counted Berendsen gives approximately 0.000000310355 despite a mean near 300 K. A correct mean is not canonical sampling. Eighteen reviewed analytic CI-NEB cases converge, while ten of twelve source-algorithm controls report convergence without passing the independent force threshold. The actual captured-source candidate has maximum projected final force 0.009773735 in nominal eV Å⁻¹; its reported zero barrier occurs because the endpoint is the maximum, not because a transition state was established. [Dynamics](../../synthapore/results/dynamics/summary.json); [source replay](../../synthapore/results/source_audit.json).

Pore work evaluates three explicit two-dimensional geometries on five midpoint grids, totaling 1636800 clearance evaluations and 90 probe masks. At 640² resolution, maximum area-fraction error is 0.000474912. Point-to-wall clearance histograms are not pore-size distributions, 3D volume or surface area. The source's probe argument is unused; its circular mask is not a honeycomb framework. Twenty synthetic adsorption points switch branches at relative pressure 0.35 with a 204.779367 cm³ g⁻¹ jump. The 152 BET fits test window/noise sensitivity, retaining unphysical constants and misfit. There is no atomistic framework, GCMC run or measured isotherm. [Geometry and adsorption summary](../../synthapore/results/pore/summary.json).

## S10. QuantumEqui original: electronic, path and thermal rejection tests

The reviewed path implementation uses established climbing-image NEB and improved-tangent formulations [12, 19]. The comparisons below assess the implementation on its stated potentials and do not introduce a new path-search theory.

The source executes unchanged, but its 31-function EHT overlap has rank eight. Canonical orthogonalization leaves capacity for 16 electrons, below the assumed 40, and therefore cannot repair the proposed electronic problem. Full-space residuals remain approximately 4.59/4.63 despite a small projected residual. The source 24385-parameter EGNN is untrained; equivariance and coordinate-derived forces cannot fix this reference failure.

An actual H₂/H STO-3G benchmark comprises 55 main SCF jobs with 42 analytic gradients, plus three pilot jobs and one gradient. Sixteen RBF fits select two surrogates using fixed 17/8/8/9 training/validation/test/OOD partitions. Joint fitting improves test energy RMSE to 5.56135×10⁻⁶ Hartree but worsens gradient RMSE to 1.63747×10⁻⁴ Hartree Å⁻¹, compared with energy-only 1.15617×10⁻⁴. Both extrapolate worse than the linear baseline. “Energy-only” fitting still uses validation gradients for hyperparameter selection. [Electronic records](../../quantumequi/results/electronic/summary.json).

Exact source-path replay uses 420 old-band evaluations plus 14 fresh final-image evaluations. Sequential in-place updates mix old forces with partly updated neighbors; 30 loops unconditionally return “converged.” The forward candidate force is 0.0183398705 in nominal kcal mol⁻¹ Å⁻¹, and endpoints are not stationary. Reviewed simultaneous-update CI-NEB instead checks endpoints and both projected/climbing forces. Forty main runs on two analytic surfaces recover known forward barriers 1.372719661 and 4 in dimensionless units; six pilot runs are separate. These are algorithm controls, not molecular barriers. [Path summary and replay](../../quantumequi/path_notes.md).

Mass-weighted Hessians project translations/rotations by SVD, retaining five rigid directions for linear and six for nonlinear molecules. The nominal kcal mol⁻¹ Å⁻² amu⁻¹ conversion factor is 108.591358535 cm⁻¹, not the source's 1302.83. Source index-based deletion discards negative modes. The reviewed source has one significant negative projected mode but fails stationarity, so no corrected source Gibbs energy is released. Eight main finite-difference source Hessians and one autograd reference separate dtype/step errors. Four artificial controls support 80 ideal-gas temperature/pressure cases and 72 low-frequency sensitivities. Translation, rotation, vibration, electronic degeneracy, pressure and symmetry are explicit; the original entropy constant and hardcoded plotted shifts are not thermal calculations. [Thermochemistry and rejection](../../quantumequi/results/thermochemistry/summary.json).

## S11. H₂ extension: electronic bias, learned error and nuclear approximation

The electronic calculations use the Psi4 software lineage [7], with the cc-pVDZ basis-family attribution given by Dunning [8]. The nuclear model uses the established Morse potential [11]; it is an approximation parameterized from the recorded electronic quantities.

RHF/UHF/FCI scans cover 25 distances and two bases. The 183 main energy drivers comprise 150 scan values, two isolated atoms, 21 minimum-search and ten curvature evaluations. Fourteen pilot calls and six recovered FCI calls give 203 actual drivers; six earlier option failures occurred before dispatch. All 50 same-basis energy comparisons satisfy FCI≤UHF≤RHF within 10⁻⁸ Hartree. UHF spin contamination becomes substantial; FCI spin diagnostics were unavailable and are not inferred from reference orbitals. FCI is exact only within the chosen finite basis.

**Table S2. FCI electronic references; no nuclear or environmental corrections.**

| Quantity | STO-3G | cc-pVDZ |
|---|---:|---:|
| Equilibrium distance / Å | 0.734865227 | 0.760893445 |
| Well depth / Hartree | 0.204142352 | 0.165116174 |
| Local curvature / Hartree Å⁻² | 1.703746695 | 1.307967354 |
| RHF−FCI at 4 Å / Hartree | 0.318301388 | 0.216407978 |
| UHF ⟨S²⟩ at 4 Å | 0.999980059 | 0.999764680 |

The new 3537-parameter distance-message network reuses only old RHF labels. Three paired seeds, two objectives and 800 epochs give six main fits and 4800 optimizer steps; two 200-epoch pilots are separate and evaluate only training/validation. Training-only scales define energy and gradient losses, while both objectives select checkpoints by the same validation energy-plus-gradient score. All selected epochs are at the budget boundary. Joint ensemble test energy/gradient RMSE is 0.000142893/0.00119298 in Hartree/Hartree Å⁻¹, and OOD values are 0.00473607/0.0127819. Each paired joint run improves, but all neural test errors exceed frozen RBF errors. Three-seed spread is uncalibrated.

At 21 shared geometries, 126 rows verify

\[
\widehat E-E_{FCI}=(\widehat E-E_{old\ RHF})+(E_{old\ RHF}-E_{new\ RHF})+(E_{new\ RHF}-E_{FCI}).
\]

Replay drift and identity residual are zero at saved precision. At 2.70 Å, joint seed 7301 has learning error 0.000423247 Hartree but total FCI error 0.254066902 Hartree, dominated by RHF bias 0.253643655 Hartree. Signed components add; RMSEs do not. FCI labels never entered training. [Correlation](../../quantumequi/results/extensions/correlation/summary.json), [learning](../../quantumequi/results/extensions/learning/summary.json), [decomposition](../../quantumequi/results/extensions/error_budget/summary.json).

Nuclear calculations parameterize V(r)=De[1−exp(−a(r−re))]² with a=√(k/2De), using bare proton/deuteron masses and J=0 second-order finite differences. They do not propagate nuclei on the actual FCI curve: the Morse representation's 25-point energy RMSE is 0.00757282/0.00420654 Hartree for STO-3G/cc-pVDZ. Main/pilot/tail studies contain 38/6/7 eigenproblems. On [0,8] Å with 1200 intervals, cc-pVDZ-derived H₂/D₂ gaps are 4117.259/2966.494 cm⁻¹, without intensity or selection-rule predictions. Coarse grids can create spurious near-threshold states, while the main H₂ domain misses the seventeenth true Morse state. Its analytic binding is 6.69637×10⁻⁷ Hartree and decay length 15.0913 Å. Extending to 96 Å and 115200 intervals recovers it but leaves 1.564% binding error. Bound-only thermal sums omit continuum, rotation and translation; they are not molecular Gibbs energies. [Main nuclear results](../../quantumequi/results/extensions/vibration/summary.json), [tail diagnostics](../../quantumequi/results/extensions/vibration/tail_followup/diagnostics.csv).

## S12. Software, source index and reproducibility

**Table S3. Recorded software, rather than a dependency installation prescription.**

| Component | Recorded version |
|---|---|
| Python | 3.12.14 |
| NumPy, principal numerical / Psi4 environment | 2.4.6 / 2.5.2 |
| SciPy / pandas | 1.18.0 / 2.3.3 |
| scikit-learn / Matplotlib | 1.9.0 / 3.11.1 |
| RDKit / PyTorch | 2026.03.5 / 2.10.0 |
| xTB / Psi4 | 6.7.1 / 1.11 |

CPU studies generally use one numerical thread. Exact seeds, architecture settings, tolerances, charges/spins and historical execution hashes are retained per module. Source code changed after some pilots, metadata finalization and tail additions; archived snapshots remain authoritative for the executed version. A small residual checks a discretized equation, a holdout score checks a specified partition, and a hash checks file identity. None substitutes for calibration or experimental replication.

**Table S4. Source and implementation index.** Links are relative to this supporting-information directory. Numerical data are linked in the corresponding sections above.

| Case | Preserved source | Reviewed implementation or method index |
|---|---|---|
| Initial | [original Python](../../data/original/simulate_all_topics.py) | [English methods](../../reports/technical_report_english.md) |
| Production | [original Python](../../production/source/run_production_pipeline_original.py) | [audit code](../../production/scripts/audit_pipeline.py) |
| Research | [original Python](../../research/source/run_research_engine_original.py) | [audit code](../../research/scripts/audit_research.py) |
| Closed loop | [original Python](../../closed_loop/source/run_closed_loop_platform.py) | [audit code](../../closed_loop/scripts/audit_closed_loop.py) |
| Toolkit | [original Python](../../toolkit/source/run_deployment_and_figures.py) | [audit code](../../toolkit/scripts/audit_toolkit.py) |
| ElectraTwin and extensions | [original Python](../../electratwin/source/electratwin_core.py) | [module/code index](../../electratwin/README.md) |
| ElectroGraph | [original Python](../../electrograph/source/electrograph_kmc_core.py) | [module/code index](../../electrograph/README.md) |
| SynthaPore | [original Python](../../synthapore/source/synthapore_engine.py) | [module/code index](../../synthapore/README.md) |
| QuantumEqui original | [original Python](../../quantumequi/source/quantum_egnn_neb_engine.py) | [module/code index](../../quantumequi/README.md) |
| QuantumEqui extensions | [executed quantum snapshot](../../quantumequi/results/extensions/correlation/executed_code.py.txt) | [extension code and data index](../../quantumequi/EXTENSIONS.md) |

No study is rerun to build this SI. From the repository root, saved-artifact validation can be invoked with an already configured interpreter:

```powershell
$chem = $env:CHEM_PYTHON
$env:OMP_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
& $chem scripts/validate_release.py
& $chem quantumequi/scripts/validate_extensions.py
```

Consult each module's index for its own validator and scientific entry points. Deliberate recalculation requires a separate checkout, appropriate existing native runtimes, and `PSI4_PYTHON`/`XTB_EXE` where applicable. Some early scripts overwrite outputs; later studies reject completed destinations. These are distinct from read-only result inspection. Convenience launchers were sometimes checked only with dry runs; the records state which individual stages actually executed. Regenerated results need new numerical review and figure QA before replacing hashes. Neither universal bitwise reproducibility nor journal acceptance is inferred from the released tests. The integrated evidence supports bounded computational comparisons and explicit failure diagnosis, while calibrated organic electrochemistry, real catalyst structures, solvent/electrode thermodynamics and experimental validation remain outside the executed scope.

## S13. Post hoc manuscript arithmetic

The manuscript adds arithmetic over frozen data, with no new quantum job, training or eigensolve: 12 MSE attribution groups, six paired comparisons and three Richardson estimates. At the eight matched test geometries, joint fitting reduces RHF learning RMSE by 60.240–73.649%, but total FCI RMSE by only 0.00838–0.02168%. At the five matched OOD geometries, the respective ranges are 62.528–97.445% and 3.403–8.214%; this subset differs from the original nine-point OOD set. The identity MSE(total)=MSE(learning)+MSE(bias)+2 mean(learning×bias) closes within 1.04083×10⁻¹⁷ Hartree². The cross term need not be positive, so these terms are not variance fractions.

At fixed 96 Å boundary, second-order Richardson extrapolation gives a weakest-binding relative deviation of −0.026882%. This is a post hoc asymptotic estimate, not another eigenproblem or certified convergence; the finest single-grid error remains 1.564%. [Analysis summary](../results/analysis_summary.json), [paired gains](../results/paired_gain_transfer.csv), [MSE attribution](../results/mse_attribution.csv), [tail extrapolation](../results/tail_richardson.csv), [residual/observable comparison](../results/residual_observable_comparison.csv).

## References

[7] D. G. A. Smith; L. A. Burns; A. C. Simmonett; et al. Psi4 1.4: Open-source software for high-throughput quantum chemistry. *The Journal of Chemical Physics* **2020, 152, 184108**. [doi:10.1063/5.0006002](https://doi.org/10.1063/5.0006002)

[8] Thom H. Dunning, Jr. Gaussian basis sets for use in correlated molecular calculations. I. The atoms boron through neon and hydrogen. *The Journal of Chemical Physics* **1989, 90, 1007–1023**. [doi:10.1063/1.456153](https://doi.org/10.1063/1.456153)

[11] Philip M. Morse. Diatomic Molecules According to the Wave Mechanics. II. Vibrational Levels. *Physical Review* **1929, 34, 57–64**. [doi:10.1103/PhysRev.34.57](https://doi.org/10.1103/PhysRev.34.57)

[12] Graeme Henkelman; Blas P. Uberuaga; Hannes Jónsson. A climbing image nudged elastic band method for finding saddle points and minimum energy paths. *The Journal of Chemical Physics* **2000, 113(22), 9901–9904**. [doi:10.1063/1.1329672](https://doi.org/10.1063/1.1329672)

[13] David L. Mobley; J. Peter Guthrie. FreeSolv: a database of experimental and calculated hydration free energies, with input files. *Journal of Computer-Aided Molecular Design* **2014, 28, 711–720**. [doi:10.1007/s10822-014-9747-x](https://doi.org/10.1007/s10822-014-9747-x)

[17] Daniel T. Gillespie. Exact stochastic simulation of coupled chemical reactions. *The Journal of Physical Chemistry* **1977, 81(25), 2340–2361**. [doi:10.1021/j100540a008](https://doi.org/10.1021/j100540a008)

[18] H. J. C. Berendsen; J. P. M. Postma; W. F. van Gunsteren; A. DiNola; J. R. Haak. Molecular dynamics with coupling to an external bath. *The Journal of Chemical Physics* **1984, 81(8), 3684–3690**. [doi:10.1063/1.448118](https://doi.org/10.1063/1.448118)

[19] Graeme Henkelman; Hannes Jónsson. Improved tangent estimate in the nudged elastic band method for finding minimum energy paths and saddle points. *The Journal of Chemical Physics* **2000, 113, 9978–9985**. [doi:10.1063/1.1323224](https://doi.org/10.1063/1.1323224)
