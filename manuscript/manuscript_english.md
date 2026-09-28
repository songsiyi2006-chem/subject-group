# Reference bias and observable dependent convergence in computational chemistry

Siyi Song

Guangxi Normal University, China

Research manuscript · 28 September 2026

## Abstract

Computational chemistry workflows increasingly combine electronic references, learned potentials, numerical solvers and downstream observables. Validation at one stage does not establish accuracy at the next. Here, a reproducible collection of molecular and model calculations is organized around two controlled chains: supervised learning of H₂ energies against alternative electronic references, and nuclear motion on FCI-parameterized Morse potentials. The electronic extension comprises 203 converged RHF, UHF and finite-basis FCI energy-driver calls; six main neural runs compare energy-only and energy–gradient objectives at paired initializations. Post hoc evaluation at eight identical test geometries shows that gradient supervision reduces energy RMSE against RHF labels by 60.24–73.65%, while reducing total error against same-basis FCI by only 0.00838–0.02168%. An exact signed decomposition retains the learning–reference cross term and identifies reference bias as the dominant limitation. In the nuclear chain, matrix residuals near machine precision coexist with an omitted shallow bound state and substantial binding-energy error. Expanding the radial domain and refining the grid recovers the state; a separately reported, fixed-domain Richardson estimate reduces its discrepancy from the analytic Morse reference to 0.02688%, without establishing ab initio spectroscopic accuracy. Supporting transport, kinetics, molecular-learning, path and analytical examples test the reach and limits of the same reasoning. The contribution is an auditable set of coupled counterexamples and an observable-specific validation procedure, rather than a new electronic-structure or learning algorithm. These results show quantitatively why improved fitting and successful internal checks can leave the requested physical conclusion largely unchanged.

Keywords: computational chemistry; machine learning potentials; reference bias; error attribution; nuclear vibration; reproducibility

## 1 Introduction

Machine learning potentials connect electronic-structure data with molecular simulation by replacing expensive evaluations with differentiable approximations. Modern models incorporate geometric symmetries and train on energies and forces, while electronic-structure calculations define the labels that are treated as targets [@unke2021; @batzner2022]. This separation creates several meanings of accuracy. A model may accurately interpolate a chosen quantum approximation, have internally consistent derivatives and remain inaccurate with respect to a better electronic reference. A numerical solver can accurately solve a discretized operator while the discretization or boundary conditions still distort a physical observable. These questions require different comparisons.

The distinction matters for organic electrosynthesis, where catalyst identity, solvent, electrode conditions, charge transfer and product selectivity must ultimately be connected. Experimentally established transformations from Haitao Tang and collaborators provide chemical motivation for computation-supported electrosynthesis [@tang1; @tang2]. They do not provide calibration data for the model studies considered here. In particular, an executed surrogate or transport solver does not establish an electrocatalytic mechanism, an oxidation potential or a preparative operating window. A numerical demonstration becomes chemically informative only when its reference, domain and target observable are specified.

Prior work has already established the importance of realistic evaluation, uncertainty assessment and physically appropriate machine learning representations [@uncertainty; @evaluation]. The present study does not claim that reference bias, distribution shift or discretization error is a newly discovered phenomenon. Its objective is to connect these familiar problems within inspectable calculation chains and quantify how much improvement at one level reaches the next. The analysis preserves unsuccessful trials and separates direct calculations from cached evaluations, synthetic labels and proposed inputs. This accounting is especially important when many small workflows are combined into a single software collection.

Two chains provide the principal evidence. In the first, paired neural fits use the same RHF/STO-3G training labels, and both are evaluated at geometries with fresh RHF and FCI references. This enables signed error attribution without interpolation between electronic datasets. In the second, FCI equilibrium quantities parameterize a Morse potential whose exact full-line spectrum provides a mathematical comparator for radial finite differences. Independent changes to grid spacing and radial extent reveal which observables remain unresolved. Earlier molecular-learning, reaction-path, thermochemical, transport, stochastic-kinetic and analytical calculations are included as supporting cases. They extend the range of failure modes but are not pooled into a universal performance score.

The study addresses three questions. How much of a reduction in fitting error survives a change to a correlated electronic reference? Can a small algebraic residual certify a near-threshold molecular-model observable? Which conclusions remain unsupported when internal consistency, conservation or successful execution is used as a substitute for the relevant external comparison? Figure 1 defines the calculation design and separates the two principal chains.

![Calculation design](figures/Fig1_Evidence_Design_english.png)

Figure 1. Two controlled calculation chains used for error attribution. Neural training uses RHF labels; FCI is introduced only for independent reference comparison. Nuclear dynamics uses an FCI-parameterized Morse approximation with its own model, boundary and discretization errors. The supporting cases retain their original data types and are not treated as additional H₂ observations.

## 2 Methods

### 2.1 Evidence set and study design

The source collection is frozen at repository commit 1ad05c243155a9f64b18e0d58c665424fde919a0. It contains the initial five-topic study and subsequent production, research, closed-loop, analytical, transport, molecular-graph, porous-material-model and quantum/path studies. Supporting Information S1–S12 documents the complete scope, including components that are not used to establish the main quantitative conclusions. Actual molecular calculations, public experimental labels, analytic controls, synthetic objectives, simulated feedback and unexecuted quantum inputs are distinguished at source-file level.

The main H₂ analysis is a small-system benchmark. Its geometries are correlated points on one dissociation curve, not independent samples of chemical space. Three random seeds probe initialization dependence, not population-level statistical uncertainty. The new error-transfer and Richardson analyses are explicitly post hoc operations on frozen data. No new electronic energy, neural training run or nuclear eigenproblem is counted for these arithmetic analyses. Model selection in the main H₂ study uses validation data; its pilot exposes only training and validation records. Some earlier supporting pilots include held-out outputs, so the complete historical collection cannot be described as an end-to-end blinded evaluation.

### 2.2 Electronic reference calculations

Psi4 1.11 [@psi4] was used for neutral H₂ in C1 symmetry at 25 internuclear separations between 0.50 and 4.00 Å. RHF, UHF and FCI were evaluated in STO-3G and cc-pVDZ [@dunning] with a common geometry at each point. The nuclei lie at ±R/2. SCF used a CORE guess, PK integrals, energy and density convergence thresholds of 10⁻¹² and a limit of 200 iterations. UHF used orbital mixing and followed internal instabilities; RHF stability was checked. FCI used the RHF reference, one root, no frozen orbitals, an energy threshold of 10⁻¹² and a residual threshold of 10⁻¹⁰. The two basis sets yield 4 and 100 determinants in the Nα=Nβ=1 sector, as verified from saved logs.

The main electronic study comprises 150 curve energies, two isolated doublet H calculations, 21 FCI minimum-search evaluations and ten FCI curvature evaluations. The search interval is 0.6–0.9 Å. Five-point curvature estimates use h=0.0025 Å; the electronic depth is Dₑ=2E(H)−E(minimum). Twenty successful pilot/recovery energy calls are counted separately, yielding 203 actual energy-driver calls. Six earlier optional-diagnostic configurations were rejected before dispatch and are not counted as energy calculations. The original study's 58 quantum jobs remain a separate record. Underlying SCF phases within FCI are not counted again.

FCI removes configuration truncation within the specified orbital basis and electronic Hamiltonian. It does not remove basis incompleteness, relativistic, nonadiabatic or environmental errors. STO-3G and cc-pVDZ are not nested contracted bases, so their observed ordering is not used as a general variational theorem. The stored FCI spin expectation is unavailable; the spin of its RHF reference orbitals is not substituted for a CI-wavefunction diagnostic.

### 2.3 Paired molecular energy learning

A 3537-parameter distance-message neural network is trained on the earlier 42-point RHF/STO-3G dataset, split into 17 training, 8 validation, 8 test and 9 stretched out-of-distribution geometries. A shared hydrogen embedding, two width-16 message/update layers with SiLU nonlinearities, two directed edges and a summed atomic-energy readout define an invariant scalar energy. Cartesian forces follow by automatic differentiation. This construction follows message-passing and continuous molecular-energy representations [@mpnn; @schnet], but is not presented as a reproduction of either cited architecture. H₂ contains only one independent distance; the model is consequently a radial regressor in this application.

Training-only energy and radial-gradient standard deviations, 0.0937046071 Hartree and 0.318790348 Hartree Å⁻¹, normalize the losses. The energy-only objective minimizes normalized energy mean-square error; the joint objective adds normalized radial-gradient mean-square error. Both select checkpoints using the same sum of normalized validation energy and gradient RMSEs. Thus energy-only training still uses validation gradients and the training gradient scale. Seeds 7301, 7302 and 7303 are paired across objectives with identical initialization hashes. Each main run uses 800 full-batch float64 Adam steps, cosine learning-rate decay from 0.003 to 0.00015, gradient clipping at norm 5 and zero weight decay. All selected checkpoints are at epoch 800; this is a budget endpoint, not proof of optimization convergence.

Frozen energy-only and joint Gaussian-RBF fits and a linear interpolation/extrapolation baseline are retained. They use the original label split. Energy invariance and force covariance are checked under rotations, reflections, translations and atom permutation, and every Cartesian force component is compared with central energy differences. These tests evaluate the learned scalar function's mathematical consistency. They do not compare its forces with correlated electronic forces.

### 2.4 Signed attribution and improvement transfer

The reference join is performed at identical saved geometries. It includes 21 separations and 126 model/seed records, comprising all eight test geometries and five of the nine original OOD geometries. Let ε be neural error relative to the old RHF label, d the old-to-fresh RHF replay difference, and b the fresh RHF–FCI difference. Then

$$
\begin{aligned}
\widehat E-E_{\rm FCI}&=\epsilon+d+b,\quad \epsilon=\widehat E-E_{\rm RHF}^{\rm old},\\
d&=E_{\rm RHF}^{\rm old}-E_{\rm RHF}^{\rm new},\quad b=E_{\rm RHF}^{\rm new}-E_{\rm FCI}.
\end{aligned}
$$

The replay difference is zero at the recorded precision for all joined geometries. The new grouped analysis therefore retains the exact two-term mean-square identity

$$
{\rm MSE}_{\rm FCI}=\langle\epsilon^2\rangle+\langle b^2\rangle+2\langle\epsilon b\rangle.
$$

The cross term can be negative. Its value is reported directly rather than interpreted as a positive variance contribution. For the same finite set, the reverse triangle inequality and triangle inequality give

$$
\left|{\rm RMSE}(b)-{\rm RMSE}(\epsilon)\right|
\leq {\rm RMSE}(\epsilon+b)
\leq {\rm RMSE}(b)+{\rm RMSE}(\epsilon).
$$

These are standard identities and bounds, not new theory. They become informative here because all terms are available at matched geometries. A relative improvement is 100[1−RMSE(joint)/RMSE(energy-only)], computed once against RHF and once against FCI for each seed and split. Test and matched-OOD results are reported separately. No FCI labels enter fitting, checkpoint selection or hyperparameter adjustment.

### 2.5 Nuclear model and numerical error separation

The electronic minimum, isolated-atom depth and local FCI curvature define a Morse approximation [@morse]

$$
V(r)=D_e[1-e^{-a(r-r_e)}]^2,\qquad a=\sqrt{k/(2D_e)}.
$$

This is a three-quantity parameterization, not a nonlinear fit to the full electronic curve. H₂ and D₂ share the same electronic potential; reduced masses use half the bare proton/deuteron masses, 1.0072764665789/2.013553212544 u. A second-order finite-difference radial Hamiltonian with Dirichlet endpoints represents J=0 nuclear motion. A symmetric tridiagonal eigensolver extracts states below Dₑ. An independently specified Morse control, Dₑ=0.2 Hartree, a=1.8 Å⁻¹ and rₑ=1 Å, checks the implementation against the exact full-line spectrum.

The principal grids use 400, 800 and 1200 intervals on [0,8] Å; domain variation and a deliberately displaced-wall control are recorded separately. Six final model/isotope spectra yield 130 main-grid levels. The main, pilot and shallow-state follow-up comprise 38, 6 and 7 unique eigenproblems, respectively. Full-line and finite radial boundary conditions are distinguished throughout. Bound-state partition sums are evaluated at 100, 298.15, 600, 1000, 2000 and 4000 K. Reported F is the bound-only vibrational Helmholtz contribution, not the full molecular Gibbs energy.

The follow-up varies radial extent and spacing independently. A new arithmetic analysis compares shallow-state binding errors with matrix residuals for the same saved eigenproblems. It also computes the standard second-order estimate Bext=(4B(h/2)−B(h))/3 at fixed right boundaries of 48 and 96 Å. This extrapolation is post hoc, assumes a leading quadratic grid error, and does not replace the archived single-grid spectrum or establish general convergence for other potentials.

## 3 Results

### 3.1 Electronic references separate correlation error from basis effects

All 50 same-basis, same-geometry comparisons satisfy EFCI≤EUHF≤ERHF within 10⁻⁸ Hartree. Twenty-eight UHF points are classified as lower spin-broken solutions and twenty-two as restricted-like. The first sampled lower branch occurs at 1.26 Å in both bases; this is not an exact location of the bifurcation. At 4 Å, RHF exceeds FCI by 0.318301388 Hartree in STO-3G and 0.216407978 Hartree in cc-pVDZ. UHF reduces this energy difference but has ⟨S²⟩ values of 0.999980059 and 0.999764680, respectively. The lower mean-field energy therefore accompanies substantial spin contamination.

Table 1. FCI quantities used to parameterize the nuclear models. Electronic depths use two isolated H atoms and exclude nuclear and thermal corrections.

| Quantity | STO-3G | cc-pVDZ |
|---|---:|---:|
| Equilibrium distance / Å | 0.734865227 | 0.760893445 |
| Minimum energy / Hartree | −1.137306051 | −1.163672981 |
| Electronic depth / Hartree | 0.204142352 | 0.165116174 |
| Curvature / Hartree Å⁻² | 1.703746695 | 1.307967354 |

The 4 Å FCI values remain 7.66273×10⁻⁶ and 4.93793×10⁻⁵ Hartree below their two-atom references. A visually flat stretched curve is therefore not treated as the asymptotic reference. Figure 2 shows the reference curves and the large discrepancy that persists when a restricted single-determinant description is retained during dissociation. These results reproduce established small-molecule physics and provide the controlled reference change needed for the learning analysis.

![Electronic reference curves](../quantumequi/reports/figures/extensions/Ext1_Dissociation_english.png)

Figure 2. Executed H₂ total energies for RHF, UHF and FCI in STO-3G and cc-pVDZ at 25 common separations. Total energies include nuclear repulsion. Lines connect saved points and do not represent additional calculations. FCI is a finite-basis electronic comparator.

### 3.2 Improved fitting transfers weakly to the correlated reference

Against the original RHF labels, joint supervision improves energy and gradient errors for all three paired seeds on both the full test and full OOD sets. Test energy RMSE decreases from 0.000226936–0.000757399 Hartree for energy-only fitting to 0.000090229–0.000216394 Hartree for joint fitting. Every neural model nevertheless has larger test energy and gradient errors than both frozen RBF models. The neural models extrapolate more accurately on the original nine-point stretched RHF set, a benefit confined to this reference and radial domain.

The matched FCI comparison changes the interpretation. RHF–FCI RMSE is 0.068653520 Hartree on the eight test geometries and 0.208255039 Hartree on the five matched OOD geometries. Joint networks have test FCI RMSEs of approximately 0.06864 Hartree even when their errors against RHF are below 0.00022 Hartree. For all three paired seeds, reducing the fitting loss yields only a small change in total test error relative to FCI (Table 2 and Figure 3).

Table 2. Relative reduction in energy RMSE after adding gradient supervision. Every comparison uses the same geometries and initialization seed; OOD contains five matched points here, not the nine points of the original learning evaluation.

| Split | Seed | Reduction against RHF / % | Reduction against FCI / % |
|---|---:|---:|---:|
| Test | 7301 | 73.64861 | 0.021679 |
| Test | 7302 | 71.42932 | 0.015494 |
| Test | 7303 | 60.24035 | 0.008384 |
| Matched OOD | 7301 | 97.44489 | 3.403104 |
| Matched OOD | 7302 | 70.12125 | 8.214032 |
| Matched OOD | 7303 | 62.52848 | 4.284502 |

![Transfer of fitting improvement](figures/Fig3_Gain_Transfer_english.png)

Figure 3. Percentage reduction in energy RMSE from energy-only to joint training, evaluated against two references at identical geometries. Connecting lines identify paired seeds. The eight-point test set and five-point matched OOD set are separate comparisons. These are observed changes in three runs, not confidence intervals or population estimates.

The MSE decomposition explains the weak transfer. For joint seed 7301 on the test set, the learning term is 1.64325×10⁻⁸ Hartree², the reference term is 4.71331×10⁻³ Hartree² and the cross term is −1.47269×10⁻⁶ Hartree². The negative cross term makes the total error slightly smaller than the RHF–FCI reference error alone; it is cancellation, not learned electron correlation. The same sign occurs for all three joint test fits. In matched OOD, the cross term is positive and increases the total error. All twelve grouped identities close within 1.05×10⁻¹⁷ Hartree². Figure 4 retains the signs, which would be lost in a chart of positive percentage contributions.

![Signed mean-square error attribution](figures/Fig4_MSE_Attribution_english.png)

Figure 4. Learning, reference, cross and total mean-square terms for the three joint models. The symmetric logarithmic axis retains negative cross terms; values are expressed in mHartree². The reference term is shared across seeds because the geometries and electronic calculations are identical. Terms are not statistically independent variance components. Panel ranges differ.

The mathematical model checks remain successful. Thirty-six transformation probes have maximum discrepancies of 2.22×10⁻¹⁶ in the corresponding energy and force units. Seventy-two full Cartesian finite-difference cases reach a maximum force discrepancy of 2.22813×10⁻¹⁰ Hartree Å⁻¹ at h=10⁻⁵ Å. These results establish consistent derivatives and symmetry of the learned function while leaving the reference bias intact. Three-seed ensemble spread is also insufficient as a calibrated uncertainty statement: only 75% of the eight test errors fall within twice the sample spread, for both energy and gradient. Eight correlated distances and three fits do not establish coverage for other molecules.

### 3.3 Nuclear accuracy depends on the requested observable

The FCI-derived Morse potentials have electronic-curve RMSEs of 0.00757282 Hartree for STO-3G and 0.00420654 Hartree for cc-pVDZ across the 25 sampled separations. This model error remains even if the nuclear equation for the Morse potential is solved exactly. The finest main grid gives 0→1 gaps of 4722.097/3397.084 cm⁻¹ for H₂/D₂ with STO-3G parameters and 4117.259/2966.494 cm⁻¹ with cc-pVDZ parameters. The corresponding anharmonic isotope ratios are 1.39004412 and 1.38792116, compared with the harmonic ratio 1.41386262. Figure 5 shows the compression of successive Morse levels relative to a harmonic approximation. No transition intensity or experimental spectral assignment is calculated.

![Isotope levels](../quantumequi/reports/figures/extensions/Ext5_Anharmonic_Levels_english.png)

Figure 5. H₂ and D₂ levels on the cc-pVDZ-derived Morse model. Solid levels are the first five bound energies from the main 1200-interval grid; dashed levels are harmonic comparisons. Energies are measured from the model minimum. This representation does not substitute for nuclear propagation on the full FCI curve.

Low levels show approximately second-order convergence, with observed ZPE orders of 2.00057–2.00093. Main-grid matrix residuals are no larger than 1.43×10⁻¹⁴ Hartree, but the maximum error among the first three levels ranges from 3.975 to 5.705 cm⁻¹ across the six spectra, relative to the analytic model at the finest main grid. Near the dissociation threshold, the cc-pVDZ-derived H₂ model has seventeen exact full-line states, whereas [0,8] Å returns sixteen. Its shallowest exact binding is only 6.69637×10⁻⁷ Hartree, with an amplitude decay length of 15.0913 Å. A domain adequate for low levels can therefore omit a state with a long tail.

The follow-up also demonstrates compensation between boundary and grid errors. At 48 Å with 7200 intervals, all seventeen states are found, yet the shallow-state binding has 125.201% relative error. The absolute binding error is about 1.72×10⁸ times the matrix residual for that same eigenproblem. At 96 Å with 115200 intervals, the single-grid binding error decreases to 1.564%, while the matrix residual grows to 1.26610×10⁻¹² Hartree. Ranking these calculations by algebraic residual alone would reverse their physical accuracy ordering (Figure 6a). The residual and binding error are not interchangeable error bounds.

Table 3. Fixed-boundary, second-order extrapolation of the shallow-state binding. Signed errors are relative to the exact full-line Morse value. These are new arithmetic estimates from archived eigenvalues, not new spectra.

| Boundary / Å | Fine intervals | Fine-grid error / % | Richardson error / % |
|---:|---:|---:|---:|
| 48 | 14400 | 25.76358 | −7.38220 |
| 48 | 28800 | 5.20144 | −1.65261 |
| 96 | 115200 | 1.56403 | −0.02688 |

![Observable dependent numerical convergence](figures/Fig6_Observable_Convergence_english.png)

Figure 6. (a) Matrix residual and absolute shallow-state binding error for the same saved eigenproblems after seventeen states are present. (b) Fine-grid and fixed-boundary Richardson errors. The 96 Å estimate improves agreement with the analytic Morse reference; it does not establish the accuracy of the electronic parameterization or a general extrapolation guarantee.

The finest 96 Å extrapolated binding is 6.69456672×10⁻⁷ Hartree, a signed discrepancy of −0.026882% from the analytic value. The weaker 48 Å extrapolations retain larger errors, consistent with unresolved boundary effects and higher-order discretization terms. Only two 96 Å grids are available; the leading-order extrapolation is not independently verified by a third fine grid. At 4000 K, the bound-state vibrational F error decreases from approximately −4.14×10⁻⁶ Hartree on the coarse long-domain grids to −6.46×10⁻⁸ Hartree on the finest grid. A relatively insensitive partition sum can therefore look satisfactory while the shallowest individual state remains inaccurate. Continuum, rotation, translation, nuclear-spin statistics and solvent contributions are absent from this comparison.

### 3.4 Supporting cases delimit the generality of internal checks

The broader collection contains several independent reasons to require target-specific comparisons. In an untrained equivariant source potential, covariance and derivatives pass numerical checks, but 1088 parameters in a final coordinate head are disconnected from its energy. The source path routine reports convergence after a fixed iteration budget; its candidate has maximum atomic force 0.0183399 in nominal units and is not stationary. A separate CI-NEB implementation follows the established climbing-image formulation [@neb] and converges on 36 analytic sweep cases plus four reversal checks. Its best barrier discrepancies are of order 10⁻¹¹ in dimensionless model units. Those results verify the controlled surfaces and do not establish a chemical transition state.

The associated thermochemical audit exposes a unit-conversion error of approximately twelvefold and a projected candidate gradient about 1834 times its stated threshold. A Gibbs correction is withheld for that nonstationary, uncalibrated potential. Explicit ideal-gas calculations on independent stationary controls instead retain translation, rotation, vibrational excitation, pressure and symmetry conventions. Seventy-two low-frequency sensitivity cases show that replacing a 0.1 cm⁻¹ mode by 100 cm⁻¹ changes G by 4.09848 kcal mol⁻¹. This intervention is not a hindered-rotor calculation and cannot be treated as a universal correction.

For experimental-label learning, three neural seeds on a 256-molecule FreeSolv subset give test RMSEs of 1.8972, 1.7468 and 1.2496 kcal mol⁻¹, compared with 1.6534 for ridge regression. Two neural seeds therefore lose to the simpler baseline. These data concern hydration free energies [@freesolv]; they do not validate oxidation potentials or reaction yields. Synthetic equivariant-denoising studies likewise include a seed that underperforms smoothing under distribution shift. Since some historical pilot outputs include held-out cases, these comparisons are descriptive retrospective evidence rather than a completely blinded model-development benchmark.

Transport and analytical examples test different failure modes. A conservative four-species electrochemical model can satisfy charge balance to numerical precision and still predict low net desired-product Faradaic efficiency because the assumed network permits overoxidation. Conservation does not imply selectivity. In sequential optimization, thousands of campaign evaluations reuse a finite pool of precomputed model values and must not be reported as new experiments. Simulated chromatographic peak recovery deteriorates when the peak width is misspecified, while source HPLC and NMR yield calculations disagree sharply. These are checks on assumptions and analytical arithmetic, not validated experimental observations. Complete parameters, source paths, counts and negative results for every earlier module are retained in the Supporting Information.

## 4 Discussion

### 4.1 What the combined analysis contributes

The main contribution is the joint quantification of two distinct barriers to interpreting improved computation. First, the matched-geometry study measures how a substantial supervised-learning gain is suppressed by an unchanged electronic reference bias. Second, the nuclear study measures how excellent discrete residuals can coexist with a missing or inaccurate near-threshold observable, and separates fixed-boundary refinement from domain enlargement. The connection is an operational validation procedure: identify the final observable, distinguish the reference and numerical approximations that affect it, and test each with an appropriate independent comparator before passing a result downstream.

This procedure is illustrated by exact signed arithmetic and controlled interventions, not by a new error identity or a new eigensolver. It is also narrower than a claim of a universally reliable AI chemistry platform. The same-basis FCI comparison provides a useful electronic target for H₂, while the full-line Morse spectrum provides a useful mathematical target for nuclear discretization. Neither is an experimental truth for general catalysis. The supporting cases show additional ways that execution, symmetry, conservation and optimization can become disconnected from their intended scientific claims; they do not estimate how frequently such failures occur in the literature.

Reference-limited improvement is particularly relevant when selecting the next calculation. If the measured reference-bias norm greatly exceeds the fitting-error norm on the target domain, further fitting alone has limited scope to change that domain's total error. The triangle bound makes this statement quantitative for the observed finite set. A better electronic reference, a changed representation of the target property, or direct calibration may then deserve priority. This is a conditional decision rule, not proof that model capacity or training data are generally unimportant. Error cancellation can also make a worse fit accidentally closer to a chosen higher-level reference; the signed cross term prevents that cancellation from being mistaken for a learned correction.

The electrosynthesis examples clarify what a chemically grounded extension would require. Tang and collaborators report a zinc single-atom cathode with anodic iodide mediation for pyridine hydroxymethylation [@tang1], an iron single-atom redox mediator for anodic organic synthesis [@tang2], and manganese-catalyzed oxygen evolution coupled to silane oxidation [@tang3]. Their tellurium-containing oxazolidinone synthesis supplies another experimentally defined transformation [@tang4]. These are distinct chemical systems and cannot be represented by relabeling the abstract A/P/B/D network or the H₂ benchmark. A future computational study would need complete structures, matching electrochemical conditions, assigned products and pathway-specific reference calculations for one chosen transformation. The papers motivate that requirement; no validation against their experiments is claimed here.

### 4.2 Relation to existing research and scope of novelty

Energy/force learning, equivariant representations, uncertainty estimation, reference-method dependence, and the distinction between numerical verification and physical validation are established research topics [@unke2021; @batzner2022; @uncertainty; @evaluation]. The potentially reusable result here is their connection through frozen, auditable examples with all intermediate quantities exposed. The H₂ learning comparison preserves initialization pairing, label-access rules and exact geometry matching; the nuclear comparison links the same eigenproblem's residual to a specific shallow-state error. The released records allow these conclusions to be checked without rerunning every scientific calculation.

An especially relevant comparison is the study by Kamath and colleagues, which already compared neural and Gaussian-process potential surfaces through both fitting quality and calculated vibrational spectra [@kamath2018]. The present neural-reference and Morse-nuclear chains are separate: the learned neural potential is not inserted into the nuclear solver. An end-to-end neural-potential-to-spectrum propagation study has therefore not been performed. The narrower distinction is the quantified attenuation of RHF fitting gains against matched FCI values, combined with a separate near-threshold grid/domain counterexample. This distinction prevents the earlier literature's contribution from being presented as new work here.

The present collection does not support claims of a novel catalyst, a new chemical reaction, chemical-space transferability, a new scalable quantum method or state-of-the-art learning accuracy. Its strongest contribution is a reproducible methodological case study. Establishing broader novelty would require multiple chemically distinct and independently assembled reference families, additional model architectures and matched compute budgets, and blinded external replication. The accompanying editorial assessment evaluates those gaps separately from the manuscript's demonstrated results.

### 4.3 Limitations and next validation steps

The primary electronic system has two electrons, the supervised input is effectively one-dimensional, and both basis sets are small. The stretched points intentionally expose a well-known mean-field limitation. Only three seeds and small correlated test sets are available. The new matched-OOD analysis includes five of nine originally designated OOD points; it does not silently extrapolate missing FCI labels. The nuclear potential is parameterized from three FCI quantities and differs from the sampled electronic curve. Its finite radial boundary problem is not identical to the analytic full-line comparator. The finest extrapolation is promising within that model but remains based on two grids.

Supporting chemistry workflows combine different evidence types and cannot establish a single end-to-end accuracy claim. In particular, no Cu reaction path has a validated quantum stationary point and connectivity analysis, no electrode/solvent model is calibrated against the cited electrosynthesis studies, and no new laboratory experiment is included. Literature examples establish relevance, not validation. A credible extension would predefine a reaction- and observable-specific reference ladder, retain a genuinely unexamined external test set, compare simple baselines and modern models under the same label/compute access, and trace uncertainties to the desired barrier, rate, selectivity or spectrum. Matched reference calculations and quantitative numerical convergence should precede mechanistic interpretation.

## 5 Conclusions

Paired H₂ fits show that large gains against RHF labels can yield very small gains against same-basis FCI. Retaining the signed learning–reference cross term distinguishes reference limitation from cancellation. A separate FCI-parameterized nuclear model shows that small matrix residuals do not establish completeness or accuracy of near-threshold bound states; grid and boundary errors must be varied separately and assessed for the requested observable. Fixed-domain extrapolation improves the shallow-state estimate in the analytic model, while leaving electronic-model limitations explicit. Together with the supporting workflow audits, these results define an inspectable route from algorithmic consistency to appropriately bounded scientific claims. The evidence supports this validation case study and its quantitative counterexamples, while broader chemical prediction remains a separate research task.

## Data and code availability

The complete underlying calculations are available in the [subject-group repository](https://github.com/songsiyi2006-chem/subject-group), with the frozen source study identified by commit 1ad05c243155a9f64b18e0d58c665424fde919a0. The manuscript directory contains the post hoc analysis script, row-level results, evidence catalog, verified references, figure provenance and manuscript sources. The Supporting Information maps every earlier module to its retained raw data, code and reports. Archived failures, pilot records and executed source snapshots are preserved. The repository integrity checks verify records and numerical identities; they do not certify experimental chemistry. Dependency and dataset licensing remain those of their respective upstream projects.

## Author and preparation statements

Siyi Song is the sole listed author, affiliated with Guangxi Normal University. AI-assisted tools were used for code implementation, orchestration, literature retrieval, translation and manuscript drafting; AI systems are not authors. Executed calculations and retained records are distinguished from proposed work throughout. This manuscript is a research draft requiring the author's final scientific approval and journal-specific disclosure review. Funding, competing-interest information and a correspondence address were not supplied and are not inferred.

## References

<!-- VERIFIED_REFERENCES -->
