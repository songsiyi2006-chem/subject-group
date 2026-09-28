# QuantumEqui extensions: correlated H₂ reference curves, supervised molecular energies and anharmonic nuclear models

## 1. Scope, evidence and calculation accounting

This extension adds executed electronic-structure calculations, supervised learning and nuclear eigenproblems to the earlier QuantumEqui audit. It asks three bounded questions: how electronic approximation changes H₂ dissociation, whether coordinate-gradient supervision improves a small neural energy model, and how isotope masses, discretization and finite boundaries affect an anharmonic nuclear model. These calculations concern H₂/D₂. They do not train or validate the original Cu-fragment network, establish an organic reaction mechanism, or provide experimental spectroscopy.

The evidence chain is deliberately explicit. Actual RHF/UHF/FCI calculations supply electronic reference curves. Neural fitting uses only the earlier frozen RHF/STO-3G data, while new FCI values enter a subsequent error decomposition. Nuclear calculations use an approximate Morse potential parameterized by FCI equilibrium quantities, rather than the complete FCI curve. Mathematical controls test the nuclear solver independently. The [extension index](../EXTENSIONS.md), [saved-evidence validation record](../results/extensions/validation.json) and [figure QA record](../results/extensions/figure_qa.json) provide the current artifact and review status.

| Workload unit | Main study | Separate pilot | Recovery or follow-up |
|---|---:|---:|---:|
| Actual electronic energy-driver calls | 183 | 14 | 6 |
| Configurations rejected before energy dispatch | 0 | 6 | 0 |
| Neural training runs | 6 | 2 | 0 |
| Neural optimizer steps | 4800 | 400 | 0 |
| Batched molecular energy evaluations | 121170 | 10050 | 0 |
| Tridiagonal nuclear eigenproblems | 38 | 6 | 7 |
| Summed nuclear interior grid nodes | 27562 | 3994 | 228593 |

Thus 209 attempted electronic configurations produced **203 actual, converged energy-driver calls**; six pilot configurations failed option validation before dispatch. The 183 main calls comprise 150 curve points, two isolated H atoms, 21 minimum-search evaluations and ten curvature evaluations. Across all stages, 87 calls explicitly request FCI; their underlying SCF phases are included in those calls, not additional quantum jobs. Neural and nuclear work introduces zero additional electronic-structure calculations. Forward evaluations, optimization steps and eigenproblems are different workload units.

All runs used one CPU thread. Recorded in-script wall times were 62.124 s for main quantum work, 5.219 s for its pilot and 2.695 s for FCI recovery; 41.792/3.973 s for main/pilot learning; and 0.419/0.088/1.976 s for main/pilot/tail nuclear calculations. These local timings exclude interpreter startup and do not predict scaling to larger molecules.

## 2. Correlated quantum curves and electronic approximation

Psi4 1.11 evaluated neutral H₂ at 25 separations from 0.50 to 4.00 Å with RHF, UHF and full configuration interaction in STO-3G and cc-pVDZ. Nuclei lie at ±R/2 in C1 symmetry. Calculations used a CORE guess, PK SCF, energy and density thresholds of 10⁻¹², and a 200-iteration SCF limit. UHF mixed its initial orbitals and followed internal instabilities; RHF stability was checked. The isolated neutral H reference was a doublet UHF calculation in each basis. These settings follow the documented [Psi4 SCF and stability workflow](https://psi4.github.io/psi4docs/master/scf.html) and [GUESS_MIX option](https://psicode.org/psi4manual/master/autodir_options_c/scf__guess_mix.html).

FCI used the RHF orbital reference, no frozen orbitals, one root, a 10⁻¹² energy threshold and a 10⁻¹⁰ residual threshold. The two bases contain two and ten spatial orbitals, giving four and 100 determinants in the Nα=Nβ=1 sector. FCI is exact only within this finite orbital space and the specified electronic Hamiltonian. It is not a complete-basis or experimental reference. [Psi4 DETCI documentation](https://psicode.org/psi4manual/master/detci.html) describes the method and implementation scope.

The pilot exposed two interface limitations. An optional S_SQUARED option was unsupported and caused six pre-driver rejections; those six FCI cases were recovered after removing that diagnostic. A determinant-count accessor returned invalid pointer-like integers in recovery metadata; counts were instead extracted from retained logs, with old values and hashes preserved. Neither issue was hidden as a successful original run. The saved FCI spin expectation is unavailable: the stored FCI-reference orbitals are RHF orbitals and cannot supply the CI wavefunction's spin expectation.

| Quantity | STO-3G | cc-pVDZ |
|---|---:|---:|
| RHF energy at 4 Å / Hartree | −0.614869974 | −0.782198208 |
| UHF energy at 4 Å / Hartree | −0.933166094 | −0.998569701 |
| FCI energy at 4 Å / Hartree | −0.933171362 | −0.998606186 |
| UHF ⟨S²⟩ at 4 Å | 0.999980059 | 0.999764680 |
| RHF − FCI at 4 Å / Hartree | 0.318301388 | 0.216407978 |
| FCI equilibrium distance / Å | 0.734865227 | 0.760893445 |
| FCI minimum energy / Hartree | −1.137306051 | −1.163672981 |
| Two isolated H energies / Hartree | −0.933163699 | −0.998556807 |
| Electronic well depth Dₑ / Hartree | 0.204142352 | 0.165116174 |
| Local FCI curvature / Hartree Å⁻² | 1.703746695 | 1.307967354 |

All 50 same-basis, same-geometry comparisons satisfy EFCI ≤ EUHF ≤ ERHF within 10⁻⁸ Hartree. UHF finds 28 lower spin-broken points and 22 restricted-like points, using ⟨S²⟩>10⁻⁴ and EUHF<ERHF−10⁻⁷ Hartree as the classification. The first sampled lower branch is at 1.26 Å in both bases; this does not locate an exact bifurcation. At long distance, the lower UHF energy accompanies substantial spin contamination. Positive sampled internal-stability eigenvalues establish local stability under that check, not a global Hartree–Fock minimum or spin purity. STO-3G and cc-pVDZ are not nested bases, so their observed cross-basis ordering is not a general variational theorem.

Bounded FCI minimum searches used [0.6,0.9] Å, followed by five energy evaluations at offsets −2h,−h,0,h,2h with h=0.0025 Å. The curvature is

\[
k=\frac{-E_{+2}+16E_{+1}-30E_0+16E_{-1}-E_{-2}}{12h^2},\qquad D_e=2E_H-E_{\min}.
\]

The two central-curvature step comparisons differ by 0.000102243 and 0.000078581 Hartree Å⁻², a numerical sensitivity diagnostic. FCI at 4 Å remains 7.66273×10⁻⁶ and 4.93793×10⁻⁵ Hartree below the respective two-atom references; 4 Å is therefore not treated as infinity. Well depths omit counterpoise, zero-point, thermal, relativistic and environmental corrections.

**Figure Ext1.** Executed total-energy dissociation curves for RHF, UHF and FCI: (a) STO-3G and (b) cc-pVDZ, each at the same 25 separations. [PNG](figures/extensions/Ext1_Dissociation_english.png) · [SVG](figures/extensions/Ext1_Dissociation_english.svg).

**Figure Ext2.** (a) Same-basis RHF and UHF energy differences from FCI, in mHartree; (b) UHF ⟨S²⟩ for both bases, with zero as the singlet reference. The plotted UHF diagnostic is not a measured FCI spin expectation. [PNG](figures/extensions/Ext2_Correlation_Spin_english.png) · [SVG](figures/extensions/Ext2_Correlation_Spin_english.svg).

## 3. Supervised distance-message energy learning

A genuine 3537-parameter PyTorch network was fitted to the frozen 42-point RHF/STO-3G reference: 17 training, eight validation, eight test and nine stretched out-of-distribution (OOD) geometries. It has a shared hydrogen embedding, two width-16 message/update layers with SiLU activation, two directed edges and a summed atomic-energy readout. Distance-only messages and scalar pooling make the energy invariant; coordinate autograd gives conservative Cartesian forces. This compact implementation follows general [message-passing](https://proceedings.mlr.press/v70/gilmer17a.html) and [continuous molecular-energy learning](https://proceedings.neurips.cc/paper/2017/hash/303ed4c69846ab36c2904d3ba8573050-Abstract.html) ideas, without claiming to reproduce SchNet. For H₂, one independent distance makes this radial regression, not evidence of chemical-space graph generalization.

\[
m_{ij}^{(l)}=\phi_l(h_i^{(l)},h_j^{(l)},\widetilde r_{ij}),\quad
h_i^{(l+1)}=h_i^{(l)}+\psi_l(h_i^{(l)},\sum_{j\ne i}m_{ij}^{(l)}),\quad
\mathbf F=-\nabla_{\mathbf X}\widehat E.
\]

All standardization uses training data only. The energy and gradient scales are sE=0.0937046071 Hartree and sg=0.318790348 Hartree Å⁻¹. The objectives are

\[
L_E=\langle[(\widehat E-E)/s_E]^2\rangle,\qquad
L_{E,g}=L_E+\langle[(\partial_R\widehat E-g)/s_g]^2\rangle.
\]

Both select checkpoints by validation RMSE(E)/sE + RMSE(g)/sg. Consequently, “energy-only” describes its fitting loss: selection still uses validation gradients and the training gradient scale. The earlier RBF comparison uses the same label-access distinction. No test/OOD labels determine parameters or hyperparameters, and no FCI labels enter training.

Two pilot runs used seed 7301 and 200 epochs, examining only training/validation. Main runs paired seeds 7301, 7302 and 7303 across objectives with identical initial-state hashes. Each used 800 full-batch float64 Adam steps, zero weight decay, cosine learning-rate decay from 0.003 to 0.00015 and a gradient-norm cap of five. All six best-validation checkpoints are at epoch 800; reaching the budget boundary is not optimization convergence. Complete learning curves and initial, periodic, selected and final states are retained.

| Model / seed | Test energy RMSE / Hartree | Test gradient RMSE / Hartree Å⁻¹ | OOD energy RMSE / Hartree | OOD gradient RMSE / Hartree Å⁻¹ |
|---|---:|---:|---:|---:|
| Neural energy-only / 7301 | 0.000486462 | 0.0121966 | 0.00767735 | 0.0147389 |
| Neural joint / 7301 | 0.000128189 | 0.00148374 | 0.000173536 | 0.00112220 |
| Neural energy-only / 7302 | 0.000757399 | 0.0121112 | 0.0277884 | 0.0571914 |
| Neural joint / 7302 | 0.000216394 | 0.00159951 | 0.00812651 | 0.0212530 |
| Neural energy-only / 7303 | 0.000226936 | 0.00509392 | 0.0159985 | 0.0406077 |
| Neural joint / 7303 | 0.0000902290 | 0.000771657 | 0.00596255 | 0.0161202 |
| Frozen RBF energy-only | 0.00000862425 | 0.000115617 | 0.175473 | 0.459351 |
| Frozen RBF joint | 0.00000556135 | 0.000163747 | 0.118352 | 0.448587 |
| Frozen linear baseline | 0.00110348 | 0.00223930 | 0.0444501 | 0.108744 |

Gradient supervision improves all paired test/OOD metrics. Nevertheless, every neural run has worse test energy and gradient errors than both frozen RBFs. Neural extrapolation is better on this particular stretched RHF set, but joint OOD energy RMSE varies from 0.000173536 to 0.00812651 Hartree across seeds. Neither selecting the best seed nor comparing one small radial dataset establishes a broadly superior architecture.

The three-seed joint ensemble gives test energy/gradient RMSE of 0.000142893/0.00119298 and OOD values of 0.00473607/0.0127819, in Hartree and Hartree Å⁻¹. Only 75% of the eight test errors lie within twice the sample spread for either quantity. The spread is uncalibrated: three correlated fits do not provide confidence intervals or guaranteed coverage.

At 0.74, 1.14 and 2.30 Å, 36 transformation checks give maximum energy-invariance and force-covariance discrepancies of 2.22×10⁻¹⁶ in their respective units. In 72 complete Cartesian finite-difference cases, the largest force error falls from 2.24512×10⁻⁴ to 2.22813×10⁻¹⁰ Hartree Å⁻¹ as the displacement decreases from 0.01 to 0.00001 Å. These checks validate derivatives of the learned scalar function, not its quantum accuracy.

**Figure Ext3.** (a) Common normalized validation score versus epoch for all paired seeds; (b) radial energy gradient, three-seed mean ± one sample standard deviation, against the RHF reference. Shading at 1.9–2.7 Å marks stretched OOD geometries; the uncertainty band is descriptive. [PNG](figures/extensions/Ext3_Learning_english.png) · [SVG](figures/extensions/Ext3_Learning_english.svg).

**Figure Ext4.** Log-scale test/OOD RMSE for (a) energy and (b) gradient. Points retain every neural seed; horizontal references are the frozen joint RBF and linear baseline, without refitting. [PNG](figures/extensions/Ext4_Generalization_english.png) · [SVG](figures/extensions/Ext4_Generalization_english.svg).

## 4. Separating learning error from reference-method bias

Matching saved predictions to new STO-3G calculations gives 126 rows at 21 identical geometries, including all eight original test points and five of the nine OOD points. No new fitting or quantum calculation is needed for this join. The exact signed identity is

\[
\widehat E-E_{\mathrm{FCI}}=
(\widehat E-E_{\mathrm{old\ RHF}})+
(E_{\mathrm{old\ RHF}}-E_{\mathrm{new\ RHF}})+
(E_{\mathrm{new\ RHF}}-E_{\mathrm{FCI}}).
\]

Both maximum reconstruction residual and old/new RHF replay drift are zero at the stored precision. At 2.70 Å, the RHF−FCI bias is 0.253643655 Hartree. For joint seed 7301, learning error is only 0.000423247 Hartree, yet total error from finite-basis FCI is 0.254066902 Hartree. A successful surrogate can therefore closely reproduce a poor dissociation reference. These signed components add; their RMSE values do not. Accidental cancellation would not demonstrate that the model learned electron correlation.

**Figure Ext7.** Three-seed mean signed learning error, RHF−FCI bias and total error from FCI at matched test/OOD geometries: (a) energy-only and (b) joint fitting, in mHartree. Connected points guide the eye; the plot does not add FCI supervision or extend evaluation to unmatched OOD points. [PNG](figures/extensions/Ext7_Error_Budget_english.png) · [SVG](figures/extensions/Ext7_Error_Budget_english.svg).

## 5. Anharmonic isotope models and the near-threshold failure analysis

For each electronic basis, the nuclear model uses

\[
V(r)=D_e[1-e^{-a(r-r_e)}]^2,\quad a=\sqrt{k/(2D_e)},\qquad
[-c_\mu\,d^2/dr^2+V(r)]u_v=E_vu_v,
\quad c_\mu=\hbar^2/(2\mu).
\]

The reduced mass uses bare proton/deuteron masses, 1.0072764665789/2.013553212544 u, from [NIST CODATA 2022](https://physics.nist.gov/cuu/Constants/Table/allascii.txt). Both isotopes share the same electronic potential. Energies are Hartree above the model minimum and distances Å; the kinetic coefficient is converted through SI units. Second-order central differences with zero endpoint wavefunctions produce a symmetric tridiagonal matrix, solved with [SciPy's selected-eigenvalue solver](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh_tridiagonal.html). This is J=0 radial finite-difference dynamics, not a high-order DVR or a full rovibrational spectrum.

An independent mathematical control sets Dₑ=0.2 Hartree, a=1.8 Å⁻¹ and rₑ=1 Å. The [Morse exact full-line spectrum](https://journals.aps.org/pr/abstract/10.1103/PhysRev.34.57) is

\[
E_v=w(v+\tfrac12)-\frac{w^2}{4D_e}(v+\tfrac12)^2,\quad
w=2a\sqrt{D_ec_\mu},\quad v+\tfrac12<2D_e/w.
\]

Here w is an energy spacing. Finite radial walls differ from full-line boundary conditions. Main grids use 400, 800 and 1200 intervals on [0,8] Å, with separate domain checks and a deliberately wall-affected control at rₑ=0.12 Å.

The FCI-derived Morse representation itself has 25-point electronic-curve RMSE of 0.00757282 Hartree for STO-3G and 0.00420654 Hartree for cc-pVDZ. Its parameters reproduce the chosen minimum, depth and curvature, not the whole FCI curve. These representation errors must be separated from nuclear grid errors.

| FCI-derived Morse model | Isotope | Harmonic wavenumber / cm⁻¹ | Numerical 0→1 gap / cm⁻¹ | Full-line Morse gap / cm⁻¹ | Numerical ZPE / Hartree | Numerical / full-line states |
|---|---|---:|---:|---:|---:|---:|
| STO-3G | H₂ | 5003.208 | 4722.097 | 4723.858 | 0.011236766 | 18 / 18 |
| STO-3G | D₂ | 3538.681 | 3397.084 | 3398.936 | 0.007979812 | 25 / 25 |
| cc-pVDZ | H₂ | 4383.737 | 4117.259 | 4118.591 | 0.009834132 | 16 / 17 |
| cc-pVDZ | D₂ | 3100.540 | 2966.494 | 2967.901 | 0.006986240 | 23 / 23 |

The table uses 1200 intervals on [0,8] Å. The harmonic H₂/D₂ ratio is 1.41386262; numerical anharmonic gap ratios are 1.39004412 and 1.38792116. These are level gaps, without intensities or selection-rule calculations, and are not assigned experimental IR lines. The model dissociation energy from its ground state is D₀=Dₑ−E₀.

Observed low-level convergence is second order, with ZPE orders 2.00057–2.00093. Main matrix residuals are at most 1.43×10⁻¹⁴ Hartree, while first-three-level errors remain 3.975–5.705 cm⁻¹ at the finest main grid. Richardson estimates reduce those discrepancies to 0.00491–0.01097 cm⁻¹, without certifying high states. Coarse discretization can also create a spurious near-threshold state: the analytic D₂ control returns 29 bound states at 400 intervals although the exact model has 28. The shifted-wall negative control loses three states for each isotope.

The main cc-pVDZ-derived H₂ spectrum instead misses its seventeenth weakly bound state. Its analytic binding is 6.69637×10⁻⁷ Hartree and amplitude decay length 15.0913 Å. Seven separate follow-up eigenproblems preserve the original main outputs:

| Right boundary / Å | Intervals | Returned / exact states | Highest-state binding relative error / % | Maximum shared-level error / cm⁻¹ |
|---:|---:|---:|---:|---:|
| 12 | 1800 | 16 / 17 | unavailable | 19.02409 |
| 24 | 3600 | 17 / 17 | 104.272 | 19.02409 |
| 48 | 7200 | 17 / 17 | 125.201 | 19.02409 |
| 48 | 14400 | 17 / 17 | 25.764 | 4.74622 |
| 48 | 28800 | 17 / 17 | 5.201 | 1.18595 |
| 96 | 57600 | 17 / 17 | 6.337 | 1.18595 |
| 96 | 115200 | 17 / 17 | 1.564 | 0.29645 |

Domain enlargement alone does not remove discretization bias. At 96 Å and 115200 intervals, all 17 states are recovered, but the weakest binding still differs by 1.564%. This remaining limit is retained rather than relabeled as complete convergence.

For 100, 298.15, 600, 1000, 2000 and 4000 K, finite bound-state sums give q=Σv exp(−Ev/kBT) and Fvib,bound=−kBT ln q. The main cc-pVDZ-derived H₂ value at 4000 K is 0.00643903 Hartree, 0.000615759 Hartree below the harmonic result. Its discrepancy from the full-line Morse bound sum is −4.09457×10⁻⁶ Hartree, reduced to −6.46456×10⁻⁸ Hartree in the finest tail calculation. Small partition-function differences do not certify individual high states. Continuum/dissociation, rotation, translation, nuclear-spin statistics and environmental terms are omitted; this bound-only F is not a complete molecular Gibbs energy.

**Figure Ext5.** The cc-pVDZ-derived Morse potential and first five numerical bound levels, with harmonic comparisons: (a) H₂ and (b) D₂. Energies are relative to the model minimum; the potential is a parameterization of FCI quantities. [PNG](figures/extensions/Ext5_Anharmonic_Levels_english.png) · [SVG](figures/extensions/Ext5_Anharmonic_Levels_english.svg).

**Figure Ext6.** (a) First-three-level maximum error versus interval count for the analytic H₂/D₂ Morse control on [0,8] Å, excluding the shifted-wall control; (b) bound-only minus harmonic vibrational F versus temperature for the cc-pVDZ-derived models. The separate long-tail study is tabulated above and is not plotted here. [PNG](figures/extensions/Ext6_Numerical_Limits_english.png) · [SVG](figures/extensions/Ext6_Numerical_Limits_english.svg).

## 6. Reproducibility, review and research limits

The [interactive explorer](quantum_explorer.html) exposes the saved numerical results; it performs no new quantum calculation. Full methods and ledgers are in the [correlation](../correlation_extension_notes.md), [learning](../learning_extension_notes.md) and [vibration](../vibration_extension_notes.md) notes. Outputs retain SHA-256 hashes, software versions, exact executed source snapshots, raw quantum logs, orbital diagnostics, checkpoints, predictions, wavefunctions and row counts. Quantum-source finalization and the later nuclear tail CLI changed the current scripts; historical snapshots identify what actually generated each dataset. A current-source hash is not substituted for an earlier executed version.

Twenty-one correlation, 14 learning, 18 vibration and four error-decomposition unit tests pass. Independent saved-data review additionally recomputed main quantum log energies and spin diagnostics, learning scales/metrics/ensembles and 107 nuclear-data checks, without repeating the science runs. Cross-model identities and hash checks support numerical consistency; they do not validate experimental chemistry. Final publication checks and visual-review scope are recorded in the linked validation and figure QA files, with separate [explorer QA](../results/extensions/explorer_qa.json) and [cross-review evidence](../results/extensions/cross_review.json).

To inspect the frozen release, use its saved-data validator. The following are portable command examples, not a claim that one wrapper executed the entire study. Existing chemistry/Psi4 interpreters and their native runtimes must already be configured; no installation is implied.

```powershell
$chem = $env:CHEM_PYTHON
$qm = $env:PSI4_PYTHON
$env:OMP_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
& $chem quantumequi/scripts/validate_extensions.py
& $chem -m unittest discover -s tests -p 'test_quantumequi_*extension.py' -v
```

Deliberate recomputation belongs in a separate checkout with an empty extension output area, while preserving the earlier frozen electronic reference and its split. Completed calculation directories reject overwrite. These current commands produce a new run; the corrected pilot need not reproduce historical pre-dispatch failures or require recovery. The old recovery ledger remains evidence of the original execution, not an extra mandatory step.

```powershell
& $qm quantumequi/scripts/correlation_extension.py --pilot
& $qm quantumequi/scripts/correlation_extension.py
& $chem quantumequi/scripts/learning_extension.py --pilot
& $chem quantumequi/scripts/learning_extension.py
& $chem quantumequi/scripts/error_budget_extension.py
& $chem quantumequi/scripts/vibration_extension.py --pilot
& $chem quantumequi/scripts/vibration_extension.py
& $chem quantumequi/scripts/vibration_extension.py --tail-followup
```

Recomputed outputs require fresh figures, validation and visual QA; cross-environment bitwise identity is not promised. This study establishes a reproducible small-system comparison with retained failures and error sources. Extension to organic electrocatalysis would still require representative electronic references, chemical-space validation, stationary-path verification and explicit environmental models. No DFT calculation, Cu-model training, electrode free energy, laboratory experiment or hardware validation is claimed here.
