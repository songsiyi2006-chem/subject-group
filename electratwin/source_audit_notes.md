# ElectraTwin submitted-program audit / 原程序数值与证据审计

This note concerns the submitted source and its one-line NumPy compatibility copy. Later reviewed implementations must be reported separately. All values below are executed model outputs or independent arithmetic/code diagnostics; none is a measured reaction result, a real device reading, or industrial qualification.

## 1. Source preservation and actual execution

- Attachment SHA-256: `152719e0fa9c4c8e1ac3316742629aa6cdcef56acaf9b784b8e53796219757ce`.
- Extracted original Python SHA-256: `41e1fd21c1b9c48b9e430f99750dd9f69ed5a97805d758b8d0249645d6c3ea35`.
- Compatibility Python SHA-256: `327bed061ee646babffff07ce593d914482fb8379ee0980254b699a1fa6fdbca`.
- Original execution: exit 1 after 6.10 s, `AttributeError` for `np.trapz`. The preserved docstring also emits an invalid-escape SyntaxWarning. The installed runtime was reused; no package installation was performed.
- Compatibility execution: exit 0 after 52.21 s on this CPU run, after replacing exactly one `np.trapz` call by `np.trapezoid`. This runtime covers the single case, twelve campaign cases, figures and exports; it is not a general speed benchmark. [NumPy's official release notes](https://numpy.org/doc/2.4/release/2.4.0-notes.html) support this API replacement.

`run_source.py` saves subprocess stdout/stderr and, after the compatibility main block returns, captures its PDE arrays and controller history with `runpy`. This observation layer does not change the model source. Logs redact local account/runtime paths. The source's institution, author, industrial status and connection claims are retained as quoted source output without endorsement.

## 2. PDE audit: a small iterate change does not establish conservation

The implemented interior stencil is

\[
u_j(C_{i,j}-C_{i-1,j})/\Delta x
=D(C_{i,j+1}-2C_{i,j}+C_{i,j-1})/\Delta y^2.
\]

The docstring's axial diffusion term is absent. Velocity is prescribed as a Poiseuille profile, temperature and diffusivity are fixed, and potential/charge transport, migration, heat transfer, pressure, gas and fouling are not solved. At fixed overpotential the boundary expression is affine in concentration; exponential dependence on a prescribed parameter does not make the coupled concentration problem nonlinear. The exchange-current density and transfer coefficient are assumed, not fitted to voltammetry. This is a simplified continuum transport model, not a first-principles molecular or validated multiphysics digital twin.

For the source case, the channel volume is 0.216 mL; Q = 450 µL/min and residence time V/Q = 28.8 s. The 100 × 45 calculation reports 973 iterations and an update norm of 9.989432 × 10⁻⁶. The independently recomputed normalized residual of its **implemented** interior stencil is 4.946831 × 10⁻⁶; the maximum dimensional residual is 0.011818608 mol m⁻³ s⁻¹. These values do not test the omitted axial-diffusion equation. The source returns `converged=True` unconditionally and clips reported conversion to 0–99.9%.

| Quantity | Reconstructed source value |
|---|---:|
| Mixed-cup conversion before rounding | 60.03536086% |
| Printed conversion | 60.04% |
| Integrated anode current before rounding | 0.02497501744 A |
| Printed current | 0.02498 A |
| Substrate removal from discrete inlet–outlet flux | 2.250163157 × 10⁻⁷ mol/s |
| Removal implied by current using source boundary n = 1 | 2.588478211 × 10⁻⁷ mol/s |
| Relative current/material-flux mismatch | +15.03513433% |
| First inlet-node trapezoidal current contribution | 0.003317589615 A |
| First inlet-node fraction of total integrated current | 13.28363283% |
| Discrete flow integration error relative to imposed Q | −0.05165289% |

The inlet row retains the uniform inlet concentration, including its anode corner, yet its current enters the trapezoidal electrode integral. That corner does not satisfy the downstream reactive boundary construction. Its contribution explains most of the observed conservation discrepancy; it is not acceptable to hide this by clipping conversion or FE. The audit deliberately retains both inconsistent flux estimates.

The mean-velocity Péclet numbers are Pe_L = 113636.36, Pe_H = 568.18 and Pe_Δx = 1147.84. A large mean axial value supports investigating a convection-dominated approximation, but it does not prove axial diffusion irrelevant at every no-slip boundary or inlet corner. [NIST's finite-volume description](https://www.ctcms.nist.gov/~wd15/fipy/documentation/numerical/discret.html) motivates auditing common face fluxes and global balances. Its [numerical-scheme documentation](https://www.ctcms.nist.gov/~wd15/fipy/documentation/numerical/scheme.html) also identifies numerical smearing from first-order upwind, so mesh refinement remains necessary.

## 3. Stoichiometry, green metrics and the fixed-grid campaign

The reactive boundary divides current by F, corresponding to one electron per substrate-consumption event in that balance. The metrics routine defaults to two electrons per product molecule. The same reaction and material basis must govern both. Using the very rounded current and conversion passed by the source controller, before its `np.clip` call, eight of twelve FE values exceed 100%.

| Flow, µL/min | η, V | Original printed FE, % | Reconstructed pre-clipping FE, % |
|---:|---:|---:|---:|
| 150 | 0.25 | 100.00 | 184.3997 |
| 150 | 0.45 | 100.00 | 151.7268 |
| 150 | 0.65 | 16.26 | 16.2607 |
| 500 | 0.25 | 100.00 | 184.7663 |
| 500 | 0.45 | 100.00 | 166.7540 |
| 500 | 0.65 | 31.04 | 31.0388 |
| 850 | 0.25 | 100.00 | 185.2421 |
| 850 | 0.45 | 100.00 | 169.8457 |
| 850 | 0.65 | 37.20 | 37.1967 |
| 1200 | 0.25 | 100.00 | 185.0425 |
| 1200 | 0.45 | 100.00 | 171.4308 |
| 1200 | 0.65 | 41.26 | 41.2596 |

The source chooses Q = 1200 µL/min, η = 0.45 V, conversion 32.03%, STY 38959.6 kg m⁻³ day⁻¹, E-factor 163.82 and SEC 0.275 kWh/kg as its compromise. These are outputs under inconsistent source assumptions; the printed 100% FE is actually 171.4308% before clipping. They must not be promoted as recommended experimental conditions or industrial performance.

The routine named EHVI executes a fixed four-by-three grid. Stub-solver probes with budgets 0, 1, 7 and 12 each produce 12 records. No surrogate, posterior uncertainty, feasibility model or acquisition optimization exists. The tested 2D hypervolume helper correctly gives area 6 for [(1,3),(2,2),(3,1),(1,1)] relative to (0,0); correctness of that rectangle union does not establish expected hypervolume improvement. The [BoTorch EHVI API](https://botorch.readthedocs.io/en/stable/acquisition.html) explicitly uses a predictive model in addition to a reference point and partitioning.

The source uses an assumed selectivity law 94 − 5η and an assumed terminal-voltage law 1.85 + η + 18.5I. No calibration data establish either law. Molecular weights are numerical constants without substrate/product structures or an atom-balanced reaction; importing RDKit does not constitute molecular-graph profiling, and no structure is processed.

The feed mass calculation omits electrolyte, workup, purification, water, cleaning, recovery and losses. It can only describe its explicitly narrow material boundary. [ACS GCI defines PMI using total inputs including water](https://www.acs.org/green-chemistry-sustainability/green-chemistry-nexus/articles/process-mass-intensity-calculation-tool.html). E = PMI − 1 follows from a common, closed mass balance when every non-product output is classed as waste; differing recycling or water conventions must be reconciled. [ACS's institutional guidance](https://www.acs.org/content/dam/acsorg/membership/white-papers/15lac-green-chemistry.pdf) defines E-factor as waste/product and explains why a single mass metric cannot establish overall environmental impact.

## 4. Instrument emulation is not a connected laboratory

The original SCPI class imports no transport driver, opens no connection and assigns its identity response from a string. `connect()` only changes a Boolean. The recorded 2.52309 V is synthetic output from 1.95 + 22.4I + Gaussian noise, which even differs from the controller's cell-voltage law. An unknown command receives an acknowledgement, an absurd current value is accepted in memory, and the voltage-limit attribute is never enforced. These are emulator diagnostics, not physical instrument tests.

[Tektronix's official reference-manual page](https://www.tek.com/tw/keithley-source-measure-units/keithley-smu-2400-series-sourcemeter-manual/model-2450-interactive-sou) confirms remote SCPI/TSP support for the 2450 family and the need to match documentation to the instrument. It does not verify this emulator's commands, the source's hardcoded identity or access to any named laboratory instrument. No actual output, pump, valve or electrochemical equipment was actuated.

## 5. Industrial rationale and evidence boundaries

Microchannel flow electrochemistry is a useful engineering approach; the attachment's assertion that it is the only industrial pathway is unsupported. A primary [OPRD study of a spinning-electrode reactor](https://pubs.acs.org/doi/10.1021/acs.oprd.3c00255) reports batch scale-up and continuous operation with a rotating-cylinder architecture. This is a counterexample to exclusivity, not proof that every batch design or this project is industrially validated. Publisher-indexed abstract/article excerpts were available; direct full-text and SI access was not obtained in this verification.

The institution, team, industrial partner, validated chemistry, calibrated transport/kinetics and real equipment remain unverified for this source. A physical digital twin would require a documented physical asset and reliable measurement linkage, identity and calibration records, model validation against independent experiments, operating limits and an explicit validation scope. Passing code tests or reproducing source figures does not supply those data.

## 6. Reproducibility and review scope

From the repository root, using a compatible scientific Python runtime:

```powershell
python electratwin/scripts/run_source.py
python electratwin/scripts/run_source.py --audit-only
```

The first command intentionally records the original NumPy failure before running the compatibility copy. The second reads saved arrays and performs small diagnostic probes; it does not launch another thirteen-case source PDE campaign. Runtime setup is documented by the repository's environment guidance. `results/source_audit.json` contains full precision, `results/sources.json` records the exact source-access scope, and original/compatibility directories retain actual logs and outputs. No source figures have been visually certified as publication-ready by this audit subtask.
