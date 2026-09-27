# Comprehensive Technical Report: Multiscale Machine Learning & Chemical Engineering for Green Organic Electrocatalysis

**Date:** 2026-09-27. **Prepared for:** the Pan–Tang research context, Guangxi Normal University. **Edition:** English; separate Chinese and section-aligned bilingual editions are supplied. The institutional affiliation in the source is a user-provided context, not an independently verified equipment inventory or authorship claim.

## 1. Abstract

The supplied three-module engine was executed unchanged and produced `research_grade_results.json`. Its 135-condition synthetic grid contains 12 nondominated records, representing 10 distinct objective vectors. The source's scalarization selects 15.71 mA/cm², 8.0 mol% catalyst and 900 rpm, with assigned yield 43.47%, Faradaic efficiency 74.49% and cell-only specific energy consumption 0.779 kWh/kg. These are model outputs, not measured performance. Eight supplied molecular structures were successfully embedded and minimized; their reported oxidation potentials are uncalibrated graph-descriptor scores and their yield labels are seeded normal draws. A first-order spherical diffusion model gives internal effectiveness factors of 0.2760 and 0.9033 at 200 µm radius for the two prescribed diffusivities. Its hard-coded statement that microporous particles above 50 µm necessarily have η below 0.40 is contradicted by its own numerical output.

Additional work comprises complete Pareto reconstruction, 231 weight combinations, three feasibility scenarios, five current-band holdout folds for a newly fitted multi-output Gaussian process, 64 converged molecular conformers, 16 actual GFN2-xTB jobs, 40 finite-volume mesh cases, 50 external-film cases and 60 transport-parameter cases. A 320-cell spherical finite-volume calculation agrees with the analytical effectiveness factor within 2.94×10⁻⁵. The η=0.40 microporous threshold is 129.379256 µm under the supplied assumptions. The contribution is a reproducible numerical and evidence audit, with a proposed one-year measurement program. No experimental substrate scope, accurate oxidation-potential predictor, validated single-atom catalyst, or journal-ready scientific discovery is established.

### 1.1 Evidence classes and scope

**[L]** primary literature or official method documentation; **[S]** executed synthetic objective or random-label calculation; **[M]** molecular mechanics/descriptors; **[Q]** executed semiempirical electronic structure; **[A]** independent code/numerical audit; **[P]** proposed experiment; **[U]** missing or unverified evidence. “Research-grade” is retained in source and output filenames as the brief's label; successful execution does not validate that claim.

The three source modules do not form a calibrated multiscale coupling. Molecular scores are not inputs to the green-optimization equations; no molecularly derived rate enters the POP model; product identity, electron stoichiometry and catalyst structure are not linked across modules. This report follows the requested five-section paper architecture while distinguishing what the equations calculate from what an experiment would need to establish.

### 1.2 Execution and provenance

The original script completed without a repair. The existing environment was reused: Python 3.12.14, NumPy 2.4.6, SciPy 1.18.0, pandas 2.3.3, scikit-learn 1.9.0, RDKit 2026.03.5 and xTB 6.7.1. Calculations were serial, CPU-only and single-threaded. The original output, logs and execution record are archived separately from all extensions. Additional scripts retain unrounded tables, identities, coordinates, charges and numerical diagnostics. Original source extraction normalizes line endings only.

```text
Extracted original script SHA-256:
@@SHA@@
```

The input supplies no measured training set, electrolysis trace, product structure, raw spectrum, catalyst synthesis or pore-transport measurement. Assigned values therefore remain assigned values even when their units are physical.

## 2. Multi-Objective Green Electrosynthesis (Module 1)

### 2.1 Model, units and actual algorithm

[S/A] The grid has 15 current densities from 5 to 35 mA/cm², three catalyst loadings (2, 5, 8 mol%) and three stirring rates (300, 600, 900 rpm). Comments describing 1–10 mol% and 200–1000 rpm are broader than the executed grid. All 135 conditions are evaluated by prescribed equations and compared pairwise. `GaussianProcessRegressor` and `Matern` are imported but never instantiated or fitted in the original script; no acquisition function, adaptive query, posterior uncertainty or outcome-feasibility constraint is implemented. The source JSON key `Module1_MultiObjective_BO` is thus a name, not evidence that Bayesian optimization occurred.

The source defines an assumed voltage and yield surface:

$$
U=2.1+0.045\ln(j+1)+0.012j,\qquad
v=\frac{j}{20}(1-e^{-\omega/350})(c/5)^{0.4},
$$

$$
Y=\mathrm{clip}\left[94\frac{v}{1+v}-0.08(j-18)^2,10,96.5\right],\quad
FE_{\%}=\mathrm{clip}[92-1.2j+0.0015\omega,15,94].
$$

Here j, c and ω are numerical values in the code's stated units; the logarithm contains an implicit current-density reference scale. Its coefficient is assigned, not fitted Tafel kinetics. With the stated resistance of 12 Ω, converting j to total current implicitly assumes 1 cm² electrode area. Alternatively, the term could be reinterpreted as an area-specific resistance, but that is not the source's stated unit.

$$
SEC=\frac{n_e F U}{3.6\times10^6 M f},\quad n_e=2,\quad
M=0.223\ \mathrm{kg\ mol^{-1}},\quad f=FE_{\%}/100.
$$

The implemented factor 3.6×10⁶ correctly converts joules to kWh with M in kg/mol. The docstring's factor 3600 would be wrong by 1000 with those same units. A check based on charge, electrolysis duration and product mass independently reproduces the implemented SEC. No product corresponding to M=223 g/mol is specified. The expression is cell electrical energy per assumed product mass at constant voltage and product-specific FE; it excludes stirring, pumps, heating, purification, solvent recovery and catalyst preparation. It is not a life-cycle greenness metric. The yield equation lacks reaction time and substrate inventory, so joint consistency of Y, FE and charge is not experimentally established.

### 2.2 Complete Pareto set and source compromise

Table 1 contains all 12 source-rounded nondominated records. A condition dominates another only when it is no worse in all three objectives and strictly better in at least one. Three records at 5 mA/cm² have identical objectives because yield is clipped at 10%; equality does not count as domination. There are therefore **10 unique objective vectors**. Independent unrounded enumeration gives the same membership and selected condition.

@@PARETO@@

The supplied “green compromise” is condition **53**, with Green_Score **4156.7141206675215** from rounded outputs. The unrounded score is 4154.82786378971. The highest-yield point instead uses 24.29 mA/cm² and returns 50.90% yield, 64.21% FE and 0.950 kWh/kg. Relative to that point, the source compromise sacrifices 7.43 yield percentage points, gains 10.28 FE points and lowers modeled SEC by about 18.0%. Neither is a laboratory recommendation.

The score $Y\,FE/SEC$ is a scalar preference, not a uniquely balanced optimum or an identified Pareto knee. Since $SEC\propto U/f$, the score is proportional to $Yf^2/U$ up to constants: FE is effectively counted twice. Catalyst loading and stirring are mostly pushed to their upper bounds because the generator rewards them without charging catalyst or mixing costs. The low-current clipped ties are the exception. The apparent green frontier inherits these modeling choices.

![Pareto objectives and scalar compromise](figures/pareto_tradeoffs@@SUFFIX@@.png)

Figure 1. Assigned response surface; color represents SEC. Stars identify the source scalarization. The connected frontier has overlapping identical objective records. Values are synthetic and uncertainty bars would not be measurement uncertainties.

### 2.3 Preference, constraints and geometry sensitivity

[A/S] A grid of 231 nonnegative weight triplets, summing to one in steps of 0.05, maximizes $w_Y\ln(Y/100)+w_F\ln(f)-w_E\ln(SEC)$ on the unrounded frontier. This dimensional scaling fixes the numerical reference units of SEC; multiplying its reference unit by a common constant does not change a given-weight ranking. Ties use the lowest condition ID. Ten conditions are selected across the grid; the source choice appears in **27/231** cases. Frequencies measure sensitivity to this arbitrary weight grid, not a probability of experimental optimality.

@@WEIGHTS@@

Three explicitly added outcome constraints illustrate how a recommendation changes when acceptable performance is stated. They are analyst-chosen scenarios, not a grant or industry standard. The third scenario is infeasible under this model; relaxing the requirements or changing chemistry would be necessary.

@@CONSTRAINTS@@

An electrode-area audit keeps 12 Ω fixed and substitutes $I=jA/1000$. It leaves the synthetic Y/FE equations unchanged solely to isolate the missing area term; this is not a scale-up prediction. Larger assumed area changes the preferred current density:

@@AREA@@

Real scale-up changes resistance, mass transfer, thermal behavior and geometry together. These results explain why area, spacing, conductivity and measured total cell voltage must be recorded before assigning a physical energy optimum.

### 2.4 Added multi-output GP and comparison with OVAT

[S/A] The extension actually fits a multi-output GP to the three unrounded synthetic objectives. Five contiguous current-band folds each hold out three of the 15 current levels: 108 training and 27 test conditions per fold. Feature scaling is fitted only on training rows. The GP uses an amplitude-scaled ARD Matérn-5/2 kernel, `alpha=1e-6`, training-target normalization, one optimizer restart and seed 42. Outputs share kernel hyperparameters; no learned cross-output covariance or mechanistic relationship is imposed. The outer folds test current interpolation/extrapolation within the same assigned surface, not new chemical reactions. [Official GP API](https://scikit-learn.org/stable/modules/generated/sklearn.gaussian_process.GaussianProcessRegressor.html).

@@GP@@

MAE/RMSE units are percentage points for Y and FE, and kWh/kg for SEC. Macro metrics average the five fold metrics. The GP improves substantially on a training-mean baseline for this smooth artificial problem. Nevertheless, stitching its out-of-fold predictions produces 15 predicted Pareto records, with precision **0.733333** and recall **0.916667** against the 12-record true model frontier. These are a retrospective diagnostic across five models, not a prospectively validated frontier from one deployed model. One kernel amplitude bound warning is archived; no tuning was performed to erase it.

OVAT can cheaply establish local operating trends but depends on its starting condition and may miss interactions. Exhaustive enumeration evaluates every modeled combination; sequential multiobjective BO would choose new measurements using predictive uncertainty and acquisition criteria. This delivery implements an added surrogate audit, not that sequential algorithm. No matched OVAT experiment was run, so no quantified experimental saving or superiority over OVAT is claimed. A future comparison should freeze equal budgets, feasible domains, reference-point choices and baseline procedures before outcomes are collected. Experimental BO literature provides precedent for that design, rather than evidence that this surface represents a reaction. [Shields et al., 2021](https://www.nature.com/articles/s41586-021-03213-y).

## 3. In Silico Substrate Scope Evaluation (Module 2)

### 3.1 Complete molecular scope table and identities

[M/S] Table 2 reproduces every source entry and its numerical output. Names are clarified where the SMILES is more specific: the benzofuran entry is the **ethyl ester**, not a generic carboxylate ion; the thiophene alcohol is 2-(thiophen-2-yl)ethanol. All eight saved structures are neutral with singlet input convention. Their inclusion does not establish that every compound is a natural product, has a specified bioactivity, or participates in one shared reaction.

@@SCOPE@@

*The source labels this number “V vs SCE,” but no reference-electrode calibration or training data support that unit assignment. **Yield values are random labels; they are neither assay nor isolated yields. Source categories are threshold labels of the same proxy, not observed reactivity. Atom indices are zero-based RDKit indices in the archived input order, not standard ring numbering. Table 2 is a descriptor and source-label inventory, not an experimental methodology-paper substrate scope.

![Eight supplied structures with atom indices](figures/substrate_atom_maps.png)

Figure 2. Red highlights mark each source-selected carbon. They do not assert product formation at that position. Geometry/identity files retain atom order and hydrogen connectivity. The image uses source names; the ethyl ester identity is clarified in Table 2.

### 3.2 Why geometry and partial charges do not validate the prediction

The source successfully embedded all eight structures with ETKDGv3 seed 42 and converged MMFF94 within 300 iterations, as independently recorded by the audit. However, the scoring path uses only calculated logP, topological polar surface area (TPSA) and the smallest Gasteiger charge among carbons bearing hydrogen:

$$
E_{proxy}=\mathrm{clip}[1.35+0.12\log P-0.008TPSA+1.5q_{min},0.75,2.30].
$$

These coefficients have no supplied fit, uncertainty, reference-electrode conversion or test set. There is no radical-cation calculation, orbital population, transition state or reaction partner in the source scoring function. Coordinate generation does not feed the formula. Repeating the charge calculation on the same molecular graph with no conformer gives a maximum charge difference of **0.0** for these inputs. All eight proxies lie inside the clipping limits, so clipping did not affect this particular table. [RDKit methods](https://www.rdkit.org/docs/RDKit_Book.html); [Gasteiger charge implementation](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/PartialCharges/GasteigerCharges.cpp).

@@LOCAL@@

Caffeine selects atom 0, an N-methyl **SP3 carbon** with positive partial charge **0.012958 e**, because it is the smallest charge in the restricted candidate set. Calling it an electron-rich aromatic SET site is unsupported. Other candidates are aromatic carbons, but that alone does not make their ranking mechanistic. SET changes the electronic state of the molecule; subsequent proton transfer, radical trapping, adsorption and competing barriers can determine site selectivity.

For a proxy below 1.45, the program draws $N(88,3^2)$; between 1.45 and 1.85 it draws $N(72,4^2)$; otherwise it draws $N(42,6^2)$. Source execution yields six “High” and two “Moderate” entries. The global seed is set by Module 1, so standalone Module 2 is not independently seeded. Reconstructing the seed-42 draws exactly reproduces all eight displayed yields. Repeating 1,000 seeds gives:

@@RANDOM@@

The different single-run yields within a category are random sampling, not substrate-specific predictions. The Gaussian generator is not intrinsically bounded to 0–100%; the stored 8,000 draws happened to remain within that interval. These distributions describe the code's noise assumption, not experimental uncertainty or chemical success rates.

### 3.3 Additional conformers and actual semiempirical calculations

[M/Q] Eight ETKDG conformers per molecule were generated with seed 271828 and optimized for up to 1,000 MMFF94 iterations. All **64/64 converged**. Duplicate conformational wells were retained; the lowest sampled MMFF energy was chosen as the starting structure for neutral GFN2-xTB/ALPB(acetonitrile) tight optimization. Each optimized neutral then underwent a +1 doublet single point at identical nuclei. All **16/16 jobs** terminated normally. Neutral/cation total atomic charges agree with 0/+1 within 10⁻⁶ e; no new dependency was installed.

@@XTB@@

Energy differences use 27.211386245988 eV/Eh. Each charged state uses its own equilibrium ALPB response; the fixed-nuclei difference is not a rigorous nonequilibrium-solvent vertical ionization energy and has no reference-electrode calibration. The source-selected caffeine carbon has a removal charge response of **−0.02830746 e**, despite the molecule losing one electron overall. Local charge redistribution can have either sign; it is not a local oxidation potential. These calculations add actual electronic-structure records, not validation of the earlier heuristic. [GFN2-xTB method](https://pubs.acs.org/doi/10.1021/acs.jctc.8b01176); [xTB options](https://xtb-docs.readthedocs.io/en/latest/commandline.html).

No DFT, charged-geometry relaxation, vibrational Hessian, thermal correction, measured Eox or competing reaction path was calculated. Neutral convergence does not certify a vibrational minimum. MMFF energies are compared only within each molecule; energy rankings across chemically different molecules have no reactivity meaning. The pilot does not establish all relevant conformations, protonation states or tautomer populations in an experimental solvent.

### 3.4 What a defensible scope predictor would require

[P/U] First define one transformation, its partners, product connectivity and measurable endpoint. Measure oxidation behavior on a common solvent/electrolyte/reference scale and retain irreversible peaks, fouling and unsuccessful runs. Separate oxidation potential prediction from product yield and regioselectivity prediction; each needs its own labels and validation. Eight structures without reactions cannot supply a universal yield predictor.

Calibrated chromatography should quantify conversion and product yield; isolated material and NMR connectivity should establish regioisomer identity. Freeze held-out substrate families before feature selection or tuning. Compare against simple baselines and report errors relative to analytical and between-day variability. A molecular or charge-response correlation is only a hypothesis until tested on independent chemistry.

## 4. POP-SAC Pore Kinetics & Thiele Modulus (Module 3)

### 4.1 Governing assumptions and analytical model

[S/A] The source assigns $k_{int}=4.2\times10^{-4}$ m³/(kg·s), $\rho=850$ kg/m³ and therefore $k_v=k_{int}\rho=0.357$ s⁻¹. It assigns $D_{eff}=1.5\times10^{-10}$ m²/s to the microporous label and $8.5\times10^{-9}$ m²/s to the hierarchical label, a ratio of 56.6667. No Cu structure, pore geometry, temperature, porosity, tortuosity or measured diffusivity is loaded.

For a uniform isothermal sphere, first-order reaction, constant effective diffusivity, zero central flux and fixed surface concentration:

$$
\frac{1}{r^2}\frac{d}{dr}\left(r^2D_{eff}\frac{dC}{dr}\right)-k_vC=0,
\qquad C(R)=C_s,\quad C'(0)=0,
$$

$$
\phi=R\sqrt{k_v/D_{eff}},\qquad
\eta=\frac{3}{\phi^2}(\phi\coth\phi-1).
$$

The source directly evaluates this analytical effectiveness factor. It imports `odeint` without using it and does not solve a new Knudsen/molecular diffusion model. η is the actual integrated rate divided by the hypothetical rate if the entire particle were at surface concentration. It is not conversion, selectivity or Faradaic efficiency. No intrinsic TOF or active-site count is supplied, so no absolute TOF is computed. The first-order spherical model is a standard reaction-engineering result. [MIT reaction/diffusion notes](https://ocw.mit.edu/courses/10-37-chemical-and-biological-reaction-engineering-spring-2007/resources/lec20_04272007_w/).

### 4.2 Full particle-size sweep and corrected design threshold

Table 3 preserves all ten source outputs. Sizes are **radii**, not diameters.

@@TRANSPORT@@

At 50 µm, microporous η is 0.7445; at 100 µm it remains 0.4890; only the 200 µm tested point lies below 0.40. Thus the hard-coded conclusion about all radii above 50 µm is false. Root solving the unrounded equation gives:

@@THRESHOLDS@@

Under these fixed assumptions, microporous particles cross η=0.40 at **129.379256 µm** radius; hierarchical particles cross η=0.85 at **260.666929 µm**. The source's hierarchical η>0.85 claim holds over the sampled 10–200 µm interval, not for arbitrarily large particles. At 200 µm, center/surface concentration is 0.00112958 for the microporous case and 0.76658717 for the hierarchical case. The model therefore predicts a strong internal gradient in the first scenario, not a demonstrated failure of a synthesized material.

![Internal utilization and concentration profiles](figures/pore_transport@@SUFFIX@@.png)

Figure 3. Left: analytical curves and source-radius markers; reference lines mark η=0.40 and its microporous threshold. Right: analytical concentration profiles and sampled 320-cell finite-volume values at R=200 µm. All parameters are prescribed rather than calibrated.

### 4.3 Independent finite-volume and external-film checks

[A] The extension discretizes the radial material balance in conservative spherical shells, with 40, 80, 160 and 320 cells for each of ten cases. Volume integration and surface flux are independently compared. Errors decrease with refinement; the maximum 320-cell η error is **2.9352×10⁻⁵**. An independent quadrature of $C/C_s=\sinh(\phi r/R)/[(r/R)\sinh\phi]$ reproduces the analytical η. This is numerical verification of the assumed equation, not validation of its parameters.

For small φ, the stable expression uses $\eta=1-\phi^2/15+2\phi^4/315-\phi^6/1575+\cdots$; for large φ, $\eta\sim3/\phi-3/\phi^2$. The source clips η below 0.01. That floor is inactive in its ten displayed cases but becomes unphysical outside them: at φ=10,000 the unclipped value is **0.00029997**, while the source returns **0.01**.

An added external-film scenario uses $Bi=k_fR/D_{eff}$ and steady flux balance:

$$
\frac{C_s}{C_b}=\left(1+\frac{\eta\phi^2}{3Bi}\right)^{-1},\qquad
\eta_{overall}=\eta\frac{C_s}{C_b}.
$$

@@EXTERNAL@@

The full file contains 50 radius/regime/Bi cases. These are hypothetical dimensionless conditions; equal Bi across different D does not imply equal physical film coefficients. Another 60 cases vary D or $k_v$ by factors of one-half/two and include $k_v\times0.4$ as an illustrative density-basis scenario. Neither set is a measured uncertainty interval. A high internal η can coexist with external transport loss.

### 4.4 Limits of a POP-SAC interpretation

The code labels 850 kg/m³ as skeletal density, whereas the volumetric reaction term requires a consistent catalyst mass per particle-volume basis. Porosity can make those densities different. The two assigned diffusivities cannot be derived from pore size alone. A gas-phase Knudsen treatment requires its own molecular-wall-collision assumptions; solvent-filled pores require consideration of liquid diffusion, partitioning, adsorption, steric hindrance and swelling. The source computes none of these. BET area does not measure liquid-phase effective diffusivity.

The approximately 15 nm pore width named in a comment is **mesoporous**; it does not by itself establish a macropore network or hierarchical connectivity. IUPAC terminology distinguishes micropores around/below 2 nm, mesopores between those and about 50 nm, and larger macropores. [IUPAC micropore](https://goldbook.iupac.org/terms/view/M03906/plain); [IUPAC mesopore](https://goldbook.iupac.org/terms/view/M03853/pdf).

An electroactive POP may also require electron/ion transport, potential distribution, site accessibility and charge-state dynamics. A uniform first-order sphere cannot identify a working SAC coordination site or separate all these limitations. A useful research question is whether measured apparent rates and concentration gradients across controlled particle sizes are consistent with a common independently constrained reaction/diffusion model.

## 5. Lab-Scale Translation & Grant Application Roadmap

### 5.1 Measurement sequence for electrosynthesis and substrate scope

[P/U] These steps are a proposed supervised program, not an executed SOP or verified local instrument inventory. A defined transformation, partner identities and a reproducible literature/group baseline are prerequisites. No target-product recipe is invented from the random scope labels.

1. Register canonical substrate/product structures, atom mapping, electron stoichiometry, batch identities and a common analytical endpoint. Establish calibration curves, internal-standard recovery, blanks and replicate preparation before optimization. Separate assay yield, isolated yield, conversion and selectivity.
2. For an initial analytical CV study, consider 1 mM analyte in MeCN with 0.10 M n-Bu4NPF6, a polished glassy-carbon working electrode, Pt counter electrode and compatible Ag/Ag+ reference. Establish the blank window first; collect 50, 100 and 200 mV/s traces, monitor fouling and document reference checks with ferrocene. These are proposed starting settings whose compatibility must be confirmed for the chosen reaction. Report peak/onset/midpoint distinctions rather than calling every number Eox.
3. Reproduce the chosen preparative baseline in independent runs. An undivided glassy-carbon/Pt cell can be a starting configuration when chemically appropriate. Measure immersed electrode area, gap, solution resistance, temperature and total cell voltage; retain current and voltage time traces. Check cathodic loss or crossover with an appropriate control before interpreting yields.
4. Calculate $FE=n_eFn_{product}/\int I\,dt$ and $SEC=\int U_{cell}I\,dt/(3.6\times10^6m_{product,kg})$ with a declared product basis. Log auxiliary power separately. Do not substitute anode overpotential for full cell voltage or an assumed 223 g/mol product for the identified product mass.
5. Freeze a feasible multiobjective study and matched OVAT/space-filling baselines under equal valid-evaluation budgets. Randomize execution where practical and interleave repeated reference conditions across days. State unacceptable potential, temperature, analytical failure and precipitation rules in advance. Failed chemistry remains data; instrument failure needs a separate status.
6. Transfer only a verified reaction to the eight-substrate panel. Keep the core procedure fixed or explicitly account for permitted changes. Determine product connectivity by 1H/13C NMR and HSQC/HMBC; use NOESY as supporting spatial evidence. HPLC-MS provides time-resolved offline composition unless a genuine online interface and latency are documented. Retain unsuccessful and mixed-regioisomer cases.

CV interpretation and electrolysis-cell design depend on the chemistry; they cannot be replaced by descriptor scores. [Electrosynthesis practice](https://www.nature.com/articles/s41570-022-00372-y).

### 5.2 POP characterization and continuous-flow translation

Prepare a defined support/metalation system before testing the abstract microporous/hierarchical labels. Measure particle-size distributions and distinguish radius, diameter and agglomerate size. Compare size fractions from the same parent batch where feasible, with independent-batch replication. Establish rate versus stirring and catalyst loading before attributing particle-size dependence solely to internal diffusion. Keep conversion low enough for an initial-rate interpretation and measure the relevant time interval rather than assuming first-order behavior.

Collect N2 physisorption at 77 K with material-compatible pretreatment, reported adsorption/desorption branches and a justified fitting range/pore model. Use ICP-MS after validated digestion for metal loading; obtain representative aberration-corrected HAADF-STEM fields for dispersion. Bright isolated features do not quantify all active atoms or identify operando coordination. Complement with XAS/XPS where accessible, leaching measurements, post-run structure and support/soluble-metal/nanoparticle controls. Gas sorption and dry microscopy cannot alone establish solvent-swollen transport or single-atom catalytic identity.

For a liquid-filled particle, estimate D by an independently designed uptake, tracer or kinetic method with partitioning accounted for. Fit $k_v$ and D only where the data can identify both; apparent rates alone can confound them. Reserve a particle size or batch for prediction. The source's 10–200 µm radii and 50 µm threshold are scenarios to test, not fabrication specifications.

For flow, measure residence-time distribution, wetted volume, pressure drop and product recovery. Start with a verified chemistry and a stable electrode/catalyst configuration; retain total cell energy and outlet composition until a prespecified steady-state criterion is met. A packed bed of POP particles is not automatically equivalent to an electrode film or the preceding microchannel model. Check external transfer, swelling, bed compaction, fouling and catalyst loss; include separation and auxiliary energy in any scale-up comparison. Demonstrate mass and charge closure before claiming a productivity advantage.

### 5.3 One-year execution plan and acceptance gates

Months refer to project start, not the current university application calendar. A freshman should prioritize one falsifiable core question and use the other modules as support. Parallel ambition does not remove instrument, supervision or coursework constraints.

![One-year project Gantt](figures/one_year_roadmap@@SUFFIX@@.png)

| Period | Core deliverable | Advancement gate |
|---|---|---|
| Months 1–3 | Literature/structure registry, training, analytical calibration, replicated baseline | Defined product, traceable records and quantification repeatability |
| Months 4–6 | Equal-budget pilot design; current/voltage balances; measured descriptor table | Metrics have consistent units and uncertainty; failures are recorded |
| Months 7–9 | Prospective optimization or controlled particle-size study; substrate validation as justified | Independent predictions improve a declared baseline or expose a useful limitation |
| Months 10–12 | Independent reproduction, negative-result analysis, grant report and manuscript draft | Raw data, code, uncertainty and claim boundaries survive internal review |

Stop or redesign if product identity is unresolved, analytical drift exceeds the proposed effect, optimizer gains vanish on repetition, transport parameters are not identifiable, or separation/electrode deterioration removes the energy benefit. A null or negative result can still satisfy a sound methodological milestone; it does not justify a chemical breakthrough claim.

### 5.4 Undergraduate innovation proposal and publication position

A defensible proposal can ask: **Can a measured multiobjective protocol lower product-specific electrical energy without sacrificing a predeclared yield/FE threshold, and can particle-scale transport explain any loss of catalyst utilization?** Choose the electrosynthesis or particle-transport question as the core after baseline feasibility; combining all eight substrate classes and all material regimes in one year is not a prerequisite.

The application dossier should contain a literature novelty comparison, a falsifiable hypothesis, measured preliminary data clearly separated from this simulation, work packages, instrument access, itemized consumables/analysis costs, responsibilities and backup plans. Use the university's current National Undergraduate Innovation and Entrepreneurship Training Program notice for eligibility, forms and deadlines; none was supplied or verified here. The schedule supports preparation and execution, not a guarantee of national approval, an award or acceptance in JACS/Angewandte/ACS Catalysis/Green Chemistry.

For a paper, computational integrity is necessary but insufficient. A publishable claim would need, for example, a prospective experimentally reproducible optimization gain, a reference-consistent substrate model with independent tests, or a constrained reaction/transport mechanism that predicts new particle-size behavior. The present release supplies reproducible software and an audit appendix; those missing chemical results remain the next work.

### 5.5 Reproduction, data availability and primary references

The release includes the original specification/script/output, the independent audit, all 135 conditions and 12 Pareto rows, per-fold predictions, eight molecular identities, 64-conformer SDF records, 16 xTB jobs, diffusion profiles, bilingual figures and three complete report editions. [Research package guide](../README.md) gives paths and environment requirements. The previous four-task and five-topic reports remain separate.

```text
python research/scripts/run_original.py
python research/scripts/audit_research.py
python research/scripts/plot_research.py --zh-font PATH_TO_CJK_FONT
python research/scripts/build_reports.py
python research/scripts/validate_research.py
python -m unittest discover -s tests -v
python scripts/validate_release.py
```

Run in a separate clone to preserve release outputs: these commands refresh saved files. Activate the existing numerical environment, set OMP/OPENBLAS/MKL threads to one, and set `XTB_EXE` when xTB is not on PATH. Each xTB job has a 180-second timeout. Public logs normalize machine-specific executable/directory paths without changing scientific values. Numerical tests validate records and equations, not chemical accuracy. No external measured training set or supporting-information experiment was reproduced.

1. Shields et al. (2021), *Bayesian reaction optimization as a tool for chemical synthesis*. [Nature](https://www.nature.com/articles/s41586-021-03213-y). Experimental design precedent only.
2. scikit-learn developers, [GaussianProcessRegressor API](https://scikit-learn.org/stable/modules/generated/sklearn.gaussian_process.GaussianProcessRegressor.html). Multi-target fitting and normalization; recorded local version governs reproduction.
3. RDKit developers, [RDKit Book](https://www.rdkit.org/docs/RDKit_Book.html) and [Gasteiger implementation](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/PartialCharges/GasteigerCharges.cpp). Conformers and graph-based descriptor provenance.
4. Bannwarth, Ehlert and Grimme (2019), GFN2-xTB. [Primary method](https://pubs.acs.org/doi/10.1021/acs.jctc.8b01176); [official xTB documentation](https://xtb-docs.readthedocs.io/en/latest/commandline.html). No target-specific calibration is implied.
5. MIT OpenCourseWare, [Reaction and Diffusion in Porous Catalysts](https://ocw.mit.edu/courses/10-37-chemical-and-biological-reaction-engineering-spring-2007/resources/lec20_04272007_w/). Official teaching source for the classical assumptions. Thiele's 1939 paper has DOI [10.1021/ie50355a027](https://doi.org/10.1021/ie50355a027); full publisher text was not retrievable in this run.
6. IUPAC Gold Book, [micropore](https://goldbook.iupac.org/terms/view/M03906/plain) and [mesopore](https://goldbook.iupac.org/terms/view/M03853/pdf). Pore terminology, not diffusivity measurements.
7. Leech and Lam (2022), [electrosynthesis methodology](https://www.nature.com/articles/s41570-022-00372-y). General measurement context; the proposed sequence above is newly designed, not copied from a group SOP.

**Conclusion.** The supplied calculations are reproducible, and their unsupported mechanistic and material-design claims can be identified quantitatively. The added models improve numerical transparency. Progress to scientific validation requires a defined reaction and independent measurements, not a stronger label for the same synthetic targets.
