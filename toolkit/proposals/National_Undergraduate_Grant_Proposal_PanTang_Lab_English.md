# National Undergraduate Innovation Training Project Proposal

**Project: AI4S and Multiscale Microenvironment Modeling for Green Organic Electrocatalytic C–H Functionalization**

**Current implementation focus: trustworthy analytical quantification, data provenance and multiobjective condition recommendation.**

Version 2026-09-28. Complete content draft for supervisor and faculty review; not an approved project, submitted application or verified current official form. Applicant identities, supervisor comments and institutional signatures must be supplied truthfully by the responsible people.

| Administrative field | Proposed content and confirmation status |
|---|---|
| University | Guangxi Normal University |
| Faculty | School of Chemistry and Pharmaceutical Sciences; confirm actual enrollment and application rules |
| Category | Proposed innovation training project; national-level recommendation depends on institutional review |
| Student lead | Name, student number, program, year and contact: to be completed by the applicant |
| Members | Names, student numbers, roles and consent: unconfirmed |
| Proposed supervisor | Prof. Haitao Tang; willingness, eligibility and signature unconfirmed |
| Intended host team | Pan–Tang group; acceptance and resource access unconfirmed |
| Formal platform name | Confirm with the faculty; public institution pages contain a historical-name marker |
| Duration | Proposed 12 months from approval, T0 to T0+12 months; subject to the award document |
| Requested funds | Planning estimate CNY 10,000; not approved and not an instrument quotation |
| Call and deadline | Next open application window unverified |

## 1. Project abstract

Organic electrosynthesis couples current, passed charge, electrode configuration, solvent, electrolyte and transport. Connecting Bayesian optimization directly to unchecked analytical numbers can turn a dilution error, peak interference or inconsistent yield definition into a false optimum. This project proposes one supervisor-selected heterocyclic C–H functionalization case, first establishing HPLC and quantitative NMR records and quality checks, then studying tradeoffs among analytical yield, Faradaic efficiency and electricity per unit product.

The existing foundation is an executable software prototype and synthetic-data audit, not an experimental discovery. The prototype distinguishes supplied values, model predictions and measured evidence; links records by hashes; flags dimensional and material-balance anomalies; produces reproducible figures; and compares small-data models with simple baselines. The supplied example gives HPLC yield 1292.00%, incompatible with its 91.39% NMR example. A synthetic deconvolution test performs well with correct peak widths but reaches mean absolute area error 14.61% under a wrong-width assumption. These observations motivate improving analytical reliability before evaluating whether active learning saves effective experiments.

After training, reaction definition and calibration, the plan allows at most 40 independent reaction runs, retaining successful, failed and unquantified outcomes. Prospective holdouts and independent repeats will evaluate the model. Paired-state DFT and multiscale diffusion modeling are auxiliary studies conditional on evidence and resources. Deliverables are a traceable dataset, analytical workflow, reproducible code and undergraduate report; journal tier, awards and patent grants are not guaranteed outcomes.

## 2. Rationale and research context

### 2.1 Chemical question and value

Electricity provides controllable electron-transfer conditions, but replacing a redox reagent does not by itself establish overall environmental benefit. Electrodes, supporting electrolyte, solvent, purification, measured cell voltage and current efficiency influence material and energy demands. Leech and Lam discuss practical electrosynthesis variables, supporting careful condition records rather than proving applicability to the proposed substrate.[1]

The first research question is how current density, supporting-electrolyte concentration and temperature affect identified product yield and charge utilization when analytical error can be estimated. The second is whether an adaptive policy finds reproducible compromise conditions more efficiently than a predeclared comparator at the same reaction budget. A high model prediction alone cannot answer either question.

### 2.2 Methodological basis and intended contribution

Shields and colleagues demonstrate Bayesian reaction optimization in chemical synthesis.[2] We propose to examine sequential design in a defined case, without transferring performance claims from unrelated reactions. Small-data evaluation must prevent preprocessing leakage, related replicates crossing train/test boundaries, selective reporting of high-yield points and confusion between measurement noise and reducible function uncertainty.

Quantitative NMR depends on a characterized internal standard, correct proton counts and suitable integration conditions. BIPM resources and ACS guidance support attention to these factors.[3,4] Chemical shift, overlap and relaxation must be assessed for the actual sample. A C3-substituted indole cannot be quantified using the C3 proton that was replaced.

The proposed contribution is an evidence-controlled case study integrating these steps. A reaction-specific novelty search remains necessary after supervisor selection; this draft makes no claim of global priority or universally accurate site-selectivity prediction.

## 3. Objectives and testable questions

| Objective | Question or hypothesis | Evaluation |
|---|---|---|
| Reliable analysis | Preparation, calibration, standard and integration records support reproducible quantification | Calibration range, independent QC, recovery, repeats and complete raw records |
| Multiobjective optimization | An adaptive policy improves search efficiency at fixed chemistry and budget | Predeclared random/designed comparator, temporal holdouts and uncertainty |
| Chemical interpretation | Descriptors or transport scenarios explain some condition effects | Comparison against measured trends, retaining unsupported hypotheses |
| Student development | Students can perform the computational, recording and communication tasks of supervised research | Reproduction scripts, research log, data dictionary and progress presentations |

Proposed pilot engineering targets are independent-QC recovery of 95–105% and quantitative-repeat RSD no greater than 5%. These are project targets for deciding whether to enter model testing, not universal regulatory limits or already achieved performance. Adapt them to concentration, detection limits and method validation, then freeze them before prospective testing. If they fail, improve analysis before using unreliable labels for optimization.

## 4. Research work packages

### 4.1 Traceable HPLC and qNMR quantification

Create independent HPLC calibration records separating absolute peak area from internal-standard area ratios. Preserve units, intercept, slope, validated concentration interval, dilution chain, volume and method version. Import pre-integrated CSV values, flag anomalies and retain untrimmed results. NMR records must identify the product resonance, proton counts, standard purity/mass, sampled fraction and relaxation conditions. Missing FIDs cannot be presented as automatically parsed spectra.

Controlled synthetic-peak tests will assess numerical area recovery with known shapes and model error with misspecified widths. Real peak identity, baseline choice and coelution require analytical review. If vendor raw chromatograms or NMR FIDs are later supported, format conversion, phase correction and baseline processing need their own validation workflow.

### 4.2 Bounded multiobjective active learning

After defining the product and reaction, begin with a limited number of variables, such as current density, electrolyte concentration and temperature. Each record links molecular identities, measured charge, assay and batch. Standardize within training partitions; compare a mean or simple-regression baseline and a fixed design with GP recommendations. Define EI/UCB and observation-noise treatment before prospective evaluation.

Yield and FE are distinct. Compute FE using the balanced reaction's electron count and identified product amount. Calculate electricity demand from measured voltage/current integration and an explicit product-mass basis. Missing product molecular weight or voltage traces preclude a justified kWh/kg value. Use exact nondominance for a Pareto front; threshold filtering remains a feasibility rule.

### 4.3 Auxiliary computation and multiscale transport

Begin with a supervised neutral/radical-cation DFT pair for a defined substrate, documenting convergence, frequencies, spin and solvation. Discuss quantitative free energies only after species, electron/proton references and transition-state evidence are established. Existing input files do not supply a computed barrier.

POP single-atom coordination, hierarchical-pore diffusion and flow electrolysis remain conditional extensions. Start with parameter-explicit transport sensitivity and order-of-magnitude checks; structure–activity relationships and TOF need actual structural, diffusion and kinetic evidence. If instrumentation or material characterization is unavailable, complete the core analytical/optimization case without manufacturing materials claims.

## 5. Technical route and experimental design

```text
Reaction/product definition -> training and calibration -> supervised baseline
           |                          |                         |
Structure/charge/transport     Raw assay and calibration    Current/voltage/T logs
           +--------------------> versioned records <------------+
                                      |
                         Baselines and prospective holdouts
                                      |
                       Budget-matched EI/UCB and comparator
                                      |
                         Independent repeats and updates
```

The proposed ceiling is 40 independent reaction runs: 24 exploratory/training runs, 6 predeclared holdout conditions, 6 independent repeats and 4 chemical controls chosen by the supervisor. This is not a count obtained after deleting failures. Reinjections are technical replicates, not independent reactions. A comparison between recommendation policies must predeclare equal budgets and allocation; the same 40 runs cannot be counted as 40 runs for each policy.

The analytical estimate is 120 HPLC injections and 20 NMR samples. Two injections for each of 40 reactions use 80 injections; seven calibration levels in triplicate use 21; the remaining 19 cover blanks, QC and rechecks. Twenty NMR samples support selected orthogonal quantification and structural confirmation. Reinjection precision does not replace independent sample preparation or reaction repetition. The final concentration range, sampling plan and selection rules will follow the analytical pilot.

Report absolute errors, independent-repeat differences and measurement uncertainty, not correlation alone. Even r = 0.9821 for six mock scope labels has no measured predictive-validation meaning. If the budget is insufficient to establish an adaptive-policy advantage, report uncertainty or a negative result rather than extending selectively until a favorable result appears.

## 6. Proposed innovation and existing foundation

The intended methodological emphasis is analytical quality before multiobjective optimization, linkage of reaction records to raw assay/calibration/electrolysis data, and retention of model failures and assumption sensitivity. Novelty requires reaction-specific comparison with prior work. At present this is a proposed method combination and undergraduate training design.

| Existing item | Inspectable evidence | Unsupported conclusion |
|---|---|---|
| Reproducible bookkeeping and feedback prototype | Unit checks, held-out GP evaluation and input generation in the repository | A completed wet-lab loop |
| Analytical audit | HPLC 1292.00% versus NMR 91.39% discrepancy identified | HPLC/NMR concordance verified |
| Synthetic deconvolution | 21 calibration rows and 80 separation/width/seed scenarios | Real instruments or peak identities validated |
| Multiobjective graphics | 23 nondominated points among 80 synthetic conditions | A real reaction optimized |
| Figure and quantum templates | English/Chinese PNG/SVG and existing DFT input files | Calculated DFT barriers, experimental scope or accepted publication |

## 7. Laboratory resources and confirmation status

The university reported a Pan–Tang electrochemistry contribution in 2023, establishing thematic relevance but not access for this project.[5] Its shared-equipment site lists a 600M NMR instrument, which establishes a public listing, not current availability, fees, operating condition or ownership by the intended group.[6] The source's blanket claims of a 500 MHz spectrometer, HPLC-MS, multichannel electrochemistry, microflow systems and ample computing nodes are not carried into this draft as confirmed resources.

The university's research-institution directory marks the historical state-key-laboratory name as former, while another 2026 university page retains the older wording. The faculty must confirm the formal host name.[7] Required resources remain individually unconfirmed: compatible electrochemical equipment, HPLC method and standards, qNMR conditions, training, sample/waste procedures and CPU access. The current statistical work runs on CPU and establishes no GPU requirement. Without actual supervisor consent, support letters or funding commitments, the corresponding application fields remain unresolved.

## 8. Schedule and milestones

| Relative month | Work | Evidence required |
|---|---|---|
| M1–M2 | Define transformation, literature search, training and resource/budget confirmation | Product/reaction definition, executable protocol and resource record |
| M3–M4 | HPLC/qNMR calibration, independent QC, baseline and repeats | Analytical records, bias assessment and decision on modeling readiness |
| M5–M7 | Budgeted exploration and recommendations; archive every outcome | Data snapshot, baseline comparison and frozen holdout plan |
| M8–M9 | Independent repeats and holdouts; optional DFT/transport interpretation | Comparisons with uncertainty and records of unsupported hypotheses |
| M10–M11 | Data review, figures, reproduction and manuscript | Auditable methods/results and recorded revisions |
| M12 | Final report and undergraduate presentation | Data dictionary, code version, expenditure record and research reflection |

The source labels September 2026–May 2028 as two years. Those month starts differ by 20 months, covering 21 named calendar months rather than exactly 24 months. This draft uses a 12-month schedule relative to approval, without inventing a future call or backdating work.

## 9. Budget and assumptions

The estimate below totals CNY 10,000 and has not been quoted or approved. Funding source, reimbursable categories, discounts and supervisor contributions require university confirmation.

| Item | Quantity | Assumed unit CNY | Total CNY | Purpose |
|---|---:|---:|---:|---|
| Reagents and solvents | 1 lot | 2500 | 2500 | Defined substrate, solvents and electrolyte |
| Electrodes and cell consumables | 1 lot | 1000 | 1000 | Pilot electrodes, connections and replacements |
| Analytical standards and column consumables | 1 lot | 2000 | 2000 | Standards, calibration and analytical supplies |
| HPLC injections | 120 | 15 | 1800 | Calibration, blanks, samples and rechecks |
| NMR samples | 20 | 100 | 2000 | Quantification checks and identification |
| CPU pilot | 1 | 200 | 200 | Approved small calculations |
| Archiving/presentation | 1 | 500 | 500 | Storage, records and presentation material |
| Total | — | — | 10000 | Unapproved planning estimate |

If funding is insufficient, retain reaction definition, calibration, independent QC and core controls, narrow the variable/substrate scope and defer POP or microflow extensions. Do not trade evidence records for nominal sample count. Revise the scope and budget through the applicable university process before exceeding the plan.

## 10. Risks, adjustments and integrity

| Risk | Response |
|---|---|
| Invalid calibration or dilution chain | Retain values and flags; check raw records and reanalyze; never clip >100% to 100% |
| Coelution, assignment or relaxation bias | Check standards and orthogonal methods; a good fit does not establish identity |
| GP below baseline performance | Report it; retain comparator and examine noise/design; no promised improvement |
| Missing DFT or instrument resources | Narrow the case, retain software/analysis results and label unperformed tasks |
| Reaction fails | Review a literature-supported model system with the supervisor; retain failures and change rationale |
| Incomplete or duplicated records | Distinguish independent and technical replicates; prevent leakage and repeated use as independent validation |

Trained personnel supervise laboratory and instrument work under local procedures. This software does not control instruments, schedule experiments or replace operator judgment. Preserve raw data, calibration, exclusions, code and figure provenance. Synthetic and illustrative figures remain labeled and cannot serve as measured manuscript results.

## 11. Deliverables, roles and evaluation

Required outputs are a provenance-aware repository, HPLC/qNMR recording and checking workflow, source-backed multiobjective figures, a supervised chemical case or a clearly documented failure, a research report, progress log and final presentation. If actual findings and novelty justify it, the supervisor may consider publication, software registration or competitions. There is no promise of one or two top-quartile papers, a granted patent or a national prize.

The student lead would manage data, reproduction and writing; appropriately trained members would record samples and analyses under supervision; the supervisor and graduate mentor would review chemistry, instrumentation, assignments and computational choices. Confirm names, commitments and consent before filling roles. No signatures may be substituted, and a prospective supervisor is not a confirmed project leader. Evaluation should emphasize questions, methods, records, failure analysis and independent communication.

## 12. Application adaptation and references

The Ministry of Education's 2019 public management document emphasizes student practice and process-oriented training; current eligibility, supervision, dates and institutional procedures still require the relevant call.[8] Guangxi Normal University published proposed 2026 project recommendations on 3 June 2026. This work did not verify a subsequent open call or obtain its official application form, so the earlier notice cannot establish that the same round remains open.[9] This is transferable content, requiring actual identities, host name, current form, supervisor comments and authentic institutional endorsements before submission.

1. Leech M. C.; Lam K. *A practical guide to electrosynthesis*. Nature Reviews Chemistry 6, 275–286 (2022). [Publisher page](https://www.nature.com/articles/s41570-022-00372-y).
2. Shields B. J. et al. *Bayesian reaction optimization as a tool for chemical synthesis*. Nature 590, 89–96 (2021). [Publisher page](https://www.nature.com/articles/s41586-021-03213-y).
3. BIPM. [qNMR resources](https://www.bipm.org/en/organic-analysis/qnmr).
4. ACS, Journal of Medicinal Chemistry. [Purity by absolute qNMR](https://pubsapp.acs.org/paragonplus/submission/jmcmar/jmcmar_purity_instructions.pdf).
5. Guangxi Normal University, 2023-05-19. [Pan–Tang electrochemistry report](https://news.gxnu.edu.cn/2023/0519/c1330a267710/page.psp).
6. University shared-equipment platform. [600M NMR listing](https://dypt.gxnu.edu.cn/genee/equipment/61).
7. University [research-institution directory](https://www.gxnu.edu.cn/1382/list.psp) and [2026 visiting-scholar announcement](https://gxttc.gxnu.edu.cn/2026/0317/c10229a337773/page.htm).
8. Ministry of Education, 2019 document 13. [National undergraduate innovation and entrepreneurship training management rules](https://www.moe.gov.cn/srcsite/A08/s5672/201907/t20190724_392132.html). Official indexed text supported the document identity and principles; direct full-text access was restricted, so comprehensive verification of current clauses is not claimed.
9. University, 2026-06-03. [Project recommendation and completion announcement](https://www.gxnu.edu.cn/2026/0603/c1441a343392/page.htm).

Web verification date: 2026-09-28. Public information does not replace equipment access, supervision consent or faculty approval. See the [technical report](../reports/deployment_toolkit_report_english.md) and [numerical results](../results/audit/audit_summary.json).
