# Editorial assessment and bounded contribution statements / 编辑评估与有界贡献陈述

Assessment date / 评估日期: 28 September 2026  
Manuscript author as supplied / 用户提供的稿件作者: Siyi Song / 宋思毅, Guangxi Normal University / 广西师范大学  
Frozen study / 冻结研究: repository commit `1ad05c243155a9f64b18e0d58c665424fde919a0`

## English assessment

### Evidence and editorial position

The present evidence supports a reproducible computational-methodology case study: two controlled H₂ calculation chains, accompanied by selected verification failures from other workflows. It does not yet establish a generally superior method, a broadly representative benchmark, an electrocatalytic mechanism or a new synthetic discovery. Its strongest candidate contribution is **quantified loss of improvement when the reference or requested observable changes**, supported by retained intermediate calculations and negative controls. High-tier publication readiness is not established by the current evidence; this is an assessment of the demonstrated scope, not a prediction of an editor's decision.

The two chains must remain distinct. The learning chain fits RHF labels and evaluates signed errors against fresh RHF and same-basis FCI at identical geometries. The nuclear chain parameterizes a Morse model from FCI equilibrium quantities and assesses its finite-difference spectrum. The trained neural potential has **not** been propagated into a nuclear spectrum. There is therefore no executed neural-potential-to-ab-initio-spectrum error budget. The FCI/Morse parameterization is also not a global fit to the electronic curve. The existing [main report](../quantumequi/reports/quantumequi_report_english.md), [extension report](../quantumequi/reports/extension_report_english.md) and [post hoc analysis](results/analysis_summary.json) preserve these boundaries.

Established work already addresses force-field evaluation beyond prediction errors, uncertainty assessment, leakage and the dependence of downstream spectra on fitted surfaces. These principles cannot be claimed as new. The useful distinction here is the combination of paired reference changes, signed attribution and observable-specific numerical counterexamples within a replayable release. [Fu et al., R03](https://openreview.net/forum?id=A8pqQipwkt), [Scalia et al., R04](https://doi.org/10.1021/acs.jcim.9b00975), [Kapoor and Narayanan, R05](https://doi.org/10.1016/j.patter.2023.100804), [Kamath et al., R07](https://doi.org/10.1063/1.5003074).

### Four falsifiable contribution statements

**C1. Reference-limited transfer of an observed learning improvement.**  
“At eight matched H₂ test geometries, adding gradient supervision to three paired neural fits reduced energy RMSE against their RHF labels by 60.240–73.649%, but reduced RMSE against same-basis FCI by only 0.00838–0.02168%. A signed decomposition identifies the unchanged RHF–FCI discrepancy and its cross term with fitting error.”

This statement can be falsified directly from [paired_gain_transfer.csv](results/paired_gain_transfer.csv) and [mse_attribution.csv](results/mse_attribution.csv). The five matched OOD geometries are a separate comparison: 62.528–97.445% reduction against RHF versus 3.403–8.214% against FCI. They are not the original nine-point OOD set. Three initializations describe these runs; they do not provide independent samples of molecular chemical space or population confidence intervals. “Gradient supervision learns correlation” would be false: no FCI label enters training. The signed identity and triangle bounds are elementary algebra; the potential contribution is the controlled quantitative example and its decision value. The exact reference is FCI **within STO-3G**, not experiment or a complete-basis limit.

The mean-square relation is
\[
\mathrm{MSE}(\epsilon+b)=\langle\epsilon^2\rangle+\langle b^2\rangle+2\langle\epsilon b\rangle.
\]
Its terms are not independent positive variance shares. A negative cross term represents cancellation; percentages of these terms must not be presented as probabilities or additive fractions of causation. Selecting which reference to improve next is a conditional conclusion for the measured geometries.

**C2. Paired mathematical success and failure of the intended target.**  
“The same trained H₂ scalar-energy functions pass the recorded symmetry and energy-derivative checks while retaining substantial error relative to the correlated electronic comparator. Separately, a source random potential passes consistency checks but cannot justify a chemically meaningful stationary point.”

The falsification test is to join the saved checkpoint identity, transformation/finite-difference records and target-error records, rather than comparing unrelated models. For the trained chain, the reported maximum transformation discrepancy is \(2.22\times10^{-16}\) in its corresponding units; the finite-difference force discrepancy is \(2.22813\times10^{-10}\) Hartree Å\(^{-1}\) at \(h=10^{-5}\) Å. Joint test FCI energy RMSE nevertheless remains approximately 0.06864 Hartree. These are logically compatible measurements. For the untrained source, force units are nominal and the reported path candidate is nonstationary; this separate audit cannot be described as trained-model performance. Equivariance, differentiation and successful linear algebra are established requirements, not sufficient chemical validation. [Unke et al., R01](https://doi.org/10.1021/acs.chemrev.0c01111), [Batzner et al., R02](https://doi.org/10.1038/s41467-022-29939-5).

**C3. Observable-dependent convergence with both omitted and spurious states.**  
“A controlled Morse nuclear benchmark separates potential representation, grid spacing and finite-domain effects. Machine-small matrix residuals coexist with a missing near-threshold state; coarse D₂ grids can also return extra below-threshold states. Improving a low-level or partition-sum metric does not certify the shallowest binding energy.”

The cc-pVDZ-parameterized H₂ model has 17 exact full-line Morse bound states; the 8 Å principal domain returns 16. Its shallowest analytic binding is \(6.6963668574\times10^{-7}\) Hartree. A saved 96 Å/115200-interval spectrum recovers the state but still has a +1.56403% binding error. Standard second-order, fixed-domain Richardson arithmetic on two saved grids yields −0.0268823%; it is a **post hoc estimate**, not another eigenproblem, a replacement for the single-grid result or an independently verified convergence order at 96 Å. The [saved comparison](results/residual_observable_comparison.csv) even contains a pair for which a larger matrix residual accompanies a smaller binding error.

The negative 400-interval D₂ controls return 29 rather than 28 states for the analytic control, 26 rather than 25 for the STO-3G model, and 24 rather than 23 for the cc-pVDZ model when compared with the corresponding full-line analytic counts. Those are discretization warnings, not new physical states. The full-line comparator and the radial Dirichlet problem have distinct boundaries. Finally, the Morse representation already differs from the 25 electronic FCI points: curve RMSE is 0.00757282 Hartree for STO-3G and 0.00420654 Hartree for cc-pVDZ. No solver refinement removes that model discrepancy. The method and analytic reference are established. [Morse, R11](https://doi.org/10.1103/PhysRev.34.57).

**C4. A traceable collection of non-equivalent validation targets.**  
“The supporting workflows expose distinct failures that would be missed by treating execution, conservation, fit quality, ensemble mean or optimization completion as a common certificate.”

Concrete examples include fixed-iteration NEB completion without force convergence; correct stochastic balances conditional on an assumed network; temperature-mean agreement without the intended variance; hydration labels mistaken for redox relevance; and analytical recovery sensitive to peak-shape misspecification. The current contribution is a documented collection and a way to map each claim to its appropriate comparator. It is **not** an estimate of failure prevalence in published software or the effectiveness of an automatic detector. Existing NEB, stochastic simulation, thermostat and database methods retain their original attribution. [Henkelman et al., R12](https://doi.org/10.1063/1.1329672), [Gillespie, R14](https://doi.org/10.1021/j100540a008), [Berendsen et al., R15](https://doi.org/10.1063/1.448118), [FreeSolv, R16](https://doi.org/10.1007/s10822-014-9747-x).

C1 and C3 are the strongest quantitative claims. C2 connects the diagnostics. C4 should remain supporting evidence unless externally assembled cases establish generality.

### Closest prior work and the remaining difference

| Prior work | What is already established | What the current data add, and what remains absent |
|---|---|---|
| [Fu et al., TMLR 2023, R03](https://openreview.net/forum?id=A8pqQipwkt) | ML force-field comparisons using molecular-simulation outcomes, beyond local prediction-error metrics. | The present release adds a small signed RHF/FCI reference-change example and numerical threshold-state controls. It lacks comparable chemical breadth and trained reactive-dynamics validation. |
| [Kamath et al., JCP 2018, R07](https://doi.org/10.1063/1.5003074) | NN and GP fitted to matched formaldehyde samples are evaluated through actual vibrational spectra on common numerical settings. | The current finite-domain/state-count controls are inspectable, but the two chains do not yet execute NN-to-spectrum propagation. “First downstream spectral evaluation” is excluded. |
| [Scalia et al., JCIM 2020, R04](https://doi.org/10.1021/acs.jcim.9b00975) | Explicit comparison of scalable molecular-prediction uncertainty methods. | Three-seed spread and observed undercoverage are negative diagnostics, not a new uncertainty method or an adequately calibrated coverage study. |
| [Joeres et al., Nature Communications 2025/2026, R06/R25](https://doi.org/10.1038/s41467-025-58606-8) | Similarity-aware splitting with a formal objective and empirical comparisons; the [addendum](https://doi.org/10.1038/s41467-025-67495-w) additionally evaluates supplied dataset splits. | Hashes and data-access records help audit this release but do not improve the splitting algorithm. The addendum's PINDER counterexample precludes a claim that one split method is universally best. |
| [dos Santos et al., JCIM 2025, R08](https://doi.org/10.1021/acs.jcim.4c01847) | Morse-based augmentation improves reactive-MLIP dissociation curves and BDEs, including a methane-combustion study. | The present work diagnoses reference and Morse-representation limits; it does not propose a competing training augmentation or demonstrate improved chemical dissociation predictions. |
| [MoleculeNet, Chemical Science 2018, R17](https://doi.org/10.1039/C7SC02664A) | Broad molecular tasks, baseline comparisons, datasets and task-dependent splits. | A 256-molecule FreeSolv subset provides supporting baseline counterexamples. It is far smaller and does not support chemical-space-wide model ranking. |
| [Zhang et al., Nature Chemical Engineering 2026, R18](https://doi.org/10.1038/s44286-026-00392-1) | Knowledge-graph/agent-supported chemical-process digital-twin construction and calibration. | Retained simulation provenance is useful, but the repository's uncalibrated process models and simulated interfaces cannot establish industrial twin validity. |

The comparison is a targeted literature assessment, not an exhaustive novelty search. “First,” “universal,” “state of the art” and “top-journal ready” are unsupported.

### Tang-group experimental context and journal classification

Four organic-electrosynthesis records were checked against publisher author lists and the affiliation of **Hai-Tao Tang at Guangxi Normal University**, not inferred from the name alone:

- [R19, JACS 2026, DOI 10.1021/jacs.6c07265](https://doi.org/10.1021/jacs.6c07265): pyridine hydroxymethylation using a zinc single-atom cathode and iodide relay. Tang's publisher entry includes the GXNU affiliation, institutional email and ORCID 0000-0001-7531-0458.
- [R20, Angewandte Chemie International Edition 2024, DOI 10.1002/anie.202404295](https://doi.org/10.1002/anie.202404295): an iron single-atom redox mediator for organic anodic transformations.
- [R21, Angewandte Chemie International Edition 2024, DOI 10.1002/anie.202315032](https://doi.org/10.1002/anie.202315032): manganese single-atom electrocatalysis of silane oxidation. Online publication is 6 December 2023; the journal volume year is 2024.
- [R22, Green Chemistry 2023, DOI 10.1039/D3GC01288C](https://doi.org/10.1039/D3GC01288C): a multicomponent electrosynthesis of tellurium-containing oxazolidinones. Tang is a GXNU coauthor; this record does not designate him as a corresponding author.

These studies establish relevant experimental chemistry, not calibration of the repository's fragment, rate constants, synthetic yields or biological activity. They do not establish an advisory relationship or collaboration with Siyi Song. The author's affiliation is user supplied.

The [reference registry](references_verified.json) separates publication verification from journal classification. No year-specific authoritative CAS major-category or JCR quartile record was verified during this assessment; affected entries remain `classification="unverified"`. This does **not** mean those journals are outside Q1. It means the present file does not certify that property. A journal's reputation, an impact factor, a Scopus quartile or a publisher list without an explicit ranking year/system must not be relabelled as CAS Zone 1 or JCR Q1. Foundational methods, the published TMLR evaluation paper and the ICML message-passing paper are included for relevance, not to manufacture an all-Q1 bibliography.

The registry includes the original [Psi4](https://doi.org/10.1063/5.0006002), [Morse](https://doi.org/10.1103/PhysRev.34.57), [CI-NEB](https://doi.org/10.1063/1.1329672), [basis-family](https://doi.org/10.1063/1.456153), [MPNN](https://proceedings.mlr.press/v70/gilmer17a.html) and [SchNet](https://doi.org/10.1063/1.5019779) attributions. Partial author lists and unavailable full texts are explicitly marked. The tangent-paper author site lists an erratum whose text was inaccessible; no claim about that erratum's content is made.

### A concrete route to a stronger methods paper

The following is **proposed work**, not an account of completed calculations. More jobs alone would not resolve the central weaknesses.

1. **Make the validation procedure independently testable.** Freeze a claim schema containing observable, units, reference Hamiltonian/basis, domain, data access, convergence criterion and allowed interpretation. Assemble both correctly specified controls and independently designed failure cases before running an automatic audit. Compare execution-only checks, symmetry/derivative checks and the full observable-specific procedure. Report detected failures, missed failures and false alarms by failure family, not a single opaque score. A blinded external maintainer should reproduce at least one complete chain from a clean environment.
2. **Test transfer across chemically different systems.** Predeclare a small, tractable set of diatomics, polyatomic vibrations and reaction paths with distinct electronic difficulties and reliable references. Hold out molecular identities or entire trajectories, not adjacent geometries alone. Compare a simple interpolator/RBF baseline and at least two learned representations with matched label access and compute budgets. Additional seeds quantify optimization variability; they do not replace new chemical systems.
3. **Close the currently missing potential-to-spectrum link.** For H₂, compute or reuse a sufficiently resolved correlated PES, assess a basis ladder and construct independently checked interpolants. Propagate the reference PES, parameterized Morse model and trained energy models through the same converged nuclear solver and boundary study. Separate electronic reference error, interpolation/representation error and nuclear discretization. First validate J=0 bound energies; only then add rotation, continuum or nonadiabatic corrections needed for a specific observable. Compare experimental spectra only after matching the Hamiltonian and corrections.
4. **Calibrate uncertainty for the chosen use case.** Use a distinct calibration set and untouched identity/domain-shift test sets. Report coverage, interval width, error–uncertainty association and a proper scoring rule against the declared reference. Retain the difference between uncertainty in fitting the label and bias of the label itself. The current eight-distance/three-seed example is insufficient for reliable general coverage estimates.
5. **Require chemistry-specific evidence before mechanistic claims.** A reaction-path study needs chemically complete endpoints, specified charge/spin, converged stationary points, a first-order saddle test and endpoint-connectivity checks on a calibrated electronic potential. Electrosynthesis additionally needs explicit electrode/solvent/reference-potential conventions and comparisons to the relevant experimental observable. The four Tang papers motivate possible targets; their reported systems cannot be replaced by the source hand fragment.
6. **Preserve prospective and retrospective separation.** The current gain-transfer and Richardson analyses are exploratory post hoc analyses. Freeze thresholds and primary outcomes before a new external study; retain negative results and compute accounting, including cached evaluations and rejected pre-driver configurations. Additional software tests support implementation reliability but cannot substitute for the above reference comparisons.

A defensible methods submission can present the current release as a transparent, limited case study now. A stronger claim of reusable scientific methodology requires demonstrable detection value, independent adoption or replication, and chemical breadth beyond H₂. No publication probability is assigned.

## 中文评估

### 证据范围与论文定位

现有结果足以支持一篇可复现的计算方法案例研究：以两条受控H₂计算链为主线，辅以其他工作流的具体验证失败。当前最值得发展的贡献是：**定量说明，在改变参考方法或最终目标物理量后，前一层计算的精度提升可能大幅衰减。** 现有证据尚未建立普适方法优势、具有代表性的广泛化学基准、电催化机理或新合成发现，也不能据此认定已达到顶刊标准。这是对证据的判断，不是对编辑决定或发表概率的预测。

两条主链需要明确区分。学习链以RHF数据训练，并在完全相同的几何位置，以新算RHF和同基组FCI评价带符号误差。核振动链由FCI平衡位置、势阱深度与曲率参数化Morse势，再研究有限差分振动谱。**训练得到的神经网络势尚未输入核振动求解器。** 因此，不能写成已完成“神经网络势到从头算振动谱”的端到端误差传播。Morse模型也不是对整条FCI势能曲线的全局拟合。

力场预测误差与模拟可靠性不等价、不确定度需要校准、数据泄漏需要约束，以及势能拟合影响振动谱，均已有发表研究。本稿应把区别落在同一批可追溯记录上的配对参考变化、保留符号的误差归因，以及针对目标物理量的数值反例。[Fu等，R03](https://openreview.net/forum?id=A8pqQipwkt)、[Scalia等，R04](https://doi.org/10.1021/acs.jcim.9b00975)、[Kapoor与Narayanan，R05](https://doi.org/10.1016/j.patter.2023.100804)、[Kamath等，R07](https://doi.org/10.1063/1.5003074)。

### 四项可证伪的贡献表述

**C1：学习改进向高阶电子参考的传递受限。** 在8个相同的H₂测试几何上，3组配对神经网络加入梯度监督后，相对RHF标签的能量RMSE下降60.240–73.649%，而相对同基组FCI的RMSE仅下降0.00838–0.02168%。5个匹配OOD几何的对应范围为62.528–97.445%和3.403–8.214%，不能混同原始9点OOD集合。读者可以直接用保存的逐行预测与参考能量推翻或重现此结论。

误差分解及三角不等式属于常规数学；潜在贡献是受控的定量结果和计算资源分配启示。3个随机种子不是3类独立化学体系；FCI只在给定有限基组和哈密顿量内精确。MSE交叉项可以为负，不能画成非负的“贡献百分比”。误差抵消也不能解释为网络学会了电子相关，因为FCI标签没有进入训练。

**C2：同一模型可通过数学一致性检查，同时未达到目标参考精度。** 保存的训练后H₂模型通过对称性及能量—力导数检查：最大变换差为相应单位下的\(2.22\times10^{-16}\)，步长\(10^{-5}\) Å处的有限差分力误差为\(2.22813\times10^{-10}\) Hartree Å\(^{-1}\)；联合训练模型的测试FCI能量RMSE仍约0.06864 Hartree。必须将同一检查点的这些记录相互对应，不能拿不同模型拼接论证。

另一个源代码随机势虽有数值一致性，路径候选点却非驻点，且能量单位未经物理校准。这是独立的源码审计案例，不能与已训练H₂模型的表现混为一谈。等变性与自动微分是已有基础条件，本稿不把这一原则宣称为首次发现。[Unke等，R01](https://doi.org/10.1021/acs.chemrev.0c01111)、[Batzner等，R02](https://doi.org/10.1038/s41467-022-29939-5)。

**C3：收敛必须针对所需物理量，且束缚态计数可出现两个方向的错误。** cc-pVDZ参数化H₂ Morse模型具有17个全实线解析束缚态；8 Å主区间仅给出16个。最浅态解析结合能仅\(6.6963668574\times10^{-7}\) Hartree。96 Å、115200区间的已有单网格结果恢复了该态，结合能仍偏高1.56403%。在两个固定边界网格上做标准二阶Richardson后处理，得到−0.0268823%的偏差；这是事后估计，没有增加本征问题，也不能替换原单网格记录或宣称已独立验证96 Å处的收敛阶。

反方向的负结果同样重要：400区间粗网格下，D₂控制势、STO-3G模型、cc-pVDZ模型分别给出29/26/24个阈下态，而其全实线解析数目为28/25/23。这些是离散误差造成的警示，不是发现了新的物理态。矩阵残差只描述离散算子的求解；径向Dirichlet边界与全实线解析问题本来也不同。即使Morse核方程完全解准，它相对25个FCI电子点的模型RMSE仍分别为0.00757282和0.00420654 Hartree。数值求解器不能消除势能表示误差。[Morse，R11](https://doi.org/10.1103/PhysRev.34.57)。

**C4：提供多类验证目标不等价的可追溯案例。** 固定迭代次数的NEB完成状态并不等于力收敛；随机网络的守恒不等于该网络的机理已经成立；恒温平均值正确不等于涨落符合目标系综；水合自由能标签不等于氧化还原性质；解析峰形假设错误可以影响回收量。此处贡献是案例记录与目标—参照的对应关系，并不是对文献失败比例或自动检测性能的估计。[CI-NEB，R12](https://doi.org/10.1063/1.1329672)、[Gillespie，R14](https://doi.org/10.1021/j100540a008)、[Berendsen等，R15](https://doi.org/10.1063/1.448118)、[FreeSolv，R16](https://doi.org/10.1007/s10822-014-9747-x)。

建议以C1和C3承担主要量化贡献，C2连接诊断逻辑，C4放入支持信息或较短的拓展部分。

### 与直接相关工作的差距

| 既有研究 | 已经完成的工作 | 本稿可保留的区别及未完成部分 |
|---|---|---|
| [Fu等，2023，R03](https://openreview.net/forum?id=A8pqQipwkt) | 超越局部预测误差，用模拟结果评价机器学习力场。 | 本稿增加配对RHF/FCI归因和阈值态数值对照，但化学广度、真实反应动力学远不及该类基准。 |
| [Kamath等，2018，R07](https://doi.org/10.1063/1.5003074) | 对相同甲醛样本训练NN/GP，并直接计算、比较振动谱。 | 本稿的域长/网格反例可复核，但尚未完成神经网络势到核谱的连接，不能宣称首次研究势面误差影响谱。 |
| [Scalia等，2020，R04](https://doi.org/10.1021/acs.jcim.9b00975) | 系统比较可扩展的分子预测不确定度方法。 | 本稿3种子和8点覆盖率只是一项负诊断，不构成新的校准算法或充分统计验证。 |
| [DataSAIL，2025及2026增补，R06/R25](https://doi.org/10.1038/s41467-025-58606-8) | 定义相似性约束拆分问题并开展多种实证比较。 | 本稿的哈希、访问记录有助审计，却不是新拆分算法。[增补文献](https://doi.org/10.1038/s41467-025-67495-w)中的PINDER反例也说明不能宣称单一拆分方法普遍最优。 |
| [dos Santos等，2025，R08](https://doi.org/10.1021/acs.jcim.4c01847) | 已用Morse数据增广改善反应机器学习势的解离曲线与键解离能。 | 本稿揭示参考与Morse表示误差，没有提出新的增广方法，也未验证更好的化学解离预测。 |
| [MoleculeNet，2018，R17](https://doi.org/10.1039/C7SC02664A) | 提供广泛任务、数据集、基线与不同拆分方式。 | 256分子FreeSolv子集可作基线负例，不能据此给出化学空间整体的模型排名。 |
| [Zhang等，2026，R18](https://doi.org/10.1038/s44286-026-00392-1) | 已发表知识图谱/代理辅助化学过程数字孪生构造与校准工作。 | 本仓库模拟过程与接口未经实验校准，不能将其写成已验证工业数字孪生。 |

本次是有明确范围的文献对照，并非穷尽检索，不支持“首次”“普适”“最先进”或“已达顶刊”的表达。

### 唐海涛相关实验文献及分区核验边界

已依据出版社作者列表和广西师范大学单位核实4篇真实有机电合成论文：R19吡啶羟甲基化、R20铁单原子介导有机阳极氧化、R21锰单原子催化硅烷氧化、R22含碲噁唑烷酮多组分电合成。上方英文部分和[文献登记表](references_verified.json)给出完整题名、DOI、作者与单位核验记录。R21在线年份为2023、卷年为2024；R22中唐海涛为共同作者，未标为通讯作者，不能自行升级其作者角色。

这些论文为课题提供实验背景，不等于本仓库复现了其催化剂、反应路径、产率或生物活性。引用也不证明宋思毅与唐海涛存在指导、合作或共同署名关系。用户提供的作者为宋思毅，单位为广西师范大学，其余关系不推断。

期刊分区和论文真实性是两项独立核验。本轮没有获得直接权威、年份明确的中科院大类分区或JCR四分位记录，因此有关条目保留`classification="unverified"`。这不是说相关期刊不属一区，而是当前文件不认证该项属性。不能把声誉、影响因子、Scopus Q1，或未写清体系/年份的目录改称中科院一区或某年JCR Q1。Psi4、Morse、NEB等原始方法，以及已出版的TMLR评价和ICML消息传递研究，因方法相关性纳入，不冒充“全部SCI一区”。

### 可执行的补强路线

以下均为**建议的后续工作，尚未执行**；不能计入当前算量。

1. **让验证程序本身成为可独立评价的方法。** 预先冻结目标物理量、单位、参考哈密顿量/基组、适用区间、数据访问、收敛条件与解释边界。由独立人员准备正确控制和失败案例，对比“仅运行成功”“对称性/导数检查”和完整目标物理量审计，分类型报告检出、漏检与误报。至少由外部人员在干净环境复现一条完整链。
2. **建立真正跨体系的参照。** 预先选择可负担、电子难点不同的双原子、少原子振动和反应路径体系；留出分子身份或完整轨迹，而非仅相邻几何。相同标签和计算预算下比较插值/RBF等简单模型及至少两类学习表示。增加随机种子不能替代增加化学体系。
3. **补上势能到振动谱的缺失连接。** 在足够分辨的相关电子势能曲线上检验基组序列和独立插值，再将参考PES、Morse模型及训练势分别输入相同、已收敛的核求解器。区分电子参考、插值/表示、径向边界和网格误差。先验证J=0束缚态，再按具体目标补旋转、连续态或非绝热项；对照实验前必须匹配哈密顿量和修正层级。
4. **针对用途校准不确定度。** 使用独立校准集与从未参与选择的身份/区间外测试集，报告覆盖率、区间宽度、误差相关性和适当评分。训练标签的拟合不确定度与标签本身的参考偏差必须分开。
5. **机理主张必须有化学专属证据。** 指定完整端点、荷电/自旋状态，完成驻点、一级鞍点及路径连通检验；电合成进一步明确电极、溶剂和参比电位，并与适当实验物理量比较。唐组文献可用于选择真实任务，却不能被源代码手工碎片替代。
6. **保持事后分析与前瞻验证的区别。** 当前改进传递率和Richardson均为冻结数据上的事后分析；下一轮外部验证前需冻结阈值和主要终点。保留负结果以及缓存评估、预派发失败等算量区别。软件测试数只能支撑实现质量，不能兑换化学验证数量。

现阶段可如实提交为有范围限定的方法案例研究。要增强为有广泛价值的方法论文，应优先增加独立检验能力、外部复现与跨化学体系证据，而不是增加图表或堆叠计算次数；不对发表概率作任何估计。
