# Audited metrics and simulation-only sequential control / 指标与模拟控制核验

This component is a reproducible numerical demonstration. It has no experimental observations, verified product structure, equipment connection or industrial qualification. The attachment's source outputs remain available separately. All reviewed metric outputs retain full numeric precision, including impossible faradaic efficiency, rather than clipping inconvenient results.

本模块属于已执行的数值演示，不含实测数据、已确认的产物结构、真实仪器连接或工业资质认定。原附件结果另行保留。计算值不会因超出预期而被截断；不存在把大于 100% 的法拉第效率修饰成 100% 的处理。

## 1. Metric boundary / 指标边界

For feed concentration in mol/L and flow in mL/min, product formation is

\[
\dot n_P=C_{in}Q\,10^{-3}(X/100)(S/100),\quad
\mathrm{FE}=\frac{nF\dot n_P}{60I}\,100\%,\quad
\mathrm{STY}=\frac{1.44\dot n_P M_P}{V_{mL}10^{-6}}.
\]

Here STY is kg/(m³ day), molecular weight is g/mol, and the product/feed stoichiometry is assumed 1:1. Electrical SEC is \(UI/(60\dot n_PM_P)\) kWh/kg. This electrical figure excludes pumping, cooling/heating, solvent recovery and downstream processing. Reactor product and recovered product are reported separately when an isolation recovery is supplied. Undefined ratios are JSON `null` with an explanatory flag; arbitrary 999 or 9999 values are not measurements.

PMI requires a declared material inventory relative to a product output. The reviewed calculator explicitly accepts incremental electrolyte, workup-solvent and water flows. When the boundary is incomplete it returns `declared_boundary_mass_intensity`, while `full_process_PMI` remains `null`. The algebraic value PMI − 1 is labelled as a declared-boundary waste ratio assuming all nonproduct mass is waste; recovered coproducts and recycle allocation require additional accounting. These distinctions follow the scope of the [ACS GCI Pharmaceutical Roundtable PMI tools](https://acsgcipr.org/tools/process-mass-intensity/), rather than claiming ACS conformance from a three-term feed estimate.

PMI 不能只用底物、偶联组分和反应溶剂的部分进料就宣称已覆盖全流程。本计算允许加入电解质、后处理溶剂和洗涤水；不完整边界只输出部分物料强度。原子经济性同样需要经核验的配平反应。附件中的 223.27、110.18 和 331.43 g/mol 仅保留为假定分子量，1:1 质量比 99.394212% 不能代替结构与反应身份确认。

The [arithmetic sensitivity record](results/control/metric_boundary_sensitivity.json) deliberately uses a hypothetical 20% conversion. The single-reaction charge-consistent case gives FE = 100%; halving the current while retaining two-electron reporting gives FE = **200%**, which is preserved and flagged. Adding assumed inputs of 0.005 g/min electrolyte, 0.5 g/min workup solvent and 0.2 g/min water, with 80% recovery, raises the declared-boundary mass intensity from **242.184624 to 893.604670**. These are sensitivity assumptions, not an optimized or measured process.

审计基准另见 [baseline_metrics.json](results/control/baseline_metrics.json)：100×48 网格、450 µL/min、0.48 V 过电位，计算转化率 **58.3872858750%**、电流 **42.2513750193 mA**。采用附件电压式 \(U=1.85+\eta+18.5I\) 和假定产物分子量，得到电压 **3.11165043786 V**、STY **29026.9472363 kg/(m³ day)**、电耗 **0.503254627151 kWh/kg**、部分物料强度 **82.9580003847**。该模型只含一个反应，FE ≈ 100% 是电荷守恒的数学结果，不是经实验确认的高选择性。

## 2. Genuine sequential recommendation / 真正的序贯推荐

The attachment's optimizer evaluates a fixed 4×3 grid; changing its `budget_iterations` argument does not allocate further trials. No surrogate is fitted and no EHVI acquisition is evaluated. The reviewed implementation actually fits posterior GPs after each observed point and selects a previously unobserved candidate using Monte Carlo expected hypervolume improvement. The foundational concept is described by [Daulton, Balandat and Bakshy (2020)](https://arxiv.org/abs/2006.05078); this small implementation is a **single-candidate finite-pool** version, not their parallel autodifferentiable qEHVI implementation. The [BoTorch acquisition documentation](https://botorch.readthedocs.io/en/stable/acquisition.html) provides the broader acquisition taxonomy.

Candidate inputs form a declared 9×9 grid: flow 100–1500 µL/min and overpotential 0.20–0.75 V. Each candidate uses the reviewed conservative transport solver at 80×32 cells. The objectives, both maximized, are **STY/50000** and **−electrical SEC**; the fixed reference is (0, −2). The single-reaction model makes FE constant and makes STY versus FE a degenerate optimization problem. We therefore do not invent side-reaction selectivity to produce a visually attractive front. The STY/SEC trade-off instead depends on the modelled current, flow and assumed voltage law; it does not establish independent physical objectives or a validated energy model.

每个种子先共享 5 个随机初始点，再进行 10 次独立序贯推荐；每种方法总预算为 15 点。8 个种子与 2 种方法产生 **16 个运行、240 条评价使用记录**。离线全池 **81 次唯一 PDE 评价**用来构造参照与缓存，另有 **1 次基准网格计算**；这些额外参考计算不冒充每种方法 15 次预算的一部分。GP 只能使用已选择点的目标值，未知点的目标不参与归一化或拟合。

Each objective has an independent normalized Matérn-5/2 GP with fixed length scales (0.3, 0.3), unit amplitude and a 10⁻⁸ numerical nugget. Inputs are normalized using declared bounds, while output centering/scaling uses observations only. No hyperparameter optimization or experimental noise fitting is performed. The 256 posterior draws per candidate use common random numbers across candidates. For each draw, exact 2D rectangle-union HVI is computed. Samples with negative STY or nonpositive SEC contribute zero acquisition; this is an assumed physical-domain mask, not a learned operational-safety constraint. Independent GP posteriors also omit objective correlations. All remaining candidates are scored; no candidate is evaluated twice within a run.

The full pool has **46 nondominated STY/SEC points**, with scaled hypervolume **1.52157817325**. Raw FE ranges from **99.9999999999689 to 100.0000000000368%**, reflecting floating-point residuals. Using the exact FE identity, only the maximum-STY candidate is nondominated in the original STY/FE formulation.

| Seed | MC-EHVI fraction of full-pool HV | Random fraction of full-pool HV | Paired HV difference |
|---:|---:|---:|---:|
| 0 | 0.985029 | 0.981589 | 0.005235 |
| 1 | 0.981783 | 0.947877 | 0.051591 |
| 2 | 0.984477 | 0.986432 | −0.002974 |
| 3 | 0.983069 | 0.984923 | −0.002821 |
| 4 | 0.975342 | 0.948635 | 0.040636 |
| 5 | 0.977870 | 0.897247 | 0.122673 |
| 6 | 0.971559 | 0.986345 | −0.022498 |
| 7 | 0.977219 | 0.943429 | 0.051414 |

Mean fractions are **0.979543496** for MC-EHVI and **0.959559623** for random selection. MC-EHVI wins **5/8**, and random selection wins **3/8**. The mean paired HV difference is **0.0304070243**, with sample SD **0.0465353274**; the 95% paired-seed t interval is **[−0.0084974830, 0.0693115316]**, which includes zero. This interval assumes independent, approximately normal seed effects and only describes this finite-pool numerical benchmark. It is not chemical uncertainty, an experimentally calibrated confidence interval, or evidence of generally superior optimization.

平均覆盖率有所改善，但种子 2、3、6 的随机对照表现更好。区间包含零，不能据此声称统计显著、普遍优越或获得工业最优条件。所有非支配候选保存在 [campaign_pareto_records.csv](results/control/campaign_pareto_records.csv)，不会用前沿的中位数自动宣称“最佳工业折中”。

## 3. Hardware simulation / 仪器状态模拟

`SimulationOnlySCPI` accepts only `SIMULATOR::ELECTRATWIN`; USB, GPIB and TCPIP resources are rejected. It imports no socket, PyVISA or instrument-vendor driver. Its identity is `ELECTRATWIN,SIMULATOR,NO-PHYSICAL-DEVICE,1.0`. The voltage response \(1.95+22.4I\) is an exact toy formula, not acquired telemetry and not a model of a specific manufacturer's instrument. The 0–0.1 A and 0.01–10 V limits are simulator bounds, not certified equipment ratings.

Unknown commands, malformed or nonfinite values and out-of-range settings fail explicitly. Output cannot start before connection. A simulated voltage-compliance violation disables output and latches a fault. Reset requires output off and zero current. Reading with output disabled reports no measurement; disconnect clears current and disables output. [The saved trace](results/control/scpi_simulator_trace.json) ends with zero current, output off and no connection. This is software state-machine verification only; real electrical interlocks and a specific instrument's command syntax require independent engineering work.

模拟器没有连接任何真实设备，也没有使用固定商业型号序列号伪装设备识别。故障锁存、关断及输入校验测试只能证明该本地状态机的行为，不能作为真实电化学平台安全性或驱动兼容性的认证。

## 4. Reproduction and tests / 复现与测试

Use the existing scientific Python environment with NumPy, SciPy and scikit-learn. Keep numerical-library threading bounded on a CPU-only host.

```text
python electratwin/scripts/metrics_control.py
python -m unittest discover -s tests -p test_electratwin_control.py -v
```

The default execution regenerates metric sensitivity, one baseline solve, the simulator trace, 81 candidate-pool solves, all 16 sequential runs, paired comparisons and summary provenance. `--skip-optimization` runs only arithmetic sensitivity and the simulator. `--budget` actually changes the per-run evaluation count, subject to initial-size and pool-size checks. Candidate-pool results are deterministic; random seeds are 0–7 by default.

**18 tests passed** on the installed environment. They independently check known molar/mass/energy units; preserved FE >100%; incomplete inventories; zero denominators; invalid inputs; rectangle-union geometry, translation/scaling and vectorized HVI versus independent recomputation; Pareto ties/directions; actual oracle-call budgets and common initial sets; and simulator resource, command, parameter, fault/compliance, measurement and disconnect states. Passing these software checks establishes neither transport-model validity nor experimental selectivity.

Inspect the machine-readable [optimization summary](results/control/optimization_summary.json), [all campaign evaluations](results/control/sequential_campaigns.csv), [paired comparisons](results/control/paired_budget_comparison.csv) and [candidate pool](results/control/candidate_pool.csv). Source hashes and runtime versions are saved in the optimization summary. Transport discretization and candidate-domain grid sensitivity are documented separately by the transport component; small algebraic residuals are not a claim of physical validation or complete grid independence.
