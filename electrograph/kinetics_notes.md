# Reviewed stochastic kinetics / 经审查的随机动力学

This module verifies an **assumed five-state independent-site model**. Its rates are inherited from the supplied demonstration; they are neither measured nor fitted. The chemical state names do not establish a POP-SAC catalytic mechanism, a spatial lattice effect, or industrial performance.

本模块对**假定的五状态独立位点模型**进行数值核验。速率沿用提供的演示代码，未经过实验测量或拟合。状态名称不能作为 POP-SAC 催化机理、空间晶格效应或工业性能的证据。

## 1. Model, corrections, and equations / 模型、修正与公式

The states, in order, are empty → substrate → radical → intermediate → adsorbed product → empty. Adsorption, SET, coupling, PCET, and desorption each advance a site by one state. SET and PCET each increment the electron counter once; desorption increments the product counter once. All five rates are irreversible pseudo-first-order rates in s⁻¹. The supplied constants are retained exactly: (F=96485.33), (R=8.314), (T=298.15), (alpha=0.5), and

\[
(k_0,k_1,k_2,k_3,k_4)=(45,120\exp[\alpha F\eta/(RT)],350,200\exp[\alpha F\eta/(RT)],80).
\]

For a rate class with state counts (n_i), event propensities are (a_i=n_i k_i). The direct stochastic simulation draws an exponential waiting time with rate (a_0^{\mathrm{tot}}=\sum_i a_i), then selects a channel with probability (a_i/a_0^{\mathrm{tot}}). Aggregating counts is exact for the independent, equivalent sites within each fixed rate class. There is no nearest-neighbor interaction, diffusion, hopping, lattice geometry, or periodic boundary operation. The algorithm follows the direct-SSA construction; “exact” means exact sampling of this specified stochastic model, not exact chemistry. [Gillespie (1977)](https://pubs.acs.org/doi/10.1021/j100540a008).

五个状态依次为空位、吸附底物、自由基、中间体和吸附产物；循环五步分别对应吸附、SET、耦合、PCET 和脱附。SET 与 PCET 各计入一个电子，脱附计入一个产物。所有速率均为不可逆拟一级速率。对于同一固定速率类别，事件倾向为 (a_i=n_i k_i)。计数型 Gillespie 抽样与这些独立等价位点的逐位点过程具有相同的总体计数分布；它不包含近邻、扩散、跳跃、空间几何或周期边界操作。

The supplied code logged pre-event populations with post-event timestamps and returned pre-event final populations. Its reported step count was the requested cap, even if execution terminated early. The reviewed implementation records post-event populations, the actual executed event count, and a final censoring row. The first scheduled event after the prescribed horizon is not executed. Occupancy is the piecewise-constant time integral over the observation window, not an unweighted average of event snapshots:

\[
\bar\theta_i={1\over N(T-b)}\int_b^T n_i(t)\,dt,\qquad
\widehat{\mathrm{TOF}}={N_{\mathrm{des}}(b,T]\over N(T-b)}.
\]

原代码将事件前占据数与事件后时间戳配对，返回的最终占据数也落后一个事件；提前终止时仍报告预设事件数。修订版记录事件后的状态、实际事件数和固定终点。超出终点的下一个事件被截尾而不执行。时间平均占据数通过分段常数轨迹积分得到，不使用事件快照的普通算术平均。

Let (Q_{ii}=-k_i), (Q_{i,i+1}=k_i), with cyclic indexing. The independent-site continuous-time Markov chain (CTMC) satisfies (p(t)=p(0)\exp(Qt)). A 10×10 augmented matrix exponential gives both probabilities and their time integrals, avoiding an inverse of singular (Q). The numerical implementation uses SciPy's matrix exponential. [SciPy `expm` documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.expm.html).

\[
\pi_i={1/k_i\over \sum_j1/k_j},\quad
\mathrm{TOF}_{\infty}={1\over\sum_j1/k_j},\quad
{\partial\log\mathrm{TOF}_{\infty}\over\partial\log k_i}=\pi_i.
\]

For a complete cycle starting from an empty site, the waiting time is the sum of five independent exponential waits. Its mean is (sum_i1/k_i), its variance is (sum_i1/k_i^2), and its CDF is computed from an acyclic five-state phase-type generator. These are conditional model identities, not experimental lifetime distributions.

CTMC 通过矩阵指数计算有限时间占据概率与积分。稳态占据数正比于各步平均停留时间 (1/k_i)，稳态 TOF 为五步平均停留时间之和的倒数。对速率的对数弹性恰好等于相应稳态占据分数。一个完整循环的等待时间是五个独立指数等待时间之和，保存的等待时间 CDF 同样只属于假定模型。

## 2. Executed calculations / 已执行计算

The full study contains **273 SSA trajectories and 9,950,504 executed events**. It includes 192 potential-scan trajectories (six potentials × 32 seeds), 48 site-scaling trajectories (64/256/1024 sites × 16 seeds), 16 two-class trajectories, one baseline with every event recorded, and 16 fixed-event diagnostics. The first 257 trajectories use \(T=1\) s and burn-in \(b=0.2\) s; the diagnostics stop at 3,500 actual events and use no burn-in. All start from empty sites. The default potential scan has 256 sites. Also saved are 105 labeled analytic rate-perturbation cases, 101 transient CTMC times, and 151 cycle-wait CDF times. The separate pilot contains five trajectories and 169,216 events; pilot and unit-test calculations are excluded from the full-study counts.

全量研究共 **273 条 SSA 轨迹和 9,950,504 个实际事件**：电势扫描 192 条、位点数扫描 48 条、两类速率位点 16 条、逐事件基准 1 条以及固定事件诊断 16 条。前 257 条采用固定终点 1 s 和预热 0.2 s，最后 16 条在实际完成 3,500 个事件后终止且不设预热。另保存 105 个带标签的解析速率扰动条件、101 个 CTMC 时间点和 151 个循环等待 CDF 时间点。独立 pilot 的 5 条轨迹和 169,216 个事件，以及单元测试计算，均不计入全量研究计数。

| η / V | Mean SSA TOF / s⁻¹ | 95% Monte Carlo interval for mean / s⁻¹ | Finite-window CTMC TOF / s⁻¹ |
|---:|---:|---:|---:|
| 0.2 | 26.413727 | 26.312386–26.515068 | 26.419159 |
| 0.3 | 26.627655 | 26.549261–26.706049 | 26.582874 |
| 0.4 | 26.660767 | 26.578592–26.742941 | 26.606421 |
| 0.5 | 26.619110 | 26.506958–26.731262 | 26.609788 |
| 0.6 | 26.592865 | 26.493277–26.692453 | 26.610268 |
| 0.7 | 26.564178 | 26.462074–26.666283 | 26.610337 |

Intervals use the Student-t distribution and the sample standard deviation across independent PRNG seeds, conditional on the fixed model. All six intervals contain their corresponding finite-window CTMC means. This is a consistency check for these draws, not proof of nominal coverage, fitted-rate uncertainty, or experimental validity. Multiplicity is not adjusted. The maximum relative mean discrepancy over this scan is 0.2043%. Individual estimates and finite-sample means can exceed an upper bound on the model's expected steady TOF.

区间使用各独立随机种子之间的样本标准差和 Student-t 分布，仅描述固定参数下均值的 Monte Carlo 抽样误差。六个区间均包含各自有限时间 CTMC 均值；本次覆盖情况不证明严格覆盖率，也不包含速率不确定性或实验误差。未作多重比较校正。扫描均值相对有限时间 CTMC 的最大偏差为 0.2043%。随机样本的 TOF 可以高于模型稳态期望的上界，二者不能混为一谈。

The original plotted list `[4.2, 12.8, 38.4, 95.1, 184.6, 260.2]` was hard-coded and is not used here. With the inherited rates, even infinitely fast electron-transfer steps leave the adsorption/coupling/desorption bound

\[
\mathrm{TOF}_{\infty}<\left({1\over45}+{1\over350}+{1\over80}\right)^{-1}=26.61034847\ \mathrm{s}^{-1}.
\]

At η=0.48 V, the analytic steady TOF is 26.60952060 s⁻¹; adsorption and desorption log elasticities are 0.59132268 and 0.33261901. SET and PCET elasticities are only (1.9444\times10^{-5}) and (1.1667\times10^{-5}). Consequently, this particular model predicts saturation, not the source's steep TOF rise. The mean cycle time is 0.03758053424 s and its variance is 0.00065824043 s².

原图中的 TOF 列表为硬编码，本次不予采用。沿用该速率时，吸附、耦合与脱附已经将稳态 TOF 期望限制在 26.61034847 s⁻¹ 以下。η=0.48 V 时，吸附和脱附的速率对数弹性分别约为 0.5913 和 0.3326，而 SET、PCET 约为 (10^{-5}) 数量级。因此该特定循环预计趋于饱和，不能支持原图持续大幅增长的趋势。

| Sites / 位点数 | Seeds / 种子数 | Sample SD of TOF / TOF样本标准差 | Long-time renewal SD approximation / 长时间更新过程近似 |
|---:|---:|---:|---:|
| 64 | 16 | 0.425537 | 0.492167 |
| 256 | 16 | 0.227059 | 0.246084 |
| 1024 | 16 | 0.115410 | 0.123042 |

The renewal approximation scales as (N^{-1/2}) for fixed observation time. Its values are asymptotic, not an exact finite-window variance benchmark. Sixteen seeds give limited precision for an estimated standard deviation. A separate static heterogeneous case assigns 192 sites reference rates and 64 sites one-quarter adsorption and one-half desorption rates. Its mean TOF is 22.15179443 s⁻¹, with a mean interval of [21.99441633, 22.30917254], compared with a finite-window weighted CTMC mean of 22.09851785 s⁻¹. There are still no spatial interactions. The fixed-3,500-event diagnostics terminate at 0.098456–0.103911 s and have mean apparent TOF 24.20722140 s⁻¹; they cannot establish steady state or be treated as unbiased fixed-time observations.

位点数研究中标准差随位点数增加而下降，与 (N^{-1/2}) 抽样缩放一致，但长时间更新过程公式只是渐近近似，16 个种子也不足以高精度估计标准差。两类速率情形将 192 个位点设为参考速率，64 个位点的吸附速率降至四分之一、脱附速率降至二分之一，得到均值 22.15179443 s⁻¹；对应有限时间解析加权均值为 22.09851785 s⁻¹。这是实际赋予不同速率的静态类别，仍没有空间相互作用。3,500 事件诊断仅运行约 0.1 s，平均表观 TOF 为 24.20722140 s⁻¹；这类路径依赖终止规则不能被描述成稳态采样。

## 3. Conservation and interpretation / 守恒与解释

For every trajectory, state populations obey (n_i(T)-n_i(b)=C_{i-1}-C_i) exactly in integer arithmetic. With state electron inventories (q=(0,0,1,1,2)), the identity is (C_\mathrm{SET}+C_\mathrm{PCET}-2C_\mathrm{des}=q\cdot[n(T)-n(b)]). All 273 trajectories satisfy these identities; the largest time-occupancy sum error is (9.33\times10^{-15}). The maximum CTMC probability-sum discrepancy in the six-potential scan is (6.36\times10^{-10}); the stiff matrix exponential is not claimed to be exact real arithmetic. The baseline's 5,363 observed products and 10,745 electrons differ from a strict 2:1 ratio by the 19-electron inventory change on occupied sites. Dividing electrons by observation time and multiplying by the elementary charge gives current for the explicitly modeled sites only; site density, electrode area, and transport are absent, so this is not a macroscopic electrode-current prediction.

所有 273 条轨迹均满足整数形式的位点守恒与电子—产物—表面库存恒等式。时间占据数和的最大误差为 (9.33\times10^{-15})。六个电势下 CTMC 概率和的最大数值偏差为 (6.36\times10^{-10})。基准观测到 5,363 个产物、10,745 个电子；与严格 2:1 比例相差的 19 个电子来自表面库存变化。电流只是这些显式位点的计数换算，没有位点密度、电极面积和传质参数，因此不能作为宏观电极电流。

## 4. Files and reproduction / 文件与复现

Run from the repository root using the existing chemistry Python environment:

```powershell
python electrograph/scripts/kinetics_reviewed.py --pilot
python electrograph/scripts/kinetics_reviewed.py
python -m unittest discover -s tests -p test_electrograph_kinetics.py -v
```

The public functions are `rate_constants`, `generator_matrix`, `stationary_solution`, `ctmc_window`, `cycle_wait_cdf`, and `simulate`. `simulate` accepts explicit seeds, rate classes, initial counts, fixed horizon, burn-in, event logging stride, and an optional diagnostic event cap. A cap reached before the observation window raises an error. Main artifacts are [summary](results/kinetics/summary.json), [all trajectory records](results/kinetics/trajectories.csv), [potential summary](results/kinetics/potential_summary.csv), [complete baseline event log](results/kinetics/baseline_trajectory.csv), [analytic rate perturbations](results/kinetics/analytic_rate_perturbations.csv), [site scaling](results/kinetics/site_scaling.csv), and [15-test log](results/kinetics/test_results.txt). Input-source and output hashes are recorded in the summary. The numerical checks verify this implementation; they do not substitute for fitted rates, a reversible thermodynamic network, neighbor interactions, transport coupling, or experimental calibration.

公开接口支持显式随机种子、速率类别、初始计数、固定终点、预热、逐事件记录和诊断用事件上限。超过事件上限而尚未进入观测窗口时会报错。结果包含源文件与 CSV 哈希。15 项单元测试覆盖生成矩阵、解析极限、守恒、终点截尾、时间加权积分、异类位点、固定事件终止与可复现性；其通过不代表真实材料或化学机理已经验证。
