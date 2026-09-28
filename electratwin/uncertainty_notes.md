# Parameter-scenario propagation and Sobol diagnostics / 参数情景传播与 Sobol 诊断

This extension executes the frozen conservative transport model under **selected independent parameter distributions**. These are not measurement error distributions, calibrated posteriors, wet-lab data or demonstrated manufacturing variability. The labels A and P retain the earlier abstract 1:1 molecular bookkeeping; molar closure does not establish elemental or full-component mass closure of a real coupling reaction. Product mass and voltage remain assumptions.

本扩展对冻结的守恒输运模型执行**所选独立参数分布**下的情景传播。这些分布不是测量误差、标定后验、湿实验数据或已验证的制造波动。A/P 仍是此前一比一分子记账的抽象标签；摩尔闭合不代表真实偶联反应的元素或全组分质量闭合。产物分子量与电压关系仍是模型假设。

## 1. Declared distributions / 声明的输入分布

| Parameter / 参数 | Distribution / 分布 | Lower / 下界 | Upper / 上界 | Units / 单位 |
|---|---|---:|---:|---|
| Flow Q / 流量 | Uniform / 均匀 | 405 | 495 | µL/min |
| Potential index η / 电位参数 | Uniform / 均匀 | 0.46 | 0.50 | V |
| Diffusivity D / 扩散系数 | Uniform / 均匀 | 0.77e-9 | 1.43e-9 | m²/s |
| Effective current scale / 有效电流尺度 | Log-uniform / 对数均匀 | 0.04 | 0.16 | A/m² |
| Channel height H / 通道高度 | Uniform / 均匀 | 0.00027 | 0.00033 | m |

The median scenario is Q = 450 µL/min, η = 0.48 V, D = 1.1 × 10⁻⁹ m²/s, current scale = 0.08 A/m² and H = 0.0003 m. For a unit coordinate u, uniform parameters use x = a + u(b−a), whereas the current scale uses x = a exp[u log(b/a)]. Thus its geometric midpoint is 0.08, not the arithmetic midpoint 0.10. Independence is imposed for this sensitivity design, not inferred from experiments; in reality, kinetic and transport parameters may covary. Changing their ranges or dependence can change every sensitivity ranking.

中位情景为 Q = 450 µL/min、η = 0.48 V、D = 1.1 × 10⁻⁹ m²/s、电流尺度 0.08 A/m²、H = 0.0003 m。均匀分布采用 x = a + u(b−a)，电流尺度采用 x = a exp[u log(b/a)]，因此其中位数为几何中点 0.08，而不是算术中点 0.10。参数独立性是分析设定，未经实验建立；改变取值范围或相关性会改变敏感性排序。

Fixed quantities are temperature 298.15 K, length 0.06 m, width 0.012 m, inlet A/P = 50/0 mM and n = 2. The assumed product molecular weight is 331.43 g/mol and the voltage law is U = 1.85 + η + 18.5I. All main uncertainty calculations use 60 × 24 cells with concentration fields omitted from returned records. The outputs are conversion, current, assumed-product STY and reactor-product electrical SEC. STY uses the reactor volume of each sampled height, rather than a fixed baseline volume. SEC excludes pumping, thermal control and downstream processing. No extra FE variation is invented: FE is already fixed by the single-reaction charge identity.

固定温度为 298.15 K，长度为 0.06 m，宽度为 0.012 m，入口 A/P 为 50/0 mM，n = 2。产物分子量假定为 331.43 g/mol，电压为 U = 1.85 + η + 18.5I。主分析采用 60 × 24 网格，不返回浓度全场。输出为转化率、电流、假定产物 STY 和按反应器产物计的电耗。每个情景的 STY 使用对应通道高度下的体积，未固定为基准体积；电耗不含泵、温控及后处理。单反应 FE 的恒等性质不被人为添加波动。

## 2. Sampling, estimators and exact accounting / 采样、估计与计数

The [official SciPy Sobol documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.qmc.Sobol.html) was opened on 28 September 2026. It describes LMS plus digital-shift scrambling and recommends power-of-two sample sizes to retain sequence balance. The implementation uses `qmc.Sobol(d=10, scramble=True, rng=default_rng(20260928))` and `random_base2(m=8)`, without skipping, thinning or optimization. A is the first five coordinates; B is the last five; ABᵢ copies A and replaces only column i by its counterpart in B. This gives 256 × (2 + 5) = **1792 main PDE solves**.

已于 2026 年 9 月 28 日打开 [SciPy Sobol 官方文档](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.qmc.Sobol.html)，核验了 LMS 加数字平移扰乱、二次幂采样及跳点或稀释会破坏平衡性质的说明。程序采用十维扰乱 Sobol、种子 20260928、`random_base2(m=8)`，不跳点、不稀释、不做后处理优化。前五维为 A，后五维为 B；ABᵢ 只将 A 的第 i 列替换为 B 对应列。主计算共 256 × (2 + 5) = **1792 次 PDE 求解**。

For each output, V is the `ddof=1` sample variance of concatenated YA and YB, and the Jansen estimates are

\[
S_i=1-\frac{\operatorname{mean}[(Y_B-Y_{AB_i})^2]}{2V},\qquad
S_{T_i}=\frac{\operatorname{mean}[(Y_A-Y_{AB_i})^2]}{2V}.
\]

Here the squared differences are averaged: `mean((YB − YAB_i)**2)` and `mean((YA − YAB_i)**2)`. Negative estimates, estimates above one and finite-sample violations of Sᵢ ≤ STᵢ are retained. Constant or numerically unresolved output variance is rejected instead of assigned an arbitrary sensitivity. Prefixes N = 64, 128 and 256 reuse saved evaluations; no extra solves are attributed to them.

每个输出使用拼接 YA、YB 的样本方差 V（`ddof=1`）。公式中的均值是**差值平方的均值**，不是均值再平方。程序不截断负估计、大于一的估计或有限样本下 Sᵢ > STᵢ 的结果；常量或数值上无法分辨的方差会明确拒绝。N = 64、128、256 前缀重用已有计算，不增加 PDE 次数。

| Analysis item / 分析项目 | PDE solves / PDE 次数 | Other records / 其他记录 |
|---|---:|---:|
| A, B and five hybrids / A、B 与五组混合矩阵 | 1792 | 1792 raw rows |
| Eight preset points on two meshes / 八个预选点双网格 | 16 | 32 output comparisons |
| Local center and two symmetric step sizes / 局部中心及两步长 | 21 | 40 derivatives |
| Cached prefix estimates / 缓存前缀估计 | 0 | 60 index rows |
| Paired-row bootstrap / 配对行自助重采样 | 0 | 500 resamples; 10000 index rows |
| Main analysis total / 主分析合计 | **1829** | No experiments / 无实验 |

All 1829 analysis solves passed numerical checks; the maximum relative material-balance error was **2.6901138973 × 10⁻¹⁴**, charge-balance error **2.6848851881 × 10⁻¹⁴** and linear residual **1.1657579284 × 10⁻¹⁴**. All four reported outputs were positive. The run took 21.48 s on this local CPU configuration, which is not a latency guarantee. Unit-test solver calls are separate software checks, not counted as part of this analysis budget.

1829 次分析求解全部通过数值检查；最大相对物料误差为 **2.6901138973 × 10⁻¹⁴**，电荷误差为 **2.6848851881 × 10⁻¹⁴**，线性残差为 **1.1657579284 × 10⁻¹⁴**，四个输出全部为正。本机 CPU 此次耗时 21.48 s，不构成时延保证。单元测试中的求解属于额外软件检查，不纳入主分析计数。

## 3. Propagated scenario distribution / 传播后的情景分布

Quantiles use only the 512 base A/B evaluations; hybrid rows are excluded to avoid counting pick-freeze reuse as additional base observations. The following central range describes the selected input scenarios, **not a 95% experimental or posterior confidence interval**. Combining the two five-dimensional projections does not claim that their union is a newly constructed 512-point balanced Sobol net.

分位数仅使用 A/B 的 512 个基础输出，不纳入混合矩阵行，避免把冻结替换的重用当作额外基础观测。以下中央范围描述所选输入情景，**不是 95% 实验区间或标定后验区间**。将两个五维投影合并，也不意味着它们构成了另一个经证明平衡的 512 点 Sobol 网。

| Output / 输出 | 2.5th percentile / 2.5% 分位 | Median / 中位数 | 97.5th percentile / 97.5% 分位 |
|---|---:|---:|---:|
| Conversion / 转化率, % | 46.65274059 | 58.21646944 | 68.29875939 |
| Current / 电流, mA | 34.04646110 | 41.76003273 | 49.17342420 |
| Assumed-product STY / 假定产物 STY, kg/(m³ day) | 22007.89154 | 28579.14634 | 36384.69734 |
| Electrical SEC / 电耗, kWh/kg | 0.4781755770 | 0.5012882243 | 0.5253083785 |

Conversion spans 44.84346541–72.41217210% across the 512 base points. That range does not cover omitted model discrepancy, unmeasured side reactions, uncertainty in molecular identity, temperature or the voltage law. The apparent narrow SEC spread partly reflects the single-reaction current/product proportionality. It should not be interpreted as a robust experimental energy-efficiency prediction.

512 个基础点的转化率最小至最大值为 44.84346541–72.41217210%。这个范围没有覆盖模型结构误差、未测副反应、分子身份、温度或电压模型的不确定性。SEC 的相对较窄分布部分来自单反应电流与产物速率成比例的结构，不能解释为真实电耗具有同样稳定性。

## 4. Sensitivity rankings, finite-sample limits and resampling / 敏感性排序、有限样本限制与重采样

The full-N total-order estimates below describe variance allocation **under the stated independent distributions**. D dominates conversion and current. D and H have similar contributions to STY because H also changes the volume in its denominator. This is a sensitivity result for the model and chosen ranges, not a universal ranking of experimental control priorities.

下表为完整 N 下的总效应估计，只适用于**所声明独立分布**。D 对转化率与电流的影响最大；在 STY 中，D 与 H 的贡献相近，因为 H 还改变了体积归一化分母。这不是对真实实验控制优先级的普遍结论。

| Parameter / 参数 | ST conversion / 转化率 | ST current / 电流 | ST STY | ST SEC |
|---|---:|---:|---:|---:|
| Q | 0.099932 | 0.067665 | 0.040956 | 0.063705 |
| η | 0.021165 | 0.021856 | 0.013584 | 0.080110 |
| D | 0.745056 | 0.767547 | 0.460676 | 0.722628 |
| Current scale / 电流尺度 | 0.053205 | 0.054670 | 0.032909 | 0.051470 |
| H | 0.078877 | 0.081341 | 0.448519 | 0.076581 |

For conversion, the first/total estimates for D are 0.740840/0.717498 at N = 64, 0.737995/0.768736 at N = 128 and 0.750875/0.745056 at N = 256. The first-order H estimate changes from 0.132477 to 0.079222 to 0.087377. These prefix changes are preserved instead of asserting complete convergence. At N = 256 the sums of first-order estimates are 1.030520, 1.040820, 1.050985 and 1.040439 for conversion, current, STY and SEC. Some first-order estimates exceed their corresponding total-order estimates. These violations of exact-index identities are finite-estimator error; they are not evidence of negative physical interactions. The code does not renormalize or clip them to make the decomposition appear exact.

以转化率为例，D 的一阶/总效应估计在 N = 64、128、256 时分别为 0.740840/0.717498、0.737995/0.768736、0.750875/0.745056；H 的一阶估计由 0.132477 变为 0.079222，再变为 0.087377。因此不能宣称估计已完全收敛。N = 256 时，转化率、电流、STY、SEC 的一阶和分别为 1.030520、1.040820、1.050985、1.040439，且部分一阶估计超过总效应估计。这些违反精确指标恒等关系的现象说明有限估计误差，不能解释成物理上的负交互作用；程序没有重新归一化或截断它们来制造精确分解的外观。

The 500 bootstrap resamples use seed 20260929 and resample matching row indices across A, B and every ABᵢ, recalculating variance each time. Their 2.5/97.5 percentiles are explicitly **paired-row resampling stability ranges**. Scrambled-QMC rows are not iid, and this resampling destroys net balance; the ranges have no rigorous QMC confidence-coverage guarantee. No independent scrambling replicate was run. For conversion D, the total-order resampling range is [0.643705, 0.873640]. For STY, D is [0.388113, 0.546007] and H is [0.381538, 0.528240], so these diagnostics do not separate their contributions decisively. Several first-order ranges cross zero; those negative values remain in the CSV files.

500 次自助重采样使用种子 20260929，对 A、B、全部 ABᵢ 同时使用相同的行索引，并每次重算方差。2.5/97.5 分位明确称为**配对行重采样稳定性范围**。扰乱 QMC 的行不是独立同分布样本，这种重采样会破坏网格平衡，不能保证严格 QMC 区间覆盖率；本次没有执行多个独立扰乱副本。转化率 D 的总效应稳定性范围为 [0.643705, 0.873640]。STY 中 D 为 [0.388113, 0.546007]、H 为 [0.381538, 0.528240]，明显重叠，不能据此精确区分二者贡献。多个一阶范围跨零，负值均原样保存在 CSV。

## 5. Mesh and local-slope diagnostics / 网格与局部斜率诊断

Eight unit-coordinate points were declared before results were examined, including two all-low/all-high corners, two alternating corners, the center and three mixed points. Each was evaluated at 60 × 24 and 180 × 72, yielding sixteen solves. The maximum absolute conversion difference is **0.1598292455 percentage points**. The largest relative coarse/refined difference is 0.2544885% for conversion, current and STY, and 0.0718705% for SEC. This is an eight-point comparison, not a domain-wide error bound; Sobol indices were not recomputed on the refined mesh. Grid error is not included in the scenario quantiles or bootstrap ranges.

八个单位坐标点在查看结果之前就已确定，包括全低和全高角点、两个交替角点、中心及三个混合点。分别采用 60 × 24、180 × 72 网格，共十六次求解。最大转化率差为 **0.1598292455 个百分点**；转化率、电流与 STY 的最大相对差为 0.2544885%，SEC 为 0.0718705%。这是八点比较，不是全域误差上界；没有在加密网格重算 Sobol 指标。网格误差未被纳入情景分位数或自助范围。

The local calculation uses the quantile center and symmetric unit-coordinate steps 0.02 and 0.01: one center plus 5 × 2 signs × 2 steps = twenty-one solves. Central derivatives dY/du are transformed to physical-coordinate derivatives using the transform Jacobian at the center; local elasticity is (x/Y)dY/dx. At step 0.01, conversion elasticities for [Q, η, D, current scale, H] are [−0.556606, 0.521692, 0.500760, 0.055849, −0.500765]. STY elasticities are [0.443395, 0.521692, 0.500760, 0.055849, −1.500772]. The largest elasticity difference between the two step sizes, over all outputs and inputs, is 0.000023071. This checks local finite-difference stability on the coarse model; it does not estimate global variance contributions or actual experimental gradients. In particular, the relatively large η elasticity does not contradict its small conversion Sobol effect, because its assigned variation range is narrow.

局部计算采用分位坐标中心和 0.02、0.01 两组对称步长，共一个中心加 5 × 2 个方向 × 2 个步长，即二十一次求解。先计算 dY/du，再用中心处变换雅可比换算为物理坐标导数；局部弹性定义为 (x/Y)dY/dx。0.01 步长下，[Q、η、D、电流尺度、H] 的转化率弹性为 [−0.556606、0.521692、0.500760、0.055849、−0.500765]，STY 弹性为 [0.443395、0.521692、0.500760、0.055849、−1.500772]。两步长在所有输入输出上的最大弹性差为 0.000023071。这只检验粗网格模型局部差分稳定性，不是全局方差贡献或实验梯度。η 的局部弹性较大，但指定变化范围很窄，因此其转化率 Sobol 效应较小并不矛盾。

## 6. Reproduction and tests / 复现与测试

Use the existing scientific environment, with OMP/BLAS/MKL threads set to one:

```text
python electratwin/scripts/uncertainty_analysis.py
python -m unittest discover -s tests -p test_electratwin_uncertainty.py -v
```

使用现有科学环境，并将 OMP/BLAS/MKL 线程设为一。默认分析包含所有上述计算；`--power` 与 `--bootstrap-replicates` 可改变规模并被完整记录。重跑前先保存交付归档，耗时、时间戳及相关文件哈希会变化。程序只写自己的不确定性结果目录，不覆盖冻结的输运、控制、源程序或此前报告。

**Thirteen new tests passed.** They check uniform/log-uniform transforms and bounds, invalid coordinates, reproducible Sobol prefixes, exact pick-freeze column replacement, analytic additive indices, the known first-order 3/7 and total-order 4/7 indices of a product of two independent U(0,1) variables, rejection of constant variance, retention of negative/above-one estimates, matched-row bootstrap, independent unit arithmetic, small-grid PDE conservation, the twenty-one-point local design and rejection of failed PDE results. These tests verify implementation properties and do not calibrate the scientific model.

**新增十三项测试全部通过。** 覆盖分布变换和边界、非法坐标、Sobol 前缀复现、冻结替换矩阵、解析加性指标、两个独立 U(0,1) 变量乘积的一阶 3/7 与总效应 4/7、常量方差拒绝、保留负或大于一的估计、配对行重采样、独立单位算术、小网格 PDE 守恒、二十一点局部设计与失败求解拒绝。这些测试验证实现性质，并未标定科学模型。

## 7. Result files and provenance / 结果文件与溯源

- [Summary and complete assumptions / 汇总与完整假设](results/uncertainty/summary.json): script and frozen-solver SHA-256, versions, seeds, counts, output hashes and limits.
- [All main design rows / 全部主设计行](results/uncertainty/design_evaluations.csv): block and paired-row identity, unit and physical inputs, outputs, residuals and evidence role.
- [N = 64/128/256 estimates / 前缀估计](results/uncertainty/sobol_prefix_estimates.csv) and [resampling ranges / 重采样范围](results/uncertainty/sobol_resampling_intervals.csv).
- [All 500 resamples / 全部重采样估计](results/uncertainty/sobol_bootstrap_replicates.csv) and [512-base-point quantiles / 基础点分位数](results/uncertainty/base_AB_quantiles.csv).
- [Grid evaluation rows / 网格评价](results/uncertainty/grid_check_evaluations.csv) and [grid differences / 网格差异](results/uncertainty/grid_check_differences.csv).
- [Local evaluation rows / 局部求解](results/uncertainty/local_difference_evaluations.csv) and [derivatives / 导数结果](results/uncertainty/local_derivatives.csv).
- [Execution log / 执行日志](results/uncertainty/run.log), [analysis script / 分析脚本](scripts/uncertainty_analysis.py) and [independent tests / 独立测试](../tests/test_electratwin_uncertainty.py).

Interpret these results as a reproducible analysis of how the declared model responds to the selected ranges. Model discrepancy, parameter identifiability, measurement uncertainty and experiment-to-model validation remain separate unresolved tasks. The current software never accesses physical hardware or asserts industrial performance.

这些结果用于可复现地说明模型如何响应所选范围。模型结构误差、参数可辨识性、测量不确定度及实验—模型验证仍是不同且尚未完成的问题。当前代码不访问真实硬件，也不宣称工业性能。
