# Supervised distance-message learning extension / 距离消息网络监督学习扩展

## Scope and methods

This extension trains a genuine small neural energy function using the frozen, actually computed H2 RHF/STO-3G reference. It performs **zero new quantum jobs**. It does not train the original random Cu-fragment EGNN or use new FCI, DFT, Cu, or electrode labels. The only species is hydrogen; the supported input is a batch of two-atom structures. Because H2 has one distinct interatomic distance, the graph calculation reduces to radial regression. No molecular-graph or chemical-space transfer claim follows.

[`scripts/learning_extension.py`](scripts/learning_extension.py) defines `HydrogenDistanceMPNN`: a shared learned hydrogen embedding, two width-16 message/update layers, SiLU nonlinearities, and a summed scalar atomic readout, totaling **3,537 parameters**. Each of the two directed edges sends a scalar-feature message based on the two node states and standardized distance. There are no self edges, coordinate-update head, cutoff, pretrained weights or hidden quantum calls. For layer l,

\[
m_{ij}^{(l)}=\phi_l(h_i^{(l)},h_j^{(l)},(r_{ij}-\bar r_{train})/s_{r,train}),\quad
h_i^{(l+1)}=h_i^{(l)}+\psi_l(h_i^{(l)},\sum_{j\ne i}m_{ij}^{(l)}).
\]

The physical-unit energy is the sum of two atomic outputs, each with offset Ētrain/2 and scale sE,train. Cartesian forces are **−∂E/∂X**, computed by coordinate autograd. For nuclei at z=±R/2, dE/dR=(F1z−F2z)/2. The gradient-supervised objective differentiates these forces with `create_graph=True`, enabling parameter training through the coordinate derivative. This is a compact original implementation of established ideas, not a SchNet reproduction; method attribution is to [Gilmer et al.](https://proceedings.mlr.press/v70/gilmer17a.html) and [Schütt et al.](https://proceedings.neurips.cc/paper/2017/hash/303ed4c69846ab36c2904d3ba8573050-Abstract.html). The [PyTorch autograd API](https://docs.pytorch.org/docs/stable/generated/torch.autograd.grad.html) documents the derivative mechanism.

The original fixed partition is preserved exactly: **17 train, 8 validation, 8 test, 9 stretched OOD points**. All preprocessing uses the training partition only: mean distance 1.14 Å, distance standard deviation 0.3919183588453085 Å, energy mean −1.0010483851947871 Hartree, energy standard deviation 0.09370460709562088 Hartree, and gradient standard deviation 0.3187903479156144 Hartree/Å. Population standard deviations use ddof=0. No displaced finite-difference labels, new quantum points or test/OOD labels enter coefficient fitting or model selection.

The paired objectives are

\[
L_E=\langle[(\widehat E-E)/s_E]^2\rangle,\qquad
L_{E,g}=L_E+\langle[(\partial_R\widehat E-g)/s_g]^2\rangle.
\]

Both select checkpoints using validation RMSE(E)/sE + RMSE(g)/sg. Thus **energy-only refers to the fitting loss**: validation gradients and the training gradient scale are still accessed for checkpoint selection. This matches the previous RBF comparison's label-access rule; it is not a comparison between complete absence and presence of force labels. No validation/test refitting occurs.

A two-run pilot used seed 7301 and 200 epochs per objective, evaluating only train/validation. Pilot validation energy/gradient RMSE was 0.0190887/0.261273 for energy-only and 0.000707781/0.0153631 for joint fitting, in Hartree and Hartree/Å. The predeclared width, optimizer and 800-epoch main plan were retained; test/OOD scores were not consulted to tune them. The main study uses seeds **7301,7302,7303**, each shared between both objectives. Initial-state hashes verify exact matching. Every run uses full batches of 17, float64, Adam with zero weight decay, cosine learning-rate decay from 0.003 to 0.00015, gradient-norm cap 5, and exactly 800 optimization steps. All six selected checkpoints are at epoch 800. This budget boundary is not evidence that optimization has converged.

## Results and retained negative findings

| Model / seed | Test E RMSE / Hartree | Test gradient RMSE / Hartree Å⁻¹ | OOD E RMSE / Hartree | OOD gradient RMSE / Hartree Å⁻¹ |
|---|---:|---:|---:|---:|
| Neural energy-only / 7301 | 0.000486462 | 0.0121966 | 0.00767735 | 0.0147389 |
| Neural joint / 7301 | 0.000128189 | 0.00148374 | 0.000173536 | 0.00112220 |
| Neural energy-only / 7302 | 0.000757399 | 0.0121112 | 0.0277884 | 0.0571914 |
| Neural joint / 7302 | 0.000216394 | 0.00159951 | 0.00812651 | 0.0212530 |
| Neural energy-only / 7303 | 0.000226936 | 0.00509392 | 0.0159985 | 0.0406077 |
| Neural joint / 7303 | 0.0000902290 | 0.000771657 | 0.00596255 | 0.0161202 |
| Frozen RBF energy-only | 0.00000862425 | 0.000115617 | 0.175473 | 0.459351 |
| Frozen RBF joint | 0.00000556135 | 0.000163747 | 0.118352 | 0.448587 |
| Frozen linear interpolation/extrapolation | 0.00110348 | 0.00223930 | 0.0444501 | 0.108744 |

All three matched seeds improve under gradient supervision on these test/OOD metrics. Nevertheless, every neural run is worse than both frozen RBFs for test energy and gradient. All energy-only runs also have worse test energy than the frozen cubic spline. The neural models extrapolate better than the frozen RBF and linear models on this particular RHF stretched set, but joint OOD energy RMSE varies from **0.000173536 to 0.00812651 Hartree** across seeds. There is no evidence that the best individual seed will generalize to another task. These errors measure approximation of the old RHF/STO-3G labels; favorable interpolation or extrapolation cannot correct the RHF dissociation approximation.

For each objective, `ensemble_predictions.csv` retains 42 per-point reference values, three-seed means and sample standard deviations (ddof=1). The joint ensemble gives test energy/gradient RMSE **0.000142893/0.00119298**, and OOD **0.00473607/0.0127819**. Mean joint test spread is 0.0000647888 Hartree and 0.000680392 Hartree/Å. Only 75% of the eight test absolute errors fall within twice the sample spread for either quantity. This is a descriptive spread-versus-error statistic, not a confidence interval, coverage guarantee or calibrated uncertainty model. OOD spread/error Spearman correlations equal one on nine ordered distances, which may reflect a shared distance trend; correlated geometries and only three seeds preclude broad calibration claims. All per-point and split summaries are retained.

Mathematical checks apply to each selected learned model at R=0.74,1.14,2.30 Å. Combined orthogonal transformations, both determinant signs, translation and atom exchange yield 36 probe rows. Maximum energy-invariance and force-covariance errors are each **2.22×10⁻¹⁶** in their respective units; net force and torque vanish to the recorded precision. Four complete Cartesian finite differences per geometry yield 72 rows. Across all models/geometries, maximum force discrepancies at h=0.01,0.001,0.0001,0.00001 Å are **2.24512×10⁻⁴, 2.24526×10⁻⁶, 2.24528×10⁻⁸, 2.22813×10⁻¹⁰ Hartree/Å**. These verify covariance and conservative derivatives of the learned function, not quantum accuracy at unobserved molecules.

## Counts, retention and reproduction

| Work | Main | Separate pilot |
|---|---:|---:|
| Training runs | 6 | 2 |
| Optimizer steps / validation epochs | 4,800 / 4,800 | 400 / 400 |
| Python model forward batches | 9,732 | 802 |
| Molecular energy evaluations within batches | 121,170 | 10,050 |
| Per-seed prediction rows | 252 | 50 |
| Per-seed metric rows | 24 | 4 |
| Covariance / Cartesian FD cases | 36 / 72 | 0 / 0 |
| New quantum jobs | 0 | 0 |

Each main run has 800 training and 800 validation forwards, covering 20,000 molecular evaluations, followed by 42 selected-checkpoint predictions and 153 geometry evaluations for mathematical checks. Thus each run has 1,622 forward batches and 20,195 molecular evaluations. A forward batch is not a quantum calculation or an optimizer epoch. For training/validation, energy-only performs 800 coordinate-gradient batches; joint fitting performs 1,600. These counts exclude unit tests and the independent earlier quantum calculations. Per-run fit times range 5.75–8.04 s; total main wall time inside the script is 41.79 s, pilot 3.97 s, both single CPU thread. Interpreter startup is additional.

All six main checkpoint archives retain initial, best-validation and final states, plus states every 200 epochs. `learning_curves.csv` contains every epoch, with training loss measured before the update and validation metrics after it. Matching initial hashes, selected hashes, code/data hashes, architecture and preprocessing travel with each checkpoint. Main code SHA256 is `a36c3945ea015eb9d9a011b10120454fed3643e063367d761ceb6a3d13aae33f`; main and pilot use matching source snapshots. All 16 main and seven pilot output hashes were checked. **14 unit tests pass**, covering frozen partitions, no preprocessing leakage, parameter training through gradient loss, complete Cartesian differences, rotation/reflection/translation/permutation covariance, force/torque balance, matched initialization, saved-state independence, sample spread and invalid geometry rejection. Tests do not rerun the study.

Key APIs: `read_reference`, `fit_preprocessing`, `HydrogenDistanceMPNN`, `energy_forces`, `radial_predictions`, `fit_one`, `evaluate_rows`, `mathematical_audit`, `ensemble_tables`. Plot/join schemas:

- `predictions.csv`: objective, seed, point_id, R_A, split, reference/predicted energy and gradient, signed errors.
- `ensemble_predictions.csv`: objective, point_id, R_A, split, reference energy/gradient, ensemble means, sample standard deviations and absolute errors.
- `learning_curves.csv`: objective, seed, epoch, pre-update losses, post-update validation metrics and common selection score.
- `metrics.csv`, `ensemble_diagnostics.csv`: split errors and uncalibrated spread diagnostics.
- `frozen_baseline_metrics.csv`: unchanged numerical values copied from the earlier baseline table, without refitting or new quantum work.

For deliberate regeneration, use a separate checkout with `CHEM_PYTHON` pointing to the existing chemistry interpreter and native runtime configured. Existing completed output is rejected to protect the frozen study. Environment, timing and hashes are recorded in `results/extensions/learning/summary.json`.

```powershell
& $env:CHEM_PYTHON quantumequi/scripts/learning_extension.py --pilot
& $env:CHEM_PYTHON quantumequi/scripts/learning_extension.py
& $env:CHEM_PYTHON -m unittest discover -s tests -p test_quantumequi_learning_extension.py -v
```

## 中文方法、结果与边界

本扩展使用此前实际计算并冻结的 H2 RHF/STO-3G 数据，训练真正的小型神经能量模型，**没有新增量子作业**。没有训练原始 Cu 片段的随机 EGNN，也没有使用新增 FCI、DFT、Cu 或电极标签。模型只支持两个氢原子组成的批量输入。H2 只有一个独立原子间距，因此图消息传递实际上退化为径向回归，不能证明跨分子图或化学空间迁移。

网络包含共享可学习氢嵌入、两层 width=16 的距离消息/节点更新、SiLU 激活和原子能求和，共 **3,537 个参数**。两条有向边使用节点状态与训练标准化后的距离生成消息，没有自边、坐标更新头、预训练参数或隐藏量子调用。能量偏移均分到两个氢原子，输出尺度来自训练能量标准差。力由坐标 autograd 给出 **−∂E/∂X**；对 z=±R/2 的两个原子，dE/dR=(F1z−F2z)/2。联合目标通过 `create_graph=True` 对坐标梯度继续求参数导数。实现借鉴上述 MPNN/SchNet 的通用思想，但不是完整复现 SchNet。

严格保留 **17 训练、8 验证、8 测试、9 拉伸分布外** 的原始分区。预处理只使用训练数据：距离均值 **1.14 Å**、标准差 **0.3919183588453085 Å**，能量均值 **−1.0010483851947871 Hartree**、标准差 **0.09370460709562088 Hartree**，梯度标准差 **0.3187903479156144 Hartree/Å**，均采用总体标准差 ddof=0。位移核验点、新量子点和测试/分布外标签均不进入参数拟合与选模。

仅能量目标使用标准化能量 MSE；联合目标再加权重为一的标准化梯度 MSE。两类模型均按“验证能量 RMSE/训练能量标准差 + 验证梯度 RMSE/训练梯度标准差”选择检查点。因此，**仅能量是指拟合损失**，其选模仍接触验证梯度及训练梯度尺度；不能宣称完全没有使用力标签。验证/测试数据没有用于重新拟合。

先导计算为种子 **7301**、两目标各 **200 轮**，只评估训练与验证。验证能量/梯度 RMSE 分别为仅能量 **0.0190887/0.261273**、联合 **0.000707781/0.0153631**，单位为 Hartree、Hartree/Å。保留预先计划的网络与 800 轮正式配置，没有根据测试或分布外得分调参。正式种子为 **7301、7302、7303**，每个种子两目标的初始状态哈希完全一致；使用全批量 17、float64、无权重衰减 Adam、从 0.003 到 0.00015 的余弦学习率以及梯度范数上限 5。六次均选中第 **800 轮**，到达预算边界不证明优化已经收敛。

上方完整表格保留了所有六次训练与冻结基线的数值。三个配对种子在联合训练后，测试与分布外误差均降低；但是**所有神经模型的测试能量和梯度仍比两个冻结 RBF 更差**。全部仅能量训练的测试能量也差于三次样条。针对这组 RHF 拉伸标签，神经网络的外推优于原 RBF 与线性基线，但联合分布外能量 RMSE 在 **0.000173536–0.00812651 Hartree** 间明显波动。不能选择表现最好的种子后宣称普适优势，拟合 RHF 的精度也不能消除 RHF 的解离近似误差。

逐点集成表给出三种子均值和样本标准差 ddof=1。联合集成测试能量/梯度 RMSE 为 **0.000142893/0.00119298**，分布外为 **0.00473607/0.0127819**。测试能量/梯度平均离散度为 **0.0000647888 Hartree、0.000680392 Hartree/Å**，两项均只有 **75%** 的八个测试误差落在两倍样本标准差内。这只是未经校准的离散度—误差统计，不是置信区间或覆盖率保证。分布外九个有序距离上 Spearman 相关均为一，可能只是随距离共同变化；三个种子及高度相关几何不足以验证不确定性校准。

每个已选模型在 **0.74、1.14、2.30 Å** 上进行数学检查。旋转/反射、平移及原子交换共 **36 条**记录，最大能量不变性与力协变误差均为 **2.22×10⁻¹⁶**。四个位移、全六个笛卡尔分量的中心差分共 **72 条**记录；位移从 0.01 降至 0.00001 Å 时，最大力差从 **2.24512×10⁻⁴** 降至 **2.22813×10⁻¹⁰ Hartree/Å**。这些验证模型的对称性和保守导数，不代表其他分子的量子精度。

正式计算为 **6 次训练、4,800 次优化、9,732 个前向批次、121,170 次批内分子能量评估**；先导另计两次训练、400 次优化、802 个批次、10,050 次分子评估。正式模型逐次评估 20,000 个训练/验证几何、42 个输出点及 153 个数学检查几何；批次、分子评估与量子作业是不同计数。主程序内总耗时 **41.79 秒**，先导 **3.97 秒**，单 CPU 线程。14 项测试全部通过，测试工作不计入正式研究。

六个主检查点归档均包含初始、验证最优、最终以及每 200 轮状态。学习曲线完整保存每轮更新前训练损失、更新后验证指标。主/先导分别 **16/7 个输出文件哈希通过**，源代码快照与数据来源可追溯。复现使用上述便携环境变量与命令，已完成结果禁止覆盖；应在独立检出目录主动重算，并重新生成 QA。所有新参考或跨方法误差分解必须在训练之外独立开展。
