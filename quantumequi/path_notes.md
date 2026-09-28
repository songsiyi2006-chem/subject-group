# QuantumEqui-NEB path audit and numerical validation / 路径审查与数值验证

## Evidence boundary / 证据边界

The supplied class was replayed without editing its algorithm, using the exact weights captured from the original run. Those weights have not been trained. Its energy and force numbers retain the source's **nominal kcal/mol and kcal/(mol Å)** convention; this label is not a calibration. The analytic surfaces below use **dimensionless coordinates, energies and forces**. Neither calculation establishes a molecular activation barrier or free energy.

原样回放使用原始运行捕获的确切权重，没有改动源 NEB 类。权重未经训练，其能量和力仅保留源代码的**名义单位**，不能视为经过校准的化学能量。下面两个解析势面使用**无量纲坐标、能量和力**。这些结果验证算法行为，不构成分子活化能、自由能或电催化机理证据。

## Exact source replay / 确切源程序回放

The source evaluates all seven image energies and forces before moving any image. It then updates interior images in increasing index order, modifying the same array. Consequently image i uses the already updated left neighbor when computing its tangent and spring term, while its true force still belongs to the old geometry. This is a sequential, order-dependent update. Source climbing begins at zero-based iteration 15. There is no force stopping rule: the method executes 30 loops and returns status “converged” unconditionally.

源程序先在旧路径上计算全部七个图像的能量与力，再按编号递增原地更新内部图像。因此第 i 个图像计算切线和弹簧力时，其左邻居已经更新，但真实力仍来自旧坐标。这会引入更新顺序依赖。攀爬从零起算第 15 轮开始；程序不检查力收敛，执行 30 轮后直接报告 converged。

The replay wrapper records the pre-update arrays from the caller frame. A return-profile hook captures the complete final band and final pre-update energy array; the previous profiler is restored afterwards. All 30 updates in each direction are reconstructed algebraically from the saved geometry, energies and forces, with exactly zero coordinate difference. The forward return dictionary, including candidate coordinates, equals the original captured dictionary. No separately reseeded model substitutes for the captured weights.

回放包装器从调用帧记录更新前坐标，通过返回事件捕获最终完整路径及陈旧能量数组，退出后恢复原有 profiler。正、反向各 30 轮更新均由保存数组代数重建，坐标差恰为零。正向返回字典（含候选坐标）与原始捕获结果完全一致，没有用重新设定随机种子的近似模型替代确切权重。

| Diagnostic / 诊断 | Forward / 正向 | Reverse / 反向 |
|---|---:|---:|
| Candidate image index / 候选编号 | 1 | 5 |
| Candidate maximum atomic true force, nominal / 候选最大原子真实力 | 0.0183398705 | 0.0183399059 |
| Final maximum projected atomic force, nominal / 最终最大投影原子力 | 0.0200995542 | 0.0200964718 |
| Maximum stale energy difference, nominal / 最大陈旧能量差 | 0.0000625849 | 0.0000625849 |
| Fresh candidate relative energy, nominal / 复算候选相对能量 | −0.0108846426 | 0.0249454975 |
| Fresh maximum over all images relative to start, nominal / 全路径相对起点最高能量 | 0 | 0.0358272791 |
| Minimum final pair distance, Å / 最终最小原子间距 | 0.4804261244 | 0.4804131927 |

Both endpoint maximum atomic forces are nonzero: 0.0180063229 and 0.0197831113. Thus neither prescribed endpoint has been demonstrated to be a local minimum. The source excludes endpoints when selecting its highest-energy candidate; this explains its negative reported forward “barrier” of −0.01. It is not evidence of a negative physical activation barrier. The candidate's nonzero gradient must be considered before interpreting any Hessian index.

两个端点最大原子力分别为 0.0180063229、0.0197831113，尚未证明为极小值。源程序选候选时排除端点，所以正向可能报告 −0.01 的“势垒”；这不是负物理活化能的证据。候选结构的梯度非零，不能仅凭 Hessian 本征值数量确认过渡态。

Reversing the endpoints and mapping the images back gives a maximum coordinate difference of **0.000273761958 Å** after the same 30 loops. This compares finite iterates, not two converged reaction barriers. A same-geometry, same-force counterfactual using simultaneous updates differs from the actual sequential step by up to **0.000049466455 Å**. Forward and reverse physical barriers would in general differ by the reaction energy; unequal barriers alone are not the order-dependence diagnostic.

端点反转并映射回原图像顺序后，30 轮终点坐标最大差为 **0.000273761958 Å**。这是有限迭代路径的比较，不是两个已收敛化学势垒的比较。在完全相同旧坐标和力下，同步更新与实际顺序更新的单步最大差为 **0.000049466455 Å**。正反向物理势垒本来就可能相差反应能，因此不能只根据势垒不同判断顺序偏差。

## Reviewed solver and new surfaces / 审查版求解器与新势面

The callback-based solver generalizes the existing [SynthaPore implementation](../synthapore/scripts/dynamics_reviewed.py). Its energy-weighted tangent follows [Henkelman and Jónsson](https://doi.org/10.1063/1.1323224), its climbing force follows [Henkelman, Uberuaga and Jónsson](https://hj.hi.is/papers/paperCI-NEB.pdf), and its inertial relaxation is a FIRE-style implementation motivated by [Bitzek et al.](https://doi.org/10.1103/PhysRevLett.97.170201). This is an explicit reuse and extension, not a new optimization method.

回调式求解器复用并泛化已有 SynthaPore 代码，采用能量加权切线、CI 力投影和 FIRE 风格惯性松弛。这是已有方法的复用与扩展，不宣称新的优化算法。

Every force is projected from one unchanged band snapshot; all images are then updated simultaneously. Fixed endpoints must satisfy the requested force threshold before optimization. Climbing starts at iteration 50. Convergence requires both the largest interior image force norm and the climbing-image true force norm to be at most the threshold. Reaching the 5,000-update cap returns a failed convergence status. Fresh energies always correspond to the saved final coordinates. The artificial optimizer step parameter is not a time step.

每轮从同一个路径快照计算全部投影力，再同步更新图像。固定端点必须先满足请求的力阈值。第 50 轮起攀爬；只有最大内部图像力和攀爬图像真实力同时低于阈值，才报告收敛。达到 5,000 次更新上限仍不收敛时返回失败。最终能量与最终坐标一致，优化步长参数不代表物理时间。

The first surface has unequal well depths:

\[
V_A(x,y)=(x^2-1)^2+0.35x+
2.5\left[y-0.55\sin(1.7x)\right]^2.
\]

Stationary x coordinates are the three real roots of \(4x^3-4x+0.35=0\), with \(y=0.55\sin(1.7x)\). The second benchmark tests a nonplanar three-dimensional path between adjacent minima of a periodic coordinate:

\[
V_B(x,y,z)=2(1-\cos x)+2(y-0.45\sin x)^2+
3\left[z-0.30(1-\cos x)\right]^2.
\]

Its endpoints are \((0,0,0)\) and \((2\pi,0,0)\), and its saddle is \((\pi,0,0.6)\). The valley functions define stationary reference points; they are not asserted to trace the entire minimum-energy path.

第一个势面两端井深不等，驻点横坐标由三次方程的三个实根给出。第二个势面具有非共面三维路径，连接周期坐标上相邻的极小值。解析谷底用于构造驻点参考，不把整条谷底曲线预先当作最小能量路径。

| Reference / 解析参考 | Tilted sine 2D / 非对称二维 | Periodic curve 3D / 三维 |
|---|---:|---:|
| Forward barrier / 正向势垒 | 1.3727196610 | 4 |
| Reverse barrier / 反向势垒 | 0.6733964499 | 4 |
| Reaction energy / 反应能差 | 0.6993232111 | 0 |
| Saddle Hessian eigenvalues / 鞍点 Hessian 本征值 | −2.487062678, 7.854002266 | −1.753167982, 4.563167982, 6 |

## Executed convergence study / 已执行的收敛计算

Each surface uses 7, 11 and 17 images, force thresholds \(10^{-3}\) and \(10^{-5}\), and seeds 11, 22 and 33: **36 sweep cases**. The initial linear band receives Gaussian interior perturbations of amplitude 0.06 multiplied by a sine envelope; endpoints remain fixed. The spring constant is 2. All 36 cases converge and all candidate analytic Hessians have index one.

每个势面采用 7／11／17 个图像、两个力阈值、三个种子，共 **36 个扫描条件**。线性初始路径的内部图像加入振幅 0.06、乘正弦包络的高斯扰动，固定两端。弹簧常数为 2。36 个条件均按力阈值收敛，候选解析 Hessian 均仅含一个负本征值。

| Surface / 势面 | Force threshold / 力阈值 | Cases / 个数 | Update range / 更新数范围 | Max barrier absolute error / 最大势垒绝对误差 | Max saddle position error / 最大鞍点位置误差 |
|---|---:|---:|---:|---:|---:|
| Tilted sine 2D | 0.001 | 9 | 147–300 | 1.65750e−7 | 3.65261e−4 |
| Tilted sine 2D | 0.00001 | 9 | 201–364 | 1.43074e−11 | 3.39351e−6 |
| Periodic curve 3D | 0.001 | 9 | 112–227 | 2.10034e−7 | 4.91792e−4 |
| Periodic curve 3D | 0.00001 | 9 | 183–290 | 1.52784e−11 | 4.19447e−6 |

Four additional reviewed runs reverse exactly the same 11-image unperturbed initial band. The mapped maximum coordinate differences are \(1.11\times10^{-16}\) and \(2.22\times10^{-16}\) for the 2D and 3D surfaces. These are deterministic numerical controls on these surfaces, not a general molecular guarantee.

另有四次审查版计算反转完全相同的 11 图像无扰动初始路径。映射坐标差分别为 \(1.11\times10^{-16}\)、\(2.22\times10^{-16}\)。这仅是两个解析势面的数值对照，不是对所有分子路径的普遍保证。

## Workload and reproduction / 工作量与复现

| Main workload / 正式工作量 | Count / 数量 |
|---|---:|
| Exact source replay directions / 原样源回放方向 | 2 |
| Source loops / 源循环 | 60 |
| Source energy–force calls / 源能量力调用 | 420 |
| Fresh final source image calls / 源最终图像复算调用 | 14 |
| Reviewed analytic sweep runs / 审查版解析扫描 | 36 |
| Sweep updates / 扫描更新 | 7,774 |
| Sweep band point energy–force evaluations / 扫描路径点评估 | 95,904 |
| Sweep endpoint checks / 扫描端点检查 | 72 |
| Reviewed reversal runs / 审查版反转计算 | 4 |
| Reversal updates / 反转更新 | 748 |
| Reversal band point evaluations / 反转路径点评估 | 8,272 |
| Reversal endpoint checks / 反转端点检查 | 8 |
| Analytic reference point energy–force evaluations / 解析参考点评估 | 6 |
| Analytic candidate Hessians plus reference Hessians / 候选与参考解析 Hessian | 36 + 6 |

Thus the reviewed main calculation has **40 runs, 8,522 updates and 104,262 point energy–force evaluations**, including endpoint and reference checks. These are batched analytic point evaluations, not neural forward passes. The separate pilot has 2 sweep cases and 4 reversal cases (11,449 point evaluations including its 12 endpoint and 6 reference checks). Pilot, original master-script work, test work and this module's two exact source replays are reported separately, and are not added to the analytic main total.

正式审查版共 **40 次求解、8,522 次更新、104,262 次点评估**，最后一项包含端点和解析参考点。这些是可批量执行的解析势点评估，不是神经网络前向次数。先导计算、原始主脚本、测试及本模块的确切源回放分别计数，不能混入正式解析总数。

The **20 tests** cover derivative identities, stationary references, unequal forward/reverse barriers, three-dimensional geometry, tangent reversal, projection, spring-free climbing, force termination, fixed endpoints, callback coordinate shapes, sequential update bias and exact saved-source reconstruction. They do not retrain a model or rerun the full study.

20 项测试覆盖导数恒等式、驻点、非对称势垒、三维几何、切线反转、力投影、攀爬去弹簧、收敛退出、固定端点、回调形状、顺序偏差及确切源记录重建；不会重新训练模型或重跑完整研究。

Reproduction, using the existing one-thread chemistry environment:

    python quantumequi/scripts/path_reviewed.py --source-audit
    python quantumequi/scripts/path_reviewed.py --pilot
    python quantumequi/scripts/path_reviewed.py
    python -m unittest discover -s tests -p test_quantumequi_path.py -v

Saved outputs: [analytic cases](results/path/analytic_cases.csv), [all sweep force histories](results/path/analytic_force_history.csv), [all final analytic images](results/path/analytic_final_images.csv), [analytic reference structures](results/path/analytic_references.json), [exact source audit](results/path/source_audit.json), [source final coordinates and forces](results/path/source_final_coordinates.json), [source final image diagnostics](results/path/source_final_images.csv), [tests](results/path/test_results.txt). The source NPZ archives contain **all 30 pre-update bands, energies and force vectors in each direction**. Source and output hashes, versions and scope are stored in the corresponding JSON summaries.

源 NPZ 文件保存两个方向各 30 轮的全部更新前路径、能量及力矢量。相关 JSON 同时记录版本、源文件与输出哈希。程序通过与数值收敛不能替代势能模型训练、分子端点优化、独立电子结构证据或实验验证。
