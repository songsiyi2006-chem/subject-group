# Numerical dynamics and CI-NEB / 数值动力学与 CI-NEB

The reviewed module tests numerical algorithms on deliberately specified analytic potentials. It contains **no chemical reaction barrier, trained neural force field, DFT calculation, hardware operation, or experimental result**. Energy units are assigned consistently to the benchmarks; they do not make the model a molecular potential.

本模块使用明确给定的解析势能函数检验数值算法，**不包含真实化学反应势垒、已训练神经力场、DFT、硬件操作或实验结果**。基准中的能量与长度单位具有一致性，但这不意味着其势能函数已经代表真实分子。

## 1. Units and Velocity Verlet / 单位与 Velocity Verlet

The internal units are Å, fs, eV, and amu. The stored mass conversion is

\[
c_m=103.64269652680505\ {\rm eV\,fs^2\,Å^{-2}\,amu^{-1}},
\quad K=\frac12c_m\sum_i m_i|\mathbf v_i|^2,
\quad \mathbf a_i=\frac{\mathbf F_i}{c_m m_i}.
\]

This follows from the explicitly declared \(1\ {\rm amu}=1.66053906660\times10^{-27}\) kg and exact electron-volt conversion \(1.602176634\times10^{-19}\) J. The value of \(k_B\) is \(8.617333262145179\times10^{-5}\) eV/K. A unit test independently reproduces the force-to-acceleration conversion through SI units. The eV-to-kcal/mol conversion is 23.06054783061903; a potential-energy difference is not relabeled a free-energy barrier.

内部单位为 Å、fs、eV 和 amu。质量换算、动能和加速度公式如上。单元测试通过独立 SI 换算核对这些关系。eV 与 kcal/mol 的换算不会将势能差转化为自由能差。

The MD benchmark has eight equal-mass particles, each 12.011 amu, with

\[
U=\frac{k}{2}\sum_i|\mathbf r_i-\bar{\mathbf r}|^2,\qquad
k=2\ {\rm eV/Å^2},\qquad
\mathbf F_i=-k(\mathbf r_i-\bar{\mathbf r}).
\]

The center-of-mass coordinate and momentum are projected initially. All subsequent forces sum to zero; Langevin noise is also projected into the zero-COM subspace. The system therefore has \(f=3N-3=21\) unconstrained internal Cartesian degrees of freedom. No bond constraints or rigid-body rotational degrees are removed. Initial positions and velocities are independent Gaussian draws at 300 K in this subspace, with explicit seeds. The spring cluster has coincident-particle minima and is a mathematical oscillator benchmark, not a chemically plausible eight-atom molecule.

八个等质量粒子构成平移不变的谐振子基准。初始质心位置与动量均被消除，后续力和噪声也保持该子空间。因此自由度为 \(3N-3=21\)，没有额外扣除刚体转动或键约束。这个模型的最低能量构型是粒子重合，故只能用于数值检验，不能称作真实八原子分子力场。

The NVE integrator uses

\[
\mathbf r_{n+1}=\mathbf r_n+h\mathbf v_n+
\frac{h^2}{2}\mathbf a_n,\qquad
\mathbf v_{n+1}=\mathbf v_n+
\frac h2(\mathbf a_n+\mathbf a_{n+1}).
\]

An exact harmonic trajectory with \(\omega=\sqrt{k/(m c_m)}\) supplies an independent reference for final positions and velocities. Timesteps violating \(\omega h<2\) are rejected. Telemetry includes the true initial state at \(t=0\), post-update states at their actual times, and the final state. The Verlet family originates in the classical molecular-dynamics integration literature. [Verlet (1967)](https://journals.aps.org/pr/abstract/10.1103/PhysRev.159.98).

NVE 计算使用上述 Velocity Verlet 更新，并与解析谐振子轨迹比较。程序拒绝超过谐振子稳定界限的时间步；记录包含真实的零时刻和最终时刻，避免把更新后的状态标记到旧时间。

The full NVE study comprises four seeds × five timesteps, each lasting 2,000 fs. The energy quantity below is the maximum absolute total-energy change divided by that trajectory's initial total energy; it is an oscillatory numerical error envelope, not a fitted linear energy-drift slope.

| Timestep / 时间步 (fs) | Largest relative energy envelope among four seeds / 四种子的最大相对能量误差 | Largest final position RMS error / 最大终点坐标均方根误差 (Å) |
|---:|---:|---:|
| 0.25 | 0.0000073361 | 0.0000414194 |
| 0.50 | 0.0000293488 | 0.0001656753 |
| 1.00 | 0.0001174680 | 0.0006626605 |
| 2.00 | 0.0004710357 | 0.0026499698 |
| 4.00 | 0.0019028858 | 0.0105877504 |

The 16 paired step-halving estimates give energy-envelope convergence orders of 1.96457–2.04409 and final-position orders of 1.99041–2.00339. This verifies second-order behavior over the tested interval and initial conditions. It does not validate arbitrary force fields or timesteps.

20 条 NVE 轨迹在相同 2,000 fs 时长下比较五个时间步。16 组相邻步长比较的能量误差收敛阶为 1.96457–2.04409，终点坐标误差阶为 1.99041–2.00339，支持本基准范围内的二阶数值收敛。

## 2. Thermostat distribution check / 恒温分布检验

Three thermostat configurations each use eight matched initial seeds, a 0.5 fs timestep, a 50 fs relaxation time, 6,000 fs total duration, and 1,000 fs burn-in. Langevin uses BAOAB splitting with an exact Ornstein–Uhlenbeck velocity substep and COM-projected Gaussian noise. Berendsen uses deterministic weak velocity rescaling; the second Berendsen case deliberately adopts the source's \(3N\) temperature formula while the saved diagnostic temperature always uses internal kinetic energy and \(3N-3\).

三种恒温设置各使用八个对应的初始种子、0.5 fs 时间步、50 fs 耦合时间、6,000 fs 总时长和 1,000 fs 预热。Langevin 使用 BAOAB 分裂及质心投影噪声；Berendsen 使用确定性的弱耦合缩放。第三种设置故意采用原代码的 \(3N\) 温度公式，作为自由度处理不一致的受控对照。

For a classical canonical harmonic model, \(2K/(k_BT)\) follows a chi-square distribution with \(f=21\) degrees of freedom. Consequently, the reference internal-temperature variance at 300 K is \(2T^2/f=8571.428571\ {\rm K^2}\). The comparison averages within-trajectory sample variances across replicas. Adjacent time points are correlated; the variance ratios are descriptive finite-trajectory diagnostics, not calibrated confidence limits.

| Algorithm / 算法 | Mean internal temperature / 内部温度均值 (K) | Mean temperature variance / 温度方差均值 (K²) | Ratio to canonical variance / 与正则方差之比 |
|---|---:|---:|---:|
| Langevin BAOAB | 302.508228 | 8423.370489 | 0.982726557 |
| Berendsen, 21 DOF | 299.999705 | 0.002660 | 0.000000310355 |
| Berendsen, source 24 DOF | 342.856740 | 0.003125 | 0.000000364530 |

The nearly exact Berendsen mean does not establish canonical sampling: its fluctuations are strongly suppressed in this integrable spring benchmark. The deliberately incorrect thermostat formula drives the internal mean near \(300\times24/21=342.857143\) K. This numerical comparison is conditional on the chosen model and coupling. BAOAB also has finite-timestep bias and finite sampling error; the observed agreement of a variance does not establish general ergodicity. [Berendsen et al. (1984)](https://pure.rug.nl/ws/files/64380902/1.448118.pdf), [Leimkuhler and Matthews](https://arxiv.org/abs/1203.5428).

Berendsen 的温度均值接近 300 K，但该可积模型中的涨落被强烈压制，因此“均值正确”不能作为正则分布证据。使用原 \(3N\) 公式时，内部温度趋近 \(300\times24/21\)。BAOAB 也存在有限步长偏差和有限采样误差，本次方差接近解析值不代表任意体系中均能充分采样。

The source MD conversion from eV/Å to newtons is dimensionally correct. Its scientific limitations are separate: the neural force field is untrained; COM velocity is removed initially but laboratory-frame kinetic energy is divided by \(3N\); an external net force may later restore bulk drift; updated coordinates are timestamped one step early; and kinetic energy is recorded before thermostat rescaling while velocities have already been rescaled. Internal thermal temperature should consistently subtract instantaneous COM drift before applying \(3N-3\), rather than blindly changing a denominator for every externally driven trajectory.

原 MD 的 eV/Å 到 N 换算本身是正确的。需要指出的问题是：力场未训练；初始消除质心速度后仍以实验室坐标系动能除以 \(3N\)；外场净力可能重新产生整体漂移；时间戳落后一步；记录动能与温度发生在缩放前，而速度已经缩放。应先一致地去除瞬时整体漂移，再讨论内部温度自由度，不能只机械修改分母。

## 3. Convergence-tested CI-NEB / 按力判据收敛的 CI-NEB

The designed two-dimensional potential is

\[
V(x,y)=(x^2-1)^2+4[y-0.65(1-x^2)]^2.
\]

Coordinates have the chosen unit Å and energy has the chosen unit eV. The minima are \((-1,0)\) and \((1,0)\); the known saddle is \((0,0.65)\) with barrier 1 eV and Hessian eigenvalues \((-4,8)\) eV/Å². The curved valley floor is not assumed to be the MEP, since the force need not be tangent to that curve. Both endpoints are first minimized from perturbed guesses; their maximum final gradient norm is \(9.6144\times10^{-12}\) eV/Å.

二维解析势的两个极小值、鞍点及势垒均已知。谷底曲线不能直接当作 MEP：在非驻点处，势能梯度未必与曲线切向一致。两个端点先从扰动位置进行最小化，再固定于 NEB；最终端点最大力范数约为 \(9.61\times10^{-12}\) eV/Å。

The implementation uses the energy-weighted improved tangent, perpendicular physical force plus parallel spring force on ordinary images, and the climbing-image force

\[
\mathbf F_{\rm CI}=\mathbf F-
2(\mathbf F\cdot\hat{\boldsymbol\tau})\hat{\boldsymbol\tau}.
\]

Climbing starts after 50 ordinary-NEB updates. FIRE-style inertial relaxation has an adaptive algorithmic step and a capped image displacement; this optimization parameter is not physical MD time. Convergence requires both the largest interior NEB-force norm and the climbing image's true-force norm to fall below the specified threshold. Final energies are evaluated at the final stored coordinates. A hard iteration cap returns an explicit nonconverged status. The tangent and climbing projection follow the primary NEB methods. [Henkelman and Jónsson (2000)](https://henkelmanlab.org/pubs/henkelman00_9978.pdf), [Henkelman, Uberuaga and Jónsson (2000)](https://doi.org/10.1063/1.1329672).

实现包含能量加权改进切向、普通图像的垂直真实力和切向弹簧力，以及爬山图像的反向切向真实力。50 次普通 NEB 更新后开启爬山。收敛同时要求整条带的最大 NEB 力和爬山图像的真实力达标。最终能量对应最终坐标；超过迭代上限不会返回“已收敛”。

The main scan tests 7, 11, and 17 images, three perturbed initial-path seeds, and force thresholds of \(10^{-3}\) and \(10^{-5}\) eV/Å. All 18 cases converged and each final saddle had exactly one negative Hessian eigenvalue.

| Force threshold / 力阈值 (eV/Å) | Cases / 条件数 | Update range / 更新次数范围 | Largest saddle-position error / 最大鞍点位置误差 (Å) | Largest absolute barrier error / 最大势垒绝对误差 (eV) |
|---:|---:|---:|---:|---:|
| 0.001 | 9 | 104–188 | 0.000249073 | 0.000000123954 |
| 0.00001 | 9 | 153–251 | 0.00000219411 | 0.00000000000962830 |

These are convergence results for an analytic potential-energy saddle, not an activation free energy or a validated transition state for the supplied indole/Cu fragment. The one-negative-mode check is two-dimensional and does not cover molecular vibrational modes.

上述结果只说明解析势能面上的数值鞍点收敛，不对应原输入吲哚/Cu 片段的真实过渡态或活化自由能。负模检查也仅限这个二维基准，不包含真实分子的完整振动模式。

## 4. Unchanged-source comparison and reproducibility / 原类对照与复现

The original ClimbingImageNEB class is extracted by AST and executed unchanged against the same analytic potential. A return-frame trace captures its post-update image coordinates without altering its method. Twelve conditions combine 7/11/17 images with 1/10/45/150 requested updates. The source reports “converged” in every case, but **10 of 12 fail** the independent \(10^{-5}\) eV/Å force check. The 7- and 11-image, 150-update cases genuinely pass this check and are retained as positive controls.

从原文件通过 AST 提取 NEB 类，保持方法代码不变，仅将势能替换成同一已知解析函数。返回时捕获更新后的坐标进行独立检查。12 个条件均被原类标记为“收敛”，但其中 **10 个不满足**独立力判据；7 和 11 图像、150 次更新的两个条件确实达标，并未被隐去。

After one update, the source reports a pre-update barrier of 2.69 eV while the final coordinates give 1.781456 eV. The maximum stale energy discrepancy is 0.908544 eV. At 45 updates, the saddle energy is already close to 1 eV, yet the largest band-force residuals are approximately \(5.00\times10^{-5}\), \(8.19\times10^{-4}\), and \(6.50\times10^{-4}\) eV/Å for 7/11/17 images. Therefore an apparently accurate rounded barrier alone does not certify whole-band convergence. The source also labels potential-energy differences as free energies without entropy or solvation calculations.

一次更新后，原方法返回的势垒仍为更新前的 2.69 eV，而实际最终坐标给出 1.781456 eV；过期能量的最大差为 0.908544 eV。45 次更新时虽然鞍点能量已经接近正确值，但整条带的力残差仍未达标。因此小数位上“看起来正确”的势垒不能替代收敛检查；势能差也不能在缺乏熵与溶剂化处理时被称为自由能差。

The full study executed 44 MD trajectories, 350,000 integration steps, 350,044 MD energy/force evaluations, 18 reviewed NEB cases with 3,129 updates and 37,337 image-point evaluations, and 12 source-class cases with 618 updates and 7,210 point evaluations. Endpoint checks and source post-update audits add 36 and 140 analytic point evaluations, respectively; two endpoint minimizations use 778 evaluations. These categories are recorded separately. Pilot and unit-test runs are excluded from these full-study counts.

全量计数为 44 条 MD 轨迹、350,000 积分步、350,044 次 MD 能量/力评估；18 个审查版 NEB 条件合计 3,129 次更新、37,337 个图像点评估；12 个原类条件合计 618 次更新、7,210 个点评估。另有 36 个 NEB 端点检查、140 个原类更新后审计点和两次端点最小化中的 778 次评估，各类别分列，pilot 与测试不计入全量研究。

The public entry point is [dynamics_reviewed.py](scripts/dynamics_reviewed.py). From the repository root:

    python synthapore/scripts/dynamics_reviewed.py --pilot
    python synthapore/scripts/dynamics_reviewed.py
    python -m unittest discover -s tests -p test_synthapore_dynamics.py -v

Saved artifacts include [summary and provenance](results/dynamics/summary.json), [all NVE summaries](results/dynamics/nve_timestep_summary.csv), [all thermostat replica summaries](results/dynamics/thermostat_replicas.csv), [all NEB force histories](results/dynamics/neb_force_history.csv), [all final bands](results/dynamics/neb_final_bands.csv), and [original-class audit](results/dynamics/source_neb_summary.csv). MD CSV telemetry is saved for the first seed at every tested timestep or thermostat configuration; all NVE initial/final states and every replica's aggregate statistics are saved. Full statistical reproduction uses the stored seeds and source. [Twenty tests](results/dynamics/test_results.txt) cover units, reversibility, COM bookkeeping, analytic gradients and Hessians, force convergence, endpoint treatment, censored iteration limits, and the source defects.

所有 NEB 收敛历史、最终图像及每条 MD 轨迹的统计汇总均已保存；MD 时间序列 CSV 保存各步长或恒温配置的第一个种子，全部 NVE 初末状态也单独保存。完整统计复现使用明确记录的种子与代码。20 项单元测试检验的是数值实现，不代替材料、分子反应或工业级有效性的证据。
