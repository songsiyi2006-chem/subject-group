# Nuclear vibration extension / 核振动扩展

This extension solves a one-dimensional, nonrotating nuclear model and compares harmonic and anharmonic isotope effects. It does not calculate experimental spectral intensities or a full molecular Gibbs energy. The numerical control uses declared mathematical parameters; the chemical models use Morse parameters determined from the separately executed H2 FCI calculations. A Morse model parameterized by FCI data is not nuclear propagation on the actual FCI potential curve.

本扩展求解一维、无转动的核运动模型，比较谐振与非谐同位素效应。不计算实验光谱强度或完整分子 Gibbs 自由能。数值对照采用预先声明的数学参数；化学模型的 Morse 参数来自独立执行的 H2 FCI 计算。由 FCI 数据参数化的 Morse 模型不等于直接在实际 FCI 势能曲线上求解核运动。

## Method and units / 方法与单位

The J = 0 radial equation is

\[
\left[-c_\mu\frac{d^2}{dr^2}+V(r)\right]u_v(r)=E_vu_v(r),\qquad
c_\mu=\frac{\hbar^2}{2\mu},\qquad
\mu=\frac{m_1m_2}{m_1+m_2}.
\]

Distances are in angstrom, energies in Hartree above the potential minimum, and wavefunctions in angstrom to the power −1/2. The implementation converts the kinetic coefficient through SI units to Hartree angstrom squared. Bare proton and deuteron masses are 1.0072764665789 and 2.013553212544 u; the corresponding homonuclear reduced masses are one half of those values. These are nuclear masses, not natural isotope-average atomic weights. Both isotopes share one Born–Oppenheimer electronic potential; diagonal Born–Oppenheimer, nonadiabatic, relativistic and QED corrections are omitted. Constants follow [NIST CODATA 2022](https://physics.nist.gov/cuu/Constants/Table/allascii.txt).

距离使用 Å，能量使用相对势能最低点的 Hartree，波函数单位为 Å^−1/2。动能系数经 SI 单位换算为 Hartree·Å²。质子和氘核质量分别为 1.0072764665789、2.013553212544 u，同核双原子分子的约化质量为相应核质量的一半。这里不使用天然同位素平均原子量。两种同位素共享 Born–Oppenheimer 电子势能面，未计入对角 Born–Oppenheimer、非绝热、相对论及 QED 修正；常数来源同上。

The radial endpoints impose u(left) = u(right) = 0. A second-order central finite difference gives diagonal V + 2cμ/Δr² and adjacent entries −cμ/Δr². M intervals contain M−1 interior unknowns, while M+1 positions include the two fixed boundary values. This is a finite-difference discretization, not a claimed high-order DVR. [SciPy's symmetric tridiagonal eigensolver](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh_tridiagonal.html) selects eigenvalues below the known model dissociation energy, with a 10^−12 Hartree exclusion margin. Quadrature normalization, orthogonality, and matrix residuals are recorded separately from discretization errors.

径向区间两端施加零边界条件。二阶中心差分产生三对角矩阵，对角元为 V + 2cμ/Δr²，紧邻非对角元为 −cμ/Δr²。M 个区间对应 M−1 个内部未知数，总坐标数 M+1 包括两个固定边界点。本程序采用有限差分，不声称使用高阶 DVR。求解器仅返回低于已知模型解离阈值的本征值，并预留 10^−12 Hartree 排除余量。波函数积分归一化、正交性和矩阵本征残差，与离散化误差分开记录。

The model potential and its exact full-real-line bound spectrum are

\[
V(r)=D_e[1-e^{-a(r-r_e)}]^2,\qquad
w=2a\sqrt{D_ec_\mu},\qquad
E_v=w(v+1/2)-\frac{w^2}{4D_e}(v+1/2)^2,
\]

with the strict bound-state condition v + 1/2 < 2De/w. Here w is an energy spacing, not an angular frequency in s^−1. Division by hc converts it to cm^−1. The harmonic zero-point energy is w/2, while the model dissociation energy from v = 0 is D0 = De−E0. These expressions follow the exactly soluble [Morse potential](https://journals.aps.org/pr/abstract/10.1103/PhysRev.34.57). A finite radial domain has different boundary conditions from the full-real-line analytical problem; agreement is conditional on negligible wall effects.

其中 w 表示能量间隔，而不是单位为 s^−1 的角频率；除以 hc 后转换为 cm^−1。谐振零点能为 w/2，振动基态到解离阈值的模型解离能为 D0 = De−E0。严格束缚条件为 v + 1/2 < 2De/w。上述解析式来自 Morse 势的精确解；有限径向区间与全实轴解析模型具有不同边界条件，二者只在壁面影响可忽略时接近。

## Predeclared numerical controls / 预先定义的数值对照

The mathematical control uses De = 0.2 Hartree, a = 1.8 Å^−1, re = 1 Å. These are assumptions, not measured or computed H2 parameters. The main grid study uses 400, 800 and 1200 intervals over [0, 8] Å. Independent domain checks use right endpoints 2, 4, 6 and 8 Å at a common spacing of 1/150 Å. The duplicate 1200-interval, 8 Å calculation is reused and counted once. First-three-state errors, observed second-order convergence, and Richardson estimates are retained, alongside all shared-state errors and the probability within 0.25 Å of either boundary.

数学对照的预设参数为 De = 0.2 Hartree、a = 1.8 Å^−1、re = 1 Å，不作为 H2 测量或电子结构计算结果。主网格研究在 [0, 8] Å 使用 400、800、1200 个区间；独立区间研究在固定 1/150 Å 步长下，将右边界设为 2、4、6、8 Å。重复的 1200 区间、8 Å 计算复用且只计数一次。保存前三能级误差、实测收敛阶、Richardson 外推误差，以及全部共同束缚态误差和距边界 0.25 Å 内的概率。

A negative boundary control keeps the same mathematical depth and range but moves re to 0.12 Å. The hard radial wall then changes the spectrum materially even when the matrix eigenproblem is solved accurately. Short-domain loss of high states and near-dissociation domain sensitivity are retained as results. An eigensolver residual alone cannot establish physical or continuum convergence.

边界负对照保持同一数学势阱深度和范围参数，将 re 移至 0.12 Å。即使矩阵本征问题求解准确，径向硬壁仍会显著改变能谱。短区间丢失高能态和近解离态的区间敏感性作为结果保留。很小的矩阵本征残差本身不能证明物理模型正确或连续问题收敛。

## FCI-derived parameterization / FCI 数据参数化

For each electronic basis, the equilibrium optimization, actual isolated H energy and local energy curvature define De = 2EH−Emin, re = re,FCI and a = sqrt(kFCI/(2De)). No nonlinear fit to the full FCI curve is performed. Every available FCI curve point is nevertheless compared with Emin + VMorse(r), so the shape error remains visible. Electronic basis error, the Morse representation error, the finite-difference grid error and the finite-domain boundary error are different quantities.

对每种电子基组，利用平衡几何优化、实际孤立 H 原子能量与局部能量曲率确定 De = 2EH−Emin、re = re,FCI、a = sqrt(kFCI/(2De))。不对完整 FCI 曲线执行非线性拟合，但逐点比较 Emin + VMorse(r) 与全部已有 FCI 能量，显式保留势能形状误差。电子基组误差、Morse 表示误差、有限差分网格误差和有限区间边界误差彼此不同。

## Bound-only thermodynamic quantities / 仅束缚态的热力学量

For each model/isotope, temperatures 100, 298.15, 600, 1000, 2000 and 4000 K are used. The finite sum qbound = sum_v exp(−Ev/kBT) gives Fvib,bound = −kBT ln qbound, Uvib = sum_v pvEv and Svib = (Uvib−Fvib,bound)/T. A second partition sum is reported with energies measured from v = 0. A harmonic oscillator comparison uses its unbounded ladder analytically. The bound-only anharmonic partition omits continuum/dissociation, rotation, translation, nuclear-spin statistics, solvent and standard-state terms. It cannot be called a complete chemical Gibbs energy or used as a solution reaction free energy. High-temperature values are conditional on retaining only bound J = 0 states, even if the occupation of the highest retained state is small.

对每组模型和同位素，在 100、298.15、600、1000、2000、4000 K 计算有限束缚态配分和 qbound = sum_v exp(−Ev/kBT)，进而得到 Fvib,bound = −kBT ln qbound、Uvib = sum_v pvEv 和 Svib = (Uvib−Fvib,bound)/T；另保存以 v = 0 为能量零点的配分和。谐振比较采用无限能级梯的解析式。非谐束缚态配分和排除连续谱/解离、转动、平动、核自旋统计、溶剂与标准态项，不能称为完整化学 Gibbs 自由能，也不能用作溶液反应自由能。即使最高保留态占据很小，高温结果仍以仅纳入束缚 J = 0 态为条件。

## Executed results / 已执行结果

The two FCI-derived Morse models retain substantial shape errors against the actual electronic curve. Their 25-point RMSE is 0.00757282 Hartree for STO-3G and 0.00420654 Hartree for cc-pVDZ; maximum absolute errors are 0.0150622 and 0.00779212 Hartree. The corresponding Morse range parameters are 2.04277699 and 1.99016299 Å^−1. These errors are not statistical uncertainty bars and are not nuclear discretization errors. The spectroscopy below belongs to those approximate Morse models.

两种 FCI 参数化 Morse 模型与实际电子能曲线仍存在明显形状误差：STO-3G、cc-pVDZ 的各 25 点 RMSE 分别为 0.00757282、0.00420654 Hartree，最大绝对误差分别为 0.0150622、0.00779212 Hartree；Morse 范围参数分别为 2.04277699、1.99016299 Å^−1。这些量不是统计不确定度，也不是核本征求解的离散化误差。下表能谱属于这些近似 Morse 模型。

| FCI-derived Morse model / 模型 | Isotope / 同位素 | Harmonic wavenumber / 谐振波数 (cm^−1) | Numerical gap 0→1 / 数值能级差 (cm^−1) | Full-line Morse gap / 全实轴解析差 (cm^−1) | Numerical ZPE / 数值零点能 (Hartree) | Numerical / full-line states 数值/解析态数 |
|---|---|---:|---:|---:|---:|---:|
| STO-3G | H2 | 5003.208 | 4722.097 | 4723.858 | 0.011236766 | 18 / 18 |
| STO-3G | D2 | 3538.681 | 3397.084 | 3398.936 | 0.007979812 | 25 / 25 |
| cc-pVDZ | H2 | 4383.737 | 4117.259 | 4118.591 | 0.009834132 | 16 / 17 |
| cc-pVDZ | D2 | 3100.540 | 2966.494 | 2967.901 | 0.006986240 | 23 / 23 |

These numerical spectra use 1200 intervals over [0, 8] Å. The harmonic H2/D2 frequency ratio is 1.41386262 for every model, as required by the chosen nuclear masses. The anharmonic numerical ratios are 1.39004412 and 1.38792116 for the STO-3G- and cc-pVDZ-derived models; exact full-line Morse ratios are 1.38980483 and 1.38771184. The isotope ratio therefore does not retain the harmonic square-root law after anharmonicity is included. These are J = 0 level gaps without an intensity or spectroscopic selection-rule calculation; they are not assigned experimental IR lines.

表中数值谱采用 [0, 8] Å 上的 1200 个区间。按所用核质量，所有模型的谐振 H2/D2 波数比均为 1.41386262；STO-3G、cc-pVDZ 参数化模型的非谐数值能级差比分别为 1.39004412、1.38792116，全实轴解析比分别为 1.38980483、1.38771184。加入非谐性后，能级差比不再严格遵循谐振质量平方根律。这里仅计算 J = 0 能级差，未计算强度和光谱选择定则，不作为已指认的实验红外谱线。

The 800→1200 grid pairs give observed ZPE convergence orders 2.00057–2.00093. At 1200 intervals, the maximum error among the first three levels ranges from 3.975 to 5.705 cm^−1 across six spectra; the two-grid Richardson estimates reduce those first-three-level discrepancies to 0.00491–0.01097 cm^−1. These extrapolations do not establish all high-state or chemical accuracy. The main maximum matrix residual is 1.43×10^−14 Hartree and maximum quadrature orthogonality error is 1.18×10^−14, demonstrating why a tiny matrix residual must not be confused with grid error. Moving the model minimum to 0.12 Å increases H2 and D2 ground energies by 0.00369867 and 0.00186043 Hartree relative to the full-line formula; it also removes three bound states in each case.

800→1200 网格对的零点能实测收敛阶为 2.00057–2.00093。六组谱在 1200 个区间下，前三能级中的最大误差为 3.975–5.705 cm^−1；双网格 Richardson 外推后的前三能级差异缩小至 0.00491–0.01097 cm^−1。这不代表所有高能态或化学精度均已验证。主研究最大矩阵本征残差为 1.43×10^−14 Hartree，积分正交误差最大为 1.18×10^−14，显示矩阵残差与网格误差必须分开。负对照将平衡点移到 0.12 Å 后，H2、D2 基态较全实轴解析值分别升高 0.00369867、0.00186043 Hartree，并各少三个束缚态。

The main cc-pVDZ-derived H2 model misses its seventeenth analytical state. Its exact binding energy is only 6.69637×10^−7 Hartree, with asymptotic amplitude decay length 15.0913 Å. A separate seven-eigenproblem follow-up preserves the original main data. It recovers 17 states at right boundary 24 Å but the highest-state binding is then wrong by 104% at the original step size. At 96 Å and 115200 intervals (step 1/1200 Å), all 17 states remain, the maximum level discrepancy is 0.29645 cm^−1 and the shallow-state binding error is still 1.564%. The 96 Å calculation at half that resolution has 6.337% binding error. This is measured improvement with a remaining precision limit, not an assertion of arbitrary convergence. The main 8 Å partition table remains unchanged and must be read as a finite-domain bound-state truncation.

主研究中，cc-pVDZ 参数化 H2 模型漏掉第 17 个解析束缚态。该态解析束缚能仅为 6.69637×10^−7 Hartree，渐近振幅衰减长度为 15.0913 Å。独立执行七个追加本征问题，保留原主数据不变。右边界扩大到 24 Å 后恢复 17 态，但在原步长下最高态束缚能仍错约 104%。使用 96 Å、115200 个区间，即 1/1200 Å 步长时，保留全部 17 态，能级最大差异为 0.29645 cm^−1，最浅态束缚能相对误差仍有 1.564%；同一区间半分辨率计算的误差为 6.337%。这是具有剩余精度限制的收敛改善，不声称任意精度收敛。主 8 Å 配分表保持不变，应按有限区间束缚态截断解释。

At 4000 K the cc-pVDZ-derived H2 main bound-only F is 0.00643903 Hartree, 0.000615759 Hartree below its harmonic value. Its difference from the full-line Morse bound partition is −4.09457×10^−6 Hartree. The finest long-domain follow-up reduces the latter discrepancy to −6.46456×10^−8 Hartree, while still omitting the physical continuum. A small partition difference does not certify the missing high state or the electronic potential shape.

4000 K 时，cc-pVDZ 参数化 H2 主模型仅束缚态的 F 为 0.00643903 Hartree，比谐振值低 0.000615759 Hartree；相对全实轴 Morse 束缚态配分结果的差为 −4.09457×10^−6 Hartree。最细长区间追加计算将后者缩小到 −6.46456×10^−8 Hartree，但仍排除物理连续谱。配分量差异较小，不等于遗漏高能态或电子势能形状已经验证。

| Workload / 计算工作量 | Pilot | Main / 主研究 | Tail follow-up / 长尾追加 |
|---|---:|---:|---:|
| Tridiagonal eigenproblems / 三对角本征问题 | 6 | 38 | 7 |
| Summed interior nodes / 累计内部节点 | 3994 | 27562 | 228593 |
| Summed returned bound states / 累计返回束缚态 | 139 | 783 | 118 |
| Main-spectrum level rows / 主谱能级行数 | 48 | 130 | 118 follow-up rows / 追加行 |
| Isotope/model spectra / 同位素模型谱 | 2 | 6 | 1 model, 7 numerical settings / 1 模型、7 设置 |
| Temperature cases / 温度案例 | 2 | 36 | 7 checks at 4000 K / 7 次 4000 K 检查 |
| Actual FCI curve residual rows / 实际 FCI 残差行 | 0 | 50 | 0 |
| New electronic-structure jobs / 新增电子结构作业 | 0 | 0 | 0 |

The 18 unit tests pass; a saved-data audit passes 107 checks covering input/output hashes, historical execution snapshots, row counts, wavefunction quadrature, energy differences, D0 accounting and partition recomputation. These checks establish the recorded numerical workflow, not experimental spectroscopy. Runtime versions are Python 3.12.14, NumPy 2.4.6 and SciPy 1.18.0. Pilot and main snapshots predate the separate tail-follow-up CLI addition; each summary identifies its own exact historical source hash. The current source matches the tail follow-up snapshot.

18 项单元测试全部通过；保存数据的只读审计通过 107 项检查，涵盖输入/输出哈希、历史执行快照、行数、波函数积分正交性、能级差、D0 能量核算与配分量复算。这些检查验证记录中的数值流程，不代表实验光谱验证。运行版本为 Python 3.12.14、NumPy 2.4.6、SciPy 1.18.0。pilot 与主研究快照早于独立长尾 CLI 入口的添加，各 summary 均引用自己的确切历史源码哈希；当前源码与长尾追加快照一致。

## Reproduction and scope / 复现与范围

Run from the repository root with the existing NumPy/SciPy chemistry environment:

    python quantumequi/scripts/vibration_extension.py --pilot
    python quantumequi/scripts/vibration_extension.py
    python quantumequi/scripts/vibration_extension.py --tail-followup
    python -m unittest discover -s tests -p test_quantumequi_vibration_extension.py -v

The pilot has no FCI input and remains in its own directory with the exact executed source snapshot. The main calculation requires the separately generated correlation curve and equilibrium reference. Source snapshots, input/output SHA-256, software versions and counts are recorded. This extension creates no new electronic-structure jobs. Eigenproblems, interior grid nodes, returned states, spectra and partition cases are separate workload units.

pilot 不读取 FCI 输入，保存在独立目录并保留当次实际执行源码快照。主计算读取独立生成的电子相关能曲线与平衡点参考。程序记录源码快照、输入/输出 SHA-256、软件版本和计算数量。本扩展不新增电子结构作业；矩阵本征问题、内部网格点、返回态数、能谱数与温度配分案例分别计数。
