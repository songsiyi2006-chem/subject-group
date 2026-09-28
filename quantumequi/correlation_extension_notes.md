# H₂ correlation and basis extension: methods and evidence

## English

### Scope and completed calculations

This extension performs real, clamped-nuclei H₂ electronic-structure calculations with Psi4 1.11, comparing RHF, UHF, and FCI in STO-3G and cc-pVDZ. FCI is a reference for this two-electron molecule in a specified finite orbital basis. It does not validate the supplied copper-containing model, an electrode, or a complete-basis experimental prediction. Earlier frozen artifacts remain unchanged.

The main grid contains 25 distances in Å: 0.50, 0.54, 0.58, 0.62, 0.66, 0.70, 0.74, 0.78, 0.82, 0.94, 1.10, 1.26, 1.30, 1.42, 1.58, 1.74, 1.80, 1.90, 2.10, 2.30, 2.50, 2.70, 3.00, 3.50, and 4.00. It includes all eight previous RHF test distances and five extrapolation distances. Thus the parent module can compare learned-model errors with RHF–FCI method differences at identical geometries. These FCI values are not new neural-training observations.

| Stage | Configuration attempts | Dispatched energy drivers | Converged drivers | Rejected before driver |
|---|---:|---:|---:|---:|
| Initial three-distance pilot | 20 | 14 | 14 | 6 |
| FCI-only pilot recovery | 6 | 6 | 6 | 0 |
| Main scan, atoms, minimum and curvature | 183 | 183 | 183 | 0 |
| Total | 209 | 203 | 203 | 6 |

The main 183 calls comprise 150 scan energies, two isolated H energies, 21 FCI minimum-search energies, and ten FCI curvature energies. Across all stages there are 87 FCI drivers and 116 stand-alone SCF drivers. Each FCI driver includes an underlying RHF calculation. The logs contain 203 SCF convergence messages; these are not additional explicit driver calls. There are zero gradient-driver and molecular Hessian-driver calls. The bond curvature comes from scalar energy differences. The actual 203 calls remain below the implementation's 206-call bound and the parent task's 250-call budget.

### Methods and spin diagnostics

H₂ has charge zero, input multiplicity one, and nuclei at z=±R/2 in Å. Calculations use C1 symmetry, the core-Hamiltonian guess, PK SCF integrals, energy and density convergence thresholds of 1e−12, a 200-iteration SCF cap, one CPU thread, and 500,000,000 bytes of requested memory. Isolated H uses charge zero, multiplicity two, and UHF; for a single electron this gives the finite-basis atomic ground-state energy.

For H₂ UHF, GUESS_MIX=True allows different α/β spatial orbitals, and STABILITY_ANALYSIS=FOLLOW requests an internal stability search and following of an unstable direction. RHF uses CHECK. These choices follow the primary [Psi4 guess-mixing documentation](https://psicode.org/psi4manual/master/autodir_options_c/scf__guess_mix.html) and [SCF/stability documentation](https://psi4.github.io/psi4docs/master/scf.html). Convergence alone does not establish symmetry breaking. We independently evaluate

\[
\langle S^2\rangle=S_z(S_z+1)+N_\beta-
\sum_{i\in\alpha,j\in\beta}|\mathbf c_i^\mathsf T S\mathbf c_j|^2,\qquad
S_z=(N_\alpha-N_\beta)/2.
\]

The overlap matrix and occupied α/β coefficients are saved. This identity applies to a single determinant. For FCI entries those matrices describe the underlying RHF reference, not the correlated CI state, and the reference spin value is never substituted for a correlated-state measurement.

A sampled point is called a lower spin-broken UHF branch only when S²>1e−4 and E_UHF<E_RHF−1e−7 Ha. There are 28 such main points and 22 restricted-like points. In both bases the first sampled lower branch occurs at 1.26 Å, while 1.10 Å is restricted-like; this does not locate an exact bifurcation. All 50 main UHF points have positive recorded lowest internal stability eigenvalues. Their grid minima are 0.04035695 and 0.04062163 in the program's orbital-Hessian convention. Local internal stability is not proof of a global Hartree–Fock minimum or spin purity.

FCI uses DETCI, RHF reference orbitals, zero frozen occupied/virtual orbitals, one root, CI_MAXITER=100, energy threshold 1e−12, residual threshold 1e−10, and one CI thread. See the primary [DETCI manual](https://psicode.org/psi4manual/master/detci.html) and [DETCI options](https://psicode.org/psi4manual/master/autodir_options_c/module__detci.html). H₂ has two or ten spatial orbitals; the complete Nα=Nβ=1 sector contains four or 100 determinants. “Full” refers only to this finite orbital space.

### Retained failures and metadata repair

Although master documentation lists optional DETCI S_SQUARED, installed Psi4 1.11 rejects it. Six pilot FCI configurations therefore failed before energy-driver dispatch; the other 14 pilot drivers succeeded. Only the six FCI calls were recovered after removing that diagnostic. Their separate pilot_fci_recovery directory and all three exact executed-code snapshots are retained. FCI S² remains null because no correlated spin expectation was returned.

The pilot also exposed invalid pointer-like integers from CIWavefunction.ndet(). These numbers never affected energies or settings. Main-run determinant counts are parsed from the numerical logs. For the recovery pilot, metadata_corrections.json preserves the original API values, original/corrected ledger hashes, and the reason for replacing them with the four/100 counts printed by DETCI. This was a metadata-only correction with no repeated quantum calculation. Each summary distinguishes the executed-source hash from the finalizer-source hash.

### Results and physical limits

All 50 same-geometry, same-basis comparisons satisfy E_FCI≤E_UHF≤E_RHF within 1e−8 Ha. Contracted STO-3G and cc-pVDZ are not nested; no general cross-basis variational theorem is asserted.

| Quantity | STO-3G | cc-pVDZ |
|---|---:|---:|
| RHF total energy at 4 Å / Ha | −0.614869973999 | −0.782198208378 |
| UHF total energy at 4 Å / Ha | −0.933166094408 | −0.998569700894 |
| FCI total energy at 4 Å / Ha | −0.933171361844 | −0.998606186141 |
| UHF ⟨S²⟩ at 4 Å | 0.999980059084 | 0.999764679904 |
| RHF−FCI at 4 Å / Ha | 0.318301387844 | 0.216407977763 |
| FCI minimum bond length / Å | 0.734865227479 | 0.760893444967 |
| FCI minimum total energy / Ha | −1.137306051222 | −1.163672981196 |
| Two isolated H energies / Ha | −0.933163699115 | −0.998556806839 |
| Electronic well depth D_e / Ha | 0.204142352108 | 0.165116174357 |
| Bond curvature / Ha Å⁻² | 1.703746695358 | 1.307967353975 |
| Difference of two central-curvature steps / Ha Å⁻² | 0.000102242570 | 0.000078580982 |

Long-bond UHF recovers much of the RHF error, but its S² near one differs substantially from a pure singlet's zero. At 4 Å, FCI differs from twice the independently computed H energy by −0.000007662729 and −0.000049379302 Ha. Therefore the endpoint is not treated as infinity: E_infinity=2E_H defines the separated-atom reference. No counterpoise correction was performed. These energies contain neither zero-point/thermal corrections nor solvent, electrode, or relativistic effects.

Bounded scalar minimization was fixed in advance on [0.60,0.90] Å with xatol=1e−7 Å and at most 12 calls per basis; it converged in ten and 11 calls. Five further energies at offsets −2h,−h,0,h,2h, with h=0.0025 Å, give

\[
k={-E_{+2}+16E_{+1}-30E_0+16E_{-1}-E_{-2}\over12h^2}.
\]

The five-point slopes are 1.45e−9 and 5.78e−8 Ha Å⁻¹. Ordinary central curvatures at h and 2h differ by the final table row. That difference is a step-sensitivity diagnostic, not a confidence interval. The equilibrium energy, curvature, and actual isolated-atom reference feed the nuclear-motion agent. A Morse curve built from them remains an approximation whose residuals to the actual FCI curve must be reported.

### Reproduction, files and checks

The executable is [correlation_extension.py](scripts/correlation_extension.py). It refuses to overwrite executed results. The current executable omits the rejected optional parameter and refuses FCI recovery when no pre-driver FCI rejection exists; each historical stage preserves its actual executed_code.py.txt. Exact historical replay uses those archived bytes with a new output location. The --finalize-only mode processes saved data without importing Psi4. A fresh successful pilot needs no recovery directory: complete-study accounting requires the main and pilot ledgers, and includes the historical recovery stage only if that directory exists. An existing recovery directory with a missing ledger is an error, rather than a silently omitted stage. Two lightweight temporary-directory tests cover both stage layouts and missing-ledger rejection; they execute no quantum driver.

Interfaces are [curve.csv](results/extensions/correlation/curve.csv), [equilibrium_reference.json](results/extensions/correlation/equilibrium_reference.json), [energy_driver_jobs.json](results/extensions/correlation/energy_driver_jobs.json), [summary.json](results/extensions/correlation/summary.json), and reference_orbitals.json. The curve table records basis, method, point_id, R_A, status, energy_Hartree, energy_total_Hartree, energy_electronic_Hartree, nuclear_repulsion_Hartree, S2, and log/job identifiers. The first energy name is an exact alias of total energy; electronic energy excludes nuclear repulsion. The equilibrium JSON is a dictionary with a references list and explicit units.

Runtime versions: Psi4 1.11, Python 3.12.14, NumPy 2.5.2, SciPy 1.18.0. Measured execution times were 5.219 s (pilot), 2.695 s (recovery), and 62.124 s (main), excluding launcher overhead. Deterministic settings/search need no random seeds. Public logs remove local absolute paths and hostnames, preserving numerical text; their original decoded-text and published-byte SHA-256 hashes are saved. Inputs are relative to quantumequi and outputs to each stage directory. Scratch files are excluded.

Twenty-one focused tests pass: analytic spin and finite-difference controls, invalid inputs, branch logic, exact grids, hashes, accounting, log energies, SCF counts, total/electronic energy identity, determinant counts, spin recomputation, large-R branches, atomic references, minimum curvature, and the CSV interface. An independent agent additionally checked all 183 main log energies (maximum discrepancy 4.88e−15 Ha), MO orthonormality (3.55e−15), and independently reconstructed reference S² (4.44e−16). No extra quantum jobs were run for these reviews. Numerical consistency does not remove finite-basis physical limitations.

## 中文

### 范围与实际计算

本扩展通过已安装的 Psi4 1.11 开展真实的固定核坐标 H₂ 计算，在 STO-3G 与 cc-pVDZ 基组内比较 RHF、UHF 和 FCI。FCI 是这一双电子小分子在指定有限轨道空间中的参考，并不意味着原始含铜模型、电极体系或完备基组实验预测得到验证。此前冻结的产物保持不变。

主网格预先固定 25 个距离，单位 Å：0.50、0.54、0.58、0.62、0.66、0.70、0.74、0.78、0.82、0.94、1.10、1.26、1.30、1.42、1.58、1.74、1.80、1.90、2.10、2.30、2.50、2.70、3.00、3.50、4.00。其中覆盖此前全部八个 RHF 测试距离及五个外推距离，使父模块能够在相同几何上比较学习误差与 RHF–FCI 方法偏差。这些 FCI 能量不是新增神经网络训练标签。

| 阶段 | 配置尝试 | 已发起能量驱动调用 | 收敛调用 | 驱动前拒绝 |
|---|---:|---:|---:|---:|
| 最初三距离先导 | 20 | 14 | 14 | 6 |
| 仅 FCI 的先导补算 | 6 | 6 | 6 | 0 |
| 主扫描、原子、极小值及曲率 | 183 | 183 | 183 | 0 |
| 合计 | 209 | 203 | 203 | 6 |

主研究 183 次调用由 150 个扫描能量、两个独立 H 原子能量、21 个 FCI 极小值搜索能量及十个 FCI 曲率能量组成。全部阶段共有 87 次 FCI 驱动与 116 次独立 SCF 驱动；每次 FCI 内部还执行 RHF。日志共出现 203 条 SCF 收敛记录，未重复计为额外显式驱动调用。梯度驱动与分子 Hessian 驱动调用均为零；键曲率使用标量能量差分。实际 203 次调用低于实现预设的 206 次上限及父任务的 250 次预算。

### 方法与自旋诊断

H₂ 电荷零、输入多重度一，两个核位于 z=±R/2，距离使用 Å。采用 C1 对称性、核 Hamiltonian 初猜、PK SCF 积分、1e−12 能量与密度收敛阈值、最多 200 次 SCF 迭代、一个 CPU 线程及 500,000,000 字节请求内存。独立 H 原子采用电荷零、多重度二和 UHF；对单电子体系，这给出该有限轨道基组内的原子基态能量。

H₂ UHF 使用 GUESS_MIX=True，允许 α/β 空间轨道不同；STABILITY_ANALYSIS=FOLLOW 搜索内部不稳定性，并尝试沿不稳定方向优化。RHF 使用 CHECK。设置依据为第一手 [Psi4 初猜文档](https://psicode.org/psi4manual/master/autodir_options_c/scf__guess_mix.html)和 [SCF/稳定性文档](https://psi4.github.io/psi4docs/master/scf.html)。单凭收敛不能证明对称性破缺，因此独立计算

\[
\langle S^2\rangle=S_z(S_z+1)+N_\beta-
\sum_{i\in\alpha,j\in\beta}|\mathbf c_i^\mathsf T S\mathbf c_j|^2,\qquad
S_z=(N_\alpha-N_\beta)/2.
\]

保存重叠矩阵及 α/β 占据轨道系数以支持独立核查。这个恒等式适用于单个行列式。FCI 条目的这些矩阵描述底层 RHF 参考，不是相关 CI 态；参考态的自旋值不能代替相关态的自旋测量。

只有同时满足 S²>1e−4 及 E_UHF<E_RHF−1e−7 Ha，才标记为采样到的较低能量破缺分支。主研究共有 28 个此类点及 22 个近似受限解。两个基组首次采样到该分支的距离均为 1.26 Å，1.10 Å 仍为近似受限解；这不等于精确分岔位置。全部 50 个 UHF 点记录的最低内部稳定性特征值均为正，两个网格最小值分别为 0.04035695 和 0.04062163，遵循软件的轨道 Hessian 约定。局部内部稳定不证明全局最优 Hartree–Fock 解或自旋纯态。

FCI 使用 DETCI、RHF 参考、零冻结占据/虚轨道、一个目标根、CI_MAXITER=100、1e−12 能量阈值、1e−10 残差阈值及一个 CI 线程。第一手来源为 [DETCI 手册](https://psicode.org/psi4manual/master/detci.html)和 [DETCI 参数文档](https://psicode.org/psi4manual/master/autodir_options_c/module__detci.html)。H₂ 分别具有两个和十个空间轨道；完整 Nα=Nβ=1 空间包含四个和 100 个行列式。“全”仅指这个有限轨道空间。

### 保留的失败与元数据修复

主版本文档虽然列出可选 DETCI S_SQUARED，本机 Psi4 1.11 却拒绝该参数。最初六个 FCI 先导配置因此在发起能量驱动之前失败，另 14 次调用全部成功。移除该诊断后，只补算六个 FCI 能量，保存在独立 pilot_fci_recovery 目录；三个阶段的精确执行源码均保留。未获得相关态自旋期望值，故 FCI S² 保持为空。

先导还发现 CIWavefunction.ndet() 返回类似指针的异常整数；这些值未参与能量计算或配置。主研究从数值日志读取行列式数。补算先导的 metadata_corrections.json 保留原始 API 值、修复前后账本哈希，以及依据 DETCI 日志改为四个/100 个的理由。这是元数据修正，没有重跑量子作业。汇总明确区分量子作业执行源码哈希与后来整理程序的源码哈希。

### 结果与物理边界

同几何、同基组的 50 组主结果均在 1e−8 Ha 容差内满足 E_FCI≤E_UHF≤E_RHF。收缩 STO-3G 与 cc-pVDZ 并非嵌套空间，不能据此宣称普遍的跨基组变分定理。

| 量 | STO-3G | cc-pVDZ |
|---|---:|---:|
| 4 Å 的 RHF 总能量 / Ha | −0.614869973999 | −0.782198208378 |
| 4 Å 的 UHF 总能量 / Ha | −0.933166094408 | −0.998569700894 |
| 4 Å 的 FCI 总能量 / Ha | −0.933171361844 | −0.998606186141 |
| 4 Å 的 UHF ⟨S²⟩ | 0.999980059084 | 0.999764679904 |
| 4 Å 的 RHF−FCI / Ha | 0.318301387844 | 0.216407977763 |
| FCI 极小值键长 / Å | 0.734865227479 | 0.760893444967 |
| FCI 极小值总能量 / Ha | −1.137306051222 | −1.163672981196 |
| 两个独立 H 原子能量 / Ha | −0.933163699115 | −0.998556806839 |
| 电子势阱深度 D_e / Ha | 0.204142352108 | 0.165116174357 |
| 键伸缩曲率 / Ha Å⁻² | 1.703746695358 | 1.307967353975 |
| 两种中心曲率步长之差 / Ha Å⁻² | 0.000102242570 | 0.000078580982 |

长键 UHF 恢复了 RHF 的大部分误差，但接近一的 S² 明显偏离纯单重态的零。4 Å 的 FCI 能量与两倍独立 H 能量之差分别为 −0.000007662729 和 −0.000049379302 Ha，故不将该端点当作无限远；分离原子参考采用 E_infinity=2E_H。未进行 counterpoise 校正，能量也不包含零点、热、溶剂、电极或相对论修正。

极小值搜索预先固定在 [0.60,0.90] Å，xatol=1e−7 Å，每基组最多 12 次调用；实际十次和 11 次收敛。再取偏移 −2h、−h、0、h、2h 的五个能量，h=0.0025 Å，通过

\[
k={-E_{+2}+16E_{+1}-30E_0+16E_{-1}-E_{-2}\over12h^2}
\]

计算曲率。五点一阶导数分别为 1.45e−9 和 5.78e−8 Ha Å⁻¹。采用 h 与 2h 的普通中心差分曲率之差见表末行；这是步长敏感性，不是统计置信区间。极小值能量、曲率和独立原子参考交给核运动代理；由它们构建的 Morse 曲线仍是近似，必须报告相对真实 FCI 曲线的残差。

### 复现、接口和验证

程序为 [correlation_extension.py](scripts/correlation_extension.py)，会拒绝覆盖已有执行产物。当前源码移除了被拒绝的可选参数，且在不存在驱动前 FCI 拒绝时拒绝重复补算；各历史阶段保存实际执行的 executed_code.py.txt，可在新的输出位置精确回放。--finalize-only 只处理保存数据，不导入 Psi4。新先导成功时不需要补算目录：完整研究账本始终要求主研究与先导账本，仅在历史补算目录存在时才纳入补算阶段；补算目录存在但账本缺失仍报错，不会默默漏计。两项轻量临时目录测试覆盖两种阶段布局及缺失账本拒绝，不执行量子驱动。

接口包括 [curve.csv](results/extensions/correlation/curve.csv)、[equilibrium_reference.json](results/extensions/correlation/equilibrium_reference.json)、[energy_driver_jobs.json](results/extensions/correlation/energy_driver_jobs.json)、[summary.json](results/extensions/correlation/summary.json) 及 reference_orbitals.json。曲线表记录 basis、method、point_id、R_A、status、energy_Hartree、energy_total_Hartree、energy_electronic_Hartree、nuclear_repulsion_Hartree、S2 及调用/日志标识。第一个能量字段是总能量的精确别名；电子能排除核排斥。平衡 JSON 顶层字典含 references 列表，单位明确。

版本为 Psi4 1.11、Python 3.12.14、NumPy 2.5.2、SciPy 1.18.0。先导、补算和主研究实测运行耗时分别为 5.219 s、2.695 s、62.124 s，不含启动器开销。确定性设置和搜索不需要随机种子。公开日志删除本地绝对路径和主机名，保留数值文本；同时保存原始解码文本及公开字节 SHA-256。输入路径相对 quantumequi，输出路径相对对应阶段目录，scratch 文件不纳入公开证据。

21 项针对性测试全部通过，覆盖解析自旋/差分对照、非法输入、分支逻辑、网格、哈希、账本、日志能量、SCF 计数、总能/电子能恒等式、行列式数、自旋重算、长键分支、原子参考、极小值曲率与 CSV 接口。独立代理还核查 183 份主日志能量，最大差 4.88e−15 Ha；MO 正交误差不超过 3.55e−15，独立重建参考 S² 的差不超过 4.44e−16。复核未新增量子作业。数值一致性不消除有限基组的物理局限。
