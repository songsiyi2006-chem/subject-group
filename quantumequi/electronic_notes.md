# Electronic-structure and surrogate evidence / 电子结构与替代模型证据

## English methods and results

### Scope and source reconstruction

This module separates three claims: diagnostics of the supplied extended-Hückel-like matrix construction; executed H2 electronic-structure calculations; and supervised interpolation of those H2 calculations. None supplies a Cu–N4 catalytic reaction energy or trains the supplied EGNN. The original eight-atom input is the inventory C2HCuN4; the original code and capture remain unchanged. The source class is extracted through its Python AST and called only to assemble its matrices, without rerunning the neural potential, NEB, or Hessian workflow.

The source constructs 31 nominal valence orbitals and counts 40 valence electrons, hence 20 doubly occupied orbitals. Its overlap expression depends on atom positions and an element-dependent radial exponent but has no angular orbital dependence. Distinct orbitals on the same atom therefore have overlap exactly one and identical rows. The Cu six-orbital block contributes 15 such pairs and the six C/N four-orbital blocks contribute 36: all 51 pairs have unit overlap. Both supplied geometries have overlap rank 8 and 23 null directions. The minimum numerical eigenvalues, −1.40e−15 and −8.08e−16, are roundoff-scale quantities; they do not remove this exact repeated-row defect.

For the source transformation, small eigenvalues are replaced by a floor η:

\[
S_\eta=U\,\mathrm{diag}(\max(s_i,\eta))U^T,\quad
X_\eta=S_\eta^{-1/2},\quad
X_\eta H X_\eta v=\epsilon v,\quad C=X_\eta v.
\]

This is an eigenproblem for the modified metric Sη. The audited quantities include both HC−SCε and HC−SηCε, and CᵀSC−I. A stable result for the modified metric does not solve the original singular problem.

| Reactant clipping floor | Source HOMO / eV | 2 sum occupied eigenvalues / eV | max abs original residual |
|---|---:|---:|---:|
| 1e−3 | 10050 | 204251.290101 | 331.898660 |
| 1e−4 | 100500 | 2045988.734345 | 1049.465152 |
| 1e−5 | 1005000 | 20463363.728759 | 3318.671454 |
| 1e−6 | 10050000 | 204637113.728200 | 10494.528399 |
| 1e−7 | 100500000 | 2046374613.728143 | 27682.126454 |
| 1e−8 | 1005000000 | 20463749613.728127 | 87538.654683 |

At every floor the maximum original-metric orthonormality error is 1. The default η=1e−5 reconstructs the captured rounded source output, including its 471896257.61 kcal/mol occupied-orbital sum. This is a reproducible numerical pathology, not a large physical binding energy. The Hamiltonian parameters and this one-electron orbital sum are also distinct from an ab initio total energy with a consistent electron-interaction and nuclear-repulsion treatment.

### Reviewed generalized eigensolver

`generalized_eigh` validates finite equal-size symmetric H and S. The strict mode rejects an indefinite, singular, or insufficiently positive metric. Indefiniteness is rejected when the negative eigenvalue exceeds the explicit roundoff tolerance 64 εmachine n max(1,|s|max). The default retained-space cutoff is max(1e−12,1e−10 smax). No angular overlap model is invented.

Canonical mode deliberately retains only the positive metric subspace:

\[
X=U_r s_r^{-1/2},\quad X^THXv=\epsilon v,\quad C=Xv.
\]

It reports projected and full-space residuals separately, as well as discarded-space Hamiltonian coupling. This is a Ritz problem in an explicitly restricted space. The distinction between positive-definite generalized eigenproblems and rank reduction is consistent with the [SciPy `eigh` specification](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html) and the orthogonalization discussion in the [Psi4 SCF manual](https://psicode.org/psi4manual/master/scf.html). These references support the methods; the saved local matrices support the following values.

| Geometry | Retained rank | Discarded directions | max abs projected residual | max abs full residual | max abs CᵀSC−I |
|---|---:|---:|---:|---:|---:|
| Reactant | 8 | 23 | 1.823481e−13 | 4.593825 | 7.018553e−15 |
| Product | 8 | 23 | 4.623372e−13 | 4.626220 | 1.400616e−14 |

The small projected residual cannot establish that the original equation is solved. Moreover, eight retained orbitals can hold at most 16 closed-shell electrons, so the 40-electron occupancy is explicitly rejected. This rank reduction is a diagnostic control, not repaired physical EHT.

Six deterministic random SPD controls, seeds 901–906, agree with the independent SciPy generalized solver to 5.95e−14 or better in their eigenvalues. A compatible rank-deficient control has a full residual below 5e−16. An incompatible rank-deficient control has a projected residual below 4e−16 but a full residual 0.282843. Twenty unit tests cover these distinctions, congruence and permutation invariance, shape/finiteness/symmetry rejection, electron capacity and parity, source overlap degeneracy, fixed dataset partitions, and analytical radial derivatives. Unit tests do not invoke Psi4 or reproduce the full study.

### Actual quantum reference

The existing Psi4 1.11 environment executed conventional PK RHF/STO-3G energies and analytical gradients for neutral singlet H2, in C1 symmetry, with `no_com` and `no_reorient`. The two nuclei lie at z=±R/2 Å. Energy and density convergence thresholds were 1e−12; maximum SCF iterations were 100. The job requested one CPU thread and 500 MB decimal memory, printed as 476.837 MiB. There were no installations, DFT jobs, experimental measurements, or Cu-cluster ab initio jobs. NumPy was 2.5.2 and SciPy 1.18.0 in this quantum execution environment.

The main scan contains 33 points R=0.50+0.04i Å, i=0,…,32, and nine stretched points R=1.90+0.10i Å, i=0,…,8. The grid and partition rule were fixed before fitting: even grid indices are training (17); odd indices 1 modulo 4 are validation (8); odd indices 3 modulo 4 are test (8); the nine longer distances are OOD. These splits assess geometry interpolation/extrapolation for one molecule. They are not independent molecular or reaction datasets.

The returned energy is the Born–Oppenheimer total energy. The separately stored electronic energy subtracts the nuclear term. Psi4's bohr-to-angstrom constant is 0.52917721067. The internal consistency identities are

\[
E_{\rm total}=E_{\rm electronic}+E_{\rm nuclear},\quad
E_{\rm nuclear}=a_0/R,\quad
\frac{dE}{dR}=\frac{g_{2z}-g_{1z}}{2},\quad F_{2z}=-g_{2z}.
\]

Gradients returned in Hartree/bohr are divided by a0 to obtain Hartree/Å. The pilot at 0.74 Å gives total energy −1.1167593073781568 Hartree, electronic energy −1.8318636461214000 Hartree, and nuclear repulsion +0.7151043387432432 Hartree. Its occupied-orbital sum is −1.1571077193736687 Hartree, explicitly different from both electronic and total energies.

The lowest point on the main grid is R=0.70 Å, Etotal=−1.1173490349988597 Hartree; dE/dR=−0.0259704221 Hartree/Å is nonzero, so this is not a claimed optimized geometry. Maximum net-force error is 1.11e−15 Hartree/Å. Maximum Fock generalized-eigen residual is 4.44e−16 and MO metric-orthonormality error 1.11e−15. AO overlap, Fock, core Hamiltonian, alpha density, MO coefficients, orbital energies, and requested analytical gradients are saved for all 55 main jobs.

Independent central energy differences use R=0.70,1.14,1.78 Å and h=0.001,0.0005 Å. The largest gradient discrepancy is 2.220864e−6 Hartree/Å at h=0.001 Å; at that same anchor halving h reduces it to 5.552152e−7. No displaced points enter surrogate training. This follows the general comparison of analytical and finite-difference gradients in the [official Psi4 gradient example](https://github.com/psi4/psi4/blob/master/samples/fd-gradient/input.dat), with the numerical results arising from the recorded local H2 jobs.

The neutral H doublet UHF/STO-3G energy is −0.46658184955727544 Hartree, giving a separated-atom reference of −0.9331636991145509 Hartree. At R=2.70 Å the restricted H2 energy remains −0.6809407604974436 Hartree. This illustrates the limitations of the restricted minimal-basis curve; it is not a benchmark of accurate dissociation or electron correlation. No FCI, correlated method, basis-convergence calculation, or H2 thermochemistry was performed.

### Supervised distance surrogate and negative results

Two deterministic radial models use 17 training distances as Gaussian centers, plus constant and linear terms:

\[
\widehat E(R)=E_0+c_0+c_1R+\sum_j c_{j+1}\exp[-(R-R_j)^2/(2\ell^2)].
\]

One fits energies; one fits energies and dE/dR. Both return forces by differentiating the fitted energy. Training energy and gradient standard deviations scale the least-squares blocks; an explicit ridge penalty regularizes all coefficients. Each family uses four prespecified lengthscales (0.06,0.12,0.24,0.48 Å) and two ridges (1e−10,1e−6): 16 linear fits total. Selection minimizes the sum of validation energy and gradient RMSE divided by their respective **training** standard deviations. Validation gradients may therefore select the energy-only model, although they are not fit targets. Both selected models use ℓ=0.24 Å and ridge=1e−10. There is no refit on validation or test data, no neural training seed, and no test-based tuning.

Baselines use the same training energies: a constant mean with zero force, piecewise linear interpolation/extrapolation, and a natural cubic spline with its analytical derivative. The linear derivative is not defined at interpolation knots; recorded non-knot test metrics avoid that ambiguity.

| Model | Test energy RMSE / Hartree | Test gradient RMSE / Hartree Å−1 | OOD energy RMSE / Hartree | OOD gradient RMSE / Hartree Å−1 |
|---|---:|---:|---:|---:|
| RBF energy only | 8.624253e−6 | 1.156169e−4 | 0.175473 | 0.459351 |
| RBF energy and force | 5.561350e−6 | 1.637466e−4 | 0.118352 | 0.448587 |
| Natural cubic spline | 2.133506e−4 | 2.974637e−3 | 0.222780 | 0.766237 |
| Piecewise linear | 1.103484e−3 | 2.239298e−3 | 0.044450 | 0.108744 |
| Training mean, zero force | 0.096071 | 0.269899 | 0.268390 | 0.161621 |

Gradient fitting improves test energy error but gives worse test gradient error than the selected energy-only RBF. Both RBFs perform much worse than the simple linear baseline in the stretched OOD region; even the zero-force baseline has lower OOD gradient error. Those outcomes are retained. Low interpolation errors do not make either radial model a transferable ML interatomic potential, and an accurate fit to RHF cannot remove RHF's physical error.

### Accounting, reproduction and evidence files

The main quantum run completed 55 SCF jobs: 42 H2 scan jobs with analytical gradients, 12 displaced-energy jobs, and one isolated H job. The successful pilot adds three SCF jobs and one gradient, giving **58 SCF jobs and 43 analytical gradients overall**, with zero failed SCF jobs. Main quantum execution took 15.81 s. The matrix audit performs two source matrix assemblies, 12 clipping controls, two canonical controls, two strict-metric rejections, two capacity rejections, and eight independent solver controls. Learning performs 16 RBF least-squares fits and three baselines, yielding 210 prediction rows and 20 metric rows. Unit-test calculations are excluded from study counts.

The first successful pilot's exact `executed_script.py.txt` is preserved. Subsequent code changes affect host-identity redaction and provenance metadata only. The pilot's historical source hash is not replaced by the final source hash. Machine paths and host identity are removed from public Psi4 logs while numerical SCF output is retained; pre-redaction and public hashes are recorded. Known Psi4 timing files are confined to this module and archived as text. An initially incorrect test-discovery working directory found zero tests; that transcript is retained under `results/electronic/environment`, followed by the final 20-test passing record.

Main and pilot `summary.json` use `inputs_sha256` relative to the `quantumequi` module and `outputs_sha256` relative to their respective result directory. `executed_source_snapshot` resolves a historical pilot source mismatch. `source_audit.json` separately hashes the source capture. The surrogate summary explicitly declares its local input-hash base. `verify_saved_energy.py` independently checks all 58 stored jobs: density-matrix energy identities, alpha-electron traces, Fock residuals, printed SCF convergence, log hashes, and printed versus stored energies. All checks pass; the largest main density-energy/Fock discrepancy is 8.88e−16 Hartree. This is numerical verification, not external accuracy validation.

Reproduction uses the installed Psi4 Python for `python quantumequi/scripts/electronic_reviewed.py --stage all`; `--stage quantum --pilot` runs the smaller quantum pilot. Stages `audit`, `quantum`, `fit`, and `summarize` are available separately. `summarize` performs no new quantum calculation. Existing outputs should be archived before a deliberate rerun; differing code is rejected if an executed snapshot already exists. Unit tests run with `python -m unittest discover -s tests -p test_quantumequi_electronic.py -v`. The saved-data energy check runs with `python quantumequi/results/electronic/verify_saved_energy.py`.

## 中文方法与结果

### 范围与原始矩阵重建

本模块严格区分三类证据：原始类扩展休克尔矩阵的诊断、实际执行的 H2 电子结构计算、以及对这些 H2 参考值的监督插值。没有计算 Cu–N4 催化反应能，也没有训练原始 EGNN。原始八原子输入的元素清单为 C2HCuN4；源码和原始捕获结果均保持不变。审计通过 Python AST 提取原始电子求解器类，仅重建矩阵，不重复原神经势、NEB 或 Hessian 流程。

源码设置31个价轨道、40个价电子，因此占据20个双占据轨道；但其重叠公式只依赖原子间距和元素径向指数，没有轨道角向依赖。同一原子上的不同轨道因此重叠为1，且整行重叠矩阵完全相同。Cu的6轨道贡献15对，6个C/N原子的4轨道各贡献6对，共51对同原子非对角元全部为1。两个原始几何的重叠秩均为8，零空间维数为23。约−1.40e−15与−8.08e−16的最小本征值属于舍入尺度，不能消除重复行造成的精确缺陷。

源码将小本征值截为η，实质上把S改为Sη，再求解新度量下的问题。必须分别检查原方程残差HC−SCε、新方程残差HC−SηCε以及CᵀSC−I。上方裁剪表给出反应物η从1e−3到1e−8的六档结果；HOMO从10050 eV增大到1005000000 eV，所有裁剪下原度量正交误差最大值均为1。默认1e−5时，原广义残差最大值为3318.671454，能够重现原始471896257.61 kcal/mol的轨道和输出。这是可重复的数值病态，不是巨大的物理结合能。一电子参数矩阵的占据轨道和，也不能直接当作具有一致电子相互作用和核排斥处理的从头算总能。

### 经审查的广义本征求解器

`generalized_eigh` 检查矩阵同阶、有限与对称；严格模式拒绝不定、奇异或正定性不足的S。不定性判据使用明确的舍入容限64 εmachine n max(1,|s|max)，默认保留阈值为max(1e−12,1e−10 smax)。程序没有虚构角向轨道重叠模型。

规范正交化模式只保留正重叠子空间，使用X=Ur sr^(−1/2)求解XᵀHX；分别报告投影残差、全空间残差和被舍弃子空间的H耦合。这是限定子空间的Ritz控制。方法背景见[ SciPy `eigh` 官方规范](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html)及[Psi4正交化说明](https://psicode.org/psi4manual/master/scf.html)；文献支持方法，具体数值来自保存矩阵。

反应物/产物的投影残差分别为1.823481e−13与4.623372e−13，但原空间残差仍为4.593825与4.626220。8个保留轨道最多容纳16个闭壳层电子，无法容纳40个电子，因此占据数检查明确拒绝。降秩操作仅是数学诊断，不能作为“已经修复的物理EHT”。

六个正定随机控制的种子固定为901–906，与独立SciPy广义求解结果的最大本征值差不超过5.95e−14。兼容的半正定控制全空间残差小于5e−16；不兼容控制投影残差小于4e−16，而全空间残差为0.282843。20项单元测试覆盖这些区别、基变换与置换不变性、错误输入拒绝、电子容量和奇偶性、源矩阵简并、固定数据分区与径向解析导数；测试不运行Psi4或完整主研究。

### 实际量子参考

现有Psi4 1.11环境以单CPU线程、500 MB十进制内存完成中性单重态H2的PK RHF/STO-3G能量和解析梯度。采用C1、`no_com`、`no_reorient`，两个H位于z=±R/2 Å；能量和密度收敛阈值均为1e−12，SCF上限100步。实际环境为NumPy 2.5.2、SciPy 1.18.0。没有新安装、DFT、实验或Cu团簇从头算。

主扫描包含R=0.50+0.04i Å、i=0,…,32的33点，以及R=1.90+0.10i Å、i=0,…,8的9个拉伸点。拟合之前固定分区：偶数索引17点训练；模4余1的8点验证；模4余3的8点测试；较长距离9点为分布外。这里只检验同一分子的几何插值/外推，不代表跨分子或反应泛化。

保存总能、电子能及核排斥，关系为Etotal=Eelectronic+Enuclear、Enuclear=a0/R。Psi4采用a0=0.52917721067 Å；Hartree/bohr梯度除以a0转为Hartree/Å，再由dE/dR=(g2z−g1z)/2、F2z=−g2z得到键长导数与力。0.74 Å试跑总能为−1.1167593073781568 Hartree，电子能−1.8318636461214000 Hartree，核排斥+0.7151043387432432 Hartree；占据轨道和−1.1571077193736687 Hartree不同于电子能和总能。

主网格最低点为R=0.70 Å，总能−1.1173490349988597 Hartree，但梯度−0.0259704221 Hartree/Å非零，因此不能称为优化平衡构型。最大合力误差1.11e−15 Hartree/Å，Fock广义残差4.44e−16，MO度量正交误差1.11e−15。55个主作业均保存AO重叠、Fock、核哈密顿量、α密度、MO系数、轨道能量及实际请求的解析梯度。

在0.70、1.14、1.78 Å，分别以h=0.001和0.0005 Å独立计算中心能量差分，共12个位移单点。最大误差为2.220864e−6 Hartree/Å；同一位置步长减半后为5.552152e−7。位移点不进入训练。方法上参考[Psi4官方解析/差分梯度实例](https://github.com/psi4/psi4/blob/master/samples/fd-gradient/input.dat)，本次数值来自本机实际作业。

中性双重态H的UHF/STO-3G能量为−0.46658184955727544 Hartree，两原子参考值−0.9331636991145509 Hartree；2.70 Å处RHF H2仍为−0.6809407604974436 Hartree，显示受限最小基组曲线的解离局限。本次没有FCI、相关方法、基组收敛或H2热化学计算，不能把这条曲线当作准确解离参考。

### 径向监督替代模型与负结果

两个确定性RBF模型均以17个训练距离为高斯中心，并加入常数和线性项；一个拟合能量，一个拟合能量及dE/dR。力来自拟合能量的解析导数。训练能量和梯度各自的标准差用于缩放最小二乘块；岭惩罚作用于全部系数。每族预先给定4个长度尺度0.06/0.12/0.24/0.48 Å和2个岭值1e−10/1e−6，共16次线性拟合。以验证能量及梯度RMSE除以对应**训练**标准差之和选择参数；所以能量模型虽不拟合梯度，也允许验证梯度参与选参。两者均选中0.24 Å、1e−10。没有在验证集或测试集重新拟合，没有神经网络训练，也没有测试集调参。

三种基线使用相同训练能量：训练均值/零力、分段线性插值和外推、自然三次样条及其解析导数。分段线性模型在结点不可微；非结点测试不受此歧义影响。上方完整误差表保留两类RBF与三类基线的相同测试和分布外指标。

联合能量/力拟合的测试能量RMSE为5.561350e−6 Hartree，低于仅能量RBF的8.624253e−6；但是测试梯度RMSE为1.637466e−4 Hartree/Å，反而高于仅能量模型的1.156169e−4。拉伸外推时，两个RBF能量RMSE为0.175473与0.118352 Hartree，均明显差于简单线性外推的0.044450 Hartree；两者的分布外梯度误差甚至高于零力基线。这些负结果全部保留。插值误差小不代表可迁移的机器学习原子势，拟合RHF更准确也不能消除RHF自身误差。

### 账本、复现与证据边界

主量子计算55个SCF作业=42个H2扫描/解析梯度+12个差分能量+1个H原子；试跑另3个SCF和1个梯度，**合计58个SCF作业、43个解析梯度，失败SCF为0**。主量子用时15.81 s。矩阵审计包含2次源矩阵构造、12个裁剪控制、2个降秩控制、2个严格模式拒绝、2个电子容量拒绝及8个独立矩阵控制。监督学习为16个RBF拟合及3个基线，保存210行预测和20行指标；单元测试工作不计入研究作业。

成功试跑的精确源码副本`executed_script.py.txt`仍然保留；之后的代码变化仅涉及主机身份脱敏与溯源元数据。试跑历史源哈希没有被最终哈希覆盖。公开Psi4日志仅移除本机路径和主机身份，保留数值和SCF信息，同时保存脱敏前与公开文件哈希。已知计时文件留在本模块，归档为文本。一次错误测试工作目录发现0个测试，记录保留于`results/electronic/environment`；随后最终20项测试全部通过。

主结果与试跑`summary.json`的输入哈希相对`quantumequi`，输出哈希相对各自结果目录；`executed_source_snapshot`说明历史代码对应关系。源审计另外绑定原始capture，拟合摘要明确其本地输入哈希基准。独立`verify_saved_energy.py`检查全部58个作业的密度矩阵能量恒等式、α电子迹、Fock残差、日志SCF收敛和哈希、打印能量与保存能量，一致性全部通过；主结果最大密度能量/Fock误差8.88e−16 Hartree。这些是数值核验，不是外部物理精度验证。

复现需使用含Psi4的现有解释器执行`python quantumequi/scripts/electronic_reviewed.py --stage all`；小试跑使用`--stage quantum --pilot`。可分别运行`audit`、`quantum`、`fit`、`summarize`，最后一项不增加量子作业。主动重跑前应归档现有结果；已有执行快照且代码不同会被拒绝。单元测试命令为`python -m unittest discover -s tests -p test_quantumequi_electronic.py -v`；保存数据检查命令为`python quantumequi/results/electronic/verify_saved_energy.py`。
