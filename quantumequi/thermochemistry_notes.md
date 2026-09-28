# Reviewed Hessian and ideal-gas thermochemistry / Hessian 与理想气体热化学复核

## English methods and evidence

The executable API is [`scripts/thermochemistry_reviewed.py`](scripts/thermochemistry_reviewed.py). It introduces unit-explicit Hessian analysis and controlled statistical-mechanical calculations. It does not calibrate the user's neural potential or provide an electrocatalytic, solution-phase, or experimental free energy. The source script and saved weights remain untouched. The source TS coordinates are read from `results/original/captured_results.json`; the NEB is not rerun by this module.

The conversion from a mass-weighted Hessian eigenvalue in kcal mol⁻¹ Å⁻² amu⁻¹ is

\[
\tilde\nu=\operatorname{sign}(\lambda)\sqrt{|\lambda|}
\frac{\sqrt{(4184/N_A)/(10^{-20}m_u)}}{2\pi c_{\rm cm}}
=\operatorname{sign}(\lambda)\sqrt{|\lambda|}\,108.59135853528242\;\mathrm{cm^{-1}}.
\]

All constants are explicit, with atomic mass constant 1.66053906892×10⁻²⁷ kg from CODATA 2022. The source factor 1302.83 is 11.9975477 times this nominal factor. It is also not the Hartree/Bohr²/amu factor, which is 5140.48714361. The source neural energy has no established kcal/mol calibration; correcting dimensional bookkeeping does not supply one.

`rigid_subspaces` centers coordinates at the mass center and constructs √m-weighted translations and infinitesimal rotations. SVD provides an orthonormal rigid basis Q and its complement B. The internal Hessian is BᵀM⁻¹ᐟ²HM⁻¹ᐟ²B. A linear molecule has five rigid directions, a nonlinear molecule six, and an isolated atom three. Coincident multi-atom geometries are rejected. Candidate columns are normalized before SVD, whose relative rank tolerance is 10⁻⁹; this is a numerical geometry convention, not a general treatment of floppy, near-linear molecules. No first-six eigenvalue deletion is performed. Raw spectra, projected spectra, asymmetry and rigid-space residuals are preserved.

`harmonic_analysis` uses the maximum atomic gradient norm to require stationarity at 10⁻⁵ in the selected energy unit per Å. A 0.01 cm⁻¹ numerical threshold distinguishes near-zero internal modes from positive/negative modes. `ideal_gas_rrho` rejects nonstationarity, internal zero modes, or a classification mismatch. A minimum must have no negative internal mode; a first-order saddle must be requested explicitly and have exactly one, which is excluded from its constrained stable-mode partition function. Projection at a nonstationary geometry is a curvature diagnostic, not a transition-state confirmation.

For positive modes, x=hcν̃/(kBT), `expm1` expressions evaluate occupations and logarithmic partition functions stably. The implementation separates ZPVE, vibrational thermal enthalpy, entropy and heat capacity. Ideal-gas translation uses volume per molecule kBT/p; classical linear/nonlinear rigid rotation uses principal moments, an explicitly supplied symmetry number, and 1 or 3/2 RT rotational energy. Translation contributes 5/2 RT to enthalpy. Ground-state electronic degeneracy g contributes R ln g to entropy, with no modeled excited electronic levels. Thus G=E+ZPVE+Hthermal−TS. Pressure is explicit in Pa; the default is 1 bar, not 1 atm. Rotational characteristic temperatures and T/max(θrot) expose the classical-rotation approximation. No solvation, electrode potential, concentration-standard-state conversion, hindered rotor or anharmonic correction is applied.

## Controls, source findings and numerical results

Four invariant pair-distance polynomial controls use V=Σ[ki(di−di,0)²/2+8(di−di,0)⁴]. Two are a linear diatomic geometry with masses 12 and 16 amu, separation 1.4 Å and k=±70 kcal mol⁻¹ Å⁻². Two are a nonlinear triangle with masses 12,16,14 amu and k=(±70,50,40). Distances and coefficients define artificial potentials, not identified molecules. At their specified stationary geometries, independently assembled analytic Cartesian Hessians agree with autograd to 3.55×10⁻¹⁵ kcal mol⁻¹ Å⁻². The diatomic frequency independently follows k(1/m1+1/m2).

| Control | Rigid rank | Projected frequencies / cm⁻¹ | Classification |
|---|---:|---|---|
| Linear minimum | 5 | 346.954630 | Minimum |
| Linear saddle | 5 | −346.954630 | First-order saddle |
| Nonlinear minimum | 6 | 233.945192, 272.910673, 378.325758 | Minimum |
| Nonlinear saddle | 6 | −338.905360, 237.572002, 300.003873 | First-order saddle |

A displaced nonlinear geometry has maximum atomic gradient 7.97240477 and is explicitly rejected. Forty FD Hessians span four stationary controls, float32/float64, and displacements 0.02,0.005,0.001,0.0002,0.00004 Å. For the nonlinear minimum, float64 maximum Hessian error falls from 0.016672862 to 6.67722×10⁻⁸, whereas float32 error eventually increases to 0.0899200. A smaller displacement is not automatically more accurate.

The exact saved source model receives one float64 autograd Hessian and eight FD Hessians at the original candidate coordinates: two dtypes and four displacements, 0.02,0.005,0.001,0.0002 Å. At the source setting float32/0.005 Å, its eigensolver, factor and selection reproduce ZPVE=0.890357736 kcal/mol and Svib=106.649406 cal mol⁻¹ K⁻¹. Both negative raw eigenvalues are in the six entries deleted by index. Its returned two-negative count does not establish two physical unstable modes: a near-zero numerical eigenvalue changes sign with numerical details. Float64 autograd has three strictly negative raw eigenvalues, but two correspond to frequencies around −5×10⁻⁸ cm⁻¹. Thresholded projection retains one significant negative mode.

The projected source spectrum contains 18 modes, minimum −1.881083387 and largest 8.744653071 nominal cm⁻¹. Its maximum atomic gradient is 0.018339870823 nominal kcal mol⁻¹ Å⁻¹, approximately 1834 times the declared stationarity threshold. All eight FD cases retain one negative projected mode while failing stationarity. The implementation therefore publishes no corrected source RRHO free energy. The float64 FD maximum Hessian discrepancy drops from 3.14173×10⁻⁷ at 0.02 Å to 3.14671×10⁻¹¹ at 0.0002 Å. Float32 instead worsens from 4.54127×10⁻⁷ to 3.21271×10⁻⁵ against the promoted-weight float64 derivative. Promoting the same stored float32 weights changes arithmetic precision, not the learned model or its training status.

Additional source defects are preserved explicitly: 78.5 cal mol⁻¹ K⁻¹ is a fixed translational/rotational entropy baseline; thermal enthalpy is omitted; the electronic reference adds a rounded untrained-neural barrier to an unrelated EHT orbital sum; the plotted “free energies” instead add hardcoded +3.2/−1.8 kcal/mol and never use the RRHO result. No source saddle-point or chemical free-energy claim is accepted.

The ideal-gas sweep evaluates four controls at T=200,250,298.15,350,500 K and p=10⁴,10⁵,101325,10⁶ Pa, with symmetry number=1 and electronic degeneracy=1. These are explicit conventions for artificial controls. At 298.15 K and 1 bar, the nonlinear minimum has ZPVE=1.265430660 kcal/mol, thermal enthalpy=3.182947980 kcal/mol, entropy=61.494904977 cal mol⁻¹ K⁻¹, and G−E=−13.886327279 kcal/mol. Its G correction changes from −8.083490738 at 200 K to −27.048086494 at 500 K. Tests verify ΔG=RT ln(p2/p1), symmetry and degeneracy effects, mass scaling, stable partition limits and invariance.

Seventy-two separate low-frequency sensitivity cases replace the lowest of three frequencies with 0.1,1,5,10,20,50,100,200 cm⁻¹ at three temperatures, optionally flooring positive modes to 50 or 100 cm⁻¹. This is an explicitly arbitrary sensitivity intervention, not quasi-RRHO or a rotor model. With [0.1,500,1000] cm⁻¹ at 298.15 K, Svib=17.91393624 cal mol⁻¹ K⁻¹; raising the lowest mode to 100 yields 4.20598579 and shifts G−E by +4.09847880 kcal/mol. A convenient cutoff does not establish a physical low-frequency correction.

## Work, files, reproduction and verification

| Workload | Main | Separate pilot |
|---|---:|---:|
| Analytic autograd Hessians, including rejected geometry | 5 | 5 |
| Analytic FD Hessians / force evaluations | 40 / 600 | 16 / 240 |
| Analytic Hessian reverse-gradient rows | 39 | 39 |
| Exact-source FD Hessians / model force evaluations | 8 / 384 | 2 / 96 |
| Exact-source autograd Hessians / reverse rows | 1 / 24 | 1 / 24 |
| Ideal-gas T/P cases | 80 | 4 |
| Low-frequency sensitivity cases | 72 | 9 |

Each autograd Hessian uses one scalar energy forward and a first-gradient calculation followed by one reverse differentiation per Cartesian component; the source forward itself also computes its force. Counts describe this algorithm and exclude pilot and tests. Pilot wall time was 3.86 s, main 6.24 s inside the script, with one CPU thread; interpreter startup is additional. Twenty unit tests passed. Main code SHA256: `f2e0402f622c62e6b8ae118539309951be3f71c20a23ee9d4c91302ac2a6c77c`. The original pilot source snapshot is `results/thermochemistry/pilot/executed_code.py.txt`, matching its own code hash. Nine data artifacts in each summary passed SHA256 rechecks. Main spectra contain 378 rows: 9 methods × (24 raw +18 projected); the pilot predates the added autograd CSV rows, preserving its execution provenance.

- `results/thermochemistry/analytic_controls.json`: geometries, masses, classifications, full spectra.
- `analytic_hessian_convergence.csv`, `source_hessian_convergence.csv`: retained dtype/step/error records.
- `source_rrho_audit.json`: original results, exact reproduction, all raw negative results and rejection rationale.
- `source_spectra.csv`: `dtype=torch.float64_autograd`, `space=projected` selects the 18 source modes for plotting.
- `ideal_gas_temperature_pressure.csv`: `pressure_Pa=100000.0` selects the four artificial-control temperature scans.
- `low_frequency_sensitivity.csv`: positive-mode sensitivity, not corrected target thermochemistry.
- `analytic_hessians.npz`, `source_hessians.npz`: arrays, including source coordinates and gradient.

```powershell
$chem=$env:CHEM_PYTHON
$env:OMP_NUM_THREADS='1'; $env:OPENBLAS_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'
& $chem quantumequi/scripts/thermochemistry_reviewed.py
& $chem -m unittest discover -s tests -p test_quantumequi_thermochemistry.py -v
```

Set `CHEM_PYTHON` to the existing chemistry interpreter and configure its native runtime first. Run deliberate regeneration in a separate checkout; it is distinct from checking published hashes and requires new provenance/QA. 先配置解释器环境变量及原生运行库；有意重新计算应使用独立检出目录，不能把重算当作既有发布哈希检查。

Primary sources consulted: [NIST CODATA 2022 constants](https://physics.nist.gov/cuu/pdf/JPCRD2022CODATA.pdf), [Gaussian thermochemistry technical note](https://gaussian.com/wp-content/uploads/dl/thermo.pdf). The latter explicitly discusses stationary geometries, ideal-gas partitions and thermal enthalpy; its older numerical constants were not copied. No Gaussian calculation, external checkpoint download or reference-spectrum fitting was performed. Source text was user supplied; new numerical implementation is original, and dependency licenses remain upstream.

## 中文方法、结果与边界

本模块提供带单位的 Hessian、刚体投影与理想气体 RRHO 数值基准，不把随机神经势校准为化学能量，也不产生溶液、电极或实验自由能。直接读取原样运行保存的权重和候选坐标；没有重跑 NEB。按照 SI/CODATA 2022 常数，kcal mol⁻¹ Å⁻² amu⁻¹ 特征值的波数换算因子为 **108.59135853528242**；原始代码 1302.83 大约高 **11.99755 倍**，而 Hartree/Bohr²/amu 对应因子是 **5140.48714361**。单位纠正仍不等于势能标定。

在质量中心构造 √m 加权的平移和转动向量，由 SVD 得到刚体空间及正交补，投影后再对角化。线性体系删除的是五维刚体子空间，非线性体系为六维，单原子为三维；不会按排序盲删六个特征值。秩阈值为 10⁻⁹，几何近线性或柔性极强时仍需专门判断。保存完整原始/投影谱、非对称残差和刚体方向残差。驻点判据为最大原子梯度不超过 10⁻⁵，内部近零模式阈值为 0.01 cm⁻¹。非驻点、内部零模或负模数不匹配时拒绝 RRHO；一阶鞍点必须显式指定，且只排除唯一不稳定方向。非驻点的投影负模只能表示局部曲率。

理想气体分配函数分别处理平移、经典刚性转动、正振动模式和电子基态简并。平移使用单分子体积 kBT/p，焓包含 5/2 RT；线性/非线性转动分别贡献 RT、3/2 RT。采用 `expm1` 稳定计算低频/高频振动项，区分 ZPVE 与热激发振动焓，最终 G=E+ZPVE+H热−TS。默认压力为 **1 bar**，并明确给定对称数和电子简并，不能与 1 atm、1 mol/L 标准态混用。没有溶剂、电子激发态、阻转子、非谐性或电极处理；输出转动特征温度以暴露经典近似的适用程度。

四个人工控制势是平移/转动不变的原子对距离多项式，包含线性双原子和非线性三角形的极小值/一阶鞍点。独立解析 Hessian 与 autograd 最大差 **3.55×10⁻¹⁵**。投影频率分别为：线性极小值 **346.954630**；线性鞍点 **−346.954630**；非线性极小值 **233.945192、272.910673、378.325758**；非线性鞍点 **−338.905360、237.572002、300.003873 cm⁻¹**。人为移动后的非驻点梯度 **7.97240477**，被明确拒绝。40 个有限差分 Hessian 比较五个位移和两种精度；双精度误差随步长缩小下降，单精度最终受舍入影响反而增大。

精确源权重审计进行了八个有限差分 Hessian 和一个双精度 autograd Hessian。在原始 float32、0.005 Å 设置下复现 ZPVE **0.890357736 kcal/mol**、振动熵 **106.649406 cal mol⁻¹ K⁻¹**。排序前六项包含两个负值，原代码全部删除。严格负特征值个数对近零数值噪声敏感：双精度原始谱有三个严格负值，其中两个仅约 **−5×10⁻⁸ cm⁻¹**。刚体投影后的 18 个模式有一个显著负模 **−1.881083387 nominal cm⁻¹**，最高模式为 **8.744653071 nominal cm⁻¹**；但最大原子梯度 **0.018339870823** 是阈值约 **1834 倍**。因此不能确认过渡态，模块不发布所谓修正后的源自由能。

源代码还使用固定的 **78.5 cal mol⁻¹ K⁻¹** 平移/转动熵基线，遗漏热焓，把 EHT 轨道和与随机神经势垒混合，绘图却另行使用 **+3.2/−1.8 kcal/mol** 固定偏移，并未调用 RRHO 自由能。所有这些不一致均保留在审计 JSON 中。

温压扫描共 **80** 个案例：四个控制体系、五个温度和四个压力，均明确设置对称数=1、电子简并=1。298.15 K、1 bar 下，非线性极小值 ZPVE=**1.265430660**、热焓=**3.182947980 kcal/mol**、总熵=**61.494904977 cal mol⁻¹ K⁻¹**、G−E=**−13.886327279 kcal/mol**。另有 **72** 个低频敏感性案例。把 [0.1,500,1000] cm⁻¹ 中最低频率强制抬升至 100 cm⁻¹，振动熵从 **17.91393624** 降至 **4.20598579 cal mol⁻¹ K⁻¹**，G 校正增加 **4.09847880 kcal/mol**。这只是任意频率截断的敏感性演示，不是准 RRHO 或阻转子修正。

主计算包含 **600** 次解析控制势有限差分力计算、**384** 次源模型有限差分力计算；autograd Hessian 的反向行数为 **39+24**。pilot 和 20 项单元测试不计入这些主任务数字。pilot/主程序内部耗时约 **3.86/6.24 秒**，均单 CPU 线程。全部 20 项测试通过，包括压力关系 ΔG=RT ln(p2/p1)、同位素缩放、旋转平移不变性、线性秩、鞍点负模、极端频率稳定性和拒绝条件。主结果和 pilot 各九个数据文件的哈希通过复核；pilot 的执行代码另存快照。英文部分列出的 API、命令、数据列选择和主来源适用于中文版。
