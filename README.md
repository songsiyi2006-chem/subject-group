# AI4S Organic Electrosynthesis & Catalyst Informatics

完整技术报告、可复现计算与证据审计。Updated 2026-09-28.

## Latest ElectroGraph-kMC / 分子图学习、构象系综与随机动力学

- **[中文完整报告](electrograph/reports/electrograph_report_chinese.md)** · **[Complete English report](electrograph/reports/electrograph_report_english.md)**
- [代码、数据与复现入口](electrograph/README.md) · [中英文 PNG/SVG 图表](electrograph/reports/figures/README.md)

三个代理完成并交叉审查新增模块：256 条实验水合自由能数据上的 4 次神经训练（240 轮）、384 次缓存搜索调用、120 个实际构象优化、2,880 次 SASA 方向积分，以及 273 条 SSA 轨迹（9,950,504 个事件）。原脚本失败、四处 API 修复及兼容运行均保留。新增 49 项单元测试。

Three neural seeds give test RMSEs of 1.8972, 1.7468 and 1.2496 kcal/mol against ridge regression at 1.6534. Learned-feature and descriptor GP searches tie at their final budgets. The supplied reaction loses C1H4, and the assumed kinetic cycle saturates near 26.61 s⁻¹ rather than the hardcoded 260.2 s⁻¹. These calculations do not validate oxidation potentials, atom reactivity, solution populations or an electrocatalytic mechanism.

## Earlier ElectraTwin-OS / 微通道传输与多目标计算平台

- **[新增计算中文报告](electratwin/reports/extension_report_chinese.md)** · **[Extended calculations: English report](electratwin/reports/extension_report_english.md)**
- [一键计算入口](electratwin/scripts/run_extensions.py) · [扩展代码、结果与复现](electratwin/README.md)
- **[最新中文报告](electratwin/reports/electratwin_report_chinese.md)**
- **[Latest English report](electratwin/reports/electratwin_report_english.md)**
- [计算、证据与复现](electratwin/README.md) · [中英文 PNG/SVG 图表](electratwin/reports/figures/README.md)

新增 39 次传输验证求解、81 点模型候选池、8 种子相同预算 MC-EHVI/随机对照、27 组水力和 9 组热量情景。修订基准转化率 58.3873%，电流 42.2514 mA；原程序物料—电流差异 15.0351%，12 点中 8 个 FE 被截断。MC-EHVI 5 胜 3 负，配对区间包含零。全部属于未经实验校准的模型；SCPI 仅内存模拟，尚未获得工业验证。

The release preserves the original failure and compatibility execution, implements a conservative model and a real sequential finite-pool acquisition calculation, and keeps numerical verification separate from physical evidence. The 240 campaign records are uses of 81 precomputed model values, not 240 experiments.

本次进一步增加 1,892 次 PDE 求解：四物种反应网络及独立对照 63 次，五参数敏感性、网格与局部导数 1,829 次；另有 64 种子 × 3 方法 × 25 次的 4,800 条缓存调用和 26 个留出拆分的 4,212 条预测记录。扩展报告与四组新图保留严重过氧化、有限样本敏感性估计误差和分块预测失败。固定 GP 在两类分块留出中均未优于二次回归；这些结果来自假设模型，未增加湿实验或实体硬件连接。

The extension implements an explicit A/P/B/D transport network, scrambled-Sobol scenario propagation and a larger cached sequential benchmark. The 1,892 new study PDE solves exclude unit-test solves; the 4,800 campaign uses and separate 60-use pilot reuse the prior pool. New calculations do not calibrate either model against experiments.

## Analytical toolkit / 分析定量、学术图表与大创申报

- **[最新中文报告](toolkit/reports/deployment_toolkit_report_chinese.md)**
- **[Latest English report](toolkit/reports/deployment_toolkit_report_english.md)**
- [完整中文大创申报书](toolkit/proposals/National_Undergraduate_Grant_Proposal_PanTang_Lab.md) · [Complete English proposal](toolkit/proposals/National_Undergraduate_Grant_Proposal_PanTang_Lab_English.md)
- [中英文 PNG/SVG 图表](toolkit/reports/figures/README.md) · [计算、证据与复现](toolkit/README.md)

原程序已执行，新增 100,000 次输入误差情景、80 组峰恢复、2,000 次 Pareto 扰动和带来源检查的 HPLC/NMR 预积分 CSV 解析。原 HPLC 产率 1292.00% 与 NMR 91.392% 明显不一致；峰宽设错使模拟面积平均绝对误差升至 14.61%。图表保留模拟/示意标识，未冒充湿实验或 DFT。申报书为完整内容稿，导师同意、平台权限和下一轮申报窗口仍待核实，未正式提交。

The reviewed toolkit adds exact dominance, analytical failure checks, editable vector graphics and separate full English/Chinese reports and proposals. Numerical success and public equipment listings do not establish assay validity, laboratory access or grant approval.

## Earlier closed-loop platform / 前一版湿实验与 DFT 输入闭环平台

- **[最新中文报告](closed_loop/reports/closed_loop_experimental_report_chinese.md)**
- **[Latest English report](closed_loop/reports/closed_loop_experimental_report_english.md)**
- [Complete bilingual report / 完整双语版](closed_loop/reports/closed_loop_experimental_report_bilingual.md)
- [Calculations, DFT inputs and reproduction / 计算、DFT 输入与复现](closed_loop/README.md)

按附件第 3 节的五部分结构完成。新增 81 组计量情景、留一验证、EI/UCB 对照、敏感性分析及 12 份配套 DFT 输入。六条反馈仍为模拟标签；标准化 GP 的 MAE 为 15.26 个百分点，未优于均值基线的 13.98。尚未执行湿实验或 DFT，SOP 缺少已定义产物、偶联当量和淬灭/分析依据，保留为设计草案。

The supplied three modules were executed and their claims audited. Reviewed quantum templates cover both the actual wet-lab substrate and the original indoline control, neutral and radical-cation states, with explicit SMD and dependent refinement steps. Input generation and numerical tests do not establish chemical validity.

## Earlier three-module research engine / 前一版三模块研究引擎

- **[最新中文报告](research/reports/research_grade_technical_report_chinese.md)**
- **[Latest English report](research/reports/research_grade_technical_report_english.md)**
- [Bilingual report / 双语合并版](research/reports/research_grade_technical_report_bilingual.md)
- [Execution, audit and reproduction / 执行、核验与复现](research/README.md)

原脚本已按原样完成；补充 231 组权重、多输出 GP 分组验证、64 个构象、16 个 xTB 作业及独立有限体积扩散核验。报告修正原始“超过 50 µm 即 η<0.40”的错误结论，并明确底物产率来自随机标签。没有湿实验、准确电位预测或生产成熟度的验证。

The latest five-section report preserves all three modules and their original outputs, with independent numerical checks and a proposed one-year measurement plan. Synthetic objective values and random yield labels remain distinct from actual molecular calculations and experiments.

## Earlier four-task pipeline / 前一版四任务流程

- **[新版中文报告](production/reports/production_technical_report_chinese.md)**
- **[New English report](production/reports/production_technical_report_english.md)**
- [Section-aligned bilingual report / 按章节对齐的双语版](production/reports/production_technical_report_bilingual.md)
- [Execution records, methods and reproduction / 执行记录与复现](production/README.md)

新附件的原脚本在 SASA 接口处失败；归档修订版已完成四任务。新增 30 种子 × 4 策略 BO 对照、32 个构象、SAC 分组留出、流动解析解与电流平衡核验、目标分子的 3 个 xTB 作业。物性嵌入优势和 C3 位点主张均未得到这些模型的支持。A/C 仍是合成目标，D 使用未校准参数；没有湿实验验证，也尚未达到生产可用性。

The new extension preserves the failed original and the explicit repair, includes negative results, and distinguishes synthetic scores, molecular calculations and proposed experiments. It follows the new four-section report specification. The earlier five-topic study below remains available with its original evidence boundaries.

## Earlier five-topic study / 原五课题报告

- **[中文完整报告](reports/technical_report_chinese.md)**
- **[Complete English report](reports/technical_report_english.md)**

两份报告均保留原规范第 3 节的三个主体部分：执行概要与实验室契合性、五课题各四部分分析、本科四阶段路线。附录 D、E 补充新计算。原始脚本与结果保持原字节，代码中的硬编码结论和合成数据标签在报告中明确审计。

Both editions preserve the requested report architecture. Original benchmark outputs remain unchanged. Appendices D and E add numerical robustness studies and real-molecule semiempirical calculations.

## What was actually computed / 实际计算范围

| Calculation | Scale | Main result and boundary |
|---|---|---|
| Paired Bayesian optimization | 30 seeds × 3 methods × 20 evaluations | Mean recommendation regret: random 3.5557, original GP 1.3570, scaled GP 1.3208 synthetic yield points; scaled vs original difference inconclusive |
| Regression generalization | 360 fits; 5,000 independent test rows | At n=60, test R²: boosting 0.9287, linear 0.9921; synthetic targets, not volts |
| SAC ranking sensitivity | 100,000 draws at each of 3 noise levels | Co–N₄ ranks first in 5.841% at original SD; generator sensitivity, not catalyst selection |
| Classification diagnosis | 30 training seeds; 5,000 held-out rows | Balanced accuracy 0.9129 overall, 0.7276 near boundaries; no real reaction labels |
| Sequential flow scenario | 191 grid points; 10,000 kinetic scenarios per flow | Best qualifying grid point 0.55 mL/min under assumptions; higher-flow rows can violate charge availability |
| Real molecular descriptors | 6 molecules; 36 GFN2-xTB jobs | 6 neutral optimizations, 30 single points; fixed-nuclei gas/ALPB comparison, charges, RDKit descriptors |

**证据边界：** 原五课题及扩展统计使用合成数据；流动参数属于未校准假设。新增六分子计算是实际执行的 GFN2-xTB 半经验计算，**不是 DFT、实测氧化电位、反应势垒、区域选择性验证或工业条件推荐**。没有湿实验记录。

**Evidence boundary:** Executed software and semiempirical calculations are distinguished from literature facts, hypotheses, and experiments. No experimental reaction, electrochemical calibration, DFT transition state, or industrial optimum is established.

## Files

```text
reports/                 Complete English and Chinese reports
reports/figures/         English/Chinese PNG and SVG figures
data/original/           Byte-preserved supplied script, original output, audit
scripts/                Additional computation, plotting, release validation
results/extended/       Per-run CSVs, test data, flow grid, summary JSON
results/molecular/      Molecular identities, XYZ, xTB logs/JSON, atomic charges
tests/                  Scientific invariants and independent source checks
provenance/             Environment, tests, publication checks, hashes
```

`data/original/benchmark_audit.json` describes the original run only. Its statements about missing repeated-seed tests do not describe the later extension. Current evidence is in `results/extended/summary.json` and `results/molecular/summary.json`.

## Reproduce

Use an existing environment when available. The recorded run used Python 3.12.14, NumPy 2.4.6, SciPy 1.18.0, pandas 2.3.3, scikit-learn 1.9.0, Matplotlib 3.11.1, RDKit 2026.03.5, and xTB 6.7.1. See [recorded package versions](provenance/requirements_recorded.txt). `requirements.txt` gives compatibility ranges, not a guarantee of identical numerical output across versions. xTB is an external executable, not supplied by the Python requirements.

Activate the environment before running. **On Windows, Conda activation must put that environment's `Library/bin` on PATH**, otherwise numerical DLLs may fail even if imports succeed. No new dependencies were installed for the delivered run. If preparing a fresh environment, install the Python requirements and obtain xTB from its [official project](https://github.com/grimme-lab/xtb).

From the repository root, set thread counts to one (PowerShell):

```powershell
$env:OMP_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
# If xtb is not on PATH: set XTB_EXE to your installed xtb executable.
python scripts/extended_benchmarks.py --output work/reproduction/extended --seeds 30
python scripts/molecular_descriptors.py --output work/reproduction/molecular --timeout 180
python -m unittest discover -s tests -v
python scripts/validate_release.py
```

The `work/` outputs are separate from the committed results. The molecular script overwrites outputs at its explicitly selected destination; use a new directory when preserving previous runs. Each xTB job is serial and limited to 180 seconds by default. The statistical suite took approximately one minute on the recorded host; runtime is machine-dependent.

Tests read committed results, independently refit representative models, check paired BO budgets, verify the original hashes, compare flow equations against an independent integrator, enforce charge/material bookkeeping, account for every reliability-bin row, and inspect all 36 xTB job records. Tests do not require rerunning xTB. Passing them verifies numerical and record integrity, not chemical prediction accuracy.

`provenance/SHA256SUMS.txt` covers the reviewed release files. The manifest itself and mutable `release_validation.json` / `test_results.txt` execution records are excluded. After intentionally changing and reviewing artifacts, stage the files, run `python scripts/update_manifest.py`, and validate again. Do not regenerate hashes merely to silence an unexplained mismatch. Archived source and calculation files preserve exact bytes through Git; authored Markdown and Python use LF.

To rerun the original benchmark without modifying its archived JSON:

```powershell
New-Item -ItemType Directory -Path work/original-reproduction -Force
Copy-Item data/original/simulate_all_topics.py work/original-reproduction/
Copy-Item data/original/audit_benchmarks.py work/original-reproduction/
Push-Location work/original-reproduction
python simulate_all_topics.py
python audit_benchmarks.py
Pop-Location
```

To regenerate figures from committed data:

```powershell
python scripts/plot_extended.py --zh-font C:/Windows/Fonts/msyh.ttc
```

On other operating systems pass an installed CJK font file. The already committed figures include embedded outlines for Chinese SVG text. Figure rendering does not rerun the calculations. Hash checks will correctly flag changed artifacts; review intentional changes before updating the manifest.

## Reading the added molecular results

Molecules: benzene, pyridine, anisole, indole, N-methylindole, benzofuran. Every molecule has explicit SMILES, atom order, charge/spin and coordinates. Neutral GFN2-xTB/ALPB(acetonitrile) geometries are reused without changing nuclei for charged and gas-phase calculations. Each charged ALPB state has its own equilibrium solvent response; the energy differences are not rigorous nonequilibrium-solvent vertical ionization energies or reference-electrode potentials. No Hessian, thermal corrections, ionic relaxation, or experimental calibration was performed.

Public xTB logs/JSON retain scientific values, with only machine-specific executable and working-directory paths normalized. Atom indices correspond to saved coordinates, not automatic chemical site names. Full methods, primary references, numerical tables, and proposed experimental validation are in each complete report.

![Added validation overview](reports/figures/extended_validation.png)
