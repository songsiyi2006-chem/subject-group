# QuantumEqui: correlation, learned forces and anharmonic nuclei

[中文扩展报告](reports/extension_report_chinese.md) · [English extension report](reports/extension_report_english.md) · [Interactive curve explorer](reports/quantum_explorer.html) · [Seven figure groups](reports/figures/extensions/README.md)

本次扩展保留原报告与计算记录，增加实际电子相关计算、能量与梯度联合神经训练、同位素非谐振动及可逐项复算的误差分解。三个协作代理分别实现模块并交叉检查；主代理整合图表、交互页面、报告和发布验证。所有计数均区分主计算、试跑、数值检查和量子作业。

This extension is a reproducible H2/H benchmark. It does not validate Cu catalysis, an organic transition state, a general-purpose neural potential or experimental spectroscopy. FCI removes configuration truncation within each finite basis; basis error remains. The nuclear calculations use an FCI-parameterized Morse approximation, rather than the full FCI potential curve.

## Additional work actually executed

| Module | Main study | Separate pilot or follow-up |
|---|---|---|
| Psi4 RHF/UHF/FCI | 183 energy-driver calls: 150 scan, 2 atomic H, 21 equilibrium search, 10 curvature | 14 successful pilot calls + 6 FCI recovery calls; 6 earlier option rejections occurred before any driver call |
| Differentiable distance MPNN | 6 runs, 3 paired seeds, 800 epochs each, 4,800 optimizer steps | 2 pilot runs, 200 epochs each, train/validation only |
| J=0 nuclear finite differences | 38 unique eigenproblems, 6 final spectra, 36 bound-state thermal cases | 6 pilot eigenproblems + 7 long-tail follow-up eigenproblems |
| Matched-geometry error accounting | 126 rows, 21 shared geometries, no interpolation | No new quantum calculation or training |
| Graphics | 7 groups × 2 languages × PNG/SVG = 28 files | 14 PNGs inspected; all 14 SVGs independently rendered and inspected |

The **203 new energy-driver calls** are additional to the original release's 58 quantum jobs. Underlying SCF iterations and convergence messages are not counted again as extra jobs. The neural training reuses the original 42 RHF/STO-3G energy/gradient labels and adds no quantum labels. The nuclear module performs 45 main/follow-up eigenproblems, with 6 pilot problems counted separately.

## Evidence and limitations / 证据与边界

- 同基组、同几何下全部 50 组比较满足 FCI ≤ UHF ≤ RHF。拉伸后 UHF 的自旋污染明显；两个基组首次采样到较低破缺分支均在 1.26 Å，这不是精确分岔点。FCI 的自旋诊断没有测得，保留为空值。
- 三个配对种子中，联合训练均改善测试集与拉伸域外集的能量和梯度误差；全部神经模型的测试误差仍高于冻结的 RBF 基线。三种子离散度不等于校准置信区间。能量单任务模型的检查点选择仍使用验证集梯度。
- RHF 标签拟合得更好，不会自动消除 RHF 相对 FCI 的偏差。126 行分解显式区分学习误差、参考重算漂移与 RHF 偏差；本次重算漂移为零。
- 核振动的 8 Å 主区间漏掉 cc-pVDZ 参数化 H2 Morse 模型最浅的一个束缚态。96 Å、115,200 区间找回该态，但束缚能相对误差仍为 1.564%。细网格误差与势能模型误差分别讨论。
- 热力学量仅来自束缚振动态配分函数，报告 Helmholtz F；未包括连续谱、转动、平动、核自旋权重或溶剂，不是完整 Gibbs 自由能。

## Methods, records and review

| Topic | Methods and implementation | Saved evidence |
|---|---|---|
| Electronic correlation | [Bilingual notes](correlation_extension_notes.md), [code](scripts/correlation_extension.py) | [Main summary](results/extensions/correlation/summary.json), [curve](results/extensions/correlation/curve.csv), [driver ledger](results/extensions/correlation/energy_driver_jobs.json) |
| Neural learning | [Bilingual notes](learning_extension_notes.md), [code](scripts/learning_extension.py) | [Summary](results/extensions/learning/summary.json), [all predictions](results/extensions/learning/predictions.csv), [metrics](results/extensions/learning/metrics.csv) |
| Nuclear motion | [Bilingual notes](vibration_extension_notes.md), [code](scripts/vibration_extension.py) | [Main summary](results/extensions/vibration/summary.json), [tail follow-up](results/extensions/vibration/tail_followup/summary.json) |
| Error decomposition | [Code](scripts/error_budget_extension.py) | [Summary](results/extensions/error_budget/summary.json), [126 rows](results/extensions/error_budget/matched_error_components.csv) |
| Publication | [Plot code](scripts/plot_extensions.py), [explorer builder](scripts/build_explorer.py) | [Figure hashes](results/extensions/figure_manifest.json), [visual QA](results/extensions/figure_qa.json), [browser QA](results/extensions/explorer_qa.json) |
| Independent checks | [Validator](scripts/validate_extensions.py) | [Cross-review](results/extensions/cross_review.json), [validation](results/extensions/validation.json) |

Each scientific run retains its executed source snapshot and version/hash records. Some finalizer and tail-follow-up code was edited after the main calculation: historical execution hashes intentionally differ from the final source hash. This preserves the actual computation history rather than implying a new scientific run.

## Inspect and reproduce

Open `reports/quantum_explorer.html` in a browser. Its 150 points are embedded locally; no server, API key or network data is required. Both basis sets, a shared energy reference, method visibility and the saved-geometry slider are selectable. The page does not compute new quantum energies.

Saved-record verification from the repository root (NumPy, SciPy and Pillow required):

```powershell
python quantumequi/scripts/validate_extensions.py
```

To repeat scientific work, use a separate checkout and existing environments. `CHEM_PYTHON` needs NumPy, SciPy, PyTorch, Matplotlib and Pillow; `PSI4_PYTHON` needs Psi4 and SciPy, with their native library runtimes activated. The electronic run used Psi4 1.11. Set interpreter paths through your environment; no machine-specific paths are embedded here.

```powershell
$chem=$env:CHEM_PYTHON
$qm=$env:PSI4_PYTHON
$env:OMP_NUM_THREADS='1'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:PYTHONIOENCODING='utf-8'
& $qm quantumequi/scripts/correlation_extension.py --pilot
& $qm quantumequi/scripts/correlation_extension.py
& $chem quantumequi/scripts/learning_extension.py --pilot
& $chem quantumequi/scripts/learning_extension.py
& $chem quantumequi/scripts/vibration_extension.py --pilot
& $chem quantumequi/scripts/vibration_extension.py
& $chem quantumequi/scripts/vibration_extension.py --tail-followup
& $chem quantumequi/scripts/error_budget_extension.py
& $chem quantumequi/scripts/plot_extensions.py
& $chem quantumequi/scripts/build_explorer.py
```

The corrected fresh pilot already executes its FCI calls: do **not** invoke `--pilot-fci-recovery` afterward. That mode exists only to recover the recorded historical pilot's six rejected options. Recomputed output has new hashes and requires renewed numerical and visual review; old QA records are not automatically valid for a new execution.
