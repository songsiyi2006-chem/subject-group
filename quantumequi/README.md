# QuantumEqui-NEB: executed computation and evidence audit

**[中文完整报告](reports/quantumequi_report_chinese.md)** · **[Complete English report](reports/quantumequi_report_english.md)** · [Figures](reports/figures/README.md)

原脚本按原样执行成功；源代码、JSON、300-DPI 原图、实际随机权重和日志全部保留。新增真实 H2/H Hartree–Fock 计算、能量/梯度代理拟合、等变性与导数核验、CI-NEB 解析势基准、质量加权 Hessian 投影和显式理想气体热化学。三个代理分工实现并交叉核对；主代理完成整合、图表和发布检查。

The source executes, but its EHT overlap has rank 8 for 31 functions and cannot represent its assumed 40 electrons. The 24,385-parameter EGNN is untrained. Its NEB candidate is nonstationary, and no valid chemical barrier or source Gibbs energy is released. The actual quantum reference is small H2/H RHF/UHF in STO-3G; analytic NEB and thermal controls verify numerical methods without validating the proposed Cu chemistry.

## What was executed

| Component | Main calculation | Separate pilot |
|---|---|---|
| Psi4 RHF/UHF | 55 SCF jobs, 42 analytic gradients | 3 SCF jobs, 1 gradient |
| Distance RBF | 16 fits, 2 selected models, 3 baselines, 210 predictions | None |
| Captured EGNN | 1,209 energy/force calls, 48 symmetry probes, 24 full Cartesian FD checks | None |
| Source NEB replay | 420 old-band evaluations plus 14 fresh final-image evaluations | None |
| Analytic CI-NEB | 40 runs, 8,522 updates, 104,262 point evaluations | 6 runs, 1,097 updates, 11,449 evaluations |
| Analytic-control FD Hessians | 40 Hessians, 600 force calls | 16 Hessians, 240 calls |
| Source FD Hessians | 8 Hessians, 384 force calls | 2 Hessians, 96 calls |
| Ideal-gas RRHO / low-frequency sensitivity | 80 / 72 cases | 4 / 9 cases |

Study counts exclude unit tests and distinguish original execution from audits. Five analytic-control and one source autograd Hessians are additional. Electronic matrix controls are not quantum jobs. “Energy-only” RBF fitting still uses validation gradients for hyperparameter selection. Joint fitting worsens test force error relative to energy-only fitting; both extrapolate worse than linear interpolation/extrapolation. Negative results remain in the reports.

## Code and saved records

- [Exact specification](source/specification.md), [source provenance](source/source_record.json), [unchanged script](source/quantum_egnn_neb_engine.py), [execution](results/original/execution.json), [requested original JSON](results/original/quantum_neb_results.json), [original panel](results/original/figures_quantum/Fig_QuantumEquiNEB_Comprehensive_Panel.png).
- [Electronic code](scripts/electronic_reviewed.py), [methods](electronic_notes.md), [summary](results/electronic/summary.json), [quantum ledger](results/electronic/quantum_job_accounting.csv), [independent saved-energy check](results/electronic/saved_energy_identity_check.json).
- [Potential audit code](scripts/potential_audit.py), [methods](potential_notes.md), [summary](results/potential/summary.json).
- [Path code](scripts/path_reviewed.py), [methods](path_notes.md), [source replay](results/path/source_audit.json), [analytic summary](results/path/summary.json), [saved-data audit](results/path/read_only_audit.json).
- [Thermochemistry code](scripts/thermochemistry_reviewed.py), [methods](thermochemistry_notes.md), [summary](results/thermochemistry/summary.json), [source rejection](results/thermochemistry/source_rrho_audit.json).
- [Cross-review record](results/cross_review.json), [figure provenance](results/figure_manifest.json), [figure QA](results/figure_qa.json), [publication validation](results/publication_validation.json).

## Verify this release

From the repository root, use a Python environment containing NumPy, SciPy, PyTorch, RDKit, Matplotlib and Pillow. Saved-record checks do not require Psi4 and do not repeat SCF jobs:

```powershell
python quantumequi/scripts/validate_quantumequi.py --require-cross-review
python -m unittest discover -s tests -p 'test_quantumequi_*.py' -v
```

The actual main environment used Python 3.12.14, NumPy 2.4.6, SciPy 1.18.0 and PyTorch 2.10.0. The electronic environment used Psi4 1.11 with NumPy 2.5.2. Exact versions and code hashes are in each calculation record. Use existing installed environments; no package installation is needed to inspect the release.

## Recompute in a separate checkout

Set `CHEM_PYTHON` and `PSI4_PYTHON` to suitable interpreter executables and activate their required native-library runtime before running. The two interpreters may be the same if all dependencies are installed. These commands reproduce the individual stages that were actually executed:

```powershell
$chem=$env:CHEM_PYTHON
$qm=$env:PSI4_PYTHON
$env:OMP_NUM_THREADS='1'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:PYTHONIOENCODING='utf-8'
& $chem quantumequi/scripts/run_source.py
& $qm quantumequi/scripts/electronic_reviewed.py --stage all
& $chem quantumequi/scripts/potential_audit.py
& $chem quantumequi/scripts/path_reviewed.py --source-audit
& $chem quantumequi/scripts/path_reviewed.py
& $chem quantumequi/scripts/thermochemistry_reviewed.py
& $chem quantumequi/scripts/plot_reviewed.py
```

The source seeds NumPy and PyTorch at 42; cross-version or cross-platform bitwise identity is not guaranteed. Recalculation overwrites stage outputs and invalidates the published hashes and visual/cross-review attestations. Preserve the release in its original checkout. Reassess regenerated results and figures before issuing new attestations. Historical pilot snapshots record the exact code executed and must not be relabeled as executions of later code. No end-to-end wrapper rerun is claimed.

The generalized solver explicitly rejects invalid metrics or electron capacity. Thermal release requires a stationary geometry with the requested saddle order. Tests and saved-data checks establish bounded numerical behavior; they do not establish chemical accuracy, model transferability, molecular saddle connectivity or publication readiness.
