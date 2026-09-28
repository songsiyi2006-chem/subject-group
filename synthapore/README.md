# SynthaPore-Omni: executable audit and numerical benchmarks

2026-09-28 · CPU computation · separate English and Chinese reports

- **[中文完整报告](reports/synthapore_report_chinese.md)** · **[Complete English report](reports/synthapore_report_english.md)**
- [Four bilingual PNG/SVG figures](reports/figures/README.md)
- [Supplied specification](source/specification.md) · [original code](source/synthapore_engine.py) · [minimal compatibility patch](source/compatibility.patch)
- [Source output JSON](results/compatibility/synthapore_flagship_results.json) · [exact source replay audit](results/source_audit.json)
- [Publication validation](results/publication_validation.json) · [agent cross-review](results/cross_review.json) · [visual inspection](results/figure_qa.json)

三个代理分工完成等变学习、动力学/NEB、几何/吸附计算，并交叉核对保存数据。原程序因 Cu 原子序数超出嵌入表而失败；保留失败日志，只将入口处 `num_species=10` 改为 `30` 后完成兼容执行。此修改会改变随机数消耗，兼容结果不是原程序的成功输出。代码全程仅作软件计算。

The original fails on Cu atomic number 29 in a ten-entry embedding. The compatibility entry point increases the embedding to 30; its extra random draws also change downstream initialization. The captured 55,394-parameter potential has never been trained on energies or forces. It describes a manually supplied eight-atom C4H2CuN fragment, not a chemically specified complete catalyst/substrate. Source outputs are preserved as claims to audit, not validated predictions.

## Executed additions

| Module | Main workload | Evidence and retained limits |
|---|---|---|
| Equivariant denoising | 368 synthetic shapes; 3 seeds × 60 epochs; 17,664 predicted-coordinate rows | Held-out and warped-shape tests, identity/PCA/fitted-smoother baselines; no reverse diffusion or molecular generation |
| Analytic MD | 44 trajectories; 350,000 integration steps | NVE timestep convergence and 24 thermostat replicas on an artificial harmonic model |
| Analytic CI-NEB | 18 reviewed cases; 12 unchanged-source comparison cases | Force stopping, optimized endpoints and Hessian checks on a known two-dimensional potential |
| Geometry | 1,636,800 point evaluations; 90 probe masks | Three artificial two-dimensional geometries, grid convergence and analytic controls |
| Synthetic adsorption | 152 BET fits; 384 perturbed observations | Window and noise sensitivity, negative constants retained; no measured isotherm |

Main counts exclude pilots, unit tests and mathematical network audits. The reports and each module summary itemize those separately. The exact captured-source NEB replay has 315 evaluations and seven final-band checks. It performs no additional MD.

Three denoisers improve held-out error against the fitted smoother, but seed 4441 is worse on warped shapes. Corrected Berendsen temperature control still collapses kinetic fluctuations; BAOAB gives a variance/reference ratio of 0.982727 on this finite benchmark. All 18 reviewed analytic NEB cases converge; 10 of 12 source comparison cases claim convergence without satisfying the independent force threshold. Source pore radius is unused, and its synthetic adsorption formula jumps at relative pressure 0.35. These negative results are retained.

## Code and detailed methods

- [Equivariant implementation](scripts/equivariant_reviewed.py), [bilingual methods](equivariant_notes.md), [summary](results/equivariant/summary.json), [actual captured-weight audit](results/equivariant/captured_source_audit.json).
- [MD and CI-NEB implementation](scripts/dynamics_reviewed.py), [bilingual methods](dynamics_notes.md), [summary](results/dynamics/summary.json).
- [Geometry and BET implementation](scripts/pore_reviewed.py), [bilingual methods](pore_notes.md), [summary](results/pore/summary.json).
- [Original failure log](results/original/stderr.log), [compatibility execution record](results/compatibility/execution.json), [source audit implementation](scripts/audit_source.py).

## Reproduction

Use an existing environment containing NumPy, SciPy, PyTorch, pandas, Matplotlib and Pillow. Recorded versions are in the module summaries. CPU threads are limited to one by the launcher. No new installation is performed. On Windows, use a Python environment whose native scientific DLL dependencies are already on `PATH`.

Run from the repository root:

```bash
python synthapore/scripts/run_suite.py --dry-run
python synthapore/scripts/run_suite.py
python -m unittest discover -s tests -p "test_synthapore_*.py" -v
python synthapore/scripts/validate_synthapore.py --require-cross-review
```

The launcher orders source execution, source replay audit, equivariant calculations, captured-weight audit, dynamics, pore calculations and figures. `--modules` accepts a subset in that fixed order; omitted prerequisites must already exist. `--timeout-seconds` bounds each stage. The individual stages were executed for this release; the convenience launcher was checked with `--dry-run`, not used to repeat the entire study.

Each reviewed calculation script supports `--pilot`. Regeneration replaces result files. Timing and serialization bytes can change across runs or environments; hashes identify the published artifacts, not a promise of universal bitwise determinism. A new run requires refreshed numerical validation, cross-review and visual inspection before publishing updated hashes. The validator does not retrain models or rerun study trajectories.

## Scientific scope

Source potential energies have only nominal units. Differentiable forces and equivariance do not establish physical accuracy. This release has no trained chemical MLIP, reverse diffusion generator, chemical free-energy barrier, atomistic POP structure, GCMC adsorption calculation or experimental calibration. The analytic 1 eV saddle is deliberately specified in a benchmark potential. The pore areas are two-dimensional probe-center accessibility, not measured material porosity or a pore-size distribution.

The deliverable is a reproducible numerical study and source audit with separate bilingual reports. Its tests support software and arithmetic integrity; they do not establish chemical mechanism, publication readiness or industrial performance.
