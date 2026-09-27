# Research engine: three-module execution and quantitative audit

**2026-09-27.** The supplied script completed unchanged. Its original results are preserved. The “research-grade” filename does not imply validated chemical predictions.

- [中文完整报告](reports/research_grade_technical_report_chinese.md)
- [Complete English report](reports/research_grade_technical_report_english.md)
- [Section-aligned bilingual report / 双语合并版](reports/research_grade_technical_report_bilingual.md)
- [Original specification](source/specification.md), [original code](source/run_research_engine_original.py), [original JSON](results/original/research_grade_results.json)
- [Additional audit summary](results/audit/audit_summary.json), [publication checks](results/publication_validation.json)

Reports follow the new brief's five sections: abstract, multiobjective electrosynthesis, complete eight-substrate table, POP particle transport, and lab/grant roadmap. Earlier report editions remain in their own directories.

| Work actually executed | Scope | Main evidence boundary |
|---|---|---|
| Original pipeline | 135 conditions, 12 Pareto records, 8 structures, 10 particle cases | Targets and transport parameters are assigned; molecular yield labels are random |
| Multiobjective audit | Rounded/unrounded frontier, 231 weight sets, 3 constraint and 3 area scenarios | Source choice depends on preferences and missing area/cost assumptions |
| Added multi-output GP | 5 contiguous current-band folds; 108 train/27 test per fold; mean baseline | Fits an artificial response surface; no sequential BO or measured chemical training set |
| Molecular pilot | 64 MMFF conformers; 8 neutral optimizations plus 8 cation GFN2-xTB/ALPB single points | Actual semiempirical jobs, not DFT, measured Eox or validated site selection |
| Random-label diagnosis | 1,000 seeds × 8 molecules | Reproduces the label generator, not substrate-specific chemistry |
| Pore-model verification | 40 mesh, 50 external-film, 60 parameter cases; threshold roots | Conservative finite-volume/analytical checks of a prescribed first-order sphere |

**关键发现：** 原脚本没有拟合 GP；“绿色折中”为合成目标评分。八个底物的产率是正态随机数，评分不依赖三维坐标。100 µm 微孔颗粒的 η=0.4890，因此硬编码的“超过 50 µm 就低于 0.40”结论错误；相应模型阈值为 129.379256 µm。不存在已完成的湿实验或保证国家级立项、获奖、论文接收的证据。

## Reproduction and files

Use the existing environment recorded in [audit_summary.json](results/audit/audit_summary.json) and [repository package record](../provenance/requirements_recorded.txt), plus xTB 6.7.1. No dependencies were installed for this delivery. On Windows activate the chemistry environment, including its `Library/bin` on PATH; otherwise numerical DLLs may fail. Set `XTB_EXE` when needed. These commands refresh generated records, so use a separate clone to preserve a release.

```powershell
$env:OMP_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
python research/scripts/run_original.py
python research/scripts/audit_research.py
python research/scripts/plot_research.py --zh-font C:/Windows/Fonts/msyh.ttc
python research/scripts/build_reports.py
python research/scripts/validate_research.py
python -m unittest discover -s tests -v
python scripts/validate_release.py
```

Use an available CJK font on other platforms. Figures are already supplied as PNG and SVG. Each xTB job is serial with a 180-second timeout; original-script execution has a 300-second timeout. Runtime varies by host. Audit functions run only under the main entry point; importing them in tests does not rerun xTB.

```text
source/                  Unmodified extracted code, brief, source hashes
results/original/        Original JSON, exit record and logs
results/audit/           Complete tables, per-fold predictions, profiles, summary
results/molecules/       Eight identities, source geometries, conformers and 16 xTB jobs
scripts/                 Reproduction, extended computations, figures and report validation
reports/                 Three report editions, shared numerical tables and templates
```

Molecular atom indices are zero-based in the supplied SMILES order. Neutral singlet and cation doublet states share the same optimized neutral nuclei. Charged ALPB states have separate equilibrium solvent responses; energy differences are not reference-electrode potentials or rigorous nonequilibrium-solvent ionization energies. No Hessian, charged geometry relaxation or experimental calibration was performed. Site charge responses are model-dependent, not reaction probabilities.

The original source hash distinguishes its unchanged scientific code from the added analysis. Public job logs normalize machine-specific executable/directory paths while preserving scientific values. Fourteen new tests independently check Pareto dominance, SEC bookkeeping, labels, grouped metrics, molecular identity, xTB convergence/charge conservation, diffusion integration, mesh convergence and film flux balance. Root tests include the preceding 29 tests. Automated checks and seven visually inspected figures establish publication/numerical integrity, not chemical validity.

The root [SHA-256 manifest](../provenance/SHA256SUMS.txt) covers this release. Fresh runs may change timings or numerical files and legitimately invalidate hashes. Review intentional changes, stage them, run `python scripts/update_manifest.py`, then validate again. Do not regenerate the manifest to hide unexplained differences.
