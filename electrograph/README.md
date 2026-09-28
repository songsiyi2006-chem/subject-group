# ElectroGraph-kMC: executed calculations and evidence audit

2026-09-28 · CPU calculations · separate complete English and Chinese reports

- **[English technical report](reports/electrograph_report_english.md)**
- **[中文完整技术报告](reports/electrograph_report_chinese.md)**
- [Eight PNG and eight editable SVG figures](reports/figures/README.md)
- [Saved-evidence publication validation](results/publication_validation.json) · [Figure inspection record](results/figure_qa.json) · [Agent cross-review](results/cross_review.json)

## What ran / 实际完成

Three agents implemented and cross-reviewed molecular graph learning, conformer/reaction analysis and independent-site kinetics. The coordinating agent preserved the supplied program, executed it, audited its claims, assembled figures and verified the release. Peer-agent review is software assistance, not independent human peer review or experimental validation.

| Main module | Actual calculation | Saved evidence |
|---|---|---|
| [Graph learning](scripts/graph_learning.py) | 256 experimental FreeSolv records; 4 neural fits, 240 epochs; 1,536 prediction rows; 384 cached search calls | [Summary](results/learning/summary.json), [notes](learning_notes.md) |
| [Conformer/CGR analysis](scripts/structure_reviewed.py) | 504 requested conformers, 120 returned and optimized, 36 pooled distinct minima; 2,880 SASA orientation evaluations | [Summary](results/structure/summary.json), [notes](structure_notes.md) |
| [Stochastic kinetics](scripts/kinetics_reviewed.py) | 273 trajectories, 9,950,504 SSA events; 105 analytic rate perturbations | [Summary](results/kinetics/summary.json), [notes](kinetics_notes.md) |
| [Source audit](scripts/audit_source.py) | 12 initialization probes plus one fixed network; 72 initialization and 192 priority forward passes | [Audit](results/source_audit.json) |

Main totals exclude pilots and unit tests. Pilots add two neural fits/four epochs/96 cached calls, two force-field optimizations and five SSA trajectories/169,216 events. The learning pilot's historical source is preserved beside its results. Calls reuse a fixed dataset; these are not new measurements. No new DFT, oxidation-potential measurements or wet experiments were performed.

学习结果保留负结果：三次 MPNN 测试 RMSE 为 1.8972、1.7468、1.2496 kcal/mol，岭回归为 1.6534；两次 MPNN 未超过基线。学习特征 GP 与描述符 GP 在八个搜索种子的最终结果均打平。水合自由能模型不支持氧化电位或原子反应性结论。

构象权重来自 MMFF94 势能与单位简并度假设，不是溶液自由能布居。反应原式 C22H21NOS → C21H17NOS 缺少 C1H4；仅在显式且唯一映射的平衡教学反应上输出键变化。

动力学实现为五状态、不可逆、独立位点模型，没有晶格邻居或传输。给定速率下的稳态期望上限为 26.61034847 s⁻¹，源码硬编码 260.2 s⁻¹ 不受模型支持。SSA 与 CTMC 一致只能证明数值实现一致，不能验证机理。

## Source preservation / 原始源码与兼容修复

[Specification](source/specification.md) and [original Python](source/electrograph_kmc_core.py) retain the supplied bytes. The [original execution](results/original/execution.json) fails at a nonexistent `HybridizationType.AROMATIC` enum. The [patch](source/compatibility.patch) makes exactly four installed-API repairs: remove that unused enum, pass ETKDG parameters through their object, use `CalcSASA(..., confIdx=...)`, and use `Minimize(maxIts=...)`.

The [compatibility execution](results/compatibility/execution.json) succeeds; an [intermediate failed attempt](results/api_repair_attempt_1/execution.json) is retained. Repairs preserve the original scientific assumptions. The source network remains untrained, source FreeSASA radii remain zero for all 30 target atoms, priority scores include random exploration, and the source potential–TOF panel is hardcoded. Its archived figure is not an endorsed reviewed research figure. [Source record](source/source_record.json) connects exact hashes and logs.

## Reproduction / 复现

Use the repository's [requirements](../requirements.txt) with Python 3.12. The executed environment reused PyTorch 2.10.0, RDKit 2026.03.5, NumPy 2.4.6 and SciPy 1.18.0, without CUDA. No automatic installation or hardware connection occurs. Set numerical threads to one for comparable runs; the launcher does this for its subprocesses.

From the repository root:

```text
python electrograph/scripts/run_suite.py --dry-run
python electrograph/scripts/run_suite.py
python electrograph/scripts/validate_electrograph.py
python -m unittest discover -s tests -p "test_electrograph_*.py"
```

The launcher runs `source → audit → learning → structure → kinetics → figures` in dependency order. `--modules learning structure kinetics figures` selects a subset, using saved prerequisites. `--timeout-seconds 1800` sets the per-stage limit. The delivered calculations ran individually; the wrapper was checked in dry-run mode. Regeneration replaces outputs and can change elapsed times, floating-point values or figure metadata. Hash validation or recorded visual QA must not be silently treated as current after regeneration.

The original source failure is retained and expected by the source runner. A compatibility-run failure stops the launcher. Archived learning data allow ordinary runs without downloading. Data preparation is separate: see [data provenance and license](data/learning/README.md). Official FreeSolv data carry CC-BY 4.0 attribution and original-source caveats; the upstream software license is distinct. Some experimental uncertainty fields are defaults/estimates. GAFF-computed values are excluded from targets and input features. Checkpoints are locally trained and should be loaded using the supplied restricted weight-loading code.

## Units, checks and limits / 单位、核验与限制

- Hydration targets, energies and model error: kcal/mol. Frozen-pool search maximizes negative hydration free energy.
- SASA: Å²; mass-weighted radius of gyration: Å; conformer temperature: K. SASA uses positive van der Waals radii, a 1.4 Å probe and 24 orientation averages; residual quadrature error remains.
- Rate constants and TOF: s⁻¹; potential index: V; trajectory times: s. Current describes only enumerated sites without electrode area/site-density conversion.
- 49 added tests: 16 learning, 18 structure and 15 kinetics. Full-repository counts are in the [release record](../provenance/release_validation.json).
- The [validator](scripts/validate_electrograph.py) checks saved hashes, independent summary arithmetic, accounting, language/table parity and local links. It does not rerun training or SSA and does not establish physical accuracy.

The ring-based split prevents overlap of defined groups, while acyclic molecules use complete canonical structures as groups; it is not a test of all forms of scaffold novelty. Eight-seed fixed-pool search intervals are descriptive. Conformer coverage, free-energy populations, electrochemical labels, calibrated kinetic rates, reverse reactions and spatial coupling remain open. The reports retain these limits rather than assigning publication probability or industrial readiness.
