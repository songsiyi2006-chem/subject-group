# ElectraTwin-OS: reviewed numerical demonstrator

完整独立的 [中文报告](reports/electratwin_report_chinese.md) · [English report](reports/electratwin_report_english.md) · [PNG/SVG figures](reports/figures/README.md)

新增独立的 [中文计算扩展报告](reports/extension_report_chinese.md) · [English calculation extension](reports/extension_report_english.md)

本目录执行并审计附件中的五个模块，补充守恒有限体积模型、真实逐步 GP MC-EHVI 推荐、绿色指标边界分析、水力和热量情景，以及仅在内存运行的 SCPI 状态机。它是可复核的数值研究演示；尚无湿实验、分子 DFT、真实仪器连接或工业验证。

The supplied program's industrial and authorship claims are preserved as source material, not endorsed. The original run failed at removed NumPy `trapz`; changing that call alone to `trapezoid` enabled a compatibility run. Its numerical and scientific defects remain visible in the source audit.

## Results and evidence

| Work completed | Inspectable record | Scope |
|---|---|---|
| Original failure and minimally repaired execution | [source record](source/source_record.json), [audit](results/source_audit.json) | 12-point fixed grid; 8 FE values above 100% were clipped by source; no original GP/EHVI |
| 39 conservative transport solves | [verification](results/transport/verification.json) | Baseline conversion 58.3872858750%, current 42.2513750193 mA; numerical checks, not calibrated chemistry |
| 81-point offline pool and one metric-baseline PDE solve | [optimization](results/control/optimization_summary.json), [pool](results/control/candidate_pool.csv) | Assumed kinetics, voltage and product mass |
| 8 seeds × 2 methods × 15 selections | [campaign records](results/control/sequential_campaigns.csv) | 240 cached evaluation uses, not 240 new PDE jobs or experiments |
| 27 hydraulic and 9 electrical-heating scenarios | [engineering](results/engineering/summary.json) | Selected fluid properties; no CFD or thermal PDE |
| Local SCPI state machine | [trace](results/control/scpi_simulator_trace.json) | Zero physical connections; final output off |
| 3 figures × 2 languages × 2 formats | [manifest](results/figure_manifest.json), [visual QA](results/figure_qa.json) | Numerical evidence labels retained |

MC-EHVI wins 5 of 8 paired seeds; random selection wins 3. Mean paired hypervolume difference is 0.0304070243, with a descriptive 95% t interval of [−0.0084974830, 0.0693115316]. This does not establish reliable superiority. FE is structurally 100% in the single-reaction model, so optimizing STY against FE degenerates; the reviewed benchmark uses STY and electrical SEC instead.

原模型物料流与电流重算差异为 15.0351%。修订模型的 A/P 摩尔守恒并不证明假设产物的完整元素或质量守恒。分子量尚无已确认结构与配平反应依据；未完整计入后处理、纯化和回收，不能称为全流程 PMI。给定电位只是有效动力学参数，原单反应模型尚未描述温度、气体、电势、迁移、副反应和寿命；本次扩展仅新增假设副反应网络。

## Extended code and executed calculations / 新增代码与已执行计算

| Module | Code and saved results | Executed scope |
|---|---|---|
| Explicit four-species reaction network | [reaction_network.py](scripts/reaction_network.py), [summary](results/reaction_network/study_summary.json), [method](network_notes.md) | 62 A/P/B/D solves + 1 frozen-model comparison; parallel and consecutive side reactions, 49-point scan and grid/analytic limits |
| Scenario propagation and sensitivity | [uncertainty_analysis.py](scripts/uncertainty_analysis.py), [summary](results/uncertainty/summary.json), [method](uncertainty_notes.md) | 1,792 Sobol-design + 16 mesh-check + 21 local-derivative PDE solves; 500 paired-row resamples use saved outputs |
| Larger sequential and prediction benchmark | [benchmark_extension.py](scripts/benchmark_extension.py), [summary](results/benchmark_extension/summary.json), [method](benchmark_extension_notes.md) | 64 seeds × 3 methods × 25 selections = 4,800 cached uses; 26 holdout splits, 4,212 prediction rows; separate 60-use pilot |
| Reproduction and publication | [run_extensions.py](scripts/run_extensions.py), [plot_extensions.py](scripts/plot_extensions.py), [saved-output validator](scripts/validate_extensions.py) | Sequential launcher, four figures in two languages and PNG/SVG; [figure hashes](results/extension_figure_manifest.json), [visual QA](results/extension_figure_qa.json), [validation](results/extension_validation.json) |

合计新增 1,892 次研究 PDE 求解，不含单元测试求解；缓存调用不重复计为 PDE 或实验。网络基准 A 转化率 58.4087%，净 P 收率 11.5438%，净 P FE 11.3871%，显示所设过氧化通道的显著影响。敏感性和优化扩展沿用原单反应模型，尚未与新网络耦合。有限 N 的 Sobol 估计出现 S>ST；未截断或冒充严格置信区间。分块留出中二次回归优于固定 GP，不能由随机拆分精度推断外推性能。

The network's abstract species, electron counts and rates are assumptions. The uncertainty module uses chosen independent distributions, not fitted experimental uncertainty. Seed comparisons and overlapping holdout predictions concern the same finite numerical pool. The extension adds 42 tests (17 network, 13 sensitivity, 12 benchmark) to the original 36 ElectraTwin tests; root release validation runs the entire repository suite.

From the repository root, preview the launcher commands or validate saved files:

```text
python electratwin/scripts/run_extensions.py --dry-run
python electratwin/scripts/validate_extensions.py
```

To regenerate all three calculation modules, or just one module, preserve the delivered files first:

```text
python electratwin/scripts/run_extensions.py
python electratwin/scripts/run_extensions.py --modules network
python electratwin/scripts/plot_extensions.py
```

The launcher uses the current interpreter, runs modules sequentially with single-thread numerical settings and stops on a failed module or timeout. Each calculation module was executed separately for this release; the launcher command construction was checked with `--dry-run`. Plot regeneration does not renew the saved visual-QA attestation. No command connects a physical instrument.

## Reproduce the original modules on an existing CPU environment

Use the repository's [requirements](../requirements.txt) and an already installed compatible Python environment. Recorded versions: Python 3.12.14, NumPy 2.4.6, SciPy 1.18.0, pandas 2.3.3, Matplotlib 3.11.1, scikit-learn 1.9.0 and RDKit 2026.3.5. The source imports RDKit but processes no molecular structures. Set OMP, OPENBLAS and MKL thread counts to one when reproducing the bounded CPU run.

From the repository root, validate the saved release without rerunning simulations:

```text
python electratwin/scripts/validate_electratwin.py
python -m unittest discover -s tests -p "test_electratwin*.py" -v
```

To regenerate calculations, first retain an unchanged release copy because these commands overwrite their result files:

```text
python electratwin/scripts/run_source.py --timeout-seconds 600
python electratwin/scripts/transport_reviewed.py
python electratwin/scripts/metrics_control.py
python electratwin/scripts/cross_review.py
python electratwin/scripts/engineering_audit.py --transport-json electratwin/results/transport/baseline_solution.json
python electratwin/scripts/plot_reviewed.py
```

`run_source.py --audit-only` recomputes the audit from saved arrays. `metrics_control.py --skip-optimization` regenerates arithmetic sensitivity and the local simulator trace only. The default optimization includes its offline pool in the saved accounting. Rerun timestamps and elapsed times can change file hashes; the saved release manifest and visual-QA attestations belong to the archived bytes. New figures require renewed visual review before updating the QA record. Validators must not be used to silently relabel regenerated files as previously reviewed.

## File and unit guide

| Location or field | Meaning |
|---|---|
| `source/` | Supplied specification, exact extracted script and compatibility-only copy |
| `results/original/`, `results/compatibility/` | Separate execution logs and original outputs; their plot titles are not scientific endorsements |
| `scripts/transport_reviewed.py` | Cell-centred conservative 2D A/P transport, half-cell Robin flux, axial upwind convection and two-direction diffusion |
| `results/transport/baseline_solution.json` | `field.A_mol_m3` and `P_mol_m3` indexed x then y; 1 mol/m³ = 1 mM; current in A |
| `grid_convergence.csv`, `screening_domain_grid_check.csv` | Conversion differences in percentage points; no global mesh-independent claim |
| `scripts/metrics_control.py` | Explicit mass/energy metrics, exact 2D hypervolume, observed-only GP campaign and simulator |
| `results/control/candidate_pool.csv` | Flow in µL/min; potential index in V; STY in kg/(m³ day); SEC in kWh/kg; FE in percent |
| `results/control/metric_boundary_sensitivity.json` | Unclipped inconsistent FE and null undefined metrics; recovered-product and reactor-product denominators separated |
| `results/engineering/` | SI hydraulic calculations and selected-property sensible-heating scenarios |
| `results/sources.json` | Primary-source links and access limitations |
| `results/publication_validation.json` | Structure, numeric consistency, links, provenance and figure checks; not physical validation |

The transport verification executes 39 solves; the control module executes 81 pool solves plus one repeated baseline. These are 121 executed PDE solves across the two workflows, with overlaps in physical conditions. Additional unit-test solves and source-model runs are separate. Three engineering tests, fifteen transport tests and eighteen metric/control tests cover numerical behavior. The earlier [analytical toolkit](../toolkit/README.md) and [Chinese grant proposal draft](../toolkit/proposals/National_Undergraduate_Grant_Proposal_PanTang_Lab.md) remain available; no grant has been submitted.
