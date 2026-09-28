# Finite-pool benchmark extension / 有限候选池基准扩展

This extension was executed on the existing **81-point deterministic transport pool**. It adds **zero PDE solves, zero experiments and zero instrument connections**. The calculation asks how selection policies and surrogate diagnostics behave on that frozen numerical surface; it does not add chemical validation, new reaction scope or industrial performance evidence.

本次扩展复用已冻结的 9×9 候选池，不把缓存读取称为新的实验或物理仿真。输入范围仍为流量 100–1500 µL/min、过电位 0.20–0.75 V；原始池由 80×32 输运网格生成，产物身份与电压关系仍为假设。完整运行共 **4800 条缓存评价使用记录**，而不是 4800 个新的独立工况。

## Execution and provenance / 执行与来源

A bounded pilot was run first: 2 seeds × 3 methods × 10 selections = 60 cached evaluation uses, with 8 held-out splits. It completed successfully and reproduced all 40 available overlapping records from the original campaign. One full run then completed in **6.748056 s** on the recorded CPU environment, with OMP, OpenBLAS and MKL each configured to one thread. The full run contains **64 seeds × 3 methods × 25 selections = 192 campaigns and 4800 uses**. The two objectives are maximized STY/50000 and negative electrical SEC, with the unchanged reference point (0, −2); the full-pool hypervolume is **1.5215781732523859**.

The original 8-seed GP/random study's first 15 evaluations were reproduced exactly: **240/240 candidate identities match**, and the maximum raw-hypervolume difference is **0**. Source hashes, the frozen pool hash, versions, timestamps, output hashes and oracle-budget checks are recorded in [summary.json](results/benchmark_extension/summary.json) and [oracle_budget_checks.json](results/benchmark_extension/oracle_budget_checks.json). The pool SHA-256 is `9a731e99a024f352ed7a51cff165aa7ae06e079bbdb99e66bf8e362c3e7a3ecf`.

原有 `metrics_control.py`、`transport_reviewed.py` 及 `results/control/` 均未修改。试运行结果保存在独立的 [pilot](results/benchmark_extension/pilot/summary.json) 目录，不能与完整运行重复累加为新的独立样本。完整运行使用 Python 3.12.14、NumPy 2.4.6、SciPy 1.18.0、scikit-learn 1.9.0；实际版本以来源记录为准。

## Sequential policy comparison / 序贯策略对照

Each seed uses the same 5 random initial points for all three policies. Existing GP MC-EHVI and uniform random selection are reused without changing their algorithm. MC-EHVI uses 256 posterior draws per remaining candidate. The added maximin policy selects the input with the largest minimum Euclidean distance from the already selected inputs in the declared, normalized two-dimensional input domain. It uses **no outcomes** to select a candidate, and exact ties resolve to the lowest remaining index. A budgeted oracle rejects duplicate calls and calls beyond 25; all 192 runs passed these checks.

三种方法分别为 GP MC-EHVI、无放回随机选择和 maximin 空间填充。每个种子共享初始点，随后分别选择自己的候选。maximin 仅使用输入坐标和已选索引，不能查看尚未评价点的目标。5、10、15、20、25 点结果均取自同一次 25 点运行的前缀，不另外执行，也不把相关前缀当成独立重复。

Mean fraction of full-pool hypervolume over 64 seeds:

| Evaluation budget | GP MC-EHVI | Random | Maximin |
|---:|---:|---:|---:|
| 5 | 0.876685929 | 0.876685929 | 0.876685929 |
| 10 | 0.973962450 | 0.941565731 | 0.962743109 |
| 15 | 0.981292266 | 0.962374910 | 0.977752064 |
| 20 | 0.987960570 | 0.972264126 | 0.985467926 |
| 25 | 0.992024355 | 0.982781182 | 0.989217617 |

At budget 25, the sample standard deviations of these fractions are **0.001733723**, **0.015598262** and **0.008467418**, respectively. The method ordering of the means does not imply that one method wins every paired run. Complete ranges and standard deviations are available in [learning_curve_summary.csv](results/benchmark_extension/learning_curve_summary.csv).

| Budget-25 paired difference in raw HV | Mean | Descriptive 95% t interval | Left wins / losses |
|---|---:|---:|---:|
| GP − random | 0.014064210 | [0.008194365, 0.019934055] | 45 / 19 |
| GP − maximin | 0.004270671 | [0.001022415, 0.007518927] | 33 / 31 |
| Maximin − random | 0.009793539 | [0.002746504, 0.016840574] | 42 / 22 |

At budget 15, **GP loses to maximin in 39 of 64 seeds**, although the average difference is positive: **0.005386695**, with interval **[−0.001078722, 0.011852111]**. At budget 20, GP also wins fewer pairs than it loses against maximin (28 versus 36). These negative results are preserved. Selection-policy performance depends on the finite surface, initialization and evaluation budget; no universal superiority or industrial optimum is established.

在 25 点预算下，GP 对随机策略 45 胜、19 负，对 maximin 33 胜、31 负；15 点预算时，对 maximin 仅 25 胜、39 负。即便平均值较高，也不能写成“稳定全面胜出”。配对区间仅描述固定候选池内的种子变化，依赖近似独立、正态种子效应假设；没有进行普遍优越性检验、p 值筛选或工业最优认定。多个预算和多组比较是相关的描述性分析，没有多重比较校正或因果结论。

The raw per-evaluation trace is [campaign_evaluations.csv](results/benchmark_extension/campaign_evaluations.csv); the 960 retained method/seed/budget checkpoints are in [budget_prefixes.csv](results/benchmark_extension/budget_prefixes.csv). All budgetwise pair comparisons are in [paired_HV_differences.csv](results/benchmark_extension/paired_HV_differences.csv).

## Held-out prediction diagnostics / 留出预测诊断

The predefined diagnostic design has **26 splits**, each containing 54 training and 27 test candidates:

- 20 random splits, with seeds 2000–2019.
- Three flow-block holdouts, each excluding three adjacent flow levels.
- Three overpotential-block holdouts, each excluding three adjacent overpotential levels.

Three fixed models are compared independently for each objective: a Matérn-5/2 GP with fixed normalized length scales (0.3, 0.3), unit amplitude and numerical nugget 10⁻⁸; a training-mean baseline; and ordinary least squares using `[1, x1, x2, x1², x1*x2, x2²]`. No test target is passed to the fitting routine and no hyperparameters are selected using a test set. All input minima/maxima and GP outcome centering/scaling are fitted using training rows only; held-out inputs receive the stored transform. Thus blocked holdout preprocessing intentionally differs from the campaign's a priori declared-domain normalization. The kernel family and fixed normalized hyperparameters are unchanged.

The separation of training-time transformations from test-time evaluation follows the [scikit-learn leakage guidance](https://scikit-learn.org/stable/common_pitfalls.html). Fixed kernels, target normalization and the numerical diagonal term are documented in the [GaussianProcessRegressor reference](https://scikit-learn.org/stable/modules/generated/sklearn.gaussian_process.GaussianProcessRegressor.html). The installed runtime version is recorded above; these documentation pages support methodological definitions, not any claim that the package validates this chemistry.

训练与测试索引、变换参数、二次回归系数和训练区间外标记保存在 [holdout_split_metadata.json](results/benchmark_extension/holdout_split_metadata.json)。低端或高端整块留出会超出**当前训练子集**的坐标范围，但仍位于原 81 点候选域内；不能据此宣称对候选域之外的新设备、新化学体系或新流量区间有预测能力。

Pooled held-out RMSEs, with the STY target scaled by 50000:

| Split type | Model | STY/50000 RMSE | −SEC RMSE (kWh/kg) |
|---|---|---:|---:|
| Random 54/27 | GP | 0.023814117 | 0.005867510 |
| Random 54/27 | Training mean | 0.317246326 | 0.091777371 |
| Random 54/27 | Quadratic | 0.073270727 | 0.015955436 |
| Flow blocks | GP | 0.230838320 | 0.059243346 |
| Flow blocks | Training mean | 0.351532797 | 0.098252572 |
| Flow blocks | Quadratic | 0.139597004 | 0.030398649 |
| Overpotential blocks | GP | 0.303768791 | 0.095824743 |
| Overpotential blocks | Training mean | 0.421253862 | 0.128631987 |
| Overpotential blocks | Quadratic | 0.184537455 | 0.040184884 |

**The quadratic baseline outperforms the fixed GP for both objectives in both blocked-holdout groups.** The GP's stronger random-split performance therefore does not demonstrate equally strong extrapolation across a missing flow or overpotential region. Hyperparameters were not retuned after observing this result.

显著的负结果是：流量分块和过电位分块留出中，简单二次回归的两个目标 RMSE 都低于固定 GP。随机留出表现较好，不能替代对跨区间预测的评价。本次保留该结果，没有根据留出误差再改核参数或挑选好看的拆分。

The GP intervals are latent-function bands `mean ± 1.96 × posterior_sd`. Standardized residuals use `(cached truth − prediction)/posterior_sd`. They assess this interpolation/extrapolation model on deterministic targets; the posterior SD is **not physical measurement uncertainty, PDE discretization error or transport-model discrepancy**.

| Split type | STY GP interval coverage | −SEC GP interval coverage | STY standardized-residual RMS | −SEC standardized-residual RMS |
|---|---:|---:|---:|---:|
| Random 54/27 | 100.0000% | 100.0000% | 0.287068 | 0.235736 |
| Flow blocks | 100.0000% | 100.0000% | 0.868403 | 0.749543 |
| Overpotential blocks | 77.7778% | 75.3086% | 1.399732 | 1.687733 |

随机留出的 100% 覆盖和过电位分块留出的约 75–78% 覆盖共同说明：名义 95% 区间的表现依赖拆分方式，并未完成概率校准。均值和二次模型未输出预测区间，相关字段保持空值。26 拆分含每目标/每模型 **702 条留出使用记录**，合计 **4212 行**；随机拆分之间重复使用候选，因此聚合统计不是 702 个新的独立样本。

Per-point means, intervals, residuals and memberships are in [holdout_predictions.csv](results/benchmark_extension/holdout_predictions.csv). The 156 split/model/target metric rows are in [holdout_split_metrics.csv](results/benchmark_extension/holdout_split_metrics.csv); the 18 split-type/model/target pooled rows are in [holdout_grouped_metrics.csv](results/benchmark_extension/holdout_grouped_metrics.csv). To convert scaled-STY RMSE to the assumed product mass basis, multiply by 50000; that unit conversion does not validate the assumed product identity.

## Reproduction and verification / 复现与核验

```text
python electratwin/scripts/benchmark_extension.py --pilot
python electratwin/scripts/benchmark_extension.py
python -m unittest discover -s tests -p test_electratwin_benchmark_extension.py -v
```

The default CLI performs the full predefined calculation. It does not need `--full` and does not modify the frozen input pool or earlier controller. Twelve tests passed after the full run, including exact budgets and no repeats; shared initialization; maximin invariance to changed outcomes; deterministic budget prefixes; disjoint complete split membership; training-only transforms; hand-computed MAE/RMSE, interval coverage and standardized residuals; recovery of an independent linear example by quadratic regression; all saved 240 original prefixes; and saved held-out-row membership.

本扩展增加的是可审计的算法对照、学习曲线与留出诊断。12 项测试证明相应软件不变量在本次环境中成立；它们不能证明电化学机理、反应选择性、产物结构或真实装置的可靠性。
