# Paired budget analysis / 配对预算分析

Post-hoc analysis of existing cached ElectraTwin optimization records, using only Python's standard library. No solver, quantum calculation, scheduler or network call is made.

本模块只分析已有数据：同一缓存候选池上的 3 种算法、64 个配对随机种子、5 个记录预算（5、10、15、20、25），共 192 条曲线。随机种子重复不是独立化学实验，也不验证底层物理模型。

From the repository root / 在仓库根目录执行：

```sh
python -m scientific_analysis.paired_budget --output work/paired-review
python -m scientific_analysis.validate_release
python -m unittest discover -s tests -p test_paired_budget_analysis.py -v
```

The output directory must not already exist. Optional arguments: `--input`, `--bootstrap-samples` (default 4000), `--seed` (default 20261012). Saved results include input and analyzer SHA256 hashes.

## Methods / 方法

- Resample whole matched seed pairs for percentile 95% confidence intervals of mean differences. Intervals are pointwise, conditional on this one candidate pool and exchangeable seeds; they are not simultaneous confidence intervals or evidence of cross-system generalization.
- Compute normalized trapezoidal AUC over recorded budgets 5–25. Intermediate values are linearly interpolated; budgets 0–4 are not invented.
- Report exact two-sided sign tests, excluding numerical ties within 1e-12. Apply Holm correction across all 12 exploratory comparisons (two baselines × five budgets plus AUC). This tests the probability of a positive difference, not the mean difference.
- Report leave-one-seed-out mean ranges to expose sensitivity to individual seeds.
- Record threshold crossing bounds for a 0.99 full-pool hypervolume fraction. Crossings between prefixes are interval-censored; initial attainment is left-censored; failure by budget 25 is right-censored. Unreached targets never receive a fictitious crossing budget.

## Recorded results / 已有数据的分析结果

Differences below are GP MC-EHVI minus baseline, in normalized hypervolume fraction units; CI = paired-bootstrap mean CI. 全部预算的比较均保存在输出中，表格只摘录最终预算和曲线面积。

| Baseline | Endpoint | Mean difference | 95% CI | Wins / 64 | Holm sign p |
|---|---|---:|---|---:|---:|
| Random | Budget 25 | 0.009243 | [0.005725, 0.013088] | 45 | 0.018753 |
| Random | Normalized AUC | 0.017908 | [0.011805, 0.024821] | 44 | 0.040590 |
| Maximin | Budget 25 | 0.002807 | [0.000931, 0.004995] | 33 | 1.000000 |
| Maximin | Normalized AUC | 0.004664 | [0.002062, 0.007427] | 31 | 1.000000 |

平均收益和逐种子胜率回答不同问题：相对 Maximin 的平均差为正，不意味着多数种子或所有预算都稳定获胜。These exploratory results do not establish universal algorithm superiority.

| Method | Reached 0.99 by budgets 5 / 10 / 15 / 20 / 25 | Unreached at 25 |
|---|---|---:|
| GP MC-EHVI | 0 / 0 / 0 / 14 / 56 | 8 |
| Random | 0 / 0 / 2 / 15 / 31 | 33 |
| Maximin | 0 / 0 / 1 / 18 / 47 | 17 |

Source: [cached budget prefixes](../electratwin/results/benchmark_extension/budget_prefixes.csv). Outputs: [complete JSON](results/analysis.json), [paired comparisons](results/comparisons.csv), [per-seed curves and crossing bounds](results/curve_metrics.csv).

Statistical references: [SciPy exact binomial test documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.binomtest.html); Holm, S. (1979), *A Simple Sequentially Rejective Multiple Test Procedure*, Scandinavian Journal of Statistics 6(2), 65–70 ([original article](https://www.jstor.org/stable/4615733)).
