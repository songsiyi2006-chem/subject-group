# ElectraTwin figures / 中英文图表

All reviewed figures represent uncalibrated numerical models. PNG files are 2160 × 1140 at 300 dpi; SVG files contain editable text and vector geometry without embedded raster images. Fonts are referenced, not embedded: Arial for English and Microsoft YaHei for Chinese. Substitute fonts require layout review.

| Figure | English | 中文 |
|---|---|---|
| 1. Concentration field and differently defined current models / 浓度场及不同电流模型 | [PNG](Fig1_Transport_Field_english.png) · [SVG](Fig1_Transport_Field_english.svg) | [PNG](Fig1_Transport_Field_chinese.png) · [SVG](Fig1_Transport_Field_chinese.svg) |
| 2. Mesh refinement and source FE inconsistency / 网格加密与原 FE 不一致 | [PNG](Fig2_Numerical_Audit_english.png) · [SVG](Fig2_Numerical_Audit_english.svg) | [PNG](Fig2_Numerical_Audit_chinese.png) · [SVG](Fig2_Numerical_Audit_chinese.svg) |
| 3. Model Pareto front and paired sequential comparison / 模型 Pareto 前沿与配对搜索对照 | [PNG](Fig3_Sequential_Optimization_english.png) · [SVG](Fig3_Sequential_Optimization_english.svg) | [PNG](Fig3_Sequential_Optimization_chinese.png) · [SVG](Fig3_Sequential_Optimization_chinese.svg) |
| 4. Reaction network and overoxidation / 反应网络与过氧化 | [PNG](Fig4_Reaction_Network_english.png) · [SVG](Fig4_Reaction_Network_english.svg) | [PNG](Fig4_Reaction_Network_chinese.png) · [SVG](Fig4_Reaction_Network_chinese.svg) |
| 5. First-order and total-effect sensitivity / 一阶与总效应敏感性 | [PNG](Fig5_Parameter_Sensitivity_english.png) · [SVG](Fig5_Parameter_Sensitivity_english.svg) | [PNG](Fig5_Parameter_Sensitivity_chinese.png) · [SVG](Fig5_Parameter_Sensitivity_chinese.svg) |
| 6. Learning curves and paired seed effects / 学习曲线与配对种子效应 | [PNG](Fig6_Optimization_Extension_english.png) · [SVG](Fig6_Optimization_Extension_english.svg) | [PNG](Fig6_Optimization_Extension_chinese.png) · [SVG](Fig6_Optimization_Extension_chinese.svg) |
| 7. Random and blocked holdout / 随机及分块留出 | [PNG](Fig7_Heldout_Prediction_english.png) · [SVG](Fig7_Heldout_Prediction_english.svg) | [PNG](Fig7_Heldout_Prediction_chinese.png) · [SVG](Fig7_Heldout_Prediction_chinese.svg) |

图 1 的原模型与修订模型使用不同的电子数和边界定义，不能当成同一物理模型的独立校准对照。图 2 的网格收敛不证明动力学真实。图 3 保留随机法胜出的三个种子；全部 81 点均为模型值。

See the full [English captions](../electratwin_report_english.md), [中文图注](../electratwin_report_chinese.md), [data hashes](../../results/figure_manifest.json) and [visual-QA scope](../../results/figure_qa.json). The source program's figures remain archived separately and were not included in this reviewed figure set.

Figures 4–7 have separate [English captions](../extension_report_english.md), [中文图注](../extension_report_chinese.md), [data and figure hashes](../../results/extension_figure_manifest.json) and [visual-QA record](../../results/extension_figure_qa.json). Figure 4 uses the new hypothetical reaction network; Figures 5–7 retain the original single-reaction model. Finite-sample S>ST estimates and blocked-prediction failures are shown without adjustment. Figure 6 intervals describe seed effects on the fixed pool, not physical uncertainty.
