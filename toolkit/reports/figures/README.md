# Reviewed figures / 经审计图表

All data are synthetic, mock or illustrative as visibly labeled. These files are not measured substrate scope, experimental optima or computed DFT free energies. Preserve labels and captions when reusing the figures. PNG is a 300-DPI preview; SVG has editable text and vector geometry without embedded raster images.

| Figure | English PNG | English SVG | 中文 PNG | 中文 SVG |
|---|---|---|---|---|
| 1. Exact Pareto front and assumed SEC contours | [PNG](Fig1_Pareto_Electrosynthesis_english.png) | [SVG](Fig1_Pareto_Electrosynthesis_english.svg) | [PNG](Fig1_Pareto_Electrosynthesis_chinese.png) | [SVG](Fig1_Pareto_Electrosynthesis_chinese.svg) |
| 2. Mock scope heatmap and descriptive correlation | [PNG](Fig2_Substrate_Scope_Heatmap_english.png) | [SVG](Fig2_Substrate_Scope_Heatmap_english.svg) | [PNG](Fig2_Substrate_Scope_Heatmap_chinese.png) | [SVG](Fig2_Substrate_Scope_Heatmap_chinese.svg) |
| 3. Hypothetical energy-level diagram | [PNG](Fig3_Reaction_Energy_Profile_english.png) | [SVG](Fig3_Reaction_Energy_Profile_english.svg) | [PNG](Fig3_Reaction_Energy_Profile_chinese.png) | [SVG](Fig3_Reaction_Energy_Profile_chinese.svg) |
| 4. Synthetic analytical assumption tests | [PNG](Fig4_Analytical_Quality_Control_english.png) | [SVG](Fig4_Analytical_Quality_Control_english.svg) | [PNG](Fig4_Analytical_Quality_Control_chinese.png) | [SVG](Fig4_Analytical_Quality_Control_chinese.svg) |

Full captions and interpretation are in the [English](../deployment_toolkit_report_english.md) and [Chinese](../deployment_toolkit_report_chinese.md) reports. Figure 1 uses the exact nondominance relation, not the original threshold. Figure 2 contains hardcoded arrays, not fitted-model predictions and experiments. Figure 3 plots supplied constants, with no DFT/TS evidence. Figure 4 tests known Gaussian shapes, their mismatch and assumed relaxation; it is not instrument validation.

English uses Arial; Chinese uses Microsoft YaHei. Editable text depends on local font availability; a missing font may be substituted. Verify final dimensions and fonts in the target journal workflow. The [manifest](../../results/figure_manifest.json) records dimensions, DPI, figure hashes and source CSV hashes; [QA](../../results/figure_qa.json) records actual visual review scope.
