# Reviewed research figures / 核验后的研究图表

Four figures, separate English and Chinese editions. Each PNG is 2400 × 1380 pixels at 300 dpi. SVG files preserve editable text and vector geometry, with no embedded raster. Font substitution remains possible on computers without Arial or Microsoft YaHei.

| Figure | English | 中文 |
|---|---|---|
| 1. Source audit: hardcoded TOF and random ranking | [PNG](Fig1_Source_Audit_english.png) · [SVG](Fig1_Source_Audit_english.svg) | [PNG](Fig1_Source_Audit_chinese.png) · [SVG](Fig1_Source_Audit_chinese.svg) |
| 2. FreeSolv regression errors and frozen-pool search | [PNG](Fig2_Graph_Benchmark_english.png) · [SVG](Fig2_Graph_Benchmark_english.svg) | [PNG](Fig2_Graph_Benchmark_chinese.png) · [SVG](Fig2_Graph_Benchmark_chinese.svg) |
| 3. Sampled target geometry and MMFF-weighted SASA | [PNG](Fig3_Conformer_Ensemble_english.png) · [SVG](Fig3_Conformer_Ensemble_english.svg) | [PNG](Fig3_Conformer_Ensemble_chinese.png) · [SVG](Fig3_Conformer_Ensemble_chinese.svg) |
| 4. SSA–CTMC numerical comparison | [PNG](Fig4_Kinetic_Validation_english.png) · [SVG](Fig4_Kinetic_Validation_english.svg) | [PNG](Fig4_Kinetic_Validation_chinese.png) · [SVG](Fig4_Kinetic_Validation_chinese.svg) |

Figure 1 uses source audit records, including deliberately preserved unsupported source values. Figure 2 reports the 51-molecule test set and eight seeds per acquisition method; search curves show means, not confidence bands. Figure 3 reads saved SDF coordinates and pooled 298.15 K descriptors; connectivity omits bond order. Figure 4a shows conditional 95% Student-t intervals for 32-seed means. Figure 4b displays only the saved 0–0.2 s analytic transient, with the stochastic trace downsampled for display (every seventh event record); complete trajectory data remain archived. The dotted right boundary is the burn-in endpoint.

All figures describe a source audit, a hydration benchmark or assumed molecular/kinetic models. They do not establish oxidation potentials, solution populations, reaction mechanisms or electrode-scale performance.

[Figure generator](../../scripts/plot_reviewed.py) · [Source data and SHA256 manifest](../../results/figure_manifest.json) · [PNG and independently rendered SVG inspection](../../results/figure_qa.json)
