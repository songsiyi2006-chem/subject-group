# QuantumEqui-NEB figures / 中英文图表

Four figure pairs are generated directly from saved calculation tables. Each language has four 2400 × 1380 PNGs at 300 DPI and four editable-text SVGs with no embedded raster images. SVG fonts are not embedded; a suitable Chinese font is needed on other systems.

| Figure / 图 | English PNG / SVG | 中文 PNG / SVG | Evidence |
|---|---|---|---|
| 1. Quantum reference | [PNG](Fig1_Quantum_Reference_english.png) · [SVG](Fig1_Quantum_Reference_english.svg) | [PNG](Fig1_Quantum_Reference_chinese.png) · [SVG](Fig1_Quantum_Reference_chinese.svg) | Source overlap clipping and 42 executed H2 RHF/STO-3G energies |
| 2. Force audit | [PNG](Fig2_Force_Audit_english.png) · [SVG](Fig2_Force_Audit_english.svg) | [PNG](Fig2_Force_Audit_chinese.png) · [SVG](Fig2_Force_Audit_chinese.svg) | Finite differences versus autograd on the captured untrained function |
| 3. Analytic CI-NEB | [PNG](Fig3_NEB_Validation_english.png) · [SVG](Fig3_NEB_Validation_english.svg) | [PNG](Fig3_NEB_Validation_chinese.png) · [SVG](Fig3_NEB_Validation_chinese.svg) | Two dimensionless designed potentials, seed 11, force tolerance 1e-5 |
| 4. Thermochemistry | [PNG](Fig4_Thermochemistry_english.png) · [SVG](Fig4_Thermochemistry_english.svg) | [PNG](Fig4_Thermochemistry_chinese.png) · [SVG](Fig4_Thermochemistry_chinese.svg) | Nonstationary source spectrum and separate ideal-gas analytic controls |

All eight PNGs and independent rasterizations of all eight SVGs were visually inspected for labels, legends, clipping and units. This is artifact QA, not experimental or chemical validation. [Exact data/figure hashes](../../results/figure_manifest.json) and [QA record](../../results/figure_qa.json) accompany the [plotting code](../../scripts/plot_reviewed.py).

原始脚本生成的 [综合面板](../../results/original/figures_quantum/Fig_QuantumEquiNEB_Comprehensive_Panel.png) 单独存档，其未校准势能、错误频率换算和硬编码自由能偏移不作为修订图表的科学结论。
