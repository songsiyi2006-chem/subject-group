# Analytical parsing, scientific graphics and grant toolkit

2026-09-28 · 执行记录与证据审计；没有新增湿实验、DFT 或正式申报。

- [中文完整技术报告](reports/deployment_toolkit_report_chinese.md)
- [Complete English technical report](reports/deployment_toolkit_report_english.md)
- [完整中文大创申报书内容稿](proposals/National_Undergraduate_Grant_Proposal_PanTang_Lab.md)
- [Complete English proposal](proposals/National_Undergraduate_Grant_Proposal_PanTang_Lab_English.md)
- [Four figures, two languages, PNG and editable SVG / 图表索引](reports/figures/README.md)
- [Numerical audit](results/audit/audit_summary.json), [publication checks](results/publication_validation.json), [figure QA](results/figure_qa.json)

The supplied program completed after its outer Python block was correctly extracted; its suggested regex truncates an internal Markdown string. The original logic and outputs remain archived in [source](source/specification.md) and [results/original](results/original/README.md). Their unsupported wet-lab/DFT/equipment labels are not endorsed.

The source HPLC arithmetic gives **1292.00%**, versus the rounded NMR example **91.392%**, so its concordance claim is false. Mass-based NMR recomputation gives **91.3252869%**. A strict dominance calculation finds **23** synthetic Pareto points, while the original threshold finds zero. Six substrate pairs are mock labels; the energy diagram contains five given constants, not DFT results.

## Additional executed calculations

| Calculation | Scale | Evidence boundary |
|---|---|---|
| Analytical input-error propagation | 100,000 draws | Assumed errors, not validated measurement uncertainty |
| qNMR purity and sample-fraction scenarios | 12 | Source integrals, no real spectrum |
| qNMR relaxation sensitivity | 6 | Assumed T1 values and ideal repeated 90-degree pulses |
| Synthetic calibration | 21 rows | Seven known levels, three observations each |
| Known-shape chromatogram recovery | 80 cases | Four separations, two widths, ten seeds; no vendor parsing |
| Pareto membership sensitivity | 2,000 draws over 80 points | Assumed objective noise, no experimental confidence probability |
| Mock scope influence diagnostics | 6 deletions | Descriptive metrics, not trained-model cross-validation |
| Proposal budget | 7 items, CNY 10,000 | Unapproved planning prices |

Peak-width mismatch raises mean absolute recovered-area error from 0.0274% to 14.6143%. The negative result remains in the report and figures. All corrected figures visibly identify their synthetic/mock/illustrative nature. The proposal replaces unsupported hardware access and guaranteed publication/award claims with a twelve-month, supervised research plan and explicit unresolved application fields.

## Reproduction

Use the existing scientific Python environment with NumPy, SciPy, pandas, Matplotlib, Pillow and RDKit. The source script imports RDKit; no new quantum calculation is performed. The archived environment is listed in [audit_summary.json](results/audit/audit_summary.json). Reuse installed packages instead of indiscriminate installation. On the recorded Windows environment the matching numerical DLL directory must be on PATH. Use one numerical thread for the bounded run.

```bash
python toolkit/scripts/run_original.py
python toolkit/scripts/audit_toolkit.py
python toolkit/scripts/plot_reviewed.py
python toolkit/scripts/validate_toolkit.py
python -m unittest discover -s tests -p test_toolkit.py
```

Run from the repository root. The first command regenerates original outputs and timestamps; it is optional when using the byte-preserved archive. Import custom **pre-integrated** HPLC or NMR CSV with `audit_toolkit.py --input-csv ... --assay hplc|nmr --role source_example|synthetic|experimental --output ...`; use a fresh output folder and the [data dictionary](data/README.md). The parser is not an instrument driver, vendor format decoder, NMR FID processor or peak-assignment engine.

For a reviewed repository release, stage intentional changes, run the module validators, then update `scripts/update_manifest.py` and run `scripts/validate_release.py`. Tests and hashes check software/document integrity, not scientific discovery. No grant application is submitted by these scripts.
