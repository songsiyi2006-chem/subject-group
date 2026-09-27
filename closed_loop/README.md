# Closed-loop AI4S: stoichiometry, DFT inputs and feedback audit

Updated 2026-09-28. All three modules in the supplied brief were executed; reports follow its five-section publication structure.

- [完整中文报告](reports/closed_loop_experimental_report_chinese.md)
- [Complete English report](reports/closed_loop_experimental_report_english.md)
- [完整章节对齐双语版 / Full bilingual report](reports/closed_loop_experimental_report_bilingual.md)
- [Original result](results/original/closed_loop_results.json), [source and extraction provenance](source/source_record.json), [execution log](results/original/execution.json)
- [Independent audit summary](results/audit/audit_summary.json), [full model results](results/audit/active_learning.json), [reviewed DFT inputs](inputs/reviewed/README.md)
- [Feedback schema and evidence rules](data/README.md), [publication validation](results/publication_validation.json)

**证据边界：** 原始反馈是六条模拟标签，没有湿实验或 DFT 结果。2.88 个百分点是给定预测与模拟标签的差异，并非交叉验证损失。新增留一验证中，标准化 GP 的 MAE 为 15.26，均值基线为 13.98；没有证明 GP 更优。SOP 缺少已定义产物、偶联当量及已核验的淬灭/分析方法，仍是待实验指导人员审阅的设计草案。

**Evidence boundary:** The updated EI ranking is conditional on mock data. The source calls UCB “EI,” uses indoline instead of the wet-lab target, and labels a CPCM ORCA input as SMD. Source outputs are preserved; reviewed inputs and analyses correct those discrepancies. File generation is not engine validation or chemical proof.

## Executed additions

| Addition | Extent | Interpretation |
|---|---|---|
| Stoichiometry sensitivity | 81 combinations | Scale, volume, area and current arithmetic; no scale-up validation |
| Partner loading | 8 scenarios | Thioanisole and thiophenol remain chemically different options |
| Charge/FE balance | 4 scenarios + current integration | Conditional on a two-electron product |
| Leave-one-out evaluation | 24 predictions over four methods | Training-only preprocessing; six mock rows |
| Candidate analysis | 270 predictions, 90 points for three GPs | 89 unseen candidates; EI and UCB retained separately |
| Acquisition sensitivity | 9 noise/xi scenarios | Different assumptions select different candidates |
| EI integration check | 200,000 normal draws | Independent check of the acquisition formula |
| DFT inputs | 4 Gaussian + 8 ORCA files | Two molecules, neutral/cation states; no DFT executed |
| Geometry reuse | 2 previously completed xTB geometries | File hashes, molecular identity, atom order and spin parity checked |

## Reproduction

Use the existing Python environment. Recorded package versions are in `results/audit/audit_summary.json`; no dependencies were installed. Set numerical thread counts to one. On Windows, the appropriate Conda environment's `Library/bin` must be on PATH for numerical DLLs.

```powershell
$env:OMP_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
python closed_loop/scripts/run_original.py
python closed_loop/scripts/audit_closed_loop.py --output work/closed-loop-audit
python closed_loop/scripts/generate_reviewed_inputs.py --output work/closed-loop-inputs
python -m unittest discover -s tests -v
python closed_loop/scripts/validate_closed_loop.py
```

The original runner refreshes its archived result and timing. To preserve release bytes, reproduce in another checkout. The audit with an explicit output directory preserves saved audit results, but regenerates the deterministic default mock CSV unless `--feedback` is supplied. The DFT generator reads prior geometries from the full repository and its SHA256 manifest; copying only `closed_loop` is insufficient. The generator does not invoke a quantum engine or scheduler. Choose one engine protocol, inspect the optimization, then follow the dependent single-point step in the input guide.

To regenerate publication assets from saved release results:

```powershell
python closed_loop/scripts/plot_closed_loop.py
python closed_loop/scripts/build_reports.py
python closed_loop/scripts/validate_closed_loop.py
```

The Chinese plots use Microsoft YaHei if available; a platform without a CJK-capable font needs an equivalent font configured before rebuilding. Validation checks file and numerical integrity, not experimental validity. The root [release validator](../scripts/validate_release.py) and [SHA256 manifest](../provenance/SHA256SUMS.txt) also cover this extension.
