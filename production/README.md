# Quantitative electrosynthesis pipeline: execution and evidence audit

**2026-09-27.** The supplied “production-grade” pipeline is an auditable computational prototype. Original execution failed at the RDKit SASA API; the separately saved repair ran successfully. Production readiness and experimental validity remain unestablished.

- [English report](reports/production_technical_report_english.md)
- [中文报告](reports/production_technical_report_chinese.md)
- [Section-aligned bilingual report / 双语合并版](reports/production_technical_report_bilingual.md)
- [Original specification](source/specification.md), [explicit repair](source/repair.patch), [repaired output](results/repaired/production_benchmark_results.json), [audit summary](results/audit/audit_summary.json)

Each report follows the new specification's four sections: scientific context, four-task quantitative breakdown, proposed wet-lab/instrumentation program, and undergraduate milestone Gantt. Separate language editions preserve the earlier user preference; the combined edition additionally supplies the filename requested in the new brief.

| Added work | Executed scope | Finding and evidence boundary |
|---|---|---|
| Bayesian optimization | 30 paired seeds × 4 strategies × 14 evaluations; 1,680 observations | Normalized physical and categorical variants have higher mean regret than random search here; artificial response surface |
| Molecular geometry | 32 ETKDG/MMFF conformers; 6 probe/algorithm combinations | Heuristic favors pendant phenyl, not C3; no regioselectivity validation |
| Target electronic structure | 3 GFN2-xTB/ALPB(MeCN) jobs | Neutral optimization plus fixed-nuclei ±1 single points; no DFT or measured Eox |
| SAC transfer | 16 artificial labels; 2 group schemes × 4 folds × 3 models | Training fit is not transferable chemical knowledge |
| Flow audit | Analytical/DOP853 comparison; current integration; 48 potential/flow and 20 parameter cases | Higher flow improves mass transfer but reduces residence time; STY is reactant disappearance |

**关键修正：** C2、C5 已被苯基／甲氧基取代，没有 C–H；Gasteiger 电荷不是自由基阳离子自旋密度；SAC 的“势垒”由公式指定；D 是固定过电位下的线性轴向 ODE，不是已解析电位、速度和浓度场的 PDE。原注释的科学主张在报告中逐项核验。

## Files and provenance

```text
source/                  Supplied specification, exact extracted code, hashes, repair diff
run_production_pipeline.py
scripts/                 Execution, additional calculations, plotting, reports, validation
results/original_attempt/ Original error and execution record
results/repaired/         Four-task JSON, unrounded intermediate tables, saved MMFF geometry
results/audit/            Per-run tables, 32-conformer SDF, atom mapping and audit summary
results/target_xtb/       Three actual xTB inputs, coordinates, charges, JSON and logs
reports/                 Separate and combined reports, templates, English/Chinese figures
```

The source extraction normalizes line endings and retains scientific code exactly. [Source hashes](source/source_record.json) distinguish the attachment from the extracted script. The SASA radius convention was changed explicitly to positive RDKit periodic-table radii; this is a model convention change, not just capitalization repair. Public logs replace local runtime/directory paths while preserving scientific values. Target atom indices are zero-based RDKit indices; they are not the one-based indices in the earlier six-molecule study.

## Reproduce

Use the [recorded environment](../provenance/requirements_recorded.txt) and an existing xTB installation. No dependencies were installed for this run. On Windows, activate the chemistry environment so its `Library/bin` is on PATH. If needed, set `XTB_EXE` to the existing xTB executable. Use a separate clone for reproduction: these scripts refresh generated files under `production/`.

Run from the repository root, with single-thread environment variables set:

```powershell
$env:OMP_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
python production/scripts/prepare_and_execute.py
python production/scripts/audit_pipeline.py --seeds 30
python production/scripts/plot_production.py --zh-font C:/Windows/Fonts/msyh.ttc
python production/scripts/build_reports.py
python production/scripts/validate_production.py
python -m unittest discover -s tests -v
python scripts/validate_release.py
```

Use an installed CJK font path on other systems. The committed PNG/SVG figures are already available; SVG text uses embedded outlines. The numerical audit took about 34 seconds on the recorded host; runtime is machine-dependent. Each xTB job is serial with a 180-second timeout. Source/repaired executions have a 300-second timeout. The original failure is expected and preserved.

`validate_production.py` checks report architecture, numerical table agreement, combined-edition coverage, execution records, data dimensions, local links and JSON integrity. The 14 additional scientific tests in [test_production.py](../tests/test_production.py) check source integrity, structure mapping, SASA bookkeeping, paired BO budgets, group holdouts and independent flow/current identities. Root discovery also runs the earlier 15 tests. These checks inspect saved xTB records and do not rerun quantum jobs.

The repository [SHA-256 manifest](../provenance/SHA256SUMS.txt) covers the release. Rerunning calculations can change timing records or numerical outputs and will legitimately invalidate hashes. Review changes before deliberately updating the manifest with `scripts/update_manifest.py`. Numerical checks do not certify chemical predictions. The report's instrument list and protocols are proposed measurements, not a verified local inventory or an executed group SOP.
