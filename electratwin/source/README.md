# Submitted source and compatibility record / 原始输入与兼容记录

`specification.md` preserves the attachment bytes. `electratwin_core.py` is the first outer Python block with LF line endings and a final newline. Its scientific logic is unchanged. Authorship, laboratory alignment, production readiness and industrial claims are copied statements, not verified affiliations or achievements.

The actual original run failed because NumPy 2.4 removed `np.trapz`. `electratwin_core_numpy_compat.py` changes exactly one call name to `np.trapezoid`; `numpy_compatibility.patch` records that change. This repair permits execution and does not repair physics, stoichiometry, metrics or the planner. SHA-256 hashes are in `source_record.json`.

原程序错误与最小修复后的计算分别保存在 `../results/original/` 与 `../results/compatibility/`。原程序输出中的“工业级”“自主优化”“实测电压”“Connected”均属于原代码输出，不能当作验证结论。SCPI 类仅改变内存变量，不打开网络、GPIB、串口或 VISA 连接。源码给出的 IP 和仪器身份仅作为原文记录；未尝试连接。

`../scripts/run_source.py` executes the source in bounded subprocesses with single-thread numerical settings. After the compatibility program returns, the wrapper saves its arrays and campaign history without changing its source logic. `--audit-only` recomputes the source audit from saved arrays and performs small in-memory planner/driver probes; those probes use a stub solver and do not count as PDE runs or experiments.
