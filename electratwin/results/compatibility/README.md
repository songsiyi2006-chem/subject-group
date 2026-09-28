# Compatibility-only source-model output / 仅兼容修复的源模型输出

These files are the actual outputs of the submitted model after the sole replacement `np.trapz` → `np.trapezoid`. The original figures, JSON payload and console claims are retained for audit, including unsupported industrial and hardware language. They are not reviewed conclusions or wet-lab evidence.

`captured_model_arrays.json` additionally records the source's complete single-case arrays and twelve campaign rows after its main block returns. `../source_audit.json` reconstructs FE before clipping, quantifies the current/material-flux mismatch, and separates a fixed grid from Bayesian optimization. `execution.json` confirms completion but does not validate the physical model. The driver is an in-memory emulator; the recorded voltage and instrument identity are synthetic.
