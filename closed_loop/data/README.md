# Feedback records and evidence roles

`mock_feedback.csv` is extracted from the supplied script. Every yield is **mock**, despite the source's HPLC comment. `reaction_id=unspecified_demo` deliberately does not assert that these labels describe the Section 2 indole reaction. This file is an example schema, not an experimental dataset.

Numerical fields: current_density_mA_cm2, electrolyte_concentration_M, temperature_K, yield_pct and prior_prediction_pct. Percentage fields use 0–100, not fractions; temperature is kelvin. `run_id` must be unique; independent replicates at the same conditions have distinct IDs. `reaction_id` identifies one common transformation.

For real data, prepare a separate CSV and set `evidence_role=experimental`. All rows in a run must share one reaction identity, including canonicalized substrate_smiles, partner_smiles and product_smiles. Supply assay=HPLC, operator and recorded_at; the time field must contain the actual laboratory record time (its format is not yet enforced by this checker). Add relative `assay_file`, `calibration_file`, `electrolysis_file` paths and each corresponding `_sha256` field. Files must resolve inside the CSV directory and have the specified SHA256. The contents, calibration response factors, integration and product assignment still need human review; a file hash is not measurement validation.

```powershell
python closed_loop/scripts/audit_closed_loop.py --feedback path/to/records.csv --role experimental --output work/real-feedback-review
```

Explicit experimental mode rejects mock records. Missing raw measurements cannot be replaced by invented data. Preserve failed runs and missing outcomes in a separate audit ledger; this minimal numerical loader does not fit records with missing yields or treat missingness as zero. The demonstration model supports three input variables only; richer reaction identity, batch effects, heteroscedastic noise and search-domain changes require deliberate model changes. Experimental mode checks links and refits the same demonstration policy; it does not make it a validated recommendation engine.

The audit writes an output snapshot and feedback-file hash. It neither edits raw assay files nor launches an experiment. Keep each round under a new versioned directory, review the model against an independently measured baseline, and rebuild reports deliberately after accepting the data. No wet-lab records were available for this release.
