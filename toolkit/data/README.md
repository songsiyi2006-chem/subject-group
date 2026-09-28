# Pre-integrated analytical input contract / 预积分数据约定

These CSVs are **source examples**, not laboratory records. Keep `evidence_role=source_example` for the supplied values. Unit conventions are deliberate: area is absolute, concentrations are mM, reaction volumes mL, standard masses mg, amounts mmol and molecular weights g/mol. Unsupported area ratios and ambiguous unit conversions are rejected. UTF-8 CSV (with optional BOM) is accepted.

## Common identity and evidence fields

| Field | Requirement |
|---|---|
| `sample_id` | Nonempty, unique within the input file |
| `evidence_role` | Must equal command-line role: `source_example`, `synthetic` or `experimental` |
| `product_identity` | Required for experimental role; a chemical identifier checked by the analyst |
| `analyst`, `recorded_at` | Required for experimental role; actual identity/time, not automatically authenticated |
| `raw_file`, `raw_sha256` | Experimental role: local file under the CSV directory and matching SHA-256 |

File existence and hashes prove identity, not chemical validity. The program does not inspect vendor metadata or authenticate an analyst. Scientific validity always remains `measurement_validated=false`; experimental-role acceptance means provenance fields passed, not the measurement passed validation.

## HPLC columns

[source_integrated_area.csv](source_integrated_area.csv) illustrates the full schema. `area_mode` must be `absolute`, `concentration_unit` must be `mM`, and `volume_unit` must be `mL`. Numerical fields are `area`, `slope_area_per_mM`, `intercept`, `dilution_factor`, `volume_mL`, `substrate_mmol`. Slope, dilution, volume and substrate amount must be positive finite values. The initial-substrate convention is mmol; one-to-one substrate/product stoichiometry is assumed.

`calibration_min_mM` and `calibration_max_mM` define the independently documented applicability range. Both may be blank for a source/synthetic example; a partially filled range is rejected. Both are required for experimental rows, with a nonnegative lower limit and larger finite upper limit. `calibration_file` and `calibration_sha256` must identify the local supporting calibration evidence for experimental rows. The source's real calibration range is unknown; the simulated 0–2 mM calibration elsewhere is not substituted for it.

The parser reports concentrations, product mmol and yield percent. It flags negative baseline-corrected concentration, out-of-range calibration and yield outside 0–100% without clipping. A numeric-range pass alone cannot establish identity, specificity, precision or recovery.

## NMR columns

[source_nmr_integrals.csv](source_nmr_integrals.csv) illustrates the full schema. `mass_unit=mg`, `amount_unit=mmol`. Fields `product_integral`, `standard_integral`, `product_H`, `standard_H`, `standard_mass_mg`, `standard_MW`, `standard_purity`, `substrate_mmol`, `sample_fraction` implement the proton-normalized amount equation. Product integral may be zero; other denominators and masses are positive finite. Purity and sampling fraction must lie in (0,1]. Proton counts represent the assigned signal's number of nuclei, not an arbitrary scaling fitted to obtain a yield.

Experimental rows additionally require `product_peak_assignment`, `standard_identity`, `relaxation_record`, and `standard_certificate_file`/`standard_certificate_sha256`. The latter file must reside under the CSV directory and match its hash. A relaxation record should describe actual acquisition conditions and their justification; the parser only checks presence, not physical adequacy. Blank example fields must not be filled with invented metadata.

No FID transform, phase/baseline correction, peak detection, multiplet assignment or purity certification is implemented. Retain unprocessed instrument data outside a public repository when it contains private information; publish only with appropriate authorization and provenance.

## Example imports

```bash
python toolkit/scripts/audit_toolkit.py --input-csv toolkit/data/source_integrated_area.csv --assay hplc --role source_example --output work/hplc-example
python toolkit/scripts/audit_toolkit.py --input-csv toolkit/data/source_nmr_integrals.csv --assay nmr --role source_example --output work/nmr-example
```

The HPLC example intentionally fails the one-to-one material balance. The NMR example retains unverified peak assignment and acquisition; absence of a numerical range flag is not validation.
