# Historical evidence snapshot

This directory contains selected exact source files from the two commits recorded
in `../audit.json`, plus explicitly named `derived_*.csv` tables and post-hoc
analysis scripts. Files named `phase*` are frozen historical source outputs, not
new calculations. The original two MIT notices are retained. Third-party source
article/workbook rights remain unchanged; selected published numerical summaries
carry their DOI and workbook hash, and no publisher figure or article is copied.

Run from any directory:

```text
python recompute_history.py
python ../plotting.py
python finalize_history.py
```

The first and third scripts use the Python standard library. Plotting requires
NumPy, Matplotlib, DejaVu Sans and Microsoft YaHei; the delivered renderer uses
the existing Windows font at `C:/Windows/Fonts/msyh.ttc`. A conda Python must be
launched with its native `Library/bin` DLL directory on PATH. No dependency is
downloaded and no original scientific program is executed. The finalizer checks
all six tables' bilingual numeric parity and all figure hashes/dimensions, then
refreshes `../audit.json`. It does not rerun training, xTB, DFT or Monte Carlo.

To recreate the exact source snapshots, `collect_history.py` accepts a directory
containing read-only clones named `aqueous-solubility-ml-benchmark` and
`ai4chem-complex-scaffolds-benchmark` at the required commits. It uses `git show`
to preserve committed bytes and refuses mismatched HEADs. Run the three commands
above afterward to regenerate post-hoc and publication metadata.

Audit depth is intentionally explicit: the complete 54-commit subject history
and latest acceptance scope for all 31 phases were inspected, while quantitative
recomputation concentrates on conformers, saved VMC block means, published trace
summaries, and matched charged-molecule energies. The eight DFT energy differences
are independently reconstructed from individual result JSON files; 72 xTB state
differences are reconstructed from the 144-job ledger. This is not an independent
reparse of every raw solver log or a complete trajectory audit. Historical running
statuses are not evidence that a process is currently active.

The solubility repository has no committed per-molecule predictions, dataset
snapshot, split IDs or fitted weights. Its reported performance is preserved in
the source README but is deliberately excluded from the recomputed performance
tables. Its single commit is also the root of the complex-scaffolds repository;
five inherited code files and four figures are byte-identical.
