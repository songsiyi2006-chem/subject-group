# FreeSolv attribution and data transformations

The experimental hydration free energies are from David L. Mobley and J. Peter
Guthrie, *FreeSolv: a database of experimental and calculated hydration free
energies, with input files*, J. Comput. Aided Mol. Des. (2014),
[doi:10.1007/s10822-014-9747-x](https://doi.org/10.1007/s10822-014-9747-x).

The [official MobleyLab repository](https://github.com/MobleyLab/FreeSolv)
identifies its data license as CC-BY 4.0, with the stated qualification that
original third-party restrictions may apply. The complete upstream license is
preserved in [FreeSolv_LICENSE.txt](FreeSolv_LICENSE.txt). Code in the upstream
repository has a separate MIT license; no upstream code is copied here.

[official_database.txt](official_database.txt) is the official v0.52 text snapshot
retrieved on 2026-09-28. Its per-compound experimental references, uncertainty
annotations, and notes are retained. A local 642-row SAMPL.csv snapshot was reused
and canonical structures and experimental labels were checked against this
official snapshot. Any mismatches are explicitly excluded and recorded in
[provenance.json](provenance.json). Calculated GAFF values are not model targets
or predictor inputs.

[selected_freesolv.csv](selected_freesolv.csv) contains 256 unique canonical
structures selected by sorted structure SHA256, without using label values to
rank them. We add structure/group/split metadata and preserve experimental
labels and original references. The split groups ring molecules by Murcko
scaffold. Acyclic molecules use the full canonical structure as their group;
this prevents exact duplicates but does not establish acyclic scaffold novelty.
Provenance records source and derivative hashes. Labels are hydration free
energies in kcal/mol, never oxidation potentials or atom-reactivity labels.
