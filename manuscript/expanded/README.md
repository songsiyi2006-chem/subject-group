# Expanded bilingual computational chemistry manuscript

This is the expanded research version assembled from six fixed public repository histories. It preserves the earlier manuscript in the parent directory.

- [English manuscript PDF](documents/manuscript_english.pdf) and [editable Word document](documents/manuscript_english.docx)
- [中文论文 PDF](documents/manuscript_chinese.pdf)与[可编辑 Word 文档](documents/manuscript_chinese.docx)
- [English source](manuscript_english.md) and [中文源稿](manuscript_chinese.md)
- [Verified references](references_verified.json), [numbered bibliography](bibliography.md) and [BibTeX](references.bib)
- [Repository inventory](../history/repository_inventory.json), [learning audit](../history/learning/audit.json), [catalysis audit](../history/catalysis/audit.json) and [pharmacology audit](../history/pharmacology/audit.json)
- [Final document checks](results/document_qa.json), [visual review](results/visual_review.json), [numeric translation parity](results/numeric_parity.json) and [release validation](results/validation.json)

The central evidence remains the matched-reference H₂ error analysis and the separately parameterized Morse nuclear model. Historical molecular, catalytic, pharmaceutical, transport and analytical records add quantitative comparisons and negative controls. Public experimental traces belong to the cited original researchers. Synthetic objectives, assigned-rate models, approximate scores, cached evaluations and actual electronic/atomistic calculations are labeled separately; no aggregate count turns these into independent experiments.

## Reproduction

Run commands from the repository root. The expansion's arithmetic and plots do not launch new quantum calculations, training, molecular dynamics, docking or experiments. NumPy/Matplotlib are needed for plots; document production also needs python-docx/Pillow and a PDF exporter.

```sh
python manuscript/expanded/scripts/analyze_multiscale.py
python manuscript/history/learning/plotting.py
python manuscript/history/catalysis/plotting.py
python manuscript/history/pharmacology/plotting.py
python manuscript/expanded/scripts/build_reference_registry.py
python manuscript/expanded/scripts/assemble_manuscripts.py
python manuscript/expanded/scripts/finalize_references.py
python manuscript/expanded/scripts/check_numeric_parity.py
python manuscript/expanded/scripts/validate_expansion.py --data-only
```

Plot renderings may acquire new renderer metadata and hashes on regeneration; the scientific input bytes and arithmetic are the reproducible objects. Rebuild documents after changing figures or text, refresh machine and visual reviews, and then run the full validator. A historical visual-review record does not certify a newly rendered PDF. Bibliographic records were verified during preparation; rerunning the registry builder records that inspection and does not perform a fresh network check.

The Word files use Letter pages, ordinary journal-style typography and high-resolution display-equation images. Equation source is retained under `assets/equations/`; the display images are not native editable Word equations. PDF export used the installed Microsoft Word application after the bundled LibreOffice renderer was found unavailable. All final PDF pages were rendered and visually inspected by AI agents; this is not human scientific peer review.

## Scope and rights

Commit histories were reviewed for lineage and important numerical results. Selected raw records, native-log excerpts and derived tables were checked; this is not an exhaustive rerun of every archived calculation. Original-source hashes and excerpt hashes are kept distinct. The original author-provided identity is Siyi Song / 宋思毅, Guangxi Normal University / 广西师范大学.

Upstream licenses and third-party data rights remain applicable. Bounded historical extracts are included for the repository owner's requested analysis; an MIT notice in one repository does not relicense publisher data or another repository. No journal-quartile status, experimental validation, priority claim or journal acceptance is certified by the software checks.
