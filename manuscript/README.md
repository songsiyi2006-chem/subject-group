# Integrated computational chemistry manuscript / 整合计算化学论文

**Reference bias and observable dependent convergence in computational chemistry**  
**计算化学中的参考偏差与观测量依赖的数值收敛**  
Siyi Song / 宋思毅 — Guangxi Normal University / 广西师范大学  
Research manuscript draft, 28 September 2026. Sole author as supplied by the user; no final author approval or journal submission is asserted.

## Read the paper / 阅读论文

| Edition | PDF | Editable Word | Linked text |
|---|---|---|---|
| English main manuscript | [PDF](documents/manuscript_english.pdf) | [DOCX](documents/manuscript_english.docx) | [Markdown](rendered/manuscript_english.md) |
| 中文正文 | [PDF](documents/manuscript_chinese.pdf) | [DOCX](documents/manuscript_chinese.docx) | [Markdown](rendered/manuscript_chinese.md) |
| English supporting information | [PDF](documents/supporting_information_english.pdf) | [DOCX](documents/supporting_information_english.docx) | [Markdown](rendered/supporting_information_english.md) |
| 中文补充信息 | [PDF](documents/supporting_information_chinese.pdf) | [DOCX](documents/supporting_information_chinese.docx) | [Markdown](rendered/supporting_information_chinese.md) |

The main text develops two controlled calculation chains. Supporting Information S1–S13 integrates the earlier five-topic benchmark, molecular calculations, production/research audits, closed-loop and analytical tools, ElectraTwin, ElectroGraph, SynthaPore and QuantumEqui. All earlier reports, source scripts and raw results remain available. This release adds **post hoc arithmetic**, not new quantum energy calls, neural training or eigenproblems.

本文以两条受控计算链组织主线；补充信息 S1–S13 覆盖此前所有模块。新增分析包括 12 组有符号 MSE 归因、6 组配对误差改善比较、3 个固定区间 Richardson 估计，以及 7 个残差与观测量误差对照。保留失败、负面结果和模型适用边界，不将合成标签、模型势或数值检查写成实验验证。

## Findings and contribution boundaries / 结论与创新边界

- At eight matched H₂ test geometries, gradient supervision improves RHF-label energy RMSE by **60.24–73.65%**, but improves total same-basis FCI RMSE by only **0.00838–0.02168%**. Signed cross terms distinguish cancellation from a change to the electronic reference.
- A matrix residual near machine precision can coexist with a missing shallow state or inaccurate binding. The finest saved single-grid binding error remains **1.564%**. New fixed-domain Richardson arithmetic gives **−0.02688%**, conditional on the second-order assumption; it is not another eigensolve or experimental spectrum.
- The neural RHF/FCI comparison and the FCI-parameterized Morse nuclear calculation are **separate chains**. No trained-neural-potential-to-vibrational-spectrum calculation has been executed.
- The identities, reference bias and finite-difference principles are established. The bounded contribution is their quantitative linkage and inspectable counterexamples. Chemical-space generalization, a new electrocatalytic mechanism and readiness for a top journal have not been demonstrated. See the [bilingual editorial assessment](editorial_assessment_bilingual.md) for prior work, falsifiable contribution statements and the missing evidence.

## References and evidence / 文献与证据

[Verified reference registry](references_verified.json) records 25 published works with source URLs and the fields actually checked; 23 have DOI records, and two use official proceedings/journal records without an invented DOI. Four papers by Haitao Tang and collaborators are identified through authors and affiliations. Nineteen distinct references are used in the main text and SI; six additional relevant records support the broader editorial assessment. See [complete bibliography](bibliography.md), [BibTeX](references.bib) and [citation audit](results/citation_audit.json).

**“中科院一区”与“JCR Q1”不混用。** 当前未取得可核查的对应年份权威分区表，因此未把登记文献统一标成“一区”；出版身份核验不等同于分区核验。唐海涛及合作者的实验论文提供研究背景，不能证明本文的小体系或合成模型完成了对应反应的验证。

The [evidence catalog](evidence_catalog.json) links 206 structured selectors to 163 hashed files across 12 modules. The [evidence audit](evidence_audit.md) separates direct calculation, public measured labels, analytic controls, synthetic targets and unexecuted proposals. Scientific source data are frozen at commit `1ad05c243155a9f64b18e0d58c665424fde919a0`.

## Reproduction / 复现

Use an existing configured scientific Python environment from the repository root:

```sh
python manuscript/scripts/analyze_error_chain.py
python manuscript/scripts/plot_manuscript.py
python manuscript/scripts/finalize_references.py
python manuscript/scripts/validate_manuscript.py --require-documents
```

The first command reads saved data and regenerates only arithmetic outputs. The plot script uses NumPy and Matplotlib; a CJK font is needed for Chinese. Validation independently recomputes numerical relationships and checks records, hashes and bilingual consistency. It is not an independent experimental replication or human peer review.

For document generation, use Python with `python-docx`, Pillow, Matplotlib, `pypdf` and `pypdfium2`, and an existing Word or LibreOffice installation. The executed equation renderer is the existing Matplotlib mathtext/STIX implementation; the attempted MathJax package download was unavailable and is not a build dependency:

```sh
python manuscript/scripts/build_documents.py --prepare-equations
python manuscript/scripts/render_equations.py
python manuscript/scripts/build_documents.py
```

On Windows, `manuscript/scripts/export_word_pdf.ps1` uses Word to export the four named files. Display equations are high-resolution rendered images in Word, with editable source notation provided; the remaining text, tables and references are editable. PDF and Word pagination can vary across renderer/font versions. Run `inspect_documents.py --render-dir <scratch-directory>` and inspect every rendered page after rebuilding. The final main papers have 14 English and 12 Chinese pages; the SI has 10 English and 9 Chinese pages. All 45 pages were rendered for review. The [document checks](results/document_qa.json), [English main visual review](results/main_english_visual_qa.json), [Chinese main visual review](results/main_chinese_visual_qa.json), [SI visual review](results/si_visual_qa.json), [figure QA](results/figure_visual_qa.json) and [scientific validation](results/validation.json) record their respective scope.
