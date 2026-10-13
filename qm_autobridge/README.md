# QM AutoBridge / 量化日志与三维学习桥接工具

The implementation is a single, standalone file: [qm_autobridge.py](../qm_autobridge.py). Its only third-party runtime dependencies are NumPy, SciPy, Matplotlib and PyTorch. It does not import cclib, PyG, DGL, a quantum engine, a scheduler or hardware drivers.

本次交付提供五个可运行模块及软件验收数据。默认执行的是**合成日志集成测试**，没有开展新的 DFT 或超算任务。它不是经大规模真实日志语料验证的通用解析器，也不是已经训练好的氧化电位预测模型。

## Run / 运行

Use Python 3.10+ with the four dependencies installed. The repository's existing scientific environment can be reused; the standalone script does not need the rest of the repository.

```sh
python qm_autobridge.py --output work/qm-demo
python qm_autobridge.py --validate-saved qm_autobridge/results
python qm_autobridge.py --parse your_single_job.out --output work/parsed-log
python -m unittest discover -s tests -p test_qm_autobridge.py -v
```

The output directory must be new; existing files are not overwritten. `--seed` controls the untrained model initialization. The sandbox writes 14 Gaussian/ORCA fixtures for substrate, radical cation, R, PC, TS, INT and P, then verifies format parity. Geometry is a fixed RDKit ETKDG embedding (seed 42) of `COc1ccc2[nH]c(-c3ccccc3)cc2c1`, C15H13NO. RDKit is not required at runtime. Energies, charge populations, frequencies and displacement vectors are invented fixtures, not quantum results; the mode vectors are not eigenvectors of a physical Hessian.

## Modules / 模块及科学边界

1. **QuantumLogParser** extracts SCF energies, ZPVE, total H/G and their corrections, temperature, charge/multiplicity, angstrom coordinates, frequencies, displacement modes, Mulliken/Hirshfeld populations, and spin-resolved frontier levels. Gaussian `D` exponents are supported. H/G already include ZPVE; it is never added twice. Missing optional quantities are `null` or empty collections with warnings.
2. **ReactionThermochemistryEngine** builds R–PC–TS–INT–P relative Gibbs profiles in kcal/mol, reports barriers from both R and PC, and fits weighted BEP slope/intercept with covariance. Explicit caller confirmation of method, solvent, standard-state and reference compatibility is required. Composition, charge, multiplicity and temperature must match. Open electrochemical reactions need separately specified proton/electron reservoirs; the module does not silently subtract neutral and cation energies.
3. **SchNetE3InvariantEngine** uses 32 Gaussian radial bins, a 5 Å cosine envelope, distance-dependent filter MLPs and residual interactions. Cutoff is applied after the filter network to eliminate bias leakage beyond the cutoff. Atomic embeddings are conditioned on charge and multiplicity. Mean-pooled scalar heads can be trained on AIP and Eox targets; atomwise barrier contributions have an additive sum. These atomwise values are model decompositions, not measurable atomic barriers. All delivered predictions are untrained arbitrary-unit outputs. AIP and electrode-referenced Eox require different, explicitly defined target data.
4. **SyntheticDFTSandbox** writes explicit synthetic fixtures, verifies the intended -342.5 cm⁻¹ TS mode and thermochemical values, and keeps neutral substrate outside the illustrative fixed-charge pathway. PC/TS/INT/P labels and BEP coordination-center labels are software examples, not calculated mechanisms or catalyst screening results.
5. **PublicationArtifactGenerator** writes an ACS-style `booktabs` LaTeX table, a standalone table document, a 300-DPI composite PNG, a vector PDF, and a structured JSON database including plotted values, input hashes, source hash and software versions. Formatting is suitable for further manuscript preparation; journal acceptance or submission compliance has not been established.

The network processes one molecule per call, with dense O(N²) pair construction and a default cap of 512 atoms. It provides differentiable tensors and works with normal PyTorch optimizers. Symmetry/gradient checks and one training-step test establish software behavior, not predictive accuracy. Charge/spin conditioning is not a substitute for electronic-structure training data.

## Supported parser scope / 日志范围

Supported patterns are Gaussian 16 standard/input orientation and low-precision `Atom AN X Y Z` mode tables, plus ORCA 5/6 angstrom coordinates, indexed frequency tables and normal-mode row matrices. Both restricted and spin-separated frontier tables are handled. The reported gap is a molecular orbital eigenvalue gap, **not** a solid-state band gap, optical excitation, ΔSCF ionization energy or oxidation potential.

Exactly one normal termination is required. Concatenated jobs, repeated thermochemistry/frequency analyses, stale thermal corrections followed by another SCF/geometry, malformed population/mode tables, fractional orbital occupations and electron-count/multiplicity inconsistencies are rejected. Split multi-job outputs into individual jobs before parsing. Unsupported formats such as Gaussian high-precision displacement tables, periodic outputs, ghost atoms and alternate population analyses require additional format-specific coverage. The parser reads one text file into memory; it is not a streaming multi-gigabyte ingestion service.

Stationary-point validation requires an explicit optimization convergence marker and a complete spectrum (Gaussian 3N−6 or linear 3N−5; ORCA 3N including rigid modes). TS classification requires exactly one strictly negative frequency and a complete nonzero displacement vector. Tiny negative modes are not silently discarded. This is a necessary diagnostic: mode inspection and IRC are still needed to establish reaction connectivity. A frequency-only log can be parsed but does not certify a separately optimized structure.

## BEP statistics / 拟合解释

Weights are inverse variances, `w_i = 1 / sigma_i²`. With supplied independent absolute standard deviations, covariance is `(XᵀWX)⁻¹` and pointwise mean confidence bands use the normal quantile. Without sigma, residual variance is estimated with n−2 degrees of freedom and Student-t bands are used. SVD checks reject singular/ill-conditioned designs. Bands describe the fitted mean, not prediction intervals. X is treated as exact; correlated errors, uncertainty in reaction energies and mechanistically distinct pathways need a richer statistical model.

## Artifacts / 文件

- [Composite figure](results/qm_autobridge_composite.png) / [vector PDF](results/qm_autobridge_composite.pdf)
- [LaTeX table](results/qm_autobridge_table.tex) / [standalone document](results/qm_autobridge_report.tex)
- [Structured database](results/qm_autobridge_database.json) / [synthetic input logs](results/synthetic_logs)
- [Tests](../tests/test_qm_autobridge.py): independent hand-written fixtures, malformed input, thermal bookkeeping, weighted-regression covariance, symmetry, cutoff behavior, finite-difference gradients and a training step.

The stored demo has a designed barrier of 18 kcal/mol from R, 21 kcal/mol from PC, and reaction energy −9 kcal/mol. These are test targets. They must not be cited as computed chemistry in the existing paper. The prior manuscript PDFs are unchanged.

## Primary references / 已核实的原始来源

- [Gaussian, Thermochemistry in Gaussian](https://gaussian.com/wp-content/uploads/dl/thermo.pdf): electronic energies and thermal corrections.
- [ORCA 6.1 manual, Thermochemistry](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/thermochemistry.html) and [Vibrational Frequencies](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/frequencies.html): output conventions. Parsing printed totals preserves the calculation's original thermochemical treatment; it does not recompute RRHO/quasi-RRHO corrections.
- Schütt, K. T.; Kindermans, P.-J.; Sauceda, H. E.; Chmiela, S.; Tkatchenko, A.; Müller, K.-R. *SchNet: A continuous-filter convolutional neural network for modeling quantum interactions*. Advances in Neural Information Processing Systems 30 (2017). [Original proceedings paper](https://papers.neurips.cc/paper/6700-schnet-a-continuous-filter-convolutional-neural-network-for-modeling-quantum-interactions.pdf). This implementation is SchNet-inspired, not a reproduction of its trained benchmark models.
