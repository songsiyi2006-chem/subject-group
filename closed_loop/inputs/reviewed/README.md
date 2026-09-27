# Reviewed DFT input templates

Twelve input files cover two molecules (target 5-methoxy-2-phenyl-1H-indole and indoline control), each as neutral singlet and radical cation doublet. These are proposed calculations, **not executed or engine-validated DFT results**. Geometry provenance and SHA256 hashes are in [input_manifest.json](input_manifest.json).

For each state, choose one engine protocol:

- Gaussian 16: `.gjf` contains B3LYP-D3BJ/def2-SVP/SMD(MeCN) optimization and frequencies, followed by an M06-2X/def2-TZVP/SMD(MeCN) single point in Link1. Open-shell states use the unrestricted variants. Inspect both steps and the first-step frequencies; Link1 can run even when a stationary point is a saddle.
- ORCA 6.1: run `*_opt.inp`, inspect convergence, frequencies and spin, then run `*_sp.inp` from that same directory. The SP input requires the **new** `*_opt.xyz` produced by the ORCA optimization. It is intentionally absent until that job runs. B3LYP/G selects Gaussian-style VWN3 correlation. Both stages explicitly enable SMD; RIJCOSX and integration grids differ from Gaussian, so numerical equality is not assumed.

Starting coordinates reuse previously completed neutral GFN2-xTB/ALPB(MeCN) optimizations from this repository. Both charge states start from that same geometry but request independent DFT optimization. They are not optimized DFT radical geometries. The two `.xyz` files describe starting structures only.

Inputs request two CPU processes: Gaussian memory 3 GB total; ORCA maxcore 1000 MB per process (nominal 2000 MB, additional overhead possible). A proposed scheduler allocation is at least 4 GB for one job; assess the actual engine's memory estimate, free RAM and local policy before execution. This is not a measured performance or memory guarantee. The original 16-process/32-GB templates are archived separately.

Required scientific review of eventual outputs: SCF and geometry convergence; frequency analysis with low-frequency inspection; wavefunction stability and doublet spin expectation; method/basis/grid/solvation sensitivity; conformer sampling; consistent thermal and standard-state conventions. Input structure/parity checks cannot replace those results. Pure-MeCN SMD is an approximation to the proposed MeCN/HFIP mixture and does not model an electrode interface or counterion explicitly.

Official syntax references: [ORCA solvation](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html), [ORCA functional definitions](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/DensityFunctionalTheory.html), [ORCA memory](https://www.faccts.de/docs/orca/6.1/tutorials/first_steps/memory.html), [Gaussian SCRF](https://gaussian.com/scrf/), [Gaussian DFT](https://gaussian.com/dft/). Checked 2026-09-28; no quantum executable was found on the checked PATH, and no job was submitted.
