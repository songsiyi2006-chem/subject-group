# Reviewed conformer and reaction-graph calculations / 构象与反应图复核

## Methods and evidence / 方法与证据

The executed study uses the installed RDKit 2026.03.5 and a single CPU thread. The seven original structures are retained. Each receives a budget of 24 ETKDGv3 conformers for each of seeds 20260928, 20260929 and 20260930. `params.pruneRmsThresh=0.15` Å sets embedding pruning through the parameter object; `MMFFOptimizeMoleculeConfs(maxIters=1000, mmffVariant="MMFF94")` returns the optimization status and energy in kcal/mol. Status 0 is required for descriptor weighting, while every returned geometry and status remains in the files. Requested conformers, returned geometries and optimization failures are separate quantities. Embedding failure counters count internal failure events, not necessarily failed requested conformers. [RDKit conformer documentation](https://www.rdkit.org/docs/RDKit_Book.html); [MMFF API](https://www.rdkit.org/docs/source/rdkit.Chem.rdForceFieldHelpers.html).

本次使用已安装的 RDKit 2026.03.5，单 CPU 线程计算。保留原始七个结构；每个结构以三个种子分别请求 24 个 ETKDGv3 构象。嵌入阶段 RMSD 剪枝阈值为 0.15 Å；随后进行最多 1000 次迭代的 MMFF94 优化。收敛状态、能量及所有返回的坐标均保留，仅状态 0 的结果参加加权。请求数、剪枝后返回数和优化失败数分开统计；嵌入内部失败事件不等于缺失构象数。

Within each seed and then across the pooled seeds, converged geometries are sorted by energy and greedily deduplicated at a symmetry-aware aligned heavy-atom RMSD of 0.35 Å. Alignment operates on copies. All hydrogens are explicit for MMFF, SASA and radius of gyration. This chosen threshold can merge hydrogen rotamers or nearby geometries; the selected representative is the lowest-energy member encountered. The implementation does not infer conformer degeneracies or basin volumes.

每个种子及合并种子的收敛结构均按能量排序，以对称性校正、对齐后的重原子 RMSD 0.35 Å 贪心去重；对齐使用副本，保存的坐标不改变。该阈值可能合并氢转子或相近构象，保留先遇到的最低能量代表。算法不推断构象简并度或势阱体积。

For each saved set of minima, weights are recomputed at 250, 298.15 and 350 K as `exp(-(E-Emin)/(R*T))`, with `R=8.31446261815324/4184` kcal mol⁻¹ K⁻¹. Effective minima count is `1/sum(w²)`. These are unit-degeneracy **force-field minima weights**, not Gibbs free-energy populations in solution. Entropy, solvent, electrode, protonation equilibria and Hessian verification are absent. A converged optimizer does not prove a true local minimum. Rg uses atomic masses and the center of mass; a geometric-center, unweighted comparison reproduces the source's mathematical definition.

同一组保存的极小值在 250、298.15 和 350 K 重加权，采用相对 MMFF 能量和上述气体常数。有效极小值数为 `1/sum(w²)`。这些是简并度取 1 的**力场极小值权重**，不是溶液 Gibbs 自由能布居；未计入熵、溶剂、电极、质子化平衡，也未进行 Hessian 检验。质量加权回转半径采用质心，同时输出原代码对应的几何中心、无质量权重值。

`classifyAtoms` returned zero radii for these non-protein molecules. SASA therefore uses positive RDKit periodic-table vdW radii and a 1.4 Å probe. Primary SASA is the mean of Lee–Richards calculations over 24 fixed pseudorandom rotations (seed 17381, QR-derived proper rotation matrices). Raw orientation values, their spread, and the difference between the first 12 and all 24 means are saved. This reduces finite surface-integration orientation artifacts; it does not establish numerical convergence. [FreeSASA API](https://www.rdkit.org/docs/source/rdkit.Chem.rdFreeSASA.html); [atomic radii API](https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html).

这批非蛋白分子由 `classifyAtoms` 得到的半径为零，因此 SASA 改用 RDKit 周期表的正值范德华半径及 1.4 Å 探针。主结果是固定种子产生的 24 个旋转下 Lee–Richards 面积平均值。保存全部方向面积、离散程度及前 12 个与全部 24 个平均值之差。这降低了离散积分的方向误差，但不能宣称积分已收敛。

## Executed results / 已执行结果

The final study requested 504 conformers in 21 replica runs, received 120 after embedding pruning, optimized all 120 to status 0, and retained 36 pooled distinct minima. It produced 84 ensemble rows, 357 weight rows, 2880 surface evaluations and 14 SDF coordinate files. The separate pilot requested 6 conformers and optimized 2; it is excluded from the main counts. No QM or wet-lab measurements were performed.

正式计算包括 21 个种子重复、504 个请求构象；剪枝后返回并优化 120 个，均收敛，合并去重后为 36 个极小值。输出 84 行系综统计、357 行权重、2880 次表面积方向积分及 14 个 SDF 坐标文件。独立先导计算请求 6 个、优化 2 个，不计入正式数量。没有量子化学或湿实验测量。

| Molecule / 分子 | Returned / 返回 | Pooled minima / 合并极小值 | Effective count at 298.15 K / 有效数 | SASA / Å² | Mass Rg / Å |
|---|---:|---:|---:|---:|---:|
| Target / 目标 | 6 | 3 | 2.751691 | 455.326579 | 3.608485 |
| Melatonin / 褪黑素 | 49 | 18 | 7.134821 | 470.631031 | 3.323095 |
| Caffeine / 咖啡因 | 3 | 1 | 1.000000 | 371.059477 | 2.464656 |
| 2-Phenylquinoline / 2-苯基喹啉 | 3 | 1 | 1.000000 | 420.532069 | 3.207699 |
| Tryptophol / 色醇 | 33 | 6 | 3.793324 | 356.289274 | 2.458393 |
| Indoline / 吲哚啉 | 8 | 1 | 1.000000 | 300.839276 | 1.961946 |
| Benzofuran ester / 苯并呋喃酯 | 18 | 6 | 3.854778 | 404.491154 | 3.005735 |

The target has formula C15H13NO. Explicit ring numbering confirms methoxy at indole C5 and phenyl at C2: source zero-based indices are N1=6, C2=7, C3=8, C3a=9, C4=10, C5=2, C6=3, C7=4 and C7a=5. Its weighted unweighted-definition Rg is 4.094013 Å, compared with mass-weighted 3.608485 Å. Melatonin's minimum energy differs by 0.574056 kcal/mol across seeds and its weighted mass Rg spans 0.069286 Å; these are incomplete-search diagnostics, not calibrated physical uncertainty. The largest first-12 versus all-24 SASA mean difference is 1.167181 Å² over the 120 optimized structures. The largest individual orientation range is 16.577209 Å², so small SASA changes must not be overinterpreted.

目标分子式为 C15H13NO。显式吲哚编号确认甲氧基在 C5、苯基在 C2，具体原子索引见上文及 JSON。无质量权重定义给出 4.094013 Å，而质量加权结果为 3.608485 Å。褪黑素三个种子的最低能量跨度为 0.574056 kcal/mol，质量加权 Rg 跨度为 0.069286 Å，说明搜索仍不完全；这些不是经过校准的物理不确定度。在全部 120 个结构中，12 与 24 个方向平均 SASA 的最大差为 1.167181 Å²，单方向范围最大为 16.577209 Å²，不能过度解释很小的 SASA 差异。

## Reaction topology / 反应拓扑

The original unmapped reaction has reactant formula C22H21NOS and product formula C21H17NOS, with products-minus-reactants differences C=−1 and H=−4. Heavy-atom counts are 25 and 24. `CSc1ccccc1` is thioanisole, while thiophenol is `Sc1ccccc1`. The source is retained as an insufficient-mapping, unbalanced input; no reaction center or mechanism is inferred. Equal heavy-atom counts would not be enough either: `CC>>CO` has two heavy atoms on each side but changes element and hydrogen inventories.

原始未映射反应左右分子式为 C22H21NOS 与 C21H17NOS，产物减反应物少 1 个 C、4 个 H，重原子数由 25 变为 24。`CSc1ccccc1` 是苯甲硫醚，苯硫酚应为 `Sc1ccccc1`。保留原始缺陷，不据此推断反应中心或机理。重原子数相同也不代表元素守恒，反例 `CC>>CO` 已纳入检验。

`mapped_bond_changes` requires a positive, unique map for every explicit atom, identical map sets, and consistent element/isotope identities. It compares bonds keyed by sorted map pairs and records order, formation, cleavage and atom-property changes. The educational, balanced ethanol→acetaldehyde+H2 example has four edits: C2–O3 single→double; C2–H4 and O3–H5 broken; H4–H5 formed. This is bookkeeping and a topology test, not an electrode mechanism. Hydrogen/charge balance includes bracket and implicit hydrogens; stereochemical reaction validation and automated atom mapping are outside scope.

映射函数要求每个显式原子具有正的唯一映射号、左右映射集合一致，且元素及同位素身份不变。函数按映射号对比较键，输出成键、断键、键级和原子属性变化。配平教学例“乙醇→乙醛+H2”恢复上述四个变化；它只是原子记账及拓扑测试，不能作为电极机理。计数包括括号氢及隐式氢；不提供自动映射或反应立体化学验证。

## Reproduction and files / 复现与文件

```bash
python electrograph/scripts/structure_reviewed.py --pilot
python electrograph/scripts/structure_reviewed.py
python -m unittest discover -s tests -p test_electrograph_structure.py -v
```

[Code / 代码](scripts/structure_reviewed.py), [summary / 汇总](results/structure/summary.json), [conformers / 构象](results/structure/conformers.csv), [ensemble statistics / 系综统计](results/structure/ensemble_statistics.csv), [weights / 权重](results/structure/minima_weights.csv), [surface quadrature / 表面方向积分](results/structure/surface_quadrature.csv), [seed variability / 种子变异](results/structure/seed_variability.csv), [reaction audit / 反应审计](results/structure/reaction_audit.json), [18 passing tests / 18 项测试通过](results/structure/test_results.txt).

The numerical/representation tests cover normalized weights and energy units, mass-weighted Rg invariance, radii and SASA bounds, coordinate preservation, convergence filtering, deterministic coordinates, unsupported inputs, formula counterexamples and mapped bond edits. They do not validate chemical accuracy. SDF coordinates are written at the standard writer precision; JSON energies and descriptor values retain floating-point precision. RDKit software uses a BSD license; its online documentation uses CC BY-SA 4.0. No dataset or checkpoint was downloaded. [License and citation source](https://www.rdkit.org/docs/Overview.html).

数值和表示测试覆盖权重、单位、质量加权 Rg 的旋转平移不变性、SASA 界限、坐标保留、收敛筛选、确定性、非法输入及映射键变化；它们不是化学准确性验证。SDF 采用标准写入精度，JSON 保留浮点数值。沿用 RDKit 软件及文档许可，没有下载数据集或模型权重。
