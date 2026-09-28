# Geometric pores and synthetic BET inversion / 几何孔道与合成 BET 反演

## Source audit / 原始程序审计

The original `ReticularPOPEngine` is executed separately through its extracted class, without running the neural/MD/NEB workflow. Its mask is `X²+Y² < 9.8²`, a circle on a square grid including both endpoints. No framework atom coordinates, bonds, graph topology, atomic masses or pore-size histogram are constructed. The `pore_probe_radius` argument is never read. Eight calls spanning two cell lengths and four probe radii confirm the unchanged answer at each cell length. The labels “hcb”, “Cu-N4”, “Connolly volume” and “PSD” are not established by that calculation.

单独提取并执行原始 `ReticularPOPEngine` 类，没有再次运行神经网络、MD 或 NEB。原掩膜是平方和小于固定半径平方，实际为方形网格中的圆；网格含两个端点。代码不构建框架原子、键、拓扑图、原子质量或孔径直方图，探针半径参数未被使用。两个胞长、四个探针半径的八次调用证实：同一胞长下输出完全不变。因此不能据此确认 hcb、Cu-N4、Connolly 体积或 PSD。

At the source's 26 Å cell length and 40×40 grid, the reported 42.5% is an endpoint-grid area estimate. The exact zero-probe circular fraction is 44.6329228%; the midpoint-grid estimate at the same resolution is 44.75%. The source's 631.2 m²/g is the prescribed monolayer capacity 145 cm³(STP)/g multiplied by 4.353 and rounded, independent of its geometry. A synthetic conversion is dimensionally possible, but the underlying capacity is not derived from material mass or adsorption evidence.

原调用胞长 26 Å、40×40 网格给出 42.5%，其实是端点网格面积估计。零探针圆形精确面积分数为 44.6329228%，相同分辨率的中点积分为 44.75%。631.2 m²/g 则是直接指定的单层容量 145 cm³(STP)/g 乘 4.353 后舍入，与几何无关；量纲换算本身可以成立，但输入容量没有原子质量或吸附证据支持。

## Reviewed geometry / 复核几何

Three analytic 2D benchmarks use a 26 Å square cell: a circular channel of radius 9.8 Å; a regular hexagonal channel of **apothem** 9.8 Å; and a periodic solid disk of radius 3 Å centered at (12.4, 2.5) Å. They are artificial geometries, not chemical frameworks. A probe center is allowed when its nearest-wall clearance exceeds the probe radius. For a circular channel, `A=π max(R−r,0)²`; for the hexagon, `A=2√3 max(a−r,0)²`. The hexagonal circumradius is `2a/√3`, not a. Periodic disk distances use `(Δ+L/2) mod L−L/2` in each coordinate. Its exact area also handles overlapping periodic exclusion disks through disk–square intersection geometry.

三个解析二维基准都使用 26 Å 方胞：半径 9.8 Å 圆孔、**内切半径** 9.8 Å 正六边形孔，以及中心位于 (12.4, 2.5) Å、半径 3 Å 的周期性实心圆盘。它们不是化学框架。探针中心到最近壁面的距离大于探针半径才可进入。圆形和六边形面积分别按上述侵蚀公式计算；六边形内切半径不等于外接圆半径。圆盘采用最小镜像距离，并以圆与最小镜像方胞交叠面积处理排斥区重叠。

Five midpoint grids (40, 80, 160, 320, 640 per axis) and explicit probe radii produce 15 clearance fields, **1,636,800 point evaluations and 90 masks**. The finest grid's maximum absolute area-fraction error is 0.000474912, or 0.0474912 percentage points, across all tested geometries/probes. A 64-point independent comparison with nine translated disk images differs by at most 3.55×10⁻¹⁵ Å. Grid errors need not decrease monotonically because of boundary alignment.

五级中点网格产生 15 个净空距离场、**1,636,800 次点评估和 90 个掩膜**。640×640 网格在全部基准中的最大面积分数误差为 0.000474912，即 0.0474912 个百分点。64 个点的最小镜像距离与显式九胞比较，最大差为 3.55×10⁻¹⁵ Å。边界与网格对齐会使误差出现非单调变化。

| Geometry / 几何 | Probe / 探针 Å | Exact center-accessible fraction / 精确中心可达分数 | 640² estimate / 网格估计 |
|---|---:|---:|---:|
| Circular channel / 圆孔 | 0 | 0.4463292285 | 0.4463378906 |
| Circular channel / 圆孔 | 1.82 | 0.2959436048 | 0.2958691406 |
| Hexagonal channel / 六边形孔 | 1.82 | 0.3263245214 | 0.3258496094 |
| Periodic disk exterior / 周期圆盘外部 | 1.82 | 0.8920314539 | 0.8920336914 |

The saved 120-bin histogram rows and 123 CDF rows describe **area-weighted point-to-wall clearance**. For circle and hexagon, the normalized analytic CDF is `1−(1−t/a)²` on `0≤t≤a`. The maximum CDF error at the 41 checked bin edges per geometry is 0.000687445. This quantity is neither a distribution of pore diameters nor an experimental PSD. Probe-center area also differs from probe-occupiable area, accessible surface area and 3D volume; terminology follows the distinctions in the [Zeo++ documentation](https://www.zeoplusplus.org/examples.html), but Zeo++ was not executed.

保存的 120 行直方图和 123 行 CDF 表示**按二维面积加权的点到壁面净空距离**。三个几何各有 41 个检验边界，最大 CDF 差为 0.000687445。这不是孔直径分布，更不是实验 PSD。探针中心可达面积也不同于探针占据面积、可达表面积和三维体积；本次没有运行 Zeo++。

## Synthetic inverse benchmark / 合成反演基准

The source's 20 prescribed points switch from BET to an empirical exponential at relative pressure 0.35. The left limit is 219.531691 cm³/g, while the right value is 424.311058 cm³/g: a **204.779367 cm³/g jump (+93.2801%)**. This is imposed by the branch definition, not evidence of capillary condensation, Type IV classification or adsorption/desorption hysteresis.

原始 20 个指定点在相对压力 0.35 处由 BET 切换至经验指数式；左右值分别为 219.531691 和 424.311058 cm³/g，**跳增 204.779367 cm³/g（93.2801%）**。该跳变来自分支定义，不能据此证明毛细凝聚、IV 型等温线或吸脱附滞后。

BET linearization fits `p/[V(1−p)] = i+s p`, with `Vm=1/(s+i)` and `C=1+s/i`. All intercepts, slopes, failed positivity checks and residuals are retained. Eight windows are fitted to each of source-rounded points, unrounded source formulas and exact BET controls: 24 deterministic fits. Selected diagnostics test positive parameters, increasing `V(1−p)`, and placement of `p_m=1/(1+√C)` inside the fitted range. These checks are not a complete standards-compliance or physical-validity assessment. The conventional nitrogen area calculation uses cross-section 0.162 nm², Avogadro constant 6.02214076×10²³ mol⁻¹, and stated STP molar volume 22414 cm³/mol (273.15 K, 1 atm), yielding factor 4.352577867 m² per cm³(STP). See [NIST SP 960-17](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication960-17.pdf).

BET 线性化采用上述公式，保留截距、斜率、非正参数及残差。舍入原始点、未舍入分段公式和精确 BET 对照分别在八个窗口拟合，共 24 次。另检查 `V(1−p)` 单调性及单层压力是否落入窗口，但不声称完全符合标准或物理适用性。氮分子截面积、阿伏伽德罗常数和 STP 体积约定均显式记录，换算系数为 4.352577867。

| Source-rounded window / 舍入点窗口 | Points / 点数 | Vm / cm³(STP) g⁻¹ | C | Diagnostic / 诊断 |
|---|---:|---:|---:|---|
| 0.01–0.30 | 6 | 145.000199 | 114.998314 | Positive; low-branch recovery / 正参数，恢复低压分支 |
| 0.10–0.35 | 5 | 145.002011 | 114.965103 | Monolayer pressure outside window / 单层压力不在窗口内 |
| 0.01–0.45 | 9 | 235.485444 | 10.566237 | `V(1−p)` fails monotonicity / 单调性失败 |
| 0.35–0.95 | 13 | 48.803456 | −1.189751 | Negative C, area withheld / C 为负，不报告面积 |
| 0.01–0.95 | 20 | 65.615805 | −3.775999 | Negative C, area withheld / C 为负，不报告面积 |

The exact BET controls recover prescribed Vm=145 and C=115 to maximum absolute errors 8.53×10⁻¹⁴ and 3.44×10⁻¹². An additional **64-seed** inverse benchmark perturbs six low-pressure points with independent mean-one lognormal noise of log standard deviation 0.01. Its 384 observations support 128 fits over two windows, for **152 total fits**. For window 0.01–0.30, fitted C has 2.5–97.5% repetition quantiles 104.439–136.165; excluding the lowest-pressure point expands this to 93.963–166.296. These are conditional synthetic-noise distributions, not confidence intervals from real adsorption data. OLS after linearization also changes the noise weighting.

精确对照恢复指定参数，最大绝对误差分别为 8.53×10⁻¹⁴ 与 3.44×10⁻¹²。额外 **64 个种子**对六个低压点加入均值乘子为 1、对数标准差为 0.01 的独立对数正态噪声，产生 384 个观测、128 次两窗口拟合，总计 **152 次拟合**。保留最低压力点时 C 的重复分位范围为 104.439–136.165，去除该点则扩为 93.963–166.296。这是条件合成噪声结果，不是实际吸附实验置信区间；线性化 OLS 也改变了噪声权重。

## Files and reproduction / 文件与复现

```bash
python synthapore/scripts/pore_reviewed.py --pilot
python synthapore/scripts/pore_reviewed.py
python -m unittest discover -s tests -p test_synthapore_pore.py -v
```

[Code / 代码](scripts/pore_reviewed.py), [summary and hashes / 汇总与哈希](results/pore/summary.json), [geometry convergence / 几何收敛](results/pore/geometry_grid_convergence.csv), [endpoint bias / 端点偏差](results/pore/endpoint_bias.csv), [clearance CDF / 净空 CDF](results/pore/clearance_CDF.csv), [source audit / 原始审计](results/pore/source_pore_audit.json), [window fits / 窗口拟合](results/pore/BET_window_fits.csv), [noise summary / 噪声汇总](results/pore/noise_summary.csv), [16 passing tests / 16 项测试通过](results/pore/test_results.txt).

The separate pilot has six fields, 24,000 point evaluations and 40 fits; these are excluded from main totals. Tests cover geometric limits, scale covariance, minimum images, exact CDF normalization, input rejection, parameter recovery and negative-fit retention. No external structures or datasets were downloaded, and no atomistic material, GCMC calculation or experiment was produced.

独立先导计算的六个距离场、24,000 次点评估和 40 次拟合不计入正式总数。测试覆盖几何极限、尺度关系、最小镜像、CDF、非法输入、参数恢复及负结果保留。本次没有下载外部结构或数据，也没有生成原子级材料、GCMC 或实验结果。
