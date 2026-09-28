# Four-species wall reaction network / 四物种壁面反应网络

This extension introduces explicit desired, parallel and consecutive pathways: **A → P, A → B, P → D**. It resolves all four concentration fields and obtains desired-product Faradaic efficiency from their integrated reaction and outlet fluxes. The rates, electron numbers and species names are transparent assumptions. The model is an executed numerical scenario, not an identified molecular mechanism, measured selectivity or validated process.

本扩展显式求解 A、P、B、D 四个物种浓度场，引入目标反应 A → P、并联副反应 A → B 和连续过氧化 P → D。目标产物 FE 从积分反应通量、出口净产物流和三条支路电流计算，因此不再结构性恒为 100%。全部物种名称、速率参数及每步二电子计量均为假设；结果不是实验选择性，也没有建立真实分子的元素配平或机理。

## 1. Equations and declared assumptions / 方程与假设

Each species obeys steady transport with fixed Poiseuille flow and equal, constant diffusivity:

$$\nabla\cdot(\mathbf u C_s-D\nabla C_s)=0,\quad s=A,P,B,D.$$

The inlet specifies total advective-diffusive flux equal to the imposed feed flux (Danckwerts); the outlet has zero axial diffusion and outward convection; the upper wall is insulating. The lower-wall outward flux vector is J = L C_w, with species order A, P, B, D and

$$L=\begin{pmatrix}k_1+k_2&0&0&0\\-k_1&k_3&0&0\\-k_2&0&0&0\\0&-k_3&0&0\end{pmatrix}.$$

Positive flux denotes removal from the liquid, while negative flux denotes generation into the liquid. If h = Δy/2, the half-cell diffusive relation C_c − C_w = h J/D gives the exact linear Robin elimination

$$C_w=(I+hL/D)^{-1}C_c,\qquad J=L(I+hL/D)^{-1}C_c.$$

The sparse system contains **four unknowns per cell**; P, B and D are not recovered by subtracting A from a prescribed total. Matrix column sums are zero because every abstract reaction converts one molecule into another. Thus total molecular amount should be conserved even though each species has its own reaction source or sink. This **does not establish elemental or mass conservation for a real organic reaction**, because no formula, molecular mass or atom mapping was supplied for A, P, B or D.

Shared finite-volume face fluxes, central diffusion and first-order upwind convection follow the standard conservative framework documented by [NIST FiPy](https://pages.nist.gov/fipy/en/latest/numerical/discret.html), checked on 2026-09-28. The implementation uses SciPy sparse matrices rather than FiPy. Both axial and transverse diffusion are present by default. Numerical upwind diffusion and unresolved local-gradient error remain; global balance closure is not a proof of mesh independence.

At D = 0 the diagnostic wall-access limit has zero wall flux and no reaction current. Wall concentrations are then undefined and exported as `null`, rather than assigning a fictitious interfacial state. Reverse flow, negative rates, negative concentrations and unsupported reverse pathways are rejected before solving. This irreversible scenario does not model cathodic reversal.

默认参数为 L = 60 mm、H = 300 µm、W = 12 mm、D = 1.1 × 10^-9 m²/s、T = 298.15 K、Q = 450 µL/min、入口 A = 50 mM，其余物种入口浓度为零。默认网格为 80 × 32，即 10,240 个标量未知数。流速在横向条带内精确积分，以保持指定体积流量。入口为总通量边界，四个场全部显式求解，壁面条件通过矩阵消元接入物料方程。

## 2. Kinetics and Faradaic accounting / 动力学与电流分账

For pathway i, an assumed potential-sensitive rate is

$$k_i(\eta)=k_{i,ref}\exp\left[\beta_i F(\eta-0.45)/(RT)\right].$$

| Pathway | Reference rate / m s^-1 | Effective exponent β | Assumed electrons |
|---|---:|---:|---:|
| A → P | 5.0 × 10^-5 | 0.50 | 2 |
| A → B | 3.0 × 10^-6 | 0.65 | 2 |
| P → D | 8.0 × 10^-6 | 0.85 | 2 |

η is an external model potential index, not a solved electrode overpotential. Distinct β values deliberately create a scenario in which the overoxidation branch grows more rapidly. They were not fitted or selected from chemical data. Each pathway remains an irreversible first-order wall reaction; no two-electron Butler–Volmer/Nernst thermodynamics is implied. Explicit rate overrides are supported and recorded for limiting tests.

At η = 0.48, the rates are k_AP = 8.9644260 × 10^-5, k_AB = 6.4082505 × 10^-6 and k_PD = 2.1583865 × 10^-5 m/s. These positive constants, surface A/P concentrations and integrated wall areas determine the extents ξ_AP, ξ_AB and ξ_PD in mol/s. Each channel carries I_i = n_i F ξ_i. With the default two electrons per step,

$$I=2F(\xi_{AP}+\xi_{AB}+\xi_{PD}),$$
$$\Delta\dot n_P=\xi_{AP}-\xi_{PD},\quad\Delta\dot n_B=\xi_{AB},\quad\Delta\dot n_D=\xi_{PD},$$
$$\mathrm{FE}_{P,net}=100\frac{2F\Delta\dot n_P}{I}.$$

An independent outlet-based electrical balance is I = 2F(Δṅ_P + Δṅ_B + 2Δṅ_D). The factor of two on D accounts for both formation of its precursor P and subsequent oxidation. Gross AP charge fraction is a different quantity: some initially formed P is later destroyed.

净 P 收率等于净 P 出口增量除以 A 进料；净 P 选择性等于净 P 增量除以 A 消耗；净 P FE 等于生成净 P 对应的二电子电荷当量除以三条支路的总电流。三者分母不同，不能混用。若入口含 P，净出口增量可以为负，此时报告的是净差量 FE，并非“绝对生成产物的 FE”；不应为了好看把它截到零。

## 3. Executed baseline and parameter pool / 基准结果与参数池

The default model gives A conversion **58.4087489913%**, net P yield **11.5438481742%**, net P selectivity **19.7639024522%**, net P FE **11.3870659398%**, and total current **73.3603391064 mA**. Outlet concentrations are A = 20.7956255044, P = 5.7719240871, B = 1.9484024642 and D = 21.4840479443 mM. Their sum remains 50 mM to floating-point precision.

| Channel | Extent / mol s^-1 | Current / mA |
|---|---:|---:|
| A → P | 2.0441979024 × 10^-7 | 39.4470227056 |
| A → B | 1.4613018482 × 10^-8 | 2.8198838829 |
| P → D | 1.6113035958 × 10^-7 | 31.0934325179 |

The gross AP charge fraction is 53.7716%, whereas net desired FE is only 11.3871%. This deliberately retained low-performance result shows the distinction between passing current through the desired first step and collecting intact product at the outlet. A conversion alone is insufficient to rank this network.

The 49-point grid spans Q = 100–1,500 µL/min and η = 0.20–0.75, with seven uniformly spaced values along each coordinate. Model net FE ranges from **0.0001740358%** at Q = 100, η = 0.75 to **98.5898713586%** at Q = 1,500, η = 0.20. The high-FE condition has only 1.0883060% A conversion and 1.0731501% net P yield; the nearly complete-conversion low-FE condition has only 0.0002970% net P yield. These outcomes are consequences of assumed rates, not experimentally discovered optima or validated operating limits.

中文结果：默认模型总电流 73.3603 mA、A 转化 58.4087%，但净 P 收率仅 11.5438%、净 P FE 仅 11.3871%，原因是连续过氧化支路消耗大量 P。49 个模型工况的 FE 非恒定，但高 FE 不等于高转化或高产率；这些负情景均保留，没有截断或人为设置选择性。

## 4. Numerical verification / 数值验证

The saved study contains **62 four-species PDE solves plus one frozen single-species comparator, 63 PDE solves in total**. The four-species count is one baseline + 49 pool evaluations + three grid cases + five limiting cases + four analytic-limit refinements. Eight invalid requests are rejected before the PDE solve. The independent unit-test solves are explicitly excluded from this study count. The solve registry records every successful study evaluation; repeated baseline/grid settings are counted because they were actually re-evaluated.

All 62 network cases passed the declared convergence and balance checks. The maximum relative species-balance error is **3.3611 × 10^-10**, total-molar-balance error **3.7649 × 10^-10**, and charge error **4.5102 × 10^-11**, normalized to feed scales. The largest errors occur in the artificial high-diffusion analytic limit, where the sparse matrix is more ill-conditioned. At the ordinary baseline the corresponding errors are around 10^-15. No clipping is applied to concentrations, yield, FE or solver status.

| Grid | Total unknowns | Net P yield / % | Net P FE / % |
|---|---:|---:|---:|
| 40 × 16 | 2,560 | 11.4970663282 | 11.3690377267 |
| 80 × 32 | 10,240 | 11.5438481742 | 11.3870659398 |
| 160 × 64 | 40,960 | 11.5672500428 | 11.3956066765 |

The last refinement changes P yield by **0.0234018686 percentage points** and FE by **0.0085407367 percentage points**. The largest solved system remains below the limit of 50,000 total species unknowns. This is a three-grid sensitivity study for the baseline; the finest tested grid is not an exact solution or an accuracy certificate for every pool point.

For no side reactions, k_AP = 5 × 10^-5 and k_AB = k_PD = 0, the network agrees with the frozen transport solver (its reverse rate also set to zero) to **1.42 × 10^-13 percentage points** in conversion and -6.94 × 10^-17 A in current. For parallel-only k_AP:k_AB = 2:1, the outlet P:B ratio is 2:1 and net FE is 66.6667%, independently of the common transport limitation. Zero kinetics preserves an arbitrary inlet mixture; zero diffusivity blocks electrode access. Product-only feed with P → D yields net P FE ≈ **-100%**, retained and labeled according to its net definition.

The independent high-mixing check uses plug velocity, D_y = 10^-4 m²/s, D_x = 0, k_AP = 2 × 10^-6, k_AB = 0 and k_PD = 4 × 10^-6 m/s. With exposure θ = LW/Q and a = k_AP + k_AB, b = k_PD, direct solution of dA/dθ = -aA and dP/dθ = k_AP A - bP gives

$$A=A_0e^{-a\theta},\qquad P=P_0e^{-b\theta}+\frac{k_{AP}A_0}{b-a}(e^{-a\theta}-e^{-b\theta}).$$

For a = b the second term is k_AP A_0 θ e^-aθ. The code evaluates the near-equal-rate expression with `expm1` to avoid cancellation. The formulas are obtained by direct integration of these stated equations; they are not literature-derived parameter data. For the tested case the exact mixed-cross-section P concentration is **7.2087720656 mM**. Axial meshes 40, 80, 160 and 320 reduce the largest outlet species discrepancy from approximately 0.043509 → 0.021844 → 0.010958 → **0.005502 mM**, consistent with first-order convergence.

Seventeen tests independently cover the sequential and parallel analytic formulas, equal/near-equal rates, all species balances, wall kinetic versus diffusion flux, charge accounting, the no-side-reaction comparator, invalid-input rejection, zero-reaction/zero-diffusion limits, and deliberate bad-solver injection. All 17 passed. Passing these tests verifies code behavior for these cases, not chemical validity.

## 5. Files and interface / 文件与接口

```text
python electratwin/scripts/reaction_network.py
python -m unittest discover -s tests -p test_electratwin_network.py -v
```

Use the existing chemistry environment with its numerical DLL path configured and OMP/BLAS/MKL thread counts set to one. No package installation, hardware connection or modification of the frozen transport/control code is required.

The callable interface is `solve_network(flow_rate_uL_min=450, potential_index_V=0.48, nx=80, ny=32, ...)`. It returns `schema_version`, `evidence_class`, `inputs`, `summary` and optional `field`. `summary` provides `conversion_A_pct`, `net_P_yield_pct`, `net_P_selectivity_pct`, `net_P_faradaic_efficiency_pct`, `current_A`, `channel_current_A`, `extent_mol_s`, `outlet_mM`, species residuals and total/charge residuals. `field.concentration_mM` contains all four [nx, ny] arrays at cell centers; `include_field=False` suppresses export of arrays without changing the computation. Optional `rate_constants_m_s` must contain exactly AP, AB and PD; `inlet_mM` is ordered A, P, B, D.

Artifacts: [source](scripts/reaction_network.py), [study summary and source hashes](results/reaction_network/study_summary.json), [baseline fields](results/reaction_network/baseline_solution.json), [49-point pool](results/reaction_network/parameter_pool.csv), [grid study](results/reaction_network/grid_convergence.csv), [limiting cases](results/reaction_network/limiting_cases.csv), [analytic comparison](results/reaction_network/analytic_sequential_limit.csv), [frozen-solver comparison](results/reaction_network/frozen_solver_comparison.json), [input rejections](results/reaction_network/input_rejections.json), [solve registry](results/reaction_network/solve_registry.csv), [tests](../tests/test_electratwin_network.py), and [test transcript](results/reaction_network/test_results.txt).

适用边界：该网络未包含真实分子结构、元素配平、离子迁移、电位分布、热场、气泡、催化剂失活或分离损失；单独的电流分账与摩尔守恒不能补足这些证据。所有数值必须保留“假设网络计算”标签，不应改写为真实实验、工业认证或自动优化出的真实配方。
