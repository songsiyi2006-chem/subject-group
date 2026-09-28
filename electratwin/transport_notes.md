# Transport model: definition, evidence and limits / 输运模型定义、证据与边界

The reviewed solver is an **executed, uncalibrated numerical model**. It is not a first-principles electrochemical calculation, an experimentally validated digital twin, or an industrial controller. The 39 transport solves and 15 independent tests below verify numerical properties under declared assumptions; they cannot establish the actual chemistry, kinetics, selectivity, equipment performance or process safety.

修订求解器属于**已执行、未标定的数值模型**。39 次输运求解及 15 项独立测试用于检验离散方程与守恒性质，不能据此证明真实反应、动力学、选择性、设备表现或工艺安全。A 与 P 是一比一转换的抽象分子物料记账标签，并非完整配平的有机反应式。

## 1. Source audit / 原程序审计

All line references below refer to the immutable [submitted source](source/electratwin_core.py).

| Source location | Observed problem | Reviewed treatment |
|---|---|---|
| Lines 204–207 | The implemented stencil has transverse diffusion and axial upwind convection, but no axial diffusion despite the displayed governing equation. | Both axial and transverse diffusive face fluxes are assembled by default. |
| Lines 213–219 versus green metrics electron count | Wall consumption uses current divided by F; downstream default metrics use two electrons. | The same specified electron count converts net molecular wall flux to current: I = nF times turnover. |
| Lines 223–224, 238 | Change between iterations is called a residual, while `converged` is always true, including iteration exhaustion. | Check the linear equation residual, global material/charge closure and concentration bounds. An injected invalid solver result produces `converged=false`. |
| Lines 228–229 | Current integration includes the inlet wall node where the inlet concentration was kept prescribed and the reactive wall condition was not updated. | Integrate exactly the same wall-face fluxes used as sinks in the material equations, over all finite-volume wall faces. |
| Lines 217–219 | Reverse rate has no explicit product concentration; it can act as a source even when no product is present. | Effective reverse rate is proportional to surface product concentration. |
| Lines 220, 235 and 102 | Concentration, conversion and FE are clipped; conversion cannot exceed 99.9% and inconsistencies can be hidden. | No scientific values are clipped; negative reverse current, trace roundoff and failed numerical checks remain visible. |
| Lines 229 and 246 | `np.trapz` is unavailable in the installed NumPy version; current is rounded before later calculations. | New finite-volume integration uses direct face summation and full-precision results. The separate source compatibility patch is not part of this new model. |

原文含轴向扩散的方程与实际离散代码不一致；电极通量使用一电子计量，绿色指标却默认二电子；入口壁节点被计入电流但未执行同一反应边界；收敛标识恒为真，截断还会掩盖不合理转化率或 FE。上述问题必须分别修正，单独让原脚本运行成功并不能解决它们。

## 2. Mathematical definition / 数学定义

For a steady, incompressible, isothermal channel with prescribed velocity u(y), each species obeys

$$\nabla\cdot(\mathbf u C_s-D\nabla C_s)=0,\qquad s\in\{A,P\}.$$

The model assumes identical diffusivities and one-to-one conversion, so the total concentration is constant: C_A + C_P = C_in,A + C_in,P = C_T. One scalar linear system is sufficient; P is reconstructed from this exact model relation. This does not model unequal diffusivities, a consumed coupling partner, proton balance, gas formation, supporting ions or side products.

At the inlet the prescribed **total** inward flux is u C_in (Danckwerts condition); this is deliberately a different boundary specification from an inlet Dirichlet concentration when axial diffusion is appreciable. At the outlet the axial diffusive flux is zero and convection leaves the domain. The upper wall has no species flux. At the lower wall,

$$J=k_f C_{A,w}-k_r C_{P,w},\qquad j=nFJ.$$

For a first cell center a distance h = Δy/2 from the wall,

$$J=\frac{(k_f+k_r)C_{A,c}-k_r C_T}{1+(k_f+k_r)h/D},\qquad C_{A,w}=C_{A,c}-Jh/D.$$

The zero-diffusivity diagnostic is defined by the D → 0 wall-access limit, giving zero flux; the uncoupled wall composition is assigned its effective equilibrium value. This is a numerical limiting case, not an experimentally accessible electrolyte model.

The default rate parameterization is

$$k_f=\frac{j_{scale}}{nFC_{ref}}\exp\left(\frac{\alpha F\eta}{RT}\right),\qquad k_r=\frac{j_{scale}}{nFC_{ref}}\exp\left(-\frac{(1-\alpha)F\eta}{RT}\right).$$

Here η is an **externally imposed potential index relative to a fixed reference state**. The exponent retains an effective one-electron transfer dependence from the attachment; n = 2 is used for molecular-to-charge accounting. This is an empirical rate law, **not a thermodynamically complete two-electron Butler–Volmer/Nernst model**. Neither electrode/electrolyte potential is solved, and none of j_scale, α or η has been calibrated to a particular molecule or electrode. This distinction matters when interpreting the otherwise exact charge closure. Concentration-dependent forward and reverse rates and reaction-current/stoichiometric-flux coupling are discussed in the [COMSOL primary modeling documentation](https://doc.comsol.com/6.3/doc/com.comsol.help.battery/battery_ug_modeling.05.08.html); that reference does not validate the numerical constants assumed here.

The finite-volume construction uses shared face fluxes, first-order upwind convection, central diffusion and a direct sparse linear solve. The Poiseuille speed is integrated over each transverse strip, so the discrete volumetric flow equals the imposed Q on every grid. The method follows the conservative cell and face-flux framework described in the [NIST FiPy documentation](https://pages.nist.gov/fipy/en/latest/numerical/discret.html); FiPy itself is not required or executed by this implementation.

中文要点：入口指定总物料通量，出口轴向扩散通量为零；壁面逆反应显式依赖产物浓度。二电子计量保证分子通量与电流的记账一致，但动力学指数仍是未标定的有效电位参数，不能当作完整二电子电化学热力学模型。没有求解温度场、电位场、迁移、压力、气泡或暂态。

## 3. Executed calculations / 已执行计算

Default inputs: L = 0.06 m, H = 0.0003 m, W = 0.012 m, D = 1.1 × 10^-9 m²/s, T = 298.15 K, Q = 450 µL/min, C_in,A = 50 mM, C_in,P = 0, η = 0.48 V as a model index, j_scale = 0.08 A/m², α = 0.5, C_ref = 50 mM, n = 2. The reactor volume is 0.216 mL, mean speed 0.00208333 m/s and nominal residence time 28.8 s. One mol/m³ equals one mM; no extra factor of 1,000 is needed for this concentration conversion.

At 100 × 48 cells the executed model gives conversion **58.3872858750%**, current **42.2513750193 mA**, outlet A concentration **20.8063570625 mM** and net turnover **2.18952322031 × 10^-7 mol/s**. The relative linear residual is 5.83 × 10^-15; feed-scaled material and charge balance errors are about 1.59 × 10^-14. These many digits expose reproducibility; they are not a statement of physical predictive precision.

| Joint grid | Cells | Model conversion / % |
|---|---:|---:|
| 50 × 12 | 600 | 58.2728316193 |
| 100 × 24 | 2,400 | 58.3900166910 |
| 200 × 48 | 9,600 | 58.4517308193 |
| 400 × 96 | 38,400 | 58.4833547496 |

The final joint refinement changes conversion by **0.0316239302 percentage points**. Axial and transverse refinements were also varied separately. The 80 × 32 screening grid gives 58.3562823302%, which differs from the finest tested grid by -0.1270724193 percentage points at the default condition. The finest tested grid is a numerical reference, not an exact solution.

First-order upwind has numerical diffusion. At the baseline centerline the estimate u Δx/2 is roughly 850 times the physical axial diffusivity on the 100-cell axial mesh, so merely including axial diffusion does not imply fine resolution of its small physical contribution. The refinement study quantifies conversion sensitivity; it does not establish all local gradients as grid independent. Removing physical axial diffusion on the 80 × 32 grid changes conversion by approximately 0.0001980 percentage points at this high-Péclet default condition.

The candidate-domain check compares 80 × 32 and 400 × 96 at the four corners Q = 100 or 1,500 µL/min and η = 0.20 or 0.75 V, plus the center Q = 800 µL/min and η = 0.475 V. The largest observed conversion difference is **0.1446419170 percentage points** at Q = 100 and η = 0.75. This covers five explicit points; it is not a bound on every point in the 81-candidate optimization pool.

There are **39 transport solves** in the verification run: one baseline, seven grid cases, 12 flow/potential sweep cases, ten solves for five domain points at two grids, five limiting cases, and four independent analytic-limit comparisons. All passed the declared numerical checks. The largest feed-scaled material residual is 3.0571 × 10^-10, occurring in the deliberately very high transverse diffusivity analytic-limit diagnostic; it is larger than the ordinary-case residual because that matrix is more ill-conditioned.

Zero wall kinetics, zero diffusion and equal A/P feed at equal forward/reverse rates preserve the expected constant state and zero net current to floating-point accuracy. A product-only feed at negative η gives a **negative current (-42.2289396144 mA on 80 × 32)**, retained in the output, and A-based conversion is `null` because the inlet A amount is zero. Trace negative concentrations around -2.2 × 10^-13 mM in a zero-reaction solve are stored without clipping and treated as floating-point roundoff against an explicit tolerance of 10^-10 C_T.

The independent plug-flow diagnostic imposes a uniform speed, very rapid transverse diffusion (10^-4 m²/s), zero axial diffusion and irreversible k_f = 10^-6 m/s. The cross-sectionally mixed continuum limit predicts X = 1 - exp(-k_f LW/Q) = **9.1535983931%**. Refining 40 → 80 → 160 → 320 axial cells reduces the conversion discrepancy approximately by half at each step; the last discrepancy is **-0.0013167311 percentage points**. This confirms the expected first-order axial discretization behavior in a separate limit, not actual reactor accuracy.

中文结果：默认算例转化率 58.3873%，电流 42.2514 mA；最细网格转化率 58.4834%，最后一次联合加密变化 0.03162 个百分点。筛选域五个检查点的最大网格差为 0.14464 个百分点。所有 39 次求解均满足所设数值检验，但不存在真实实验验证；单反应模型的 FE 恒为 100% 是结构性假设，不能解释为优化成功。

## 4. Reproducibility and result schema / 复现与接口

From the repository root, using the existing chemistry environment with its numerical DLL path configured and OMP/BLAS/MKL thread counts set to one:

```text
python electratwin/scripts/transport_reviewed.py
python -m unittest discover -s tests -p test_electratwin_transport.py -v
```

Import `solve_transport` from [the reviewed module](scripts/transport_reviewed.py). The returned object uses `schema_version = electratwin.transport.v1` and contains `model`, `evidence_class`, `assumptions`, `inputs`, `summary` and optional `field` (`include_field=False` omits the arrays). `summary` includes `conversion_pct`, signed `current_A`, `reacted_mol_s`, `outlet_A_mM`, `outlet_P_mM`, `material_balance_relative_error`, `charge_balance_relative_error`, `linear_residual_relative`, `nonnegative`, `converged` and elapsed time. Geometry is in `inputs`; volume and residence time are in `summary`. `field.A_mol_m3` and `field.P_mol_m3` have shape [nx, ny] at cell centers, with x coordinates in mm and y coordinates in µm. Wall profiles have nx entries. No endpoint interpolation or inlet-node trapezoid is used for current.

Results: [baseline field](results/transport/baseline_solution.json), [verification summary](results/transport/verification.json), [grid study](results/transport/grid_convergence.csv), [parameter sweep](results/transport/parameter_sweep.csv), [domain grid check](results/transport/screening_domain_grid_check.csv), [limiting cases](results/transport/limiting_cases.csv), [analytic plug-flow comparison](results/transport/analytic_plug_limit.csv), [15 independent tests](../tests/test_electratwin_transport.py), and [test transcript](results/transport/test_results.txt).

Each parameter-sweep row is an evaluated model condition with full-precision summary fields; it is not an experimental observation. Runtime is an observed property of this local CPU run and is not a guaranteed industrial latency. The physical model and reported constants require external calibration and validation before any quantitative real-reactor interpretation.
