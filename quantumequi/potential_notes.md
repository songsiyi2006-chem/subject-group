# Captured EGNN mathematical audit / 实际保存 EGNN 的数学审计

The unchanged source creates a 24,385-parameter, three-layer EGNN on eight supplied atoms (C2HCuN4). No model fitting, energy-force training set, chemical charge state or multiplicity is supplied. Its energy/force labels are nominal kcal/mol and kcal/(mol Å); the labels do not calibrate the neural function.

原始代码构造含 24,385 个参数的三层 EGNN，输入为八个手工坐标（C2HCuN4）。没有拟合、能量—力训练集、化学电荷态或多重度。源码的 kcal/mol 与 kcal/(mol Å) 是赋予输出的名义单位，不能自动校准随机神经函数。

## Implementation and probes

`scripts/potential_audit.py` loads `results/original/untrained_source_weights.pt` with `weights_only=True`. Float32 tests reproduce the captured reactant energy exactly. Float64 tests convert these same quantized parameter values to float64; they do not reinitialize the model or recover precision absent from the original weights.

Three geometries are checked: original reactant, original product, and the source-selected post-update candidate. For each of two dtypes, four deterministic QR transformations (seeds 101, 202, 303, 404) are applied with both determinant signs, translation (3.1, -2.2, 0.7) Å and deterministic atom permutations. For row-vector coordinates, the expected transformations are X'=XQᵀ+t and F'=FQᵀ, with the same atomic permutation applied to Z, X and F. The tested total energy must be invariant.

脚本只读取实际保存权重。单精度反应物能量与原始输出完全相同。双精度把同一组单精度权重的数值提升为 float64，不重新随机初始化。三个几何、两种精度、四个 QR 随机种子、正反两种行列式，加上平移与原子置换，共进行 48 次组合变换检验。总力和关于几何中心的总力矩也独立记录。

| Probe | Float32 | Float64 |
|---|---:|---:|
| Maximum energy covariance error | 4.76837158203125e-7 | 4.440892098500626e-16 |
| Maximum force covariance error | 1.3102842418705185e-8 | 3.469446951953614e-17 |

For every coordinate at each geometry, central energy differences use steps 10⁻², 10⁻³, 10⁻⁴ and 10⁻⁵ Å. Each full check calls the potential 48 times. At 10⁻⁴ Å, the largest float64 force discrepancy across the three geometries is 1.3121292247175731e-11; float32 cancellation instead reaches 0.0030195973813533783 in nominal force units. At 10⁻⁵ Å the largest float32 discrepancy reaches 0.016164911445230246. These tests compare the derivative of the same neural function, not a quantum-mechanical reference force.

对每个几何的全部 24 个笛卡尔分量，分别采用四个位移步长进行中心差分；每项完整检查需 48 次势函数调用。在 10⁻⁴ Å 步长下，双精度最大力差为 1.3121292247175731e-11；单精度抵消误差最大为 0.0030195973813533783。步长继续降至 10⁻⁵ Å 时，单精度最大误差增至 0.016164911445230246。这是同一函数的微分验证，不是量子力准确性验证。

## Disconnected final coordinate head

The final layer updates coordinates that are never consumed by the energy head. Automatic differentiation marks its three `phi_x` parameter tensors (1,024 + 32 + 32 = 1,088 parameters) disconnected from total energy. Adding 7 to every one of those parameters changes neither the reactant energy nor its forces, exactly in the double-precision check. Earlier coordinate heads do affect subsequent messages. This finding identifies unused energy-model parameters; it does not invalidate the tested E(3) covariance.

最后一层更新后的坐标没有进入能量头，因此最后一个 `phi_x` 的 1,088 个参数与能量、力无连接。将这些参数全部加 7 后，能量和力均保持不变。前两层坐标更新仍能影响后续消息。此结果说明存在冗余参数，并不否定已经检验的等变性。

## Workload and reproduction

Main audit: 48 transformation probes, 24 full Cartesian finite-difference checks, 1,209 energy/force evaluations, zero training epochs. The latter count includes six geometry/dtype baselines, one parameter-connectivity evaluation and two ablation evaluations. Four unit tests exercise covariance, force/torque balance, differentiation and atom relabeling; their calls are excluded from the study count.

主要审计包含 48 次变换、24 次全坐标差分及 1,209 次能量/力调用；没有训练。四项测试的调用另计，不纳入主计算量。

```bash
python quantumequi/scripts/potential_audit.py
python -m unittest discover -s tests -p "test_quantumequi_potential.py" -v
```

CSV files retain every transformation error, difference step and baseline force component. `results/potential/summary.json` links source, checkpoint, script and table hashes. The mathematical architecture background is [Satorras et al., 2021](https://proceedings.mlr.press/v139/satorras21a.html). Source absence of a cutoff and chemical calibration remains outside these passing mathematical checks.
