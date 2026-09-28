# Equivariant geometry audit and supervised denoising / 等变几何审计与监督去噪

## English

This module separates two questions: whether energy/force transformations are
mathematically consistent, and whether a network actually learns a specified
denoising task. Neither the source potential nor the reviewed charge-constrained
potential is trained on molecular energies or forces. Their outputs have no
validated eV or chemical-force interpretation. The trained model instead predicts
clean **synthetic eight-node cycle coordinates**, in arbitrary length units.
There is no learned time-dependent score, forward diffusion schedule, or reverse
sampling chain; this is not a molecular diffusion generator. The primary EGNN
paper defines the equivariant construction, while DDPM makes clear the additional
probabilistic training and reverse-process structure absent here.
[EGNN](https://proceedings.mlr.press/v139/satorras21a.html),
[DDPM](https://proceedings.neurips.cc/paper/2020/hash/4c5bcfec8584af0d967f1ab10179ca4b-Abstract.html).

### Source defects and the exact captured-weight audit

The main diagnostic uses an explicitly recorded small source model
(`num_species=30`, hidden width 16, two layers, seed 991), alongside a separately
initialized reviewed model of width 24, two layers and prescribed total charge
zero. It is an architecture probe, not the exact weights from the supplied main
execution. A separate supplemental audit loads the **actual captured 55,394
source parameters**, embedding size 30, width 48 and three layers, together with
the eight-node initial coordinates and field from the compatibility run. Changing
embedding size from 10 to 30 changes RNG consumption; matching seeds alone would
not recover those weights. The supplemental weight, input, source and replay
script hashes are saved in [captured_source_audit.json](results/equivariant/captured_source_audit.json).

For the exact captured model, joint rotations/reflections of coordinates and the
external field preserve energy to 2.38419e-7 and transform forces to 5.70254e-9
maximum absolute error in float32. Permutation energy error is zero and force
error 3.72529e-9. Rotating coordinates while holding a laboratory field fixed is
a different physical configuration, so it is not the covariance test.

The captured raw charges sum to 2.80738759 without a fixed molecular-charge
constraint. Translating coordinates by (4, −3, 7), with the source field
(0, 0, 0.05), shifts energy by −0.982585669 and changes individual forces by
up to 0.001329373. With geometry-dependent total charge Q(R), the source obeys

\[
U(R+t)-U(R)=-(t\cdot E)Q(R),\qquad
F(R+t)-F(R)=(t\cdot E)\nabla_R Q(R).
\]

The force-identity residual is 9.34828e-9. The energy shift alone would be expected
for a *fixed* charged system; the geometry-dependent charge and consequent
origin-dependent internal forces are the defect. The actual chemical charge of
the incomplete typed fragment is not established by this test.

The source cutoff only zeros radial inputs, leaving node features and MLP biases
active. At a separation of 8, beyond the cutoff 6, the captured first layer still
changes node features by 0.037263751 and coordinates by 0.019687176 compared with
removing the edge. `.double()` fails because an aggregation tensor is created as
float32. Captured-source force finite-difference maximum errors at steps
0.01/0.001/0.0001 are 1.92267e-5/2.25190e-4/1.650996e-3: smaller steps aggravate
float32 subtraction. This is numerical differentiation evidence, not a force
accuracy estimate.

### Reviewed mathematical implementation

The reviewed layer multiplies **complete scalar messages and coordinate
coefficients** by the cosine cutoff, preserves dtype, and registers both radial
centers and widths as buffers. Scalar messages depend on invariant node features
and distances; vector increments are scalar-weighted relative coordinates. The
cutoff value and its first derivative vanish at the boundary. Edges beyond the
cutoff have zero effect in the numerical test. No claim of independent fragments
sharing a global molecular-charge constraint is made: the potential interface
represents one specified charged object.

The charge projection is

\[
q_i=\widetilde q_i-\frac1N\sum_j\widetilde q_j+\frac{Q_{\rm fixed}}N.
\]

The neutral numerical example uses Q_fixed=0. Its energy translation error is
1.11022e-16 and force error 1.64799e-17. For fixed nonzero charge, a separate test
checks the expected constant energy shift and total force Q_fixed E, while
individual forces remain origin invariant. Projection fixes the mathematical
charge constraint, not the accuracy of partial charges. The small reviewed
double-precision model has finite-difference maximum error 5.80916e-12 at step
1e-4. Its architecture differs from the source; this is not a controlled
same-weights precision comparison. All original defects and error values remain
archived.

### Actual synthetic training

The dataset contains 368 distinct clean shapes: 192 training, 48 validation, 64
test and 64 held-out warped cycles. The first three partitions are planar
ellipses with radius uniformly 0.8–1.6, axis ratio 0.75–1.25, random phase, proper
rotations or reflections, and translations. The OOD set adds a second-harmonic
out-of-plane displacement of amplitude 0.25–0.45. Generator seeds are
33001/33002/33003/33004. Each shape has eight ordered cycle vertices. Independent
Gaussian noise has nominal sigma 0.08/0.16/0.24; its graph mean is subtracted,
so centroid noise is absent and per-component variance is reduced by 7/8.
This deliberate simplification makes centroid recovery unnecessary. A shape
appears once; no augmented copy crosses partitions.

The 13,176-parameter network has width 24, three equivariant layers, and inputs
the noise level and known cycle connectivity. Coordinate updates are centered
per graph. Arbitrary node relabeling must also relabel adjacency; the public
interface accepts an explicit adjacency for this purpose. Centering updates does
not turn denoising into physical dynamics. The network uses Adam, learning rate
0.003, weight decay 1e-5, batches of 32, gradient norm cap 5. Seeds 4441/4442/4443
each run 60 epochs; validation loss alone selects epochs 57/60/60. Test and OOD
outcomes do not choose checkpoints.

Baselines are unchanged noisy coordinates; per-input PCA projection onto a best
plane; and a cycle-neighbor smoother with one training-fitted scalar per noise
level. The latter has coefficients 0.130937054/0.303102639/0.427779964. PCA has a
strong planarity prior; it is deliberately tested on warped shapes where that
prior is wrong. Metrics are raw coordinate RMSE without post-hoc alignment and
RMSE over the 28 unordered inter-node distances per shape.

| Model | Test coordinate RMSE | OOD coordinate RMSE | Test distance RMSE | OOD distance RMSE |
|---|---:|---:|---:|---:|
| Identity | 0.162532 | 0.164768 | 0.250556 | 0.241461 |
| PCA plane | 0.146258 | 0.203743 | 0.250818 | 0.252519 |
| Training-fitted cycle smoother | 0.126434 | 0.132234 | 0.233569 | 0.232804 |
| EGNN seed 4441 | 0.093198 | 0.138630 | 0.149020 | 0.152178 |
| EGNN seed 4442 | 0.096786 | 0.126767 | 0.156405 | 0.155963 |
| EGNN seed 4443 | 0.097257 | 0.124032 | 0.158868 | 0.159351 |

The network improves all reported in-distribution coordinate errors, but seed
4441 is worse than the simple smoother on OOD coordinate recovery. PCA becomes
worse than doing nothing on that OOD set. Different metrics need not rank models
identically; a lower distance error does not imply better absolute-coordinate
recovery. The three learned models' joint reflection/rotation/translation checks
on four test shapes, evaluated after a double-precision copy, have maximum errors
1.77636e-15/8.88178e-16/8.88178e-16. These are numerical covariance checks, not
empirical confidence intervals or molecular generative performance.

### Costs, tests and reproduction

The frozen main study trains three models for 180 epochs, has 2,944 noisy nodes,
24 aggregate metric rows, 2,208 per-shape metric rows and 17,664 predicted node
coordinate rows. The same synthetic dataset is reused by a separate two-epoch,
one-model pilot; it is not 736 unique shapes. The pilot source predates the
optional-adjacency API and extra input checks; its exact hash-matched bytes are
retained in [the pilot snapshot](results/equivariant/pilot/equivariant_pilot_source.py.txt).
Pilot and main use identical default model computation; the pilot is excluded
from reported final metrics. Main local runtime was 25.75 s, with PyTorch and
BLAS/OpenMP threads fixed to one; this is not a portable performance guarantee.
Versions are Python 3.12.14, NumPy 2.4.6 and PyTorch 2.10.0.

The main architecture audit separately makes 123 source and 121 reviewed
energy/force evaluations, 124 source charge-only forwards, four extra layer
forwards and one expected dtype-failure probe. The pilot repeats those audit
counts. The exact-captured-weight supplement adds 157 energy/force evaluations,
one charge-only forward, two layer forwards and one expected dtype failure,
with **zero additional training, MD or NEB runs**. These categories describe
public function calls, not all internal tensor operations. Unit tests are
excluded. The supplement has its own manifest so the primary frozen training
source and result hashes remain intact.

Eighteen fast tests pass, covering cutoff and derivative, whole-message gating,
dtype/gradients, O(3) and translation covariance, graph relabeling, charge/force
translation laws, finite differences, input rejection, shape uniqueness,
centroid convention, baselines, checkpoint replay and output hashes. No full
training occurs during tests. Reproduce from the repository root:

```text
python synthapore/scripts/equivariant_reviewed.py --pilot
python synthapore/scripts/equivariant_reviewed.py
python synthapore/results/equivariant/audit_captured_source.py
python -m unittest discover -s tests -p test_synthapore_equivariant.py -v
```

The supplemental audit requires the separately captured source weights and
initial inputs. Current `--pilot` uses the current optional-adjacency interface;
the historical snapshot is retained to explain its original source hash.
See [summary](results/equivariant/summary.json),
[metrics](results/equivariant/denoising_metrics.csv),
[source/reviewed audit](results/equivariant/source_and_reviewed_audit.json), and
[test log](results/equivariant/test_results.txt). No download, DFT, instrument or
new experimental measurement was used in this module.

## 中文

本模块分别回答两个问题：能量和力的变换是否在数学上自洽，以及网络是否真正学会了
指定的去噪任务。源势函数和审查版电荷约束势均未用分子能量或力训练，其数值不能作为
经验证的 eV 能量或化学力。实际训练的模型预测**合成八节点环的无噪坐标**，单位为
任意长度。没有时变分数函数、正向扩散日程或反向采样链，因此不称为分子扩散生成器。
EGNN 原论文支持等变构造；DDPM 原论文则说明本代码尚未实现的概率训练和反向过程要求。

### 原模型缺陷与精确捕获权重审计

主诊断使用显式记录的小型源模型（num_species=30、宽度16、两层、种子991），以及
另行初始化的宽度24、两层、总电荷指定为零的审查版模型。这是架构探针，不是原主程序
的同一组权重。补充审计另加载兼容执行中**实际保存的55,394个源参数**：嵌入大小30、
宽度48、三层，并使用其八节点初始坐标和外场。嵌入大小由10改为30会改变随机数消耗，
仅重新设置相同种子不能保证恢复这些权重。权重、输入、源文件和重放脚本哈希另存。

精确捕获模型的坐标与外场共同旋转或反射时，float32 能量误差最大2.38419e-7，
力协变误差最大5.70254e-9；置换能量误差为零、力误差3.72529e-9。若只旋转坐标，
保持实验室外场方向不变，则物理配置改变，不能要求同样的协变关系。

捕获模型的原始电荷总和为2.80738759，没有固定分子总电荷约束。保持原外场
(0,0,0.05)，平移(4,−3,7)使能量变化−0.982585669，单原子力最大改变0.001329373。
当Q依赖几何时，上述能量平移式和力平移式给出附加项，后者数值残差9.34828e-9。
固定带电体系发生常数能量平移本身并非错误；这里的问题是总电荷依赖几何，进而使
内部力依赖坐标原点。本检验没有确定这个不完整元素片段的真实化学电荷。

原截止函数仅令径向输入归零，节点特征和MLP偏置仍起作用。距离8大于截止6时，
捕获模型第一层相对移除该边，仍改变标量特征0.037263751、坐标0.019687176。
double转换因聚合张量硬编码float32而失败。精确捕获模型在差分步长
0.01/0.001/0.0001下，力最大误差分别为1.92267e-5/2.25190e-4/1.650996e-3，
缩小步长反而放大float32相减误差。这是数值微分检查，不是化学力准确度评估。

### 审查版数学实现

审查版对**完整标量消息和坐标系数**施加余弦截止，保留dtype，并把径向中心与宽度
均注册为buffer。消息仅依赖标量节点特征和距离；向量更新为相对坐标乘标量。
截止边界值和一阶导数为零，超截止边在检验中无影响。电荷投影遵循上式：原始电荷
减去均值，再加Q_fixed/N。接口描述一个总电荷指定的对象，不声称在共享全局电荷
约束时可以把多个互不关联片段当作独立体系。

中性数值例取Q_fixed=0，平移能量误差1.11022e-16、力误差1.64799e-17。
另对固定非零电荷检查预期常数能量平移和总力Q_fixed E，单原子力仍应与原点无关。
投影修复总电荷数学约束，不验证局部电荷。小型审查版双精度模型在步长1e-4下的
最大力差分误差为5.80916e-12。它与源模型架构不同，不能称为同权重条件下的纯精度
比较。原始缺陷与结果均保留。

### 实际合成数据训练

数据含368个不同无噪形状：训练192、验证48、测试64、分布外扭曲环64。
前三部分是平面椭圆，半径0.8–1.6、轴比0.75–1.25，并加入随机相位、旋转/反射和平移。
分布外环加入幅度0.25–0.45的二次谐波面外位移。生成器种子为
33001/33002/33003/33004。每个形状有八个按环顺序排列的节点；高斯噪声名义标准差
0.08/0.16/0.24，随后减去图内均值，因此质心噪声为零、单分量方差减为原来的7/8。
该简化免去了恢复未知质心的困难。每个干净形状只出现一次，没有增强副本跨拆分。

网络宽度24、三层，共13,176参数，输入已知噪声级和环连接关系。每个图的坐标更新
去均值。任意节点重标号必须同时重标号邻接；公开接口允许传入对应邻接矩阵。
去中心化更新不等于物理动力学。Adam学习率0.003、权重衰减1e-5、批量32、
梯度范数上限5。4441/4442/4443三个种子各训练60轮，仅依据验证损失选取
57/60/60轮，测试和分布外表现不参与选择。

基线包括原样保留噪声坐标、按输入做最佳平面PCA投影，以及每个噪声级仅用训练数据
拟合一个系数的环邻居平滑器。三个平滑系数为0.130937054/0.303102639/0.427779964。
PCA使用很强的平面先验，故意在违反该先验的扭曲形状上检验。指标为未作事后对齐的
原坐标RMSE，以及每个形状28对无序节点间距离的RMSE；完整数值见上表。

网络改善了所有列出的同分布坐标误差，但种子4441在分布外坐标恢复上不如简单平滑器。
PCA在分布外甚至不如原样保留噪声。坐标和距离指标可能有不同排名，较低距离误差
不必然代表绝对坐标恢复更好。把三个训练模型转为双精度后，对四个测试形状作共同
反射/旋转/平移，最大误差为1.77636e-15/8.88178e-16/8.88178e-16。这些是数值
协变性检查，不是经验置信区间或分子生成性能。

### 计数、测试与复现

主研究实际训练三个模型、180轮，有2,944个带噪节点、24行汇总指标、2,208行逐形状
指标、17,664行预测节点坐标。独立pilot是一模型两轮，复用同一数据，不能说有736个
唯一形状。pilot早于可选邻接API和额外输入校验，精确源字节及匹配哈希已保存。
默认模型计算相同，pilot不进入最终结果比较。主计算本地耗时25.75秒；PyTorch与
BLAS/OpenMP均为单线程，不构成跨硬件性能保证。版本为Python3.12.14、NumPy2.4.6、
PyTorch2.10.0。

主架构审计另有123次源模型、121次审查版能量/力调用，124次源电荷专用前向，
四次额外层调用及一次预期dtype失败。pilot重复这些审计数。精确捕获权重的补充审计
再增加157次能量/力调用、一次电荷前向、两次层调用及一次预期dtype失败，
**没有新增训练、MD轨迹或NEB运行**。计数指公开函数调用，不是所有内部张量运算，
也不含单元测试。补充审计使用独立清单，保留主训练源代码和结果哈希。

18项快速测试通过，覆盖截止与导数、完整消息屏蔽、dtype/梯度、O(3)和平移协变、
图重标号、电荷/力平移律、有限差分、输入拒绝、形状唯一性、质心约定、基线、权重
重放及哈希。测试期间不训练。复现命令见英文部分；补充审计要求已存在的源权重与初始
输入。当前pilot使用当前邻接接口，历史快照解释原始哈希。本模块没有下载、DFT、
仪器操作或新实验测量。
