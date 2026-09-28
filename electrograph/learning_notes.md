# Reviewed graph learning and frozen-pool search / 图学习与冻结候选池检验

## English

The supplied script has no supervised training loop, checkpoint, oxidation
labels, or atom-reactivity labels. Its output named oxidation potential comes
from random network weights. Its purported uncertainty is an independent
uniform random number rather than a GP posterior. Those outputs cannot validate
oxidation potentials, Fukui indices, or substrate reactivity. The separate
[reviewed implementation](scripts/graph_learning.py) replaces that claim with an
actual, explicitly different supervised **hydration free-energy benchmark**.

The local 642-row FreeSolv/SAMPL snapshot was matched by canonical isomeric
SMILES and experimental values to the official FreeSolv v0.52 snapshot. All
642 records agreed within 1e-8 kcal/mol; no discrepancies were excluded. The
first 256 canonical-structure SHA256 values define a label-blind subset.
Ring molecules are grouped by nonchiral Murcko scaffold; acyclic molecules are
grouped by full canonical structure. Whole groups are allocated in decreasing
size, with hashed tie order, to the largest remaining 60/20/20 deficit. This
gives 154 training, 51 validation, and 51 test molecules with disjoint structures
and group keys. There are 38/42/45 acyclic molecules in those partitions, and
62/50/50 groups. Acyclic close analogues may cross partitions: this is **not a
claim of general scaffold novelty**. Official experimental uncertainty and
reference fields, attribution, data transformations, full upstream CC-BY 4.0
license and its original-source caveat are retained in the
[data directory](data/learning/README.md).

The heavy-atom graph has 24 atom features: 13 element indicators including an
unknown bucket; charge, degree, H count, aromaticity, ring membership, mass and
isotope; and four absolute CIP categories. The CIP feature avoids the atom-order
dependence of raw clockwise/anticlockwise tags. There are 11 bond features:
five bond-type indicators including unknown, conjugation, ring membership, and
four stereo categories. Every bond produces two directed edges. Each of two
layers constructs a full 24x24 matrix from the edge features and sends
`W(e_ji) h_j / sqrt(24)` to the destination; destination messages are summed and
update the atom state through a GRU. Isolated atoms receive a zero-message
update. Sum and mean pooling are concatenated into a 48-dimensional molecular
embedding. The scalar head has 29,369 total trainable parameters. No atom-level
property head is retained without atom-level labels.

Targets are experimental hydration free energies in kcal/mol, standardized by
training-only mean −4.6901298701 and population SD 4.2638749810. Three independent
initialization/minibatch seeds (20260928–20260930) each run 60 Adam epochs with
learning rate 0.003, weight decay 1e-4, batch size 32 and gradient norm cap 5.
Checkpoints minimize validation MSE, at epochs 39, 44 and 60. Test labels are not
used for fitting, normalization or checkpoint choice. A separate control shuffles
only training labels (seed 20261001); it still uses the true validation labels
for checkpoint selection and selects epoch 2. This is a bounded negative control,
not a permutation significance test. A ten-descriptor ridge baseline scales
features on training rows only and selects alpha from {0.1, 1, 10, 100} by validation
RMSE; alpha 0.1 is selected. A training-mean baseline is also included.

| Model | Test RMSE (kcal/mol) | Test MAE (kcal/mol) |
|---|---:|---:|
| MPNN seed 20260928 | 1.897239 | 1.426569 |
| MPNN seed 20260929 | 1.746786 | 1.269662 |
| MPNN seed 20260930 | 1.249601 | 1.016765 |
| Descriptor ridge | 1.653425 | 1.202677 |
| Training mean | 4.805322 | 3.503397 |
| Shuffled-training-label MPNN | 4.852263 | 3.589107 |

The three MPNN RMSE values have mean 1.631209 and sample SD 0.338936. This SD
describes training-seed variation on the same split, not a generalization
confidence interval. Two of three MPNNs are worse than ridge; stable neural
superiority is not established. All 1,536 prediction rows (256 structures x six
models), 240 epoch-loss rows, exact split memberships, scalers, and four neural
weight checkpoints are saved.

The retrospective search maximizes **negative experimental hydration free
energy**, solely over the 51 test structures. It does not optimize yield or
oxidation. The learned feature encoder uses the predeclared seed 20260928,
independent of which seed performs best on the test set. Its feature scaler is
fit to training embeddings only. A competing GP uses the ten training-scaled
descriptors. Both GPs have a fixed unit-amplitude Matérn-5/2 kernel, Euclidean
distance divided by sqrt(feature count), length scale one, and a standardized
noise variance of 1e-5. At every step outcome mean/SD and posterior conditioning
use only already observed pool labels. UCB is `mean + 1.96 * latent_sd`; ties
choose the lowest available index. The SD is a genuine conditional GP posterior
quantity but is not calibrated experimental or physical uncertainty. Kernel
hyperparameters are not fit. This is a frozen neural embedding plus GP, not
jointly trained end-to-end deep kernel learning.

Eight seeds (0–7), three methods and 16 distinct queries each give 24 campaigns
and **384 cached label calls**. All methods share the same four initial points
within each seed. Encoder training consumes 154 historical labels, plus 51
validation labels for checkpoint selection; this additional historical label
cost is disclosed separately from test-pool calls. Both GP methods find the
pool minimum, −23.62 kcal/mol for polyhydroxy molecule mobley_4587267, within
nine calls for all eight seeds. At call 16 both have zero simple regret and tie
in every seed. Random search finds that optimum for 1/8 seeds and has mean regret
12.60625 kcal/mol. The paired GP-minus-random reduction is 12.60625 with descriptive
95% t interval [7.581861, 17.630639], seven wins and one tie. This single pool has
a readily distinguished polyhydroxy extreme; neither general search superiority
nor a neural-feature advantage over descriptors follows.

The full run took 45.60 s on the local CPU with one numerical thread, a local
observation rather than a hardware-general performance claim. A separate
two-epoch pilot used two neural fits, four epochs, six campaigns and 96 calls;
its random-order chirality feature preceded the final absolute-CIP correction.
Its original source bytes are preserved and hash-verified in
[the pilot snapshot](results/learning/pilot/graph_learning_pilot_snapshot.py.txt).
It is excluded from final comparisons. Across pilot and final run the actual
accounting is six neural fits, 244 epochs and 480 cached calls. No DFT or new
experiment was performed. Experimental label uncertainties are recorded but
are not modeled in the unweighted squared-error loss or in the GP noise term.

Sixteen fast tests pass: reciprocal edges and shapes; unknown/isolated atoms;
invalid input; atom permutation; graph batching; equivalent chiral SMILES;
finite edge gradients; canonical/group identity; no split overlap; label-blind
allocation; finite descriptors; GP conditioning against dense algebra; affine
output behavior; invalid GP state; query budgets/shared starts; first adaptive
selection independent of unobserved labels; and saved checkpoint reload.
Some tests cover multiple invariants. These are software checks, not chemical
validation. Use the [JSON summary](results/learning/summary.json),
[metrics](results/learning/regression_metrics.csv),
[query rows](results/learning/campaign_evaluations.csv), and
[test log](results/learning/test_results.txt) for exact values and provenance.

Reproduction after installing the repository's already declared scientific
dependencies: `python electrograph/scripts/graph_learning.py`. It uses the
packaged selected data without network access. The optional `--pilot` writes to
a separate pilot directory; it executes the current corrected implementation,
so it does not reproduce the archived pre-correction pilot byte-for-byte.
Rebuilding the selection from a local SAMPL file is supported with
`--prepare-from PATH_TO_SAMPL.csv`. No hardware interfaces are present.

Primary methodological references: [Gilmer et al., MPNN, 2017](https://proceedings.mlr.press/v70/gilmer17a.html),
[PyTorch GRUCell](https://docs.pytorch.org/docs/stable/generated/torch.nn.GRUCell.html),
and [scikit-learn data leakage guidance](https://scikit-learn.org/stable/common_pitfalls.html).

## 中文

原脚本没有监督训练循环、检查点、氧化电位标签或原子反应性标签；所谓氧化电位来自随机网络，
所谓不确定度来自独立均匀随机数。它们不能验证氧化电位、Fukui 指数或底物反应性。
本次另建的审查版实际训练了边条件 MPNN，但明确以 **水合自由能** 为监督目标，
没有把另一个性质的成功冒充为电化学验证。

复用本地 FreeSolv/SAMPL 的 642 条记录，并按规范异构 SMILES 与实验值逐条核对官方
v0.52 快照。642 条均在 1e-8 kcal/mol 容差内一致，无不一致记录被排除。按照规范结构
SHA256 排序取前 256 个分子，不按标签数值筛选。环状分子按无手性 Murcko 骨架分组，
无环分子按完整规范结构分组；组按大小递减、哈希破同序，分配到 60/20/20 目标中缺额最大
的部分，最终训练/验证/测试为 154/51/51，结构与组键均无交叉。三部分无环分子分别为
38/42/45 个，组数为 62/50/50。无环近似物仍可能跨集合，**不能宣称普遍骨架新颖性**。
原实验误差、文献字段、作者署名、完整 CC-BY 4.0 许可及上游原始来源条款均已保留。

模型用重原子图、24 维原子特征和 11 维键特征，含显式未知类别。手性用绝对 CIP 标签，
避免直接使用随原子顺序改变的顺/逆时针标签。每条键双向传递消息；两层边 MLP 各产生
24x24 变换矩阵，消息为 `W(e_ji) h_j / sqrt(24)`，求和后通过 GRU 更新。孤立原子接收
零消息。求和与均值池化组成 48 维表示；全模型共 29,369 个参数。没有原子标签，因此
不保留声称可预测原子反应性或 Fukui 指数的输出头。

实验水合自由能单位为 kcal/mol，仅用训练集均值 −4.6901298701、标准差 4.2638749810
做标准化。三个固定种子 20260928–20260930 各训练 60 轮，Adam 学习率 0.003、
权重衰减 1e-4、批量 32、梯度范数上限 5。按验证 MSE 最小选择第 39/44/60 轮检查点，
测试标签不参与拟合、标准化或检查点选择。种子 20261001 的负对照仅打乱训练标签，
仍用真实验证标签选择检查点，选择第 2 轮；它不是置换显著性检验。
十个 RDKit 描述符的岭回归只在训练集拟合缩放，验证集从 0.1/1/10/100 中选择 alpha=0.1。

三个 MPNN 的测试 RMSE 分别为 1.897239、1.746786、1.249601 kcal/mol；
岭回归 1.653425，训练均值 4.805322，置乱标签网络 4.852263。上表同时给出对应 MAE。
MPNN 的种子均值为 1.631209、样本标准差 0.338936；这反映同一拆分上的训练随机性，
不等于泛化置信区间。三个种子中有两个不如岭回归，因此没有证明神经网络稳定占优。
已保存 1,536 条预测、240 条轮次损失、逐分子拆分、缩放参数及四份权重检查点。

回顾性候选池搜索仅在 51 个测试分子上 **最大化负实验水合自由能**，不优化产率或氧化。
学习表征预先指定种子 20260928，不根据测试优劣挑选编码器；表征缩放仅在训练集拟合。
对照 GP 使用十个训练集缩放后的描述符。两者使用固定 Matérn-5/2 核，特征距离除以
特征维数的平方根，幅度和长度尺度均为 1，标准化噪声方差 1e-5。每一步仅根据已查询
标签更新输出均值/标准差与后验，UCB 为 `后验均值 + 1.96 × 潜在函数后验标准差`。
这是真实条件 GP 后验，不是随机构造的误差条；它也不是经校准的实验或物理不确定度。
核参数固定，编码器与 GP 没有端到端联合训练。

八个种子、三种方法、各 16 个不同查询，共 24 组搜索、**384 次缓存标签调用**。
同一种子三个方法共享四个初始点。学习表征还消耗 154 个历史训练标签和 51 个验证标签，
该历史成本与候选池查询数分开披露。两种 GP 均在九次调用内、八个种子全部找到最优
多羟基分子 mobley_4587267，其水合自由能 −23.62 kcal/mol。到第 16 次调用，两种
GP 的简单遗憾均为零、八次全部平局；随机策略仅 1/8 找到最优，平均遗憾为 12.60625
kcal/mol。GP 相对随机的配对改善为 12.60625，描述性 95% t 区间为
[7.581861, 17.630639]，七胜一平。这个固定池存在容易通过描述符区分的多羟基极值，
结果既不能证明普遍搜索优越性，也没有证明学习表征胜过描述符。

主运行在本地单线程 CPU 用时 45.60 秒，不是跨硬件性能保证。另存的两轮 pilot
包括两次网络训练、四轮、六组搜索和 96 次缓存调用；它发生在绝对 CIP 特征修复前，
原始源码字节已归档并核对哈希，不纳入最终比较。合计实际执行六次训练、244 轮和
480 次缓存调用。没有新实验或 DFT。原始实验误差被保留，但不参与无权平方误差损失，
也不等于 GP 中固定的噪声项。16 项快速测试覆盖图不变性、边梯度、拆分隔离、
GP 代数、预算、观测边界和权重重载；软件测试不等于化学验证。

复现命令：`python electrograph/scripts/graph_learning.py`。它读取随仓库分发的选定数据，
无需联网。`--pilot` 使用当前修正代码写入独立目录，不逐字复现旧版 pilot；
`--prepare-from PATH_TO_SAMPL.csv` 可从本地原表重建选样。精确参数、版本、数据与源码
哈希见 JSON 汇总，逐步查询与全部负结果保留在 CSV。
