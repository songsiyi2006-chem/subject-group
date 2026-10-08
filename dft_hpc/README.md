# DFT supplement: preparation, execution and evidence collection

## 中文使用说明

本模块把论文补算建议落实为四个固定几何单点任务：Q02（3-甲氧基吡啶母体）、Q03（3-氰基吡啶母体），各计算环氮质子化阳离子（+1，单重态）与加一个电子后的中性自由基（0，双重态）。方法为 PBE0/def2-TZVP，气相。这里不是未质子化吡啶的普通中性/阴离子配对。

两份 XYZ 来自已冻结的上游提交，字节哈希与原论文数据一致，来源见 [geometry_sources.json](geometry_sources.json)。原项目许可证及数据权利继续适用；本目录仅保存本次研究所需的两个原始坐标。Q01 已有较大基组比较，不重复列入本批。喹唑啉酮阴离子任务保留为待补全原始输入，尚未生成执行输入。

本版只支持上述四项电子单点：没有实现几何优化、频率、TS、IRC、MECP或溶剂电位预测。未知的波函数稳定性保持 `not_checked`；诊断能量配对不等于化学验收通过。

### 本地准备（不需要安装 Psi4）

从仓库根目录运行，目标目录必须不存在：

```sh
python -m dft_hpc.workflow prepare work/dft-campaign --cores 8 --memory-mb 16000 --hours 4 --concurrency 2
python -m dft_hpc.workflow collect work/dft-campaign/campaign.json --output work/dft-before-run.json
```

第一条生成可独立复制的目录，包括原始坐标、冻结任务表、运行器和 Slurm 脚本；第二条在尚未执行时应报告两个不完整配对，能量为 null。准备操作不会联网、登录或提交超算，不会覆盖已有目录。

### 超算执行

将整个生成目录上传至计算中心。请先根据该中心文档加载带 Psi4 的 Python 环境，并确认 `python -c "import psi4; print(psi4.__version__)"` 成功。在生成目录中执行：

```sh
# 先按本中心要求填写 account / partition；资源只是起始配置，不是耗时承诺。
sbatch --account=YOUR_ACCOUNT --partition=YOUR_CPU_PARTITION submit.slurm
```

任务数组默认最多同时运行两项。`--cpus-per-task` 与 Psi4 线程数一致，Psi4 内存预算取调度内存的80%，其余留给解释器和运行时。仅使用单节点 CPU；GPU、多节点MPI和软件环境加载需由计算中心另行配置。`--time` 仅由 Slurm 执行，直接运行 Python 时不会自动强制此时限。

每次执行写入新的 `runs/<job_id>/<attempt_uuid>/`，保留输入定义、原始输出及独立结果。成功返回只表示本次单点调用完成；被调度器强制终止的任务可能保持 started，收集器不会将其当成成功。没有自动重试。

```sh
python workflow.py collect campaign.json --output collected.json
```

收集器验证坐标、运行器、原始日志及任务身份。不同几何、方法、基组或环境不能混合配对；缺失、非有限能量和未完成SCF拒绝配对；同一状态有多个成功尝试时报告 ambiguous，不能任意选择能量较低的一次。应保留原目录，重新建立具有明确单次选择的新研究批次，而不是删除不喜欢的结果。

输出的差值为 **E(charge=0) − E(charge=+1)**，单位同时保存 Hartree 和 eV。它不是实验还原电位，也不是自由能。收集器检查记录一致性，不是独立重算原生输出或证明日志真实性；自旋诊断及波函数稳定性仍需科研审查。

### 测试与扩展

```sh
python -m unittest discover -s tests -p test_dft_hpc.py -v
```

测试里的数值是明确标记的合成样例，只验证软件行为。修改任务化学定义需要修改并复核源码、重新生成批次；本版故意拒绝手工篡改冻结任务清单。下一阶段可独立增加已验证的优化/频率流程，而不能把单点记录重新命名为优化结果。

## Saved local execution evidence / 已执行本地试跑

2026-10-08 使用已安装 Psi4 1.11、2 个 CPU 线程、2000 MB 调度预算（引擎使用80%）执行了 **1 项** Q02 阳离子 PBE0/def2-TZVP 固定几何单点。原生输出和返回值一致：电子能量为 **−362.89229397789865 Hartree**。UTC起止时间记录在 [结果文件](results/local_pilot_result.json)。其余三个任务未执行，完整电荷配对数为零，远程提交数为零。该试跑不验证优化、频率、波函数稳定性或化学预测准确度。

[原生日志可公开副本](results/local_pilot_engine_portable.out)仅移除了本机路径与主机名；它不是原始日志的逐字节副本。[来源记录](results/local_pilot_provenance.json)分别保留原始日志与公开副本的哈希，原始日志仍保留本机。运行 `python -m dft_hpc.validate_release` 可独立提取保存日志中的最终能量并验证来源与缺失配对状态。单元测试使用的合成样例不计入这项真实试跑。

## English scope

This bounded module prepares four archived-geometry PBE0/def2-TZVP gas-phase single points for Q02/Q03 cation/radical pairs. It does not submit to a scheduler. Preparation and collection use the Python standard library; explicit execution requires a site-provided Psi4 environment (including NumPy). Keep the complete prepared directory together. Geometry bytes and provenance are frozen, attempts are append-only, duplicate successful attempts remain ambiguous, and unavailable evidence stays unknown. Collected energy differences are diagnostic electronic quantities, not validated minima, free energies or electrode potentials. No new manuscript result is implied by generated inputs or synthetic tests.

## Documentation consulted

- [Slurm job arrays and bounded concurrency](https://slurm.schedmd.com/job_array.html)
- [Slurm sbatch resource and time options](https://slurm.schedmd.com/sbatch.html)
- [Psi4 runtime memory and threading](https://psi4.github.io/psi4docs/master/psithoninput.html)

These are implementation references, not evidence that this workflow has run on a remote cluster.
