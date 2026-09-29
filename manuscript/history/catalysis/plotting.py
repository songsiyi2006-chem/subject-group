"""Five bilingual historical-audit figures; reads frozen extracts, no science jobs."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
SOURCES = HERE / "sources"
OUT = HERE / "figures"
BLUE, ORANGE, GREY = "#176b8e", "#bd6c16", "#6f7480"


def rows(name):
    with (SOURCES / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "text.color": "black", "axes.labelcolor": "black", "xtick.color": "black", "ytick.color": "black", "figure.facecolor": "white", "axes.facecolor": "white", "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none", "axes.unicode_minus": True})
    manifest = []
    gates = rows("posthoc_spin_gates.csv")
    spin = [r for r in rows("C04_spin_states.csv") if r["status"] == "converged"]
    ml = rows("posthoc_egnn_metrics.csv")
    neb = rows("posthoc_proton_wire.csv")
    xtb_values = json.loads((SOURCES / "C20_excerpt.json").read_text(encoding="utf-8"))["selected_values"]
    dft_values = json.loads((SOURCES / "C21_excerpt.json").read_text(encoding="utf-8"))["selected_values"]
    frequencies = []
    for run_id in ("P009-neutral-gfn2-gas-hess", "R001-rotor-hess"):
        prefix = next(k.removesuffix("/id") for k, v in xtb_values.items() if k.endswith("/id") and v == run_id)
        frequencies.append(xtb_values[prefix + "/minimum_vibrational_frequency_cm1"])
    dft_counts = [sum(r["status"] != "resource_preflight_not_launched" for r in dft_values["/runs"]), dft_values["/completed_quantum_runs"], dft_values["/failed_or_terminated_runs"], dft_values["/preflight_rejections_not_counted_as_quantum_runs"], len(dft_values["/energy_pairs"])]

    for language in ("english", "chinese"):
        cn = language == "chinese"
        plt.rcParams["font.family"] = ["Microsoft YaHei", "DejaVu Sans"] if cn else ["DejaVu Sans", "Microsoft YaHei"]

        def save(fig, stem, inputs):
            for ext in ("png", "svg"):
                path = OUT / f"{stem}_{language}.{ext}"
                fig.savefig(path, dpi=300, facecolor="white")
                manifest.append({"file": path.relative_to(HERE).as_posix(), "sha256": sha(path), "source_sha256": {f"sources/{p}": sha(SOURCES / p) for p in inputs}})
            plt.close(fig)

        fig, ax = plt.subplots(figsize=(9, 4.2), layout="constrained")
        labels = ["目标槽位", "实际状态计算", "SCF 收敛", "自旋与结构筛选通过", "可计算的原始垂直能隙", "合格诊断能隙", "接受的 MECP"] if cn else ["Target slots", "Physical state attempts", "SCF converged", "Spin + identity eligible", "Raw vertical gaps", "Eligible diagnostic gaps", "Accepted MECP"]
        y = np.arange(len(gates))
        vals = [int(r["count"]) for r in gates]
        ax.barh(y, vals, color=BLUE, height=.65)
        for i, val in enumerate(vals):
            ax.text(val + .7, i, str(val), va="center")
        ax.set_yticks(y, labels); ax.invert_yaxis(); ax.set_xlim(0, 60)
        ax.set_xlabel("记录数；阶段不是独立样本" if cn else "Record count; stages are not independent samples")
        ax.set_title("自旋筛选记录：有收敛不等于有可用交叉曲面" if cn else "Spin screening: convergence does not establish a crossing")
        save(fig, "Figure_C1_Spin_Gates", ["posthoc_spin_gates.csv"])

        fig, ax = plt.subplots(figsize=(9, 4.2), layout="constrained")
        names = []
        for r in spin:
            metal, rest = r["catalyst_id"].split("_", 1)
            family, sub = rest.rsplit("_", 1)
            names.append(f"{metal} {family.replace('_', '-')}/{sub}, M={r['multiplicity']}")
        vals = [float(r["spin_contamination"]) for r in spin]
        colors = [ORANGE if r["spin_contamination_flag"] == "True" else BLUE for r in spin]
        y = np.arange(len(spin))
        ax.barh(y, vals, color=colors, height=.65)
        for i, v in enumerate(vals):
            ax.text(max(v, 0) + .015, i, f"{v:.3f}" if abs(v) > 1e-9 else "~0", va="center", fontsize=9)
        ax.set_yticks(y, names, fontsize=8.5); ax.invert_yaxis(); ax.set_xlim(-.025, .8)
        ax.set_xlabel("〈S²〉−S(S+1)；橙色为原记录自旋筛选未通过" if cn else "〈S²〉−S(S+1); orange marks recorded spin-screen failures")
        ax.set_title("八个收敛状态中的自旋污染" if cn else "Spin contamination among eight converged states")
        save(fig, "Figure_C2_Spin_Quality", ["C04_spin_states.csv"])

        fig, ax = plt.subplots(figsize=(9, 4.2), layout="constrained")
        x = np.arange(len(ml)); width = .24
        for offset, key, label, color, hatch in [(-1, "model_MAE_eV", "辅助 EGNN" if cn else "Auxiliary EGNN", BLUE, None), (0, "equal_reference_MAE_eV", "等参考能基线" if cn else "Equal-reference baseline", ORANGE, "//"), (1, "training_mean_MAE_eV", "训练均值基线" if cn else "Training-mean baseline", GREY, "..")]:
            ax.bar(x + offset*width, [float(r[key]) for r in ml], width, label=label, color=color, hatch=hatch)
        ax.set_xticks(x, ["训练：374 对", "验证：80 对", "测试：86 对"] if cn else ["Train: 374 pairs", "Validation: 80 pairs", "Test: 86 pairs"])
        ax.set_ylabel("相对构象电子能 MAE / eV" if cn else "Relative conformer-energy MAE / eV")
        ax.set_ylim(0, .41); ax.legend(frameon=False, ncol=3, fontsize=8.7, loc="upper center")
        ax.set_title("辅助模型与简单基线；无势垒或 MECP 训练标签" if cn else "Auxiliary model versus baselines; no barrier or MECP labels")
        save(fig, "Figure_C3_EGNN_Baselines", ["posthoc_egnn_metrics.csv"])

        fig, axs = plt.subplots(1, 2, figsize=(9, 4.2), layout="constrained")
        labels = ["起点 1", "起点 2", "起点 1 续算"] if cn else ["Start 1", "Start 2", "Start 1\ncontinuation"]
        x = np.arange(3)
        ratios = [float(r["residual_over_target"]) for r in neb]
        axs[0].bar(x, ratios, color=ORANGE)
        axs[0].axhline(1, color="black", linestyle="--", linewidth=1)
        axs[0].set_ylabel("末次 NEB 力残差 / 目标" if cn else "Last NEB force residual / target")
        axs[0].set_title("a  未达到力收敛门槛" if cn else "a  Force threshold not reached")
        maxima = [float(r["band_maximum_eV"]) for r in neb]
        axs[1].bar(x, maxima, color=BLUE)
        axs[1].set_ylabel("采样路径最大势能差 / eV" if cn else "Sampled band maximum above reactant / eV")
        axs[1].set_title("b  不是活化自由能" if cn else "b  Not activation free energies")
        for ax, values in zip(axs, [ratios, maxima]):
            ax.set_xticks(x, labels)
            ax.set_ylim(0, max(values)*1.2)
            for i, val in enumerate(values): ax.text(i, val + max(values)*.025, f"{val:.2f}", ha="center", fontsize=9)
        save(fig, "Figure_C4_Unaccepted_Paths", ["posthoc_proton_wire.csv"])

        fig, axs = plt.subplots(1, 2, figsize=(9, 4.2), layout="constrained")
        axs[0].bar([0, 1], frequencies, color=[ORANGE, BLUE])
        axs[0].axhline(0, color="black", linewidth=.8)
        axs[0].set_xticks([0, 1], ["初始", "转子修复后"] if cn else ["Initial", "After rotor repair"])
        axs[0].set_ylim(-80, 110)
        axs[0].set_ylabel("最低振动频率 / cm⁻¹" if cn else "Lowest vibrational frequency / cm⁻¹")
        axs[0].set_title("a  xTB 分子驻点诊断" if cn else "a  xTB stationary-point check")
        for i, val in enumerate(frequencies): axs[0].text(i, val + (5 if val > 0 else -12), f"{val:.2f}", ha="center")
        vals = dft_counts
        labs = ["原生启动", "SCF 完成", "超时", "预检未启动", "完整电荷对"] if cn else ["Native launches", "SCF completed", "Timeouts", "Preflight, not launched", "Complete charge pairs"]
        axs[1].barh(np.arange(5), vals, color=BLUE, height=.65)
        axs[1].set_yticks(np.arange(5), labs, fontsize=8.5); axs[1].invert_yaxis(); axs[1].set_xlim(0, 4.8)
        for i, val in enumerate(vals): axs[1].text(val + .1, i, str(val), va="center")
        axs[1].set_xlabel("记录数；预检不计入原生启动" if cn else "Count; preflights excluded from launches")
        axs[1].set_title("b  定核 DFT；模型不含 Cu" if cn else "b  Fixed-nuclei DFT; no Cu")
        save(fig, "Figure_C5_Molecular_Preflight", ["C20_excerpt.json", "C21_excerpt.json", "C23_excerpt.json"])
    (OUT / "manifest.json").write_text(json.dumps({"generator_sha256": sha(Path(__file__)), "figure_groups": 5, "canvas_inches": [9, 4.2], "png_dpi": 300, "files": manifest, "scope": "Plots of fixed-commit historical records and lightweight arithmetic only."}, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"figure_files": len(manifest)}))


if __name__ == "__main__":
    main()
