"""Check saved ElectroGraph provenance and publication files without reruns.

Passing checks establish artifact consistency, not physical validation.
No models are trained and no chemistry or stochastic simulations are executed.
"""
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import sys
from urllib.parse import unquote
import xml.etree.ElementTree as ET
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
RESULT = ROOT / "results/publication_validation.json"
STATES = ("empty", "substrate", "radical", "intermediate", "product")
STEMS = ("Fig1_Source_Audit", "Fig2_Graph_Benchmark", "Fig3_Conformer_Ensemble", "Fig4_Kinetic_Validation")


def main():
    checks = []

    def check(name, ok, detail=None):
        item = {"check": name, "passed": bool(ok)}
        if detail is not None:
            value = str(detail)
            for path, label in ((ROOT, "electrograph"), (REPO, "repository")):
                value = value.replace(str(path).replace("\\", "\\\\"), label).replace(str(path), label)
            item["detail"] = value
        checks.append(item)

    def guarded(name, function):
        try:
            function()
        except Exception as exc:
            check(name + " readable and valid schema", False, f"{type(exc).__name__}: {exc}")

    def read(name):
        return json.loads((ROOT / name).read_text(encoding="utf-8-sig"))

    def rows(name, expected=None):
        with (ROOT / name).open(encoding="utf-8-sig", newline="") as handle:
            values = list(csv.DictReader(handle))
        if expected is not None:
            check(str(name) + " row count", len(values) == expected)
        return values

    def digest(path):
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()

    def hashes(mapping, base, label):
        check(label + " nonempty hash inventory", bool(mapping))
        for name, expected in mapping.items():
            path = base / name
            check(label + " " + name, path.is_file() and digest(path) == expected)

    def near(a, b, tol=1e-10):
        return math.isclose(float(a), float(b), rel_tol=tol, abs_tol=tol)

    def module_provenance():
        for module in ("learning", "kinetics", "structure"):
            for prefix in ("", "pilot/"):
                name = f"results/{module}/{prefix}summary.json"
                def validate_summary(name=name, module=module, prefix=prefix):
                    summary = read(name)
                    if module == "learning" and prefix == "pilot/" and "source_snapshot" in summary:
                        snapshot = summary["source_snapshot"]
                        path = ROOT / f"results/{module}/{prefix}" / snapshot["path"]
                        check(name + " archived pilot source identity",
                              summary["pilot"] and len(summary["source_sha256"]) == 1 and
                              set(summary["source_sha256"].values()) == {snapshot["sha256"]} and
                              path.is_file() and digest(path) == snapshot["sha256"] and bool(snapshot["reason"]))
                    else:
                        hashes(summary["source_sha256"], REPO, name + " source")
                    hashes(summary["output_sha256"], ROOT / f"results/{module}/{prefix}", name + " output")
                guarded(name, validate_summary)

    def source():
        record = read("source/source_record.json")
        original = read("results/original/execution.json")
        compatibility = read("results/compatibility/execution.json")
        attempt = read("results/api_repair_attempt_1/execution.json")
        audit = read("results/source_audit.json")
        check("archived specification bytes", digest(ROOT / "source/specification.md") == record["specification_sha256"])
        check("original source execution and audit hashes", digest(ROOT / "source/electrograph_kmc_core.py") ==
              record["original_python_sha256"] == original["source_sha256"] == audit["source_sha256"])
        check("compatibility source execution and audit hashes", digest(ROOT / "source/electrograph_kmc_core_compat.py") ==
              record["compatibility_python_sha256"] == compatibility["source_sha256"] == audit["compatibility_sha256"])
        check("failed intermediate source and execution retained", digest(ROOT / "source/electrograph_kmc_core_api_attempt_1.py") ==
              record["archived_failed_intermediate_source_sha256"] == attempt["source_sha256"] and attempt["exit_code"] == 1)
        check("original failure and repaired completion retained", original["exit_code"] == 1 and compatibility["exit_code"] == 0
              and not original["timed_out"] and not compatibility["timed_out"])
        check("no source physical instrument connections", all(not r["physical_instrument_connected"] for r in (original, compatibility, attempt)))
        first = (ROOT / "source/electrograph_kmc_core.py").read_bytes()
        second = (ROOT / "source/electrograph_kmc_core_compat.py").read_bytes()
        newline = "\r\n" if b"\r\n" in first else "\n"
        expected = first.decode("utf-8")
        substitutions = [
            ("Chem.rdchem.HybridizationType.SP3," + newline + "                      Chem.rdchem.HybridizationType.AROMATIC]",
             "Chem.rdchem.HybridizationType.SP3]"),
            ("        cids = AllChem.EmbedMultipleConfs(" + newline +
             "            mol, numConfs=self.num_confs," + newline +
             "            params=AllChem.ETKDGv3()," + newline +
             "            pruneRmsdThreshold=self.rmsd_thresh" + newline + "        )",
             "        params = AllChem.ETKDGv3()" + newline +
             "        params.pruneRmsThresh = self.rmsd_thresh" + newline +
             "        cids = AllChem.EmbedMultipleConfs(mol, self.num_confs, params)"),
            ("rdFreeSASA.calcSASA(mol, radii_table, confId=cid)", "rdFreeSASA.CalcSASA(mol, radii_table, confIdx=cid)"),
            ("ff.Minimize(maxIters=400)", "ff.Minimize(maxIts=400)"),
        ]
        for index, (before, after) in enumerate(substitutions):
            check(f"compatibility repair {index + 1} unique source anchor", expected.count(before) == 1)
            expected = expected.replace(before, after)
        check("exactly four API repairs and no scientific replacement", len(record["repairs"]) == 4 and expected.encode("utf-8") == second)
        check("audit tied to captured execution values", digest(ROOT / "results/compatibility/captured_results.json") == audit["captured_sha256"])
        check("audit script bytes", digest(ROOT / "scripts/audit_source.py") == audit["audit_script_sha256"])
        hashes(audit["output_sha256"], ROOT / "results", "source audit output")
        init = rows("results/untrained_initialization_probes.csv", 72)
        priority = rows("results/random_priority_probes.csv", 192)
        curve = rows("results/hardcoded_curve_audit.csv", 6)
        counts = audit["counts"]
        check("source audit 264 forwards and 13 model initializations",
              counts["initialization_forward_passes"] == 72 and counts["priority_forward_passes"] == 192 and
              counts["sensitivity_model_initializations"] == 12 and counts["fixed_weight_model_initializations"] == 1 and
              counts["total_model_initializations"] == 13 and counts["fixed_weight_priority_runs"] == 32)
        check("audit performs no additional kMC or conformer search", counts["new_kMC_trajectories"] == counts["new_conformer_embeddings"] == 0)
        check("12 initialization seeds each cover six molecules", Counter(r["seed"] for r in init) == {str(i): 6 for i in range(12)})
        check("32 priority seeds each rank six molecules", Counter(r["uniform_seed"] for r in priority) == {str(i): 6 for i in range(32)})
        neural = audit["source_neural_model"]
        check("untrained source and random uncertainty are not relabeled as fitted", not any(neural[k] for k in (
            "training_labels_supplied", "optimizer_or_backpropagation", "trained_checkpoint_loaded",
            "measured_reference_electrode_calibration", "GP_model_fitted", "posterior_uncertainty_computed")))
        check("non-carbon source atom argmax retained", any(r["actual_element"] != "C" for r in neural["original_selected_atom_elements"]))
        check("source curve fixed numerical entries retained", [float(r["source_hardcoded_TOF_s_1"]) for r in curve] == [4.2, 12.8, 38.4, 95.1, 184.6, 260.2])
        check("source reaction claim explicitly identified", audit["source_reaction"]["reaction_class_is_hardcoded"]
              and not audit["source_reaction"]["atom_mapping_supplied"] and not audit["source_reaction"]["bond_changes_computed"])

    def learning():
        summary = read("results/learning/summary.json")
        provenance = read("data/learning/provenance.json")
        data = rows("data/learning/selected_freesolv.csv", 256)
        check("selected data hashes and provenance", digest(ROOT / "data/learning/selected_freesolv.csv") ==
              provenance["selected_sha256"] == summary["selected_data_sha256"] and
              digest(ROOT / "data/learning/provenance.json") == summary["dataset_provenance_sha256"])
        check("official FreeSolv snapshot and license hashes", digest(ROOT / "data/learning/official_database.txt") ==
              provenance["official_database_sha256"] and digest(ROOT / "data/learning/FreeSolv_LICENSE.txt") == provenance["license_sha256"])
        check("hydration target and attribution retained", "hydration" in provenance["dataset"].lower() and
              provenance["target_unit"] == "kcal/mol" and "Mobley" in provenance["attribution"] and "CC-BY-4.0" in provenance["license"])
        official = {}
        for line in (ROOT / "data/learning/official_database.txt").read_text(encoding="utf-8").splitlines():
            if line.strip() and not line.startswith("#"):
                values = [value.strip() for value in line.split(";")]
                official[values[0]] = values
        check("selected labels and uncertainty fields match archived official data",
              all(near(r["expt_kcal_mol"], official[r["molecule_id"]][3]) and
                  near(r["experimental_uncertainty_kcal_mol"], official[r["molecule_id"]][4]) for r in data))
        check("no duplicated selected molecule IDs or canonical structures", len({r["molecule_id"] for r in data}) ==
              len({r["smiles"] for r in data}) == 256)
        splits = read("results/learning/splits.json")
        check("154 51 51 split counts", summary["split_counts"] == {"train": 154, "validation": 51, "test": 51} and
              Counter(r["split"] for r in data) == summary["split_counts"])
        for split in ("train", "validation", "test"):
            indices = [i for i, r in enumerate(data) if r["split"] == split]
            check(split + " recorded IDs and groups", splits[split]["row_indices"] == indices and
                  splits[split]["molecule_ids"] == [data[i]["molecule_id"] for i in indices] and
                  splits[split]["groups"] == sorted({data[i]["group_key"] for i in indices}))
        check("declared grouping has no split leakage", all(not set(splits[a]["groups"]) & set(splits[b]["groups"])
              for a, b in [("train", "validation"), ("train", "test"), ("validation", "test")]))
        check("acyclic split boundary disclosed", "does NOT establish acyclic scaffold novelty" in provenance["acyclic_policy"])
        training_values = [float(r["expt_kcal_mol"]) for r in data if r["split"] == "train"]
        fits = summary["training"]["fits"]
        check("four actual training runs and 240 epochs", len(fits) == 4 and summary["counts"]["neural_training_runs"] == 4 and
              summary["counts"]["training_epochs_executed"] == sum(r["epochs_executed"] for r in fits) == 240)
        check("training-only target normalization", all(near(r["target_mean"], statistics.mean(training_values)) and
              near(r["target_scale"], statistics.pstdev(training_values)) for r in fits))
        losses = rows("results/learning/training_losses.csv", 240)
        for fit in fits:
            group = [r for r in losses if int(r["seed"]) == fit["seed"]]
            selected = min(group, key=lambda r: float(r["validation_standardized_MSE"]))
            check(fit["model"] + " validation checkpoint selection", int(selected["epoch"]) == fit["best_epoch"] and
                  near(selected["validation_standardized_MSE"], fit["validation_standardized_MSE"]))
        embedding = rows("results/learning/learned_embeddings.csv", 256)
        scaler = read("results/learning/encoder_scaler.json")
        train_ids = splits["train"]["row_indices"]
        check("latent scaler uses training IDs only", scaler["fit_ids"] == train_ids and scaler["predeclared_seed"] == 20260928)
        check("latent scaler mean and scale reconstructed", all(
            near(scaler["mean"][j], statistics.mean(float(embedding[i][f"z{j}"]) for i in train_ids)) and
            near(scaler["scale"][j], statistics.pstdev(float(embedding[i][f"z{j}"]) for i in train_ids))
            for j in range(48)))
        predictions = rows("results/learning/predictions.csv", 1536)
        metrics = rows("results/learning/regression_metrics.csv", 18)
        by_id = {r["molecule_id"]: r for r in data}
        check("prediction truth and split attribution", all(r["split"] == by_id[r["molecule_id"]]["split"] and
              near(r["observed_hydration_kcal_mol"], by_id[r["molecule_id"]]["expt_kcal_mol"]) for r in predictions))
        metric_lookup = {}
        for metric in metrics:
            group = [r for r in predictions if (r["model"], r["split"]) == (metric["model"], metric["split"])]
            residual = [float(r["predicted_hydration_kcal_mol"]) - float(r["observed_hydration_kcal_mol"]) for r in group]
            check(metric["model"] + " " + metric["split"] + " metrics reconstructed",
                  len(group) == int(metric["count"]) and near(metric["RMSE_kcal_mol"], math.sqrt(statistics.mean(x*x for x in residual))) and
                  near(metric["MAE_kcal_mol"], statistics.mean(abs(x) for x in residual)))
            metric_lookup[(metric["model"], metric["split"])] = float(metric["RMSE_kcal_mol"])
        check("mixed neural versus ridge result retained",
              metric_lookup[("mpnn_seed20260930", "test")] < metric_lookup[("descriptor_ridge", "test")] <
              min(metric_lookup[(name, "test")] for name in ("mpnn_seed20260928", "mpnn_seed20260929")))
        check("shuffled-label negative control retained", metric_lookup[("shuffled_train_mpnn_seed20261001", "test")] >
              max(metric_lookup[(f"mpnn_seed{seed}", "test")] for seed in (20260928, 20260929, 20260930)))
        campaign = rows("results/learning/campaign_evaluations.csv", 384)
        budget = read("results/learning/campaign_checks.json")
        check("24 campaigns with 16 actual unique cached calls", len(budget) == 24 and
              all(r["oracle_calls"] == r["unique_calls"] == r["budget"] == 16 and r["test_pool_only"] for r in budget))
        test_ids = splits["test"]["row_indices"]
        optimum = max(-float(data[i]["expt_kcal_mol"]) for i in test_ids)
        grouped = defaultdict(list)
        for row in campaign:
            grouped[(int(row["seed"]), row["method"])].append(row)
        methods = ("latent_gp_ucb", "descriptor_gp_ucb", "random")
        check("campaign identities", set(grouped) == {(seed, method) for seed in range(8) for method in methods})
        for seed in range(8):
            starts = []
            for method in methods:
                group = grouped[(seed, method)]
                starts.append([r["candidate_index"] for r in group[:4]])
                check(f"campaign {seed} {method} budget and order", len(group) == len({r["candidate_index"] for r in group}) == 16 and
                      [int(r["step"]) for r in group] == list(range(1, 17)) and
                      [int(r["observed_count_before_selection"]) for r in group] == list(range(16)))
                best = -math.inf
                for row in group:
                    candidate = data[test_ids[int(row["candidate_index"])]]
                    value = -float(candidate["expt_kcal_mol"])
                    best = max(best, value)
                    check(f"campaign {seed} {method} step {row['step']} observed target and regret",
                          candidate["molecule_id"] == row["molecule_id"] and near(value, row["objective_negative_hydration_kcal_mol"]) and
                          near(best, row["best_observed_objective"]) and near(optimum - best, row["simple_regret_kcal_mol"]))
                    if method != "random" and int(row["step"]) > 4:
                        check(f"campaign {seed} {method} step {row['step']} recorded posterior acquisition",
                              math.isfinite(float(row["posterior_sd_selected"])) and float(row["posterior_sd_selected"]) >= 0 and
                              near(row["UCB_selected"], float(row["posterior_mean_selected"]) + 1.96 * float(row["posterior_sd_selected"])))
            check(f"campaign {seed} shared four initial queries", starts[0] == starts[1] == starts[2])
        paired = rows("results/learning/campaign_paired_comparisons.csv", 3)
        check("learned GP versus descriptor GP terminal tie retained",
              any(r["left"] == "latent_gp_ucb" and r["right"] == "descriptor_gp_ucb" and
                  int(r["ties"]) == 8 and near(r["mean_regret_reduction_left_vs_right"], 0) for r in paired))
        check("retrospective benchmark has no new experiments or DFT", summary["counts"]["new_experiments"] == summary["counts"]["new_DFT"] == 0)
        check("frozen-encoder and conditional GP uncertainty boundaries", "not jointly trained" in " ".join(summary["limitations"]) and
              "neither experimental uncertainty" in summary["active_learning"]["uncertainty"] and
              "not counted as pool calls" in summary["active_learning"]["extra_pretraining_cost"])

    def kinetics():
        summary = read("results/kinetics/summary.json")
        cases = rows("results/kinetics/trajectories.csv", 273)
        check("273 study trajectories and 9950504 actual events", summary["counts"]["study_SSA_trajectories"] == len(cases) == 273 and
              sum(int(r["total_events"]) for r in cases) == summary["counts"]["executed_SSA_events"] == 9950504)
        check("kinetics family accounting", Counter(r["family"] for r in cases) ==
              {"potential_scan": 192, "baseline": 1, "site_scaling": 48, "two_rate_classes": 16, "fixed_event_transient": 16})
        check("distinct study seeds", len({r["seed"] for r in cases}) == len(cases))
        for row in cases:
            final, burn, events = (json.loads(row[key]) for key in
                                   ("final_counts_by_class", "burn_in_counts_by_class", "reaction_counts_observed_by_class"))
            balanced = all(final[c][i] - burn[c][i] == events[c][(i-1) % 5] - events[c][i]
                           for c in range(len(final)) for i in range(5))
            charge = all(events[c][1] + events[c][3] - 2 * events[c][4] ==
                         sum((final[c][i] - burn[c][i]) * (0, 0, 1, 1, 2)[i] for i in range(5)) for c in range(len(final)))
            positive = all(isinstance(x, int) and x >= 0 for matrix in (final, burn, events) for cls in matrix for x in cls)
            check("kinetics seed " + row["seed"] + " exact class balances", balanced and charge and positive and
                  sum(sum(cls) for cls in final) == int(row["sites"]) and float(row["occupancy_sum_error"]) < 1e-12)
            duration, sites = float(row["observation_duration_s"]), int(row["sites"])
            check("kinetics seed " + row["seed"] + " TOF and current from counts",
                  near(row["TOF_s_1"], int(row["products_observed"]) / (sites * duration)) and
                  near(row["electrons_per_site_s"], int(row["electrons_observed"]) / (sites * duration)) and
                  math.isclose(float(row["current_for_explicit_sites_A"]), int(row["electrons_observed"]) * 1.602176634e-19 / duration,
                               rel_tol=1e-12, abs_tol=1e-30))
        check("257 prescribed-horizon trajectories and 16 event-cap diagnostics",
              all(float(r["simulated_time_s"]) == 1.0 and float(r["burn_in_s"]) == .2 for r in cases if r["family"] != "fixed_event_transient") and
              all(int(r["total_events"]) == 3500 and r["stop_reason"] == "event_limit" for r in cases if r["family"] == "fixed_event_transient"))
        potential = rows("results/kinetics/potential_summary.csv", 6)
        for row in potential:
            eta = float(row["eta_V"])
            factor = math.exp(.5 * 96485.33 * eta / (8.314 * 298.15))
            reference = 1 / (1/45 + 1/(120*factor) + 1/350 + 1/(200*factor) + 1/80)
            check("CTMC and closed-form rate identity eta " + row["eta_V"], near(row["steady_TOF_s_1"], reference) and
                  abs(float(row["finite_window_CTMC_TOF_s_1"]) - reference) < 1e-7 and
                  float(row["CTMC_probability_sum_error"]) < 1e-8 and float(row["CTMC_window_occupancy_sum_error"]) < 1e-8)
        check("analytic saturation below original hardcoded 260.2 retained", near(summary["rate_limited_TOF_upper_bound_s_1"],
              1/(1/45 + 1/350 + 1/80)) and summary["source_hardcoded_TOFs_not_used"][-1] == 260.2 and
              not summary["model"]["inter_site_interactions"] and not summary["model"]["spatial_periodic_boundary_operations"])
        for name, expected in (("analytic_rate_perturbations", 105), ("analytic_transient", 101),
                               ("cycle_wait_distribution", 151), ("site_scaling", 3), ("fixed_event_transients", 16)):
            rows("results/kinetics/" + name + ".csv", expected)
        trajectory = rows("results/kinetics/baseline_trajectory.csv", 33834)
        baseline = read("results/kinetics/baseline.json")
        integrated = [0.] * 5
        monotonic = True
        for left, right in zip(trajectory, trajectory[1:]):
            t0, t1 = float(left["time_s"]), float(right["time_s"])
            monotonic = monotonic and t1 > t0
            dt = max(0., t1 - max(t0, .2))
            for i, state in enumerate(STATES):
                integrated[i] += dt * int(left["n_" + state])
        check("full baseline post-event log reconstructs time occupancy", monotonic and
              all(near(integrated[i] / (256*.8), baseline["time_average_coverage"][i]) for i in range(5)) and
              [int(trajectory[-1]["n_" + state]) for state in STATES] == baseline["final_counts"])
        check("prior read-only independent kMC audit passed", read("results/kinetics/audit.json")["passed"])

    def structure():
        summary = read("results/structure/summary.json")
        c = summary["counts"]
        check("structure study count partition", [c[k] for k in ("molecules", "replica_runs", "requested_conformers",
              "returned_after_embedding_pruning", "MMFF_minimizations", "converged", "nonconverged",
              "pooled_unique_minima", "SASA_orientation_evaluations", "coordinate_files")] ==
              [7, 21, 504, 120, 120, 120, 0, 36, 2880, 14])
        conformers = rows("results/structure/conformers.csv", 120)
        weights = rows("results/structure/minima_weights.csv", 357)
        ensembles = rows("results/structure/ensemble_statistics.csv", 84)
        quadrature = rows("results/structure/surface_quadrature.csv", 2880)
        rows("results/structure/seed_variability.csv", 21)
        check("36 pooled distinct minima reproduced from ensemble rows", sum(int(r["unique_minima"]) for r in ensembles
              if r["scope"] == "pooled" and near(r["temperature_K"], 298.15)) == 36)
        check("all recorded MMFF optimizations converged", all(int(r["optimization_status"]) == 0 for r in conformers))
        conformer_index = {(r["molecule_id"], r["conformer_id"]): r for r in conformers}
        check("unique conformer identities", len(conformer_index) == 120)
        weight_groups = defaultdict(list)
        for row in weights:
            weight_groups[(row["molecule_id"], row["scope"], row["temperature_K"])].append(row)
        check("84 ensemble weight groups", len(weight_groups) == 84)
        for ensemble in ensembles:
            key = (ensemble["molecule_id"], ensemble["scope"], ensemble["temperature_K"])
            group = weight_groups[key]
            T, R = float(ensemble["temperature_K"]), float(summary["method"]["R_kcal_mol_K"])
            factors = [math.exp(-float(r["relative_energy_kcal_mol"]) / (R*T)) for r in group]
            normalization = sum(factors)
            good = all(near(r["forcefield_minima_weight"], f/normalization) for r, f in zip(group, factors))
            check("Boltzmann minima weights " + " ".join(key), good and
                  near(sum(float(r["forcefield_minima_weight"]) for r in group), 1) and
                  int(ensemble["unique_minima"]) == len(group))
            for column, descriptor in (("weighted_sasa_A2", "sasa_A2"),
                                       ("weighted_mass_weighted_Rg_A", "mass_weighted_Rg_A"),
                                       ("weighted_unweighted_Rg_A", "unweighted_Rg_A")):
                predicted = sum(float(r["forcefield_minima_weight"]) *
                                float(conformer_index[(r["molecule_id"], r["conformer_id"])][descriptor]) for r in group)
                check("ensemble " + " ".join(key) + " " + descriptor, near(predicted, ensemble[column]))
        quads = defaultdict(list)
        for row in quadrature:
            quads[(row["molecule_id"], row["conformer_id"])].append(row)
        check("24 SASA orientations per conformer", len(quads) == 120 and
              all(len(group) == 24 and {int(r["orientation_index"]) for r in group} == set(range(24)) for group in quads.values()))
        check("SASA average reconstructed from finite orientations", all(near(
              statistics.mean(float(r["sasa_A2"]) for r in group), conformer_index[key]["sasa_A2"]) for key, group in quads.items()))
        sdf_files = list((ROOT / "results/structure").glob("*.sdf"))
        check("14 coordinate files and 240 embedded plus optimized records", len(sdf_files) == 14 and
              sum(p.read_text(encoding="utf-8").count("$$$$") for p in sdf_files) == 240)
        for molecule in summary["molecule_summary"]:
            name = molecule["molecule_id"]
            record = read("results/structure/" + name + ".json")
            hashes(record["coordinate_sha256"], ROOT / "results/structure", name + " coordinate snapshot")
            check(name + " formula count and coordinate correspondence", record["formula"] == molecule["formula"] and
                  len(record["records"]) == molecule["returned_after_embedding_pruning"] and
                  all((ROOT / "results/structure" / (name + suffix)).read_text(encoding="utf-8").count("$$$$") ==
                      molecule["returned_after_embedding_pruning"] for suffix in ("_embedded.sdf", "_optimized.sdf")))
        reaction = read("results/structure/reaction_audit.json")
        source_rxn = reaction["source_unmapped_reaction"]
        check("source formula and absent mapping defect retained", source_rxn["status"] == "insufficient_atom_mapping" and
              source_rxn["inventory"]["reactants"]["formula"] == "C22H21NOS" and
              source_rxn["inventory"]["products"]["formula"] == "C21H17NOS" and
              source_rxn["inventory"]["element_delta_products_minus_reactants"] == {"C": -1, "H": -4, "N": 0, "O": 0, "S": 0} and
              source_rxn["reaction_center"] is None and source_rxn["mechanism"] is None)
        example = reaction["educational_balanced_example"]
        check("mapped educational CGR computes four bond changes without mechanism", example["inventory"]["elementally_balanced"] and
              len(example["bond_changes"]) == 4 and example["reaction_center_atom_maps"] == [2, 3, 4, 5] and
              example["mechanistic_or_reaction_class_claim"] is None)
        counterexample = reaction["equal_heavy_count_counterexample"]
        check("heavy atom equality is not elemental conservation", counterexample["heavy_atom_count_equal"] and
              not counterexample["elementally_balanced"])
        check("target graph identity retained", reaction["target_identity"]["formula"] == "C15H13NO" and
              all(reaction["target_identity"]["checks"].values()))
        rejected = read("results/structure/input_rejections.json")
        check("three unsupported or invalid inputs retained", len(rejected) == 3 and
              all(r["status"] in ("unsupported", "rejected") for r in rejected))
        check("force-field and search limits not recast as quantum or experiment", c["QM_calculations"] == c["experimental_measurements"] == 0 and
              "unit degeneracy" in summary["method"]["weights"] and "not buried volume" in " ".join(summary["limitations"]))

    def figures():
        manifest = read("results/figure_manifest.json")["figures"]
        expected = {f"{stem}_{language}.{format}" for stem in STEMS for language in ("english", "chinese") for format in ("png", "svg")}
        check("16 reviewed figures with specified names", len(manifest) == 16 and {Path(r["file"]).name for r in manifest} == expected and
              Counter(r["format"] for r in manifest) == {"png": 8, "svg": 8})
        for figure in manifest:
            path = ROOT / figure["file"]
            if not path.is_file():
                path = ROOT / "reports/figures" / figure["file"]
            name = Path(figure["file"]).name
            check("figure current bytes " + name, path.is_file() and digest(path) == figure["sha256"])
            hashes(figure["sources"], ROOT, name + " source")
            if figure["format"] == "png":
                with Image.open(path) as im:
                    check("PNG pixels and resolution " + name, list(im.size) == figure["pixels"] == [2400, 1380] and
                          all(abs(value - 300) < .1 for value in im.info.get("dpi", (0, 0))))
            else:
                tree = ET.parse(path)
                ns = {"s": "http://www.w3.org/2000/svg"}
                nodes = tree.findall(".//s:text", ns)
                labels = " ".join("".join(node.itertext()) for node in nodes)
                check("editable SVG without embedded raster " + name, len(nodes) > 10 and
                      not tree.findall(".//s:image", ns) and figure["editable_text"] and not figure["embedded_raster"])
                check("visible figure evidence boundary " + name, bool(re.search(
                    r"UNTRAINED|HYDRATION|FORCE.FIELD|MMFF94|ASSUMED|MODEL|NUMERICAL|AUDIT|未训练|水合|力场|假定|假设|模型|数值|审计", labels, re.I)))
        qa = read("results/figure_qa.json")
        actual = {Path(r["file"]).name: r["sha256"] for r in manifest}
        recorded = {Path(name).name: value for name, value in qa["reviewed_figure_hashes"].items()}
        check("visual QA covers all current eight PNG and eight SVG", actual == recorded and
              qa["png_files_visually_reviewed"] == 8 and qa["svg_files_rendered_and_visually_reviewed"] == 8)
        check("visual QA identifies current generator", qa["generator_sha256"] == digest(ROOT / "scripts/plot_reviewed.py"))

    def cross_review():
        record = read("results/cross_review.json")
        check("three peer-agent reviews and coordinator review recorded", len(record["reviews"]) == 4 and
              {r["reviewer"] for r in record["reviews"]} == {
                  "electrograph_structure", "electrograph_kinetics", "electrograph_learning", "coordinating_assistant"})
        check("cross-review records have no unresolved blocking findings",
              all(r["blocking_findings"] == [] for r in record["reviews"]))
        for group, mapping in record["reviewed_artifact_sha256"].items():
            hashes(mapping, ROOT, "cross-review current evidence " + group)

    def reports():
        texts = {}
        expected = [f"{section}.{index}" for section, count in [(1, 3), (2, 4), (3, 3), (4, 4), (5, 3)]
                    for index in range(1, count + 1)]
        for language in ("english", "chinese"):
            path = ROOT / "reports" / f"electrograph_report_{language}.md"
            text = path.read_text(encoding="utf-8")
            texts[language] = text
            check(language + " five main sections", re.findall(r"^## (\d+)\.", text, re.M) == ["1", "2", "3", "4", "5"])
            check(language + " seventeen specified ordered subsections", re.findall(r"^### (\d+\.\d+) ", text, re.M) == expected)
            for number, stem in enumerate(STEMS, 1):
                check(language + f" figure {number} caption and PNG SVG links",
                      ("**" + ("Figure " if language == "english" else "图 ") + str(number)) in text and
                      all(stem + "_" + language + "." + format in text for format in ("png", "svg")))
            check(language + " no unfinished placeholders and balanced fences",
                  not re.search(r"@@[A-Z_]+@@", text) and text.count(chr(96)*3) % 2 == 0 and
                  text.count("\\[") == text.count("\\]") and text.count("\\(") == text.count("\\)"))
        numbers = lambda text: [re.findall(r"(?<![A-Za-z0-9.])[-−]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", line)
                                for line in text.splitlines() if line.startswith("|")]
        check("separate English report contains no Chinese prose", not re.search(r"[\u4e00-\u9fff]", texts["english"]))
        check("English Chinese numerical table parity", numbers(texts["english"]) == numbers(texts["chinese"]))

    def local_links():
        check("module README present", (ROOT / "README.md").is_file())
        paths = list(ROOT.glob("*.md")) + list((ROOT / "reports").rglob("*.md"))
        check("authored documentation available", bool(paths))
        for path in sorted(paths):
            text = path.read_text(encoding="utf-8")
            for link in re.findall(r"\]\(([^)]+)\)", text):
                if link.startswith(("https://", "http://", "mailto:", "#", "data:")):
                    continue
                name = unquote(link.split("#", 1)[0].strip("<>"))
                target = path.parent / name
                check(path.relative_to(ROOT).as_posix() + " local link " + link,
                      target.exists() or target.resolve() == RESULT.resolve())

    for name, function in [("module provenance", module_provenance), ("source", source), ("learning", learning),
                           ("kinetics", kinetics), ("structure", structure), ("figures", figures),
                           ("cross-review", cross_review), ("reports", reports), ("local links", local_links)]:
        guarded(name, function)
    result = {
        "passed": all(item["passed"] for item in checks),
        "checks": len(checks),
        "failed": [item for item in checks if not item["passed"]],
        "scope": "Saved source provenance, numerical accounting, negative results, bilingual report parity, local links and publication QA records. No model retraining, stochastic reruns, new chemistry or physical validation.",
        "validator_sha256": digest(Path(__file__)),
        "items": checks,
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({key: result[key] for key in ("passed", "checks", "failed")}, indent=2, ensure_ascii=False))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
