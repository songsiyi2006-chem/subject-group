"""Validate frozen manuscript evidence without executing any scientific solver.

Only the Python standard library is required.  Run from any working directory:
    python manuscript/scripts/validate_manuscript.py
    python manuscript/scripts/validate_manuscript.py --require-documents

The latter also requires all eight publication documents, rendered Markdown,
and the root author's document QA record.  File checks are not visual review,
literature re-verification, or experimental validation.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import struct
import unicodedata
import xml.etree.ElementTree as ET
import zipfile
import zlib

ROOT = Path(__file__).resolve().parents[2]
HERE = ROOT / "manuscript"
OUTPUT = HERE / "results/validation.json"
FIGURES = ("Fig1_Evidence_Design", "Fig3_Gain_Transfer",
           "Fig4_MSE_Attribution", "Fig6_Observable_Convergence")
LANGUAGES = ("english", "chinese")


class Audit:
    def __init__(self, require_documents=False):
        self.require_documents = require_documents
        self.checks = []
        self.inputs = {}
        self.maxima = {}
        self.details = {}

    def check(self, name, condition, detail=None):
        item = {"id": name, "status": "pass" if condition else "fail"}
        if detail is not None:
            item["detail"] = detail
        self.checks.append(item)
        return bool(condition)

    def pending(self, name, detail):
        self.checks.append({"id": name, "status": "fail" if self.require_documents else "pending",
                            "detail": detail})

    def read(self, relative):
        path = ROOT / relative
        raw = path.read_bytes()
        self.inputs[path.relative_to(ROOT).as_posix()] = hashlib.sha256(raw).hexdigest()
        return raw

    def text(self, relative):
        return self.read(relative).decode("utf-8-sig")

    def js(self, relative):
        return json.loads(self.text(relative))

    def csv(self, relative):
        return list(csv.DictReader(self.text(relative).splitlines()))

    def hash(self, name, relative, expected):
        raw = self.read(relative)
        return self.check(name, hashlib.sha256(raw).hexdigest() == expected,
                          {"path": str(relative), "expected_sha256": expected})

    def close(self, name, actual, expected, *, atol=1e-14, rtol=1e-12, category=None):
        a, e = float(actual), float(expected)
        difference = abs(a - e)
        if category:
            self.maxima[category] = max(self.maxima.get(category, 0.0), difference)
        return self.check(name, math.isfinite(a) and math.isfinite(e)
                          and difference <= atol + rtol * abs(e),
                          {"actual": a, "expected": e, "absolute_difference": difference,
                           "absolute_tolerance": atol, "relative_tolerance": rtol})

    def section(self, name, function):
        try:
            function(self)
        except Exception as exc:
            self.check(name + ".unhandled_input_error", False,
                       {"type": type(exc).__name__, "message": str(exc)})


def pointer(value, key):
    for token in key.strip("/").split("/") if key else []:
        token = token.replace("~1", "/").replace("~0", "~")
        value = value[int(token)] if isinstance(value, list) else value[token]
    return value


def check_catalog(a):
    data = a.js("manuscript/evidence_catalog.json")
    a.check("catalog.module_count", len(data["modules"]) == 12)
    a.check("catalog.no_unresolved_selectors", not data["catalog_internal_verification"]["unresolved_sources"])
    a.check("catalog.author", data["author"] == {"name": "宋思毅", "affiliation": "广西师范大学"})
    files = data["artifact_inventory"]
    a.check("catalog.unique_paths", len({x["path"] for x in files}) == len(files))
    for item in files:
        a.hash("catalog.hash." + item["path"], item["path"], item["sha256"])
    count = 0
    for module in data["modules"]:
        for entry in (module["workload_records"] + module["quantitative_evidence"]
                      + module["separate_pilot_records"]):
            src = entry["source"]
            key = "catalog.selector." + module["id"] + "." + str(count)
            if src["type"] == "json_pointer":
                actual = pointer(a.js(src["path"]), src["pointer"])
                valid = actual == entry["value"]
            elif src["type"] == "csv_row_count":
                valid = len(a.csv(src["path"])) == entry["value"]
            elif src["type"] == "csv_filter":
                actual = [{"csv_line": i + 2, "fields": row}
                          for i, row in enumerate(a.csv(src["path"]))
                          if all(row.get(k) == str(v) for k, v in src["where"].items())]
                valid = actual == entry["records"] and len(actual) == entry["selected_rows"]
            else:
                valid = False
            a.check(key, valid, src)
            count += 1
        for path in module["report_paths"] + module["raw_and_audit_paths"]:
            a.check("catalog.link." + path, (ROOT / path).is_file())
    a.check("catalog.selector_count", count == data["catalog_internal_verification"]["selected_records_checked"])
    accounting = data["homogeneous_job_accounting"]
    xtb = [p for group in accounting["xTB"]["parts"] for p in group["records"]]
    a.check("accounting.xTB55", len(xtb) == len(set(xtb)) == 55)
    a.check("accounting.xTB_success", all(a.js(p)["returncode"] == 0 for p in xtb))
    old = sum(a.js("quantumequi/results/electronic/" + sub + "reference_summary.json")["counts"]["total_SCF_jobs"]
              for sub in ("", "pilot/"))
    correlation = a.js("quantumequi/results/extensions/correlation/summary.json")["complete_study_accounting"]
    a.check("accounting.Psi4_energy_drivers", old + correlation["energy_driver_calls_total"] == 261)
    a.check("accounting.FCI_drivers", correlation["explicit_FCI_drivers_total"] == 87)
    a.check("accounting.no_heterogeneous_grand_total", accounting["grand_total"] is None)
    a.details["catalog"] = {"modules": 12, "selectors": count, "artifact_hashes": len(files)}


def check_analysis(a):
    prefix = "manuscript/results/"
    summary = a.js(prefix + "analysis_summary.json")
    for path, digest in summary["inputs_sha256"].items():
        a.hash("analysis.input." + path, path, digest)
    for name, digest in summary["outputs_sha256"].items():
        a.hash("analysis.output." + name, prefix + name, digest)
    a.hash("analysis.script", "manuscript/scripts/analyze_error_chain.py", summary["script_sha256"])
    source = a.csv("quantumequi/results/extensions/error_budget/matched_error_components.csv")
    attribution = a.csv(prefix + "mse_attribution.csv")
    gains = a.csv(prefix + "paired_gain_transfer.csv")
    richardson = a.csv(prefix + "tail_richardson.csv")
    comparison = a.csv(prefix + "residual_observable_comparison.csv")
    tail = a.csv("quantumequi/results/extensions/vibration/tail_followup/diagnostics.csv")
    a.check("analysis.source_rows", len(source) == 126)
    a.check("analysis.unique_geometry_count", len({float(x["R_A"]) for x in source}) == 21)
    a.check("analysis.unique_model_geometry_rows", len({(x["seed"], x["objective"], x["R_A"]) for x in source}) == 126)
    a.check("analysis.zero_replay_drift", all(float(x["reference_replay_drift_Hartree"]) == 0 for x in source))
    a.check("analysis.output_group_counts", [len(attribution), len(gains), len(richardson), len(comparison)] == [12, 6, 3, 7])
    truth = {}
    expected_keys = {(s, str(seed), obj) for s in ("test", "ood") for seed in (7301, 7302, 7303)
                     for obj in ("energy_only", "energy_gradient")}
    a.check("analysis.group_keys", {(x["split"], x["seed"], x["objective"]) for x in attribution} == expected_keys)
    a.check("analysis.paired_group_keys", {(x["split"], x["seed"]) for x in gains}
            == {(split, str(seed)) for split in ("test", "ood") for seed in (7301, 7302, 7303)})
    geometry_sets = {}
    negative_cross = 0
    for row in attribution:
        key = (row["split"], row["seed"], row["objective"])
        selected = [x for x in source if x["seed"] == key[1] and x["objective"] == key[2]
                    and (x["split"] == "test" if key[0] == "test" else "ood" in x["split"].lower())]
        distances = tuple(sorted(float(x["R_A"]) for x in selected))
        geometry_sets.setdefault(key[0], set()).add(distances)
        a.check("analysis.n." + str(key), len(selected) == (8 if key[0] == "test" else 5))
        learning = [float(x["learning_error_Hartree"]) for x in selected]
        bias = [float(x["RHF_minus_FCI_Hartree"]) for x in selected]
        total = [float(x["total_error_vs_FCI_Hartree"]) for x in selected]
        # Independent compensated summation directly over source rows; do not
        # import the manuscript analysis or infer total MSE from its components.
        n = len(selected)
        ml = math.fsum(x*x for x in learning) / n
        mb = math.fsum(x*x for x in bias) / n
        mt = math.fsum(x*x for x in total) / n
        cross = 2 * math.fsum(x*y for x, y in zip(learning, bias)) / n
        rl, rb, rt = math.sqrt(ml), math.sqrt(mb), math.sqrt(mt)
        values = {"learning_MSE_Hartree2": ml, "reference_bias_MSE_Hartree2": mb,
                  "total_MSE_Hartree2": mt, "cross_term_Hartree2": cross,
                  "learning_RMSE_Hartree": rl, "reference_bias_RMSE_Hartree": rb,
                  "total_FCI_RMSE_Hartree": rt,
                  "MSE_identity_residual_Hartree2": mt-ml-mb-cross,
                  "triangle_lower_bound_Hartree": abs(rb-rl),
                  "triangle_upper_bound_Hartree": rb+rl,
                  "n_geometries": n,
                  "signed_cancellation_points": sum(x*y < 0 for x, y in zip(learning, bias))}
        for field, value in values.items():
            a.close("analysis.recalculate." + str(key) + "." + field, row[field], value,
                    atol=2e-16, rtol=2e-13, category=field)
        a.close("analysis.signed_MSE_identity." + str(key), mt, ml+mb+cross, atol=2e-16, rtol=2e-13)
        a.check("analysis.triangle." + str(key), abs(rb-rl)-1e-15 <= rt <= rb+rl+1e-15)
        negative_cross += cross < 0
        truth[key] = values
    a.check("analysis.same_geometries_across_models", all(len(x) == 1 for x in geometry_sets.values()))
    a.check("analysis.OOD_distances", next(iter(geometry_sets["ood"])) == (1.9, 2.1, 2.3, 2.5, 2.7))
    a.check("analysis.negative_cross_retained", negative_cross == 5)
    for row in gains:
        key = (row["split"], row["seed"])
        initial, joint = truth[key + ("energy_only",)], truth[key + ("energy_gradient",)]
        for column, metric in [("learning_RMSE_reduction_percent", "learning_RMSE_Hartree"),
                               ("total_FCI_RMSE_reduction_percent", "total_FCI_RMSE_Hartree")]:
            a.close("analysis.paired." + str(key) + "." + column, row[column],
                    100*(1-joint[metric]/initial[metric]), atol=2e-11, rtol=2e-13, category=column)
        a.close("analysis.paired_initial." + str(key), row["energy_only_FCI_RMSE_Hartree"], initial["total_FCI_RMSE_Hartree"])
        a.close("analysis.paired_joint." + str(key), row["joint_FCI_RMSE_Hartree"], joint["total_FCI_RMSE_Hartree"])
    grid = {(float(x["right_A"]), int(x["intervals"])): x for x in tail}
    a.check("analysis.unique_tail_grids", len(grid) == len(tail) == 7)
    a.check("analysis.residual_comparison_grid_keys",
            {(float(x["right_A"]), int(x["intervals"])) for x in comparison} == set(grid))
    expected_pairs = {(48.0, 7200, 14400), (48.0, 14400, 28800), (96.0, 57600, 115200)}
    a.check("analysis.Richardson_pairing", {(float(x["right_A"]), int(x["coarse_intervals"]), int(x["fine_intervals"])) for x in richardson} == expected_pairs)
    for row in richardson:
        boundary = float(row["right_A"])
        coarse = grid[(boundary, int(row["coarse_intervals"]))]
        fine = grid[(boundary, int(row["fine_intervals"]))]
        a.check("analysis.Richardson_same_state." + str(boundary) + row["fine_intervals"],
                int(coarse["bound_states_returned"]) == int(fine["bound_states_returned"]) == 17)
        a.close("analysis.Richardson_grid_ratio." + row["fine_intervals"],
                float(coarse["grid_spacing_A"])/float(fine["grid_spacing_A"]), 2)
        bc, bf = float(coarse["last_state_binding_Hartree"]), float(fine["last_state_binding_Hartree"])
        exact = float(fine["last_full_line_exact_binding_Hartree"])
        extrapolated = bf + (bf-bc)/3
        a.close("analysis.Richardson_value." + row["fine_intervals"], row["extrapolated_binding_Hartree"], extrapolated,
                atol=1e-20, rtol=1e-13, category="Richardson_binding_Hartree")
        a.close("analysis.Richardson_error." + row["fine_intervals"], row["extrapolated_signed_relative_error_percent"],
                100*(extrapolated/exact-1), atol=2e-11, rtol=1e-13, category="Richardson_relative_percentage_points")
        a.close("analysis.Richardson_fine_error." + row["fine_intervals"], row["fine_signed_relative_error_percent"],
                100*(bf/exact-1), atol=2e-11, rtol=1e-13)
    for row in comparison:
        key = (float(row["right_A"]), int(row["intervals"]))
        stored = grid[key]
        residual = float(stored["eigen_residual_max_Hartree"])
        a.close("analysis.residual." + str(key), row["eigen_residual_Hartree"], residual, atol=0, rtol=0)
        a.close("analysis.partition." + str(key), row["F4000K_error_Hartree"], stored["F4000K_numerical_minus_full_line_Morse_Hartree"], atol=0, rtol=0)
        if int(stored["bound_states_returned"]) == 17:
            error = abs(float(stored["last_state_binding_Hartree"]) - float(stored["last_full_line_exact_binding_Hartree"]))
            a.close("analysis.observable_error." + str(key), row["same_state_binding_absolute_error_Hartree"], error, atol=0, rtol=0)
            a.close("analysis.error_residual_ratio." + str(key), row["error_to_matrix_residual_ratio"], error/residual, atol=0, rtol=1e-14, category="residual_ratio")
        else:
            a.check("analysis.missing_state_not_misidentified." + str(key), row["same_state_binding_absolute_error_Hartree"] == row["error_to_matrix_residual_ratio"] == "")
    for name in ("new_quantum_jobs", "new_training_runs", "new_eigensolves"):
        a.check("analysis.no_new_computation." + name, summary["counts"][name] == 0)
    a.details["independent_recalculation"] = {
        "method": "Compensated stdlib sums from frozen point-level CSV, independent of analyze_error_chain.py functions.",
        "groups": 12, "paired_gains": 6, "Richardson_estimates": 3,
        "residual_observable_rows": 7, "matched_test_points": 8, "matched_OOD_points": 5,
        "negative_cross_terms": negative_cross,
        "interpretation": "Cross terms are signed cancellation/reinforcement, not positive variance percentages. Richardson uses saved fixed-domain grids, not new eigensolves."}


def check_references(a):
    data = a.js("manuscript/references_verified.json")
    refs = data["references"]
    keys = {x["citation_key"] for x in refs}
    a.check("references.count25", len(refs) == len(keys) == data["reference_count"] == 25)
    a.check("references.unique_ids", len({x["id"] for x in refs}) == 25)
    doi = [x["doi"].lower() for x in refs if x.get("doi")]
    a.check("references.unique_nonnull_DOI", len(doi) == len(set(doi)))
    for ref in refs:
        name = "references." + ref["citation_key"]
        a.check(name + ".published", ref["publication_status"] == "published")
        a.check(name + ".source_metadata", bool(ref["verified_source"]) and all(x.get("url") and x.get("access_level") and x.get("verified_fields") for x in ref["verified_source"]))
        a.check(name + ".claim_boundary", bool(ref.get("does_not_support")) and all(ref["supported_claim"].get(lang) for lang in LANGUAGES))
        if not ref["author_list_complete"]:
            a.check(name + ".incomplete_authors_disclosed", bool(ref.get("author_note")))
    tang = [x for x in refs if x.get("identity_verification", {}).get("author") == "Hai-Tao Tang"]
    a.check("references.Tang_count4", len(tang) == data["tang_gxnu_author_verified_count"] == 4)
    for ref in tang:
        identity = ref["identity_verification"]
        a.check("references.Tang_identity." + ref["citation_key"],
                "Hai-Tao Tang" in ref["authors"] and "Guangxi Normal University" in identity["affiliation"]
                and identity.get("no_inferred_connection_to_manuscript_author") is True)
    citations = {}
    for lang in LANGUAGES:
        text = a.text(f"manuscript/manuscript_{lang}.md")
        citations[lang] = Counter(re.findall(r"@([A-Za-z0-9_-]+)", text))
        a.check("references.resolved_body." + lang, set(citations[lang]) <= keys,
                {"used_keys": sorted(citations[lang]), "unknown": sorted(set(citations[lang])-keys)})
    a.check("references.bilingual_citation_multiplicity", citations["english"] == citations["chinese"])
    bib = HERE / "references.bib"
    if bib.exists():
        bkeys = re.findall(r"@[A-Za-z]+\s*\{\s*([^,\s]+)", a.text(bib.relative_to(ROOT)))
        a.check("references.BibTeX25", len(bkeys) == 25 and set(bkeys) == keys)
    else:
        a.pending("references.BibTeX", "references.bib has not yet been produced; registry itself is checked.")
    bibliography = HERE / "bibliography.md"
    if bibliography.exists():
        text = a.text(bibliography.relative_to(ROOT))
        entries = re.findall(r"(?m)^\[(\d+)\] (.+)$", text)
        a.check("references.bibliography25", [int(x[0]) for x in entries] == list(range(1, 26)))
        for ref in refs:
            a.check("references.bibliography_title." + ref["citation_key"], ref["title"] in text)
    else:
        a.pending("references.bibliography", "Complete bibliography not yet generated.")
    citation_path = HERE / "results/citation_audit.json"
    if citation_path.exists():
        record = a.js(citation_path.relative_to(ROOT))
        a.check("references.citation_renderer_passed", record["passed"] and not record["failures"])
        a.hash("references.registry_render_hash", "manuscript/references_verified.json", record["registry_sha256"])
        a.hash("references.renderer_hash", "manuscript/scripts/finalize_references.py", record["script_sha256"])
        for name, digest in record["output_sha256"].items():
            a.hash("references.rendered_output." + name, name, digest)
        a.check("references.global_numbering25", {x["key"] for x in record["global_numbering"]} == keys
                and [x["number"] for x in record["global_numbering"]] == list(range(1, 26)))
    else:
        a.pending("references.citation_audit", "Citation renderer's record not yet available.")
    a.details["reference_scope"] = "Checks recorded primary-source metadata and local citation consistency. It does not re-browse sources, certify journal quartiles, or require padding the body with all 25 references."


def png_details(raw):
    if raw[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError("Invalid PNG signature")
    pos = 8
    dimensions = None
    dpi = None
    end = False
    while pos + 12 <= len(raw):
        size = struct.unpack(">I", raw[pos:pos+4])[0]
        kind = raw[pos+4:pos+8]
        payload = raw[pos+8:pos+8+size]
        crc = struct.unpack(">I", raw[pos+8+size:pos+12+size])[0]
        if zlib.crc32(kind+payload) & 0xffffffff != crc:
            raise ValueError("PNG CRC mismatch")
        if kind == b"IHDR":
            dimensions = struct.unpack(">II", payload[:8])
        if kind == b"pHYs" and len(payload) == 9 and payload[-1] == 1:
            dpi = [x*0.0254 for x in struct.unpack(">II", payload[:8])]
        pos += 12 + size
        if kind == b"IEND":
            end = True
            break
    if not end or pos != len(raw):
        raise ValueError("Incomplete PNG or unexpected trailing bytes")
    return {"dimensions": dimensions, "dpi": dpi}


def check_figures(a):
    data = a.js("manuscript/results/figure_manifest.json")
    files = data["files"]
    expected = {f"figures/{stem}_{lang}.{ext}" for stem in FIGURES for lang in LANGUAGES for ext in ("png", "svg")}
    a.check("figures.four_bilingual_groups", data["figure_groups"] == 4 and len(files) == 16
            and {x["file"] for x in files} == expected)
    a.hash("figures.generator_hash", "manuscript/scripts/plot_manuscript.py", data["generator_sha256"])
    image_info = {}
    for file in files:
        path = "manuscript/" + file["file"]
        a.hash("figures.hash." + file["file"], path, file["sha256"])
        a.check("figures.provenance_present." + file["file"], bool(file["sources"]))
        for source, digest in file["sources"].items():
            a.hash("figures.data." + file["file"] + "." + source, source, digest)
        raw = a.read(path)
        if path.endswith(".png"):
            info = png_details(raw)
            image_info[file["file"]] = info
            a.check("figures.PNG_integrity." + file["file"], info["dimensions"] == (2880, 1472))
            a.check("figures.PNG_dpi." + file["file"], bool(info["dpi"]) and all(abs(x-320) < .1 for x in info["dpi"]))
        else:
            tree = ET.fromstring(raw)
            text_elements = [x for x in tree.iter() if x.tag.rsplit("}", 1)[-1] == "text"]
            a.check("figures.editable_SVG." + file["file"], tree.tag.endswith("svg") and bool(text_elements))
    a.details["figures"] = {"PNG_metadata": image_info,
                           "scope": "Integrity, editable text and data provenance; visual layout requires separate inspection."}
    qa_path = HERE / "results/figure_visual_qa.json"
    if qa_path.exists():
        record = a.js(qa_path.relative_to(ROOT))
        a.hash("figures.review_manifest_hash", "manuscript/results/figure_manifest.json", record["figure_manifest_sha256"])
        a.hash("figures.review_generator_hash", "manuscript/scripts/plot_manuscript.py", record["generator_sha256"])
        a.check("figures.review_scope", record["counts"]["png_individually_viewed"] == 8
                and record["counts"]["svg_structurally_parsed"] == 8)
        for item in record["png_visual_review"]:
            a.hash("figures.review_image." + item["file"], "manuscript/" + item["file"], item["sha256"])
        a.details["figures"]["separate_visual_QA"] = qa_path.relative_to(ROOT).as_posix()


def tables(text):
    groups, current = [], []
    for line in text.splitlines() + [""]:
        if line.strip().startswith("|"):
            current.append([x.strip() for x in line.strip().strip("|").split("|")])
        elif current:
            groups.append(current)
            current = []
    return groups


def number(value):
    return float(unicodedata.normalize("NFKC", value).replace("−", "-"))


def rounded_check(a, name, shown, actual):
    shown = unicodedata.normalize("NFKC", shown).replace("−", "-")
    digits = len(shown.split(".", 1)[1]) if "." in shown else 0
    a.close(name, number(shown), actual, atol=.5001*10**(-digits), rtol=0)


def check_bilingual(a):
    texts = {lang: a.text(f"manuscript/manuscript_{lang}.md") for lang in LANGUAGES}
    a.check("language.English_no_CJK", re.search(r"[\u3400-\u9fff]", texts["english"]) is None)
    a.check("language.Chinese_has_CJK", len(re.findall(r"[\u3400-\u9fff]", texts["chinese"])) > 3000)
    tab = {lang: tables(t) for lang, t in texts.items()}
    a.check("tables.count3", len(tab["english"]) == len(tab["chinese"]) == 3)
    for index in range(3):
        en, cn = tab["english"][index], tab["chinese"][index]
        start = 1 if index < 2 else 0
        a.check("tables.numeric_parity." + str(index+1),
                [row[start:] for row in en[2:]] == [row[start:] for row in cn[2:]])
    eq = a.js("quantumequi/results/extensions/correlation/equilibrium_reference.json")["references"]
    for row, field in zip(tab["english"][0][2:], ("R_e_A", "E_min_Hartree", "D_e_Hartree", "curvature_Hartree_A2")):
        for index, ref in enumerate(eq):
            rounded_check(a, "tables.equilibrium." + field + ref["basis"], row[index+1], ref[field])
    gain = {(x["split"], x["seed"]): x for x in a.csv("manuscript/results/paired_gain_transfer.csv")}
    for row in tab["english"][1][2:]:
        record = gain[("test" if row[0] == "Test" else "ood", row[1])]
        for col, key in [(2, "learning_RMSE_reduction_percent"), (3, "total_FCI_RMSE_reduction_percent")]:
            rounded_check(a, "tables.gain." + row[0] + row[1] + key, row[col], record[key])
    rich = {(int(float(x["right_A"])), int(x["fine_intervals"])): x for x in a.csv("manuscript/results/tail_richardson.csv")}
    for row in tab["english"][2][2:]:
        record = rich[(int(row[0]), int(row[1]))]
        rounded_check(a, "tables.Richardson.fine." + row[1], row[2], record["fine_signed_relative_error_percent"])
        rounded_check(a, "tables.Richardson.extrapolated." + row[1], row[3], record["extrapolated_signed_relative_error_percent"])
    decimals = {}
    for lang, text in texts.items():
        body = re.split(r"(?m)^## (?:References|参考文献)\s*$", text)[0]
        body = re.sub(r"!?\[[^\]]*\]\([^)]*\)", "", body)
        body = re.sub(r"\[@[^\]]+\]", "", body)
        body = unicodedata.normalize("NFKC", body).replace("−", "-")
        decimals[lang] = Counter(re.findall(r"(?<![A-Za-z\d])[-+]?\d+\.\d+(?:[eE][-+]?\d+)?", body))
        images = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", text)
        a.check("figures.body_image_count." + lang, len(images) == 6)
        a.check("figures.body_image_language." + lang, all(path.endswith("_"+lang+".png") for path in images))
        for link in re.findall(r"!?\[[^\]]*\]\(([^)]+)\)", text):
            if not re.match(r"[a-z]+://|#", link):
                a.check("manuscript.local_link." + lang + "." + link, (HERE/link.split("#")[0]).is_file())
    a.check("language.body_decimal_multiset", decimals["english"] == decimals["chinese"],
            {"english_count": sum(decimals["english"].values()), "chinese_count": sum(decimals["chinese"].values()),
             "English_only": dict(decimals["english"]-decimals["chinese"]),
             "Chinese_only": dict(decimals["chinese"]-decimals["english"])})
    a.details["bilingual_scope"] = "Main tables are checked against frozen data, citations by multiplicity, decimals by multiset. This is not a proof of full semantic translation equivalence."


def check_documents(a):
    documents = []
    for kind in ("manuscript", "supporting_information"):
        for lang in LANGUAGES:
            for ext in ("docx", "pdf"):
                path = HERE / "documents" / f"{kind}_{lang}.{ext}"
                if not path.exists():
                    a.pending("documents.missing." + path.name, "Publication document not yet generated.")
                    continue
                raw = a.read(path.relative_to(ROOT))
                documents.append({"file": path.relative_to(HERE).as_posix(),
                                  "sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)})
                if ext == "docx":
                    with zipfile.ZipFile(path) as archive:
                        a.check("documents.DOCX_zip." + path.name, archive.testzip() is None)
                        tree = ET.fromstring(archive.read("word/document.xml"))
                        text = "".join(x.text or "" for x in tree.iter() if x.tag.endswith("}t"))
                        a.check("documents.DOCX_content." + path.name, len(text) > 5000)
                        a.check("documents.DOCX_placeholders." + path.name,
                                "VERIFIED_REFERENCES" not in text and not re.search(r"\[@[A-Za-z]", text))
                else:
                    a.check("documents.PDF_envelope." + path.name, raw.startswith(b"%PDF-") and b"%%EOF" in raw[-2048:])
    for kind in ("manuscript", "supporting_information"):
        for lang in LANGUAGES:
            path = HERE / "rendered" / f"{kind}_{lang}.md"
            if path.exists():
                text = a.text(path.relative_to(ROOT))
                a.check("documents.rendered_placeholder." + path.name,
                        "VERIFIED_REFERENCES" not in text and not re.search(r"\[@[A-Za-z]", text))
            else:
                a.pending("documents.rendered_missing." + path.name, "Rendered Markdown not yet generated.")
    candidates = [HERE / "results/document_qa.json", HERE / "document_qa.json"]
    qa_path = next((p for p in candidates if p.exists()), None)
    if qa_path is None:
        a.pending("documents.QA_record", "Root author's document_qa.json is not yet available; no visual/page review is inferred.")
    else:
        qa = a.js(qa_path.relative_to(ROOT))
        a.check("documents.QA_object", isinstance(qa, dict) and bool(qa))
        a.check("documents.QA_passed", qa.get("passed") is True)
        expected_names = {kind+"_"+lang for kind in ("manuscript", "supporting_information") for lang in LANGUAGES}
        a.check("documents.QA_all_four", len(qa["documents"]) == 4
                and {x["name"] for x in qa["documents"]} == expected_names)
        for document in qa["documents"]:
            name = document["name"]
            for extension in ("pdf", "docx"):
                a.hash("documents.QA_hash." + name + "." + extension,
                       document[extension], document[extension+"_sha256"])
            a.check("documents.QA_page_records." + name,
                    len(document["pages"]) == document["page_count"] > 0
                    and [p["page"] for p in document["pages"]] == list(range(1, document["page_count"]+1)))
            a.check("documents.QA_page_metadata." + name,
                    all(p["text_characters"] > 0 and len(p["image_sha256"]) == 64
                        and all(n > 0 for n in p["rendered_pixels"]) for p in document["pages"]))
            a.check("documents.QA_declared_checks." + name, all(document["checks"].values()))
        a.details["document_QA_record"] = {"path": qa_path.relative_to(ROOT).as_posix(),
                                          "sha256": a.inputs[qa_path.relative_to(ROOT).as_posix()],
                                          "top_level_keys": list(qa),
                                          "scope": "Document hashes and page-record consistency checked. No fresh PDF text extraction or visual inspection occurs in this validator; see the separate review's scope."}
    a.details["documents"] = documents


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-documents", action="store_true")
    args = parser.parse_args()
    audit = Audit(args.require_documents)
    for name, fn in [("catalog", check_catalog), ("analysis", check_analysis),
                     ("references", check_references), ("figures", check_figures),
                     ("bilingual", check_bilingual), ("documents", check_documents)]:
        audit.section(name, fn)
    counts = Counter(x["status"] for x in audit.checks)
    status = "fail" if counts["fail"] else ("pass_with_pending_documents" if counts["pending"] else "pass")
    result = {"schema_version": "manuscript.validation.v1",
              "timestamp_policy": "No execution timestamp: identical source files and CLI flags produce identical report bytes.",
              "status": status, "require_documents": args.require_documents,
              "counts": {k: counts[k] for k in ("pass", "fail", "pending")},
              "scope": "Saved-data scientific identities, source/citation integrity, bilingual tables/numbers and artifact checks. No fresh scientific calculation or external source lookup.",
              "new_quantum_jobs": 0, "new_training_runs": 0, "new_eigenproblems": 0,
              "maximum_absolute_recalculation_differences": audit.maxima,
              "details": audit.details, "failures": [x for x in audit.checks if x["status"] == "fail"],
              "pending": [x for x in audit.checks if x["status"] == "pending"],
              "checks": audit.checks, "inputs_sha256": dict(sorted(audit.inputs.items())),
              "input_hash_base": "repository root",
              "validator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "limitations": ["Automated checks do not establish experimental chemical accuracy or broad generalization.",
                              "Recorded literature metadata is checked locally, not independently re-browsed.",
                              "PNG/SVG integrity and provenance do not replace visual figure review.",
                              "The bilingual numeric audit does not certify every translated sentence."]}
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"status": status, "counts": result["counts"], "failures": result["failures"],
                      "pending": len(result["pending"])}, ensure_ascii=True, indent=2))
    raise SystemExit(1 if counts["fail"] else 0)


if __name__ == "__main__":
    main()
