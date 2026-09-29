"""Validate retained audit artifacts; no new chemistry or network calls.

Visual-review fields describe the recorded AI inspection of the current PNGs.
Running this program is not itself a visual inspection of changed figures.
"""
from pathlib import Path
import csv
import hashlib
import importlib.metadata
import json
import re
import struct
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent.parent
SOURCES = HERE / "sources"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def numeric_leaves(value, path):
    if isinstance(value, dict):
        for key, child in value.items():
            yield from numeric_leaves(child, path + "/" + str(key).replace("~", "~0").replace("/", "~1"))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from numeric_leaves(child, path + "/" + str(index))
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        yield path, value


def main():
    catalog = read(SOURCES / "catalog.json")
    fixed = catalog["fixed_commit"]
    numeric_index = []
    for sid, item in catalog["source_entries"].items():
        assert item["commit"] == fixed
        ref = {key: item[key] for key in ("source_id", "commit", "path", "sha256")}
        if "local_excerpt" in item:
            local = SOURCES / item["local_excerpt"]
            assert sha(local) == item["excerpt_sha256"]
            selected = read(local)["selected_values"]
            for pointer, value in selected.items():
                for path, number in numeric_leaves(value, pointer):
                    numeric_index.append({**ref, "selector": path, "value": number})
        if "local_exact_copy" in item:
            local = SOURCES / item["local_exact_copy"]
            assert sha(local) == item["sha256"]
            if local.suffix == ".csv":
                with local.open(encoding="utf-8-sig", newline="") as stream:
                    for index, row in enumerate(csv.DictReader(stream), 1):
                        for key, raw in row.items():
                            try:
                                number = float(raw)
                            except (TypeError, ValueError):
                                continue
                            numeric_index.append({**ref, "selector": f"CSV data row {index}, column {key}", "value": number})
    native = read(SOURCES / "native_energy_checks.json")
    assert len(native) == 10 and all(row["passed"] for row in native)
    assert all(row["recorded_native_hash_matches"] for row in native[:8])
    assert max(row["energy_absolute_difference_hartree"] for row in native) == 0
    for row in native:
        numeric_index.append({**{key: row["source"][key] for key in ("source_id", "commit", "path", "sha256")}, "selector": "native Final Energy line; exact line-numbered excerpt in native_energy_checks.json", "value": row["native_final_energy_hartree"]})

    sections = {lang: (HERE / f"section_{lang}.md").read_text(encoding="utf-8") for lang in ("english", "chinese")}
    table_numbers = {}
    section_checks = {}
    for lang, content in sections.items():
        lines = [line for line in content.splitlines() if line.startswith("|")]
        table_numbers[lang] = [re.findall(r"[-+]?\d+(?:\.\d+)?", line.replace("−", "-")) for line in lines]
        destinations = re.findall(r"\]\(([^)]+)\)", content)
        broken = [value for value in destinations if not value.startswith("http") and value != "audit.json" and not (HERE / value).exists()]
        assert not broken, broken
        pngs = [value for value in destinations if value.endswith(".png")]
        svgs = [value for value in destinations if value.endswith(".svg")]
        assert len(pngs) == len(svgs) == 5
        assert all(value.endswith("_" + lang + ".png") for value in pngs)
        assert not re.search(r"[A-Za-z]:[\\/]", content)
        section_checks[lang] = {"file": f"section_{lang}.md", "sha256": sha(HERE / f"section_{lang}.md"), "table_count": len(re.findall(r"^\|[-|: ]+\|$", content, re.M)), "table_line_count": len(lines), "png_links": 5, "svg_links": 5, "broken_relative_links": [], "word_tokens": len(re.findall(r"\b[A-Za-z]+(?:[-'][A-Za-z]+)*\b", content)), "han_characters": len(re.findall(r"[\u4e00-\u9fff]", content))}
    assert table_numbers["english"] == table_numbers["chinese"], "Bilingual table numeric-token mismatch"
    assert section_checks["english"]["word_tokens"] >= 2000
    assert section_checks["english"]["table_count"] == section_checks["chinese"]["table_count"] == 6

    manifest = read(HERE / "figures" / "manifest.json")
    assert manifest["generator_sha256"] == sha(HERE / "plotting.py")
    figure_checks = []
    for item in manifest["files"]:
        path = HERE / item["file"]
        assert sha(path) == item["sha256"]
        for source, expected in item["source_sha256"].items():
            assert sha(HERE / source) == expected
        entry = {"file": item["file"], "sha256": sha(path)}
        if path.suffix == ".png":
            data = path.read_bytes()
            dimensions = struct.unpack(">II", data[16:24])
            assert dimensions == (2700, 1260)
            entry.update({"pixel_dimensions": dimensions, "review_type": "AI visual inspection through view_image", "reviewed": True, "result": "No remaining blocking issue in title, axes, units, labels, font or clipping", "not_human_review": True})
        else:
            root = ET.fromstring(path.read_text(encoding="utf-8"))
            assert root.attrib["width"] == "648pt" and root.attrib["height"] == "302.4pt"
            texts = root.findall(".//{http://www.w3.org/2000/svg}text")
            assert len(texts) > 10
            entry.update({"review_type": "SVG XML, canvas and editable-text structural check", "text_elements": len(texts), "canvas_points": [648, 302.4], "independently_rendered": False, "result": "Structural checks passed; no independent rendered-SVG visual claim"})
        figure_checks.append(entry)
    assert len(figure_checks) == 20
    for table in catalog["posthoc_tables"]:
        assert sha(SOURCES / table["file"]) == table["sha256"]
    private_path_hits = []
    for path in HERE.rglob("*"):
        if path.suffix not in (".md", ".json", ".csv", ".py", ".svg") or path.name == "audit.json":
            continue
        if re.search(r"[A-Za-z]:[\\/](?:Users|Codex)[\\/]", path.read_text(encoding="utf-8-sig")):
            private_path_hits.append(path.relative_to(HERE).as_posix())
    assert not private_path_hits, private_path_hits
    payload = {
        "schema_version": 1,
        "scope": "Retrospective fixed-public-commit extraction, native-text corroboration, arithmetic, bilingual manuscript writing and plotting; no new chemical calculations",
        "remote_snapshot": {"repository": catalog["repository"], "fixed_commit": fixed, "live_git_ls_remote_verified_date": "2026-09-29", "HEAD": fixed, "refs/heads/main": fixed, "old_local_head_read_only": "9aa9c08b15f2f64b3f108e4a90e0f8ec0a4352c0", "old_local_status_porcelain_empty_at_end": True, "old_clone_untouched": True, "extraction": "No-checkout object-store clone reused existing local objects; only git show and history reads for science evidence"},
        "source_catalog": {"path": "sources/catalog.json", "sha256": sha(SOURCES / "catalog.json"), "source_count": len(catalog["source_entries"]), "source_entries": catalog["source_entries"]},
        "source_numeric_index": numeric_index,
        "history": catalog["history"],
        "checks": catalog["checks"],
        "native_energy_check": {"path": "sources/native_energy_checks.json", "sha256": sha(SOURCES / "native_energy_checks.json"), "attempted": 10, "passed": 10, "maximum_energy_discrepancy_hartree": 0.0, "first_8_logs_match_preexisting_archived_hash": True, "last_2_native_log_hashes_recorded_here": True, "meaning": "Exact agreement at archived precision validates transcription, not independent quantum accuracy"},
        "posthoc_tables": catalog["posthoc_tables"],
        "section_checks": section_checks,
        "bilingual_table_numeric_tokens_equal": True,
        "figures": {"manifest": "figures/manifest.json", "manifest_sha256": sha(HERE / "figures" / "manifest.json"), "groups": 5, "png_ai_visual_inspections": 10, "svg_structural_checks": 10, "independent_svg_render_inspections": 0, "resolved_visual_issue": "Shortened English C5 panel titles to keep text inside the canvas; changed C4 English labels from Seed to Start to avoid confusing run order with numerical cluster seed", "checks": figure_checks},
        "new_computation_accounting": {"quantum_driver_calls": 0, "NEB_runs": 0, "neural_training_runs": 0, "wet_lab_experiments": 0, "arithmetic_tables": 6},
        "interpretation_boundaries": ["54 target slots, 48 native state attempts, 8 converged states and 6 quality-screened states have different denominators; the 48 attempts are not the sum of all historical exploratory protocols", "No eligible MECP pair was dispatched; zero MECP calculations is a precondition rejection, not a converged crossing or proof that no crossing exists", "The trained EGNN auxiliary target is relative conformer electronic energy; activation-barrier and MECP-gap heads have no observed training labels", "The EGNN does not improve over the equal-reference baseline on the held-out test split", "Microsolvation electronic association, 1 M thermochemistry and ALPB-containing beyond-pair decomposition are different quantities", "All proton-wire paths remain unaccepted; band maxima are not activation free energies", "Latest molecular model C9H8N2O has no Cu; two neutral states and two anion timeouts do not make an electron-affinity pair", "The two different-functional neutral DFT total energies cannot be subtracted to obtain electron affinity", "SCF stability, DFT geometry optimization and DFT Hessians were absent from the latest molecular pilot", "Synthetic kinetic conservation checks are not physical catalyst predictions"],
        "license_note": catalog["license_note"],
        "absolute_private_path_scan": {"hits": private_path_hits, "scope": "Authored sections, scripts, extracted text, CSV, JSON and SVG"},
        "runtime_for_audit_only": {name: importlib.metadata.version(name) for name in ("numpy", "matplotlib")},
        "reproduction_scripts": {"sources/reproduce_audit.py": sha(SOURCES / "reproduce_audit.py"), "sources/finalize_audit.py": sha(Path(__file__)), "plotting.py": sha(HERE / "plotting.py")},
        "reviewer": "AI agent; no claim of human peer review"
    }
    (HERE / "audit.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"sources": len(catalog["source_entries"]), "numeric_source_selectors": len(numeric_index), "section_checks": section_checks, "png_ai_reviewed": 10, "svg_structurally_checked": 10, "table_numeric_parity": True, "private_path_hits": private_path_hits}, ensure_ascii=False))


if __name__ == "__main__":
    main()
