"""Render numbered, bilingual-consistent citations from a verified registry.

Standard-library-only publishing operation; no scientific calculations, network
requests, source-manuscript mutations, or journal-tier inferences are performed.
Run from any directory. --check compares expected output bytes without replacing
them; the citation audit is refreshed in either mode.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
from urllib.parse import unquote

BASE = Path(__file__).resolve().parents[1]
REPO = BASE.parent.parent
SOURCES = (
    "manuscript_english.md",
    "manuscript_chinese.md",
)
REGISTRY = BASE / "references_verified.json"
AUDIT = BASE / "results" / "citation_audit.json"
CITE = re.compile(r"\[\s*(@[^\]\n]+)\]")
KEY = re.compile(r"@[A-Za-z0-9_:-]+")
LINK = re.compile(r"(!?\[[^\]\n]*\])\(([^\)\n]+)\)")
PLACEHOLDER = "<!-- VERIFIED_REFERENCES -->"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def relative(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def groups(text: str) -> list[list[str]]:
    """Accept the deliberately small source syntax [@key; @key]."""
    result = []
    for match in CITE.finditer(text):
        fields = [part.strip() for part in match.group(1).split(";")]
        if not all(KEY.fullmatch(part) for part in fields):
            raise ValueError(f"Unsupported citation syntax: {match.group(0)}")
        result.append([part[1:] for part in fields])
    if "[@" in CITE.sub("", text):
        raise ValueError("Malformed or unmatched [@ citation marker")
    return result


def flatten(items: list[list[str]]) -> list[str]:
    return [key for group in items for key in group]


def distinct(items: list[str]) -> list[str]:
    return list(dict.fromkeys(items))


def entry(ref: dict) -> str:
    authors = "; ".join(ref["authors"])
    if not ref.get("author_list_complete", False):
        authors += "; et al."
    location = str(ref["year"])
    if ref.get("volume"):
        location += ", " + str(ref["volume"])
    if ref.get("issue"):
        location += "(" + str(ref["issue"]) + ")"
    if ref.get("pages_or_article"):
        location += ", " + str(ref["pages_or_article"])
    link = (
        f"[doi:{ref['doi']}](https://doi.org/{ref['doi']})"
        if ref.get("doi") else f"[Publisher record]({ref['url']})"
    )
    return f"{authors.rstrip('.')}. {ref['title']}. *{ref['journal']}* **{location}**. {link}"


def bibliography(refs: dict[str, dict], keys: list[str], numbers: dict[str, int]) -> str:
    return "\n\n".join(f"[{numbers[key]}] {entry(refs[key])}" for key in keys)


def tex_escape(value: object) -> str:
    value = str(value)
    # UTF-8 names are retained, avoiding invented transliterations.
    return value.replace("&", r"\&").replace("%", r"\%").replace("_", r"\_")


def bibtex(ref: dict) -> str:
    conference = ref.get("publication_type") == "peer_reviewed_conference_paper"
    kind = "inproceedings" if conference else "article"
    suffix_names = {
        "Thom H. Dunning, Jr.": "Dunning, Jr., Thom H.",
        "Tucker Carrington, Jr.": "Carrington, Jr., Tucker",
    }
    author_names = [suffix_names.get(name, name) for name in ref["authors"]]
    if not ref.get("author_list_complete", False):
        author_names.append("others")
    fields = {
        "author": " and ".join(author_names),
        "title": "{" + ref["title"] + "}",
        "booktitle" if conference else "journal": ref["journal"],
        "year": str(ref["year"]),
    }
    for original, target in (("volume", "volume"), ("issue", "number"),
                             ("pages_or_article", "pages"), ("doi", "doi"), ("url", "url")):
        if ref.get(original):
            fields[target] = str(ref[original]).replace("–", "--") if target == "pages" else ref[original]
    notes = []
    if not ref.get("author_list_complete", False):
        notes.append("Author list deliberately abbreviated to the verified names; et al. retained")
    if not ref.get("doi"):
        notes.append("No DOI provided in the inspected official publication record")
    if notes:
        fields["note"] = "; ".join(notes)
    body = ",\n".join(f"  {name} = {{{tex_escape(value)}}}" for name, value in fields.items())
    return f"@{kind}{{{ref['citation_key']},\n{body}\n}}"


def rebase_links(text: str, destination: Path) -> tuple[str, list[dict]]:
    checked = []

    def convert(match: re.Match) -> str:
        label, target = match.groups()
        angle = target.startswith("<") and target.endswith(">")
        path_text = target[1:-1] if angle else target
        if re.match(r"^[a-zA-Z][a-zA-Z0-9+.-]*:", path_text) or path_text.startswith(("#", "//")):
            return match.group(0)
        # Existing manuscripts use simple local targets without optional titles.
        bare, marker, anchor = path_text.partition("#")
        resolved = (BASE / unquote(bare)).resolve()
        exists = resolved.exists()
        inside = resolved.is_relative_to(REPO)
        checked.append({"source_target": path_text, "exists": exists, "inside_repository": inside})
        rebased = Path(os.path.relpath(resolved, destination.parent)).as_posix()
        if marker:
            rebased += "#" + anchor
        return f"{label}(<{rebased}>)" if angle else f"{label}({rebased})"

    return LINK.sub(convert, text), checked


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="verify rendered bytes without replacing them")
    args = parser.parse_args()
    failures: list[str] = []
    warnings: list[str] = []
    audit: dict = {
        "schema_version": 1,
        "scope": "Bibliographic identity and rendering integrity, not scientific or journal-tier validation.",
        "numbering_rule": "First appearance in expanded English main text; any Chinese-only keys are diagnosed. Uncited registry entries follow registry order in the complete bibliography only.",
        "no_new_scientific_calculations": True,
        "failures": failures,
        "warnings": warnings,
    }
    try:
        raw_registry = REGISTRY.read_bytes()
        registry = json.loads(raw_registry.decode("utf-8-sig", errors="strict"))
        records = registry["references"]
        refs = {r["citation_key"]: r for r in records}
        if len(refs) != len(records):
            failures.append("Duplicate citation keys in reference registry")
        ids = [r["id"] for r in records]
        if len(set(ids)) != len(ids):
            failures.append("Duplicate registry IDs")
        dois = [r["doi"].casefold() for r in records if r.get("doi")]
        if len(set(dois)) != len(dois):
            failures.append("Duplicate non-null DOIs")
        if registry.get("reference_count") != len(records):
            failures.append("Registry reference_count disagrees with its rows")
        bibliographic_unverified = [r["citation_key"] for r in records
                                   if r.get("publication_status") != "published" or not r.get("verified_source")]
        if bibliographic_unverified:
            failures.append("Unverified publication records: " + ", ".join(bibliographic_unverified))

        source_bytes = {name: (BASE / name).read_bytes() for name in SOURCES}
        source_text = {name: value.decode("utf-8-sig", errors="strict") for name, value in source_bytes.items()}
        parsed = {name: groups(value) for name, value in source_text.items()}
        all_keys = distinct([key for name in SOURCES for key in flatten(parsed[name])])
        dangling = sorted(set(all_keys) - set(refs))
        if dangling:
            raise ValueError("Dangling citation keys: " + ", ".join(dangling))
        duplicated_groups = {name: [group for group in values if len(set(group)) != len(group)]
                             for name, values in parsed.items()}
        if any(duplicated_groups.values()):
            failures.append("Repeated key inside a citation group")
        parity = {}
        for en, cn in ((SOURCES[0], SOURCES[1]),):
            parity[en] = {"chinese_source": cn, "same_group_order": parsed[en] == parsed[cn],
                          "same_unique_keys": set(flatten(parsed[en])) == set(flatten(parsed[cn]))}
            if not all((parity[en]["same_group_order"], parity[en]["same_unique_keys"])):
                failures.append("English/Chinese citation mismatch: " + en)

        cited = distinct(flatten(parsed[SOURCES[0]]))
        chinese_only = [key for key in all_keys if key not in cited]
        cited += chinese_only
        unused = [r["citation_key"] for r in records if r["citation_key"] not in cited]
        ordered = cited + unused
        numbers = {key: i + 1 for i, key in enumerate(ordered)}
        outputs: dict[Path, str] = {}
        documents = {}
        for name, text in source_text.items():
            dest = BASE / "rendered" / name
            rebased, local_links = rebase_links(text, dest)

            def replace_citation(match: re.Match) -> str:
                keys = groups(match.group(0))[0]
                return "[" + ", ".join(str(numbers[key]) for key in keys) + "]"

            rendered = CITE.sub(replace_citation, rebased)
            doc_keys = sorted(set(flatten(parsed[name])), key=numbers.__getitem__)
            ref_text = bibliography(refs, doc_keys, numbers)
            marker_count = rendered.count(PLACEHOLDER)
            if marker_count > 1:
                failures.append("Duplicate reference placeholder: " + name)
            if marker_count:
                rendered = rendered.replace(PLACEHOLDER, ref_text)
            else:
                heading = "## 参考文献" if "chinese" in name else "## References"
                rendered = rendered.rstrip() + "\n\n" + heading + "\n\n" + ref_text + "\n"
            rendered = rendered.replace("\r\n", "\n").rstrip() + "\n"
            if "[@" in rendered or PLACEHOLDER in rendered:
                failures.append("Unresolved citation or bibliography marker: " + name)
            missing_links = [link for link in local_links if not link["exists"] or not link["inside_repository"]]
            if missing_links:
                failures.append("Missing or outside-repository local targets: " + name)
            outputs[dest] = rendered
            documents[name] = {
                "source_sha256": sha(source_bytes[name]),
                "citation_groups": len(parsed[name]),
                "citation_occurrences": len(flatten(parsed[name])),
                "unique_cited_references": len(doc_keys),
                "reference_keys": doc_keys,
                "reference_numbers": [numbers[key] for key in doc_keys],
                "rendered_file": relative(dest),
                "local_link_count": len(local_links),
                "local_link_failures": missing_links,
            }

        cited_text = bibliography(refs, cited, numbers)
        uncited_text = bibliography(refs, unused, numbers)
        outputs[BASE / "bibliography.md"] = (
            "# Verified publication bibliography\n\n"
            "Numbers match the English/Chinese expanded main manuscripts. "
            "Only cited entries appear in each rendered document. Publication verification is "
            "separate from journal-tier verification: no CAS/JCR classification is asserted. "
            "Incomplete author lists retain et al.; absent DOIs are not invented.\n\n"
            "## Cited in the expanded manuscripts\n\n" + cited_text +
            "\n\n## Verified registry entries not cited in these documents\n\n" +
            (uncited_text or "None.") + "\n"
        )
        outputs[BASE / "references.bib"] = (
            "% All verified registry entries; numbering belongs to rendered Markdown.\n"
            "% No journal classification is inferred; incomplete author lists use and others.\n\n"
            + "\n\n".join(bibtex(ref) for ref in records) + "\n"
        )
        utf8_checks = {}
        for name, text in {**source_text, **{relative(p): t for p, t in outputs.items()},
                           "references_verified.json": raw_registry.decode("utf-8-sig")}.items():
            sentinels = [x for x in ("\ufffd", "â€", "Ã") if x in text]
            utf8_checks[name] = {"strict_utf8_decode": True, "suspicious_sequences": sentinels}
            if sentinels:
                failures.append("Encoding warning in " + name)

        output_hashes = {}
        for path, text in outputs.items():
            expected = text.encode("utf-8")
            if args.check:
                if not path.exists() or path.read_bytes() != expected:
                    failures.append("Missing or stale rendered output: " + relative(path))
            else:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(expected)
            output_hashes[relative(path)] = sha(expected)
        source_unchanged = all((BASE / name).read_bytes() == value for name, value in source_bytes.items())
        if not source_unchanged:
            failures.append("Source manuscript changed during rendering")
        audit.update({
            "registry_sha256": sha(raw_registry),
            "script_sha256": sha(Path(__file__).read_bytes()),
            "registry_references": len(records), "globally_cited_references": len(cited),
            "unused_verified_keys": unused, "dangling_keys": dangling,
            "duplicate_citation_groups": duplicated_groups,
            "unverified_publication_keys": bibliographic_unverified,
            "journal_classifications": {r["citation_key"]: r.get("journal_classification", {}).get("classification") for r in records},
            "partial_author_list_keys": [r["citation_key"] for r in records if not r.get("author_list_complete", False)],
            "null_doi_keys": [r["citation_key"] for r in records if not r.get("doi")],
            "global_numbering": [{"number": numbers[key], "key": key, "registry_id": refs[key]["id"], "cited": key in cited} for key in ordered],
            "language_parity": parity, "documents": documents,
            "source_markdown_unchanged_by_renderer": source_unchanged,
            "utf8_checks": utf8_checks,
            "output_sha256": output_hashes,
        })
    except (OSError, ValueError, KeyError, TypeError) as exc:
        failures.append(f"{type(exc).__name__}: {exc}")
    audit["passed"] = not failures
    audit["failure_count"] = len(failures)
    AUDIT.parent.mkdir(parents=True, exist_ok=True)
    AUDIT.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"passed": audit["passed"], "failure_count": len(failures),
                      "registry_references": audit.get("registry_references"),
                      "globally_cited_references": audit.get("globally_cited_references"),
                      "failures": failures}, ensure_ascii=False))
    return 0 if audit["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
