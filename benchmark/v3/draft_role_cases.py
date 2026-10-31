#!/usr/bin/env python3
"""Draft pending-review role-truth registry cases from negation coverage.

Negative truth comes only from in-scope negative statements retained in the
coverage input.  Positive truth comes only from experimental Swiss-Prot GO or
EC annotations.  The result is a review queue, not accepted training truth.
"""
from __future__ import annotations

import argparse
from collections import defaultdict, deque
from datetime import date
import hashlib
import json
from pathlib import Path
import re
import sys

import evidence_codes as EV
import mine_negative_function_statements as MINER
import prepare_role_truth_loci as REGISTRY


SCHEMA = "gspa.role-truth-drafts/1"
DEFAULT_CLASSES = ("activity_negative", "gaf_not")
# HTP is deliberately absent: the drafting contract names the individual GO
# codes accepted here rather than adopting every high-throughput code.
EXPERIMENTAL_GO = frozenset({
    "EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HDA", "HMP", "HGI", "HEP",
})
EXPERIMENTAL_EC = EV.EXPERIMENTAL_ECO
PUBMED_RE = re.compile(r"(?:PubMed|PMID):(\d+)", re.I)
OUTCOMES = ("drafted_cases", "drafted_null_groups", "drafted_negative_only")
REJECTED_OUTCOME = "drafted_rejected"
FILENAMES = {name: name + ".jsonl" for name in (*OUTCOMES, REJECTED_OUTCOME)}
COMPLETE_EC = re.compile(r"^EC:\d+\.\d+\.\d+\.\d+$")
EC2GO_LINE = re.compile(r"^(EC:\S+)\s+>\s+GO:.*?;\s*(GO:\d+)\s*$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def file_info(path: Path) -> dict:
    path = Path(path)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def prepare_output(path: Path) -> None:
    path = Path(path)
    if path.exists():
        if not path.is_dir():
            raise FileExistsError(f"output path is not a directory: {path}")
        if any(path.iterdir()):
            raise FileExistsError(f"refusing to overwrite non-empty directory: {path}")
    else:
        path.mkdir(parents=True)


def read_coverage(path: Path, file_index: int) -> list[tuple[tuple[int, int], dict]]:
    records = []
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError as error:
                raise ValueError(f"invalid JSON on {path}:{line_number}: {error}") from error
            if not isinstance(row, dict):
                raise ValueError(f"expected an object on {path}:{line_number}")
            records.append(((file_index, line_number), row))
    return records


def load_go_parents(path: Path | None) -> dict[str, set[str]]:
    """Load transitive-is_a-ready child-to-parent links from a small OBO parser."""
    if path is None:
        return {}
    parents: dict[str, set[str]] = defaultdict(set)
    current = None
    in_term = False

    def flush():
        if current and current.get("id"):
            parents[current["id"]].update(current.get("parents", ()))

    with Path(path).open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith("["):
                flush()
                current = None
                in_term = line == "[Term]"
            elif not in_term or not line or line.startswith("!"):
                continue
            elif line.startswith("id:"):
                current = {"id": line[3:].strip(), "parents": []}
            elif current is not None and line.startswith("is_a:"):
                parent = line[5:].split("!", 1)[0].strip()
                if parent:
                    current["parents"].append(parent)
    flush()
    return {child: set(values) for child, values in parents.items()}


def load_ec2go(path: Path | None) -> dict[str, set[str]]:
    """Load GO mappings for complete four-field EC identifiers only."""
    if path is None:
        return {}
    mappings = defaultdict(set)
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith("!"):
                continue
            match = EC2GO_LINE.fullmatch(line)
            if match is None:
                raise ValueError(f"invalid ec2go line on {path}:{line_number}")
            ec, go_id = match.groups()
            if COMPLETE_EC.fullmatch(ec):
                mappings[ec].add(go_id)
    return {ec: set(go_ids) for ec, go_ids in mappings.items()}


def is_descendant(term: str, ancestor: str, parents: dict[str, set[str]]) -> bool:
    if term == ancestor:
        return False
    pending = deque(parents.get(term, ()))
    seen = set()
    while pending:
        current = pending.popleft()
        if current == ancestor:
            return True
        if current in seen:
            continue
        seen.add(current)
        pending.extend(parents.get(current, ()))
    return False


def positive_evidence(entry: dict | None, function_type: str, function_id: str,
                      go_parents: dict[str, set[str]],
                      ec2go: dict[str, set[str]] | None = None) -> dict | None:
    """Return one deterministic experimental annotation, preferring exact GO."""
    if entry is None:
        return None
    accession = entry["accessions"][0]
    if function_type == "EC":
        matches = sorted({(ec, evidence) for ec, evidence in entry.get("ec", [])
                          if ec == function_id and evidence in EXPERIMENTAL_EC})
        if matches:
            ec, evidence = matches[0]
            return {"basis": f"Swiss-Prot {accession} {ec} {evidence}",
                    "evidence_code": evidence}
        return None

    exact = set()
    descendants = set()
    for item in entry.get("go_mf", []):
        go_id, evidence = item.get("go_id"), item.get("evidence")
        if evidence not in EXPERIMENTAL_GO:
            continue
        if go_id == function_id:
            exact.add((go_id, evidence))
        elif go_parents and is_descendant(go_id, function_id, go_parents):
            descendants.add((go_id, evidence))
    if exact:
        go_id, evidence = sorted(exact)[0]
        basis = f"Swiss-Prot {accession} {go_id} {evidence}"
    elif descendants:
        go_id, evidence = sorted(descendants)[0]
        basis = (f"Swiss-Prot {accession} {go_id} {evidence} "
                 f"(descendant of {function_id})")
    else:
        mapped_exact = set()
        mapped_descendants = set()
        for ec, evidence in entry.get("ec", []):
            if evidence not in EXPERIMENTAL_EC:
                continue
            ec_number = ec[3:] if ec.startswith("EC:") else ec
            for go_id in (ec2go or {}).get("EC:" + ec_number, ()):
                item = (ec_number, evidence, go_id)
                if go_id == function_id:
                    mapped_exact.add(item)
                elif go_parents and is_descendant(go_id, function_id, go_parents):
                    mapped_descendants.add(item)
        if mapped_exact:
            ec_number, evidence, go_id = sorted(mapped_exact)[0]
            suffix = f"(ec2go → {go_id})"
        elif mapped_descendants:
            ec_number, evidence, go_id = sorted(mapped_descendants)[0]
            suffix = f"(ec2go → {go_id}, descendant of {function_id})"
        else:
            return None
        return {
            "basis": f"Swiss-Prot {accession} EC {ec_number} {evidence} {suffix}",
            "evidence_code": evidence,
            "positive_source": "ec2go",
        }
    return {"basis": basis, "evidence_code": evidence}


def numeric_key(value: str):
    try:
        return (0, int(value), value)
    except (TypeError, ValueError):
        return (1, 0, str(value))


def statement_references(row: dict) -> list[str]:
    """Collect normalized source ids available on a coverage row."""
    pmids = {str(value) for value in row.get("pubmed_ids", [])}
    pmids.update(PUBMED_RE.findall(row.get("sentence", "") or ""))
    references = []
    for raw in row.get("references", []) or []:
        value = str(raw).strip()
        match = PUBMED_RE.fullmatch(value)
        if match:
            pmids.add(match.group(1))
        elif value.startswith(("DOI:", "PMC")):
            references.append(value)
        else:
            pmids.update(PUBMED_RE.findall(value))
    result = ["PMID:" + value for value in sorted(pmids, key=numeric_key)]
    result.extend(sorted(set(references)))
    return result


def negative_basis(row: dict) -> str:
    sentence = re.sub(r"\s+", " ", row.get("sentence", "") or "").strip()
    pmids = [ref for ref in statement_references(row) if ref.startswith("PMID:")]
    return sentence + (" [" + ", ".join(pmids) + "]" if pmids else "")


def skip_description(ref: tuple[int, int], row: dict, coverage_paths: list[Path], reason: str) -> dict:
    return {
        "coverage": str(coverage_paths[ref[0]]),
        "line": ref[1],
        "statement_index": row.get("statement_index"),
        "genome_id": row.get("genome_id", ""),
        "function_id": row.get("function_id", ""),
        "accession": row.get("accession", ""),
        "protein_id": row.get("protein_id", ""),
        "reason": reason,
    }


def case_id(genome_id: str, function_id: str) -> str:
    return f"{genome_id}-{function_id.replace(':', '-').replace('.', '-')}".lower()


def markdown_text(value: str) -> str:
    return re.sub(r"\s+", " ", value or "").replace("|", "\\|").strip()


def truncate_sentence(value: str, limit: int = 120) -> str:
    value = markdown_text(value)
    return value if len(value) <= limit else value[:limit - 3].rstrip() + "..."


def draft_case(group: dict, entries: dict[str, dict], go_parents: dict[str, set[str]],
               added: str, ec2go: dict[str, set[str]] | None = None
               ) -> tuple[dict | None, str | None]:
    universe = group["universe"]
    negative_rows = group["negative_rows"]
    negative_members = set(negative_rows)
    negative_entries = [entries[member] for member in sorted(negative_members) if member in entries]
    if not negative_entries:
        return None, "negative_member_swissprot_entry_missing"
    source_entry = negative_entries[0]
    if not source_entry.get("taxid") or not source_entry.get("organism"):
        return None, "negative_member_taxid_or_organism_missing"

    references = sorted({ref for rows in negative_rows.values() for _, row in rows
                         for ref in statement_references(row)})
    if not any(ref.startswith(("PMID:", "DOI:", "PMC")) for ref in references):
        return None, "no_primary_reference"

    function_id = group["function_id"]
    members = []
    conflict = False
    sentences = []
    for protein in universe:
        entry = entries.get(protein)
        positive = positive_evidence(
            entry, group["function_type"], function_id, go_parents, ec2go=ec2go)
        bases = []
        if protein in negative_rows:
            unique = {}
            for _, row in negative_rows[protein]:
                sentence = row.get("sentence", "") or ""
                sentences.append(sentence)
                basis = negative_basis(row)
                unique[(sentence, tuple(statement_references(row)))] = basis
            bases.extend(unique[key] for key in sorted(unique))
        if positive is not None:
            bases.append(positive["basis"])
        if protein in negative_rows and positive is not None:
            label, role_truth, conflict = None, "unknown", True
        elif protein in negative_rows:
            label, role_truth = 0, "inactive_homolog"
        elif positive is not None:
            label, role_truth = 1, "catalytic_enzyme"
        else:
            label, role_truth = None, "unknown"
            bases = [f"no experimental annotation or negative statement located for {function_id}"]
        member = {
            "protein_id": protein,
            "swissprot_accession": protein,
            "label": label,
            "role_truth": role_truth,
            "basis": " | ".join(bases),
        }
        if label == 1 and positive is not None:
            member["evidence_code"] = positive["evidence_code"]
        if positive is not None and positive.get("positive_source"):
            member["positive_source"] = positive["positive_source"]
        members.append(member)

    positives = sum(member["label"] == 1 for member in members)
    negatives = sum(member["label"] == 0 for member in members)
    unknowns = sum(member["label"] is None for member in members)
    if negatives == 0:
        return None, "no_clean_negative_after_conflict_resolution"
    if positives:
        outcome = "drafted_cases"
    elif unknowns == 0:
        outcome = "drafted_null_groups"
    else:
        outcome = "drafted_negative_only"

    unique_sentences = list(dict.fromkeys(sentences))
    case = {
        "case_id": case_id(group["genome_id"], function_id),
        "status": "pending_review",
        "genome_id": group["genome_id"],
        "taxid": str(source_entry["taxid"]),
        "organism": source_entry["organism"],
        "in_training_cohort": False,
        "function_type": group["function_type"],
        "function_id": function_id,
        "role": "UNSPECIFIED" if group["function_type"] == "GO" else "CATALYTIC",
        "members": members,
        "references": references,
        "source_statement": "Negative-statement draft: " + " / ".join(
            negative_basis(row) for rows in negative_rows.values() for _, row in rows
        ),
        "reference_files": [],
        "universe_note": (f"Coverage-derived universe of {len(universe)} candidates for {function_id}; "
                          "automated evidence draft, completeness and labels require human review"),
        "added": added,
        "reviewer": f"draft_role_cases.py {added}; human review pending",
    }
    if positives > 1:
        case["allow_multi_positive"] = True
    if conflict:
        case["conflict"] = True
    if outcome == "drafted_null_groups":
        case["null_group_candidate"] = True
    return (case, outcome), " / ".join(unique_sentences)


def write_jsonl(path: Path, cases: list[dict]) -> None:
    with Path(path).open("w", encoding="utf-8") as handle:
        for case in cases:
            handle.write(json.dumps(case, sort_keys=True) + "\n")


def run(args: argparse.Namespace) -> dict:
    coverage_paths = [Path(path) for path in args.coverage]
    swissprot = Path(args.swissprot)
    out = Path(args.out)
    if args.max_universe < 2:
        raise ValueError("--max-universe must be at least 2")
    classes = {value.strip() for value in args.classes.split(",") if value.strip()}
    if not classes:
        raise ValueError("--classes must select at least one statement class")
    prepare_output(out)

    records = []
    for file_index, path in enumerate(coverage_paths):
        records.extend(read_coverage(path, file_index))
    skipped = {}
    groups = {}
    genomes = set()

    def skip(ref, row, reason):
        skipped.setdefault(ref, (row, reason))

    for ref, row in records:
        genome = str(row.get("genome_id") or "")
        if genome:
            genomes.add(genome)
        if row.get("statement_class") not in classes:
            skip(ref, row, "statement_class_not_selected")
            continue
        is_gaf = row.get("kind") == "gaf_not" or row.get("statement_class") == "gaf_not"
        if is_gaf and not args.include_bp_cc and row.get("aspect") != "F":
            skip(ref, row, "gaf_aspect_not_F")
            continue
        if not row.get("in_negation_scope"):
            skip(ref, row, "not_in_negation_scope")
            continue
        if not row.get("protein_in_universe"):
            skip(ref, row, "negative_protein_not_in_universe")
            continue
        members = tuple(sorted({str(member) for member in row.get("universe_members", [])}))
        try:
            declared_size = int(row.get("universe_size", len(members)))
        except (TypeError, ValueError):
            skip(ref, row, "invalid_universe_size")
            continue
        if not genome or row.get("function_type") not in ("GO", "EC") or not row.get("function_id"):
            skip(ref, row, "missing_genome_or_function_identity")
            continue
        if row.get("truncated"):
            skip(ref, row, "universe_truncated")
            continue
        if declared_size > args.max_universe or len(members) > args.max_universe:
            skip(ref, row, "universe_exceeds_max_universe")
            continue
        if declared_size != len(members) or len(members) < 2:
            skip(ref, row, "incomplete_or_singleton_universe")
            continue
        protein = str(row.get("protein_id") or "")
        if protein not in members:
            skip(ref, row, "negative_protein_absent_from_universe_members")
            continue
        key = (genome, row["function_type"], row["function_id"])
        group = groups.setdefault(key, {
            "genome_id": genome,
            "function_type": row["function_type"],
            "function_id": row["function_id"],
            "universe": members,
            "negative_rows": defaultdict(list),
            "rows": [],
            "inconsistent_universe": False,
        })
        if group["universe"] != members:
            group["inconsistent_universe"] = True
        group["negative_rows"][protein].append((ref, row))
        group["rows"].append((ref, row))

    # With taxid selection, retain selected taxa.  Otherwise retain exactly the
    # primary accessions occurring in a kept universe, in one streaming pass.
    requested_taxids = MINER.read_taxids(args.taxids) if args.taxids else None
    selected_accessions = set()
    if requested_taxids is None:
        selected_accessions = {member for group in groups.values() for member in group["universe"]}
    entries = {}
    entries_scanned = 0
    for entry in MINER.iter_entries(swissprot):
        entries_scanned += 1
        accession = entry["accessions"][0]
        selected = (str(entry.get("taxid") or "") in requested_taxids
                    if requested_taxids is not None else accession in selected_accessions)
        if selected:
            if accession in entries:
                raise ValueError(f"duplicate Swiss-Prot primary accession: {accession}")
            entries[accession] = entry

    go_obo = Path(args.go_obo) if args.go_obo else None
    go_parents = load_go_parents(go_obo)
    ec2go_path = Path(args.ec2go) if args.ec2go else None
    ec2go = load_ec2go(ec2go_path)
    added = date.today().isoformat()
    routed = {name: [] for name in FILENAMES}
    review_rows = []
    rejected_review_rows = []
    for key in sorted(groups):
        group = groups[key]
        if group["inconsistent_universe"]:
            for ref, row in group["rows"]:
                skip(ref, row, "inconsistent_universe_for_case")
            continue
        result, detail = draft_case(group, entries, go_parents, added, ec2go=ec2go)
        if result is None:
            for ref, row in group["rows"]:
                skip(ref, row, detail)
            continue
        drafted, outcome = result
        if any(row.get("proposed_function", {}).get("mapping_too_broad")
               for rows in group["negative_rows"].values() for _, row in rows):
            outcome = REJECTED_OUTCOME
            drafted["reason"] = "mapping_too_broad"
        routed[outcome].append(drafted)
        if outcome == REJECTED_OUTCOME:
            rejected_review_rows.append((drafted, detail))
        else:
            review_rows.append((drafted, detail))

    for outcome in FILENAMES:
        routed[outcome].sort(key=lambda case: case["case_id"])
        write_jsonl(out / FILENAMES[outcome], routed[outcome])

    by_genome = {}
    for genome in sorted(genomes):
        counts = {outcome: sum(case["genome_id"] == genome for case in routed[outcome])
                  for outcome in OUTCOMES}
        rejected = sum(case["genome_id"] == genome for case in routed[REJECTED_OUTCOME])
        if rejected:
            counts[REJECTED_OUTCOME] = rejected
        counts["total"] = sum(counts.values())
        by_genome[genome] = counts
    by_outcome = {outcome: len(routed[outcome]) for outcome in OUTCOMES}
    if routed[REJECTED_OUTCOME]:
        by_outcome[REJECTED_OUTCOME] = len(routed[REJECTED_OUTCOME])
    by_outcome["total"] = sum(by_outcome.values())

    with (out / "draft_review.md").open("w", encoding="utf-8") as handle:
        handle.write("# Draft role-truth cases for human review\n")
        for genome in sorted(genomes):
            handle.write(f"\n## {markdown_text(genome)}\n\n")
            handle.write("| case_id | function | members (+/−/?) | negative sentence |\n")
            handle.write("|---|---|---|---|\n")
            for drafted, sentences in sorted(review_rows, key=lambda item: item[0]["case_id"]):
                if drafted["genome_id"] != genome:
                    continue
                symbols = " ".join(
                    f"{markdown_text(member['protein_id'])}("
                    f"{ '+' if member['label'] == 1 else '−' if member['label'] == 0 else '?' })"
                    for member in drafted["members"]
                )
                handle.write(
                    f"| {markdown_text(drafted['case_id'])} | {markdown_text(drafted['function_id'])}"
                    f" | {symbols} | {truncate_sentence(sentences)} |\n"
                )
        if rejected_review_rows:
            handle.write("\n## Rejected\n\n")
            handle.write("| case_id | function | members (+/−/?) | negative sentence | reason |\n")
            handle.write("|---|---|---|---|---|\n")
            for drafted, sentences in sorted(
                    rejected_review_rows, key=lambda item: item[0]["case_id"]):
                symbols = " ".join(
                    f"{markdown_text(member['protein_id'])}("
                    f"{ '+' if member['label'] == 1 else '−' if member['label'] == 0 else '?' })"
                    for member in drafted["members"]
                )
                handle.write(
                    f"| {markdown_text(drafted['case_id'])} | {markdown_text(drafted['function_id'])}"
                    f" | {symbols} | {truncate_sentence(sentences)} | {drafted['reason']} |\n"
                )

    # Validate the actual files, using the registry's existing negative-only
    # relaxation for both files that are allowed to lack a positive.
    validation = {
        FILENAMES["drafted_cases"]: REGISTRY.validate_registry(
            out / FILENAMES["drafted_cases"], go_obo=go_obo),
        FILENAMES["drafted_null_groups"]: REGISTRY.validate_registry(
            out / FILENAMES["drafted_null_groups"], go_obo=go_obo, require_positive=False),
        FILENAMES["drafted_negative_only"]: REGISTRY.validate_registry(
            out / FILENAMES["drafted_negative_only"], go_obo=go_obo, require_positive=False),
    }
    inputs = {
        "coverage": [file_info(path) for path in coverage_paths],
        "swissprot": file_info(swissprot),
    }
    if args.taxids:
        inputs["taxids"] = file_info(args.taxids)
    if go_obo:
        inputs["go_obo"] = file_info(go_obo)
    if ec2go_path:
        inputs["ec2go"] = file_info(ec2go_path)
    cli_args = {
        "coverage": [str(path) for path in coverage_paths],
        "swissprot": str(swissprot),
        "taxids": str(args.taxids) if args.taxids else None,
        "accessions_from_coverage": requested_taxids is None,
        "go_obo": str(go_obo) if go_obo else None,
        "ec2go": str(ec2go_path) if ec2go_path else None,
        "max_universe": args.max_universe,
        "classes": sorted(classes),
        "include_bp_cc": bool(args.include_bp_cc),
        "out": str(out),
    }
    skipped_rows = [skip_description(ref, skipped[ref][0], coverage_paths, skipped[ref][1])
                    for ref, _ in records if ref in skipped]
    report = {
        "schema": SCHEMA,
        "status": "verified",
        "date": added,
        "inputs": inputs,
        "cli_args": cli_args,
        "counts": {
            "coverage_records": len(records),
            "kept_statements": len(records) - len(skipped_rows),
            "skipped_statements": len(skipped_rows),
            "swissprot_entries_scanned": entries_scanned,
            "swissprot_entries_retained": len(entries),
            "by_genome": by_genome,
            "by_outcome": by_outcome,
        },
        "skipped_statements": skipped_rows,
        "validation": validation,
        "output_sha256": {
            name: sha256(out / name) for name in [*FILENAMES.values(), "draft_review.md"]
        },
    }
    (out / "draft_report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                                           encoding="utf-8")
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coverage", type=Path, action="append", required=True)
    parser.add_argument("--swissprot", type=Path, required=True)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--taxids", type=Path)
    selection.add_argument("--accessions-from-coverage", action="store_true",
                           help="retain only universe-member primary accessions (default)")
    parser.add_argument("--go-obo", type=Path)
    parser.add_argument("--ec2go", type=Path)
    parser.add_argument("--max-universe", type=int, default=8)
    parser.add_argument("--classes", default=",".join(DEFAULT_CLASSES))
    parser.add_argument("--include-bp-cc", action="store_true")
    parser.add_argument("--out", type=Path, required=True)
    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        report = run(args)
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError,
            REGISTRY.RoleTruthError) as error:
        parser.exit(2, f"error: {error}\n")
    print(json.dumps(report["counts"]["by_outcome"], sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
