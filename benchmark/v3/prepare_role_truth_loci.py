#!/usr/bin/env python3
"""Validate reviewed role-truth cases and augment a training cohort with them.

These are complete, reviewed candidate loci with experimental role labels; they are
not a sufficient null-model training set. ``build-loci`` preserves every original
training example and split, adds only declared positive/negative truth, and certifies the
whole exported candidate universe. Evaluation-only cases are never built here.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import io
import json
import re
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "scripts"))
import audit_gaf_negatives as A  # noqa: E402
import match_swissprot_labels as M  # noqa: E402
import prepare_training_data as P  # noqa: E402

T = importlib.import_module("train_context_model")  # noqa: E402

SCHEMA = "gspa.role-truth-loci/1"
STATUSES = ("accepted", "pending_review", "rejected")
FUNCTION_TYPES = ("GO", "EC")
SPLITS = ("train", "calibration", "test")
GO_ID = re.compile(r"GO:\d{7}")
EC_ID = re.compile(r"EC:\d+\.\d+\.\d+\.\d+")
SHA256 = re.compile(r"[0-9a-f]{64}")


class RoleTruthError(ValueError):
    """A registry or cohort violates the reviewed-locus protocol."""


def sha256(path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_json(path):
    try:
        return json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise RoleTruthError(f"cannot read JSON {path}: {exc}") from exc


def _read_jsonl(path):
    rows = []
    try:
        for line in Path(path).read_text().splitlines():
            if line.strip():
                rows.append(json.loads(line))
    except (OSError, json.JSONDecodeError) as exc:
        raise RoleTruthError(f"cannot read JSONL {path}: {exc}") from exc
    return rows


def load_registry(path) -> list[dict]:
    """Read JSONL without ever consulting the pinned Swiss-Prot corpus."""
    cases = []
    with Path(path).open() as fh:
        for line_number, line in enumerate(fh, 1):
            if not line.strip():
                continue
            try:
                case = json.loads(line)
            except json.JSONDecodeError as exc:
                raise RoleTruthError(f"{path}:{line_number}: invalid JSON: {exc}") from exc
            if not isinstance(case, dict):
                raise RoleTruthError(f"{path}:{line_number}: case must be a JSON object")
            cases.append(case)
    return cases


def _text(value, field, where, *, empty=False):
    if not isinstance(value, str) or (not empty and not value.strip()):
        raise RoleTruthError(f"{where}: {field} must be {'a string' if empty else 'a non-empty string'}")
    if "\t" in value or "\n" in value or "\r" in value:
        raise RoleTruthError(f"{where}: {field} must not contain tabs or newlines")


def _go_primary_ids(path):
    """Return non-obsolete primary GO ids; alt_ids are deliberately excluded."""
    ids, current, obsolete, in_term = set(), None, False, False
    with Path(path).open() as fh:
        for raw in fh:
            line = raw.strip()
            if line == "[Term]":
                if in_term and current and not obsolete:
                    ids.add(current)
                current, obsolete, in_term = None, False, True
            elif line.startswith("["):
                if in_term and current and not obsolete:
                    ids.add(current)
                current, obsolete, in_term = None, False, False
            elif in_term and line.startswith("id: "):
                current = line[4:].strip()
            elif in_term and line == "is_obsolete: true":
                obsolete = True
    if in_term and current and not obsolete:
        ids.add(current)
    return ids


def verify_reference_files(case, references):
    """Verify every archived source declared by a case; return resolved paths."""
    paths = []
    for artifact in case.get("reference_files", []):
        relative = Path(artifact["path"])
        path = references / relative
        if not path.is_file():
            raise RoleTruthError(f"{case['case_id']}: reference artifact is missing: {path}")
        actual = sha256(path)
        if actual != artifact["sha256"]:
            raise RoleTruthError(f"{case['case_id']}: reference artifact SHA-256 differs: {artifact['path']}")
        paths.append(path)
    return paths


def validate_registry(registry, references=None, go_obo=None, *, require_positive=True):
    """Schema-check cases and optionally verify source hashes and GO identities.

    ``registry`` may be a JSONL path or an already loaded list. Negative-only
    observations use ``require_positive=False``; they otherwise follow the same
    schema. A successful return is a deterministic JSON-compatible summary; any
    violation raises RoleTruthError.
    """
    cases = load_registry(registry) if isinstance(registry, (str, Path)) else list(registry)
    references = Path(references) if references is not None else None
    primary_go = _go_primary_ids(go_obo) if go_obo is not None else None
    seen, multi = set(), []
    by_status, by_cohort, by_type = Counter(), Counter(), Counter()
    required = ("case_id", "status", "genome_id", "taxid", "organism", "in_training_cohort",
                "function_type", "function_id", "role", "members", "references", "source_statement",
                "reference_files", "universe_note", "added", "reviewer")
    for index, case in enumerate(cases, 1):
        if not isinstance(case, dict):
            raise RoleTruthError(f"case {index}: must be an object")
        where = f"case {index} ({case.get('case_id', 'no id')})"
        for field in required:
            if field not in case:
                raise RoleTruthError(f"{where}: missing field {field}")
        for field in ("case_id", "genome_id", "taxid", "organism", "function_type", "function_id",
                      "role", "source_statement", "universe_note", "added", "reviewer"):
            _text(case[field], field, where)
        if case["case_id"] in seen:
            raise RoleTruthError(f"{where}: duplicate case_id {case['case_id']!r}")
        seen.add(case["case_id"])
        if case["status"] not in STATUSES:
            raise RoleTruthError(f"{where}: invalid status {case['status']!r}")
        if not isinstance(case["in_training_cohort"], bool):
            raise RoleTruthError(f"{where}: in_training_cohort must be boolean")
        if case["function_type"] not in FUNCTION_TYPES:
            raise RoleTruthError(f"{where}: function_type must be GO or EC")
        pattern = GO_ID if case["function_type"] == "GO" else EC_ID
        if not pattern.fullmatch(case["function_id"]):
            raise RoleTruthError(f"{where}: invalid {case['function_type']} function_id")
        if primary_go is not None and case["function_type"] == "GO" and case["function_id"] not in primary_go:
            raise RoleTruthError(f"{where}: function_id is not a non-obsolete primary GO id")
        try:
            datetime.strptime(case["added"], "%Y-%m-%d")
        except ValueError as exc:
            raise RoleTruthError(f"{where}: added must be an ISO date") from exc
        if "allow_multi_positive" in case and not isinstance(case["allow_multi_positive"], bool):
            raise RoleTruthError(f"{where}: allow_multi_positive must be boolean")
        members = case["members"]
        if not isinstance(members, list) or len(members) < 2:
            raise RoleTruthError(f"{where}: at least two members are required")
        proteins, positives, negatives = set(), 0, 0
        for member_index, member in enumerate(members, 1):
            mwhere = f"{where} member {member_index}"
            if not isinstance(member, dict):
                raise RoleTruthError(f"{mwhere}: must be an object")
            for field in ("protein_id", "swissprot_accession", "label", "role_truth", "basis"):
                if field not in member:
                    raise RoleTruthError(f"{mwhere}: missing field {field}")
            for field in ("protein_id", "role_truth", "basis"):
                _text(member[field], field, mwhere)
            label = member["label"]
            if label is not None and (isinstance(label, bool) or label not in (0, 1)):
                raise RoleTruthError(f"{mwhere}: label must be integer 0, integer 1, or null")
            if label is None and member["role_truth"] != "unknown":
                raise RoleTruthError(f"{mwhere}: a null label requires role_truth 'unknown'")
            _text(member["swissprot_accession"], "swissprot_accession", mwhere, empty=label is None)
            if not member["swissprot_accession"] and label is not None:
                raise RoleTruthError(f"{mwhere}: swissprot_accession may be empty only for a null-label member")
            if "evidence_code" in member:
                _text(member["evidence_code"], "evidence_code", mwhere)
            if member["protein_id"] in proteins:
                raise RoleTruthError(f"{mwhere}: duplicate protein_id {member['protein_id']!r}")
            proteins.add(member["protein_id"])
            positives += label == 1
            negatives += label == 0
        if require_positive and positives != 1 and not (case.get("allow_multi_positive") and positives >= 1):
            raise RoleTruthError(f"{where}: exactly one label-1 member is required unless allow_multi_positive is set")
        if not require_positive and positives > 1 and not case.get("allow_multi_positive"):
            raise RoleTruthError(f"{where}: at most one label-1 member is allowed unless allow_multi_positive is set")
        if require_positive and positives == 0:
            raise RoleTruthError(f"{where}: at least one label-1 member is required")
        if negatives == 0:
            raise RoleTruthError(f"{where}: at least one label-0 member is required")
        if positives > 1:
            multi.append(case["case_id"])
        refs = case["references"]
        if not isinstance(refs, list) or not refs or any(not isinstance(r, str) or not r.strip() for r in refs):
            raise RoleTruthError(f"{where}: references must be a non-empty string list")
        if not any(r.startswith(("PMID:", "DOI:", "PMC")) for r in refs):
            raise RoleTruthError(f"{where}: at least one PMID:, DOI:, or PMC reference is required")
        artifacts = case["reference_files"]
        if not isinstance(artifacts, list):
            raise RoleTruthError(f"{where}: reference_files must be a list")
        for artifact in artifacts:
            if not isinstance(artifact, dict) or set(artifact) != {"path", "sha256"}:
                raise RoleTruthError(f"{where}: each reference_files entry needs exactly path and sha256")
            _text(artifact["path"], "reference_files.path", where)
            relative = Path(artifact["path"])
            if relative.is_absolute() or ".." in relative.parts:
                raise RoleTruthError(f"{where}: reference artifact paths must be relative and stay below --references")
            if not isinstance(artifact["sha256"], str) or not SHA256.fullmatch(artifact["sha256"]):
                raise RoleTruthError(f"{where}: reference_files.sha256 must be 64 lowercase hexadecimal characters")
        if "gaf_assertion" in case:
            gaf = case["gaf_assertion"]
            required_gaf = {"inventory", "accession", "relation", "evidence", "reference", "candidate_status"}
            if not isinstance(gaf, dict) or not required_gaf <= set(gaf):
                raise RoleTruthError(f"{where}: incomplete gaf_assertion")
            for field in required_gaf:
                _text(gaf[field], f"gaf_assertion.{field}", where)
        if references is not None:
            verify_reference_files(case, references)
        by_status[case["status"]] += 1
        by_cohort[str(case["in_training_cohort"]).lower()] += 1
        by_type[case["function_type"]] += 1
    return {
        "schema": "gspa.role-truth-registry-validation/1", "status": "valid", "cases": len(cases),
        "cases_by_status": {key: by_status[key] for key in STATUSES if by_status[key]},
        "cases_by_in_training_cohort": {key: by_cohort[key] for key in ("true", "false") if by_cohort[key]},
        "cases_by_function_type": {key: by_type[key] for key in FUNCTION_TYPES if by_type[key]},
        "allow_multi_positive_cases": multi,
    }


def verify_entry(entry_path, accession, protein, taxid, proteins, snapshot_date=None):
    """Check one UniProt record against an acquired protein, without a corpus scan."""
    entries = list(M.iter_entries(entry_path))
    if len(entries) != 1:
        raise RoleTruthError("expected one UniProt record")
    entry = entries[0]
    if entry["accessions"][0] != accession or entry["taxid"] != str(taxid):
        raise RoleTruthError("UniProt identity/taxon differs")
    if protein not in entry["xrefs"] or proteins.get(protein) != entry["sequence"] or not entry["sequence"]:
        raise RoleTruthError("full sequence and versioned accession must both match")
    dates = [line for line in entry_path.read_text().splitlines()
             if line.startswith("DT   ") and "sequence version" in line]
    if len(dates) != 1:
        raise RoleTruthError("expected one UniProt sequence-version date")
    if snapshot_date is not None and datetime.strptime(dates[0][5:16], "%d-%b-%Y").date() > snapshot_date:
        raise RoleTruthError("sequence version postdates the annotation snapshot")
    basis = "versioned_accession_and_full_sequence"
    if snapshot_date is not None:
        basis += "; last sequence update predates GAF"
    return {"sequence_check": "verified", "accession": accession, "protein_id": protein,
            "taxid": str(taxid), "sequence_sha256": hashlib.sha256(entry["sequence"].encode()).hexdigest(),
            "length": len(entry["sequence"]), "sequence_version": dates[0], "basis": basis}


def _gaf_cross_check(case, references):
    """Reproduce the predecessor's inventory check only for a GAF-derived case."""
    gaf = case["gaf_assertion"]
    inventory_path = Path(gaf["inventory"])
    if not inventory_path.is_absolute():
        inventory_path = references / inventory_path
    inventory_path = inventory_path.resolve()
    inventory = _read_json(inventory_path)
    if inventory.get("schema") != "gspa.gaf-negative-inventory/1" or inventory.get("status") != "verified":
        raise RoleTruthError(f"{case['case_id']}: complete GAF inventory required")
    candidates = inventory_path.parent / "cohort-candidates.jsonl"
    if sha256(candidates) != inventory["output_sha256"]["cohort-candidates.jsonl"]:
        raise RoleTruthError(f"{case['case_id']}: GAF cohort records changed")
    negative = next(m for m in case["members"] if m["swissprot_accession"] == gaf["accession"])
    found = [r for r in _read_jsonl(candidates)
             if r["accession"] == gaf["accession"] and r["go_id"] == case["function_id"]
             and r["relation"] == gaf["relation"] and r["evidence"] == gaf["evidence"]
             and gaf["reference"] in r["references"]]
    if len(found) != 1 or not any(j["genome_id"] == case["genome_id"] and j["protein_id"] == negative["protein_id"]
                                and j["status"] == gaf["candidate_status"] for j in found[0]["candidate_joins"]):
        raise RoleTruthError(f"{case['case_id']}: expected unqualified exact GAF assertion absent")
    snapshots = [s for s in inventory["headers"] if s.startswith("!date-generated: ")]
    if len(snapshots) != 1:
        raise RoleTruthError(f"{case['case_id']}: GAF inventory needs one snapshot date")
    snapshot = datetime.strptime(snapshots[0][17:27], "%Y-%m-%d").date()
    return inventory_path, candidates, snapshot


def _label_function_id(case, value):
    value = value.strip().replace("_", ":")
    if case["function_type"] == "EC" and not value.upper().startswith("EC:"):
        value = "EC:" + value
    return value


def _primary_reference(case):
    return next((r for r in case["references"] if r.startswith(("PMID:", "DOI:", "PMC"))), case["references"][0])


def _acquired_proteins(case, provenance):
    organism = provenance.get("organisms", {}).get(str(case["taxid"]))
    if organism is None:
        raise RoleTruthError(f"{case['case_id']}: taxon absent from acquisition provenance")
    assembly = organism["assembly"]
    if assembly["used_accession"] != case["genome_id"]:
        raise RoleTruthError(f"{case['case_id']}: assembly identity differs")
    fasta = Path(assembly["directory"]) / "protein.faa"
    if sha256(fasta) != assembly["member_sha256"]["protein.faa"]:
        raise RoleTruthError(f"{case['case_id']}: acquired genome proteins changed")
    return fasta, M.PL.fasta_sequences(fasta)


def _check_original_rows(before, after, expected_new):
    """Original examples are immutable; only declared positive/negative ids may be new."""
    keys = ("example_id", "features", "label", "split", "group")
    before_ids = {r["example_id"] for r in before}
    after_by_id = {r["example_id"]: r for r in after}
    if len(after_by_id) != len(after):
        raise RoleTruthError("augmented table has duplicate example identities")
    for old in before:
        new = after_by_id.get(old["example_id"])
        if new is None or any(old[k] != new[k] for k in keys):
            raise RoleTruthError("candidate features, labels, identities or splits changed")
    for new in after:
        if new["example_id"] in before_ids:
            continue
        if new["example_id"] not in expected_new or new["label"] != expected_new[new["example_id"]]:
            raise RoleTruthError("an undeclared training example was added")
    return True


def build(registry, original_manifest, provenance_path, references, out, statuses=("accepted",)):
    """Build complete reviewed loci into a copy of an original training cohort."""
    registry = Path(registry).resolve()
    original_manifest = Path(original_manifest).resolve()
    provenance_path = Path(provenance_path).resolve()
    references, out = Path(references).resolve(), Path(out).resolve()
    all_cases = load_registry(registry)
    validate_registry(all_cases, references=references)
    selected = [c for c in all_cases if c["in_training_cohort"] and c["status"] in set(statuses)]
    if not selected:
        raise RoleTruthError("no registry cases selected")
    old = _read_json(original_manifest)
    inputs = old["inputs"]
    provenance = _read_json(provenance_path)
    original_table = Path(old["table"])
    if sha256(original_table) != old["table_sha256"]:
        raise RoleTruthError("original training table changed")
    for key in ("contract", "labels", "splits", "label_policy", "registry", "go_obo", "example_clusters"):
        if inputs.get(key):
            path = Path(inputs[key])
            if inputs.get(key + "_sha256") and sha256(path) != inputs[key + "_sha256"]:
                raise RoleTruthError("original input hash differs: " + key)
    with Path(inputs["labels"]).open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        original_fields = list(reader.fieldnames or [])
        labels = list(reader)
    if not labels:
        raise RoleTruthError("this augmentation requires original positive-only labels")
    with Path(inputs["splits"]).open(newline="") as fh:
        split_rows = list(csv.DictReader(fh, delimiter="\t"))
    split_by_genome = {r["genome_id"]: r["split"] for r in split_rows}
    source_by_genome = {r["genome_id"]: r["source"] for r in split_rows}

    # Resolve all declared truth before writing anything. A reviewed positive may
    # be absent from Swiss-Prot-derived labels; only that declared gap is filled.
    # Negatives copy the active row's evidence and other source fields.
    active_rows, expected_new = {}, {}
    positive_additions, declared_positive_additions = [], []
    for case in selected:
        if case["genome_id"] not in source_by_genome:
            raise RoleTruthError(f"{case['case_id']}: genome absent from original splits")
        for member in case["members"]:
            if member["label"] is None:
                continue
            matches = [r for r in labels if r["genome_id"] == case["genome_id"]
                       and r["protein_id"] == member["protein_id"]
                       and _label_function_id(case, r["function_id"]) == case["function_id"]]
            original_positives = [r for r in matches if r.get("label", "1") in ("", "1")
                                  and not r.get("qualifier", "")]
            original_negatives = [r for r in matches if r.get("label") == "0" or r.get("qualifier") == "NOT"]
            example_id = (f"{case['genome_id']}|{member['protein_id']}|{case['function_type']}|"
                          f"{case['function_id']}|{case['role']}")
            if member["label"] == 1:
                if original_negatives:
                    raise RoleTruthError(f"{case['case_id']}: conflicting original truth (label-0 row) for "
                                         f"{member['protein_id']}")
                if len(original_positives) > 1:
                    raise RoleTruthError(f"{case['case_id']}: expected at most one original positive for "
                                         f"{member['protein_id']}")
                if original_positives:
                    active_rows[(case["case_id"], member["protein_id"])] = original_positives[0]
                else:
                    evidence = member.get("evidence_code", "IDA")
                    reference = _primary_reference(case)
                    row = {field: "" for field in original_fields}
                    row.update({"genome_id": case["genome_id"], "protein_id": member["protein_id"],
                                "function_id": case["function_id"], "evidence": evidence, "tier": "experimental",
                                "reference": reference, "swissprot_accession": member["swissprot_accession"],
                                "label": "1"})
                    if "function_type" in row:
                        row["function_type"] = case["function_type"]
                    if "role" in row:
                        row["role"] = case["role"]
                    active_rows[(case["case_id"], member["protein_id"])] = row
                    positive_additions.append(row)
                    expected_new[example_id] = 1
                    declared_positive_additions.append({
                        "case_id": case["case_id"], "protein_id": member["protein_id"],
                        "swissprot_accession": member["swissprot_accession"],
                        "function_type": case["function_type"], "function_id": case["function_id"],
                        "role": case["role"], "example_id": example_id, "evidence": evidence,
                        "tier": "experimental", "reference": reference,
                    })
            else:
                if original_positives:
                    raise RoleTruthError(f"{case['case_id']}: negative member already has an original positive for "
                                         f"{member['protein_id']}")
                if matches:
                    raise RoleTruthError(f"{case['case_id']}: negative member already has original truth")
                expected_new[example_id] = 0
    if any(r.get("label", "1") != "1" or r.get("qualifier", "") for r in labels):
        raise RoleTruthError("this augmentation requires original positive-only labels")
    for row in labels:
        row["label"] = "1"
        row["reference"] = ""
    labels.extend(positive_additions)
    for case in selected:
        reference = _primary_reference(case)
        for member in case["members"]:
            if member["label"] == 1:
                active_rows[(case["case_id"], member["protein_id"])]["reference"] = reference
        for member in case["members"]:
            if member["label"] == 0:
                active = next(m for m in case["members"] if m["label"] == 1)
                row = dict(active_rows[(case["case_id"], active["protein_id"])],
                           protein_id=member["protein_id"], swissprot_accession=member["swissprot_accession"],
                           label="0", reference=reference)
                labels.append(row)
    fields = list(original_fields)
    for field in ("label", "reference", "swissprot_accession"):
        if field not in fields:
            fields.append(field)
    if declared_positive_additions and "tier" not in fields:
        fields.append("tier")

    policy = _read_json(Path(inputs["label_policy"]))
    per_policy = {}
    certifiable_cases = [c for c in selected if not any(m["label"] is None for m in c["members"])]
    for case in selected:
        source = source_by_genome[case["genome_id"]]
        key = (source, case["function_type"])
        per_policy.setdefault(key, []).append(case)
    for (source, kind), cases in per_policy.items():
        try:
            target = policy["datasets"][source][kind]
        except KeyError as exc:
            raise RoleTruthError(f"label policy has no dataset/function entry for {source}/{kind}") from exc
        target["explicit_negative_truth"] = (cases[0]["source_statement"] if len(cases) == 1 else
            "\n".join(f"{c['case_id']}: {c['source_statement']}" for c in cases))
    # Newlines in a multi-case policy justification are intentional JSON content;
    # the registry itself restricts source statements to one line each.
    policy["locus_truth_complete"] = (certifiable_cases[0]["source_statement"] if len(certifiable_cases) == 1 else
        "\n".join(f"{c['case_id']}: {c['source_statement']}" for c in certifiable_cases))

    # Verify source artifacts, acquired sequences, and optional GAF provenance.
    files = [registry, original_manifest, provenance_path, original_table, Path(__file__),
             Path(A.__file__ or ""), Path(M.__file__ or ""), Path(M.PL.__file__ or ""),
             Path(P.__file__ or ""), Path(T.__file__ or ""), Path(M.EV.__file__ or "")]
    for key in ("contract", "labels", "splits", "label_policy", "registry", "go_obo", "example_clusters"):
        if inputs.get(key):
            files.append(Path(inputs[key]))
    files.extend(p for p in Path(inputs["feature_dir"]).rglob("*") if p.is_file())
    sequence_checks = []
    protein_cache = {}
    for case in selected:
        files.extend(verify_reference_files(case, references))
        if case["genome_id"] not in protein_cache:
            protein_cache[case["genome_id"]] = _acquired_proteins(case, provenance)
        fasta, proteins = protein_cache[case["genome_id"]]
        files.append(fasta)
        snapshot = None
        if "gaf_assertion" in case:
            inventory, candidates, snapshot = _gaf_cross_check(case, references)
            files.extend((inventory, candidates))
        for member in case["members"]:
            accession = member["swissprot_accession"]
            entry = references / (accession + ".txt")
            if not accession or not entry.is_file():
                sequence_checks.append({"case_id": case["case_id"], "protein_id": member["protein_id"],
                                        "accession": accession, "sequence_check": "skipped",
                                        "reason": ("no Swiss-Prot accession declared" if not accession else
                                                   "UniProt .txt record not present"),
                                        "expected_path": str(entry) if accession else None})
                continue
            files.append(entry)
            check = verify_entry(entry, member["swissprot_accession"], member["protein_id"], case["taxid"],
                                 proteins, snapshot)
            check["case_id"] = case["case_id"]
            sequence_checks.append(check)
    identity = {str(p.resolve()): sha256(p) for p in files}

    out.mkdir(parents=True, exist_ok=False)
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(labels)
    (out / "labels.tsv").write_text(stream.getvalue())
    (out / "label_policy.json").write_text(json.dumps(policy, indent=2) + "\n")
    locus = out / "locus_groups.tsv"
    with locus.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(("genome_id", "locus_group_id", "candidate_key", "truth_complete"))
        for case in certifiable_cases:
            for member in case["members"]:
                candidate = (f"{member['protein_id']}|{case['function_type']}|{case['function_id']}|{case['role']}")
                writer.writerow((case["genome_id"], case["case_id"], candidate, "true"))
    result = P.build_table(Path(inputs["feature_dir"]), Path(inputs["contract"]), out / "labels.tsv",
        Path(inputs["splits"]), out / "label_policy.json", out / "table",
        registry_path=Path(inputs["registry"]), go_obo=inputs.get("go_obo"),
        example_clusters=inputs.get("example_clusters"), locus_groups=locus,
        disjoint=tuple(old["split_audit"]["disjoint_keys"]))
    before = _read_jsonl(original_table)
    after = _read_jsonl(result["table"])
    _check_original_rows(before, after, expected_new)
    contract = T.load_contract(Path(inputs["contract"]))
    rows = T.read_table(Path(result["table"]), contract)
    T.validate_locus_truth(rows)

    audits = {g["locus_group_id"]: g for g in result["locus_audit"]["groups"]}
    after_ids = {r["example_id"] for r in after}
    for case in selected:
        if not any(m["label"] is None for m in case["members"]):
            continue
        declared = [f"{m['protein_id']}|{case['function_type']}|{case['function_id']}|{case['role']}"
                    for m in case["members"]]
        retained = [candidate for m, candidate in zip(case["members"], declared)
                    if f"{case['genome_id']}|{candidate}" in after_ids]
        audits[case["case_id"]] = {
            "genome_id": case["genome_id"], "locus_group_id": case["case_id"], "eligible": False,
            "reasons": ["candidate_universe_incomplete_by_declaration"],
            "declared_members": sorted(declared), "retained_members": sorted(retained),
        }
    locus_audit = dict(result["locus_audit"])
    locus_audit["groups"] = sorted(audits.values(), key=lambda g: (g["genome_id"], g["locus_group_id"]))
    locus_audit["declared_groups"] = len(audits)
    locus_audit["eligible_groups"] = sum(g["eligible"] for g in audits.values())
    case_reports = []
    complete: dict[str, int] = dict.fromkeys(SPLITS, 0)
    for case in selected:
        audit = audits.get(case["case_id"])
        if audit is None:
            raise RoleTruthError(f"{case['case_id']}: case silently missing from locus certification")
        split = split_by_genome[case["genome_id"]]
        if audit["eligible"]:
            complete[split] += 1
        case_reports.append({"case_id": case["case_id"], "split": split, **audit})
    for path, digest in identity.items():
        if sha256(path) != digest:
            raise RoleTruthError("source changed: " + path)
    report = {
        "schema": SCHEMA, "status": "verified", "registry_sha256": sha256(registry),
        "selected_case_ids": [c["case_id"] for c in selected], "cases": case_reports,
        "sources": {c["case_id"]: c["source_statement"] for c in selected},
        "sequence_checks": sequence_checks,
        "skipped_sequence_checks": [c for c in sequence_checks if c["sequence_check"] == "skipped"],
        "declared_positive_additions": declared_positive_additions,
        "identity_sha256": identity, "rows": len(rows), "original_features_labels_splits_unchanged": True,
        "locus_audit": locus_audit,
        "explicit_negative_rows": sum(r.get("label_basis") == "explicit_negative_truth" for r in rows),
        "complete_groups_by_split": complete,
        "scope": "Reviewed complete training loci only; not a sufficient null-model training set or independent evaluation.",
        "output_sha256": {str(p.relative_to(out)): sha256(p) for p in out.rglob("*") if p.is_file()},
    }
    (out / "verification.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


# Descriptive alias for callers that prefer the subcommand name.
build_loci = build


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    validate = sub.add_parser("validate")
    validate.add_argument("--registry", type=Path, required=True)
    validate.add_argument("--references", type=Path)
    validate.add_argument("--go-obo", type=Path)
    validate.add_argument("--negative-only", action="store_true",
                          help="validate negative-only observations (no positive member required)")
    builder = sub.add_parser("build-loci")
    for name in ("registry", "original-manifest", "provenance", "references", "out"):
        builder.add_argument("--" + name, type=Path, required=True)
    builder.add_argument("--statuses", default="accepted")
    args = parser.parse_args(argv)
    try:
        if args.command == "validate":
            report = validate_registry(args.registry, references=args.references, go_obo=args.go_obo,
                                       require_positive=not args.negative_only)
        else:
            statuses = tuple(s.strip() for s in args.statuses.split(",") if s.strip())
            unknown = set(statuses) - set(STATUSES)
            if unknown:
                raise RoleTruthError(f"unknown statuses: {sorted(unknown)}")
            report = build(args.registry, args.original_manifest, args.provenance, args.references, args.out, statuses)
            report = {k: report[k] for k in ("status", "rows", "explicit_negative_rows", "complete_groups_by_split")}
    except (OSError, KeyError, json.JSONDecodeError, RoleTruthError, P.PreparationError, T.ProvenanceError) as exc:
        parser.exit(2, f"error: {exc}\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
