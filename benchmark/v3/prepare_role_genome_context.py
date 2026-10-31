#!/usr/bin/env python3
"""Acquire role-genome pathway observations and prepare joint 21-genome arms.

The three role genomes are integrated jointly with the 18 source cohort
members.  Gap records are a stdlib-only replay of PathwayGapBuilder's observed
-table rule; cohort replay parity is mandatory before role gaps are emitted.
No role truth or target functional annotation is read by this program.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shlex
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prepare_context_modes as CONTEXT  # noqa: E402
import replay_cohort_pathways as REPLAY  # noqa: E402
import run_source_arm_context_comparison as SOURCE  # noqa: E402

DATA = Path.home() / ".cache" / "gspa-v3-data"
DEFAULT_COHORT_ROOTS = (
    DATA / "cohort-pathway-search-v1",
    DATA / "cohort-pathway-search-archaea-v2",
)
DEFAULT_COHORT_REPLAY = DATA / "cohort-pathway-replay-v2"
CONTEXTS = ("SEQONLY", "UNIVERSAL", "PATHWAY", "ORGANISM")
KINDS = ("none", "logistic", "gbdt")
OLD_CLAIMS_PREFIX = "/home/leechuck/tmp/source-context-results/training/"
NEW_CLAIMS_PREFIX = "/data/hohndor/gspa-v3/source-context-full-v1/jobs/10406/checkpoint/training/"
ACQUISITION_SCHEMA = "gspa.role-genome-pathway-acquisition/1"
CONTEXT_SCHEMA = "gspa.role-genome-context-inputs/1"
LAYOUT_FIELDS = ["protein_id", "contig", "start", "end", "strand"]


class RoleContextError(ValueError):
    """An input, parity, or plan contract was not satisfied."""


class PendingGenomes(RoleContextError):
    """One or more role-genome searches have not completed."""


def require(condition, message):
    if not condition:
        raise RoleContextError(message)


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            value.update(block)
    return value.hexdigest()


def save_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def read_json(path):
    return json.loads(Path(path).read_text())


def read_tsv(path):
    return REPLAY.rows(Path(path))


def write_tsv(path, fields, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def file_hashes(paths):
    return {str(Path(path).resolve()): digest(path) for path in paths}


def role_rows(path):
    rows = read_tsv(path)
    required = {"genome_id", "taxid", "proteome_sha256", "claims_path"}
    require(rows and required <= set(rows[0]), "role_genomes.tsv lacks required columns")
    ids = [row["genome_id"] for row in rows]
    require(len(ids) == len(set(ids)), "duplicate genome in role_genomes.tsv")
    return rows


def role_layout_rows(path, proteomes=False):
    rows = role_rows(path)
    required = {"key", "proteome"} if proteomes else {"key"}
    for row in rows:
        missing = sorted(column for column in required if not row.get(column))
        require(not missing,
                f"role_genomes.tsv lacks {', '.join(missing)} for {row['genome_id']}")
    return rows


def fasta_identifiers(path, expected_sha256, genome_id):
    path = Path(path).resolve()
    require(path.is_file(), f"missing proteome FASTA for {genome_id}: {path}")
    actual_sha256 = digest(path)
    require(actual_sha256 == expected_sha256,
            f"proteome sha256 differs for {genome_id}: {path}")
    identifiers = set()
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.startswith(">"):
                continue
            fields = line[1:].split()
            require(fields, f"empty FASTA identifier for {genome_id} at line {line_number}")
            protein_id = fields[0]
            require(protein_id not in identifiers,
                    f"duplicate FASTA protein id for {genome_id}: {protein_id}")
            identifiers.add(protein_id)
    return path, identifiers


def layout_identifiers(path, fasta_ids, contigs, genome_id):
    path = Path(path).resolve()
    require(path.is_file(), f"missing layout for {genome_id}: {path}")
    require(isinstance(contigs, dict), f"metadata contigs missing for {genome_id}")
    identifiers = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames == LAYOUT_FIELDS,
                f"layout header must be {' '.join(LAYOUT_FIELDS)} for {genome_id}: {path}")
        for line_number, row in enumerate(reader, 2):
            require(None not in row and all(row.get(field) is not None for field in LAYOUT_FIELDS),
                    f"malformed layout row for {genome_id} at line {line_number}")
            protein_id = row["protein_id"]
            require(protein_id, f"empty layout protein id for {genome_id} at line {line_number}")
            require(protein_id not in identifiers,
                    f"duplicate layout protein id for {genome_id}: {protein_id}")
            require(protein_id in fasta_ids,
                    f"layout protein id absent from proteome FASTA for {genome_id}: {protein_id}")
            require(row["contig"] in contigs,
                    f"layout contig absent from metadata for {genome_id}: {row['contig']}")
            identifiers.add(protein_id)
    return identifiers


def safe_output(directory, relative):
    path = Path(relative)
    require(not path.is_absolute() and ".." not in path.parts,
            "observation output path is not relative to its genome directory: " + relative)
    return directory / path


def acquire(search_root, role_genomes, cohort_reference_hashes, out):
    """Verify role observations, writing only a pending or fully verified report."""
    search_root = Path(search_root).resolve()
    cohort_reference_hashes = Path(cohort_reference_hashes).resolve()
    out = Path(out)
    records = {}
    pending = []
    for role in role_rows(role_genomes):
        gid = role["genome_id"]
        directory = search_root / "observations" / gid
        observation_path = directory / "observation.json"
        observation = read_json(observation_path) if observation_path.exists() else {}
        if observation.get("status") != "observed":
            pending.append(gid)
            records[gid] = {
                "status": "pending",
                "search_root": str(search_root),
                "observation_path": str(observation_path),
                "observed_status": observation.get("status", "missing"),
            }
            continue
        require(observation.get("proteome_sha256") == role["proteome_sha256"],
                f"proteome sha256 differs for {gid}")
        reference_path = directory / "reference-hashes.json"
        require(reference_path.is_file(), f"missing reference-hashes.json for {gid}")
        require(reference_path.read_bytes() == cohort_reference_hashes.read_bytes(),
                f"reference hashes differ from cohort for {gid}")
        output_hashes = observation.get("output_hashes")
        require(isinstance(output_hashes, dict) and output_hashes,
                f"observation has no output hashes for {gid}")
        verified_outputs = {}
        for relative, expected in sorted(output_hashes.items()):
            path = safe_output(directory, relative)
            require(path.is_file(), f"missing observed output for {gid}: {relative}")
            actual = digest(path)
            require(actual == expected, f"output sha256 differs for {gid}: {relative}")
            verified_outputs[relative] = actual
        reactions = directory / "outputs" / f"{gid}-all-Reactions.tbl"
        pathways = directory / "outputs" / f"{gid}-all-Pathways.tbl"
        require(reactions.is_file() and pathways.is_file(),
                f"reaction/pathway tables are missing for {gid}")
        reaction_rows = read_tsv(reactions)
        pathway_rows = read_tsv(pathways)
        records[gid] = {
            "status": "observed",
            "search_root": str(search_root),
            "observation_path": str(observation_path),
            "observation_sha256": digest(observation_path),
            "reference_hashes_path": str(reference_path),
            "reference_hashes_sha256": digest(reference_path),
            "proteome_sha256": observation["proteome_sha256"],
            "output_hashes": verified_outputs,
            "reactions_path": str(reactions),
            "pathways_path": str(pathways),
            "reaction_rows": len(reaction_rows),
            "pathway_rows": len(pathway_rows),
            "predicted_pathways": sum(row["Prediction"].strip().upper() == "TRUE"
                                      for row in pathway_rows),
        }
    report = {
        "schema": ACQUISITION_SCHEMA,
        "status": "pending" if pending else "verified",
        "genomes": records,
        "pending_genomes": pending,
        "inputs": file_hashes((role_genomes, cohort_reference_hashes)),
        "scope": "Observed search integrity only; gap records are hypotheses, not evidence of biological absence.",
    }
    save_json(out / "acquisition.json", report)
    if pending:
        raise PendingGenomes("role genome searches pending: " + ", ".join(pending))
    return report


def reaction_flag(value):
    require(isinstance(value, str) and value.strip().upper() in ("TRUE", "FALSE"),
            "invalid pathway/reaction boolean: " + repr(value))
    return value.strip().upper() == "TRUE"


def identifier(value, label):
    require(isinstance(value, str), f"missing {label}")
    result = value.strip().strip("|")
    require(bool(result), f"empty or invalid {label}")
    return result


def ec_numbers(value):
    """PathwayGapBuilder.ecNumbers, including preliminary/partial EC rules."""
    require(isinstance(value, str), "missing reaction EC field")
    numbers = set()
    for token in re.split(r"[/,;\s]+", value.strip()):
        ec = re.sub(r"^EC:", "", token, count=1)
        if not ec or ec in ("NA", "-"):
            continue
        if re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+\.(?!n[0-9]+$)[A-Za-z][A-Za-z0-9]*", ec):
            continue
        if not (re.fullmatch(r"[0-9]+(?:\.[0-9]+){0,3}", ec) or
                re.fullmatch(r"[0-9]+\.(?:[0-9]+|-)\.(?:[0-9]+|-)\.(?:[0-9]+|n[0-9]+|-)", ec)):
            raise RoleContextError("invalid reaction EC: " + token)
        numbers.add(ec)
    return numbers


def load_ec2go(path):
    """PathwayLoader.loadEc2GoAlternatives, retaining every alternative."""
    result = defaultdict(list)
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("!") or line.startswith("#"):
                continue
            standard = re.search(r"^(EC:\S+)\s*>\s*GO:[^;]+;\s*(GO:\d{7})", line)
            if standard:
                ec, go = standard.group(1)[3:], standard.group(2)
            else:
                fields = line.strip().split("\t")
                if len(fields) < 2 or not re.fullmatch(r"GO:\d{7}", fields[1].strip()):
                    continue
                ec = re.sub(r"^EC:", "", fields[0].strip(), count=1)
                go = fields[1].strip()
            if go not in result[ec]:
                result[ec].append(go)
    return {ec: sorted(terms) for ec, terms in sorted(result.items())}


def pathway_predictions(rows):
    predictions = {}
    for row in rows:
        pid = identifier(row.get("ID"), "pathway ID")
        predicted = reaction_flag(row.get("Prediction"))
        require(pid not in predictions or predictions[pid] == predicted,
                "conflicting pathway predictions: " + pid)
        predictions[pid] = predicted
    return predictions


def build_gap_records(reactions, pathways, ec2go):
    """Pure-Python equivalent of PathwayGapBuilder.build for gapsmith/proteome."""
    predictions = pathway_predictions(pathways)
    grouped = defaultdict(list)
    for row in reactions:
        pathway = identifier(row.get("pathway"), "reaction pathway ID")
        reaction = identifier(row.get("rxn"), "reaction ID")
        require(row.get("status") in ("good_blast", "bad_blast", "no_blast",
                                     "no_seq_data", "spontaneous"),
                "unknown reaction status: " + repr(row.get("status")))
        reaction_flag(row.get("spont"))
        reaction_flag(row.get("is_complex"))
        ec_numbers(row.get("ec"))
        grouped[(pathway, reaction)].append(row)
    records = []
    for pathway, reaction in sorted(grouped):
        rows = grouped[(pathway, reaction)]
        excluded = (predictions.get(pathway) is not True or
                    any(row["status"] == "good_blast" for row in rows) or
                    any(row["status"] == "no_seq_data" for row in rows) or
                    any(reaction_flag(row["spont"]) or row["status"] == "spontaneous"
                        for row in rows) or
                    any(reaction_flag(row["is_complex"]) for row in rows))
        if excluded:
            continue
        ecs = sorted(set().union(*(ec_numbers(row["ec"]) for row in rows)))
        guessed = any(row["status"] == "bad_blast" for row in rows)
        for ec in ecs or [None]:
            gos = sorted(set(ec2go.get(ec, []))) or [None]
            for go in gos:
                records.append({
                    "pathway_id": "gapsmith:" + pathway,
                    "reaction_id": reaction,
                    "ec_number": ec,
                    "go_term": go,
                    "gapseq_guessed": guessed,
                })
    return records


def derive_gap_records(reactions, pathways, ec2go):
    return build_gap_records(read_tsv(reactions), read_tsv(pathways), load_ec2go(ec2go))


def write_jsonl(path, records):
    with Path(path).open("w") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")


def cohort_parity(roots, replay_root, ec2go):
    """Compare canonical derived and JVM replay records for every available genome."""
    mapping = load_ec2go(ec2go)
    reports = {}
    for root in roots:
        root = Path(root).resolve()
        observations = root / "observations"
        if not observations.is_dir():
            continue
        for directory in sorted(observations.iterdir()):
            if not directory.is_dir():
                continue
            gid = directory.name
            replay_file = (Path(replay_root) / gid / "annotation" /
                           f"{gid}_integrated.tsv.gaps.jsonl")
            observation_path = directory / "observation.json"
            if not replay_file.is_file() or not observation_path.is_file():
                continue
            if read_json(observation_path).get("status") != "observed":
                continue
            require(gid not in reports, "multiple observed cohort searches for " + gid)
            reactions = directory / "outputs" / f"{gid}-all-Reactions.tbl"
            pathways = directory / "outputs" / f"{gid}-all-Pathways.tbl"
            derived = build_gap_records(read_tsv(reactions), read_tsv(pathways), mapping)
            replayed = REPLAY.jsonl(replay_file)
            derived_canonical = set(REPLAY.canonical(derived))
            replay_canonical = set(REPLAY.canonical(replayed))
            missing = sorted(replay_canonical - derived_canonical)
            unexpected = sorted(derived_canonical - replay_canonical)
            reports[gid] = {
                "search_root": str(root),
                "equal": not missing and not unexpected,
                "derived_records": len(derived),
                "replay_records": len(replayed),
                "equal_records": len(derived_canonical & replay_canonical),
                "unequal_records": len(missing) + len(unexpected),
                "missing_records": [json.loads(value) for value in missing],
                "unexpected_records": [json.loads(value) for value in unexpected],
            }
    require(reports, "no cohort replay gap files found for parity")
    failed = sorted(gid for gid, report in reports.items() if not report["equal"])
    return {
        "schema": "gspa.role-genome-gap-parity/1",
        "status": "failed" if failed else "verified",
        "genomes": reports,
        "genomes_compared": len(reports),
        "genomes_equal": len(reports) - len(failed),
        "genomes_unequal": len(failed),
        "failed_genomes": failed,
    }


def prepare_gaps(acquisition, ec2go, out, cohort_roots=DEFAULT_COHORT_ROOTS,
                 cohort_replay=DEFAULT_COHORT_REPLAY, parity_only=False):
    out = Path(out)
    parity = cohort_parity(cohort_roots, cohort_replay, ec2go)
    save_json(out / "parity.json", parity)
    if parity["status"] != "verified":
        details = "; ".join(
            f"{gid}: missing={len(parity['genomes'][gid]['missing_records'])}, "
            f"unexpected={len(parity['genomes'][gid]['unexpected_records'])}"
            for gid in parity["failed_genomes"])
        raise RoleContextError("cohort gap parity failed; see parity.json: " + details)
    if parity_only:
        return parity
    acquired = read_json(acquisition)
    require(acquired.get("schema") == ACQUISITION_SCHEMA and acquired.get("status") == "verified",
            "verified role-genome acquisition required")
    counts = {}
    outputs = {}
    for gid, observed in sorted(acquired["genomes"].items()):
        require(observed.get("status") == "observed", "role genome is not observed: " + gid)
        reactions, pathways = Path(observed["reactions_path"]), Path(observed["pathways_path"])
        raw, pwy = read_tsv(reactions), read_tsv(pathways)
        gaps = build_gap_records(raw, pwy, load_ec2go(ec2go))
        REPLAY.verify_gap_groups(gaps, REPLAY.eligible_groups(raw, pwy))
        target = out / f"{gid}.gaps.jsonl"
        write_jsonl(target, gaps)
        outputs[target.name] = digest(target)
        counts[gid] = {
            "gap_records": len(gaps),
            "records_with_go": sum(record["go_term"] is not None for record in gaps),
            "gapseq_guessed_records": sum(record["gapseq_guessed"] for record in gaps),
            "eligible_reaction_groups": len(REPLAY.eligible_groups(raw, pwy)),
            "reaction_rows": len(raw),
            "pathway_rows": len(pwy),
        }
    receipt = {
        "schema": "gspa.role-genome-gaps/1",
        "status": "verified",
        "counts": counts,
        "parity": {key: parity[key] for key in ("genomes_compared", "genomes_equal", "genomes_unequal")},
        "inputs": file_hashes((acquisition, ec2go)),
        "outputs": outputs,
    }
    save_json(out / "receipt.json", receipt)
    return receipt


def _group_keys(rows):
    return {tuple(sorted(row.items())) for row in rows}


def check_context_parity(old_metadata, new_metadata, old_groups, new_groups):
    """Field-exact cohort metadata equality and inclusion of every old group row."""
    old_genomes = old_metadata["genomes"]
    new_genomes = new_metadata["genomes"]
    changed = sorted(gid for gid in old_genomes if new_genomes.get(gid) != old_genomes[gid])
    missing = sorted(_group_keys(old_groups) - _group_keys(new_groups))
    return {
        "status": "failed" if changed or missing else "verified",
        "cohort_genomes": len(old_genomes),
        "changed_cohort_records": changed,
        "old_pathway_group_rows": len(old_groups),
        "new_pathway_group_rows": len(new_groups),
        "missing_pathway_group_rows": [dict(items) for items in missing],
        "new_genomes": sorted(set(new_genomes) - set(old_genomes)),
    }


def context_source_inputs(cohort_inputs):
    receipt = read_json(Path(cohort_inputs) / "receipt.json")
    require(receipt.get("status") == "verified", "verified cohort context inputs required")
    graphs = sorted(path for path in receipt["inputs"] if path.endswith("/seed_reactions.tsv"))
    manifests = sorted(path for path in receipt["inputs"] if path.endswith("/manifest.json"))
    require(len(graphs) == 1, "cohort receipt does not identify exactly one reaction graph")
    require(len(manifests) == 2, "cohort receipt does not identify exactly two search roots")
    return Path(graphs[0]), [Path(path).parent for path in manifests], receipt


def prepare_disjoint_roots(cohort_roots, role_root, graph, metadata, out):
    """Run the unchanged preparer over disjoint cohort/role manifest universes.

    audit_cohort_pathways.audit quite properly rejects changed genomes within a
    cohort's attempt roots.  A role root is a new partition, not a changed
    cohort, so audit each partition independently and supply their verified
    union while prepare() performs its usual mapping and receipt work.
    """
    # The 18-genome bacteria manifest must precede its archaeal retry subset.
    cohort_roots = sorted(cohort_roots,
                          key=lambda root: len(read_json(root / "manifest.json")["genomes"]),
                          reverse=True)
    original_audit = CONTEXT.audit
    cohort = original_audit(cohort_roots)
    role = original_audit([Path(role_root).resolve()])
    require(cohort["status"] == role["status"] == "verified",
            "complete verified cohort and role acquisitions required")
    require(not set(cohort["observations"]) & set(role["observations"]),
            "cohort and role acquisitions overlap")
    require(cohort["reference_tables"] == role["reference_tables"],
            "cohort and role searches use different reference tables")
    joint = dict(cohort)
    joint["expected_genomes"] = cohort["expected_genomes"] + role["expected_genomes"]
    joint["observed_genomes"] = cohort["observed_genomes"] + role["observed_genomes"]
    joint["missing_genomes"] = []
    joint["attempts"] = cohort["attempts"] + role["attempts"]
    joint["observations"] = {**cohort["observations"], **role["observations"]}
    try:
        CONTEXT.audit = lambda unused_roots: joint
        return CONTEXT.prepare([*cohort_roots, Path(role_root).resolve()], graph, metadata, out)
    finally:
        CONTEXT.audit = original_audit


def prepare_context(cohort_inputs, metadata, role_root, layouts, out):
    cohort_inputs = Path(cohort_inputs).resolve()
    metadata = Path(metadata).resolve()
    role_root = Path(role_root).resolve()
    layouts = Path(layouts).resolve()
    out = Path(out)
    # The role arm receipt writes role_genomes.tsv beside the families directory
    # containing the metadata supplied to this step.
    role_genomes = metadata.parent.parent / "role_genomes.tsv"
    require(role_genomes.is_file(), f"missing role genome table: {role_genomes}")
    roles = role_layout_rows(role_genomes, proteomes=True)
    require(len(roles) == 3, "joint context must have exactly three role genome rows")
    graph, roots, cohort_receipt = context_source_inputs(cohort_inputs)
    result = prepare_disjoint_roots(roots, role_root, graph, metadata, out)
    old_metadata = read_json(cohort_inputs / "organism-metadata.json")
    new_metadata = read_json(out / "organism-metadata.json")
    role_records = {}
    layout_input_paths = []
    fasta_input_paths = []
    for role in roles:
        gid = role["genome_id"]
        require(gid in new_metadata["genomes"], f"joint context metadata lacks role genome {gid}")
        fasta_path, fasta_ids = fasta_identifiers(
            role["proteome"], role["proteome_sha256"], gid)
        layout_path = layouts / f"{role['key']}.layout.tsv"
        layout_ids = layout_identifiers(
            layout_path, fasta_ids, new_metadata["genomes"][gid].get("contigs"), gid)
        unplaced = sorted(fasta_ids - layout_ids)
        new_metadata["genomes"][gid]["unplacedProteins"] = unplaced
        layout_path = layout_path.resolve()
        role_records[gid] = {
            "layout_path": str(layout_path),
            "layout_sha256": digest(layout_path),
            "layout_proteins": len(layout_ids),
            "fasta_proteins": len(fasta_ids),
            "unplaced_proteins": len(unplaced),
        }
        layout_input_paths.append(layout_path)
        fasta_input_paths.append(fasta_path)
    save_json(out / "organism-metadata.json", new_metadata)
    old_groups = read_tsv(cohort_inputs / "pathway-groups.tsv")
    new_groups = read_tsv(out / "pathway-groups.tsv")
    parity = check_context_parity(old_metadata, new_metadata, old_groups, new_groups)
    save_json(out / "context-parity.json", parity)
    if parity["status"] != "verified":
        raise RoleContextError("cohort context parity failed; see context-parity.json")
    require(len(parity["new_genomes"]) == 3 and set(parity["new_genomes"]) == set(role_records),
            "joint context must add exactly the three role genomes")
    metadata_sha256 = digest(out / "organism-metadata.json")
    preparer_receipt_path = out / "receipt.json"
    preparer_receipt = read_json(preparer_receipt_path)
    require(isinstance(preparer_receipt.get("outputs"), dict),
            "context preparer receipt lacks output hashes")
    preparer_receipt["outputs"]["organism-metadata.json"] = metadata_sha256
    save_json(preparer_receipt_path, preparer_receipt)
    mapping = read_json(out / "mapping-audit.json")["genomes"]
    new_counts = {gid: mapping[gid] for gid in parity["new_genomes"]}
    receipt = {
        "schema": CONTEXT_SCHEMA,
        "status": "verified",
        "genomes": result["genomes"],
        "groups": result["groups"],
        "role_genomes": role_records,
        "cohort_parity": parity,
        "new_genome_counts": new_counts,
        "graph": str(graph),
        "cohort_roots": [str(root) for root in roots],
        "role_root": str(role_root),
        "inputs": file_hashes((cohort_inputs / "receipt.json", metadata,
                               role_root / "manifest.json", role_genomes,
                               *layout_input_paths, *fasta_input_paths)),
        "outputs": {"organism-metadata.json": metadata_sha256},
        "preparer_receipt_sha256": digest(preparer_receipt_path),
        "cohort_receipt_sha256": digest(cohort_inputs / "receipt.json"),
    }
    save_json(out / "role-context-receipt.json", receipt)
    return receipt


def replace_flag(command, flag, value):
    positions = [index for index, token in enumerate(command) if token == flag]
    require(len(positions) == 1 and positions[0] + 1 < len(command),
            "command must contain exactly one valued flag: " + flag)
    command[positions[0] + 1] = str(value)


def parse_path_maps(values):
    result = []
    for value in values:
        require("=" in value, "invalid --path-map (expected LOCAL=REMOTE): " + value)
        local, remote = value.split("=", 1)
        require(local and remote, "invalid empty --path-map endpoint: " + value)
        result.append((local, remote))
    return result


def map_path(value, path_maps):
    """Rewrite a path prefix in a scalar or in a --flag=value command token."""
    text = str(value)
    prefix = ""
    if text.startswith("--") and "=" in text:
        prefix, text = text.split("=", 1)
        prefix += "="
    for local, remote in sorted(path_maps, key=lambda pair: len(pair[0]), reverse=True):
        if text == local or text.startswith(local.rstrip("/") + "/"):
            return prefix + remote.rstrip("/") + text[len(local.rstrip("/")):]
    return prefix + text


def make_plans(plans_10406, role_genomes, gaps, context, families, layouts, models_root,
               jar, out, path_maps=(), cohort_claims_dir=None):
    plans_10406 = Path(plans_10406).resolve()
    role_genomes = Path(role_genomes).resolve()
    gaps, context, families = Path(gaps).resolve(), Path(context).resolve(), Path(families).resolve()
    layouts = Path(layouts).resolve()
    out = Path(out)
    cohort_claims = None
    if cohort_claims_dir is not None:
        cohort_claims_dir = Path(cohort_claims_dir).resolve()
        require(cohort_claims_dir.is_dir(),
                f"cohort claims directory does not exist: {cohort_claims_dir}")
        cohort_claims = {"directory": str(cohort_claims_dir), "files": {}}
    plans_root = out / "plans"
    plans_root.mkdir(parents=True, exist_ok=False)
    roles = role_layout_rows(role_genomes)
    role_ids = {row["genome_id"] for row in roles}
    context_receipt_path = context / "role-context-receipt.json"
    context_receipt = read_json(context_receipt_path)
    require(context_receipt.get("schema") == CONTEXT_SCHEMA and
            context_receipt.get("status") == "verified",
            "verified role-genome context receipt required")
    context_roles = context_receipt.get("role_genomes")
    require(isinstance(context_roles, dict) and set(context_roles) == role_ids,
            "context receipt role genomes differ from role_genomes.tsv")
    role_layouts = {}
    for role in roles:
        gid = role["genome_id"]
        layout_path = (layouts / f"{role['key']}.layout.tsv").resolve()
        require(layout_path.is_file(), f"missing layout for {gid}: {layout_path}")
        actual_sha256 = digest(layout_path)
        expected_sha256 = context_roles[gid].get("layout_sha256")
        require(actual_sha256 == expected_sha256,
                f"layout sha256 differs from context receipt for {gid}: {layout_path}")
        role_layouts[gid] = layout_path
    input_paths = [role_genomes, families / "cohort_orthogroups.tsv",
                   families / "cohort_metadata.json", context / "organism-metadata.json",
                   context / "pathway-groups.tsv", context_receipt_path,
                   Path(__file__).resolve()]
    input_paths.extend(role_layouts[gid] for gid in sorted(role_ids))
    input_paths.extend(gaps / f"{gid}.gaps.jsonl" for gid in sorted(role_ids))
    arms = []
    plan_index = {}
    external_inputs = {"jar": str(jar), "models": {}}
    for ctx in CONTEXTS:
        source_manifest = plans_10406 / ctx / "manifest.tsv"
        source_command = plans_10406 / ctx / "command.json"
        input_paths.extend((source_manifest, source_command))
        with source_manifest.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields, cohort_rows = reader.fieldnames, list(reader)
        require(fields == SOURCE.MANIFEST_COLUMNS, f"unexpected {ctx} manifest columns")
        require(len(cohort_rows) == 18, f"{ctx} source manifest must have 18 cohort rows")
        cohort_ids = {row["genome_id"] for row in cohort_rows}
        require(len(cohort_ids) == 18 and not cohort_ids & role_ids,
                f"{ctx} cohort/role genome population is invalid")
        base_command = read_json(source_command)
        require(isinstance(base_command, list) and all(isinstance(v, str) for v in base_command),
                f"{ctx} source command is not a string argv")
        require("--learned-model" not in base_command, f"{ctx} source command is already learned")
        for kind in KINDS:
            arm = plans_root / ctx / kind
            inference = arm / "inference"
            inference.mkdir(parents=True)
            rows = []
            for original in cohort_rows:
                claims = original["claims_path"]
                require(claims.startswith(OLD_CLAIMS_PREFIX),
                        f"unexpected cohort claims path in {ctx}: {claims}")
                gid = original["genome_id"]
                if cohort_claims is None:
                    claims = NEW_CLAIMS_PREFIX + claims[len(OLD_CLAIMS_PREFIX):]
                else:
                    claims = str((cohort_claims_dir / f"{gid}.claims.jsonl").resolve())
                    require(Path(claims).is_file(),
                            f"missing role-fenced cohort claims for {gid}: {claims}")
                    if gid not in cohort_claims["files"]:
                        cohort_claims["files"][gid] = {
                            "path": claims, "sha256": digest(claims)}
                        input_paths.append(claims)
                rows.append(dict(original, claims_path=claims,
                                 out_path=str(inference / f"{gid}.integrated.tsv")))
            template = cohort_rows[0]
            for role in roles:
                gid = role["genome_id"]
                row = dict(template)
                row.update({
                    "genome_id": gid,
                    "claims_path": role["claims_path"],
                    "out_path": str(inference / f"{gid}.integrated.tsv"),
                    "taxon": "NCBITaxon_" + role["taxid"],
                    "essential_profile": "bacteria",
                    "gaps_path": str(gaps / f"{gid}.gaps.jsonl"),
                    "genome_layout_path": str(role_layouts[gid]),
                })
                rows.append(row)
            command = list(base_command)
            if ctx == "SEQONLY":
                fields, rows = SOURCE.sequence_only_manifest(fields, rows)
                command = SOURCE.strip_context(command)
                # strip_context adds a global --enable-priors none safeguard.
                # Every SEQONLY row already carries enable_priors=none, and the
                # arm contract is stricter: is_context_flag must reject every
                # context token, including that otherwise redundant global.
                while "--enable-priors" in command:
                    index = command.index("--enable-priors")
                    del command[index:index + 2]
            replace_flag(command, "-jar", jar)
            replace_flag(command, "--batch-manifest", arm / "manifest.tsv")
            replace_flag(command, "--feature-export-dir", arm / "features")
            replace_flag(command, "--batch-stats-out", inference / "batch_stats.tsv")
            for flag, name in (("--features-out", "rlc-features.tsv"),
                               ("--suggestions-out", "rlc-suggestions.tsv")):
                if flag in command:
                    replace_flag(command, flag, inference / name)
            for flag, path in (("--cohort-metadata", context / "organism-metadata.json"),
                               ("--cohort-orthogroups", families / "cohort_orthogroups.tsv"),
                               ("--cohort-pathways", context / "pathway-groups.tsv")):
                if flag in command:
                    replace_flag(command, flag, path)
            model = "-"
            if kind != "none":
                model = str(Path(models_root) / ctx / "models" / f"model_{kind}.json")
                command.extend(("--learned-model", model))
                external_inputs["models"][f"{ctx}/{kind}"] = model
            if ctx == "SEQONLY":
                forbidden = [token for token in command
                             if token.split("=", 1)[0] != "--learned-model"
                             and SOURCE.is_context_flag(token)]
                require(not forbidden,
                        "SEQONLY command contains context flags: " + repr(forbidden))
            require(len(rows) == 21, f"{ctx}/{kind} manifest does not have 21 rows")
            local_rows = [dict(row) for row in rows]
            cluster_rows = [{column: map_path(value, path_maps) for column, value in row.items()}
                            for row in rows]
            cluster_command = [map_path(token, path_maps) for token in command]
            write_tsv(arm / "manifest.local.tsv", fields, local_rows)
            write_tsv(arm / "manifest.tsv", fields, cluster_rows)
            save_json(arm / "command.json", cluster_command)
            manifest_cluster = map_path(arm / "manifest.tsv", path_maps)
            command_cluster = map_path(arm / "command.json", path_maps)
            inference_cluster = map_path(inference, path_maps)
            arms.append({"ctx": ctx, "kind": kind, "manifest": manifest_cluster,
                         "command": command_cluster, "inference_dir": inference_cluster,
                         "model": model})
            plan_index[f"{ctx}/{kind}"] = {
                "manifest": manifest_cluster,
                "manifest_local": map_path(arm / "manifest.local.tsv", path_maps),
                "command": command_cluster,
                "inference_dir": inference_cluster,
                "model": model,
                "argv": cluster_command,
            }
    write_tsv(plans_root / "arms.tsv", ["ctx", "kind", "manifest", "command", "inference_dir", "model"], arms)
    if cohort_claims is not None:
        require(len(cohort_claims["files"]) == 18,
                "role-fenced cohort claims must cover all 18 cohort genomes")
    output_paths = sorted(path for path in plans_root.rglob("*") if path.is_file())
    report = {
        "schema": "gspa.role-genome-arm-plans/1",
        "status": "verified",
        "arms": len(arms),
        "genomes_per_manifest": 21,
        "path_maps": [{"local": local, "remote": remote} for local, remote in path_maps],
        "plans": plan_index,
        "inputs": file_hashes(input_paths),
        "external_inputs": external_inputs,
        "outputs": {str(path.relative_to(plans_root)): digest(path) for path in output_paths},
    }
    if cohort_claims is not None:
        report["cohort_claims"] = cohort_claims
    save_json(plans_root / "plan_report.json", report)
    return report


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    steps = parser.add_subparsers(dest="step", required=True)
    acquire_parser = steps.add_parser("acquire")
    acquire_parser.add_argument("--search-root", type=Path, required=True)
    acquire_parser.add_argument("--role-genomes", type=Path, required=True)
    acquire_parser.add_argument("--cohort-reference-hashes", type=Path, required=True)
    acquire_parser.add_argument("--out", type=Path, required=True)
    gaps_parser = steps.add_parser("gaps")
    gaps_parser.add_argument("--acquisition", type=Path)
    gaps_parser.add_argument("--ec2go", type=Path, required=True)
    gaps_parser.add_argument("--out", type=Path, required=True)
    gaps_parser.add_argument("--cohort-root", type=Path, action="append")
    gaps_parser.add_argument("--cohort-replay", type=Path, default=DEFAULT_COHORT_REPLAY)
    gaps_parser.add_argument("--parity-only", action="store_true")
    context_parser = steps.add_parser("context")
    context_parser.add_argument("--cohort-inputs", type=Path, required=True)
    context_parser.add_argument("--metadata", type=Path, required=True)
    context_parser.add_argument("--role-root", type=Path, required=True)
    context_parser.add_argument("--layouts", type=Path, required=True)
    context_parser.add_argument("--out", type=Path, required=True)
    plan_parser = steps.add_parser("plan")
    for name in ("plans-10406", "role-genomes", "gaps", "context", "families", "layouts",
                 "models-root", "jar", "out"):
        plan_parser.add_argument("--" + name, type=Path, required=True)
    plan_parser.add_argument("--path-map", action="append", default=[], metavar="LOCAL=REMOTE")
    plan_parser.add_argument("--cohort-claims-dir", type=Path,
                             help="role-fenced DIR/<genome_id>.claims.jsonl for all 18 cohort rows")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        if args.step == "acquire":
            result = acquire(args.search_root, args.role_genomes,
                             args.cohort_reference_hashes, args.out)
            summary = {"status": result["status"], "genomes": len(result["genomes"])}
        elif args.step == "gaps":
            require(args.parity_only or args.acquisition is not None,
                    "--acquisition is required unless --parity-only is used")
            result = prepare_gaps(args.acquisition, args.ec2go, args.out,
                                  args.cohort_root or DEFAULT_COHORT_ROOTS,
                                  args.cohort_replay, args.parity_only)
            if args.parity_only:
                summary = {key: result[key] for key in ("status", "genomes_compared",
                                                        "genomes_equal", "genomes_unequal")}
            else:
                summary = {"status": result["status"], "counts": result["counts"],
                           "parity": result["parity"]}
        elif args.step == "context":
            result = prepare_context(args.cohort_inputs, args.metadata, args.role_root,
                                     args.layouts, args.out)
            summary = {"status": result["status"], "genomes": result["genomes"],
                       "groups": result["groups"], "new_genome_counts": result["new_genome_counts"]}
        else:
            result = make_plans(args.plans_10406, args.role_genomes, args.gaps,
                                args.context, args.families, args.layouts, args.models_root,
                                args.jar, args.out, parse_path_maps(args.path_map),
                                cohort_claims_dir=args.cohort_claims_dir)
            summary = {"status": result["status"], "arms": result["arms"],
                       "genomes_per_manifest": result["genomes_per_manifest"]}
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0
    except (RoleContextError, FileNotFoundError, KeyError, json.JSONDecodeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
