#!/usr/bin/env python3
"""Prepare label-blind joint-cohort batches for the WP4 family blocks.

Gap records reuse the role-genome replay and mandatory cohort parity check.
Context inputs reuse the disjoint-root cohort/role preparer.  Plans preserve
the 10406 command templates except for the explicitly declared WP4 flags.
No functional truth labels are read by this module.

Two plan modes exist.  ``plan --block B --target T`` is the Task L single-block
mode (base arms of one block, 10406 cohort gaps).  ``plan --blocks all|LIST``
is the full design (Task R): one shared ``cohort`` block holds R0/R0-SC and the
``ref10406`` D4 arms; every NCBI family block holds R1/R2 (variant ``base``),
add-one (regime ADD1), S1/S2 shuffle, duplicate, removal and fragmentation arms
built from ``prepare_wp4_perturbations.py`` outputs.  In that mode every
genome's gap records and ``plausibleReactions`` must come from fenced gapsmith
searches (``gaps``/``context`` with ``--search-root SET=ROOT``), bound to the
same observation per genome and to one reference per superkingdom.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import audit_cohort_pathways as AUDIT  # noqa: E402
import prepare_role_genome_context as ROLE  # noqa: E402
import run_source_arm_context_comparison as SOURCE  # noqa: E402

DATA = Path.home() / ".cache" / "gspa-v3-data"
DEFAULT_COHORT_ROOTS = ROLE.DEFAULT_COHORT_ROOTS
DEFAULT_COHORT_REPLAY = ROLE.DEFAULT_COHORT_REPLAY
CONTEXTS = ("UNIVERSAL", "PATHWAY", "ORGANISM")
KINDS = ("none", "logistic", "gbdt")
REGIMES = ("R0", "R0-SC", "R1", "R2")
HOLDOUTS = {"R0": "CLADE", "R0-SC": "STRAIN_CLUSTER",
            "R1": "CLADE", "R2": "STRAIN_CLUSTER"}
PEER_REGIMES = {"R1", "R2"}
LAYOUT_FIELDS = ROLE.LAYOUT_FIELDS
GAPS_SCHEMA = "gspa.wp4-peer-gaps/1"
CONTEXT_SCHEMA = "gspa.wp4-context-inputs/1"
PLAN_SCHEMA = "gspa.wp4-arm-plans/1"
# arms.tsv: one row per analysed (block, regime, ctx, kind, variant); rows of the
# fragmentation family share one integration run (run_id) per target and arm.
ARM_FIELDS = ["block", "target", "regime", "ctx", "kind", "variant", "sealed_genomes",
              "tier", "holdout", "genomes", "manifest", "command", "inference_dir",
              "run_dir", "run_id", "model"]
RUN_FIELDS = ["run_id", "block", "regime", "ctx", "kind", "run_variant", "holdout",
              "genomes", "manifest", "command", "inference_dir", "model", "variants"]
# Output paths deliberately limited to the flags named by the WP4 contract.
ALLOWED_COMMAND_FLAGS = (
    "--batch-manifest", "--feature-export-dir", "--batch-stats-out",
    "--features-out", "--suggestions-out",
    "--cohort-metadata", "--cohort-orthogroups", "--cohort-pathways",
    "--cohort-context-mode", "--cohort-holdout", "--learned-model",
)

# Re-export the sibling's small, well-tested I/O and mapping helpers.
digest = ROLE.digest
save_json = ROLE.save_json
read_json = ROLE.read_json
read_tsv = ROLE.read_tsv
write_tsv = ROLE.write_tsv
file_hashes = ROLE.file_hashes
replace_flag = ROLE.replace_flag
parse_path_maps = ROLE.parse_path_maps
map_path = ROLE.map_path


class WP4BatchError(ValueError):
    """An acquisition, parity, coverage, or plan contract was not satisfied."""


def require(condition, message):
    if not condition:
        raise WP4BatchError(message)


def csv_values(value, allowed=None, label="value"):
    if isinstance(value, str):
        values = [item.strip() for item in value.split(",") if item.strip()]
    else:
        values = list(value)
    require(len(values) == len(set(values)), f"duplicate {label}: {values}")
    if allowed is not None:
        unknown = sorted(set(values) - set(allowed))
        require(not unknown, f"unknown {label}: {', '.join(unknown)}")
    return values


def fasta_identifiers(path, genome_id):
    path = Path(path).resolve()
    require(path.is_file(), f"missing proteome FASTA for {genome_id}: {path}")
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
    require(identifiers, f"empty proteome FASTA for {genome_id}: {path}")
    return identifiers


def read_unplaced(path, genome_id):
    path = Path(path).resolve()
    require(path.is_file(), f"missing unplaced.txt for {genome_id}: {path}")
    values = []
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        protein_id = line.strip()
        require(protein_id and not any(c.isspace() for c in protein_id),
                f"invalid unplaced protein for {genome_id} at line {line_number}")
        values.append(protein_id)
    require(len(values) == len(set(values)), f"duplicate unplaced protein for {genome_id}")
    return values


def verify_layout(path, fasta_ids, contigs, unplaced, genome_id):
    """Validate placed IDs and require layout + unplaced to partition the FASTA."""
    path = Path(path).resolve()
    require(path.is_file(), f"missing layout for {genome_id}: {path}")
    require(isinstance(contigs, dict) and contigs,
            f"metadata contigs missing for {genome_id}")
    unplaced = set(unplaced)
    require(unplaced <= fasta_ids,
            f"unplaced protein absent from FASTA for {genome_id}: "
            + ", ".join(sorted(unplaced - fasta_ids)))
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
    require(not identifiers & unplaced,
            f"protein is both placed and unplaced for {genome_id}: "
            + ", ".join(sorted(identifiers & unplaced)))
    require(identifiers | unplaced == fasta_ids,
            f"layout/unplaced ids do not equal FASTA ids for {genome_id}; missing: "
            + ", ".join(sorted(fasta_ids - identifiers - unplaced)))
    return identifiers


def manifest_rows(peer_root, peers, exact=False):
    peer_root = Path(peer_root).resolve()
    manifest_path = peer_root / "manifest.json"
    manifest = read_json(manifest_path)
    rows = manifest.get("genomes")
    require(isinstance(rows, list), f"peer search manifest has no genomes list: {manifest_path}")
    result = {}
    for row in rows:
        gid = row.get("genome_id")
        require(isinstance(gid, str) and gid and gid not in result,
                f"invalid or duplicate genome in {manifest_path}")
        result[gid] = row
    require(set(peers) <= set(result),
            "peer search manifest lacks: " + ", ".join(sorted(set(peers) - set(result))))
    if exact:
        require(set(result) == set(peers),
                "peer root/metadata population differs; use one disjoint root per family block")
    return result


def normalize_superkingdom(value, genome_id):
    if isinstance(value, dict):
        value = value.get("name")
    require(isinstance(value, str), f"taxonomy superkingdom missing for {genome_id}")
    normalized = value.strip().lower()
    require(normalized in ("bacteria", "archaea"),
            f"invalid taxonomy superkingdom for {genome_id}: {value!r}")
    return normalized


def essential_profile(taxonomy, genome_id, manifest_row=None, fallback=None):
    value = None
    for key in ("superkingdom", "taxonomy", "kingdom"):
        if taxonomy.get(key) is not None:
            value = taxonomy[key]
            break
    if value is None and manifest_row is not None:
        # Current pilot acquisition files predate the explicit superkingdom
        # field; the observed search manifest records the same NCBI partition.
        value = manifest_row.get("taxonomy")
    if value is None:
        value = fallback
    profile = normalize_superkingdom(value, genome_id)
    if fallback is not None:
        require(profile == fallback,
                f"taxonomy superkingdom differs from context receipt for {genome_id}")
    return profile


def assembly_status(taxonomy, genome_id):
    explicit = taxonomy.get("assemblyStatus", taxonomy.get("assembly_status"))
    if explicit is not None:
        status = str(explicit).upper()
        require(status in ("COMPLETE", "DRAFT"),
                f"invalid assembly status for {genome_id}: {explicit!r}")
        return status
    level = taxonomy.get("assembly_level", taxonomy.get("assemblyLevel"))
    require(level is not None, f"taxonomy assembly level/status missing for {genome_id}")
    return "COMPLETE" if str(level).strip().lower() in ("complete", "complete genome") else "DRAFT"


def peer_record(peers_dir, genome_id, manifest_row=None, essential_fallback=None):
    directory = Path(peers_dir).resolve() / genome_id
    names = ("protein.faa", "layout.tsv", "unplaced.txt", "contigs.json",
             "taxonomy.json", "provenance.json")
    paths = {name: directory / name for name in names}
    for name, path in paths.items():
        require(path.is_file(), f"missing acquisition input for {genome_id}: {path}")
    provenance = read_json(paths["provenance.json"])
    if provenance.get("protein_id_scheme") is not None:
        require(provenance["protein_id_scheme"] == "accession__raw",
                f"unexpected protein id scheme for {genome_id}")
    recorded_hashes = provenance.get("file_sha256s", {})
    for name in names[:-1]:
        if name in recorded_hashes:
            require(digest(paths[name]) == recorded_hashes[name],
                    f"acquisition sha256 differs for {genome_id}/{name}")
    fasta_ids = fasta_identifiers(paths["protein.faa"], genome_id)
    taxonomy = read_json(paths["taxonomy.json"])
    contigs = read_json(paths["contigs.json"])
    for contig, record in contigs.items() if isinstance(contigs, dict) else ():
        require(isinstance(record, dict) and type(record.get("length")) is int and
                record["length"] > 0 and type(record.get("closed")) is bool,
                f"invalid contig record for {genome_id}:{contig}")
    unplaced = read_unplaced(paths["unplaced.txt"], genome_id)
    layout_ids = verify_layout(paths["layout.tsv"], fasta_ids, contigs, unplaced, genome_id)
    clade, strain = taxonomy.get("clade"), taxonomy.get("strainCluster")
    require(isinstance(clade, str) and clade.startswith("ncbi_genus_") and
            isinstance(strain, str) and strain.startswith("ncbi_species_"),
            f"taxonomy clade/strainCluster missing or invalid for {genome_id}")
    species = taxonomy.get("species_taxid")
    if species is None and isinstance(taxonomy.get("species"), dict):
        species = taxonomy["species"].get("taxid")
    require(isinstance(species, (str, int)) and str(species).isdigit(),
            f"taxonomy species taxid missing for {genome_id}")
    if manifest_row is not None:
        if manifest_row.get("sha256") is not None:
            require(manifest_row["sha256"] == digest(paths["protein.faa"]),
                    f"search manifest proteome sha256 differs for {genome_id}")
        if manifest_row.get("proteins") is not None:
            require(manifest_row["proteins"] == len(fasta_ids),
                    f"search manifest protein count differs for {genome_id}")
    return {
        "genome_id": genome_id,
        "fasta_path": str(paths["protein.faa"]),
        "fasta_sha256": digest(paths["protein.faa"]),
        "fasta_ids": fasta_ids,
        "layout_path": str(paths["layout.tsv"].resolve()),
        "layout_sha256": digest(paths["layout.tsv"]),
        "layout_ids": layout_ids,
        "unplacedProteins": sorted(unplaced),
        "contigs": contigs,
        "taxonomy": taxonomy,
        "clade": clade,
        "strainCluster": strain,
        "assemblyStatus": assembly_status(taxonomy, genome_id),
        "species_taxid": str(species),
        "essential_profile": essential_profile(
            taxonomy, genome_id, manifest_row, essential_fallback),
        "provenance_path": str(paths["provenance.json"].resolve()),
    }


def public_peer_record(record):
    return {
        "layout_path": record["layout_path"],
        "layout_sha256": record["layout_sha256"],
        "layout_proteins": len(record["layout_ids"]),
        "fasta_path": record["fasta_path"],
        "fasta_sha256": record["fasta_sha256"],
        "fasta_proteins": len(record["fasta_ids"]),
        "unplaced_proteins": len(record["unplacedProteins"]),
        "unplacedProteins": record["unplacedProteins"],
        "species_taxid": record["species_taxid"],
        "essential_profile": record["essential_profile"],
        "clade": record["clade"],
        "strainCluster": record["strainCluster"],
        "assemblyStatus": record["assemblyStatus"],
    }


def verify_observation(peer_root, record):
    gid = record["genome_id"]
    directory = Path(peer_root).resolve() / "observations" / gid
    observation_path = directory / "observation.json"
    require(observation_path.is_file(), f"missing gapsmith observation for {gid}")
    observation = read_json(observation_path)
    require(observation.get("status") == "observed", f"gapsmith search is not observed for {gid}")
    require(observation.get("proteome_sha256") == record["fasta_sha256"],
            f"observation proteome sha256 differs for {gid}")
    output_hashes = observation.get("output_hashes")
    require(isinstance(output_hashes, dict) and output_hashes,
            f"observation has no output hashes for {gid}")
    for relative, expected in sorted(output_hashes.items()):
        path = ROLE.safe_output(directory, relative)
        require(path.is_file(), f"missing observed output for {gid}: {relative}")
        require(digest(path) == expected, f"output sha256 differs for {gid}: {relative}")
    reactions = directory / "outputs" / f"{gid}-all-Reactions.tbl"
    pathways = directory / "outputs" / f"{gid}-all-Pathways.tbl"
    require(reactions.is_file() and pathways.is_file(),
            f"reaction/pathway tables are missing for {gid}")
    reference_path = directory / "reference-hashes.json"
    require(reference_path.is_file(), f"missing reference-hashes.json for {gid}")
    return observation_path, reactions, pathways, reference_path


def prepare_gaps(peer_root, peers_dir, peers, ec2go, out,
                 cohort_roots=DEFAULT_COHORT_ROOTS,
                 cohort_replay=DEFAULT_COHORT_REPLAY, parity_only=False):
    peers = csv_values(peers, label="peer")
    out = Path(out)
    parity = ROLE.cohort_parity(cohort_roots, cohort_replay, ec2go)
    save_json(out / "parity.json", parity)
    if parity["status"] != "verified":
        raise WP4BatchError("cohort gap parity failed; see parity.json")
    if parity_only:
        return parity
    require(peers, "--peers is required unless --parity-only is used")
    peer_root = Path(peer_root).resolve()
    rows = manifest_rows(peer_root, peers)
    mapping = ROLE.load_ec2go(ec2go)
    counts, gap_records, input_paths = {}, {}, [Path(ec2go).resolve(), peer_root / "manifest.json"]
    references = {}
    for gid in peers:
        record = peer_record(peers_dir, gid, rows[gid])
        observation_path, reactions, pathways, reference_path = verify_observation(peer_root, record)
        reference_bytes = reference_path.read_bytes()
        taxonomy = record["essential_profile"]
        if taxonomy in references:
            require(references[taxonomy] == reference_bytes,
                    f"reference hashes differ within {taxonomy} peers")
        references[taxonomy] = reference_bytes
        raw, pwy = read_tsv(reactions), read_tsv(pathways)
        gaps = ROLE.build_gap_records(raw, pwy, mapping)
        ROLE.REPLAY.verify_gap_groups(gaps, ROLE.REPLAY.eligible_groups(raw, pwy))
        target = out / f"{gid}.gaps.jsonl"
        ROLE.write_jsonl(target, gaps)
        gap_records[gid] = {
            "path": str(target.resolve()), "sha256": digest(target),
            "essential_profile": record["essential_profile"],
            "observation_sha256": digest(observation_path),
            "reference_hashes_sha256": digest(reference_path),
        }
        counts[gid] = {
            "gap_records": len(gaps),
            "records_with_go": sum(item["go_term"] is not None for item in gaps),
            "gapseq_guessed_records": sum(item["gapseq_guessed"] for item in gaps),
            "eligible_reaction_groups": len(ROLE.REPLAY.eligible_groups(raw, pwy)),
            "reaction_rows": len(raw), "pathway_rows": len(pwy),
        }
        input_paths.extend((observation_path, reactions, pathways, reference_path,
                            Path(record["fasta_path"]), Path(record["provenance_path"])))
    receipt = {
        "schema": GAPS_SCHEMA, "status": "verified", "counts": counts,
        "peers": gap_records,
        "parity": {key: parity[key] for key in
                   ("genomes_compared", "genomes_equal", "genomes_unequal")},
        "inputs": file_hashes(input_paths),
        "outputs": {f"{gid}.gaps.jsonl": value["sha256"]
                    for gid, value in gap_records.items()},
    }
    save_json(out / "receipt.json", receipt)
    return receipt


def family_coverage(sequence, peers, peer_records):
    path = Path(sequence).resolve() / "cohort_orthogroups.tsv"
    require(path.is_file(), f"missing joint orthogroups: {path}")
    selected = set(peers)
    proteins = {gid: set() for gid in peers}
    seen = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames is not None and
                {"genome_id", "protein_id"} <= set(reader.fieldnames),
                f"joint orthogroups lacks genome_id/protein_id: {path}")
        for line_number, row in enumerate(reader, 2):
            gid, protein = row["genome_id"], row["protein_id"]
            if gid not in selected:
                continue
            require(protein, f"empty family protein id at {path}:{line_number}")
            key = (gid, protein)
            require(key not in seen, f"duplicate peer family membership: {gid}/{protein}")
            seen.add(key)
            require(protein in peer_records[gid]["fasta_ids"],
                    f"family protein absent from FASTA for {gid}: {protein}")
            proteins[gid].add(protein)
    result = {}
    for gid in peers:
        record = peer_records[gid]
        covered = proteins[gid]
        require(covered <= record["layout_ids"] | set(record["unplacedProteins"]),
                f"family protein not in layout or unplacedProteins for {gid}: "
                + ", ".join(sorted(covered - record["layout_ids"] - set(record["unplacedProteins"]))))
        result[gid] = {
            "proteins": covered,
            "familySearchComplete": covered == record["fasta_ids"],
            "familySearchCoverage": len(covered) / len(record["fasta_ids"]),
        }
    return result


def read_claim_proteins(path, genome_id, contract=None):
    proteins = set()
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError as error:
                raise WP4BatchError(f"invalid claim JSON for {genome_id} at line {line_number}") from error
            protein = row.get("protein_id")
            require(isinstance(protein, str) and protein,
                    f"claim lacks protein_id for {genome_id} at line {line_number}")
            if contract is not None:
                require(row.get("metadata", {}).get("predictor_contract") == contract,
                        f"claim predictor contract differs from claims_report for {genome_id} "
                        f"at line {line_number}")
            proteins.add(protein)
    return proteins


def verify_claims(sequence, peers, peer_records):
    sequence = Path(sequence).resolve()
    report_path = sequence / "claims_report.json"
    report = read_json(report_path)
    require(report.get("schema") == "gspa.wp4-peer-claims/1",
            "WP4 peer claims report schema required")
    require(report.get("contract_differs_from_cohort") == [] and
            report.get("contract_gate", {}).get("passed") is True,
            "peer claim contract differs from cohort contract")
    cohort_contract = report.get("cohort_contract")
    if cohort_contract:
        require(report.get("contract") == cohort_contract,
                "claims_report contract and cohort_contract differ")
    genomes = report.get("genomes")
    require(isinstance(genomes, dict), "claims_report has no genomes object")
    result = {}
    for gid in peers:
        recorded = genomes.get(gid)
        require(isinstance(recorded, dict), f"claims_report lacks peer {gid}")
        path = (sequence / "claims" / f"{gid}.claims.jsonl").resolve()
        require(path.is_file(), f"missing peer claims for {gid}: {path}")
        if recorded.get("file"):
            # The report records the path where the bundle was written; a relocated bundle (staged to a
            # cluster mirror) keeps the layout, so compare the bundle-relative path. sha256 binds the bytes.
            recorded_path = Path(recorded["file"])
            require(recorded_path.resolve() == path or
                    recorded_path.parts[-2:] == ("claims", f"{gid}.claims.jsonl"),
                    f"claims_report path differs for {gid}")
        actual_sha = digest(path)
        require(recorded.get("sha256") == actual_sha,
                f"peer claim sha256 differs from claims_report for {gid}")
        proteins = read_claim_proteins(path, gid, report.get("contract"))
        record = peer_records[gid]
        covered = record["layout_ids"] | set(record["unplacedProteins"])
        require(proteins <= covered,
                f"claim protein not in layout or unplacedProteins for {gid}: "
                + ", ".join(sorted(proteins - covered)))
        require(proteins <= record["fasta_ids"],
                f"claim protein absent from FASTA for {gid}")
        with path.open() as handle:
            records = sum(1 for line in handle if line.strip())
        result[gid] = {"path": str(path), "sha256": actual_sha,
                       "proteins": proteins, "records": records}
    return {"report": report, "report_path": report_path, "claims": result}


def prepare_context(cohort_inputs, sequence, peer_root, peers_dir, peers, out):
    peers = csv_values(peers, label="peer")
    require(peers, "context requires at least one peer")
    cohort_inputs = Path(cohort_inputs).resolve()
    sequence = Path(sequence).resolve()
    peer_root = Path(peer_root).resolve()
    out = Path(out)
    rows = manifest_rows(peer_root, peers, exact=True)
    peer_records = {gid: peer_record(peers_dir, gid, rows[gid]) for gid in peers}
    families = family_coverage(sequence, peers, peer_records)
    graph, roots, unused_cohort_receipt = ROLE.context_source_inputs(cohort_inputs)
    old_metadata = read_json(cohort_inputs / "organism-metadata.json")
    metadata_path = sequence / "cohort_metadata.json"
    metadata = read_json(metadata_path)
    old_genomes = old_metadata.get("genomes", {})
    require(set(metadata.get("genomes", {})) == set(old_genomes) | set(peers),
            "sequence metadata population is not cohort plus selected peers")
    for gid, record in old_genomes.items():
        require(metadata["genomes"].get(gid) == record,
                f"sequence metadata changed cohort entry before context preparation: {gid}")
    for gid in peers:
        peer = peer_records[gid]
        expected = {
            "assemblyStatus": peer["assemblyStatus"], "contigs": peer["contigs"],
            "unplacedProteins": peer["unplacedProteins"], "clade": peer["clade"],
            "strainCluster": peer["strainCluster"],
            "familySearchComplete": families[gid]["familySearchComplete"],
            "familySearchCoverage": families[gid]["familySearchCoverage"],
        }
        actual = metadata["genomes"][gid]
        for key, value in expected.items():
            require(actual.get(key) == value,
                    f"sequence metadata {key} differs from acquisition for {gid}")
    result = ROLE.prepare_disjoint_roots(roots, peer_root, graph, metadata_path, out)
    new_metadata = read_json(out / "organism-metadata.json")
    for gid in peers:
        # Reassert the acquisition file, not a derived FASTA-minus-layout guess.
        new_metadata["genomes"][gid]["unplacedProteins"] = peer_records[gid]["unplacedProteins"]
    save_json(out / "organism-metadata.json", new_metadata)
    old_groups = read_tsv(cohort_inputs / "pathway-groups.tsv")
    new_groups = read_tsv(out / "pathway-groups.tsv")
    parity = ROLE.check_context_parity(old_metadata, new_metadata, old_groups, new_groups)
    save_json(out / "context-parity.json", parity)
    require(parity["status"] == "verified", "cohort context parity failed; see context-parity.json")
    require(set(parity["new_genomes"]) == set(peers),
            "joint context must add exactly the selected peers")
    metadata_sha = digest(out / "organism-metadata.json")
    preparer_receipt_path = out / "receipt.json"
    preparer_receipt = read_json(preparer_receipt_path)
    preparer_receipt["outputs"]["organism-metadata.json"] = metadata_sha
    save_json(preparer_receipt_path, preparer_receipt)
    mapping = read_json(out / "mapping-audit.json")["genomes"]
    peer_receipts = {}
    for gid in peers:
        peer_receipts[gid] = public_peer_record(peer_records[gid])
        peer_receipts[gid].update({
            "familySearchComplete": families[gid]["familySearchComplete"],
            "familySearchCoverage": families[gid]["familySearchCoverage"],
        })
    input_paths = [cohort_inputs / "receipt.json", metadata_path, peer_root / "manifest.json",
                   sequence / "cohort_orthogroups.tsv"]
    for gid in peers:
        directory = Path(peers_dir).resolve() / gid
        input_paths.extend(directory / name for name in
                           ("protein.faa", "layout.tsv", "unplaced.txt", "contigs.json",
                            "taxonomy.json", "provenance.json"))
    receipt = {
        "schema": CONTEXT_SCHEMA, "status": "verified", "genomes": result["genomes"],
        "groups": result["groups"], "peer_genomes": peer_receipts,
        "cohort_parity": parity,
        "new_genome_counts": {gid: mapping[gid] for gid in peers},
        "graph": str(graph), "cohort_roots": [str(root) for root in roots],
        "peer_root": str(peer_root), "inputs": file_hashes(input_paths),
        "outputs": {"organism-metadata.json": metadata_sha},
        "preparer_receipt_sha256": digest(preparer_receipt_path),
        "cohort_receipt_sha256": digest(cohort_inputs / "receipt.json"),
    }
    save_json(out / "wp4-context-receipt.json", receipt)
    return receipt


def load_cohort_proteins(layouts, cohort_ids):
    layouts = Path(layouts).resolve()
    combined = layouts / "all_proteins.faa"
    result = {gid: set() for gid in cohort_ids}
    if combined.is_file():
        with combined.open() as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.startswith(">"):
                    continue
                fields = line[1:].split()
                require(fields, f"empty cohort FASTA identifier at line {line_number}")
                namespaced = fields[0]
                require("__" in namespaced,
                        f"cohort FASTA id is not namespaced at line {line_number}: {namespaced}")
                gid, protein = namespaced.split("__", 1)
                if gid in result:
                    require(protein not in result[gid],
                            f"duplicate cohort FASTA protein: {namespaced}")
                    result[gid].add(protein)
    else:
        for gid in cohort_ids:
            fasta = layouts.parent / "genomes" / gid / "protein.faa"
            result[gid] = fasta_identifiers(fasta, gid)
    for gid, proteins in result.items():
        require(proteins, f"cohort FASTA has no proteins for {gid}")
    return result


def argv_segments(command):
    segments, current = [], []
    for token in command:
        if token.startswith("-") and current:
            segments.append(current)
            current = [token]
        else:
            current.append(token)
    if current:
        segments.append(current)
    return segments


def checked_command_diff(template, command):
    """Return valued-flag changes and reject every undeclared argv change."""
    before = argv_segments(template)
    after = argv_segments(command)
    learned = [segment for segment in after if segment[0] == "--learned-model"]
    require(len(learned) <= 1, "command contains repeated --learned-model")
    after_without_learned = [segment for segment in after if segment[0] != "--learned-model"]
    require([segment[0] for segment in before] ==
            [segment[0] for segment in after_without_learned],
            "command diff adds, removes, or reorders a non-learned flag")
    differences = []
    for old, new in zip(before, after_without_learned):
        if old[1:] != new[1:]:
            differences.append({"flag": old[0], "before": old[1:], "after": new[1:]})
    for segment in learned:
        differences.append({"flag": "--learned-model", "before": [], "after": segment[1:]})
    forbidden = sorted({item["flag"] for item in differences} - set(ALLOWED_COMMAND_FLAGS))
    require(not forbidden, "command diff contains disallowed flags: " + ", ".join(forbidden))
    return differences


def strip_seqonly(command):
    """Retain the role-genome SEQONLY stripping helper for later WP4 controls."""
    result = SOURCE.strip_context(list(command))
    while "--enable-priors" in result:
        index = result.index("--enable-priors")
        del result[index:index + 2]
    return result


# Descriptive alias matching the sibling helper's name.
strip_context = strip_seqonly


def restrict_metadata(source, genome_ids, out_path):
    """Write the arm's organism metadata restricted to exactly its manifest genomes.

    ``CohortContextBuilder.build`` rejects metadata that describes a genome which
    is not a row of the arm's manifest (R0 on joint cohort+peers metadata failed
    on Shaheen).  Pathway groups carry no genome ids and orthogroup rows of
    genomes outside the cohort are skipped by the engine, so only the metadata
    is restricted.
    """
    source = Path(source)
    document = read_json(source)
    require(set(document) == {"genomes"},
            f"cohort metadata with keys other than 'genomes' cannot be restricted safely: {source}")
    genomes = document["genomes"]
    ids = list(genome_ids)
    missing = sorted(set(ids) - set(genomes))
    require(not missing, f"cohort metadata {source} lacks manifest genomes: " + ", ".join(missing))
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    save_json(out_path, {"genomes": {gid: genomes[gid] for gid in sorted(ids)}})
    check_arm_metadata(out_path, ids)
    return {"path": str(out_path.resolve()), "sha256": digest(out_path),
            "source": str(source.resolve()), "source_sha256": digest(source), "genomes": len(ids)}


def check_arm_metadata(metadata_path, genome_ids):
    """Fail when an arm's metadata names a genome outside its manifest (or misses one)."""
    genomes = set(read_json(metadata_path).get("genomes", {}))
    extra = sorted(genomes - set(genome_ids))
    require(not extra, f"arm metadata {metadata_path} describes genomes outside its manifest: "
            + ", ".join(extra))
    missing = sorted(set(genome_ids) - genomes)
    require(not missing, f"arm metadata {metadata_path} lacks manifest genomes: " + ", ".join(missing))


def make_plans(plans_10406, peers_dir, gaps, context, sequence, layouts, models_root, out,
               block, target, peers, regimes=REGIMES, contexts=CONTEXTS, kinds=KINDS,
               path_maps=(), cohort_claims_dir=None):
    peers = csv_values(peers, label="peer")
    regimes = csv_values(regimes, REGIMES, "regime")
    contexts = csv_values(contexts, CONTEXTS, "context")
    kinds = csv_values(kinds, KINDS, "kind")
    require(regimes and contexts and kinds, "at least one regime, context, and kind is required")
    require(isinstance(block, str) and block and block not in (".", "..") and
            "/" not in block and "\\" not in block, "invalid --block name")
    require(target and target not in peers, "target must be a cohort genome, not a peer")
    plans_10406 = Path(plans_10406).resolve()
    peers_dir, gaps = Path(peers_dir).resolve(), Path(gaps).resolve()
    context, sequence = Path(context).resolve(), Path(sequence).resolve()
    layouts, models_root = Path(layouts).resolve(), Path(models_root).resolve()
    out = Path(out).resolve()
    block_root = out / block

    context_receipt_path = context / "wp4-context-receipt.json"
    context_receipt = read_json(context_receipt_path)
    require(context_receipt.get("schema") == CONTEXT_SCHEMA and
            context_receipt.get("status") == "verified", "verified WP4 context receipt required")
    context_peers = context_receipt.get("peer_genomes")
    require(isinstance(context_peers, dict) and set(context_peers) == set(peers),
            "context receipt peer population differs from --peers")
    gaps_receipt_path = gaps / "receipt.json"
    gaps_receipt = read_json(gaps_receipt_path)
    require(gaps_receipt.get("schema") == GAPS_SCHEMA and gaps_receipt.get("status") == "verified" and
            set(gaps_receipt.get("peers", {})) == set(peers),
            "verified WP4 gaps for exactly the selected peers required")

    peer_records = {}
    for gid in peers:
        receipt = context_peers[gid]
        peer_records[gid] = peer_record(
            peers_dir, gid, essential_fallback=receipt.get("essential_profile"))
        require(receipt.get("layout_sha256") == peer_records[gid]["layout_sha256"] and
                receipt.get("fasta_sha256") == peer_records[gid]["fasta_sha256"] and
                receipt.get("unplacedProteins") == peer_records[gid]["unplacedProteins"],
                f"acquisition differs from context receipt for {gid}")
    families = family_coverage(sequence, peers, peer_records)
    claims = verify_claims(sequence, peers, peer_records)
    for gid in peers:
        metadata_record = read_json(context / "organism-metadata.json")["genomes"][gid]
        require(metadata_record.get("unplacedProteins") == peer_records[gid]["unplacedProteins"],
                f"context unplacedProteins differs from acquisition for {gid}")
        require(metadata_record.get("familySearchComplete") == families[gid]["familySearchComplete"] and
                metadata_record.get("familySearchCoverage") == families[gid]["familySearchCoverage"],
                f"context family coverage differs from orthogroups for {gid}")
        gap_path = (gaps / f"{gid}.gaps.jsonl").resolve()
        require(gap_path.is_file() and digest(gap_path) == gaps_receipt["peers"][gid]["sha256"],
                f"peer gap file differs from gaps receipt for {gid}")

    families_report_path = sequence / "families_report.json"
    families_report = read_json(families_report_path)
    require(families_report.get("schema") == "gspa.wp4-peer-families/1" and
            families_report.get("cohort_preservation", {}).get("byte_identical_prefix") is True,
            "families report does not certify byte-identical cohort rows")

    templates = {}
    cohort_ids = None
    input_paths = [context_receipt_path, gaps_receipt_path,
                   context / "organism-metadata.json", context / "pathway-groups.tsv",
                   sequence / "cohort_metadata.json", sequence / "cohort_orthogroups.tsv",
                   claims["report_path"], families_report_path, Path(__file__).resolve()]
    cohort_claims = None
    if cohort_claims_dir is not None:
        cohort_claims_dir = Path(cohort_claims_dir).resolve()
        require(cohort_claims_dir.is_dir(),
                f"cohort claims directory does not exist: {cohort_claims_dir}")
        cohort_claims = {"directory": str(cohort_claims_dir), "files": {}}
    for ctx in contexts:
        source_manifest = plans_10406 / ctx / "manifest.tsv"
        source_command = plans_10406 / ctx / "command.json"
        input_paths.extend((source_manifest, source_command))
        with source_manifest.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields, rows = reader.fieldnames, list(reader)
        require(fields == SOURCE.MANIFEST_COLUMNS, f"unexpected {ctx} manifest columns")
        require(rows, f"{ctx} source manifest is empty")
        ids = [row["genome_id"] for row in rows]
        require(len(ids) == len(set(ids)) and target in ids and not set(ids) & set(peers),
                f"{ctx} source cohort population is invalid")
        if cohort_ids is None:
            cohort_ids = ids
        else:
            require(set(ids) == set(cohort_ids), f"{ctx} cohort population differs")
        command = read_json(source_command)
        require(isinstance(command, list) and command and
                all(isinstance(value, str) for value in command),
                f"{ctx} source command is not a string argv")
        require("--learned-model" not in command, f"{ctx} source command is already learned")
        templates[ctx] = {"fields": fields, "rows": rows, "command": command}

    cohort_id_set = set(cohort_ids)
    joint_metadata = read_json(context / "organism-metadata.json")
    require(set(joint_metadata.get("genomes", {})) == cohort_id_set | set(peers),
            "context metadata population is not template cohort plus selected peers")
    cohort_proteins = load_cohort_proteins(layouts, cohort_id_set)
    cohort_layouts, cohort_gaps = {}, {}
    for gid in cohort_ids:
        layout = (layouts / f"{gid}.layout.tsv").resolve()
        metadata_record = joint_metadata["genomes"][gid]
        verify_layout(layout, cohort_proteins[gid], metadata_record.get("contigs"),
                      metadata_record.get("unplacedProteins", []), gid)
        cohort_layouts[gid] = layout
        input_paths.append(layout)
    for original in templates[contexts[0]]["rows"]:
        gid = original["genome_id"]
        gap_path = Path(original["gaps_path"]).resolve()
        require(gap_path.is_file(), f"missing 10406 gap records for {gid}: {gap_path}")
        cohort_gaps[gid] = gap_path
        input_paths.append(gap_path)
    for ctx in contexts:
        for row in templates[ctx]["rows"]:
            require(Path(row["gaps_path"]).resolve() == cohort_gaps[row["genome_id"]],
                    f"10406 gap path differs across contexts for {row['genome_id']}")

    target_sha256 = None
    target_claim_path = None
    for ctx in contexts:
        for original in templates[ctx]["rows"]:
            gid = original["genome_id"]
            if cohort_claims_dir is None:
                claims_path = original["claims_path"]
                require(claims_path.startswith(ROLE.OLD_CLAIMS_PREFIX),
                        f"unexpected cohort claims path in {ctx}: {claims_path}")
                claims_path = ROLE.NEW_CLAIMS_PREFIX + claims_path[len(ROLE.OLD_CLAIMS_PREFIX):]
                if gid != target:
                    continue
            else:
                claims_path = str((cohort_claims_dir / f"{gid}.claims.jsonl").resolve())
                require(Path(claims_path).is_file(),
                        f"missing cohort claims for {gid}: {claims_path}")
                if gid not in cohort_claims["files"]:
                    cohort_claims["files"][gid] = {"path": claims_path, "sha256": digest(claims_path)}
                    input_paths.append(claims_path)
            if gid == target:
                path = Path(claims_path)
                require(path.is_file(),
                        f"target claims must be locally readable for sha256; use --cohort-claims-dir: {path}")
                actual = digest(path)
                if target_sha256 is None:
                    target_sha256, target_claim_path = actual, path.resolve()
                    input_paths.append(path)
                else:
                    require(actual == target_sha256,
                            f"target claim files are not byte-identical across contexts: {path}")
    require(target_sha256 is not None, "target claims were not found")
    if cohort_claims is not None:
        require(set(cohort_claims["files"]) == cohort_id_set,
                "cohort claims directory must cover every cohort genome")
    for gid in peers:
        input_paths.append(Path(claims["claims"][gid]["path"]))
        input_paths.append(gaps / f"{gid}.gaps.jsonl")
        directory = peers_dir / gid
        input_paths.extend(directory / name for name in
                           ("protein.faa", "layout.tsv", "unplaced.txt", "contigs.json",
                            "taxonomy.json", "provenance.json"))

    block_root.mkdir(parents=True, exist_ok=False)
    arms, plan_index = [], {}
    external_models = {}
    for regime in regimes:
        holdout = HOLDOUTS[regime]
        for ctx in contexts:
            template = templates[ctx]
            for kind in kinds:
                arm = block_root / regime / ctx / kind
                inference = arm / "inference"
                inference.mkdir(parents=True)
                rows = []
                for original in template["rows"]:
                    gid = original["genome_id"]
                    if cohort_claims_dir is None:
                        claims_path = ROLE.NEW_CLAIMS_PREFIX + original["claims_path"][len(ROLE.OLD_CLAIMS_PREFIX):]
                    else:
                        claims_path = str((cohort_claims_dir / f"{gid}.claims.jsonl").resolve())
                    rows.append(dict(original, claims_path=claims_path,
                                     out_path=str(inference / f"{gid}.integrated.tsv"),
                                     gaps_path=str(cohort_gaps[gid]),
                                     genome_layout_path=str(cohort_layouts[gid])))
                if regime in PEER_REGIMES:
                    peer_template = next(row for row in template["rows"] if row["genome_id"] == target)
                    for gid in peers:
                        peer = peer_records[gid]
                        row = dict(peer_template)
                        row.update({
                            "genome_id": gid,
                            "claims_path": claims["claims"][gid]["path"],
                            "out_path": str(inference / f"{gid}.integrated.tsv"),
                            "taxon": "NCBITaxon_" + peer["species_taxid"],
                            "essential_profile": peer["essential_profile"],
                            "gaps_path": str((gaps / f"{gid}.gaps.jsonl").resolve()),
                            "genome_layout_path": peer["layout_path"],
                        })
                        rows.append(row)
                command = list(template["command"])
                replace_flag(command, "--batch-manifest", arm / "manifest.tsv")
                replace_flag(command, "--feature-export-dir", arm / "features")
                replace_flag(command, "--batch-stats-out", inference / "batch_stats.tsv")
                replace_flag(command, "--cohort-metadata", context / "organism-metadata.json")
                replace_flag(command, "--cohort-orthogroups", sequence / "cohort_orthogroups.tsv")
                replace_flag(command, "--cohort-pathways", context / "pathway-groups.tsv")
                replace_flag(command, "--cohort-context-mode", ctx)
                replace_flag(command, "--cohort-holdout", holdout)
                for flag, name in (("--features-out", "rlc-features.tsv"),
                                   ("--suggestions-out", "rlc-suggestions.tsv")):
                    if flag in command:
                        replace_flag(command, flag, inference / name)
                model = "-"
                if kind != "none":
                    model = str(models_root / ctx / "models" / f"model_{kind}.json")
                    command.extend(("--learned-model", model))
                    external_models[f"{ctx}/{kind}"] = model
                arm_context = restrict_metadata(context / "organism-metadata.json",
                                                [row["genome_id"] for row in rows],
                                                arm / "context" / "organism-metadata.json")
                replace_flag(command, "--cohort-metadata", arm_context["path"])
                differences = checked_command_diff(template["command"], command)
                expected_rows = len(cohort_ids) + (len(peers) if regime in PEER_REGIMES else 0)
                require(len(rows) == expected_rows,
                        f"{regime}/{ctx}/{kind} manifest has wrong genome count")
                local_rows = [dict(row) for row in rows]
                cluster_rows = [{column: map_path(value, path_maps) for column, value in row.items()}
                                for row in rows]
                cluster_command = [map_path(token, path_maps) for token in command]
                write_tsv(arm / "manifest.local.tsv", template["fields"], local_rows)
                write_tsv(arm / "manifest.tsv", template["fields"], cluster_rows)
                save_json(arm / "command.json", cluster_command)
                key = f"{regime}/{ctx}/{kind}"
                run_dir = f"{block}/{regime}/{ctx}/{kind}"
                arm_record = {
                    "block": block, "target": target, "regime": regime, "ctx": ctx,
                    "kind": kind, "variant": "base", "sealed_genomes": target, "tier": "",
                    "manifest": map_path(arm / "manifest.tsv", path_maps),
                    "command": map_path(arm / "command.json", path_maps),
                    "inference_dir": map_path(inference, path_maps),
                    "run_dir": run_dir, "run_id": run_dir,
                    "model": map_path(model, path_maps), "holdout": holdout,
                    "genomes": len(rows),
                }
                arms.append(arm_record)
                plan_index[key] = dict(arm_record, argv=cluster_command,
                                       target_claim_sha256=target_sha256,
                                       command_diff=differences,
                                       context_files={"organism-metadata.json": arm_context})
    write_tsv(block_root / "arms.tsv", ARM_FIELDS, arms)
    output_paths = sorted(path for path in block_root.rglob("*") if path.is_file())
    report = {
        "schema": PLAN_SCHEMA, "status": "verified", "block": block, "target": target,
        "peers": peers, "regimes": regimes, "contexts": contexts, "kinds": kinds,
        "arms": len(arms), "cohort_genomes": len(cohort_ids),
        "genomes_per_regime": {regime: len(cohort_ids) +
                               (len(peers) if regime in PEER_REGIMES else 0)
                               for regime in regimes},
        "target_claims": {"path": str(target_claim_path), "sha256": target_sha256},
        "peer_claim_contract": {
            "contract": claims["report"].get("contract"),
            "contract_differs_from_cohort": claims["report"]["contract_differs_from_cohort"],
            "gate_passed": True,
        },
        "command_diff_policy": list(ALLOWED_COMMAND_FLAGS),
        "path_maps": [{"local": local, "remote": remote} for local, remote in path_maps],
        "plans": plan_index, "inputs": file_hashes(input_paths),
        "external_inputs": {"models": external_models},
        "outputs": {str(path.relative_to(block_root)): digest(path) for path in output_paths},
    }
    if cohort_claims is not None:
        report["cohort_claims"] = cohort_claims
    save_json(block_root / "plan_report.json", report)
    return report


# ---------------------------------------------------------------------------
# Task R: fenced search sets and the full WP4 design (all blocks, all arms)
# ---------------------------------------------------------------------------

FENCED_CONTEXT_SCHEMA = "gspa.wp4-fenced-context-inputs/1"
DESIGN_SCHEMA = "gspa.wp4-design-plans/1"
PERTURBATION_SCHEMA = "gspa.wp4-perturbation/1"
COHORT_BLOCK = "cohort"
ALL_CONTEXTS = ("SEQONLY",) + CONTEXTS
LEARNED_KINDS = ("logistic", "gbdt")
DESIGN_REGIMES = REGIMES + ("ADD1",)
DESIGN_HOLDOUTS = dict(HOLDOUTS, ADD1="STRAIN_CLUSTER")
FAMILIES = ("base", "ref10406", "addone", "S1", "S2", "dup", "removal", "fragment")
SEQONLY_FORBIDDEN_FLAGS = ("--cohort", "--cohort-metadata", "--cohort-orthogroups",
                           "--cohort-pathways", "--cohort-context-mode", "--cohort-holdout",
                           "--features-out", "--suggestions-out", "--rlc-suggester")
VARIANT_TSV_FIELDS = ["variant", "target", "sealed_genomes", "path"]
DUP_ROMAN = ("i", "ii", "iii", "iv", "v")
VARIANT_PATTERNS = {
    "removal": re.compile(r"^removal$"),
    "fragment": re.compile(r"^fragment:[0-9]+(k|bp):[0-9]+$"),
    "S1": re.compile(r"^S1:[0-9]+$"),
    "S2": re.compile(r"^S2:[0-9]+$"),
    "dup": re.compile(r"^dup:(i|ii|iii|iv|v)$"),
}
MODE_FAMILIES = {"removal": ("removal",), "fragment": ("fragment",),
                 "shuffle": ("S1", "S2"), "duplicate": ("dup",)}

# Design §2: NCBI family block -> (target accession, target taxid, 10406 split).
# Sulfolobaceae is one block with three targets (one batch of 18 + 14 genomes).
BLOCK_TABLE = {
    "methanosarcinaceae": [("GCA_000007065.1", "192952", "test")],
    "clostridiaceae": [("GCA_000008765.1", "272562", "test")],
    "pectobacteriaceae": [("GCA_000147055.1", "198628", "test")],
    "magnetospirillaceae": [("GCA_000009985.1", "342108", "test")],
    "pasteurellaceae": [("GCA_000027305.1", "71421", "calibration")],
    "campylobacteraceae": [("GCA_000009085.1", "192222", "calibration")],
    "methanobacteriaceae": [("GCA_000145295.1", "79929", "calibration")],
    "sulfolobaceae": [("GCA_000007005.1", "273057", "train"),
                      ("GCA_000012285.1", "330779", "train"),
                      ("GCA_000011205.1", "273063", "train")],
    "thermococcaceae": [("GCA_000009965.1", "69014", "train")],
    "haloferacaceae": [("GCA_000025685.1", "309800", "train")],
    "synechococcaceae": [("GCA_000012525.1", "1140", "train")],
    "methanococcaceae": [("GCA_000011585.1", "267377", "train")],
    "moraxellaceae": [("GCA_000046845.1", "62977", "train")],
    "acidilobaceae": [("GCA_000011125.1", "272557", "train")],
    "thermoplasmataceae": [("GCA_000195915.1", "273075", "train")],
    "legionellaceae": [("GCA_000008485.1", "272624", "train")],
}


def parse_named_paths(values, label):
    """``NAME=PATH`` values (repeatable flags) -> ordered {name: resolved path}."""
    result = {}
    for value in values or ():
        if isinstance(value, (tuple, list)):
            name, path = value
        else:
            require("=" in str(value), f"invalid {label} (expected NAME=PATH): {value}")
            name, path = str(value).split("=", 1)
        require(name and str(path) and name not in result,
                f"invalid or duplicate {label} name: {name!r}")
        result[name] = Path(path).resolve()
    return result


def parse_perturbation_args(values):
    """``DIR`` or ``DIR=CONTEXT_DIR`` values -> [(dir, context or None)]."""
    result = []
    for value in values or ():
        if isinstance(value, (tuple, list)):
            directory, context = value
        elif "=" in str(value):
            directory, context = str(value).split("=", 1)
        else:
            directory, context = value, None
        result.append((Path(directory).resolve(), Path(context).resolve() if context else None))
    require(len({item[0] for item in result}) == len(result), "duplicate --perturbations directory")
    return result


def strip_plausible(record):
    return {key: value for key, value in record.items() if key != "plausibleReactions"}


def cohort_manifest_index(cohort_roots):
    """genome -> 10406 cohort search proteome sha256 (the bytes every later search must reuse)."""
    index = {}
    for root in cohort_roots:
        path = Path(root).resolve() / "manifest.json"
        if not path.is_file():
            continue
        for row in read_json(path).get("genomes", []):
            gid, sha = row.get("genome_id"), row.get("sha256")
            require(isinstance(gid, str) and isinstance(sha, str), f"invalid cohort manifest row in {path}")
            require(index.get(gid, {}).get("sha256", sha) == sha,
                    f"cohort proteome sha256 differs between 10406 search roots for {gid}")
            index[gid] = {"sha256": sha, "manifest": str(path)}
    return index


def perturbation_proteomes(directories):
    """Copies listed in ``gapsmith_inputs.tsv`` -> expected proteome sha256 (file re-hashed)."""
    index = {}
    for directory in directories:
        path = Path(directory) / "gapsmith_inputs.tsv"
        if not path.is_file():
            continue
        for row in read_tsv(path):
            gid = row.get("genome_id")
            proteome = (Path(directory) / row.get("proteome", "")).resolve()
            require(gid and gid not in index, f"invalid or duplicate gapsmith input genome: {gid}")
            require(proteome.is_file() and digest(proteome) == row.get("sha256"),
                    f"perturbation proteome differs from gapsmith_inputs.tsv for {gid}")
            index[gid] = {"sha256": row["sha256"], "path": str(proteome),
                          "source_genome": row.get("source_genome"),
                          "perturbation": row.get("perturbation")}
    return index


def prepare_set_gaps(search_roots, ec2go, out, peers_dir=None, perturbations=(),
                     cohort_roots=DEFAULT_COHORT_ROOTS, cohort_replay=DEFAULT_COHORT_REPLAY):
    """Gap records for every genome of every named search root (cohort, peers, copies).

    The gap derivation is first proven against the JVM replay of the 10406 cohort
    searches (mandatory parity).  Each observed genome must have searched exactly
    the proteome bytes of its source: the 10406 cohort search manifest (cohort),
    the acquisition (peers) or ``gapsmith_inputs.tsv`` (perturbation copies).
    All genomes of one superkingdom must share one reference manifest.
    """
    search_roots = parse_named_paths(search_roots, "--search-root")
    require(search_roots, "at least one --search-root SET=ROOT is required")
    out = Path(out)
    parity = ROLE.cohort_parity(cohort_roots, cohort_replay, ec2go)
    save_json(out / "parity.json", parity)
    if parity["status"] != "verified":
        raise WP4BatchError("cohort gap parity failed; see parity.json")
    cohort_index = cohort_manifest_index(cohort_roots)
    copy_index = perturbation_proteomes([Path(p).resolve() for p in perturbations or ()])
    mapping = ROLE.load_ec2go(ec2go)
    peers_dir = Path(peers_dir).resolve() if peers_dir else None
    genomes, counts, references, sets = {}, {}, {}, {}
    input_paths = [Path(ec2go).resolve()]
    for name, root in search_roots.items():
        rows = manifest_rows(root, [])
        require(rows, f"search root {name} lists no genomes: {root}")
        input_paths.append(root / "manifest.json")
        sets[name] = {"root": str(root), "genomes": sorted(rows)}
        for gid, row in rows.items():
            require(gid not in genomes, f"genome {gid} appears in two search roots")
            if gid in copy_index:
                kind, expected = "copy", copy_index[gid]["sha256"]
            elif peers_dir is not None and (peers_dir / gid).is_dir():
                kind = "peer"
                expected = peer_record(peers_dir, gid, row)["fasta_sha256"]
            elif gid in cohort_index:
                kind, expected = "cohort", cohort_index[gid]["sha256"]
            else:
                raise WP4BatchError(f"no known proteome for searched genome {gid} in set {name}")
            require(row.get("sha256") == expected,
                    f"search manifest proteome differs from the {kind} proteome for {gid}")
            superkingdom = normalize_superkingdom(row.get("taxonomy"), gid)
            observation_path, reactions, pathways, reference_path = verify_observation(
                root, {"genome_id": gid, "fasta_sha256": expected})
            reference_bytes = reference_path.read_bytes()
            require(references.setdefault(superkingdom, reference_bytes) == reference_bytes,
                    f"reference hashes differ within {superkingdom} genomes across search roots")
            raw, pwy = read_tsv(reactions), read_tsv(pathways)
            gaps = ROLE.build_gap_records(raw, pwy, mapping)
            ROLE.REPLAY.verify_gap_groups(gaps, ROLE.REPLAY.eligible_groups(raw, pwy))
            target = out / f"{gid}.gaps.jsonl"
            ROLE.write_jsonl(target, gaps)
            genomes[gid] = {
                "set": name, "kind": kind, "search_root": str(root),
                "path": str(target.resolve()), "sha256": digest(target),
                "essential_profile": superkingdom, "proteome_sha256": expected,
                "observation_sha256": digest(observation_path),
                "reference_hashes_sha256": digest(reference_path),
            }
            counts[gid] = {
                "gap_records": len(gaps),
                "records_with_go": sum(item["go_term"] is not None for item in gaps),
                "eligible_reaction_groups": len(ROLE.REPLAY.eligible_groups(raw, pwy)),
                "reaction_rows": len(raw), "pathway_rows": len(pwy),
            }
            input_paths.extend((observation_path, reactions, pathways, reference_path))
    receipt = {
        "schema": GAPS_SCHEMA, "status": "verified", "mode": "search-sets",
        "sets": sets, "genomes": genomes, "counts": counts,
        "reference_hashes_sha256": {key: _sha256_bytes(value)
                                    for key, value in sorted(references.items())},
        "parity": {key: parity[key] for key in
                   ("genomes_compared", "genomes_equal", "genomes_unequal")},
        "inputs": file_hashes(input_paths),
        "outputs": {f"{gid}.gaps.jsonl": value["sha256"] for gid, value in genomes.items()},
    }
    save_json(out / "receipt.json", receipt)
    return receipt


def _sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def audit_subset(root, genome_ids):
    """``audit_cohort_pathways.audit`` of one search root over a subset of its manifest rows.

    ``staging-hashes.json`` is verified completely; ``check`` runs only for the
    manifest rows named in ``genome_ids``.  The result has the shape of
    ``audit`` with ``expected_genomes`` = the subset size and status
    ``verified`` iff every subset genome is observed.
    """
    root = Path(root).resolve()
    for path, sha in read_json(root / "staging-hashes.json").items():
        if digest(root / path) != sha:
            raise ValueError("staged input changed: " + str(root / path))
    manifest = read_json(root / "manifest.json")
    wanted = set(genome_ids)
    rows = [row for row in manifest["genomes"] if row["genome_id"] in wanted]
    universe = {row["genome_id"] for row in rows}
    attempts, observed, reference_sets, tables = [], {}, {}, None
    for row in rows:
        gid = row["genome_id"]
        directory = root / "observations" / gid
        if not directory.exists():
            continue
        path = directory / "execution.json"
        status = read_json(path).get("status") if path.exists() else "preparing"
        attempts.append(dict(genome_id=gid, root=str(root), status=status))
        if status != "observed":
            continue
        result, refs = AUDIT.check(root, manifest, row)
        tax = row["taxonomy"]
        if tax in reference_sets and reference_sets[tax] != refs:
            raise ValueError("reference bytes differ within taxonomy: " + tax)
        reference_sets[tax] = refs
        shared = {table: refs[table] for table in AUDIT.TABLES}
        if tables is not None and tables != shared:
            raise ValueError("reference tables differ across genomes")
        tables = shared
        if gid in observed:
            raise ValueError("multiple observed runs require explicit comparison: " + gid)
        observed[gid] = result
    missing = sorted(universe - set(observed))
    return dict(schema="gspa.cohort-pathway-audit/1", status="incomplete" if missing else "verified",
                expected_genomes=len(universe), observed_genomes=len(observed), missing_genomes=missing,
                attempts=attempts, observations=observed, reference_tables=tables,
                scope="Acquisition integrity only; no biological accuracy or reference-overlap "
                      "exclusion claim.")


def prepare_partitions(cohort_roots, partition_roots, graph, metadata, out, restrict_to=None):
    """Run the unchanged context preparer over N disjoint, separately audited partitions.

    Generalises ``prepare_role_genome_context.prepare_disjoint_roots``: the cohort
    roots (disjoint universes, largest manifest first) are audited one by one,
    every other root is its own partition; observations must be disjoint and all partitions must use the
    same reference tables.  With ``restrict_to`` (genome ids), every non-cohort
    partition is audited and used only for those genomes (``audit_subset``);
    the cohort roots stay unrestricted.
    """
    cohort_roots = sorted((Path(root).resolve() for root in cohort_roots),
                          key=lambda root: len(read_json(root / "manifest.json")["genomes"]),
                          reverse=True)
    partition_roots = [Path(root).resolve() for root in partition_roots]
    original_audit = ROLE.CONTEXT.audit
    if restrict_to is None:
        partition_audits = [original_audit([root]) for root in partition_roots]
    else:
        wanted = set(restrict_to)
        available = {row["genome_id"] for root in cohort_roots
                     for row in read_json(root / "manifest.json")["genomes"]}
        subsets = []
        for root in partition_roots:
            ids = {row["genome_id"] for row in read_json(root / "manifest.json")["genomes"]}
            require(ids & wanted, f"search root contributes no --metadata genome: {root}")
            subsets.append(ids & wanted)
            available |= ids
        missing = sorted(wanted - available)
        require(not missing, "metadata genomes in no search root: " + ", ".join(missing))
        partition_audits = [audit_subset(root, ids) for root, ids in zip(partition_roots, subsets)]
    # Cohort roots are disjoint universes (per superkingdom), not repeated attempts of one
    # universe, so each is audited alone and merged like any other partition.
    audits = [original_audit([root]) for root in cohort_roots] + partition_audits
    require(all(item["status"] == "verified" for item in audits),
            "complete verified acquisitions are required for every search root")
    seen = set()
    for item in audits:
        overlap = seen & set(item["observations"])
        require(not overlap, "search roots overlap: " + ", ".join(sorted(overlap)))
        seen |= set(item["observations"])
        require(item["reference_tables"] == audits[0]["reference_tables"],
                "search roots use different reference tables")
    joint = dict(audits[0])
    joint["expected_genomes"] = sum(item["expected_genomes"] for item in audits)
    joint["observed_genomes"] = sum(item["observed_genomes"] for item in audits)
    joint["missing_genomes"] = []
    joint["attempts"] = [attempt for item in audits for attempt in item["attempts"]]
    joint["observations"] = {gid: value for item in audits
                             for gid, value in item["observations"].items()}
    try:
        ROLE.CONTEXT.audit = lambda unused_roots: joint
        return ROLE.CONTEXT.prepare([*cohort_roots, *partition_roots], graph, metadata, out)
    finally:
        ROLE.CONTEXT.audit = original_audit


def prepare_fenced_context(cohort_inputs, metadata, cohort_roots, search_roots, out,
                           restrict_to_metadata=False):
    """Context inputs (metadata with plausibleReactions, pathway groups) from fenced searches.

    ``metadata`` is the joint metadata of one genome set (sequence
    ``cohort_metadata.json`` for cohort + peers, a perturbation output's
    ``cohort_metadata.json`` for its copies); its population must equal the union
    of the fenced cohort roots and the named search roots.  Only
    ``plausibleReactions`` may change; cohort records must equal the 10406 records
    apart from ``plausibleReactions``, whose differences are reported (D4).
    With ``restrict_to_metadata`` each search root is used only for the
    genomes of ``metadata`` (one perturbation target's copies from a shared
    per-superkingdom root); the receipt records the restriction.
    """
    cohort_inputs = Path(cohort_inputs).resolve()
    metadata = Path(metadata).resolve()
    search_roots = parse_named_paths(search_roots, "--search-root")
    cohort_roots = [Path(root).resolve() for root in cohort_roots or ()]
    require(cohort_roots, "fenced context requires --cohort-root (the fenced cohort search)")
    require(search_roots, "fenced context requires at least one --search-root SET=ROOT")
    graph, old_roots, unused_receipt = ROLE.context_source_inputs(cohort_inputs)
    require(not {root.resolve() for root in old_roots} & set(cohort_roots),
            "fenced context must not reuse the 10406 cohort search roots")
    out = Path(out)
    input_metadata = read_json(metadata)
    old_metadata = read_json(cohort_inputs / "organism-metadata.json")
    if restrict_to_metadata:
        result = prepare_partitions(cohort_roots, list(search_roots.values()), graph, metadata, out,
                                    restrict_to=set(input_metadata["genomes"]))
    else:
        result = prepare_partitions(cohort_roots, list(search_roots.values()), graph, metadata, out)
    new_metadata = read_json(out / "organism-metadata.json")
    require(set(new_metadata["genomes"]) == set(input_metadata["genomes"]),
            "context population differs from the input metadata")
    for gid, record in input_metadata["genomes"].items():
        require(strip_plausible(new_metadata["genomes"][gid]) == strip_plausible(record),
                f"context preparation changed a non-plausibleReactions field for {gid}")
        require(isinstance(new_metadata["genomes"][gid].get("plausibleReactions"), list),
                f"no fenced plausibleReactions for {gid}")
    changed, d4 = [], {}
    for gid, record in old_metadata["genomes"].items():
        require(gid in new_metadata["genomes"], f"fenced context lacks cohort genome {gid}")
        if strip_plausible(new_metadata["genomes"][gid]) != strip_plausible(record):
            changed.append(gid)
        old = set(record.get("plausibleReactions") or ())
        new = set(new_metadata["genomes"][gid]["plausibleReactions"])
        d4[gid] = {"plausible_10406": len(old), "plausible_fenced": len(new),
                   "only_10406": len(old - new), "only_fenced": len(new - old)}
    require(not changed, "fenced context changed cohort conventions: " + ", ".join(sorted(changed)))
    old_groups = ROLE._group_keys(read_tsv(cohort_inputs / "pathway-groups.tsv"))
    new_groups = ROLE._group_keys(read_tsv(out / "pathway-groups.tsv"))
    audit = read_json(out / "acquisition-audit.json")
    observations, partitions = {}, {"cohort": {"roots": [str(r) for r in cohort_roots], "genomes": []}}
    root_names = {str(root): name for name, root in search_roots.items()}
    for name, root in search_roots.items():
        partitions[name] = {"roots": [str(root)], "genomes": []}
    for gid, value in sorted(audit["observations"].items()):
        name = root_names.get(str(Path(value["search_root"]).resolve()), "cohort")
        partitions[name]["genomes"].append(gid)
        observations[gid] = {"set": name, "search_root": value["search_root"],
                             "observation_sha256": value["observation_sha256"],
                             "reference_manifest_sha256": value["reference_manifest_sha256"]}
    metadata_sha = digest(out / "organism-metadata.json")
    receipt = {
        "schema": FENCED_CONTEXT_SCHEMA, "status": "verified", "mode": "fenced",
        "genomes": result["genomes"], "groups": result["groups"],
        "metadata_input": str(metadata), "partitions": partitions,
        "observations": observations,
        "reference_tables": audit.get("reference_tables"),
        "d4": {"cohort_plausible_reactions": d4,
               "pathway_group_rows": {"10406": len(old_groups), "fenced": len(new_groups),
                                      "only_10406": len(old_groups - new_groups),
                                      "only_fenced": len(new_groups - old_groups)},
               "note": "D4: reported, not interpreted"},
        "graph": str(graph),
        "inputs": file_hashes([cohort_inputs / "receipt.json", cohort_inputs / "organism-metadata.json",
                               cohort_inputs / "pathway-groups.tsv", metadata,
                               *[root / "manifest.json" for root in cohort_roots],
                               *[root / "manifest.json" for root in search_roots.values()]]),
        "outputs": {"organism-metadata.json": metadata_sha,
                    "pathway-groups.tsv": digest(out / "pathway-groups.tsv")},
        "preparer_receipt_sha256": digest(out / "receipt.json"),
    }
    if restrict_to_metadata:
        receipt["restriction"] = {"mode": "metadata", "roots": {
            name: {"root_genomes": len(read_json(root / "manifest.json")["genomes"]),
                   "used_genomes": len(partitions[name]["genomes"])}
            for name, root in search_roots.items()}}
    save_json(out / "wp4-context-receipt.json", receipt)
    return receipt


# ------------------------------------------------------------------ design rules
def design_rules(block, info):
    """Design §4/§9 arm matrix of one block: [(family, regime, ctx, kind, dup_variants)]."""
    rules = []
    if block == COHORT_BLOCK:
        rules += [("base", "R0", ctx, kind, None) for ctx in ALL_CONTEXTS for kind in KINDS]
        rules += [("base", "R0-SC", ctx, "none", None) for ctx in CONTEXTS]
        rules += [("ref10406", "R0", ctx, "none", None) for ctx in CONTEXTS]
        return rules
    split, has_peers = info["split"], bool(info["peers"])
    test, calibration = split == "test", split == "calibration"
    if has_peers:
        rules += [("base", "R2", ctx, "none", None) for ctx in CONTEXTS]
        if test:
            rules += [("base", "R2", ctx, kind, None) for ctx in CONTEXTS for kind in LEARNED_KINDS]
            rules += [("base", "R2", "SEQONLY", "none", None), ("base", "R1", "UNIVERSAL", "none", None)]
        if test or calibration:
            rules += [("addone", "ADD1", ctx, "none", None) for ctx in CONTEXTS]
            rules += [("S1", "R2", ctx, "none", None) for ctx in CONTEXTS]
            rules += [("S2", "R2", "PATHWAY", "none", None)]
            rules += [("dup", "R2", "UNIVERSAL", "none", None),
                      ("dup", "R2", "UNIVERSAL", "logistic", ("dup:i",))]
    if test:
        regimes = ("R0", "R2") if has_peers else ("R0",)
        for family in ("removal", "fragment"):
            rules += [(family, regime, ctx, kind, None) for regime in regimes
                      for ctx in CONTEXTS for kind in ("none", "logistic")]
    return rules


def peer_tier(role, peer_clade, target_clade):
    """ND = near_duplicate role; C = same genus (clade) as the target; F otherwise."""
    if role == "near_duplicate":
        return "ND"
    return "C" if peer_clade == target_clade else "F"


def load_blocks(peers_tsv, block_table=None, selected="all"):
    """Family blocks with their targets and the unique peers of their peers.tsv rows."""
    table = block_table if block_table is not None else BLOCK_TABLE
    rows = read_tsv(peers_tsv)
    require(rows and {"target_accession", "peer_accession", "role"} <= set(rows[0]),
            "peers.tsv lacks target_accession/peer_accession/role")
    targets = {target: block for block, items in table.items() for target, _taxid, _split in items}
    unknown = sorted({row["target_accession"] for row in rows} - set(targets))
    require(not unknown, "peers.tsv targets outside the block table: " + ", ".join(unknown))
    names = list(table) if selected in (None, "all", ["all"]) else csv_values(selected, table, "block")
    require(names, "no blocks selected")
    require(COHORT_BLOCK not in table, f"block name {COHORT_BLOCK!r} is reserved")
    blocks = {}
    for name in names:
        require(name and "/" not in name and name not in (".", ".."), f"invalid block name {name!r}")
        items = table[name]
        splits = {split for _target, _taxid, split in items}
        require(len(splits) == 1, f"block {name} mixes 10406 splits")
        peers, roles = [], {}
        for row in rows:
            if targets[row["target_accession"]] != name:
                continue
            peer = row["peer_accession"]
            if peer not in peers:
                peers.append(peer)
            roles[(row["target_accession"], peer)] = row["role"]
        blocks[name] = {"targets": [target for target, _taxid, _split in items],
                        "split": splits.pop(), "peers": peers, "roles": roles}
    return blocks


# ------------------------------------------------------------------ design planner
class DesignPlanner:
    """Plans every arm of the WP4 design for the selected blocks (fail-closed)."""

    def __init__(self, plans_10406, peers_dir, peers_tsv, sequence, layouts, models_root, out,
                 context, gaps, perturbations=(), cohort_claims_dir=None, blocks="all",
                 block_table=None, regimes=DESIGN_REGIMES, contexts=ALL_CONTEXTS, kinds=KINDS,
                 families=FAMILIES, path_maps=(), expect_references=None):
        self.plans_10406 = Path(plans_10406).resolve()
        self.peers_dir = Path(peers_dir).resolve()
        self.peers_tsv = Path(peers_tsv).resolve()
        self.sequence = Path(sequence).resolve()
        self.layouts = Path(layouts).resolve()
        self.models_root = Path(models_root).resolve()
        self.out = Path(out).resolve()
        self.context = Path(context).resolve()
        self.gap_dirs = [Path(path).resolve() for path in
                         ([gaps] if isinstance(gaps, (str, Path)) else gaps)]
        self.perturbation_args = parse_perturbation_args(perturbations)
        self.cohort_claims_dir = Path(cohort_claims_dir).resolve() if cohort_claims_dir else None
        self.regimes = csv_values(regimes, DESIGN_REGIMES, "regime")
        self.contexts = csv_values(contexts, ALL_CONTEXTS, "context")
        self.kinds = csv_values(kinds, KINDS, "kind")
        self.families = csv_values(families, FAMILIES, "family")
        require(self.regimes and self.contexts and self.kinds and self.families,
                "at least one regime, context, kind and family is required")
        self.path_maps = list(path_maps)
        self.expect_references = dict(expect_references or {})
        self.block_table = block_table
        self.selected_blocks = blocks
        self.inputs = [Path(__file__).resolve(), self.peers_tsv]
        self.cache = {}

    # -------------------------------------------------------------- loading
    def load_templates(self):
        self.templates, cohort_ids = {}, None
        for ctx in self.contexts:
            manifest = self.plans_10406 / ctx / "manifest.tsv"
            command_path = self.plans_10406 / ctx / "command.json"
            self.inputs.extend((manifest, command_path))
            with manifest.open(newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                fields, rows = reader.fieldnames, list(reader)
            require(fields == SOURCE.MANIFEST_COLUMNS, f"unexpected {ctx} manifest columns")
            ids = [row["genome_id"] for row in rows]
            require(rows and len(ids) == len(set(ids)), f"{ctx} source cohort population is invalid")
            if cohort_ids is None:
                cohort_ids = ids
            require(set(ids) == set(cohort_ids), f"{ctx} cohort population differs")
            command = read_json(command_path)
            require(isinstance(command, list) and command and
                    all(isinstance(value, str) for value in command),
                    f"{ctx} source command is not a string argv")
            require("--learned-model" not in command, f"{ctx} source command is already learned")
            if ctx == "SEQONLY":
                present = sorted(set(command) & set(SEQONLY_FORBIDDEN_FLAGS))
                require(not present, "SEQONLY template carries context flags: " + ", ".join(present))
            else:
                for flag in ("--cohort-metadata", "--cohort-orthogroups", "--cohort-pathways",
                             "--cohort-context-mode", "--cohort-holdout"):
                    require(command.count(flag) == 1, f"{ctx} template lacks exactly one {flag}")
                require(command[command.index("--cohort-context-mode") + 1] == ctx,
                        f"{ctx} template context mode differs")
            self.templates[ctx] = {"fields": fields, "rows": {row["genome_id"]: row for row in rows},
                                   "order": ids, "command": command}
        self.cohort_ids = cohort_ids
        # Claims are shared by every template; gap records only by the context templates (the SEQONLY
        # template carries no gaps, written as "-").
        first = self.templates[self.contexts[0]]["rows"]
        with_gaps = [ctx for ctx in self.contexts if ctx != "SEQONLY"]
        for ctx in self.contexts:
            for gid, row in self.templates[ctx]["rows"].items():
                require(row["claims_path"] == first[gid]["claims_path"],
                        f"10406 claims path differs across contexts for {gid}")
                if ctx != "SEQONLY":
                    require(row["gaps_path"] == self.templates[with_gaps[0]]["rows"][gid]["gaps_path"],
                            f"10406 gaps path differs across contexts for {gid}")

    def load_blocks(self):
        self.blocks = load_blocks(self.peers_tsv, self.block_table, self.selected_blocks)
        self.targets = [target for info in self.blocks.values() for target in info["targets"]]
        missing = sorted(set(self.targets) - set(self.cohort_ids))
        require(not missing, "block targets missing from the 10406 cohort: " + ", ".join(missing))
        self.peers = []
        for info in self.blocks.values():
            for peer in info["peers"]:
                require(peer not in self.cohort_ids, f"peer {peer} is a cohort genome")
                require(peer not in self.peers, f"peer {peer} belongs to two blocks")
                self.peers.append(peer)

    def load_context(self):
        receipt_path = self.context / "wp4-context-receipt.json"
        receipt = read_json(receipt_path)
        require(receipt.get("schema") == FENCED_CONTEXT_SCHEMA and receipt.get("status") == "verified"
                and receipt.get("mode") == "fenced",
                "a verified fenced WP4 context receipt is required (context --search-root)")
        self.base_metadata_path = (self.context / "organism-metadata.json").resolve()
        self.base_pathways_path = (self.context / "pathway-groups.tsv").resolve()
        require(digest(self.base_metadata_path) == receipt["outputs"]["organism-metadata.json"] and
                digest(self.base_pathways_path) == receipt["outputs"]["pathway-groups.tsv"],
                "base context outputs differ from their receipt")
        self.base_metadata = read_json(self.base_metadata_path)["genomes"]
        needed = set(self.cohort_ids) | set(self.peers)
        require(needed <= set(self.base_metadata),
                "base context lacks: " + ", ".join(sorted(needed - set(self.base_metadata))))
        self.observations = dict(receipt.get("observations") or {})
        self.inputs.extend((receipt_path, self.base_metadata_path, self.base_pathways_path))
        self.base_context_receipt = receipt

    def load_gaps(self):
        self.gap_index = {}
        for directory in self.gap_dirs:
            receipt_path = directory / "receipt.json"
            receipt = read_json(receipt_path)
            require(receipt.get("schema") == GAPS_SCHEMA and receipt.get("status") == "verified",
                    f"verified WP4 gaps receipt required: {receipt_path}")
            entries = receipt.get("genomes") or receipt.get("peers") or {}
            require(isinstance(entries, dict) and entries, f"gaps receipt lists no genomes: {receipt_path}")
            self.inputs.append(receipt_path)
            for gid, entry in entries.items():
                require(gid not in self.gap_index, f"gap records for {gid} in two --gaps directories")
                path = (directory / f"{gid}.gaps.jsonl").resolve()
                require(path.is_file() and digest(path) == entry.get("sha256"),
                        f"gap file differs from its receipt for {gid}")
                self.gap_index[gid] = dict(entry, path=str(path))

    def load_sequence(self):
        records = {}
        for gid in self.peers:
            # Acquisition taxonomy.json files predate the explicit superkingdom field; the verified gaps
            # receipt records the superkingdom of the fenced search root each peer was observed in.
            # Used only when taxonomy.json has no superkingdom (as the pilot's search-manifest fallback).
            observed = getattr(self, "gap_index", {}).get(gid, {}).get("essential_profile")
            records[gid] = peer_record(self.peers_dir, gid, {"taxonomy": observed} if observed else None)
            self.inputs.extend(self.peers_dir / gid / name for name in
                               ("protein.faa", "layout.tsv", "unplaced.txt", "contigs.json",
                                "taxonomy.json", "provenance.json"))
        self.peer_records = records
        families = family_coverage(self.sequence, self.peers, records) if self.peers else {}
        self.claims = verify_claims(self.sequence, self.peers, records) if self.peers else \
            {"report": {}, "report_path": None, "claims": {}}
        for gid in self.peers:
            record = self.base_metadata[gid]
            require(record.get("unplacedProteins") == records[gid]["unplacedProteins"],
                    f"context unplacedProteins differs from acquisition for {gid}")
            require(record.get("familySearchComplete") == families[gid]["familySearchComplete"] and
                    record.get("familySearchCoverage") == families[gid]["familySearchCoverage"],
                    f"context family coverage differs from orthogroups for {gid}")
            require(record.get("clade") == records[gid]["clade"] and
                    record.get("strainCluster") == records[gid]["strainCluster"],
                    f"context clade/strainCluster differs from acquisition for {gid}")
        families_report = read_json(self.sequence / "families_report.json")
        require(families_report.get("schema") == "gspa.wp4-peer-families/1" and
                families_report.get("cohort_preservation", {}).get("byte_identical_prefix") is True,
                "families report does not certify byte-identical cohort rows")
        self.base_orthogroups_path = (self.sequence / "cohort_orthogroups.tsv").resolve()
        self.inputs.extend((self.sequence / "families_report.json", self.base_orthogroups_path,
                            self.sequence / "cohort_metadata.json"))
        if self.claims["report_path"] is not None:
            self.inputs.append(self.claims["report_path"])
        self.sequence_metadata = read_json(self.sequence / "cohort_metadata.json")["genomes"]
        # Tiers (label-free): ND by role, C by genus clade, F otherwise.
        self.tiers = {}
        for block, info in self.blocks.items():
            for (target, peer), role in info["roles"].items():
                self.tiers[(target, peer)] = peer_tier(
                    role, records[peer]["clade"], self.base_metadata[target].get("clade"))

    def load_cohort(self):
        self.cohort_proteins = load_cohort_proteins(self.layouts, set(self.cohort_ids))
        self.genomes = {}
        if self.cohort_claims_dir is not None:
            require(self.cohort_claims_dir.is_dir(),
                    f"cohort claims directory does not exist: {self.cohort_claims_dir}")
        self.target_claims = {}
        first = self.templates[self.contexts[0]]["rows"]
        gap_rows = self.templates[next((c for c in self.contexts if c != "SEQONLY"), self.contexts[0])]["rows"]
        for gid in self.cohort_ids:
            layout = (self.layouts / f"{gid}.layout.tsv").resolve()
            if self.cohort_claims_dir is not None:
                claims = (self.cohort_claims_dir / f"{gid}.claims.jsonl").resolve()
                require(claims.is_file(), f"missing cohort claims for {gid}: {claims}")
            else:
                original = first[gid]["claims_path"]
                require(original.startswith(ROLE.OLD_CLAIMS_PREFIX),
                        f"unexpected cohort claims path: {original}")
                claims = Path(ROLE.NEW_CLAIMS_PREFIX + original[len(ROLE.OLD_CLAIMS_PREFIX):])
            readable = claims.is_file()
            if gid in self.targets:
                require(readable, f"target claims must be locally readable for sha256; "
                                  f"use --cohort-claims-dir: {claims}")
                self.target_claims[gid] = {"path": str(claims), "sha256": digest(claims)}
            if readable:
                self.inputs.append(claims)
            self.genomes[gid] = {
                "kind": "cohort", "fasta_ids": self.cohort_proteins[gid], "layout": str(layout),
                "claims": str(claims), "claims_readable": readable, "template": gid,
                "source": None, "family": gid, "gaps_fenced": self.gap_index.get(gid, {}).get("path"),
                "gaps_10406": str(Path(gap_rows[gid]["gaps_path"]).resolve()),
            }
            self.inputs.append(layout)
        for block, info in self.blocks.items():
            for gid in info["peers"]:
                record = self.peer_records[gid]
                self.genomes[gid] = {
                    "kind": "peer", "fasta_ids": record["fasta_ids"], "layout": record["layout_path"],
                    "claims": self.claims["claims"][gid]["path"], "claims_readable": True,
                    "template": info["targets"][0], "source": None, "family": gid,
                    "taxon": "NCBITaxon_" + record["species_taxid"],
                    "essential_profile": record["essential_profile"],
                    "unplaced": record["unplacedProteins"], "contigs": record["contigs"],
                    "fasta_path": record["fasta_path"],
                    "gaps_fenced": self.gap_index.get(gid, {}).get("path"),
                }

    # -------------------------------------------------------------- perturbations
    def load_perturbations(self):
        self.perturbations, self.ledger = [], []
        for directory, context in self.perturbation_args:
            receipt_path = directory / "receipt.json"
            receipt = read_json(receipt_path)
            require(receipt.get("schema") == PERTURBATION_SCHEMA and receipt.get("mode") in MODE_FAMILIES,
                    f"not a WP4 perturbation output: {directory}")
            actual = {path.relative_to(directory).as_posix(): digest(path)
                      for path in sorted(directory.rglob("*")) if path.is_file()}
            actual.pop("receipt.json", None)
            require(actual == receipt.get("outputs"),
                    f"perturbation outputs differ from their receipt: {directory}")
            with (directory / "variants.tsv").open(newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                require(reader.fieldnames == VARIANT_TSV_FIELDS,
                        f"variants.tsv header must be {' '.join(VARIANT_TSV_FIELDS)}: {directory}")
                variants = list(reader)
            mode = receipt["mode"]
            for line, row in enumerate(variants, 2):
                family = next((name for name in MODE_FAMILIES[mode]
                               if VARIANT_PATTERNS[name].match(row["variant"] or "")), None)
                require(family is not None,
                        f"{directory}/variants.tsv:{line}: variant {row['variant']!r} is not a {mode} variant")
                row["family"] = family
            item = {"dir": directory, "mode": mode, "receipt": receipt, "variants": variants,
                    "context": context, "metadata": None}
            self.inputs.append(receipt_path)
            if mode in ("removal", "fragment"):
                require(context is not None,
                        f"{mode} perturbations need their fenced context: --perturbations DIR=CONTEXT")
                item["metadata"] = self.perturbation_context(directory, context)
            self.perturbations.append(item)
        self.copy_metadata_overrides = {}

    def perturbation_context(self, directory, context):
        receipt_path = context / "wp4-context-receipt.json"
        receipt = read_json(receipt_path)
        require(receipt.get("schema") == FENCED_CONTEXT_SCHEMA and receipt.get("status") == "verified"
                and receipt.get("mode") == "fenced", f"verified fenced context required: {context}")
        require(Path(receipt.get("metadata_input", "")).resolve() ==
                (directory / "cohort_metadata.json").resolve(),
                f"context {context} was not prepared from {directory}/cohort_metadata.json")
        metadata_path = (context / "organism-metadata.json").resolve()
        require(digest(metadata_path) == receipt["outputs"]["organism-metadata.json"],
                f"perturbation context metadata differs from its receipt: {context}")
        genomes = read_json(metadata_path)["genomes"]
        source = read_json(directory / "cohort_metadata.json")["genomes"]
        for gid, record in genomes.items():
            if gid in self.base_metadata:
                require(record == self.base_metadata[gid],
                        f"perturbation context record differs from the base context for {gid}")
            else:
                require(strip_plausible(record) == strip_plausible(source.get(gid, {})),
                        f"perturbation context changed copy metadata for {gid}")
        for gid, value in (receipt.get("observations") or {}).items():
            if gid in self.observations:
                require(self.observations[gid]["observation_sha256"] == value["observation_sha256"],
                        f"perturbation context used a different search of {gid}")
            else:
                self.observations[gid] = value
        self.inputs.extend((receipt_path, metadata_path))
        groups_equal = digest(context / "pathway-groups.tsv") == digest(self.base_pathways_path)
        return {"path": str(metadata_path), "genomes": genomes, "pathway_groups_equal_base": groups_equal}

    def register_copy(self, item, copy_id, relative):
        """A perturbation copy genome: files, source, claims subset, exact-copy checks."""
        directory = item["dir"]
        base = (directory / "genomes" / relative / copy_id).resolve()
        claims = (directory / "claims" / relative / f"{copy_id}.claims.jsonl").resolve()
        for name in ("protein.faa", "layout.tsv", "unplaced.txt", "contigs.json", "taxonomy.json"):
            require((base / name).is_file(), f"missing perturbation file for {copy_id}: {base / name}")
        require(claims.is_file(), f"missing perturbation claims for {copy_id}: {claims}")
        key = (str(base), str(claims))
        if copy_id in self.genomes:
            existing = self.genomes[copy_id]
            require(existing["kind"] == "copy" and existing["files_key"] == key,
                    f"copy id {copy_id} collides with another genome")
            return existing
        taxonomy = read_json(base / "taxonomy.json")
        source = taxonomy.get("source_genome")
        require(source in self.genomes and self.genomes[source]["kind"] != "copy",
                f"copy {copy_id} names an unknown source genome {source!r}")
        fasta = dict(parse_fasta_file(base / "protein.faa"))
        fasta_ids = set(fasta)
        source_record = self.genomes[source]
        unplaced = read_unplaced(base / "unplaced.txt", copy_id)
        contigs = read_json(base / "contigs.json")
        # Claims of a copy are a sub-multiset of its source's claim lines (same contract).
        copy_lines = claim_lines(claims)
        source_lines = claim_lines(source_record["claims"]) if source_record["claims_readable"] else None
        if source_lines is not None:
            extra = copy_lines - source_lines
            require(not extra, f"copy {copy_id} carries claims absent from its source {source}")
        proteins = {json.loads(line).get("protein_id") for line in copy_lines.elements()}
        require(proteins <= fasta_ids, f"copy {copy_id} claims proteins outside its FASTA")
        source_fasta = self.source_sequences(source)
        require(all(source_fasta.get(pid) == seq for pid, seq in fasta.items()),
                f"copy {copy_id} has proteins that are not exact source proteins of {source}")
        exact = fasta_ids == set(source_fasta)
        record = {
            "kind": "copy", "fasta_ids": fasta_ids, "layout": str(base / "layout.tsv"),
            "claims": str(claims), "claims_readable": True, "template": source_record["template"],
            "source": source, "family": source, "unplaced": unplaced, "contigs": contigs,
            "perturbation": item["mode"], "exact_copy": exact, "files_key": key,
            "gaps_fenced": self.gap_index.get(copy_id, {}).get("path"),
            "copy_lines_equal_source": source_lines is not None and copy_lines == source_lines,
        }
        if source_record["kind"] == "peer":
            record["taxon"] = source_record["taxon"]
            record["essential_profile"] = source_record["essential_profile"]
        self.genomes[copy_id] = record
        self.inputs.extend((base / "protein.faa", base / "layout.tsv", claims))
        return record

    def source_sequences(self, gid):
        key = ("sequences", gid)
        if key not in self.cache:
            record = self.genomes[gid]
            if record["kind"] == "peer":
                self.cache[key] = dict(parse_fasta_file(record["fasta_path"]))
            else:
                self.cache[key] = cohort_sequences(self.layouts, gid)
        return self.cache[key]

    # -------------------------------------------------------------- arm specs
    def arm_specs(self):
        """Every (block, regime, ctx, kind) run of the design with its rows and inputs."""
        specs, self.not_planned = [], []
        cohort_info = {"targets": list(self.targets), "split": "mixed", "peers": []}
        blocks = [(COHORT_BLOCK, cohort_info)] + list(self.blocks.items())
        for block, info in blocks:
            for family, regime, ctx, kind, only in design_rules(block, info):
                if family in ("addone",):
                    for target in info["targets"]:
                        for peer in info["peers"]:
                            if self.tiers.get((target, peer)) != "C":
                                continue
                            specs.append(self.spec(block, info, family, regime, ctx, kind,
                                                   variant=f"addone:{peer}", tier="C",
                                                   population=self.cohort_ids + [peer],
                                                   sealed=[(target, target)]))
                elif family in ("base", "ref10406"):
                    peers = info["peers"] if regime in PEER_REGIMES else []
                    variant = "base" if family == "base" else "ref10406"
                    specs.append(self.spec(block, info, family, regime, ctx, kind, variant=variant,
                                           population=self.cohort_ids + peers,
                                           sealed=[(t, t) for t in info["targets"]]))
                else:
                    specs.extend(self.perturbation_specs(block, info, family, regime, ctx, kind, only))
        return [spec for spec in specs if spec is not None]

    def spec(self, block, info, family, regime, ctx, kind, variant, population, sealed,
             tier="", metadata=None, orthogroups=None, pathways=None, row_variants=None,
             run_variant=None, ledger=None):
        return {"block": block, "targets": info["targets"], "family": family, "regime": regime,
                "ctx": ctx, "kind": kind, "variant": variant, "tier": tier,
                "population": list(population), "sealed": sealed,
                "metadata": metadata, "orthogroups": orthogroups, "pathways": pathways,
                "row_variants": row_variants or [(variant, sealed, tier)],
                "run_variant": run_variant if run_variant is not None else variant,
                "ledger": ledger or []}

    def perturbation_specs(self, block, info, family, regime, ctx, kind, only):
        specs = []
        for item in self.perturbations:
            if family not in MODE_FAMILIES[item["mode"]]:
                continue
            rows = [row for row in item["variants"] if row["family"] == family]
            if only:
                rows = [row for row in rows if row["variant"] in only]
            if family == "fragment":
                by_target = defaultdict(list)
                for row in rows:
                    if row["target"] in info["targets"]:
                        by_target[row["target"]].append(row)
                for target, group in sorted(by_target.items()):
                    specs.append(self.copy_spec(block, info, item, group, regime, ctx, kind,
                                                run_variant="fragment"))
                continue
            for row in rows:
                if family in ("S1", "S2"):
                    require(row["target"] == "*" and not row["sealed_genomes"],
                            f"{family} variant rows must have target '*' and no sealed genomes")
                    specs.append(self.shuffle_spec(block, info, item, row, regime, ctx, kind))
                elif row["target"] in info["targets"]:
                    specs.append(self.copy_spec(block, info, item, [row], regime, ctx, kind))
        return specs

    def shuffle_spec(self, block, info, item, row, regime, ctx, kind):
        path = (item["dir"] / row["path"]).resolve()
        require(path.is_file() and str(path).startswith(str(item["dir"])),
                f"{row['variant']} file is missing or outside its output: {path}")
        shuffled = set(item["receipt"].get("parameters", {}).get("peers") or ())
        missing = sorted(set(info["peers"]) - shuffled)
        sealed = [(t, t) for t in info["targets"]]
        ledger_key = (str(item["dir"]), row["variant"], row["target"])
        if row["family"] == "S1" and missing:
            self.ledger_skip(ledger_key, block, "S1 did not shuffle block peers: " + ", ".join(missing))
            return None
        population = self.cohort_ids + info["peers"]
        if row["family"] == "S1":
            self.check_s1(path, set(info["peers"]))
            return self.spec(block, info, "S1", regime, ctx, kind, row["variant"], population, sealed,
                             orthogroups=str(path), ledger=[ledger_key])
        self.check_s2(path)
        return self.spec(block, info, "S2", regime, ctx, kind, row["variant"], population, sealed,
                         pathways=str(path), ledger=[ledger_key])

    def copy_spec(self, block, info, item, rows, regime, ctx, kind, run_variant=None):
        """Removal, fragment (all fragment variants of one target in one run) and duplicate arms."""
        family = rows[0]["family"]
        target = rows[0]["target"]
        relative = Path(rows[0]["path"] or ".")
        require(not relative.is_absolute() and ".." not in relative.parts,
                f"invalid variant path {relative}")
        row_variants, copies, ledger = [], [], []
        for row in rows:
            require(row["target"] == target and Path(row["path"] or ".") == relative,
                    "variants sharing one run must share target and path")
            sealed = []
            for entry in csv_values(row["sealed_genomes"] or target, label="sealed genome"):
                genome, _, mapped = entry.partition("=")
                mapped = mapped or genome
                require(mapped == target, f"sealed genome {entry} does not map to target {target}")
                sealed.append((genome, mapped))
            require((target, target) in sealed, f"{row['variant']} does not seal its target {target}")
            for genome, _ in sealed:
                if genome != target:
                    copies.append(genome)
            tier = row["variant"].split(":")[1] if family == "fragment" else ""
            row_variants.append((row["variant"], sealed, tier))
            ledger.append((str(item["dir"]), row["variant"], row["target"]))
        if family == "dup":
            population, metadata, orthogroups = self.duplicate_inputs(block, info, item, rows[0], relative)
            for gid in population:
                if self.genomes.get(gid, {}).get("kind") == "copy":
                    copies.append(gid)
        else:
            for copy_id in copies:
                self.register_copy(item, copy_id, relative)
            population = self.cohort_ids + (info["peers"] if regime in PEER_REGIMES else []) + \
                sorted(set(copies))
            metadata = item["metadata"]["path"]
            orthogroups = str((item["dir"] / "cohort_orthogroups.tsv").resolve())
            listed = {r["genome_id"] for r in read_tsv(item["dir"] / "gapsmith_inputs.tsv")}
            missing = sorted(set(copies) - listed)
            require(not missing, "copies absent from gapsmith_inputs.tsv: " + ", ".join(missing))
        copies = sorted(set(copies))
        self.check_derived_orthogroups(orthogroups, {gid: self.genomes[gid]["source"] for gid in copies},
                                       population)
        variant = rows[0]["variant"] if len(rows) == 1 else run_variant
        return self.spec(block, info, family, regime, ctx, kind, variant, population,
                         row_variants[0][1], tier=row_variants[0][2], metadata=metadata,
                         orthogroups=orthogroups, row_variants=row_variants,
                         run_variant=run_variant or rows[0]["variant"], ledger=ledger)

    def duplicate_inputs(self, block, info, item, row, relative):
        directory = (item["dir"] / relative).resolve()
        population = [line.strip() for line in (directory / "population.txt").read_text().splitlines()
                      if line.strip()]
        variant_metadata = read_json(directory / "cohort_metadata.json")["genomes"]
        base_population = set(self.cohort_ids) | set(info["peers"])
        copies_dir = item["dir"] / "genomes" / relative
        copy_ids = sorted(path.name for path in copies_dir.iterdir()) if copies_dir.is_dir() else []
        for copy_id in copy_ids:
            self.register_copy(item, copy_id, relative)
        expected = set(base_population) | set(copy_ids)
        if row["variant"] == "dup:v":
            near = [peer for peer in info["peers"] if any(
                self.tiers.get((t, peer)) == "ND" for t in info["targets"])]
            require(len(near) == 1, f"dup:v needs exactly one near-duplicate in block {block}")
            expected -= set(near)
        require(set(population) == expected,
                f"{row['variant']} population differs from the block population of {block}")
        require(set(variant_metadata) == set(population),
                f"{row['variant']} metadata population differs from population.txt")
        composed = {}
        for gid in population:
            record = variant_metadata[gid]
            if gid in copy_ids:
                source = self.genomes[gid]["source"]
                require(self.genomes[gid]["exact_copy"],
                        f"duplicate copy {gid} is not an exact copy of {source}")
                base = self.base_metadata[source]
                require({k: v for k, v in strip_plausible(record).items() if k != "strainCluster"} ==
                        {k: v for k, v in strip_plausible(self.sequence_metadata[source]).items()
                         if k != "strainCluster"},
                        f"duplicate copy {gid} metadata differs from its source {source}")
                composed[gid] = dict(base, strainCluster=record["strainCluster"])
                self.genomes[gid]["gaps_fenced"] = self.genomes[source]["gaps_fenced"]
                self.genomes[gid]["gaps_reused_from"] = source
            else:
                expected_record = {k: v for k, v in strip_plausible(self.sequence_metadata[gid]).items()
                                   if k != "strainCluster"}
                require({k: v for k, v in strip_plausible(record).items() if k != "strainCluster"} ==
                        expected_record, f"{row['variant']} changed metadata of {gid}")
                if record.get("strainCluster") != self.sequence_metadata[gid].get("strainCluster"):
                    require(row["variant"] == "dup:iii",
                            f"{row['variant']} relabels {gid}; only dup:iii may relabel")
                composed[gid] = dict(self.base_metadata[gid], strainCluster=record["strainCluster"])
        slug = variant_slug(row["variant"])
        target_dir = self.out / block / "_inputs" / f"{info['targets'][0]}-{slug}"
        path = target_dir / "organism-metadata.json"
        if not path.exists():
            target_dir.mkdir(parents=True, exist_ok=True)
            save_json(path, {"genomes": {gid: composed[gid] for gid in sorted(composed)}})
        else:
            require(read_json(path)["genomes"] == composed, f"composed metadata changed: {path}")
        orthogroups = str((directory / "cohort_orthogroups.tsv").resolve())
        ordered = [gid for gid in self.cohort_ids if gid in population] + \
            [gid for gid in info["peers"] if gid in population] + \
            sorted(gid for gid in population if gid in copy_ids)
        return ordered, str(path.resolve()), orthogroups

    # -------------------------------------------------------------- derived input checks
    def orthogroup_index(self, path):
        """genome -> [(protein, payload, raw row)]; only the base file stays cached (memory)."""
        key = ("ortho", str(path))
        if key not in self.cache:
            index = defaultdict(list)
            with Path(path).open("rb") as handle:
                header = handle.readline()
                require(header.decode().rstrip("\n").split("\t")[:3] ==
                        ["genome_id", "protein_id", "orthogroup_id"],
                        f"unexpected orthogroup header: {path}")
                for number, raw in enumerate(handle, 2):
                    fields = raw.decode().rstrip("\n").split("\t")
                    require(len(fields) >= 3 and all(fields[:3]),
                            f"invalid orthogroup row at {path}:{number}")
                    index[fields[0]].append((fields[1], tuple(fields[2:]), raw))
            for gid, rows in index.items():
                proteins = [protein for protein, _, _ in rows]
                require(len(proteins) == len(set(proteins)), f"protein has two family rows in {gid}: {path}")
            index = dict(index)
            self.cache[("famprot", str(path))] = {
                gid: frozenset(protein for protein, _, _ in rows) for gid, rows in index.items()}
            if Path(path).resolve() != self.base_orthogroups_path:
                return index
            self.cache[key] = index
        return self.cache[key]

    def family_proteins(self, path):
        key = ("famprot", str(path))
        if key not in self.cache:
            self.orthogroup_index(path)
        return self.cache[key]

    def check_derived_orthogroups(self, path, copies, population):
        """Non-copy rows byte-identical to the base file; copy rows = source rows of kept proteins."""
        key = ("ortho-check", str(path), tuple(sorted(copies.items())), tuple(sorted(population)))
        if key in self.cache:
            return
        derived, base = self.orthogroup_index(path), self.orthogroup_index(self.base_orthogroups_path)
        for gid in population:
            if gid in copies:
                source = {protein: payload for protein, payload, _ in base.get(copies[gid], [])}
                kept = self.genomes[gid]["fasta_ids"]
                rows = {protein: payload for protein, payload, _ in derived.get(gid, [])}
                expected = {protein: payload for protein, payload in source.items() if protein in kept}
                require(rows == expected, f"family rows of copy {gid} differ from its source rows: {path}")
            else:
                require([raw for _, _, raw in derived.get(gid, [])] ==
                        [raw for _, _, raw in base.get(gid, [])],
                        f"family rows of {gid} differ from the base orthogroups: {path}")
        self.cache[key] = True

    def check_s1(self, path, shuffled_peers):
        key = ("s1", str(path), tuple(sorted(shuffled_peers)))
        if key in self.cache:
            return
        derived, base = self.orthogroup_index(path), self.orthogroup_index(self.base_orthogroups_path)
        require(set(derived) == set(base), f"S1 genome set differs from the base orthogroups: {path}")
        for gid, rows in base.items():
            if gid in self.cohort_ids or gid not in shuffled_peers:
                if gid in self.cohort_ids:
                    require([raw for _, _, raw in derived[gid]] == [raw for _, _, raw in rows],
                            f"S1 changed cohort family rows of {gid}: {path}")
                continue
            require({protein for protein, _, _ in derived[gid]} == {protein for protein, _, _ in rows},
                    f"S1 changed the protein set of {gid}: {path}")
            require(Counter(payload for _, payload, _ in derived[gid]) ==
                    Counter(payload for _, payload, _ in rows),
                    f"S1 changed the family multiset of {gid}: {path}")
        self.cache[key] = True

    def check_s2(self, path):
        key = ("s2", str(path))
        if key in self.cache:
            return
        base, derived = read_tsv(self.base_pathways_path), read_tsv(path)
        group = lambda rows: Counter((r["pathway_id"], r["group_id"]) for r in rows)  # noqa: E731
        reaction = lambda rows: Counter(r["reaction_id"] for r in rows)  # noqa: E731
        require(group(derived) == group(base) and reaction(derived) == reaction(base),
                f"S2 groups are not a degree-preserving shuffle of the base pathway groups: {path}")
        self.cache[key] = True

    def ledger_skip(self, key, block, reason):
        self.skips.append({"dir": key[0], "variant": key[1], "target": key[2], "block": block,
                           "reason": reason})

    # -------------------------------------------------------------- manifests
    def manifest_rows(self, spec, inference):
        template = self.templates[spec["ctx"]]
        fenced = spec["family"] != "ref10406"
        rows = []
        for gid in spec["population"]:
            genome = self.genomes[gid]
            base_row = template["rows"][genome["template"]]
            if genome["kind"] == "cohort":
                row = dict(base_row)
            else:
                row = dict(base_row, genome_id=gid)
                if "taxon" in genome:
                    row["taxon"] = genome["taxon"]
                    row["essential_profile"] = genome["essential_profile"]
            if fenced:
                gaps = genome["gaps_fenced"]
                require(gaps, f"no fenced gap records for {gid} (--gaps)")
            else:
                require(genome["kind"] == "cohort", "ref10406 arms hold cohort genomes only")
                gaps = genome["gaps_10406"]
            row.update(claims_path=genome["claims"], out_path=str(inference / f"{gid}.integrated.tsv"),
                       gaps_path=str(gaps), genome_layout_path=genome["layout"])
            rows.append(row)
        return rows

    def arm_inputs(self, spec, ctx_command):
        """(metadata, orthogroups, pathways) paths of one arm."""
        if spec["family"] == "ref10406":
            command = self.templates[spec["ctx"]]["command"]
            if spec["ctx"] == "SEQONLY":
                return None, None, None
            return tuple(command[command.index(flag) + 1] for flag in
                         ("--cohort-metadata", "--cohort-orthogroups", "--cohort-pathways"))
        return (spec["metadata"] or str(self.base_metadata_path),
                spec["orthogroups"] or str(self.base_orthogroups_path),
                spec["pathways"] or str(self.base_pathways_path))

    def metadata_genomes(self, path):
        key = ("meta", str(path))
        if key not in self.cache:
            require(Path(path).is_file(), f"missing cohort metadata: {path}")
            self.cache[key] = read_json(path)["genomes"]
            self.inputs.append(Path(path))
        return self.cache[key]

    def verify_manifest(self, spec, rows, metadata_path, orthogroups_path):
        ids = [row["genome_id"] for row in rows]
        require(len(ids) == len(set(ids)), f"duplicate genome in {spec['block']} manifest")
        for variant, sealed, _tier in spec["row_variants"]:
            absent = sorted({genome for genome, _ in sealed} - set(ids))
            require(not absent, f"sealed genomes absent from the {variant} manifest: " + ", ".join(absent))
        genomes = self.metadata_genomes(metadata_path) if metadata_path else None
        families = self.family_proteins(orthogroups_path) if orthogroups_path else None
        fenced = spec["family"] != "ref10406"
        for row in rows:
            gid = row["genome_id"]
            genome = self.genomes[gid]
            if genome["kind"] == "cohort" and gid in self.target_claims:
                require(row["claims_path"] == self.target_claims[gid]["path"],
                        f"target claims path differs across arms for {gid}")
            if genomes is not None:
                record = genomes.get(gid)
                require(isinstance(record, dict), f"arm metadata lacks {gid}: {metadata_path}")
                if fenced:
                    require(isinstance(record.get("plausibleReactions"), list),
                            f"no fenced plausibleReactions for {gid} in {metadata_path}")
                    self.check_reference(gid)
                key = ("layout", gid, genome["layout"], json.dumps(
                    [record.get("contigs"), record.get("unplacedProteins")], sort_keys=True))
                if key not in self.cache:
                    placed = verify_layout(genome["layout"], genome["fasta_ids"], record.get("contigs"),
                                           record.get("unplacedProteins", []), gid)
                    if "unplaced" in genome:
                        require(sorted(record.get("unplacedProteins", [])) == sorted(genome["unplaced"]),
                                f"metadata unplacedProteins differs from unplaced.txt for {gid}")
                        require(record.get("contigs") == genome["contigs"],
                                f"metadata contigs differ from contigs.json for {gid}")
                    self.cache[key] = placed | set(record.get("unplacedProteins", []))
                covered = self.cache[key]
                if genome["claims_readable"]:
                    claimed = self.claim_proteins(genome["claims"])
                    require(claimed <= covered, f"claim protein not in layout or unplacedProteins "
                                                f"for {gid}: " + ", ".join(sorted(claimed - covered)[:5]))
                if families is not None:
                    family_proteins = families.get(gid, frozenset())
                    require(family_proteins <= covered,
                            f"family protein not in layout or unplacedProteins for {gid}")
        self.check_collisions(spec["block"], ids)

    def check_collisions(self, block, ids):
        """Protein ids may repeat across genomes only inside one copy family (copies keep ids)."""
        key = ("collisions", tuple(sorted(ids)))
        if key in self.cache:
            return
        owners = {}
        for gid in ids:
            genome = self.genomes[gid]
            for protein in genome["fasta_ids"]:
                family = owners.setdefault(protein, genome["family"])
                require(family == genome["family"],
                        f"protein id {protein} collides between unrelated genomes {gid} and the "
                        f"{family} copy family in {block}")
        self.cache[key] = True

    def claim_proteins(self, path):
        key = ("claims", str(path))
        if key not in self.cache:
            self.cache[key] = read_claim_proteins(path, Path(path).name)
        return self.cache[key]

    def check_reference(self, gid):
        """Gap records and plausibleReactions of a genome come from the same fenced search."""
        genome = self.genomes[gid]
        source = genome.get("gaps_reused_from") or gid
        observation = self.observations.get(source)
        gaps = self.gap_index.get(source)
        require(observation is not None, f"no fenced context observation for {source}")
        require(gaps is not None, f"no fenced gap records for {source}")
        require(gaps.get("observation_sha256") == observation["observation_sha256"],
                f"gap records and plausibleReactions of {source} come from different searches")
        reference = observation["reference_manifest_sha256"]
        require(gaps.get("reference_hashes_sha256") in (None, reference),
                f"gap and context reference manifests differ for {source}")
        profile = gaps.get("essential_profile") or "unknown"
        known = self.references.setdefault(profile, reference)
        require(known == reference, f"{profile} genomes were searched with different references")
        expected = self.expect_references.get(profile)
        require(expected in (None, reference),
                f"{profile} reference manifest {reference} differs from the expected fenced reference")

    # -------------------------------------------------------------- writing
    def command(self, spec, run, inference, metadata, orthogroups, pathways):
        ctx, kind = spec["ctx"], spec["kind"]
        template = self.templates[ctx]["command"]
        command = list(template)
        replace_flag(command, "--batch-manifest", run / "manifest.tsv")
        replace_flag(command, "--feature-export-dir", run / "features")
        replace_flag(command, "--batch-stats-out", inference / "batch_stats.tsv")
        if ctx != "SEQONLY":
            replace_flag(command, "--cohort-metadata", metadata)
            replace_flag(command, "--cohort-orthogroups", orthogroups)
            replace_flag(command, "--cohort-pathways", pathways)
            replace_flag(command, "--cohort-context-mode", ctx)
            replace_flag(command, "--cohort-holdout", DESIGN_HOLDOUTS[spec["regime"]])
            for flag, name in (("--features-out", "rlc-features.tsv"),
                               ("--suggestions-out", "rlc-suggestions.tsv")):
                if flag in command:
                    replace_flag(command, flag, inference / name)
        model = "-"
        if kind != "none":
            model = str(self.models_root / ctx / "models" / f"model_{kind}.json")
            command.extend(("--learned-model", model))
        differences = checked_command_diff(template, command)
        if ctx == "SEQONLY":
            present = sorted(set(command) & set(SEQONLY_FORBIDDEN_FLAGS))
            require(not present, "SEQONLY command carries context flags: " + ", ".join(present))
        return command, model, differences

    def selected(self, spec):
        return (spec["regime"] in self.regimes and spec["ctx"] in self.contexts and
                spec["kind"] in self.kinds and spec["family"] in self.families)

    def run(self):
        self.skips, self.references = [], {}
        self.load_templates()
        self.load_blocks()
        self.load_context()
        self.load_gaps()
        self.load_sequence()
        self.load_cohort()
        self.load_perturbations()
        require(not self.out.exists() or not any(self.out.iterdir()),
                f"refusing to plan into a non-empty output directory: {self.out}")
        self.out.mkdir(parents=True, exist_ok=True)
        specs = self.arm_specs()
        arms, runs, plans = [], [], {}
        planned_ledger = defaultdict(list)
        by_block = defaultdict(list)
        external_models = {}
        for spec in specs:
            if not self.selected(spec):
                continue
            slug = "" if spec["run_variant"] == "base" else "variant-" + variant_slug(spec["run_variant"])
            relative = Path(spec["block"], spec["regime"], spec["ctx"], spec["kind"], slug) if slug else \
                Path(spec["block"], spec["regime"], spec["ctx"], spec["kind"])
            run = self.out / relative
            require(not run.exists(), f"two arms plan the same run directory: {relative}")
            inference = run / "inference"
            inference.mkdir(parents=True)
            rows = self.manifest_rows(spec, inference)
            metadata, orthogroups, pathways = self.arm_inputs(spec, None)
            self.verify_manifest(spec, rows, metadata, orthogroups)
            context_files = {}
            if spec["ctx"] != "SEQONLY":
                context_files["organism-metadata.json"] = restrict_metadata(
                    metadata, [row["genome_id"] for row in rows], run / "context" / "organism-metadata.json")
                metadata = context_files["organism-metadata.json"]["path"]
            command, model, differences = self.command(spec, run, inference, metadata, orthogroups, pathways)
            if model != "-":
                external_models[f"{spec['ctx']}/{spec['kind']}"] = model
            cluster_rows = [{column: map_path(value, self.path_maps) for column, value in row.items()}
                            for row in rows]
            fields = self.templates[spec["ctx"]]["fields"]
            write_tsv(run / "manifest.local.tsv", fields, rows)
            write_tsv(run / "manifest.tsv", fields, cluster_rows)
            cluster_command = [map_path(token, self.path_maps) for token in command]
            save_json(run / "command.json", cluster_command)
            run_id = relative.as_posix()
            common = {
                "block": spec["block"], "regime": spec["regime"], "ctx": spec["ctx"],
                "kind": spec["kind"], "holdout": DESIGN_HOLDOUTS[spec["regime"]]
                if spec["ctx"] != "SEQONLY" else "NONE",
                "genomes": len(rows),
                "manifest": map_path(run / "manifest.tsv", self.path_maps),
                "command": map_path(run / "command.json", self.path_maps),
                "inference_dir": map_path(inference, self.path_maps),
                "model": map_path(model, self.path_maps),
            }
            runs.append(dict(common, run_id=run_id, run_variant=spec["run_variant"],
                             variants=",".join(v for v, _, _ in spec["row_variants"])))
            for variant, sealed, tier in spec["row_variants"]:
                row = dict(common, target=",".join(spec["targets"]), variant=variant,
                           sealed_genomes=",".join(g if g == t else f"{g}={t}" for g, t in sealed),
                           tier=tier, run_dir=run_id, run_id=run_id)
                arms.append(row)
                by_block[spec["block"]].append(row)
                arm_id = f"{spec['block']}/{spec['regime']}/{spec['ctx']}/{spec['kind']}/{variant}"
                require(arm_id not in plans, f"duplicate arm id {arm_id}")
                plans[arm_id] = dict(row, family=spec["family"], argv=cluster_command,
                                     command_diff=differences,
                                     inputs={"metadata": metadata, "orthogroups": orthogroups,
                                             "pathways": pathways},
                                     context_files=context_files,
                                     target_claim_sha256={t: self.target_claims[t]["sha256"]
                                                          for t in spec["targets"]})
            for key in spec["ledger"]:
                planned_ledger[key].append(run_id)
        self.check_sealed_keys(arms)
        ledger = self.variant_ledger(planned_ledger)
        self.write_outputs(arms, runs, plans, by_block, ledger, external_models)
        return self.report

    def check_sealed_keys(self, arms):
        keys = set()
        for row in arms:
            for entry in row["sealed_genomes"].split(","):
                genome = entry.split("=", 1)[0]
                key = (genome, row["regime"], row["ctx"], row["kind"], row["variant"])
                require(key not in keys, f"genome/arm key sealed twice: {key}")
                keys.add(key)

    def variant_ledger(self, planned):
        skipped = {(item["dir"], item["variant"], item["target"]): item for item in self.skips}
        ledger = []
        for item in self.perturbations:
            for row in item["variants"]:
                key = (str(item["dir"]), row["variant"], row["target"])
                entry = {"dir": key[0], "mode": item["mode"], "variant": row["variant"],
                         "target": row["target"], "sealed_genomes": row["sealed_genomes"]}
                if planned.get(key):
                    entry.update(status="planned", runs=sorted(set(planned[key])))
                elif key in skipped:
                    entry.update(status="skipped", reason=skipped[key]["reason"])
                else:
                    entry.update(status="skipped", reason=self.skip_reason(item, row))
                ledger.append(entry)
            if item["mode"] == "duplicate":
                for target, summary in sorted(item["receipt"].get("counts", {}).get("targets", {}).items()):
                    for slug, state in sorted(summary.get("variants", {}).items()):
                        if str(state.get("status", "")).startswith("skipped"):
                            ledger.append({"dir": str(item["dir"]), "mode": "duplicate",
                                           "variant": "dup:" + slug.split("_", 1)[0], "target": target,
                                           "status": "skipped",
                                           "reason": "not produced by prepare_wp4_perturbations: "
                                                     + state["status"]})
        for entry in ledger:
            require(entry["status"] == "planned" or entry.get("reason"),
                    f"variant {entry['variant']} of {entry['target']} yields no arm and no reason")
        return ledger

    def skip_reason(self, item, row):
        target = row["target"]
        if target != "*" and target not in self.targets:
            return "target not in the selected blocks"
        blocks = [block for block, info in self.blocks.items()
                  if target == "*" or target in info["targets"]]
        designed = any(rule[0] == row["family"] and (rule[4] is None or row["variant"] in rule[4])
                       for block in blocks for rule in design_rules(block, self.blocks[block]))
        if not designed:
            return "design §4 plans no " + row["family"] + " arms for this target's block"
        return ("all arms filtered out by --regimes/--contexts/--kinds/--families "
                f"(regimes={','.join(self.regimes)}; contexts={','.join(self.contexts)}; "
                f"kinds={','.join(self.kinds)}; families={','.join(self.families)})")

    def write_outputs(self, arms, runs, plans, by_block, ledger, external_models):
        tiers_path = self.out / "peer_tiers.tsv"
        write_tsv(tiers_path, ["target_accession", "peer_accession", "tier"],
                  [{"target_accession": t, "peer_accession": p, "tier": tier}
                   for (t, p), tier in sorted(self.tiers.items())])
        for block, rows in sorted(by_block.items()):
            targets = sorted({t for row in rows for t in row["target"].split(",")})
            block_report = {
                "schema": DESIGN_SCHEMA, "status": "verified", "block": block,
                "target": ",".join(targets),
                "target_claims": (dict(self.target_claims[targets[0]]) if len(targets) == 1 else {}),
                "target_claims_by_genome": {t: self.target_claims[t] for t in targets},
                "peers": [] if block == COHORT_BLOCK else self.blocks[block]["peers"],
                "arms": len(rows),
                "plans": {key: value for key, value in plans.items() if value["block"] == block},
            }
            write_tsv(self.out / block / "arms.tsv", ARM_FIELDS, rows)
            save_json(self.out / block / "plan_report.json", block_report)
        write_tsv(self.out / "arms.tsv", ARM_FIELDS, arms)
        write_tsv(self.out / "runs.tsv", RUN_FIELDS, runs)
        counts = defaultdict(lambda: defaultdict(int))
        run_counts = defaultdict(int)
        for row in arms:
            counts[plans[f"{row['block']}/{row['regime']}/{row['ctx']}/{row['kind']}/{row['variant']}"]
                   ["family"]][row["regime"]] += 1
        for row in runs:
            run_counts[row["regime"]] += 1
        planned_families = {plan["family"] for plan in plans.values()}
        not_planned = [family for family in FAMILIES if family not in planned_families]
        output_paths = sorted(path for path in self.out.rglob("*") if path.is_file())
        self.report = {
            "schema": DESIGN_SCHEMA, "status": "verified", "blocks": sorted(self.blocks),
            "targets": self.targets, "peers": len(self.peers), "cohort_genomes": len(self.cohort_ids),
            "regimes": self.regimes, "contexts": self.contexts, "kinds": self.kinds,
            "families": self.families, "arms": len(arms), "runs": len(runs),
            "arms_by_family_regime": {family: dict(sorted(values.items()))
                                      for family, values in sorted(counts.items())},
            "runs_by_regime": dict(sorted(run_counts.items())),
            "families_without_arms": not_planned,
            "target_claims_by_genome": self.target_claims,
            "peer_claim_contract": {
                "contract": self.claims["report"].get("contract"),
                "contract_differs_from_cohort": self.claims["report"].get("contract_differs_from_cohort"),
                "gate_passed": bool(self.peers)},
            "references": self.references,
            "d4": self.base_context_receipt.get("d4"),
            "perturbation_contexts": {str(item["dir"]): {
                "context": str(item["context"]) if item["context"] else None,
                "pathway_groups_equal_base": (item["metadata"] or {}).get("pathway_groups_equal_base")}
                for item in self.perturbations},
            "variant_ledger": ledger,
            "command_diff_policy": list(ALLOWED_COMMAND_FLAGS),
            "path_maps": [{"local": local, "remote": remote} for local, remote in self.path_maps],
            "external_inputs": {"models": external_models},
            "plans": plans, "inputs": file_hashes(sorted({Path(p) for p in self.inputs})),
            "outputs": {str(path.relative_to(self.out)): digest(path) for path in output_paths},
        }
        save_json(self.out / "plan_report.json", self.report)


def variant_slug(variant):
    return re.sub(r"[^A-Za-z0-9_.-]", "-", variant)


def parse_fasta_file(path):
    entries, identifier, chunks = [], None, []
    with Path(path).open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if identifier is not None:
                    entries.append((identifier, "".join(chunks)))
                fields = line[1:].split()
                require(fields, f"empty FASTA header in {path}")
                identifier, chunks = fields[0], []
            elif line.strip():
                chunks.append(line.strip())
    if identifier is not None:
        entries.append((identifier, "".join(chunks)))
    identifiers = [item[0] for item in entries]
    require(len(identifiers) == len(set(identifiers)), f"duplicate FASTA protein id in {path}")
    return entries


def cohort_sequences(layouts, gid):
    combined = Path(layouts) / "all_proteins.faa"
    result = {}
    for identifier, sequence in parse_fasta_file(combined):
        genome, _, protein = identifier.partition("__")
        if genome == gid:
            result[protein] = sequence
    require(result, f"cohort FASTA has no proteins for {gid}")
    return result


def claim_lines(path):
    with Path(path).open() as handle:
        return Counter(line if line.endswith("\n") else line + "\n" for line in handle if line.strip())


def make_design_plans(**kwargs):
    return DesignPlanner(**kwargs).run()


def design_counts(peers_tsv, peers_dir, sequence, block_table=None, blocks="all",
                  fragment_variants=6, s1_seeds=5, s2_seeds=3):
    """Arm rows and integration runs of the full design without perturbation outputs.

    Assumes the Task P conventions: one removal variant per target, ``fragment_variants``
    fragment variants per target sharing one run, S1/S2 seeds, and dup:i..v unless a C
    peer (i, ii) or the near-duplicate (iii, v) is missing.
    """
    infos = load_blocks(peers_tsv, block_table, blocks)
    metadata = read_json(Path(sequence) / "cohort_metadata.json")["genomes"]
    targets = [t for info in infos.values() for t in info["targets"]]
    rows, runs = defaultdict(lambda: defaultdict(int)), defaultdict(lambda: defaultdict(int))
    cohort = {"targets": targets, "split": "mixed", "peers": []}
    for block, info in [(COHORT_BLOCK, cohort)] + list(infos.items()):
        tiers = {}
        for (target, peer), role in info.get("roles", {}).items():
            clade = read_json(Path(peers_dir) / peer / "taxonomy.json").get("clade")
            tiers[(target, peer)] = peer_tier(role, clade, metadata[target].get("clade"))
        for family, regime, ctx, kind, only in design_rules(block, info):
            if family == "addone":
                n_rows = n_runs = sum(1 for tier in tiers.values() if tier == "C")
            elif family in ("base", "ref10406"):
                n_rows = n_runs = 1
            elif family == "S1":
                n_rows = n_runs = s1_seeds
            elif family == "S2":
                n_rows = n_runs = s2_seeds
            elif family == "removal":
                n_rows = n_runs = len(info["targets"])
            elif family == "fragment":
                n_runs, n_rows = len(info["targets"]), len(info["targets"]) * fragment_variants
            else:  # dup
                n_rows = 0
                for target in info["targets"]:
                    has_c = any(t == target and tier == "C" for (t, _p), tier in tiers.items())
                    has_nd = any(t == target and tier == "ND" for (t, _p), tier in tiers.items())
                    for roman in DUP_ROMAN:
                        if only and f"dup:{roman}" not in only:
                            continue
                        if (roman in ("i", "ii") and not has_c) or (roman in ("iii", "v") and not has_nd):
                            continue
                        n_rows += 1
                n_runs = n_rows
            rows[family][regime] += n_rows
            runs[family][regime] += n_runs
    total_rows = sum(sum(v.values()) for v in rows.values())
    total_runs = sum(sum(v.values()) for v in runs.values())
    return {"blocks": len(infos), "targets": len(targets),
            "arms_by_family_regime": {k: dict(sorted(v.items())) for k, v in sorted(rows.items())},
            "runs_by_family_regime": {k: dict(sorted(v.items())) for k, v in sorted(runs.items())},
            "arms": total_rows, "runs": total_runs}


def read_block_table(path):
    """``block, target_accession, target_taxid, split`` TSV -> block table (overrides §2)."""
    table = {}
    for row in read_tsv(path):
        require(row.get("block") and row.get("target_accession") and row.get("split"),
                f"invalid block table row in {path}")
        table.setdefault(row["block"], []).append(
            (row["target_accession"], row.get("target_taxid", ""), row["split"]))
    require(table, f"empty block table: {path}")
    return table


def parse_expect_references(values):
    result = {}
    for value in values or ():
        require("=" in value, f"invalid --expect-reference (expected SUPERKINGDOM=SHA256): {value}")
        name, sha = value.split("=", 1)
        result[normalize_superkingdom(name, "--expect-reference")] = sha
    return result


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    steps = parser.add_subparsers(dest="step", required=True)

    gaps_parser = steps.add_parser("gaps", help="gap records from gapsmith observations")
    gaps_parser.add_argument("--peer-root", type=Path, help="Task L: one peer search root")
    gaps_parser.add_argument("--peers-dir", type=Path)
    gaps_parser.add_argument("--peers", default="", help="comma-separated peer accessions")
    gaps_parser.add_argument("--search-root", action="append", default=[], metavar="SET=ROOT",
                             help="fenced search root of one genome set (cohort, peers, removal, ...)")
    gaps_parser.add_argument("--perturbations", type=Path, action="append", default=[],
                             help="perturbation output whose gapsmith_inputs.tsv names copy proteomes")
    gaps_parser.add_argument("--ec2go", type=Path, required=True)
    gaps_parser.add_argument("--out", type=Path, required=True)
    gaps_parser.add_argument("--cohort-root", type=Path, action="append",
                             help="10406 cohort search roots (parity and cohort proteome sha256)")
    gaps_parser.add_argument("--cohort-replay", type=Path, default=DEFAULT_COHORT_REPLAY)
    gaps_parser.add_argument("--parity-only", action="store_true")

    context_parser = steps.add_parser("context", help="organism metadata and pathway groups")
    context_parser.add_argument("--cohort-inputs", type=Path, required=True)
    context_parser.add_argument("--sequence", type=Path)
    context_parser.add_argument("--peer-root", type=Path)
    context_parser.add_argument("--peers-dir", type=Path)
    context_parser.add_argument("--peers", default="")
    context_parser.add_argument("--search-root", action="append", default=[], metavar="SET=ROOT",
                                help="fenced mode: search root of one non-cohort genome set")
    context_parser.add_argument("--cohort-root", type=Path, action="append", default=[],
                                help="fenced mode: the fenced cohort search root(s)")
    context_parser.add_argument("--metadata", type=Path,
                                help="fenced mode: joint metadata of the genome set "
                                     "(default SEQUENCE/cohort_metadata.json)")
    context_parser.add_argument("--restrict-to-metadata", action="store_true",
                                help="fenced mode: use each --search-root only for the --metadata "
                                     "genomes (cohort roots stay unrestricted)")
    context_parser.add_argument("--out", type=Path, required=True)

    plan_parser = steps.add_parser("plan", help="arm manifests, commands, arms.tsv, plan_report.json")
    for name in ("plans-10406", "peers-dir", "context", "sequence", "layouts", "models-root", "out"):
        plan_parser.add_argument("--" + name, type=Path, required=True)
    plan_parser.add_argument("--gaps", type=Path, action="append", required=True,
                             help="gaps output directory (repeatable in design mode)")
    plan_parser.add_argument("--block", help="Task L single-block mode: block name")
    plan_parser.add_argument("--target", help="Task L single-block mode: target accession")
    plan_parser.add_argument("--peers", default="")
    plan_parser.add_argument("--blocks", help="design mode: 'all' or comma-separated family blocks")
    plan_parser.add_argument("--block-table", type=Path, help="design mode: override the §2 block table")
    plan_parser.add_argument("--peers-tsv", type=Path, help="design mode: peer-selection peers.tsv")
    plan_parser.add_argument("--perturbations", action="append", default=[], metavar="DIR[=CONTEXT]",
                             help="design mode: prepare_wp4_perturbations.py output (removal/fragment "
                                  "need =CONTEXT, their fenced context)")
    plan_parser.add_argument("--families", default=",".join(FAMILIES))
    plan_parser.add_argument("--expect-reference", action="append", default=[],
                             metavar="SUPERKINGDOM=SHA256")
    plan_parser.add_argument("--regimes")
    plan_parser.add_argument("--contexts")
    plan_parser.add_argument("--kinds", default=",".join(KINDS))
    plan_parser.add_argument("--path-map", action="append", default=[], metavar="LOCAL=REMOTE")
    plan_parser.add_argument("--cohort-claims-dir", type=Path)

    design_parser = steps.add_parser("design", help="arm/run counts of the full design (no files written)")
    design_parser.add_argument("--peers-tsv", type=Path, required=True)
    design_parser.add_argument("--peers-dir", type=Path, required=True)
    design_parser.add_argument("--sequence", type=Path, required=True)
    design_parser.add_argument("--blocks", default="all")
    design_parser.add_argument("--block-table", type=Path)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        if args.step == "gaps":
            roots = args.cohort_root or DEFAULT_COHORT_ROOTS
            if args.search_root:
                require(not args.parity_only and args.peer_root is None,
                        "--search-root cannot be combined with --peer-root or --parity-only")
                result = prepare_set_gaps(args.search_root, args.ec2go, args.out, args.peers_dir,
                                          args.perturbations, roots, args.cohort_replay)
                summary = {"status": result["status"], "sets": {k: len(v["genomes"])
                                                                 for k, v in result["sets"].items()},
                           "parity": result["parity"]}
            else:
                require(args.parity_only or (args.peer_root is not None and args.peers_dir is not None),
                        "--peer-root and --peers-dir are required unless --parity-only is used")
                result = prepare_gaps(args.peer_root, args.peers_dir, args.peers, args.ec2go,
                                      args.out, roots, args.cohort_replay, args.parity_only)
                summary = ({"status": result["status"], "genomes_compared": result["genomes_compared"]}
                           if args.parity_only else
                           {"status": result["status"], "counts": result["counts"],
                            "parity": result["parity"]})
        elif args.step == "context":
            if args.search_root:
                metadata = args.metadata or (args.sequence / "cohort_metadata.json"
                                             if args.sequence else None)
                require(metadata is not None, "fenced context needs --metadata or --sequence")
                if args.restrict_to_metadata:
                    result = prepare_fenced_context(args.cohort_inputs, metadata, args.cohort_root,
                                                    args.search_root, args.out, restrict_to_metadata=True)
                else:
                    result = prepare_fenced_context(args.cohort_inputs, metadata, args.cohort_root,
                                                    args.search_root, args.out)
                summary = {"status": result["status"], "genomes": result["genomes"],
                           "groups": result["groups"],
                           "partitions": {k: len(v["genomes"]) for k, v in result["partitions"].items()}}
                if args.restrict_to_metadata:
                    summary["restriction"] = result["restriction"]
            else:
                require(not args.restrict_to_metadata,
                        "--restrict-to-metadata requires fenced mode (--search-root)")
                require(args.sequence and args.peer_root and args.peers_dir and args.peers,
                        "context needs --sequence, --peer-root, --peers-dir and --peers "
                        "(or --search-root for fenced mode)")
                result = prepare_context(args.cohort_inputs, args.sequence, args.peer_root,
                                         args.peers_dir, args.peers, args.out)
                summary = {"status": result["status"], "genomes": result["genomes"],
                           "groups": result["groups"], "new_genome_counts": result["new_genome_counts"]}
        elif args.step == "design":
            table = read_block_table(args.block_table) if args.block_table else None
            summary = design_counts(args.peers_tsv, args.peers_dir, args.sequence, table, args.blocks)
        elif args.blocks:
            require(args.peers_tsv is not None, "design mode needs --peers-tsv")
            require(args.block is None and args.target is None and not args.peers,
                    "--blocks cannot be combined with --block/--target/--peers")
            table = read_block_table(args.block_table) if args.block_table else None
            result = make_design_plans(
                plans_10406=args.plans_10406, peers_dir=args.peers_dir, peers_tsv=args.peers_tsv,
                sequence=args.sequence, layouts=args.layouts, models_root=args.models_root,
                out=args.out, context=args.context, gaps=args.gaps,
                perturbations=args.perturbations, cohort_claims_dir=args.cohort_claims_dir,
                blocks=args.blocks, block_table=table,
                regimes=args.regimes or ",".join(DESIGN_REGIMES),
                contexts=args.contexts or ",".join(ALL_CONTEXTS), kinds=args.kinds,
                families=args.families, path_maps=parse_path_maps(args.path_map),
                expect_references=parse_expect_references(args.expect_reference))
            summary = {"status": result["status"], "arms": result["arms"], "runs": result["runs"],
                       "arms_by_family_regime": result["arms_by_family_regime"],
                       "skipped_variants": sum(e["status"] == "skipped"
                                               for e in result["variant_ledger"])}
        else:
            require(args.block and args.target, "plan needs --blocks (design) or --block and --target")
            require(len(args.gaps) == 1, "single-block mode takes exactly one --gaps directory")
            result = make_plans(args.plans_10406, args.peers_dir, args.gaps[0], args.context,
                                args.sequence, args.layouts, args.models_root, args.out,
                                args.block, args.target, args.peers,
                                args.regimes or ",".join(REGIMES),
                                args.contexts or ",".join(CONTEXTS), args.kinds,
                                parse_path_maps(args.path_map), args.cohort_claims_dir)
            summary = {"status": result["status"], "arms": result["arms"],
                       "genomes_per_regime": result["genomes_per_regime"]}
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0
    except (WP4BatchError, ROLE.RoleContextError, FileNotFoundError, KeyError,
            json.JSONDecodeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
