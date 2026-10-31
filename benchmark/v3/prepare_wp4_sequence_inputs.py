#!/usr/bin/env python3
"""Prepare label-free WP4 peer claims and cohort-compatible family memberships.

Claims reuse the role-genome arm's DIAMOND search, claim writer and registered
predictor-contract calculation.  Families use a two-stage protocol instead of
the role arm's joint reclustering: peers first search the *original* cohort
representatives, then only proteins without a cohort hit are clustered.  Thus
all cohort orthogroup bytes remain an unchanged prefix of the output.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prepare_price_cohort as P  # noqa: E402
import prepare_price_layout as PL  # noqa: E402
import prepare_role_genome_arm as ARM  # noqa: E402

COHORT_CONTRACTS = ARM.COHORT_CONTRACTS
COHORT_REPORT = (Path.home() / "tmp" / "gspa-role-genomes-v1" / "arm-inputs" /
                 "families" / "families_report.json")
ACQUISITION_FILES = ("protein.faa", "layout.tsv", "contigs.json",
                     "taxonomy.json", "provenance.json")


class SequenceInputError(ValueError):
    pass


def _sha(path):
    return ARM.H.sha256_file(Path(path))


def _record(path):
    return ARM.file_record(path)


def _save(path, value):
    ARM.save_json(Path(path), value)


def _flatten_accessions(values):
    if not values:
        return None
    result = [item for group in values for item in
              (group if isinstance(group, list) else [group])]
    if len(set(result)) != len(result):
        raise SequenceInputError("an accession was requested more than once")
    return result


def select_accessions(peers_dir, requested=None):
    peers_dir = Path(peers_dir)
    requested = _flatten_accessions(requested)
    if requested is None:
        if not peers_dir.is_dir():
            raise SequenceInputError(f"peers directory does not exist: {peers_dir}")
        accessions = sorted(p.name for p in peers_dir.iterdir()
                            if p.is_dir() and (p / "protein.faa").is_file())
    else:
        accessions = sorted(requested)
    if not accessions:
        raise SequenceInputError("no peer accessions selected")
    for accession in accessions:
        if accession in (".", "..") or Path(accession).name != accession:
            raise SequenceInputError(f"invalid peer accession: {accession!r}")
        directory = peers_dir / accession
        if not directory.is_dir() or not (directory / "protein.faa").is_file():
            raise SequenceInputError(f"peer proteome is missing: {directory / 'protein.faa'}")
    return accessions


def _check_bare_fasta(path):
    """Acquisition must already have removed product descriptions from headers."""
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith(">") and (len(line.rstrip("\n\r")) == 1 or
                                          any(c.isspace() for c in line.rstrip("\n\r")[1:])):
                raise SequenceInputError(f"non-bare FASTA header in {path}: {line.strip()!r}")


def _output_map(out, paths):
    out = Path(out).resolve()
    result = {}
    for path in paths:
        path = Path(path).resolve()
        result[path.relative_to(out).as_posix()] = _sha(path)
    return dict(sorted(result.items()))


def enforce_contract_gate(report, expected=None):
    """Permit only an exact contract, or a declared reference-only difference."""
    if "contract_differs_from_cohort" not in report:
        raise SequenceInputError("claims report lacks contract_differs_from_cohort")
    differences = list(report["contract_differs_from_cohort"])
    exact = not differences
    declared_reference = expected == "reference_sha256" and differences == ["reference_sha256"]
    if expected not in (None, "reference_sha256"):
        raise SequenceInputError(f"unsupported expected contract difference: {expected!r}")
    if not (exact or declared_reference):
        expectation = "no difference" if expected is None else "only reference_sha256"
        raise SequenceInputError(
            f"claims predictor contract gate failed: expected {expectation}; "
            f"actual differences={differences}")
    return differences


def run_claims(peers_dir, reference, out, accessions=None, threads=4,
               executable="diamond", cohort_contracts=None,
               expect_contract_difference=None):
    """Search selected peer proteomes and write gated, contracted claim files."""
    if threads < 1:
        raise SequenceInputError("threads must be positive")
    peers_dir, reference, out = Path(peers_dir), Path(reference), Path(out)
    accessions = select_accessions(peers_dir, accessions)
    out.mkdir(parents=True, exist_ok=True)
    staged = out / "proteomes"
    staged.mkdir(exist_ok=True)
    paths = []
    inputs = {"reference_faa": _record(reference / "reference.faa"),
              "reference_labels": _record(reference / "reference_labels.tsv"),
              "reference_report": _record(reference / "reference_report.json")}
    for accession in accessions:
        source = peers_dir / accession / "protein.faa"
        _check_bare_fasta(source)
        target = staged / f"{accession}.faa"
        shutil.copyfile(source, target)
        if _sha(source) != _sha(target):
            raise SequenceInputError(f"staged proteome hash differs for {accession}")
        paths.append(target)
        inputs[f"{accession}_proteome"] = _record(source)

    search = ARM.run_search(reference, paths, out / "search", executable, threads)
    contracts = cohort_contracts if cohort_contracts is not None else COHORT_CONTRACTS
    report = ARM.make_claims(out / "search", reference, out / "claims", contracts)
    enforce_contract_gate(report, expect_contract_difference)
    if set(report["genomes"]) != set(accessions):
        raise SequenceInputError("claim genome set differs from selected peer accessions")

    nested_report = out / "claims" / "claims_report.json"
    report.update({
        "schema": "gspa.wp4-peer-claims/1",
        "accessions": accessions,
        "expected_contract_difference": expect_contract_difference,
        "contract_gate": {"passed": True, "allowed_difference": expect_contract_difference,
                          "actual_differences": report["contract_differs_from_cohort"]},
    })
    report["inputs"].update(inputs)
    output_paths = [out / "claims" / "proteins.tsv", out / "search" / "search_report.json"]
    output_paths.extend(Path(rec["file"]) for rec in report["genomes"].values())
    output_paths.extend(out / "search" / name for name in search["outputs"])
    output_paths.extend(paths)
    report["outputs"] = _output_map(out, output_paths)
    final_report = out / "claims_report.json"
    nested_report.unlink()
    _save(final_report, report)
    return report


def _read_cohort(path):
    """Recover original representatives from their stable fam30 identifiers."""
    path = Path(path)
    raw = path.read_bytes()
    if not raw.endswith(b"\n"):
        raise SequenceInputError("cohort orthogroups must end with a newline")
    membership = {}
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if next(reader, None) != P.HEADER:
            raise SequenceInputError(f"unexpected cohort orthogroup header in {path}")
        for line_number, row in enumerate(reader, 2):
            if len(row) != 5 or not row[0] or not row[1] or not row[2]:
                raise SequenceInputError(f"invalid cohort orthogroup row at {path}:{line_number}")
            key = f"{row[0]}__{row[1]}"
            if key in membership:
                raise SequenceInputError(f"cohort protein occurs twice: {key}")
            membership[key] = row[2]
    if not membership:
        raise SequenceInputError("cohort orthogroups contain no proteins")
    return raw, membership


def _representatives(cohort_fasta, membership):
    sequences = PL.fasta_sequences(cohort_fasta)
    if set(sequences) != set(membership):
        raise SequenceInputError("cohort FASTA and cohort orthogroup proteins differ")
    representatives = defaultdict(list)
    for namespaced, family in membership.items():
        if P.stable_id("fam30", namespaced) == family:
            representatives[family].append(namespaced)
    bad = {family: reps for family, reps in representatives.items() if len(reps) != 1}
    missing = set(membership.values()) - set(representatives)
    if bad or missing:
        raise SequenceInputError("could not recover exactly one original representative per cohort family")
    rep_for_member = {}
    for namespaced, family in membership.items():
        rep_for_member[namespaced] = representatives[family][0]
    return sequences, rep_for_member, {family: reps[0] for family, reps in representatives.items()}


def _recorded_protocol(cohort_report):
    path = Path(cohort_report)
    document = json.loads(path.read_text())
    recorded = document["mmseqs"]
    # This is the same recorded-protocol consistency check used by mmseqs_command.
    ARM.mmseqs_command(recorded, "mmseqs", "old.faa", "old-prefix", "old-tmp", 4)
    parameters = recorded.get("parameters", {})
    if parameters.get("min-seq-id") != "0.3" or parameters.get("c") != "0.8":
        raise SequenceInputError("cohort MMseqs2 protocol is not 30% identity / 80% coverage")
    return recorded


def search_command(recorded, executable, query, target, output, tmp, threads):
    """Build easy-search with the cohort's identity/coverage/sensitivity settings."""
    parsed = P.parse_mmseqs(None, "recorded", recorded["command"])["parameters"]
    if parsed != recorded.get("parameters"):
        raise SequenceInputError("recorded MMseqs2 command and parameters disagree")
    tokens = recorded["command"].split()
    wanted = {"--min-seq-id", "-c", "--cov-mode", "-s", "--max-seqs", "--threads"}
    options, index = [], 5
    while index < len(tokens):
        token = tokens[index]
        if token in wanted:
            if index + 1 >= len(tokens):
                raise SequenceInputError("truncated recorded MMseqs2 option")
            if token != "--threads":
                options.extend((token, tokens[index + 1]))
            index += 2
        else:
            index += 1
    return [str(executable), "easy-search", str(query), str(target), str(output), str(tmp),
            *options, "--threads", str(threads),
            "--format-output", "query,target,fident,qcov,evalue,bits"]


def _write_fasta(path, sequences, order=None):
    with Path(path).open("w") as handle:
        for namespaced in sorted(sequences) if order is None else order:
            handle.write(f">{namespaced}\n{sequences[namespaced]}\n")


def _best_cohort_hits(path, peer_sequences, representative_family):
    best, rows = {}, 0
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6:
                raise SequenceInputError(f"invalid MMseqs2 search row at {path}:{line_number}")
            query, target = fields[:2]
            if query not in peer_sequences or target not in representative_family:
                raise SequenceInputError(f"MMseqs2 search contains a foreign identifier at line {line_number}")
            try:
                identity, coverage, evalue, bitscore = (float(fields[i]) for i in range(2, 6))
            except ValueError as error:
                raise SequenceInputError(f"non-numeric MMseqs2 search row at line {line_number}") from error
            if not all(math.isfinite(x) for x in (identity, coverage, evalue, bitscore)):
                raise SequenceInputError(f"non-finite MMseqs2 search score at line {line_number}")
            identity = identity / 100 if identity > 1 else identity
            coverage = coverage / 100 if coverage > 1 else coverage
            rows += 1
            if identity < 0.3 or coverage < 0.8:
                continue
            rank = (-bitscore, evalue, -identity, -coverage, target)
            if query not in best or rank < best[query][0]:
                best[query] = (rank, target)
    return {query: value[1] for query, value in best.items()}, rows


def _assembly_status(taxonomy):
    explicit = taxonomy.get("assemblyStatus", taxonomy.get("assembly_status"))
    if explicit is not None:
        status = str(explicit).upper()
        if status not in ("COMPLETE", "DRAFT"):
            raise SequenceInputError(f"invalid peer assemblyStatus: {explicit!r}")
        return status
    level = taxonomy.get("assembly_level", taxonomy.get("assemblyLevel"))
    if level is None:
        raise SequenceInputError("peer taxonomy lacks assembly level/status")
    return "COMPLETE" if str(level).strip().lower() in ("complete", "complete genome") else "DRAFT"


def _peer_metadata(accession, directory, proteins):
    taxonomy = json.loads((directory / "taxonomy.json").read_text())
    contigs = json.loads((directory / "contigs.json").read_text())
    if not isinstance(contigs, dict) or not contigs:
        raise SequenceInputError(f"invalid contigs.json for {accession}")
    for contig, record in contigs.items():
        if (not contig or not isinstance(record, dict) or
                type(record.get("length")) is not int or record["length"] < 1 or
                type(record.get("closed")) is not bool):
            raise SequenceInputError(f"invalid contig record for {accession}:{contig}")
    placed = ARM.read_role_layout(directory / "layout.tsv", set(proteins), contigs)
    clade, strain = taxonomy.get("clade"), taxonomy.get("strainCluster")
    if not isinstance(clade, str) or not clade or not isinstance(strain, str) or not strain:
        raise SequenceInputError(f"peer taxonomy lacks clade/strainCluster for {accession}")
    return {
        "assemblyStatus": _assembly_status(taxonomy),
        "contigs": contigs,
        "unplacedProteins": sorted(set(proteins) - placed),
        "clade": clade,
        "strainCluster": strain,
        "familySearchCoverage": 1.0,
        "familySearchComplete": True,
    }


def run_families(peers_dir, cohort_fasta, cohort_orthogroups, cohort_metadata, out,
                 accessions=None, threads=4, executable="mmseqs",
                 cohort_report=None):
    """Assign peers to frozen cohort families and cluster only the remainder."""
    if threads < 1:
        raise SequenceInputError("threads must be positive")
    peers_dir, out = Path(peers_dir), Path(out)
    cohort_fasta = Path(cohort_fasta)
    cohort_orthogroups, cohort_metadata = Path(cohort_orthogroups), Path(cohort_metadata)
    cohort_report = Path(cohort_report) if cohort_report is not None else COHORT_REPORT
    accessions = select_accessions(peers_dir, accessions)
    out.mkdir(parents=True, exist_ok=True)
    recorded = _recorded_protocol(cohort_report)
    cohort_bytes, cohort_membership = _read_cohort(cohort_orthogroups)
    cohort_hash = hashlib.sha256(cohort_bytes).hexdigest()
    if cohort_hash != _sha(cohort_orthogroups):
        raise SequenceInputError("cohort orthogroups changed while being read")
    cohort_sequences, cohort_owner, representative_family = _representatives(
        cohort_fasta, cohort_membership)

    metadata_document = json.loads(cohort_metadata.read_text())
    if not isinstance(metadata_document.get("genomes"), dict):
        raise SequenceInputError("cohort metadata has no genomes object")
    genomes = metadata_document["genomes"]
    cohort_genomes = {P.split_namespace(key)[0] for key in cohort_sequences}
    if set(genomes) != cohort_genomes:
        raise SequenceInputError("cohort metadata and FASTA genome sets differ")
    if set(genomes) & set(accessions):
        raise SequenceInputError("a selected peer is already present in cohort metadata")

    peer_sequences, proteins_by_genome = {}, {}
    inputs = {
        "cohort_fasta": _record(cohort_fasta),
        "cohort_orthogroups": _record(cohort_orthogroups),
        "cohort_metadata": _record(cohort_metadata),
        "cohort_report": _record(cohort_report),
    }
    peer_metadata = {}
    for accession in accessions:
        directory = peers_dir / accession
        for name in ACQUISITION_FILES:
            if not (directory / name).is_file():
                raise SequenceInputError(f"peer acquisition file is missing: {directory / name}")
            inputs[f"{accession}_{name}"] = _record(directory / name)
        fasta = directory / "protein.faa"
        _check_bare_fasta(fasta)
        found = PL.fasta_sequences(fasta)
        proteins_by_genome[accession] = set(found)
        peer_metadata[accession] = _peer_metadata(accession, directory, found)
        for protein, sequence in found.items():
            namespaced = f"{accession}__{protein}"
            if namespaced in cohort_sequences or namespaced in peer_sequences:
                raise SequenceInputError(f"duplicate namespaced peer protein: {namespaced}")
            peer_sequences[namespaced] = sequence

    representative_sequences = {rep: cohort_sequences[rep]
                                for rep in representative_family.values()}
    peer_fasta = out / "peer_proteins.faa"
    representative_fasta = out / "cohort_representatives.faa"
    _write_fasta(peer_fasta, peer_sequences)
    _write_fasta(representative_fasta, representative_sequences)
    search_dir = out / "search"
    search_dir.mkdir(exist_ok=True)
    hits_path = search_dir / "hits.tsv"
    command = search_command(recorded, executable, peer_fasta, representative_fasta,
                             hits_path, search_dir / "tmp", threads)
    search_protocol = ARM.run_mmseqs(command, search_dir, executable)
    assigned, hit_rows = _best_cohort_hits(
        hits_path, peer_sequences, set(representative_family.values()))
    unplaced = sorted(set(peer_sequences) - set(assigned))
    unassigned_fasta = out / "unassigned.faa"
    _write_fasta(unassigned_fasta, {key: peer_sequences[key] for key in unplaced}, unplaced)

    cluster_dir = out / "cluster"
    cluster_dir.mkdir(exist_ok=True)
    new_pairs, cluster_protocol = [], None
    cluster_prefix = cluster_dir / "clusters"
    cluster_path = cluster_dir / "clusters_cluster.tsv"
    if unplaced:
        cluster_command = ARM.mmseqs_command(recorded, executable, unassigned_fasta,
                                             cluster_prefix, cluster_dir / "tmp", threads)
        cluster_protocol = ARM.run_mmseqs(cluster_command, cluster_dir, executable)
        try:
            new_pairs = P.read_clusters(cluster_path)
        except P.CohortError as error:
            raise SequenceInputError(f"invalid unassigned-protein clusters: {error}") from error
    else:
        cluster_path.write_text("")

    all_proteins = {genome: set() for genome in cohort_genomes}
    for namespaced in cohort_sequences:
        genome, protein = P.split_namespace(namespaced)
        all_proteins[genome].add(protein)
    all_proteins.update(proteins_by_genome)
    pairs = [(cohort_owner[member], member) for member in cohort_sequences]
    pairs.extend((target, query) for query, target in assigned.items())
    pairs.extend(new_pairs)
    owner, audit = P.audit_membership(pairs, all_proteins)
    if not audit["exactly_once"]:
        raise SequenceInputError(f"peer/cohort family membership is not exact-once: {audit}")
    audit.pop("_missing")
    family = P.family_ids(owner, "fam30")
    for member, expected in cohort_membership.items():
        if family[member] != expected:
            raise SequenceInputError(f"cohort family changed for {member}")

    output_orthogroups = out / "cohort_orthogroups.tsv"
    peer_reports = {}
    with output_orthogroups.open("wb") as handle:
        handle.write(cohort_bytes)
        for accession in accessions:
            assigned_count, new_families, new_protein_count = 0, set(), 0
            rows = []
            for protein in sorted(proteins_by_genome[accession]):
                namespaced = f"{accession}__{protein}"
                family_id = family[namespaced]
                rows.append(f"{accession}\t{protein}\t{family_id}\t-\t1\n".encode())
                if namespaced in assigned:
                    assigned_count += 1
                else:
                    new_families.add(family_id)
                    new_protein_count += 1
            for row in rows:
                handle.write(row)
            peer_reports[accession] = {
                "proteins": len(proteins_by_genome[accession]),
                "assigned_to_cohort_clusters": assigned_count,
                "new_clusters": len(new_families),
                "proteins_in_new_clusters": new_protein_count,
                "unassigned": 0,
                "unplacedProteins": len(peer_metadata[accession]["unplacedProteins"]),
            }
    written = output_orthogroups.read_bytes()
    if (hashlib.sha256(written[:len(cohort_bytes)]).hexdigest() != cohort_hash or
            written[:len(cohort_bytes)] != cohort_bytes):
        raise SequenceInputError("cohort orthogroup byte-prefix verification failed")

    extended = {genome: dict(record) for genome, record in genomes.items()}
    extended.update(peer_metadata)
    output_metadata = out / "cohort_metadata.json"
    output_metadata.write_text(json.dumps(
        {"genomes": {key: extended[key] for key in sorted(extended)}}, indent=2) + "\n")

    output_paths = [output_orthogroups, output_metadata, peer_fasta,
                    representative_fasta, unassigned_fasta, hits_path,
                    search_dir / "mmseqs.log", cluster_path]
    if cluster_protocol is not None:
        output_paths.append(cluster_dir / "mmseqs.log")
    report = {
        "schema": "gspa.wp4-peer-families/1",
        "certified_family_search": True,
        "accessions": accessions,
        "peers": peer_reports,
        "proteins": len(peer_sequences),
        "assigned_to_cohort_clusters": len(assigned),
        "new_clusters": len({owner[key] for key in unplaced}),
        "unassigned": 0,
        "search_hit_rows": hit_rows,
        "cluster_audit": audit,
        "mmseqs": {"recorded": recorded, "search": search_protocol,
                   "unassigned_cluster": cluster_protocol},
        "cohort_preservation": {"sha256": cohort_hash, "bytes": len(cohort_bytes),
                                "byte_identical_prefix": True},
        "inputs": inputs,
        "outputs": _output_map(out, output_paths),
        "interpretation": P.INTERPRETATION +
                          " Peer-only extension uses assignment to frozen representatives, not joint reclustering.",
    }
    _save(out / "families_report.json", report)
    return report


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    steps = parser.add_subparsers(dest="step", required=True)
    claims = steps.add_parser("claims")
    claims.add_argument("--peers-dir", type=Path, required=True)
    claims.add_argument("--reference", type=Path, required=True)
    claims.add_argument("--out", type=Path, required=True)
    claims.add_argument("--accession", action="append", nargs="+")
    claims.add_argument("--threads", type=int, default=4)
    claims.add_argument("--diamond", default="diamond")
    claims.add_argument("--expect-contract-difference", choices=["reference_sha256"])
    families = steps.add_parser("families")
    families.add_argument("--peers-dir", type=Path, required=True)
    families.add_argument("--cohort-fasta", type=Path, required=True)
    families.add_argument("--cohort-orthogroups", type=Path, required=True)
    families.add_argument("--cohort-metadata", type=Path, required=True)
    families.add_argument("--out", type=Path, required=True)
    families.add_argument("--accession", action="append", nargs="+")
    families.add_argument("--threads", type=int, default=4)
    families.add_argument("--mmseqs", default="mmseqs")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.step == "claims":
        result = run_claims(args.peers_dir, args.reference, args.out, args.accession,
                            args.threads, args.diamond,
                            expect_contract_difference=args.expect_contract_difference)
        summary = {"accessions": result["accessions"],
                   "contract_differs_from_cohort": result["contract_differs_from_cohort"]}
    else:
        result = run_families(args.peers_dir, args.cohort_fasta, args.cohort_orthogroups,
                              args.cohort_metadata, args.out, args.accession,
                              args.threads, args.mmseqs)
        summary = {"accessions": result["accessions"], "proteins": result["proteins"],
                   "assigned_to_cohort_clusters": result["assigned_to_cohort_clusters"],
                   "new_clusters": result["new_clusters"], "unassigned": result["unassigned"]}
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
