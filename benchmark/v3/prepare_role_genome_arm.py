#!/usr/bin/env python3
"""Prepare sequence, claim, family and layout inputs for the three role genomes.

The three genomes were examined during development and are not strict held-out
controls.  This command reuses the training cohort's DIAMOND and MMseqs2
protocols; it never reads role truth or target functional annotations.

Subcommands are deliberately separate and each writes a JSON receipt.  NCBI
taxonomy responses are cached by URL so a families rerun can be offline.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import shlex
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import acquire_training_cohort as A  # noqa: E402
import homology_baseline as H  # noqa: E402
import prepare_price_cohort as P  # noqa: E402
import prepare_price_layout as PL  # noqa: E402
import replay_predictor_contract as R  # noqa: E402
import select_peer_genomes as S  # noqa: E402

DATA = Path.home() / ".cache" / "gspa-v3-data"
COHORT = DATA / "training-cohort"
PROVENANCE = COHORT / "provenance.json"
SELECTION = HERE / "results" / "label-rich-selection" / "selection_manifest.json"
SYNONYMS = HERE / "registry" / "evaluation_exclusions.json"
COHORT_CONTRACTS = DATA / "predictor-contract-v1" / "contracts.json"
SCHEMA = "gspa.role-genome-arm-inputs/1"
ROLE_GENOMES = {
    "mtb": {"genome_id": "GCF_000195955.2", "taxid": "83332"},
    "ecoli": {"genome_id": "GCF_000005845.2", "taxid": "511145"},
    "bsubtilis": {"genome_id": "GCF_000009045.1", "taxid": "224308"},
}
DEFAULT_EXCLUDE_TAXIDS = ["83332", "1773", "83333", "511145", "224308"]
DEFAULT_EXCLUDE_GENERA = ["mycobacterium", "mycolicibacterium", "mycobacteroides",
                          "escherichia", "shigella", "bacillus"]


class RoleArmError(ValueError):
    pass


def save_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def file_record(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": H.sha256_file(path)}


def output_records(out, names):
    return {name: H.sha256_file(Path(out) / name) for name in names}


def count_lines(path):
    with Path(path).open() as handle:
        return sum(1 for _ in handle)


def expand_genus_synonyms(genera, synonyms):
    """Case-insensitive closure over both keys and values of the registry map."""
    adjacency = defaultdict(set)
    for name, values in synonyms.items():
        name = name.lower()
        for value in values:
            adjacency[name].add(value.lower())
            adjacency[value.lower()].add(name)
    found = {str(g).lower() for g in genera}
    pending = list(found)
    while pending:
        for value in adjacency[pending.pop()] - found:
            found.add(value)
            pending.append(value)
    return found


def build_role_reference(swissprot, provenance, selection, synonyms, out_dir,
                         extra_taxids=DEFAULT_EXCLUDE_TAXIDS,
                         extra_genera=DEFAULT_EXCLUDE_GENERA,
                         provenance_path=None, selection_path=None, synonyms_path=None):
    """Apply the baseline fence, then filter its strict-label output by role taxa/genera.

    The second pass deliberately consumes ``H.build_reference`` output rather
    than reproducing its evidence-tier logic.  Extra-rule counts therefore
    describe entries which survived the original cohort fence and label rule.
    """
    out_dir = Path(out_dir)
    base = H.build_reference(swissprot, provenance, selection, synonyms, out_dir)
    base_sequences = PL.fasta_sequences(out_dir / "reference.faa")
    extra_taxids = {str(t) for t in extra_taxids}
    extra_genera = expand_genus_synonyms(extra_genera, synonyms)
    dropped = {}
    counts = {"dropped_extra_taxid": 0, "dropped_extra_genus": 0}
    for entry in H.iter_entries(swissprot):
        accession = entry["accessions"][0]
        if accession not in base_sequences:
            continue
        if str(entry["taxid"]) in extra_taxids:
            dropped[accession] = "taxid"
            counts["dropped_extra_taxid"] += 1
        elif H.genus_of(entry["organism"]) in extra_genera:
            dropped[accession] = "genus"
            counts["dropped_extra_genus"] += 1
    kept = set(base_sequences) - set(dropped)
    with (out_dir / "reference.faa").open("w") as handle:
        for accession, sequence in base_sequences.items():
            if accession in kept:
                handle.write(f">{accession}\n{sequence}\n")
    label_lines = (out_dir / "reference_labels.tsv").read_text().splitlines()
    with (out_dir / "reference_labels.tsv").open("w") as handle:
        handle.write(label_lines[0] + "\n")
        for line in label_lines[1:]:
            if line.split("\t", 1)[0] in kept:
                handle.write(line + "\n")
    report = dict(base)
    report.update(counts)
    report.update({
        "schema": "gspa.role-genome-reference/1",
        "kept": len(kept),
        "excluded_cohort_taxids": sorted(provenance["organisms"]),
        "excluded_extra_taxids": sorted(extra_taxids, key=int),
        "excluded_extra_genera": sorted(extra_genera),
        "reference_faa_sha256": H.sha256_file(out_dir / "reference.faa"),
        "reference_labels_sha256": H.sha256_file(out_dir / "reference_labels.tsv"),
        "inputs": {"swissprot": file_record(swissprot)},
        "outputs": output_records(out_dir, ("reference.faa", "reference_labels.tsv")),
    })
    for key, path in (("provenance", provenance_path), ("selection", selection_path),
                      ("synonyms", synonyms_path)):
        if path is not None:
            report["inputs"][key] = file_record(path)
    save_json(out_dir / "reference_report.json", report)
    return report


def proteome_record(path):
    path = Path(path).resolve()
    key = path.name[:-4] if path.name.endswith(".faa") else path.stem
    role = ROLE_GENOMES.get(key, {})
    sequences = PL.fasta_sequences(path)
    return {"key": key, "genome_id": role.get("genome_id", key),
            "path": str(path), "sha256": H.sha256_file(path),
            "n_proteins": len(sequences)}, sequences


def combine_proteomes(paths, output):
    records, sequences, owners = [], {}, {}
    for path in paths:
        record, found = proteome_record(path)
        if any(r["genome_id"] == record["genome_id"] for r in records):
            raise RoleArmError(f"duplicate genome {record['genome_id']}")
        for accession, sequence in found.items():
            if accession in owners:
                raise RoleArmError(f"query accession {accession} occurs in both {owners[accession]} and {record['genome_id']}")
            owners[accession] = record["genome_id"]
            sequences[accession] = sequence
        records.append(record)
    with Path(output).open("w") as handle:
        for record in records:
            _, found = proteome_record(record["path"])
            for accession, sequence in found.items():
                handle.write(f">{accession}\n{sequence}\n")
    return records, sequences, owners


def diamond_commands(executable, query, reference, out_dir, threads=4):
    out_dir = Path(out_dir)
    makedb = [executable, "makedb", "--in", str(reference), "-d", str(out_dir / "reference"),
              "--threads", str(threads)]
    blastp = [executable, "blastp", "-q", str(query), "-d", str(out_dir / "reference.dmnd"),
              "-o", str(out_dir / "hits.tsv"), "--threads", str(min(threads, 4)),
              "--evalue", "1e-5", "--max-target-seqs", "10", "--outfmt", "6",
              "qseqid", "sseqid", "pident", "qcovhsp", "evalue", "bitscore", "--sensitive"]
    return makedb, blastp


def search_parameters(command):
    """Parse a role search command with the cohort's strict contract parser."""
    tokens = shlex.split(command)
    if not tokens:
        raise RoleArmError("empty DIAMOND command")
    tokens[0] = "diamond"
    return R.search_parameters(" ".join(tokens))


def run_search(reference_dir, proteomes, out_dir, executable="diamond", threads=4):
    reference_dir, out_dir = Path(reference_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    records, _, _ = combine_proteomes(proteomes, out_dir / "proteomes.faa")
    makedb, blastp = diamond_commands(executable, out_dir / "proteomes.faa",
                                      reference_dir / "reference.faa", out_dir, threads)
    made = subprocess.run(makedb, check=True, capture_output=True, text=True)
    (out_dir / "makedb.log").write_text(made.stdout + made.stderr)
    searched = subprocess.run(blastp, capture_output=True, text=True)
    (out_dir / "diamond.log").write_text(searched.stdout + searched.stderr)
    if searched.returncode:
        raise RoleArmError(f"diamond failed ({searched.returncode}); see {out_dir / 'diamond.log'}")
    version_result = subprocess.run([executable, "version"], check=True, capture_output=True, text=True)
    version = version_result.stdout.strip()
    executable_path = shutil.which(executable)
    if executable_path is None:
        raise RoleArmError(f"DIAMOND executable unavailable: {executable}")
    command = " ".join(blastp)
    report = {
        "schema": "gspa.role-genome-search/1",
        "command": command,
        "makedb_command": " ".join(makedb),
        "parameters": search_parameters(command),
        "version": version,
        "executable": str(Path(executable_path).resolve()),
        "tool_sha256": H.sha256_file(executable_path),
        "threads": threads,
        "proteomes": records,
        "hit_rows": count_lines(out_dir / "hits.tsv"),
        "hits_sha256": H.sha256_file(out_dir / "hits.tsv"),
        "inputs": {"reference": file_record(reference_dir / "reference.faa"),
                   **{r["key"]: {"path": r["path"], "sha256": r["sha256"]} for r in records}},
        "outputs": output_records(out_dir, ("proteomes.faa", "reference.dmnd", "hits.tsv",
                                             "makedb.log", "diamond.log")),
    }
    save_json(out_dir / "search_report.json", report)
    return report


def make_claims(search_dir, reference_dir, out_dir, cohort_contracts=COHORT_CONTRACTS):
    search_dir, reference_dir, out_dir = Path(search_dir), Path(reference_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    search = json.loads((search_dir / "search_report.json").read_text())
    reference_report = json.loads((reference_dir / "reference_report.json").read_text())
    with (out_dir / "proteins.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("genome_id", "protein_id", "kind", "genome_protein_id"))
        for proteome in sorted(search["proteomes"], key=lambda r: r["genome_id"]):
            sequences = PL.fasta_sequences(proteome["path"])
            for accession in sorted(sequences):
                writer.writerow((proteome["genome_id"], accession, "target", accession))
    labels = H.read_reference_labels(reference_dir / "reference_labels.tsv")
    claims = H.claims_from_hits(search_dir / "hits.tsv", labels, "homology")
    files = H.write_claims_jsonl(claims, out_dir / "proteins.tsv", out_dir)
    references = {"sequence_sha256": reference_report["reference_faa_sha256"],
                  "annotation_sha256": reference_report["reference_labels_sha256"]}
    parameters = search_parameters(search["command"])
    identity, contract_parameters = R.contract(
        {"version": search["version"]}, search["tool_sha256"], references,
        parameters, "homology", H.sha256_file(H.__file__))
    for genome in sorted(files):
        path = Path(files[genome]["file"])
        contracted = path.with_name(path.name + ".contracted")
        with path.open() as source, contracted.open("w") as target:
            for line in source:
                row = json.loads(line)
                row["metadata"]["predictor_contract"] = identity
                target.write(json.dumps(row, sort_keys=True) + "\n")
        contracted.replace(path)
        files[genome]["sha256"] = H.sha256_file(path)
        files[genome]["file"] = str(path.resolve())
    cohort_identity = None
    if cohort_contracts is not None:
        cohort_doc = json.loads(Path(cohort_contracts).read_text())
        cohort_identity = cohort_doc["homology"]["identity"]
    differences = sorted(key for key in identity
                         if cohort_identity is None or identity[key] != cohort_identity.get(key))
    report = {
        "schema": "gspa.role-genome-claims/1",
        "mode": "homology",
        "genomes": files,
        "contract": identity,
        "contract_parameters": contract_parameters,
        "reference": references,
        "contract_differs_from_cohort": differences,
        "cohort_contract": {key: cohort_identity.get(key) for key in differences} if cohort_identity else None,
        "inputs": {"hits": file_record(search_dir / "hits.tsv"),
                   "reference_labels": file_record(reference_dir / "reference_labels.tsv"),
                   "search_report": file_record(search_dir / "search_report.json"),
                   "reference_report": file_record(reference_dir / "reference_report.json"),
                   "adapter": file_record(H.__file__)},
        "outputs": {**{g: r["sha256"] for g, r in files.items()},
                    "proteins.tsv": H.sha256_file(out_dir / "proteins.tsv")},
    }
    save_json(out_dir / "claims_report.json", report)
    return report


def mmseqs_command(recorded, executable, fasta, prefix, tmp, threads=4):
    """Retain the recorded easy-cluster argv, changing only executable/paths/threads."""
    tokens = shlex.split(recorded["command"])
    if len(tokens) < 5 or tokens[1] != "easy-cluster":
        raise RoleArmError("cohort report does not record mmseqs easy-cluster")
    parsed = P.parse_mmseqs(None, "recorded", recorded["command"])["parameters"]
    if parsed != recorded.get("parameters"):
        raise RoleArmError("cohort MMseqs2 command and parameter report disagree")
    tokens[0] = executable
    tokens[2:5] = [str(fasta), str(prefix), str(tmp)]
    if "--threads" in tokens:
        tokens[tokens.index("--threads") + 1] = str(threads)
    else:
        tokens.extend(("--threads", str(threads)))
    return tokens


def run_mmseqs(command, out_dir, executable="mmseqs"):
    out_dir = Path(out_dir)
    result = subprocess.run(command, capture_output=True, text=True)
    (out_dir / "mmseqs.log").write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RoleArmError(f"mmseqs failed ({result.returncode}); see {out_dir / 'mmseqs.log'}")
    version = subprocess.run([executable, "version"], check=True, capture_output=True,
                             text=True).stdout.strip()
    protocol = P.parse_mmseqs(None, version, " ".join(command))
    return protocol


def read_gff_contigs(path, genome_id):
    contigs, declared = {}, None
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#!genome-build-accession "):
                declared = line.strip().split(":")[-1]
            elif line.startswith("##sequence-region "):
                fields = line.split()
                if len(fields) != 4:
                    raise RoleArmError(f"invalid sequence-region in {path}: {line.strip()}")
                _, contig, low, high = fields
                if contig in contigs or int(low) != 1 or int(high) < 1:
                    raise RoleArmError(f"invalid or repeated sequence-region in {path}: {line.strip()}")
                contigs[contig] = {"length": int(high), "closed": True}
    if declared is not None and declared != genome_id:
        raise RoleArmError(f"GFF assembly {declared!r} != {genome_id!r}")
    if not contigs:
        raise RoleArmError(f"no sequence-region lines in {path}")
    return contigs


def read_role_layout(path, protein_ids, contigs):
    placed = set()
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["protein_id", "contig", "start", "end", "strand"]:
            raise RoleArmError(f"unexpected layout columns in {path}: {reader.fieldnames}")
        for row in reader:
            protein = row["protein_id"]
            if protein not in protein_ids:
                raise RoleArmError(f"layout protein {protein} is absent from its proteome")
            if protein in placed:
                raise RoleArmError(f"layout repeats protein {protein}")
            if row["contig"] not in contigs:
                raise RoleArmError(f"layout protein {protein} names an undeclared contig")
            placed.add(protein)
    return placed


def parse_assignment(text, parts, kind):
    fields = text.split("=", parts - 1)
    if len(fields) != parts or not all(fields):
        raise RoleArmError(f"invalid --{kind} assignment: {text}")
    return fields


def collect_new_genomes(proteome_args, gff_args, layout_args, taxid_args):
    proteomes = dict((p[0], p[1:]) for p in (parse_assignment(x, 3, "proteome") for x in proteome_args))
    gffs = dict(parse_assignment(x, 2, "gff") for x in gff_args)
    layouts = dict(parse_assignment(x, 2, "layout") for x in layout_args)
    taxids = dict(parse_assignment(x, 2, "taxid") for x in taxid_args)
    if not proteomes or not (set(proteomes) == set(gffs) == set(layouts) == set(taxids)):
        raise RoleArmError("proteome, GFF, layout and taxid keys must match")
    records = []
    for key in sorted(proteomes):
        genome_id, fasta = proteomes[key]
        expected = ROLE_GENOMES.get(key)
        if expected and (genome_id != expected["genome_id"] or taxids[key] != expected["taxid"]):
            raise RoleArmError(f"{key} does not match its pinned genome/taxid")
        records.append({"key": key, "genome_id": genome_id, "taxid": taxids[key],
                        "proteome": str(Path(fasta).resolve()), "gff": str(Path(gffs[key]).resolve()),
                        "layout": str(Path(layouts[key]).resolve())})
    return records


def fetch_role_taxonomy(new_genomes, out_dir, cache_dir, opener=A.default_opener):
    out_dir = Path(out_dir)
    cache = S.ResponseCache(Path(cache_dir), opener)
    taxids = sorted({r["taxid"] for r in new_genomes})
    url = A.NCBI_TAXONOMY_XML.format(ids=",".join(taxids))
    cached = cache.get(url, ".xml")
    A.fetch_taxonomy(taxids, out_dir, lambda _: cached["body"])
    # The raw body is on disk in the cache; the report keeps only its identity.
    return out_dir / "selected_taxonomy.xml", {k: v for k, v in cached.items() if k != "body"}


def write_family_outputs(combined_fasta, cohort_inputs, new_genomes, clusters, taxonomy,
                         out_dir, mmseqs, taxonomy_cache=None):
    """Write cohort-compatible families/metadata using prepare_price_cohort primitives."""
    out_dir, cohort_inputs = Path(out_dir), Path(cohort_inputs)
    metadata_in = json.loads((cohort_inputs / "cohort_metadata.json").read_text())
    sequences = PL.fasta_sequences(combined_fasta)
    proteins = defaultdict(set)
    for namespaced in sequences:
        genome, protein = P.split_namespace(namespaced)
        proteins[genome].add(protein)
    new_ids = {r["genome_id"] for r in new_genomes}
    if set(proteins) != set(metadata_in["genomes"]) | new_ids:
        raise RoleArmError("combined FASTA and cohort/new genome catalogues differ")
    pairs = P.read_clusters(clusters)
    owner, audit = P.audit_membership(pairs, proteins)
    if not audit["exactly_once"]:
        raise RoleArmError(f"MMseqs2 membership is not exact-once: {audit}")
    audit.pop("_missing")
    family = P.family_ids(owner, "fam30")
    genomes = {g: copy.deepcopy(rec) for g, rec in metadata_in["genomes"].items()}
    taxa, absent = P.read_taxonomy(taxonomy, {int(r["taxid"]) for r in new_genomes})
    if absent:
        raise RoleArmError(f"taxonomy is missing requested taxids: {absent}")
    new_report = {}
    for record in sorted(new_genomes, key=lambda r: r["genome_id"]):
        genome, taxid = record["genome_id"], int(record["taxid"])
        info = taxa[taxid]
        if info["genus"] is None or info["species"] is None:
            raise RoleArmError(f"taxonomy {taxid} lacks a unique genus or species rank")
        found = PL.fasta_sequences(record["proteome"])
        if set(found) != proteins[genome]:
            raise RoleArmError(f"combined FASTA differs from proteome for {genome}")
        contigs = read_gff_contigs(record["gff"], genome)
        placed = read_role_layout(record["layout"], set(found), contigs)
        covered = sum(f"{genome}__{p}" in owner for p in found)
        coverage = covered / len(found)
        genomes[genome] = {
            "assemblyStatus": "COMPLETE", "contigs": contigs,
            "unplacedProteins": sorted(set(found) - placed),
            "clade": f"ncbi_genus_{info['genus']}",
            "strainCluster": f"ncbi_species_{info['species']}",
            "familySearchCoverage": round(coverage, 6) if coverage < 1 else 1.0,
            "familySearchComplete": bool(audit["exactly_once"]) and covered == len(found),
        }
        new_report[genome] = {**record, "proteome_sha256": H.sha256_file(record["proteome"]),
                              "proteins": len(found), "in_clusters": covered,
                              "unplaced_proteins": len(set(found) - placed),
                              "taxonomy": info}
    rows = 0
    with (out_dir / "cohort_orthogroups.tsv").open("w") as handle:
        handle.write("\t".join(P.HEADER) + "\n")
        for genome in sorted(proteins):
            for protein in sorted(proteins[genome]):
                namespaced = f"{genome}__{protein}"
                handle.write(f"{genome}\t{protein}\t{family[namespaced]}\t-\t1\n")
                rows += 1
    (out_dir / "cohort_metadata.json").write_text(json.dumps({"genomes": dict(sorted(genomes.items()))}, indent=2) + "\n")
    sizes = Counter(owner.values())
    genomes_per_family, per_genome_family = defaultdict(set), Counter()
    for member, representative in owner.items():
        genome = P.split_namespace(member)[0]
        genomes_per_family[representative].add(genome)
        per_genome_family[(representative, genome)] += 1
    paralogous = {rep for (rep, _), count in per_genome_family.items() if count > 1}
    stats = {"families": len(sizes),
             "singleton_families": sum(1 for n in sizes.values() if n == 1),
             "largest_family": max(sizes.values()),
             "families_in_more_than_one_genome": sum(1 for gs in genomes_per_family.values() if len(gs) > 1),
             "families_with_in_genome_paralogs": len(paralogous)}
    with (out_dir / "families.tsv").open("w") as handle:
        handle.write("orthogroup_id\trepresentative\tproteins\tgenomes\n")
        for representative, size in sorted(sizes.items()):
            handle.write(f"{P.stable_id('fam30', representative)}\t{representative}\t{size}\t"
                         f"{len(genomes_per_family[representative])}\n")
    cohort_representatives = {rep for member, rep in owner.items()
                              if P.split_namespace(member)[0] not in new_ids}
    shared_total = 0
    for genome in sorted(new_ids):
        shared = sum(1 for protein in proteins[genome]
                     if owner[f"{genome}__{protein}"] in cohort_representatives)
        new_report[genome]["proteins_sharing_cohort_family"] = shared
        shared_total += shared
    inputs = {"combined_fasta": file_record(combined_fasta), "clusters": file_record(clusters),
              "taxonomy": file_record(taxonomy),
              "cohort_orthogroups": file_record(cohort_inputs / "cohort_orthogroups.tsv"),
              "cohort_metadata": file_record(cohort_inputs / "cohort_metadata.json"),
              "cohort_report": file_record(cohort_inputs / "report.json")}
    for record in new_genomes:
        for kind in ("proteome", "gff", "layout"):
            inputs[f"{record['key']}_{kind}"] = file_record(record[kind])
    report = {
        "schema": "gspa.role-genome-families/1", "certified_family_search": True,
        "mmseqs": mmseqs, "cluster_audit": audit, "family_stats": stats,
        "new_genome_proteins_sharing_cohort_family": shared_total,
        "new_genomes": new_report,
        "taxonomy": {"requested_taxids": sorted(int(r["taxid"]) for r in new_genomes),
                     "resolved": len(taxa), "raw_response_cache": taxonomy_cache},
        "output_rows": rows, "inputs": inputs,
        "outputs": output_records(out_dir, ("cohort_orthogroups.tsv", "cohort_metadata.json", "families.tsv")),
        "interpretation": P.INTERPRETATION,
    }
    save_json(out_dir / "families_report.json", report)
    return report


def run_families(cohort_fasta, cohort_inputs, new_genomes, out_dir, executable="mmseqs",
                 threads=4, offline_cache=None, opener=A.default_opener):
    out_dir, cohort_inputs = Path(out_dir), Path(cohort_inputs)
    out_dir.mkdir(parents=True, exist_ok=True)
    recorded = json.loads((cohort_inputs / "report.json").read_text())["mmseqs"]
    cohort_sequences = PL.fasta_sequences(cohort_fasta)
    all_sequences = dict(cohort_sequences)
    for record in new_genomes:
        for accession, sequence in PL.fasta_sequences(record["proteome"]).items():
            namespaced = f"{record['genome_id']}__{accession}"
            if namespaced in all_sequences:
                raise RoleArmError(f"duplicate namespaced protein {namespaced}")
            all_sequences[namespaced] = sequence
    with (out_dir / "all_proteins.faa").open("w") as handle:
        for namespaced in sorted(all_sequences):
            handle.write(f">{namespaced}\n{all_sequences[namespaced]}\n")
    work = out_dir / "mmseqs"
    work.mkdir(exist_ok=True)
    command = mmseqs_command(recorded, executable, out_dir / "all_proteins.faa",
                             work / "clusters", work / "tmp", threads)
    protocol = run_mmseqs(command, out_dir, executable)
    clusters = work / "clusters_cluster.tsv"
    cache_dir = offline_cache or out_dir / "ncbi-cache"
    taxonomy, cached = fetch_role_taxonomy(new_genomes, out_dir, cache_dir, opener)
    report = write_family_outputs(out_dir / "all_proteins.faa", cohort_inputs, new_genomes,
                                  clusters, taxonomy, out_dir, protocol, cached)
    report["inputs"]["cohort_fasta"] = file_record(cohort_fasta)
    report["outputs"].update(output_records(out_dir, ("all_proteins.faa", "selected_taxonomy.xml", "mmseqs.log")))
    save_json(out_dir / "families_report.json", report)
    return report


def make_receipt(root):
    root = Path(root)
    reports = {
        "reference": json.loads((root / "reference" / "reference_report.json").read_text()),
        "search": json.loads((root / "search" / "search_report.json").read_text()),
        "claims": json.loads((root / "claims" / "claims_report.json").read_text()),
        "families": json.loads((root / "families" / "families_report.json").read_text()),
    }
    search_by_genome = {r["genome_id"]: r for r in reports["search"]["proteomes"]}
    rows = []
    for genome, family in sorted(reports["families"]["new_genomes"].items(),
                                 key=lambda item: item[1]["key"]):
        proteome, claim = search_by_genome[genome], reports["claims"]["genomes"][genome]
        if proteome["n_proteins"] != family["proteins"]:
            raise RoleArmError(f"protein count differs between search and families for {genome}")
        rows.append({"genome_id": genome, "key": family["key"], "taxid": family["taxid"],
                     "proteome": proteome["path"], "proteome_sha256": proteome["sha256"],
                     "n_proteins": proteome["n_proteins"], "claims_path": claim["file"],
                     "claims_sha256": claim["sha256"], "n_claims": claim["records"],
                     "layout_path": family["layout"], "gff_path": family["gff"]})
    columns = ["genome_id", "key", "taxid", "proteome", "proteome_sha256", "n_proteins",
               "claims_path", "claims_sha256", "n_claims", "layout_path", "gff_path"]
    with (root / "role_genomes.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    receipt = {
        "schema": SCHEMA, "steps": reports,
        "step_receipts": {name: file_record(root / name / f"{name}_report.json") for name in reports},
        "role_genomes": rows,
        "outputs": {"role_genomes.tsv": H.sha256_file(root / "role_genomes.tsv")},
        "limitations": [
            "These role genomes were examined during development; they are not held-out in the strict sense.",
            "The fenced Swiss-Prot reference differs from the cohort reference by the additional role taxid and genus exclusions.",
            "Claims use DIAMOND-only evidence; no CLEAN predictions are supplied for these genomes.",
            f"Tool versions: {reports['search']['version']}; MMseqs2 {reports['families']['mmseqs']['version']}.",
        ],
    }
    save_json(root / "receipt.json", receipt)
    return receipt


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    steps = parser.add_subparsers(dest="step", required=True)
    ref = steps.add_parser("reference")
    ref.add_argument("--swissprot", type=Path, required=True)
    ref.add_argument("--out", type=Path, required=True)
    ref.add_argument("--exclude-taxid", action="append")
    ref.add_argument("--exclude-genus", action="append")
    ref.add_argument("--provenance", type=Path, default=PROVENANCE)
    ref.add_argument("--selection", type=Path, default=SELECTION)
    ref.add_argument("--synonyms", type=Path, default=SYNONYMS)
    search = steps.add_parser("search")
    search.add_argument("--reference", type=Path, required=True)
    search.add_argument("--proteome", type=Path, nargs="+", required=True)
    search.add_argument("--out", type=Path, required=True)
    search.add_argument("--diamond", default="diamond")
    search.add_argument("--threads", type=int, default=4)
    claims = steps.add_parser("claims")
    claims.add_argument("--search", type=Path, required=True)
    claims.add_argument("--reference", type=Path, required=True)
    claims.add_argument("--out", type=Path, required=True)
    claims.add_argument("--cohort-contracts", type=Path, default=COHORT_CONTRACTS)
    families = steps.add_parser("families")
    families.add_argument("--cohort-fasta", type=Path, required=True)
    families.add_argument("--cohort-inputs", type=Path, required=True)
    families.add_argument("--proteome", action="append", required=True, metavar="KEY=GENOME_ID=FAA")
    families.add_argument("--gff", action="append", required=True, metavar="KEY=GFF")
    families.add_argument("--layout", action="append", required=True, metavar="KEY=TSV")
    families.add_argument("--taxid", action="append", required=True, metavar="KEY=TAXID")
    families.add_argument("--out", type=Path, required=True)
    families.add_argument("--mmseqs", default="mmseqs")
    families.add_argument("--threads", type=int, default=4)
    families.add_argument("--offline-cache", type=Path)
    receipt = steps.add_parser("receipt")
    receipt.add_argument("--out", type=Path, required=True)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.step == "reference":
        provenance = json.loads(args.provenance.read_text())
        selection = json.loads(args.selection.read_text())
        synonyms = json.loads(args.synonyms.read_text()).get("genus_synonyms", {})
        result = build_role_reference(
            args.swissprot, provenance, selection, synonyms, args.out,
            args.exclude_taxid or DEFAULT_EXCLUDE_TAXIDS,
            args.exclude_genus or DEFAULT_EXCLUDE_GENERA,
            args.provenance, args.selection, args.synonyms)
        summary = {key: result[key] for key in ("entries_read", "kept", "dropped_cohort_taxid",
                   "dropped_cohort_genus", "dropped_extra_taxid", "dropped_extra_genus", "no_strict_label")}
    elif args.step == "search":
        result = run_search(args.reference, args.proteome, args.out, args.diamond, args.threads)
        summary = {"proteomes": len(result["proteomes"]), "hit_rows": result["hit_rows"],
                   "version": result["version"]}
    elif args.step == "claims":
        result = make_claims(args.search, args.reference, args.out, args.cohort_contracts)
        summary = {"claims": {g: r["records"] for g, r in result["genomes"].items()},
                   "contract_differs_from_cohort": result["contract_differs_from_cohort"]}
    elif args.step == "families":
        new_genomes = collect_new_genomes(args.proteome, args.gff, args.layout, args.taxid)
        result = run_families(args.cohort_fasta, args.cohort_inputs, new_genomes, args.out,
                              args.mmseqs, args.threads, args.offline_cache)
        summary = {**result["family_stats"],
                   "new_genome_proteins_sharing_cohort_family": result["new_genome_proteins_sharing_cohort_family"],
                   "taxonomy": result["taxonomy"]}
    else:
        result = make_receipt(args.out)
        summary = {"schema": result["schema"], "genomes": len(result["role_genomes"])}
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
