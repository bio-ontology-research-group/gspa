#!/usr/bin/env python3
"""Acquire label-free, exact-version NCBI packages for the WP4 peer genomes.

``acquire`` reads the frozen peer table/selection, excludes cohort assemblies and
GCA/GCF twins, and writes one directory per unique accession.  The original NCBI
package members remain below ``package/``; the derived FASTA at the accession
root has genome-unique ``<accession>__<raw_id>`` protein identifiers so product
names cannot leak into predictors.  Coordinates come only from GFF CDS
``protein_id`` attributes.  Proteins without
one usable locus are retained in the FASTA and listed in ``unplaced.txt``.

NCBI Datasets report pages and Taxonomy batches recorded by
``select_peer_genomes.py`` are reused through their raw-response cache.  Only a
missing/stale genome package normally requires a download.  ``receipt`` performs
no network access and verifies all requested output hashes.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import acquire_training_cohort as A  # noqa: E402
import prepare_price_layout as PL  # noqa: E402
import select_peer_genomes as S  # noqa: E402

SCHEMA = "gspa.wp4-peer-genome-acquisition/1"
PACKAGE_SCHEMA = "gspa.wp4-peer-genome-package/1"
RECEIPT_SCHEMA = "gspa.wp4-peer-genome-receipt/1"
LAYOUT_FIELDS = ["protein_id", "contig", "start", "end", "strand"]
PACKAGE_FILES = ("protein.faa", "genomic.gff", "genomic.fna", "assembly_data_report.json")
OUTPUT_FILES = ("protein.faa", "genomic.gff", "layout.tsv", "unplaced.txt",
                "contigs.json", "taxonomy.json")
PROTEIN_ID_SCHEME = "accession__raw"


class AcquisitionError(ValueError):
    pass


def save_json(path: Path, value: dict) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def file_record(path: Path) -> dict:
    return {"path": str(Path(path).resolve()), "sha256": A.sha256_file(path)}


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text())


def load_peer_inputs(peers_path: Path, selection_path: Path) -> tuple:
    """Return unique accession order and their target/role rows, without labels."""
    selection = read_json(selection_path)
    if selection.get("schema") != S.SCHEMA:
        raise AcquisitionError(f"unexpected peer selection schema {selection.get('schema')!r}")
    details = {}
    for target in selection.get("targets", []):
        for peer in target.get("peers", []):
            details[(str(target.get("taxid")), peer.get("accession"))] = peer
    with Path(peers_path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        if not {"target_taxid", "role"} <= fields or not fields & {"accession", "peer_accession"}:
            raise AcquisitionError("peers.tsv lacks target_taxid, accession/peer_accession or role")
        rows = []
        for line_number, row in enumerate(reader, 2):
            accession = row.get("accession") or row.get("peer_accession")
            if not accession:
                raise AcquisitionError(f"peers.tsv:{line_number}: empty accession")
            key = (str(row["target_taxid"]), accession)
            detail = details.get(key, {})
            rows.append({
                "target_taxid": str(row["target_taxid"]),
                "target_accession": row.get("target_accession") or None,
                "accession": accession,
                "role": row["role"],
                "peer_species_taxid": (row.get("peer_species_taxid") or
                                       detail.get("species_taxid")),
                "peer_organism": row.get("peer_organism") or detail.get("organism"),
                "rank_reasons": detail.get("rank_reasons") or (
                    json.loads(row["rank_reasons"]) if row.get("rank_reasons", "").startswith("[")
                    else ([row["rank_reasons"]] if row.get("rank_reasons") else [])),
            })
    if not rows:
        raise AcquisitionError("peers.tsv contains no peers")
    by_accession = {}
    for row in rows:
        by_accession.setdefault(row["accession"], []).append(row)
    for accession, records in by_accession.items():
        records.sort(key=lambda r: (r["target_taxid"], r["target_accession"] or "", r["role"]))
    return sorted(by_accession), by_accession, selection


def choose_only(accessions: list, only) -> list:
    if not only:
        return list(accessions)
    requested = sorted(set(only))
    unknown = sorted(set(requested) - set(accessions))
    if unknown:
        raise AcquisitionError("--only names accessions absent from peers.tsv: " + ", ".join(unknown))
    return requested


def resolve_input_path(value: str, selection_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute() or path.exists():
        return path
    beside = selection_path.resolve().parent / path
    return beside if beside.exists() else path


def cohort_accessions(selection: dict, selection_path: Path) -> set:
    """Include used, declared and paired cohort assemblies in the twin fence."""
    source = (selection.get("inputs") or {}).get("provenance")
    if not source or not source.get("path"):
        raise AcquisitionError("peer_selection.json does not identify cohort provenance")
    path = resolve_input_path(source["path"], Path(selection_path))
    if not path.is_file():
        raise AcquisitionError(f"cohort provenance is missing: {path}")
    if source.get("sha256") and A.sha256_file(path) != source["sha256"]:
        raise AcquisitionError(f"cohort provenance sha256 differs from peer selection: {path}")
    provenance = read_json(path)
    found = set()
    for record in provenance.get("organisms", {}).values():
        assembly = record.get("assembly") or {}
        for key in ("used_accession", "declared_accession", "paired_accession"):
            if assembly.get(key):
                found.add(assembly[key])
    if not found:
        raise AcquisitionError("cohort provenance contains no assembly accessions")
    return found


def cached_source(record: dict, path: Path) -> dict:
    return {key: record.get(key) for key in ("url", "url_sha256", "sha256", "size_bytes")
            if record.get(key) is not None} | {"cache_file": str(path.resolve())}


class SelectionCache:
    """Read the exact raw responses recorded in peer_selection.json."""

    def __init__(self, selection: dict, cache_dir: Path, opener=A.default_opener):
        self.selection = selection
        self.cache = S.ResponseCache(Path(cache_dir), opener)

    def body(self, record: dict, suffix: str) -> tuple:
        url = record.get("url")
        if not url:
            raise AcquisitionError("cached NCBI record has no URL")
        expected = self.cache.path_for(url, suffix)
        original = Path(record.get("cache_file") or "")
        path = expected if expected.exists() else original
        if path.is_file():
            body = path.read_bytes()
        else:
            fetched = self.cache.get(url, suffix)
            body, path = fetched["body"], Path(fetched["cache_file"])
        if record.get("sha256") and sha256_bytes(body) != record["sha256"]:
            raise AcquisitionError(f"cached response sha256 differs for {url}")
        if record.get("size_bytes") is not None and len(body) != record["size_bytes"]:
            raise AcquisitionError(f"cached response size differs for {url}")
        return body, cached_source(record, path)

    def reports(self, accessions: list) -> tuple:
        wanted = set(accessions)
        found, sources = {}, {}
        for page in self.selection.get("ncbi_pages", []):
            body, source = self.body(page, ".json")
            document = json.loads(body)
            for report in document.get("reports") or []:
                accession = report.get("accession")
                if accession not in wanted:
                    continue
                if accession in found and found[accession] != report:
                    raise AcquisitionError(f"conflicting cached dataset reports for {accession}")
                found[accession], sources[accession] = report, source
        for accession in sorted(wanted - set(found)):
            url = A.NCBI_REPORT.format(acc=accession)
            fetched = self.cache.get(url, ".json")
            reports = json.loads(fetched["body"]).get("reports") or []
            if len(reports) != 1 or reports[0].get("accession") != accession:
                raise AcquisitionError(f"{accession}: dataset_report did not return exactly that assembly")
            found[accession] = reports[0]
            sources[accession] = {key: fetched[key] for key in
                                  ("url", "url_sha256", "sha256", "size_bytes", "cache_file")}
        return found, sources

    def taxonomy(self) -> tuple:
        classifications, sources = {}, []
        batches = self.selection.get("ncbi_taxonomy_batches") or []
        if not batches:
            raise AcquisitionError("peer_selection.json has no cached taxonomy batches")
        for batch in batches:
            body, source = self.body(batch, ".xml")
            root = ET.fromstring(body)
            for taxon in root.findall("Taxon"):
                own_id = taxon.findtext("TaxId")
                if not own_id:
                    raise AcquisitionError("NCBI Taxonomy record lacks TaxId")
                ranks = {"species": [], "genus": []}
                for parent in taxon.findall("LineageEx/Taxon"):
                    rank = parent.findtext("Rank")
                    if rank in ranks:
                        ranks[rank].append((parent.findtext("TaxId"), parent.findtext("ScientificName")))
                own_rank = taxon.findtext("Rank")
                own = (own_id, taxon.findtext("ScientificName"))
                if own_rank in ranks and all(entry[0] != own_id for entry in ranks[own_rank]):
                    ranks[own_rank].append(own)
                resolved = {}
                for rank, entries in ranks.items():
                    unique = list(dict.fromkeys(entries))
                    if len(unique) == 1 and all(unique[0]):
                        resolved[rank] = {"taxid": str(unique[0][0]), "name": unique[0][1]}
                if own_id in classifications and classifications[own_id] != resolved:
                    raise AcquisitionError(f"conflicting cached taxonomy for {own_id}")
                classifications[str(own_id)] = resolved
            sources.append(source)
        return classifications, sources


def nested(record: dict, *names: str) -> dict:
    for name in names:
        value = record.get(name)
        if isinstance(value, dict):
            return value
    return {}


def assembly_info(report: dict) -> dict:
    return nested(report, "assembly_info", "assemblyInfo")


def assembly_level(report: dict) -> str:
    info = assembly_info(report)
    return info.get("assembly_level") or info.get("assemblyLevel") or ""


def paired_accession(report: dict):
    pair = report.get("paired_accession") or report.get("pairedAccession")
    if pair:
        return pair
    info = assembly_info(report)
    paired = nested(info, "paired_assembly", "pairedAssembly")
    return paired.get("accession")


def organism_taxid(report: dict) -> str:
    organism = nested(report, "organism")
    taxid = organism.get("tax_id") if "tax_id" in organism else organism.get("taxId")
    if taxid is None:
        raise AcquisitionError(f"{report.get('accession')}: assembly report has no organism taxid")
    return str(taxid)


def organism_name(report: dict) -> str:
    organism = nested(report, "organism")
    return organism.get("organism_name") or organism.get("organismName") or ""


def complete_assembly(report: dict) -> bool:
    return assembly_level(report).strip().lower().replace("_", " ") in {
        "complete", "complete genome"}


def already_acquired(out_dir: Path, accession: str) -> bool:
    directory = Path(out_dir) / accession
    # A provenance file is an acquisition marker even when later verification
    # finds damage; refusing the twin is safer than silently replacing its role.
    return ((directory / "provenance.json").exists() or
            all((directory / name).is_file() for name in OUTPUT_FILES))


def enforce_exclusions(accession: str, paired: str, cohort: set, out_dir: Path) -> None:
    if accession in cohort:
        raise AcquisitionError(f"{accession}: accession is a cohort assembly")
    if paired:
        if paired == accession:
            raise AcquisitionError(f"{accession}: assembly reports itself as its paired accession")
        if paired in cohort:
            raise AcquisitionError(f"{accession}: paired accession {paired} is a cohort assembly")
        if already_acquired(out_dir, paired):
            raise AcquisitionError(f"{accession}: paired accession {paired} is already acquired")


def package_valid(package_dir: Path, accession: str) -> bool:
    source_path = package_dir / "package.json"
    if not source_path.is_file():
        return False
    try:
        source = read_json(source_path)
    except (OSError, json.JSONDecodeError):
        return False
    if (source.get("schema") != PACKAGE_SCHEMA or source.get("accession") != accession or
            not (source.get("package") or {}).get("sha256")):
        return False
    hashes = source.get("member_sha256s") or {}
    for name in PACKAGE_FILES:
        path = package_dir / name
        if not hashes.get(name) or not path.is_file() or A.sha256_file(path) != hashes[name]:
            return False
    return True


def ensure_package(accession: str, out_dir: Path, opener=A.default_opener) -> dict:
    package_dir = Path(out_dir) / accession / "package"
    source_path = package_dir / "package.json"
    if package_valid(package_dir, accession):
        return read_json(source_path)
    fetched = A.fetch(A.NCBI_PACKAGE.format(acc=accession), opener=opener)
    package_dir.mkdir(parents=True, exist_ok=True)
    members = A.extract_package(fetched["body"], accession, package_dir)
    missing = [name for name in PACKAGE_FILES if not members.get(name)]
    if missing:
        raise AcquisitionError(f"{accession}: NCBI package lacks " + ", ".join(missing))
    source = {"schema": PACKAGE_SCHEMA, "accession": accession,
              "package": A.provenance(fetched), "member_sha256s": members}
    save_json(source_path, source)
    if not package_valid(package_dir, accession):
        raise AcquisitionError(f"{accession}: extracted package failed hash verification")
    return source


def audit_cds_joins(gff_path: Path, protein_ids: set) -> dict:
    """Reject CDS protein joins absent from FASTA; count unjoinable CDS rows."""
    total, without, unknown = 0, 0, set()
    with Path(gff_path).open() as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise AcquisitionError(f"{gff_path}:{line_number}: expected 9 GFF fields")
            if fields[2] != "CDS":
                continue
            total += 1
            protein_id = PL.attributes(fields[8]).get("protein_id")
            if not protein_id:
                without += 1
            elif protein_id not in protein_ids:
                unknown.add(protein_id)
    if unknown:
        raise AcquisitionError(f"{gff_path}: CDS protein_id absent from FASTA: " +
                               ", ".join(sorted(unknown)))
    return {"gff_cds": total, "gff_cds_without_protein_id": without}


def qualified_protein_id(accession: str, raw_id: str) -> str:
    return f"{accession}__{raw_id}"


def write_genome_unique_fasta(source: Path, destination: Path, accession: str) -> dict:
    sequences = PL.fasta_sequences(source)
    with Path(destination).open("w") as handle:
        for raw_id in sorted(sequences):
            protein_id = qualified_protein_id(accession, raw_id)
            handle.write(f">{protein_id}\n{sequences[raw_id]}\n")
    # Reparse the actual output: all derived files use exactly these IDs.
    return PL.fasta_sequences(destination)


def write_layout(gff_path: Path, accession: str, protein_ids: set,
                 layout_path: Path, unplaced_path: Path) -> tuple:
    # GFF joins use the original package IDs; only the derived rows are qualified.
    raw_rows, unavailable, contigs = PL.read_coordinates(gff_path, accession, protein_ids)
    if not contigs:
        raise AcquisitionError(f"{accession}: GFF declares no sequence regions")
    raw_unplaced = sorted(protein_ids - {row[0] for row in raw_rows})
    if raw_unplaced != sorted(unavailable):
        raise AcquisitionError(f"{accession}: unplaced set is not FASTA minus layout")
    rows = [(qualified_protein_id(accession, row[0]), *row[1:]) for row in raw_rows]
    unplaced = [qualified_protein_id(accession, raw_id) for raw_id in raw_unplaced]
    with Path(layout_path).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(LAYOUT_FIELDS)
        writer.writerows(rows)
    Path(unplaced_path).write_text("".join(protein + "\n" for protein in unplaced))
    return rows, unplaced, contigs


def build_outputs(accession: str, report: dict, taxonomy: dict, taxonomy_sources: list,
                  report_source: dict, package: dict, selection_rows: list,
                  out_dir: Path) -> dict:
    directory = Path(out_dir) / accession
    package_dir = directory / "package"
    directory.mkdir(parents=True, exist_ok=True)
    # Root GFF is a byte-for-byte copy; only the derived FASTA is rewritten.
    shutil.copyfile(package_dir / "genomic.gff", directory / "genomic.gff")
    raw_sequences = PL.fasta_sequences(package_dir / "protein.faa")
    sequences = write_genome_unique_fasta(
        package_dir / "protein.faa", directory / "protein.faa", accession)
    joins = audit_cds_joins(directory / "genomic.gff", set(raw_sequences))
    rows, unplaced, contigs = write_layout(
        directory / "genomic.gff", accession, set(raw_sequences),
        directory / "layout.tsv", directory / "unplaced.txt")
    closed = complete_assembly(report)
    contig_doc = {contig: {"length": values["length"], "closed": closed}
                  for contig, values in sorted(contigs.items())}
    save_json(directory / "contigs.json", contig_doc)
    taxid = organism_taxid(report)
    classification = taxonomy.get(taxid)
    if not classification or not classification.get("species") or not classification.get("genus"):
        raise AcquisitionError(f"{accession}: cached taxonomy lacks species/genus for organism {taxid}")
    expected_species = {str(row["peer_species_taxid"]) for row in selection_rows
                        if row.get("peer_species_taxid")}
    species = classification["species"]
    if expected_species and expected_species != {species["taxid"]}:
        raise AcquisitionError(f"{accession}: taxonomy species {species['taxid']} differs from selection " +
                               "/".join(sorted(expected_species)))
    genus = classification["genus"]
    taxonomy_doc = {
        "accession": accession,
        "organism_taxid": taxid,
        "organism_name": organism_name(report),
        "species_taxid": species["taxid"], "species_name": species["name"],
        "genus_taxid": genus["taxid"], "genus_name": genus["name"],
        "species": species, "genus": genus,
        "assembly_level": assembly_level(report),
        "clade": f"ncbi_genus_{genus['taxid']}",
        "strainCluster": f"ncbi_species_{species['taxid']}",
    }
    save_json(directory / "taxonomy.json", taxonomy_doc)
    counts = {"proteins": len(sequences), "layout_rows": len(rows),
              "unplaced": len(unplaced), "contigs": len(contig_doc), **joins}
    provenance = {
        "schema": SCHEMA, "status": "acquired", "accession": accession,
        "protein_id_scheme": PROTEIN_ID_SCHEME,
        "paired_accession": paired_accession(report),
        "selection_rows": selection_rows,
        "urls": {"package": package["package"]["url"],
                 "dataset_report": report_source.get("url"),
                 "taxonomy_batches": [source.get("url") for source in taxonomy_sources]},
        "package_sha256": package["package"]["sha256"],
        "package_source": package["package"],
        "package_member_sha256s": package["member_sha256s"],
        "dataset_report_source": report_source,
        "taxonomy_sources": taxonomy_sources,
        "file_sha256s": {name: A.sha256_file(directory / name) for name in OUTPUT_FILES},
        "counts": counts,
    }
    # This is deliberately written last: its presence denotes completed outputs.
    save_json(directory / "provenance.json", provenance)
    return provenance


def verify_protein_ids(directory: Path, accession: str, provenance: dict,
                       problems: list) -> None:
    """Check the shared derived ID space independently of recorded hashes."""
    if provenance and provenance.get("protein_id_scheme") != PROTEIN_ID_SCHEME:
        problems.append("invalid:protein_id_scheme")

    fasta_ids = None
    fasta_path = directory / "protein.faa"
    if fasta_path.is_file():
        try:
            fasta_ids = set(PL.fasta_sequences(fasta_path))
        except (OSError, ValueError) as error:
            problems.append(f"invalid:protein.faa:{error}")
        else:
            prefix = accession + "__"
            if any(not protein_id.startswith(prefix) for protein_id in fasta_ids):
                problems.append("protein_id_scheme:protein.faa")
            recorded_count = ((provenance or {}).get("counts") or {}).get("proteins")
            if recorded_count is not None and recorded_count != len(fasta_ids):
                problems.append("protein_id_count:protein.faa")

    layout_ids = None
    layout_path = directory / "layout.tsv"
    if layout_path.is_file():
        try:
            with layout_path.open(newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames != LAYOUT_FIELDS:
                    raise ValueError("unexpected layout columns")
                layout_rows = list(reader)
            layout_ids = [row["protein_id"] for row in layout_rows]
            if any(not protein_id for protein_id in layout_ids):
                raise ValueError("empty protein_id")
            if len(layout_ids) != len(set(layout_ids)):
                raise ValueError("duplicate protein_id")
            layout_ids = set(layout_ids)
        except (OSError, ValueError, csv.Error) as error:
            problems.append(f"invalid:layout.tsv:{error}")
            layout_ids = None
        if layout_ids is not None and fasta_ids is not None and not layout_ids <= fasta_ids:
            problems.append("protein_id_mismatch:layout.tsv")

    unplaced_ids = None
    unplaced_path = directory / "unplaced.txt"
    if unplaced_path.is_file():
        try:
            unplaced_rows = unplaced_path.read_text().splitlines()
            if any(not protein_id for protein_id in unplaced_rows):
                raise ValueError("empty protein_id")
            if len(unplaced_rows) != len(set(unplaced_rows)):
                raise ValueError("duplicate protein_id")
            unplaced_ids = set(unplaced_rows)
        except (OSError, ValueError) as error:
            problems.append(f"invalid:unplaced.txt:{error}")
    if fasta_ids is not None and layout_ids is not None and unplaced_ids is not None:
        if layout_ids & unplaced_ids or layout_ids | unplaced_ids != fasta_ids:
            problems.append("protein_id_partition:layout.tsv+unplaced.txt")


def verify_genome(out_dir: Path, accession: str) -> dict:
    directory = Path(out_dir) / accession
    problems = []
    provenance_path = directory / "provenance.json"
    provenance = None
    if not provenance_path.is_file():
        problems.append("missing:provenance.json")
    else:
        try:
            provenance = read_json(provenance_path)
        except (OSError, json.JSONDecodeError) as error:
            problems.append(f"invalid:provenance.json:{error}")
    if provenance and (provenance.get("schema") != SCHEMA or
                       provenance.get("accession") != accession):
        problems.append("invalid:provenance.json")
    hashes = (provenance or {}).get("file_sha256s") or {}
    for name in OUTPUT_FILES:
        path = directory / name
        if not path.is_file():
            problems.append("missing:" + name)
        elif not hashes.get(name):
            problems.append("unrecorded:" + name)
        elif A.sha256_file(path) != hashes[name]:
            problems.append("sha256:" + name)
    package_hashes = (provenance or {}).get("package_member_sha256s") or {}
    if not (directory / "package" / "package.json").is_file():
        problems.append("missing:package/package.json")
    for name in PACKAGE_FILES:
        path = directory / "package" / name
        if not path.is_file():
            problems.append("missing:package/" + name)
        elif not package_hashes.get(name):
            problems.append("unrecorded:package/" + name)
        elif A.sha256_file(path) != package_hashes[name]:
            problems.append("sha256:package/" + name)
    verify_protein_ids(directory, accession, provenance, problems)
    record = {"accession": accession, "directory": str(directory.resolve()),
              "status": "verified" if not problems else "incomplete",
              "problems": problems}
    if provenance:
        record.update({"counts": provenance.get("counts"),
                       "package_sha256": provenance.get("package_sha256"),
                       "provenance_sha256": A.sha256_file(provenance_path)})
    return record


def make_receipt(peers_path: Path, selection_path: Path, out_dir: Path, only=None) -> dict:
    accessions, _, _ = load_peer_inputs(peers_path, selection_path)
    requested = choose_only(accessions, only)
    genomes = {accession: verify_genome(out_dir, accession) for accession in requested}
    status = "verified" if all(record["status"] == "verified" for record in genomes.values()) else "incomplete"
    receipt = {
        "schema": RECEIPT_SCHEMA, "status": status,
        "requested_accessions": requested, "genomes": genomes,
        "inputs": {"peers": file_record(peers_path), "selection": file_record(selection_path)},
    }
    save_json(Path(out_dir) / "acquisition.json", receipt)
    return receipt


def acquire(peers_path: Path, selection_path: Path, out_dir: Path, only=None,
            cache_dir=None, opener=A.default_opener) -> dict:
    peers_path, selection_path, out_dir = Path(peers_path), Path(selection_path), Path(out_dir)
    accessions, by_accession, selection = load_peer_inputs(peers_path, selection_path)
    requested = choose_only(accessions, only)
    cohort = cohort_accessions(selection, selection_path)
    configured_cache = cache_dir or (selection.get("cli_args") or {}).get("offline_cache")
    cache = SelectionCache(selection, Path(configured_cache) if configured_cache else out_dir / "ncbi-cache",
                           opener)
    reports, report_sources = cache.reports(requested)
    # Complete the full twin/cohort preflight before the first package download.
    for accession in requested:
        if reports[accession].get("accession") != accession:
            raise AcquisitionError(f"cached report accession mismatch for {accession}")
        enforce_exclusions(accession, paired_accession(reports[accession]), cohort, out_dir)
        pair = paired_accession(reports[accession])
        if pair in requested:
            raise AcquisitionError(f"{accession}: paired accession {pair} is also requested")
    taxonomy, taxonomy_sources = cache.taxonomy()
    out_dir.mkdir(parents=True, exist_ok=True)
    for accession in requested:
        cached_report = reports[accession]
        package = ensure_package(accession, out_dir, opener)
        package_report = read_json(out_dir / accession / "package" / "assembly_data_report.json")
        if organism_taxid(package_report) != organism_taxid(cached_report):
            raise AcquisitionError(f"{accession}: package and cached report organism taxids differ")
        cached_pair, package_pair = paired_accession(cached_report), paired_accession(package_report)
        if cached_pair and package_pair and cached_pair != package_pair:
            raise AcquisitionError(f"{accession}: package and cached report paired accessions differ")
        # A pair absent from the selection page still receives the same fence.
        enforce_exclusions(accession, package_pair or cached_pair, cohort, out_dir)
        build_outputs(accession, package_report, taxonomy, taxonomy_sources,
                      report_sources[accession], package, by_accession[accession], out_dir)
    return make_receipt(peers_path, selection_path, out_dir, requested)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    steps = parser.add_subparsers(dest="step", required=True)
    for name in ("acquire", "receipt"):
        step = steps.add_parser(name)
        step.add_argument("--peers", type=Path, required=True)
        step.add_argument("--selection", type=Path, required=True)
        step.add_argument("--out", type=Path, required=True)
        step.add_argument("--only", nargs="+", metavar="ACCESSION")
        if name == "acquire":
            step.add_argument("--cache", type=Path, help="raw NCBI Datasets/Taxonomy response cache")
    return parser


def main(argv=None, opener=A.default_opener) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.step == "acquire":
            result = acquire(args.peers, args.selection, args.out, args.only, args.cache, opener)
        else:
            result = make_receipt(args.peers, args.selection, args.out, args.only)
        print(json.dumps({"status": result["status"], "requested": len(result["requested_accessions"]),
                          "verified": sum(record["status"] == "verified"
                                          for record in result["genomes"].values())},
                         indent=2, sort_keys=True))
        return 0
    except (AcquisitionError, FileNotFoundError, KeyError, ET.ParseError,
            json.JSONDecodeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
