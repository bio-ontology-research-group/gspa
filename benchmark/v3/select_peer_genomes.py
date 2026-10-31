#!/usr/bin/env python3
"""Select pinned, unlabelled close-relative assemblies for the label-rich cohort.

For every selected target, NCBI Datasets v2 is queried for complete assemblies in the
target's family.  Close peers come from other species (at most a configurable number per
species); same-species assemblies are selected separately as near-duplicate controls.
Only selection manifests and reports are written: no genome package is downloaded.

Raw NCBI responses are cached by the SHA-256 of their URL.  A rerun with the same output
and cache directories therefore needs no network and ranks exactly the same reports.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import acquire_training_cohort as cohort  # noqa: E402

SCHEMA = "gspa.peer-genome-selection/1"
NCBI_FAMILY_REPORT = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{taxid}/dataset_report"
TSV_COLUMNS = ("target_taxid", "target_accession", "peer_accession", "peer_species_taxid",
               "peer_organism", "role")


class PeerSelectionError(ValueError):
    pass


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def family_report_url(family_taxid, source="refseq", assembly_level="complete", page_token=None) -> str:
    level = {"complete": "complete_genome"}.get(assembly_level, assembly_level)
    query = [("filters.assembly_level", level),
             ("filters.assembly_source", source),
             ("filters.exclude_atypical", "true")]
    if page_token:
        query.append(("page_token", page_token))
    return NCBI_FAMILY_REPORT.format(taxid=family_taxid) + "?" + urllib.parse.urlencode(query)


class ResponseCache:
    """Raw response cache with cohort retry/pause semantics and a 3-request/s limit."""

    def __init__(self, directory: Path, opener=cohort.default_opener, min_interval=0.34):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.opener = opener
        self.min_interval = min_interval
        self.last_request = None

    def path_for(self, url: str, suffix: str) -> Path:
        return self.directory / (sha256_bytes(url.encode()) + suffix)

    def _network_open(self, url: str) -> bytes:
        if self.last_request is not None:
            delay = self.min_interval - (time.monotonic() - self.last_request)
            if delay > 0:
                time.sleep(delay)
        try:
            return self.opener(url)
        finally:
            self.last_request = time.monotonic()

    def get(self, url: str, suffix=".json") -> dict:
        path = self.path_for(url, suffix)
        if path.exists():
            body = path.read_bytes()
        else:
            fetched = cohort.fetch(url, opener=self._network_open)
            body = fetched["body"]
            path.write_bytes(body)
        return {"body": body, "url": url, "url_sha256": sha256_bytes(url.encode()),
                "sha256": sha256_bytes(body), "size_bytes": len(body),
                "cache_file": str(path)}


def load_inputs(selection_path: Path, provenance_path: Path) -> tuple:
    selection = json.loads(Path(selection_path).read_text())
    if selection.get("schema") != "gspa.label-rich-selection/1":
        raise PeerSelectionError(f"unexpected selection schema {selection.get('schema')!r}")
    provenance = json.loads(Path(provenance_path).read_text())
    targets, skipped = [], []
    for organism in selection.get("organisms", []):
        taxid = str(organism["taxid"])
        assembly = (provenance.get("organisms", {}).get(taxid) or {}).get("assembly") or {}
        accession = assembly.get("used_accession")
        if not accession:
            # An organism whose assembly never resolved (e.g. taxid 431944, NCBI taxid
            # mismatch) is not a cohort genome: record it, never invent a target.
            skipped.append({"taxid": taxid, "status": assembly.get("status", "unresolved"),
                            "reason": assembly.get("reason", "no used_accession in provenance")})
            continue
        targets.append({"selection": organism, "taxid": taxid, "target_accession": accession})
    if not targets:
        raise PeerSelectionError("selection contains no resolved organisms")
    # Exclude both members of every cohort assembly pair: the GenBank/RefSeq twin of the
    # target is the same sequence, not a near-duplicate strain.
    cohort_accessions = {acc for record in provenance.get("organisms", {}).values()
                         for a in [record.get("assembly") or {}]
                         for acc in (a.get("used_accession"), a.get("paired_accession")) if acc}
    inputs = {
        "selection": {"path": str(selection_path), "sha256": cohort.sha256_file(selection_path)},
        "provenance": {"path": str(provenance_path), "sha256": cohort.sha256_file(provenance_path)},
        "skipped_unresolved_targets": skipped,
    }
    return targets, cohort_accessions, inputs


def species_from_taxonomy(root: ET.Element) -> dict:
    """Map each taxon to its species rank, or to itself when no species rank exists."""
    found = {}
    for taxon in root.findall("Taxon"):
        own_id = taxon.findtext("TaxId")
        species = None
        if taxon.findtext("Rank") == "species":
            species = (own_id, taxon.findtext("ScientificName"))
        if species is None:
            for parent in taxon.findall("LineageEx/Taxon"):
                if parent.findtext("Rank") == "species":
                    species = (parent.findtext("TaxId"), parent.findtext("ScientificName"))
        if own_id:
            found[str(own_id)] = {"taxid": str(species[0]) if species else str(own_id),
                                  "name": species[1] if species else taxon.findtext("ScientificName"),
                                  "source": "ncbi_taxonomy" if species else "no_species_rank"}
    return found


def target_species(path: Path) -> dict:
    """Read species taxids from the taxonomy XML already fetched by fetch_taxonomy."""
    return species_from_taxonomy(ET.parse(path).getroot())


def resolve_species(taxids, cache: ResponseCache) -> tuple:
    """Resolve organism taxids through cached NCBI Taxonomy requests of at most 200 ids."""
    ordered = sorted({str(taxid) for taxid in taxids if taxid is not None})
    found, batches = {}, []
    for start in range(0, len(ordered), 200):
        batch = ordered[start:start + 200]
        url = cohort.NCBI_TAXONOMY_XML.format(ids=",".join(batch))
        fetched = cache.get(url, ".xml")
        batches.append({k: fetched[k] for k in
                        ("url", "url_sha256", "sha256", "size_bytes", "cache_file")})
        found.update(species_from_taxonomy(ET.fromstring(fetched["body"])))
    missing = sorted(set(ordered) - set(found))
    if missing:
        raise PeerSelectionError(f"NCBI Taxonomy returned no record for candidate taxids: {','.join(missing)}")
    return found, batches


def organism_genus(organism: dict) -> str:
    if organism.get("genus"):
        return str(organism["genus"])
    name = organism.get("organism_name") or ""
    tokens = name.replace("Candidatus ", "").split()
    return tokens[0].strip("[]") if tokens else ""


def candidate_from_report(report: dict) -> dict:
    organism = report.get("organism") or {}
    stats = report.get("assembly_stats") or {}
    organism_taxid = organism.get("tax_id")
    species_taxid = (organism.get("species_tax_id") or organism.get("species_taxid")
                     or organism_taxid)
    return {
        "accession": report.get("accession"),
        "organism": organism.get("organism_name"),
        "organism_taxid": str(organism_taxid) if organism_taxid is not None else None,
        "species_taxid": str(species_taxid) if species_taxid is not None else None,
        "species_name": organism.get("species_name"),
        "species_source": "reported_species" if organism.get("species_tax_id") or organism.get("species_taxid") else "unresolved",
        "genus": organism_genus(organism),
        "assembly_level": (report.get("assembly_info") or {}).get("assembly_level"),
        "refseq_category": report.get("refseq_category"),
        "release_date": report.get("release_date") or "",
        "contig_n50": stats.get("contig_n50"),
        "n_contigs": stats.get("number_of_contigs"),
    }


def preferred_category(candidate: dict) -> bool:
    category = (candidate.get("refseq_category") or "").lower()
    return "reference" in category or "representative" in category


def rank_candidates(candidates: list, target_genus: str, refseq_category_available=True) -> list:
    """Apply each available ranking criterion in order, using stable sorts and accession ties."""
    genus = target_genus.lower()
    ordered = sorted(candidates, key=lambda c: c["accession"] or "")
    ordered.sort(key=lambda c: c.get("contig_n50") or 0, reverse=True)
    ordered.sort(key=lambda c: c.get("release_date") or "", reverse=True)
    if refseq_category_available:
        ordered.sort(key=preferred_category, reverse=True)
    ordered.sort(key=lambda c: (c.get("genus") or "").lower() == genus, reverse=True)
    return ordered


def rank_reasons(candidate: dict, target_genus: str, refseq_category_available=True) -> list:
    same_genus = (candidate.get("genus") or "").lower() == target_genus.lower()
    reasons = [same_genus and "same_genus" or "same_family_only"]
    if refseq_category_available:
        reasons.append("reference_or_representative" if preferred_category(candidate)
                       else "not_reference_or_representative")
    reasons.extend([f"release_date:{candidate.get('release_date') or 'missing'}",
                    f"contig_n50:{candidate.get('contig_n50') if candidate.get('contig_n50') is not None else 'missing'}",
                    "accession_tiebreak"])
    return reasons


def select_for_target(reports: list, target: dict, cohort_accessions: set, *, peers_per_target=6,
                      same_species=1, max_per_species=1, species_by_taxid=None) -> dict:
    """Normalize, exclude, rank and select one target's close and near-duplicate peers."""
    refseq_category_available = any("refseq_category" in report for report in reports)
    by_accession = {}
    for report in reports:
        candidate = candidate_from_report(report)
        if species_by_taxid is not None:
            resolved = species_by_taxid.get(candidate.get("organism_taxid"))
            if resolved is not None:
                candidate.update(species_taxid=resolved["taxid"], species_name=resolved["name"],
                                 species_source=resolved["source"])
        accession = candidate.get("accession")
        if accession and accession not in by_accession:
            by_accession[accession] = candidate
    candidates = [c for c in by_accession.values()
                  if c["accession"] not in cohort_accessions | {target["target_accession"]}]
    species = str(target["species_taxid"])
    near_pool = [c for c in candidates if c.get("species_taxid") == species]
    close_pool = [c for c in candidates if c.get("species_taxid") != species]

    near = rank_candidates(near_pool, target["genus"], refseq_category_available)[:same_species]
    close, per_species = [], {}
    for candidate in rank_candidates(close_pool, target["genus"], refseq_category_available):
        if len(close) >= peers_per_target:
            break
        key = candidate.get("species_taxid") or "accession:" + candidate["accession"]
        if per_species.get(key, 0) >= max_per_species:
            continue
        close.append(candidate)
        per_species[key] = per_species.get(key, 0) + 1

    peers = []
    for role, chosen in (("close_peer", close), ("near_duplicate", near)):
        for candidate in chosen:
            peers.append({**candidate, "role": role,
                          "rank_reasons": rank_reasons(candidate, target["genus"], refseq_category_available)})
    return {"candidates_found": len(candidates), "peers": peers,
            "refseq_category_available": refseq_category_available,
            "shortfalls": {"close_peers": max(0, peers_per_target - len(close)),
                           "near_duplicates": max(0, same_species - len(near))}}


def fetch_family_reports(family_taxid: str, cache: ResponseCache, source: str, assembly_level: str) -> tuple:
    reports, pages = [], []
    token = None
    seen_tokens = set()
    while True:
        url = family_report_url(family_taxid, source, assembly_level, token)
        fetched = cache.get(url, ".json")
        page = json.loads(fetched["body"])
        reports.extend(page.get("reports") or [])
        pages.append({k: fetched[k] for k in ("url", "url_sha256", "sha256", "size_bytes", "cache_file")})
        token = page.get("next_page_token")
        if not token:
            break
        if token in seen_tokens:
            raise PeerSelectionError(f"family {family_taxid}: repeated next_page_token {token!r}")
        seen_tokens.add(token)
    return reports, pages


def cli_args_dict(args: argparse.Namespace) -> dict:
    return {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()}


def write_outputs(out_dir: Path, manifest: dict) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "peer_selection.json").write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
    with open(out_dir / "peers.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter="\t", lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader()
        for target in manifest["targets"]:
            for peer in target["peers"]:
                writer.writerow({"target_taxid": target["taxid"],
                                 "target_accession": target["target_accession"],
                                 "peer_accession": peer["accession"],
                                 "peer_species_taxid": peer["species_taxid"],
                                 "peer_organism": peer["organism"], "role": peer["role"]})
    lines = ["# Close-relative peer genome selection", "",
             "| Target taxid | Family | Candidates found | Close peers | Near-duplicates | Shortfall |",
             "|---|---|---:|---:|---:|---|"]
    for target in manifest["targets"]:
        short = target["shortfalls"]
        lines.append("| {taxid} | {family} | {candidates_found} | {close} | {near} | close={close_short}; near_duplicate={near_short} |".format(
            taxid=target["taxid"], family=target["family"], candidates_found=target["candidates_found"],
            close=sum(p["role"] == "close_peer" for p in target["peers"]),
            near=sum(p["role"] == "near_duplicate" for p in target["peers"]),
            close_short=short["close_peers"], near_short=short["near_duplicates"]))
    drifted = [target for target in manifest["targets"]
               if not target["taxonomy_verification"]["matches"]]
    if drifted:
        lines.extend(["", "## Taxonomy notes", ""])
        for target in drifted:
            verification = target["taxonomy_verification"]
            selection, ncbi = verification["selection"], verification["ncbi"]
            lines.append(
                f"- Target {target['taxid']}: selection {selection.get('genus')}/{selection.get('family')}; "
                f"NCBI {ncbi.get('genus')}/{ncbi.get('family')}. "
                f"Policy: `{verification['policy']}`.")
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n")


def run(args: argparse.Namespace, opener=cohort.default_opener) -> dict:
    targets, cohort_accessions, inputs = load_inputs(args.selection, args.provenance)
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)  # fetch_taxonomy writes here before write_outputs
    cache = ResponseCache(args.offline_cache, opener)
    taxids = [target["taxid"] for target in targets]
    taxonomy_url = cohort.NCBI_TAXONOMY_XML.format(ids=",".join(sorted(taxids)))
    taxonomy_cached = cache.get(taxonomy_url, ".xml")
    # Reuse the acquisition helper unchanged; its read is served from the raw-response cache.
    classifications = cohort.fetch_taxonomy(taxids, args.out_dir, lambda url: cache.get(url, ".xml")["body"])
    species = target_species(args.out_dir / "selected_taxonomy.xml")

    manifest_targets, all_pages, all_taxonomy_batches = [], {}, {}
    for target in targets:
        selected = target["selection"]
        taxid = target["taxid"]
        ncbi = classifications[taxid]
        genus, family = selected.get("genus"), selected.get("family")
        matches = genus == ncbi.get("genus") and family == ncbi.get("family")
        # The selection recorded the taxonomy at acquisition time; NCBI may have moved a
        # family since (Aeropyrum: Desulfurococcaceae -> Acidilobaceae). The family query
        # must use NCBI's current placement, so the disagreement is recorded, not fatal.
        if not matches:
            print(json.dumps({"taxonomy_mismatch": taxid, "selection": [genus, family],
                              "ncbi": [ncbi.get("genus"), ncbi.get("family")]}), file=sys.stderr)
        if not ncbi.get("family_taxid"):
            raise PeerSelectionError(f"target {taxid}: NCBI taxonomy has no family taxid")
        target_record = {**target, "genus": ncbi.get("genus") or genus, "family": ncbi.get("family") or family,
                         "species_taxid": species[taxid]["taxid"]}
        reports, pages = fetch_family_reports(str(ncbi["family_taxid"]), cache, args.source, args.assembly_level)
        organism_taxids = [(report.get("organism") or {}).get("tax_id") for report in reports]
        candidate_species, taxonomy_batches = resolve_species(organism_taxids, cache)
        result = select_for_target(reports, target_record, cohort_accessions,
                                   peers_per_target=args.peers_per_target, same_species=args.same_species,
                                   max_per_species=args.max_per_species, species_by_taxid=candidate_species)
        verification = {"selection": {"genus": genus, "family": family},
                        "ncbi": {"genus": ncbi.get("genus"), "genus_taxid": ncbi.get("genus_taxid"),
                                 "family": ncbi.get("family"), "family_taxid": ncbi.get("family_taxid")},
                        "matches": matches}
        if not matches:
            verification["policy"] = "ncbi_current_placement_used"
        manifest_targets.append({
            "taxid": taxid, "genus": target_record["genus"], "family": target_record["family"],
            "species_taxid": target_record["species_taxid"],
            "target_accession": target["target_accession"], "proposed_split": selected.get("proposed_split"),
            "taxonomy_verification": verification,
            "ncbi_page_urls": [p["url"] for p in pages], **result})
        for page in pages:
            all_pages[page["url"]] = page
        for batch in taxonomy_batches:
            all_taxonomy_batches[batch["url"]] = batch

    manifest = {"schema": SCHEMA, "created_utc": cohort.utcnow(), "cli_args": cli_args_dict(args),
                "inputs": inputs,
                "ncbi_taxonomy": {k: taxonomy_cached[k] for k in
                                  ("url", "url_sha256", "sha256", "size_bytes", "cache_file")},
                "ncbi_pages": [all_pages[url] for url in sorted(all_pages)],
                "ncbi_taxonomy_batches": [all_taxonomy_batches[url] for url in sorted(all_taxonomy_batches)],
                "targets": manifest_targets}
    write_outputs(args.out_dir, manifest)
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--peers-per-target", type=int, default=6)
    parser.add_argument("--same-species", type=int, default=1)
    parser.add_argument("--max-per-species", type=int, default=1)
    parser.add_argument("--assembly-level", default="complete")
    parser.add_argument("--source", choices=("refseq", "genbank"), default="refseq")
    parser.add_argument("--offline-cache", type=Path)
    return parser


def main(argv=None, opener=cohort.default_opener) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.offline_cache is None:
        args.offline_cache = args.out_dir / "ncbi-cache"
    if min(args.peers_per_target, args.same_species, args.max_per_species) < 0:
        parser.error("peer counts and species cap must be non-negative")
    manifest = run(args, opener)
    print(json.dumps({"targets": len(manifest["targets"]),
                      "peers": sum(len(t["peers"]) for t in manifest["targets"]),
                      "shortfalls": sum(sum(t["shortfalls"].values()) for t in manifest["targets"])},
                     indent=1, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
