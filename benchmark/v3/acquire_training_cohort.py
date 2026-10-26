#!/usr/bin/env python3
"""Acquire the label-rich training cohort from public sources and audit it against Price.

The selection (``results/label-rich-selection/selection_manifest.json``) is read and never edited. This
tool writes a NEW provenance manifest. Subcommands, each restartable (files already on disk are reused
only when their recorded hash still matches):

``proteomes``   fetch the public UniProt proteome JSON of every selected organism's dominant proteome and
                check it against the selection: proteome id, NCBI taxid, organism name/strain, proteome type,
                and a well-formed assembly accession (``GCA_`` / ``GCF_`` with an explicit version).
``assemblies``  resolve each declared assembly through NCBI Datasets (``dataset_report``), check accession,
                version currency and taxid, then download the exact package (protein FASTA, GFF3, genomic DNA)
                and extract it in the layout ``prepare_price_layout.prepare`` expects. When the declared
                GenBank assembly carries no protein FASTA, its RefSeq pair is used ONLY if NCBI reports that
                pair as ``current``, and the substitution is recorded; no other genome is ever substituted.
``exclusions``  audit every selected organism against ALL Price organisms: the NCBI taxonomy dump of the Price
                cohort, the committed Price acquisition manifests, the Price host list and the evaluation
                exclusions registry. A shared genus (by name, synonym or NCBI genus taxid) removes the organism;
                a shared family is reported.

Every network read records URL, UTC timestamp, size and SHA-256. Nothing here guesses a mapping: an
unresolvable organism stays in the manifest with ``status`` and ``reason``.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import sys
import time
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prepare_price_layout as PL  # noqa: E402  (FASTA/GFF reading and layout writing are reused unchanged)
from select_label_rich_genomes import load_exclusions  # noqa: E402

REPO = HERE.parents[1]
SELECTION_DEFAULT = HERE / "results" / "label-rich-selection" / "selection_manifest.json"
DATA_DEFAULT = Path.home() / ".cache" / "gspa-v3-data" / "training-cohort"
PRICE_FAMILIES_DEFAULT = Path.home() / ".cache" / "gspa-v3-data" / "price-families"
SCHEMA = "gspa.training-cohort-provenance/1"
UNIPROT_PROTEOME = "https://rest.uniprot.org/proteomes/{id}"
NCBI_REPORT = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/dataset_report"
NCBI_PACKAGE = ("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download"
                "?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=PROT_FASTA")
NCBI_TAXONOMY_XML = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&retmode=xml&id={ids}"
ASSEMBLY_RE = re.compile(r"^GC[AF]_\d{9}\.\d+$")
USER_AGENT = "gspa-v3-training-cohort/1 (public read-only acquisition)"


class AcquisitionError(ValueError):
    pass


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path) -> str:
    return PL.digest(path)


def default_opener(url: str, timeout: int = 300) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "*/*"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def fetch(url: str, opener=default_opener, retries: int = 3, pause: float = 2.0) -> dict:
    """``{body, url, retrieved_utc, sha256, size_bytes}``; transient failures are retried, 404 is not."""
    last = None
    for attempt in range(retries):
        try:
            body = opener(url)
            return {"body": body, "url": url, "retrieved_utc": utcnow(), "sha256": sha256_bytes(body), "size_bytes": len(body)}
        except urllib.error.HTTPError as exc:
            last = exc
            if exc.code in (400, 404):
                break
        except (urllib.error.URLError, OSError) as exc:
            last = exc
        time.sleep(pause * (attempt + 1))
    raise AcquisitionError(f"cannot fetch {url}: {last}")


def provenance(f: dict) -> dict:
    return {k: f[k] for k in ("url", "retrieved_utc", "sha256", "size_bytes")}


def load_selection(path=SELECTION_DEFAULT) -> dict:
    sel = json.loads(Path(path).read_text())
    if sel.get("schema") != "gspa.label-rich-selection/1":
        raise AcquisitionError(f"unexpected selection schema {sel.get('schema')!r}")
    sel["_sha256"] = sha256_file(path)
    return sel


def read_manifest(path: Path) -> dict:
    if path.exists():
        return json.loads(path.read_text())
    return {"schema": SCHEMA, "organisms": {}}


def write_manifest(path: Path, manifest: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")


# ----------------------------------------------------------------------------- proteomes

def check_proteome(doc: dict, org: dict) -> dict:
    """Compare one UniProt proteome JSON with the selection record. Returns a dict of verdicts."""
    wanted = org["dominant_uniprot_proteome"]["id"]
    problems = []
    if doc.get("id") != wanted:
        problems.append(f"proteome_id_mismatch:{doc.get('id')}")
    taxon = doc.get("taxonomy") or {}
    if str(taxon.get("taxonId")) != str(org["taxid"]):
        problems.append(f"taxid_mismatch:selection={org['taxid']},proteome={taxon.get('taxonId')}")
    name = taxon.get("scientificName", "")
    if not name or not org["organism"].startswith(name):
        problems.append(f"organism_name_mismatch:{name!r}")
    if doc.get("proteomeType") not in ("Reference proteome", "Reference and representative proteome",
                                       "Representative proteome", "Other proteome"):
        problems.append(f"proteome_type:{doc.get('proteomeType')}")
    if doc.get("redundantTo"):
        problems.append(f"redundant_to:{doc['redundantTo']}")
    asm = (doc.get("genomeAssembly") or {}).get("assemblyId", "")
    if not ASSEMBLY_RE.match(asm or ""):
        problems.append(f"assembly_accession_missing_or_unversioned:{asm!r}")
    return {"ok": not problems, "problems": problems, "assembly_id": asm or None,
            "taxon_id": taxon.get("taxonId"), "strain": doc.get("strain"), "scientific_name": name,
            "proteome_type": doc.get("proteomeType"), "protein_count": doc.get("proteinCount"),
            "assembly_level": (doc.get("genomeAssembly") or {}).get("level"),
            "assembly_source": (doc.get("genomeAssembly") or {}).get("source"),
            "genome_accessions": sorted(x["id"] for c in doc.get("components", [])
                                        for x in c.get("proteomeCrossReferences", []) if x.get("database") == "GenomeAccession")}


def acquire_proteomes(selection: dict, out_dir: Path, opener=default_opener) -> dict:
    manifest = read_manifest(out_dir / "provenance.json")
    (out_dir / "proteomes").mkdir(parents=True, exist_ok=True)
    for org in selection["organisms"]:
        pid = org["dominant_uniprot_proteome"]["id"]
        path = out_dir / "proteomes" / f"{pid}.json"
        rec = manifest["organisms"].setdefault(org["taxid"], {})
        f = fetch(UNIPROT_PROTEOME.format(id=pid), opener)
        path.write_bytes(f["body"])
        verdict = check_proteome(json.loads(f["body"]), org)
        rec.update({"taxid": org["taxid"], "organism": org["organism"], "proteome_id": pid,
                    "declared_strain": verdict["strain"], "proteome_file": str(path),
                    "proteome_source": provenance(f), "proteome_check": verdict})
    manifest["selection"] = {"file": str(SELECTION_DEFAULT.relative_to(REPO)), "sha256": selection["_sha256"]}
    write_manifest(out_dir / "provenance.json", manifest)
    return manifest


# ----------------------------------------------------------------------------- assemblies

def report_check(report: dict, accession: str, taxid: str) -> list:
    problems = []
    if report.get("accession") != accession:
        problems.append(f"report_accession:{report.get('accession')}")
    if report.get("current_accession") != accession:
        problems.append(f"version_superseded:{report.get('current_accession')}")
    status = (report.get("assembly_info") or {}).get("assembly_status")
    if status != "current":
        problems.append(f"assembly_status:{status}")
    got = str((report.get("organism") or {}).get("tax_id"))
    if got != str(taxid):
        problems.append(f"ncbi_taxid_mismatch:selection={taxid},ncbi={got}")
    return problems


def fetch_report(accession: str, opener) -> tuple:
    f = fetch(NCBI_REPORT.format(acc=accession), opener)
    reports = json.loads(f["body"]).get("reports") or []
    if len(reports) != 1:
        raise AcquisitionError(f"{accession}: dataset_report returned {len(reports)} reports")
    return reports[0], f


def extract_package(zip_bytes: bytes, accession: str, dest: Path) -> dict:
    """Write protein.faa / genomic.gff / genomic.fna / assembly_data_report.json for exactly ``accession``."""
    zf = zipfile.ZipFile(io.BytesIO(zip_bytes))
    names = zf.namelist()
    base = f"ncbi_dataset/data/{accession}/"
    wanted = {"protein.faa": [n for n in names if n == base + "protein.faa"],
              "genomic.gff": [n for n in names if n == base + "genomic.gff"],
              "genomic.fna": [n for n in names if n.startswith(base) and n.endswith("_genomic.fna")]}
    reports = [n for n in names if n.endswith("assembly_data_report.jsonl")]
    if len(reports) != 1:
        raise AcquisitionError(f"{accession}: expected one assembly_data_report.jsonl, found {len(reports)}")
    lines = [l for l in zf.read(reports[0]).decode().splitlines() if l.strip()]
    if len(lines) != 1 or json.loads(lines[0]).get("accession") != accession:
        raise AcquisitionError(f"{accession}: package report names a different assembly")
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "assembly_data_report.json").write_text(lines[0] + "\n")
    out = {"assembly_data_report.json": sha256_file(dest / "assembly_data_report.json")}
    for local, members in wanted.items():
        if len(members) != 1:
            out[local] = None
            continue
        (dest / local).write_bytes(zf.read(members[0]))
        out[local] = sha256_file(dest / local)
    return out


def acquire_assemblies(out_dir: Path, opener=default_opener) -> dict:
    manifest = read_manifest(out_dir / "provenance.json")
    for taxid, rec in sorted(manifest["organisms"].items()):
        check = rec.get("proteome_check") or {}
        if not check.get("ok"):
            rec["assembly"] = {"status": "unresolved", "reason": "proteome_check_failed:" + ";".join(check.get("problems", ["not_run"]))}
            continue
        declared = check["assembly_id"]
        asm = {"declared_accession": declared, "declared_by": "UniProt proteome " + rec["proteome_id"]}
        report, f = fetch_report(declared, opener)
        asm["declared_report_source"] = provenance(f)
        problems = report_check(report, declared, taxid)
        asm["declared_report_problems"] = problems
        asm["ncbi_organism"] = report.get("organism")
        asm["assembly_level"] = (report.get("assembly_info") or {}).get("assembly_level")
        asm["paired_accession"] = report.get("paired_accession")
        if problems:
            rec["assembly"] = {**asm, "status": "unresolved", "reason": ";".join(problems)}
            continue
        used, used_report = declared, report
        dest_root = out_dir / "genomes"
        pkg = fetch(NCBI_PACKAGE.format(acc=used), opener)
        info = extract_package(pkg["body"], used, dest_root / used)
        asm["used_package_source"] = provenance(pkg)
        if info["protein.faa"] is None or info["genomic.gff"] is None:
            pair = report.get("paired_accession")
            pair_report, pf = fetch_report(pair, opener) if pair else (None, None)
            pair_status = ((report.get("assembly_info") or {}).get("paired_assembly") or {}).get("status")   # as NCBI states it for the declared assembly
            if not pair or report_check(pair_report, pair, taxid) or pair_status != "current":
                rec["assembly"] = {**asm, "status": "unresolved", "reason": "declared_assembly_has_no_protein_annotation_and_no_current_pair"}
                continue
            pkg = fetch(NCBI_PACKAGE.format(acc=pair), opener)
            info = extract_package(pkg["body"], pair, dest_root / pair)
            asm.update({"substituted_with_pair": True, "pair_report_source": provenance(pf), "used_package_source": provenance(pkg)})
            used, used_report = pair, pair_report
            if info["protein.faa"] is None or info["genomic.gff"] is None:
                rec["assembly"] = {**asm, "status": "unresolved", "reason": "pair_has_no_protein_annotation"}
                continue
        seqs = PL.fasta_sequences(dest_root / used / "protein.faa")
        asm.update({"status": "resolved", "used_accession": used, "directory": str(dest_root / used), "member_sha256": info,
                    "proteins": len(seqs), "annotation": (used_report.get("annotation_info") or {}).get("name"),
                    "annotation_release": (used_report.get("annotation_info") or {}).get("release_date"),
                    "ncbi_organism_used": used_report.get("organism")})
        rec["assembly"] = asm
    write_manifest(out_dir / "provenance.json", manifest)
    return manifest


# ----------------------------------------------------------------------------- Price exclusion audit

def read_price_taxonomy(path: Path) -> list:
    """Every top-level taxon of the Price taxonomy dump with its genus/family (name and taxid) and lineage ids."""
    out = []
    for t in ET.parse(path).getroot().findall("Taxon"):
        ranks = {}
        for x in t.findall("LineageEx/Taxon"):
            ranks.setdefault(x.findtext("Rank"), (x.findtext("ScientificName"), int(x.findtext("TaxId"))))
        own = (t.findtext("Rank"), t.findtext("ScientificName"), int(t.findtext("TaxId")))
        if own[0] in ("genus", "family"):
            ranks.setdefault(own[0], (own[1], own[2]))
        out.append({"taxid": own[2], "name": own[1], "rank": own[0],
                    "genus": ranks.get("genus"), "family": ranks.get("family"),
                    "lineage_ids": [int(x.findtext("TaxId")) for x in t.findall("LineageEx/Taxon")]})
    return out


def read_price_acquisitions(results_dir: Path) -> dict:
    """``{accession: {taxid, name}}`` from every committed Price acquisition manifest (jsonl with organism)."""
    found = {}
    for path in sorted(results_dir.glob("development-*/price-*acquisition.jsonl")):
        for line in path.read_text().splitlines():
            if not line.strip():
                continue
            rec = json.loads(line)
            org = rec.get("organism")
            if org and rec.get("accession_requested"):
                found[rec["accession_requested"]] = {"taxid": org.get("taxId"), "name": org.get("organismName"),
                                                     "manifest": str(path.relative_to(REPO))}
    return found


def read_price_hosts(manifest_path: Path) -> list:
    with open(manifest_path, newline="") as fh:
        return sorted({r["host_assembly"] for r in csv.DictReader(fh, delimiter="\t")})


def fetch_taxonomy(taxids, out_dir: Path, opener) -> dict:
    """NCBI Taxonomy (efetch XML, the format of the Price dump) for the selected taxids: ``{taxid: classification}``.

    Genus and family come from NCBI lineage ranks, independently of the Swiss-Prot lineage text; an organism
    whose genus cannot be read from NCBI stops the audit instead of passing it.
    """
    f = fetch(NCBI_TAXONOMY_XML.format(ids=",".join(sorted(taxids))), opener)
    path = out_dir / "selected_taxonomy.xml"
    path.write_bytes(f["body"])
    found = {}
    for t in read_price_taxonomy(path):          # same parser as the Price dump
        found[str(t["taxid"])] = {"genus": t["genus"][0] if t["genus"] else None, "genus_taxid": t["genus"][1] if t["genus"] else None,
                                  "family": t["family"][0] if t["family"] else None, "family_taxid": t["family"][1] if t["family"] else None,
                                  "ncbi_name": t["name"], "lineage_ids": t["lineage_ids"], "source": provenance(f)}
    missing = sorted(set(map(str, taxids)) - set(found))
    if missing:
        raise AcquisitionError(f"NCBI Taxonomy returned no record for {missing}")
    nogenus = sorted(t for t, v in found.items() if not v["genus_taxid"])
    if nogenus:
        raise AcquisitionError(f"NCBI Taxonomy gives no genus for {nogenus}: the Price audit cannot be run")
    return found


def audit_price_conflicts(org: dict, ncbi: dict, price_taxa: list, synonyms: dict, registry_genera: dict) -> dict:
    """Genus-level conflicts exclude; family-level overlaps are reported. Matching is by name, registry synonym and NCBI taxid."""
    names = {n for n in (org.get("genus"), ncbi.get("genus")) if n}
    for n in list(names):
        names.update(synonyms.get(n, []))
    names = {n.lower() for n in names}
    genus_ids = {i for i in (ncbi.get("genus_taxid"),) if i}
    family_ids = {i for i in (ncbi.get("family_taxid"),) if i}
    family_names = {n.lower() for n in (org.get("family"), ncbi.get("family")) if n}
    genus_hits, family_hits = [], []
    for t in price_taxa:
        pg, pf = t.get("genus"), t.get("family")
        pg_names = {pg[0].lower()} | {s.lower() for s in synonyms.get(pg[0], [])} if pg else set()
        if pg and (pg_names & names or pg[1] in genus_ids or ncbi.get("genus_taxid") in t["lineage_ids"] or t["taxid"] in (ncbi.get("genus_taxid"),)):
            genus_hits.append({"price_taxid": t["taxid"], "price_name": t["name"], "price_genus": pg[0], "basis": "genus"})
        elif pf and (pf[0].lower() in family_names or pf[1] in family_ids):
            family_hits.append({"price_taxid": t["taxid"], "price_name": t["name"], "price_family": pf[0]})
    registry_hit = next((registry_genera[n] for n in names if n in registry_genera), None)
    return {"genus_conflicts": genus_hits, "family_overlaps": family_hits, "registry_genus_conflict": registry_hit,
            "excluded": bool(genus_hits or registry_hit),
            "checked_names": sorted(names), "checked_genus_taxid": ncbi.get("genus_taxid")}


def audit_exclusions(selection: dict, out_dir: Path, price_families: Path, opener=default_opener, results_dir=HERE / "results") -> dict:
    manifest = read_manifest(out_dir / "provenance.json")
    tax_path = price_families / "taxonomy.xml"
    price_taxa = read_price_taxonomy(tax_path)
    acquisitions = read_price_acquisitions(results_dir)
    hosts = read_price_hosts(REPO / "benchmark" / "manifests" / "suite6_price149_manifest.tsv")
    excl = load_exclusions()
    synonyms = json.loads((HERE / "registry" / "evaluation_exclusions.json").read_text()).get("genus_synonyms", {})
    # Price taxa may be known only from the acquisition manifests (not the taxonomy dump): their names count too.
    dump_ids = {t["taxid"] for t in price_taxa}
    acquisition_only = {a: v for a, v in acquisitions.items() if int(v["taxid"]) not in dump_ids}
    unresolved_hosts = sorted(set(hosts) - set(acquisitions))
    classification = fetch_taxonomy(manifest["organisms"], out_dir, opener)
    for taxid, rec in sorted(manifest["organisms"].items()):
        org = next(o for o in selection["organisms"] if o["taxid"] == taxid)
        ncbi = classification[taxid]
        verdict = audit_price_conflicts(org, ncbi, price_taxa, synonyms, excl["genera"])
        rec["price_audit"] = {**verdict, "ncbi_classification": ncbi,
                              "selection_price_genus_check_was": org["price_genus_check"],
                              "family_shared_with_excluded_panel_organism": org["family_shared_with_excluded_panel_organism"]}
    manifest["price_audit_inputs"] = {
        "taxonomy_xml": {"path": str(tax_path), "sha256": sha256_file(tax_path), "taxa": len(price_taxa),
                         "genera": sorted({t["genus"][0] for t in price_taxa if t["genus"]}),
                         "families": sorted({t["family"][0] for t in price_taxa if t["family"]})},
        "acquisition_manifests": {"assemblies": len(acquisitions), "taxids": sorted({int(v["taxid"]) for v in acquisitions.values()}),
                                  "not_in_taxonomy_dump": acquisition_only},
        "price_host_assemblies": len(hosts), "price_hosts_without_organism_record": unresolved_hosts,
        "evaluation_exclusions_registry_sha256": excl["sha256"]}
    manifest["excluded_by_price_audit"] = sorted(t for t, r in manifest["organisms"].items() if r["price_audit"]["excluded"])
    write_manifest(out_dir / "provenance.json", manifest)
    return manifest


# ----------------------------------------------------------------------------- layouts

def build_layouts(out_dir: Path, layout_dir: Path) -> dict:
    """Reuse ``prepare_price_layout.prepare`` on the resolved assembly directories (coordinates only, no labels)."""
    manifest = read_manifest(out_dir / "provenance.json")
    keep = {r["assembly"]["used_accession"] for r in manifest["organisms"].values() if r.get("assembly", {}).get("status") == "resolved"}
    import shutil
    staging = out_dir / "layout_input"
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)
    for acc in sorted(keep):
        (staging / acc).mkdir()
        for name in ("protein.faa", "genomic.gff", "assembly_data_report.json"):
            (staging / acc / name).symlink_to((out_dir / "genomes" / acc / name).resolve())
    if layout_dir.exists():
        shutil.rmtree(layout_dir)
    result = PL.prepare(staging, layout_dir)
    manifest["layout"] = {"dir": str(layout_dir), "report_sha256": sha256_file(layout_dir / "report.json"),
                          "all_proteins_sha256": result["all_proteins_sha256"]}
    write_manifest(out_dir / "provenance.json", manifest)
    return result


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=("proteomes", "assemblies", "exclusions", "layouts"))
    ap.add_argument("--selection", type=Path, default=SELECTION_DEFAULT)
    ap.add_argument("--data-dir", type=Path, default=DATA_DEFAULT)
    ap.add_argument("--price-families", type=Path, default=PRICE_FAMILIES_DEFAULT)
    ap.add_argument("--layout-dir", type=Path)
    ap.add_argument("--manifest-copy", type=Path, help="also write the provenance manifest here (committed copy)")
    args = ap.parse_args(argv)
    selection = load_selection(args.selection)
    if args.step == "proteomes":
        manifest = acquire_proteomes(selection, args.data_dir)
    elif args.step == "assemblies":
        manifest = acquire_assemblies(args.data_dir)
    elif args.step == "exclusions":
        manifest = audit_exclusions(selection, args.data_dir, args.price_families)
    else:
        build_layouts(args.data_dir, args.layout_dir or args.data_dir / "layout")
        manifest = read_manifest(args.data_dir / "provenance.json")
    if args.manifest_copy:
        write_manifest(args.manifest_copy, manifest)
    orgs = manifest["organisms"].values()
    print(json.dumps({"organisms": len(manifest["organisms"]),
                      "proteome_ok": sum(1 for r in orgs if r.get("proteome_check", {}).get("ok")),
                      "assembly_resolved": sum(1 for r in orgs if r.get("assembly", {}).get("status") == "resolved"),
                      "price_excluded": manifest.get("excluded_by_price_audit", [])}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
