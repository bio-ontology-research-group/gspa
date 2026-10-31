#!/usr/bin/env python3
"""Build per-batch manifests for Empty Quarter (17,055) and benchmark (145) genomes.

Produces:
  - <work_root>/manifest_all.tsv: full genome table with resolved input locations
  - <work_root>/batches_index.tsv: summary of each batch (batch_id, n_genomes, datasets, n_tars)
  - <work_root>/batches/batch_NNN/batch_manifest.tsv: per-batch genome manifest (~500 genomes/batch)

Grouping rules:
  - Soil metagenome bins (dataset == "soil_bin") are grouped by sample so that each
    site tar (<trip>/<site>/<sample>.tar) is read by exactly one batch.
  - Culture genomes (RAK-PR catalogue) are packed into subsequent batches of up to batch_size.
  - Benchmark genomes (22 Price matrix + 18 recurrence cohort + 105 recurrence peers) use
    the exact same batch_manifest.tsv schema so the same staging and gapsmith jobs run on them.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


MANIFEST_COLUMNS = [
    "batch_id",
    "mag_id",
    "dataset",
    "sample",
    "bin",
    "domain",
    "protein_source",
    "expected_n_proteins",
    "gapseq_source_ibex",
    "source_kind",
    "tar_path",
    "bin_fasta_member",
    "bakta_prefix",
    "direct_faa_path",
    "direct_gff3_path",
    "direct_fna_path",
    "expected_faa_sha256",
]

BATCH_INDEX_COLUMNS = [
    "batch_id",
    "batch_index",
    "n_genomes",
    "n_samples",
    "n_tars",
    "n_bakta",
    "n_pyrodigal",
    "n_ncbi_proteome",
    "n_archaea",
    "n_bacteria",
    "datasets",
    "manifest_path",
]


def count_fasta_records(path: Path) -> int:
    count = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                count += 1
    return count


def load_sites_map(sites_tsv: Path) -> Dict[str, Tuple[str, str]]:
    """Map sample_id -> (trip, site) from /project/k10018/empty-quarter/meta/sites.tsv."""
    mapping: Dict[str, Tuple[str, str]] = {}
    with open(sites_tsv, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample_id = (row.get("sample_id") or row["sample"]).strip()   # real sites.tsv column is "sample"
            trip = row["trip"].strip()
            site = row["site"].strip()
            mapping[sample_id] = (trip, site)
    return mapping


def load_rakpr_inputs_map(rakpr_dir: Path) -> Dict[str, str]:
    """Map culture sample_id -> raw assembly FASTA path from RAK-PR inputs.tsv if present."""
    inputs_tsv = rakpr_dir / "inputs.tsv"
    mapping: Dict[str, str] = {}
    if not inputs_tsv.is_file():
        return mapping
    with open(inputs_tsv, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2 or parts[0] in ("sample", "id", "sample_id"):
                continue
            sample_id, raw_path = parts[0].strip(), parts[1].strip()
            staged_candidate = rakpr_dir / "assemblies" / os.path.basename(raw_path)
            if staged_candidate.is_file():
                mapping[sample_id] = str(staged_candidate)
            else:
                mapping[sample_id] = raw_path
    return mapping


def build_eq_records(
    eq_per_bin_tsv: Path,
    sites_tsv: Path,
    sites_tar_dir: Path,
    rakpr_dir: Path,
) -> List[Dict[str, str]]:
    """Build manifest records for all 17,055 Empty Quarter genomes."""
    sites_map = load_sites_map(sites_tsv)
    rakpr_inputs = load_rakpr_inputs_map(rakpr_dir)

    records: List[Dict[str, str]] = []
    with open(eq_per_bin_tsv, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            mag_id = row["mag_id"].strip()
            dataset = row["dataset"].strip()
            sample = row["sample"].strip()
            bin_name = row["bin"].strip()
            domain = row["domain"].strip()
            protein_source = row["protein_source"].strip()
            expected_n_proteins = str(int(row["n_proteins"].strip()))
            gapseq_source = row.get("gapseq_source", "").strip()

            if protein_source not in ("bakta", "pyrodigal"):
                raise ValueError(
                    "Unsupported protein_source for " + mag_id + ": " + protein_source
                )
            if domain not in ("Bacteria", "Archaea"):
                raise ValueError("Unsupported domain for " + mag_id + ": " + domain)

            rec: Dict[str, str] = {
                "batch_id": "",
                "mag_id": mag_id,
                "dataset": dataset,
                "sample": sample,
                "bin": bin_name,
                "domain": domain,
                "protein_source": protein_source,
                "expected_n_proteins": expected_n_proteins,
                "gapseq_source_ibex": gapseq_source,
                "source_kind": "",
                "tar_path": "",
                "bin_fasta_member": "",
                "bakta_prefix": "",
                "direct_faa_path": "",
                "direct_gff3_path": "",
                "direct_fna_path": "",
                "expected_faa_sha256": "",
            }

            if dataset == "soil_bin":
                if sample not in sites_map:
                    raise KeyError(
                        "Soil sample " + sample + " (mag " + mag_id + ") not found in " + str(sites_tsv)
                    )
                trip, site = sites_map[sample]
                rel_sample = trip + "/" + site + "/" + sample
                rec["source_kind"] = "site_tar"
                rec["tar_path"] = str(sites_tar_dir / trip / site / (sample + ".tar"))
                rec["bin_fasta_member"] = (
                    rel_sample + "/wgs/binning/single_sample/output_bins/" + bin_name + ".fa.gz"
                )
                if protein_source == "bakta":
                    rec["bakta_prefix"] = (
                        rel_sample + "/wgs/annotation/single_sample/bakta/" + bin_name + "/" + bin_name
                    )
            else:
                rec["source_kind"] = "rakpr"
                if protein_source == "bakta":
                    bdir = rakpr_dir / "results" / "annotation" / "bakta" / bin_name
                    rec["direct_faa_path"] = str(bdir / (bin_name + ".faa"))
                    rec["direct_gff3_path"] = str(bdir / (bin_name + ".gff3"))
                    rec["direct_fna_path"] = str(bdir / (bin_name + ".fna"))
                else:
                    # RAK-PR inputs.tsv is keyed by genome_id (= bin name, e.g. C-1.1), not by sample (C-1).
                    asm_path = rakpr_inputs.get(
                        bin_name,
                        str(rakpr_dir / "assemblies" / (bin_name + ".fa")),
                    )
                    rec["direct_fna_path"] = asm_path

            records.append(rec)

    return records


def assign_eq_batches(
    records: List[Dict[str, str]],
    batch_size: int = 500,
) -> List[List[Dict[str, str]]]:
    """Group soil_bin genomes by sample into batches (~batch_size), then pack culture genomes."""
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")

    soil_by_sample: "OrderedDict[str, List[Dict[str, str]]]" = OrderedDict()
    culture_records: List[Dict[str, str]] = []

    for rec in records:
        if rec["dataset"] == "soil_bin":
            sample = rec["sample"]
            if sample not in soil_by_sample:
                soil_by_sample[sample] = []
            soil_by_sample[sample].append(rec)
        else:
            culture_records.append(rec)

    batches: List[List[Dict[str, str]]] = []
    current_batch: List[Dict[str, str]] = []

    for sample, sample_recs in soil_by_sample.items():
        if current_batch and (len(current_batch) + len(sample_recs) > batch_size):
            batches.append(current_batch)
            current_batch = []
        current_batch.extend(sample_recs)

    if current_batch:
        batches.append(current_batch)
        current_batch = []

    for rec in culture_records:
        current_batch.append(rec)
        if len(current_batch) >= batch_size:
            batches.append(current_batch)
            current_batch = []

    if current_batch:
        batches.append(current_batch)

    for idx, batch in enumerate(batches):
        batch_id = "batch_{0:03d}".format(idx)
        for rec in batch:
            rec["batch_id"] = batch_id

    return batches


def build_benchmark_records(
    price_targets_json: Path,
    price_genomes_dir: Path,
    wp4_search_dir: Path,
    peers_tsv: Optional[Path] = None,
    price_layout_report: Optional[Path] = None,
) -> List[Dict[str, str]]:
    """Build manifest records for the 22 Price + 18 cohort + 105 peer benchmark genomes."""
    records: List[Dict[str, str]] = []

    # 1. Price 22 genomes
    if not price_targets_json.is_file():
        raise FileNotFoundError("Price targets.json not found: " + str(price_targets_json))
    with open(price_targets_json, "r", encoding="utf-8") as fh:
        targets_data = json.load(fh)

    if "aliases" in targets_data and isinstance(targets_data["aliases"], dict):
        price_assemblies = sorted(targets_data["aliases"].keys())
    elif "alias_resolution" in targets_data and isinstance(targets_data["alias_resolution"], dict):
        price_assemblies = sorted(targets_data["alias_resolution"].keys())
    elif "targets" in targets_data and isinstance(targets_data["targets"], dict):
        price_assemblies = sorted(targets_data["targets"].keys())
    else:
        raise ValueError("Cannot extract Price assembly IDs from " + str(price_targets_json))

    layout_info: Dict[str, dict] = {}
    if price_layout_report is None:
        candidate = price_targets_json.parent / "layout" / "report.json"
        if candidate.is_file():
            price_layout_report = candidate
    if price_layout_report is not None and price_layout_report.is_file():
        with open(price_layout_report, "r", encoding="utf-8") as fh:
            rep = json.load(fh)
            layout_info = rep.get("assemblies", {})

    for asm in price_assemblies:
        info = layout_info.get(asm, {})
        faa_path = price_genomes_dir / asm / "protein.faa"
        gff_path = price_genomes_dir / asm / "genomic.gff"
        fna_path = price_genomes_dir / asm / "genomic.fna"
        n_prot = info.get("proteins", 0)
        if not n_prot and faa_path.is_file():
            n_prot = count_fasta_records(faa_path)
        sha256 = info.get("input_hashes", {}).get("protein.faa", "")

        records.append(
            {
                "batch_id": "",
                "mag_id": asm,
                "dataset": "price_matrix",
                "sample": asm,
                "bin": asm,
                "domain": "Bacteria",
                "protein_source": "ncbi_proteome",
                "expected_n_proteins": str(int(n_prot)),
                "gapseq_source_ibex": "benchmark_target",
                "source_kind": "direct_files",
                "tar_path": "",
                "bin_fasta_member": "",
                "bakta_prefix": "",
                "direct_faa_path": str(faa_path),
                "direct_gff3_path": str(gff_path) if gff_path.is_file() else str(gff_path),
                "direct_fna_path": str(fna_path) if fna_path.is_file() else "",
                "expected_faa_sha256": sha256,
            }
        )

    # 2. Recurrence cohort (18) and peers (105)
    expected_peers: Optional[set] = None
    if peers_tsv is not None and peers_tsv.is_file():
        expected_peers = set()
        with open(peers_tsv, "r", encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                expected_peers.add(row["peer_accession"].strip())

    seen_cohort: set = set()
    seen_peers: set = set()

    for set_name, dataset_label, seen_set in (
        ("cohort", "recurrence_cohort", seen_cohort),
        ("peers", "recurrence_peer", seen_peers),
    ):
        for sk_dir, domain_label in (("archaea", "Archaea"), ("bacteria", "Bacteria")):
            manifest_path = wp4_search_dir / set_name / sk_dir / "manifest.json"
            proteomes_dir = wp4_search_dir / set_name / sk_dir / "proteomes"
            if manifest_path.is_file():
                with open(manifest_path, "r", encoding="utf-8") as fh:
                    mdata = json.load(fh)
                for g in sorted(mdata.get("genomes", []), key=lambda x: x["genome_id"]):
                    gid = g["genome_id"].strip()
                    if gid in seen_set:
                        continue
                    seen_set.add(gid)
                    local_faa = proteomes_dir / (gid + ".faa")
                    faa_str = str(local_faa) if local_faa.is_file() else g.get("proteome", str(local_faa))
                    records.append(
                        {
                            "batch_id": "",
                            "mag_id": gid,
                            "dataset": dataset_label,
                            "sample": set_name,
                            "bin": gid,
                            "domain": g.get("taxonomy", domain_label),
                            "protein_source": "ncbi_proteome",
                            "expected_n_proteins": str(int(g.get("proteins", 0))),
                            "gapseq_source_ibex": "benchmark_target",
                            "source_kind": "direct_files",
                            "tar_path": "",
                            "bin_fasta_member": "",
                            "bakta_prefix": "",
                            "direct_faa_path": faa_str,
                            "direct_gff3_path": "",
                            "direct_fna_path": "",
                            "expected_faa_sha256": g.get("sha256", ""),
                        }
                    )
            elif proteomes_dir.is_dir():
                for faa_file in sorted(proteomes_dir.glob("*.faa")):
                    gid = faa_file.stem
                    if gid in seen_set:
                        continue
                    seen_set.add(gid)
                    n_prot = count_fasta_records(faa_file)
                    records.append(
                        {
                            "batch_id": "",
                            "mag_id": gid,
                            "dataset": dataset_label,
                            "sample": set_name,
                            "bin": gid,
                            "domain": domain_label,
                            "protein_source": "ncbi_proteome",
                            "expected_n_proteins": str(int(n_prot)),
                            "gapseq_source_ibex": "benchmark_target",
                            "source_kind": "direct_files",
                            "tar_path": "",
                            "bin_fasta_member": "",
                            "bakta_prefix": "",
                            "direct_faa_path": str(faa_file),
                            "direct_gff3_path": "",
                            "direct_fna_path": "",
                            "expected_faa_sha256": "",
                        }
                    )

    if expected_peers is not None and seen_peers and seen_peers != expected_peers:
        missing = sorted(expected_peers - seen_peers)
        extra = sorted(seen_peers - expected_peers)
        raise ValueError(
            "Peer set mismatch against peers.tsv: missing="
            + str(missing)
            + " extra="
            + str(extra)
        )

    return records


def assign_fixed_batches(
    records: List[Dict[str, str]],
    batch_size: int = 500,
) -> List[List[Dict[str, str]]]:
    """Pack records sequentially into batches of size up to batch_size."""
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    batches: List[List[Dict[str, str]]] = []
    for i in range(0, len(records), batch_size):
        chunk = records[i : i + batch_size]
        batch_id = "batch_{0:03d}".format(len(batches))
        for rec in chunk:
            rec["batch_id"] = batch_id
        batches.append(chunk)
    return batches


def write_manifests(
    batches: List[List[Dict[str, str]]],
    work_root: Path,
) -> Tuple[Path, Path]:
    """Write manifest_all.tsv, batches_index.tsv, and per-batch batch_manifest.tsv files."""
    work_root.mkdir(parents=True, exist_ok=True)
    batches_dir = work_root / "batches"
    batches_dir.mkdir(parents=True, exist_ok=True)

    all_manifest_path = work_root / "manifest_all.tsv"
    index_path = work_root / "batches_index.tsv"

    with open(all_manifest_path, "w", encoding="utf-8", newline="") as all_fh, open(
        index_path, "w", encoding="utf-8", newline=""
    ) as idx_fh:
        all_writer = csv.DictWriter(all_fh, fieldnames=MANIFEST_COLUMNS, delimiter="\t")
        all_writer.writeheader()

        idx_writer = csv.DictWriter(idx_fh, fieldnames=BATCH_INDEX_COLUMNS, delimiter="\t")
        idx_writer.writeheader()

        for batch_idx, batch in enumerate(batches):
            batch_id = "batch_{0:03d}".format(batch_idx)
            bdir = batches_dir / batch_id
            bdir.mkdir(parents=True, exist_ok=True)
            bmanifest = bdir / "batch_manifest.tsv"

            with open(bmanifest, "w", encoding="utf-8", newline="") as b_fh:
                b_writer = csv.DictWriter(b_fh, fieldnames=MANIFEST_COLUMNS, delimiter="\t")
                b_writer.writeheader()
                for rec in batch:
                    b_writer.writerow(rec)
                    all_writer.writerow(rec)

            samples = sorted(set(r["sample"] for r in batch))
            tars = sorted(set(r["tar_path"] for r in batch if r["tar_path"]))
            datasets = sorted(set(r["dataset"] for r in batch))
            idx_writer.writerow(
                {
                    "batch_id": batch_id,
                    "batch_index": str(batch_idx),
                    "n_genomes": str(len(batch)),
                    "n_samples": str(len(samples)),
                    "n_tars": str(len(tars)),
                    "n_bakta": str(sum(1 for r in batch if r["protein_source"] == "bakta")),
                    "n_pyrodigal": str(sum(1 for r in batch if r["protein_source"] == "pyrodigal")),
                    "n_ncbi_proteome": str(
                        sum(1 for r in batch if r["protein_source"] == "ncbi_proteome")
                    ),
                    "n_archaea": str(sum(1 for r in batch if r["domain"] == "Archaea")),
                    "n_bacteria": str(sum(1 for r in batch if r["domain"] == "Bacteria")),
                    "datasets": ",".join(datasets),
                    "manifest_path": str(bmanifest),
                }
            )

    return all_manifest_path, index_path


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Build per-batch manifests for Empty Quarter or benchmark genomes on Shaheen."
    )
    ap.add_argument(
        "--mode",
        choices=("eq", "benchmark"),
        required=True,
        help="Manifest mode: 'eq' (17,055 Empty Quarter genomes) or 'benchmark' (145 benchmark genomes).",
    )
    ap.add_argument(
        "--work-root",
        type=Path,
        required=True,
        help="Output workspace directory (holds manifest_all.tsv, batches_index.tsv, batches/batch_NNN/).",
    )
    ap.add_argument(
        "--batch-size",
        type=int,
        default=int(os.environ.get("BATCH_SIZE", "500")),
        help="Target number of genomes per batch (default: 500).",
    )
    # EQ mode arguments
    ap.add_argument(
        "--eq-per-bin-tsv",
        type=Path,
        default=Path(
            os.environ.get(
                "EQ_PER_BIN_TSV",
                str(Path(__file__).resolve().parent.parent / "results" / "eq_mags_per_bin.tsv"),
            )
        ),
    )
    ap.add_argument(
        "--sites-tsv",
        type=Path,
        default=Path(os.environ.get("EQ_SITES_META", "/project/k10018/empty-quarter/meta/sites.tsv")),
    )
    ap.add_argument(
        "--sites-tar-dir",
        type=Path,
        default=Path(os.environ.get("EQ_SITES_TAR_DIR", "/project/k10018/empty-quarter/data/sites-tar")),
    )
    ap.add_argument(
        "--rakpr-dir",
        type=Path,
        default=Path(os.environ.get("EQ_RAKPR_DIR", "/project/k10018/empty-quarter/data/rakpr")),
    )
    # Benchmark mode arguments
    ap.add_argument(
        "--price-targets-json",
        type=Path,
        default=Path(
            os.environ.get(
                "PRICE_TARGETS_JSON",
                "/project/k10018/gspa-v3/price-full-matrix-v2/targets.json",
            )
        ),
    )
    ap.add_argument(
        "--price-layout-report",
        type=Path,
        default=Path(
            os.environ.get(
                "PRICE_LAYOUT_REPORT",
                "/project/k10018/gspa-v3/price-full-matrix-v2/layout/report.json",
            )
        ),
    )
    ap.add_argument(
        "--price-genomes-dir",
        type=Path,
        default=Path(
            os.environ.get(
                "PRICE_GENOMES_DIR",
                "/project/k10018/gspa-v3/price-raw-sources/genomes",
            )
        ),
    )
    ap.add_argument(
        "--wp4-search-dir",
        type=Path,
        default=Path(
            os.environ.get(
                "WP4_SEARCH_DIR",
                "/scratch/project/k10018/gspa-v3/out/wp4-search-v1",
            )
        ),
    )
    ap.add_argument(
        "--peers-tsv",
        type=Path,
        default=Path(
            os.environ.get(
                "PEERS_TSV",
                str(
                    Path(__file__).resolve().parent.parent.parent
                    / "v3"
                    / "results"
                    / "peer-selection-v1"
                    / "peers.tsv"
                ),
            )
        ),
    )
    return ap.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    if args.mode == "eq":
        records = build_eq_records(
            eq_per_bin_tsv=args.eq_per_bin_tsv,
            sites_tsv=args.sites_tsv,
            sites_tar_dir=args.sites_tar_dir,
            rakpr_dir=args.rakpr_dir,
        )
        batches = assign_eq_batches(records, batch_size=args.batch_size)
    else:
        records = build_benchmark_records(
            price_targets_json=args.price_targets_json,
            price_genomes_dir=args.price_genomes_dir,
            wp4_search_dir=args.wp4_search_dir,
            peers_tsv=args.peers_tsv,
            price_layout_report=args.price_layout_report,
        )
        batches = assign_fixed_batches(records, batch_size=args.batch_size)

    all_path, idx_path = write_manifests(batches, args.work_root)
    print(
        "[build_manifest] mode={0} total_genomes={1} n_batches={2} index={3}".format(
            args.mode, len(records), len(batches), idx_path
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
