#!/usr/bin/env python3
"""Sidecar for genome-level FOSS predictors (phage / prophage / plasmid /
viral-contig detection).

Output schema is the new 6-column genomic-region TSV::

    contig_id<TAB>region_start<TAB>region_end<TAB>region_type<TAB>score<TAB>attributes

Coordinates are 1-based inclusive on the contig. ``attributes`` is a
free-form ``key=val|key=val`` string (e.g.
``completeness=98|taxonomy=Caudoviricetes``).

Manifest format (TSV with header)::

    tag<TAB>genome_fasta<TAB>gff_path<TAB>output_dir

``gff_path`` is required for PhiSpy (per-CDS features); optional for
geNomad / CheckV. Use the literal string ``-`` if omitted.

Predictors
----------
- ``genomad`` — Camargo et al. 2024 (BSD-3-Clause). Convolutional NN +
  marker HMMs. Detects viral and plasmid contigs / regions, with
  per-gene functional annotation.
- ``checkv``  — Nayfach et al. 2021 (BSD-3-Clause). Estimates viral
  contig completeness + contamination via marker-gene HMMs.
- ``phispy``  — Akhter et al. 2012, McNair et al. 2018 (MIT). Random
  forest on codon-usage / strand-bias features for prophage region
  calls.

Usage
-----
::

    run_genomic_predictors.py --predictor genomad \\
        --manifest manifest.tsv --db-path /path/to/genomad_db

    run_genomic_predictors.py --predictor checkv \\
        --manifest manifest.tsv --db-path /path/to/checkv-db

    run_genomic_predictors.py --predictor phispy \\
        --manifest manifest.tsv --trainset data/trainSets/Mycobacterium

Each runner shells out to the upstream tool (typically via Singularity
biocontainer) and rewrites its output into the 6-col TSV.
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterator, Optional

LOG = logging.getLogger("run_genomic_predictors")

OUTPUT_HEADER = ["contig_id", "region_start", "region_end",
                 "region_type", "score", "attributes"]


@dataclass
class ManifestRow:
    tag: str
    genome_fasta: Path
    gff_path: Optional[Path]
    output_dir: Path


def read_manifest(path: Path) -> list[ManifestRow]:
    rows: list[ManifestRow] = []
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_tag = header.index("tag")
            i_fa  = header.index("genome_fasta")
            i_gff = header.index("gff_path")
            i_out = header.index("output_dir")
        except ValueError as exc:
            raise SystemExit(f"manifest missing required column: {exc}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i_tag, i_fa, i_gff, i_out):
                continue
            gff = None if f[i_gff] in ("", "-") else Path(f[i_gff])
            rows.append(ManifestRow(
                tag=f[i_tag],
                genome_fasta=Path(f[i_fa]),
                gff_path=gff,
                output_dir=Path(f[i_out]),
            ))
    return rows


def output_path(row: ManifestRow, predictor: str) -> Path:
    row.output_dir.mkdir(parents=True, exist_ok=True)
    return row.output_dir / f"{row.tag}.{predictor}.genomic.tsv"


def open_output(path: Path):
    fh = path.open("w", newline="")
    writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
    writer.writerow(OUTPUT_HEADER)
    return fh, writer


def fmt_attrs(d: dict[str, str]) -> str:
    return "|".join(f"{k}={v}" for k, v in sorted(d.items()) if v not in (None, ""))


# ---------- geNomad -------------------------------------------------------

def run_genomad(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Run geNomad end-to-end mode and parse its provirus + plasmid summaries.

    Native output (under ``<out>/<tag>_summary/``):
      - ``<tag>_virus_summary.tsv`` — viral contigs/proviruses with
        coordinates, length, and topology
      - ``<tag>_plasmid_summary.tsv`` — plasmid contigs

    For each row, we map:
      - virus_summary → region_type=prophage when has parent contig coords;
        region_type=viral_contig otherwise
      - plasmid_summary → region_type=plasmid
    """
    if not args.db_path:
        raise SystemExit("genomad requires --db-path")
    if not args.genomad_sif and not shutil.which("genomad"):
        raise SystemExit("genomad not on PATH and no --genomad-sif provided")

    for row in rows:
        out = output_path(row, "genomad")
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            if args.genomad_sif:
                in_dir = row.genome_fasta.parent.resolve()
                cmd = ["singularity", "exec",
                       "--bind", f"{in_dir}:/in",
                       "--bind", f"{tmp_path}:/out",
                       "--bind", f"{args.db_path}:/db",
                       args.genomad_sif,
                       "genomad", "end-to-end",
                       "--cleanup",
                       f"/in/{row.genome_fasta.name}",
                       "/out", "/db"]
            else:
                cmd = ["genomad", "end-to-end", "--cleanup",
                       str(row.genome_fasta), str(tmp_path), str(args.db_path)]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                LOG.warning("%s genomad exit %d: %s",
                            row.tag, r.returncode, r.stderr[-500:])
                continue

            stem = row.genome_fasta.stem
            summary = tmp_path / f"{stem}_summary"
            virus_tsv = summary / f"{stem}_virus_summary.tsv"
            plas_tsv  = summary / f"{stem}_plasmid_summary.tsv"

            fh, w = open_output(out)
            n = 0
            if virus_tsv.exists():
                with virus_tsv.open() as fin:
                    header = fin.readline().rstrip("\n").split("\t")
                    h = {c: i for i, c in enumerate(header)}
                    for line in fin:
                        f = line.rstrip("\n").split("\t")
                        if len(f) < len(header):
                            continue
                        seqid = f[h.get("seq_name", 0)]
                        score = float(f[h["virus_score"]]) if "virus_score" in h else 0.0
                        topo = f[h.get("topology", -1)] if "topology" in h else ""
                        coords = f[h.get("coordinates", -1)] if "coordinates" in h else ""
                        length = f[h.get("length", -1)] if "length" in h else ""
                        taxonomy = f[h.get("taxonomy", -1)] if "taxonomy" in h else ""
                        # If ``coordinates`` is set we have an integrated provirus
                        # (e.g. "12345-67890"); otherwise the whole contig is viral.
                        if coords and "-" in coords:
                            try:
                                start, end = (int(x) for x in coords.split("-", 1))
                                contig = seqid.rsplit("|", 1)[0]
                            except ValueError:
                                continue
                            region_type = "prophage"
                        else:
                            try:
                                start = 1
                                end = int(length) if length else 0
                            except ValueError:
                                continue
                            contig = seqid
                            region_type = "viral_contig"
                        attrs = fmt_attrs({
                            "topology": topo,
                            "taxonomy": taxonomy,
                            "length": length,
                        })
                        if score >= args.min_score:
                            w.writerow([contig, start, end, region_type,
                                        f"{score:.4f}", attrs])
                            n += 1
            if plas_tsv.exists():
                with plas_tsv.open() as fin:
                    header = fin.readline().rstrip("\n").split("\t")
                    h = {c: i for i, c in enumerate(header)}
                    for line in fin:
                        f = line.rstrip("\n").split("\t")
                        if len(f) < len(header):
                            continue
                        seqid = f[h.get("seq_name", 0)]
                        score = float(f[h["plasmid_score"]]) if "plasmid_score" in h else 0.0
                        topo = f[h.get("topology", -1)] if "topology" in h else ""
                        length = f[h.get("length", -1)] if "length" in h else ""
                        if score < args.min_score:
                            continue
                        try:
                            end = int(length) if length else 0
                        except ValueError:
                            end = 0
                        attrs = fmt_attrs({"topology": topo, "length": length})
                        w.writerow([seqid, 1, end, "plasmid",
                                    f"{score:.4f}", attrs])
                        n += 1
            fh.close()
        LOG.info("%s genomad: %d genomic regions -> %s", row.tag, n, out)


# ---------- CheckV --------------------------------------------------------

def run_checkv(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Run CheckV end-to-end and parse quality_summary.tsv.

    CheckV emits (in ``<out>/quality_summary.tsv``):
      contig_id, contig_length, provirus, proviral_length, gene_count,
      viral_genes, host_genes, checkv_quality, miuvig_quality,
      completeness, completeness_method, contamination, kmer_freq, warnings

    We map to:
      region_type='viral_contig' (whole-contig classification)
      score = completeness / 100
      attributes = checkv_quality, miuvig_quality, completeness, contamination
    """
    if not args.db_path:
        raise SystemExit("checkv requires --db-path")
    if not args.checkv_sif and not shutil.which("checkv"):
        raise SystemExit("checkv not on PATH and no --checkv-sif provided")

    for row in rows:
        out = output_path(row, "checkv")
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            if args.checkv_sif:
                in_dir = row.genome_fasta.parent.resolve()
                cmd = ["singularity", "exec",
                       "--bind", f"{in_dir}:/in",
                       "--bind", f"{tmp_path}:/out",
                       "--bind", f"{args.db_path}:/db",
                       args.checkv_sif,
                       "checkv", "end_to_end",
                       f"/in/{row.genome_fasta.name}",
                       "/out", "-d", "/db", "-t", str(args.threads)]
            else:
                cmd = ["checkv", "end_to_end",
                       str(row.genome_fasta), str(tmp_path),
                       "-d", str(args.db_path), "-t", str(args.threads)]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                LOG.warning("%s checkv exit %d: %s",
                            row.tag, r.returncode, r.stderr[-500:])
                continue

            qs = tmp_path / "quality_summary.tsv"
            if not qs.exists():
                LOG.warning("%s: no checkv quality_summary.tsv", row.tag)
                continue

            fh, w = open_output(out)
            n = 0
            with qs.open() as fin:
                header = fin.readline().rstrip("\n").split("\t")
                h = {c: i for i, c in enumerate(header)}
                for line in fin:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < len(header):
                        continue
                    contig = f[h.get("contig_id", 0)]
                    length = f[h.get("contig_length", -1)]
                    completeness = f[h.get("completeness", -1)] or "0"
                    try:
                        end = int(length)
                        comp = float(completeness) if completeness != "NA" else 0.0
                    except ValueError:
                        continue
                    norm = comp / 100.0
                    if norm < args.min_score:
                        continue
                    attrs = fmt_attrs({
                        "checkv_quality": f[h.get("checkv_quality", -1)],
                        "miuvig_quality": f[h.get("miuvig_quality", -1)],
                        "completeness":   completeness,
                        "contamination":  f[h.get("contamination", -1)],
                        "viral_genes":    f[h.get("viral_genes", -1)],
                    })
                    w.writerow([contig, 1, end, "viral_contig",
                                f"{norm:.4f}", attrs])
                    n += 1
            fh.close()
        LOG.info("%s checkv: %d viral contigs -> %s", row.tag, n, out)


# ---------- PhiSpy --------------------------------------------------------

def run_phispy(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Run PhiSpy and parse its prophage_coordinates.tsv.

    Native output (in ``<out>/prophage_coordinates.tsv``):
      pp_number, contig, start, end, ...

    PhiSpy needs a GenBank file as input (not raw FASTA). The sidecar
    expects ``<row.genome_fasta>`` to point at a GBK file when the user
    enables PhiSpy. (geNomad and CheckV take FASTA; PhiSpy takes GBK.)
    """
    if not args.phispy_sif and not shutil.which("PhiSpy.py"):
        raise SystemExit("PhiSpy.py not on PATH and no --phispy-sif provided")

    trainset = args.phispy_trainset or "data/trainSets/genericAll.txt"

    for row in rows:
        out = output_path(row, "phispy")
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            if args.phispy_sif:
                in_dir = row.genome_fasta.parent.resolve()
                cmd = ["singularity", "exec",
                       "--bind", f"{in_dir}:/in",
                       "--bind", f"{tmp_path}:/out",
                       args.phispy_sif,
                       "PhiSpy.py",
                       "-o", "/out",
                       f"/in/{row.genome_fasta.name}"]
            else:
                cmd = ["PhiSpy.py", "-o", str(tmp_path), str(row.genome_fasta)]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                LOG.warning("%s phispy exit %d: %s",
                            row.tag, r.returncode, r.stderr[-500:])
                continue

            coord_tsv = tmp_path / "prophage_coordinates.tsv"
            if not coord_tsv.exists():
                LOG.info("%s: no PhiSpy prophage_coordinates.tsv (no prophages)", row.tag)
                # still write empty header so downstream is consistent
                fh, _ = open_output(out)
                fh.close()
                continue

            fh, w = open_output(out)
            n = 0
            with coord_tsv.open() as fin:
                for line in fin:
                    if not line.strip() or line.startswith("#"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 4:
                        continue
                    # PhiSpy doesn't emit header; columns:
                    # pp_number, contig, start, end, [attL_start..]
                    try:
                        contig = f[1]
                        start = int(f[2])
                        end = int(f[3])
                    except ValueError:
                        continue
                    score = 1.0  # PhiSpy doesn't emit a confidence score; use 1
                    if score < args.min_score:
                        continue
                    w.writerow([contig, start, end, "prophage",
                                f"{score:.4f}", f"pp_number={f[0]}"])
                    n += 1
            fh.close()
        LOG.info("%s phispy: %d prophage regions -> %s", row.tag, n, out)


# ---------- registry + CLI ------------------------------------------------

RUNNERS: dict[str, Callable] = {
    "genomad": run_genomad,
    "checkv":  run_checkv,
    "phispy":  run_phispy,
}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--predictor", required=True, choices=sorted(RUNNERS.keys()))
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--min-score", type=float, default=0.5)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--db-path", type=Path,
                    help="genomad / checkv: path to the model DB directory")
    ap.add_argument("--genomad-sif",
                    help="Singularity image for geNomad (optional; biocontainer recommended)")
    ap.add_argument("--checkv-sif",
                    help="Singularity image for CheckV (optional)")
    ap.add_argument("--phispy-sif",
                    help="Singularity image for PhiSpy (optional)")
    ap.add_argument("--phispy-trainset",
                    help="PhiSpy trainset file path (default: genericAll)")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    rows = read_manifest(args.manifest)
    LOG.info("manifest: %d rows", len(rows))
    LOG.info("predictor=%s", args.predictor)
    RUNNERS[args.predictor](rows, args)
    LOG.info("done")


if __name__ == "__main__":
    main()
