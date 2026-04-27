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
                       "--threads", str(args.threads),
                       f"/in/{row.genome_fasta.name}",
                       "/out", "/db"]
            else:
                cmd = ["genomad", "end-to-end", "--cleanup",
                       "--threads", str(args.threads),
                       str(row.genome_fasta), str(tmp_path), str(args.db_path)]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                LOG.warning("%s genomad exit %d: %s",
                            row.tag, r.returncode, r.stderr[-3000:])
                LOG.warning("%s genomad stdout tail: %s",
                            row.tag, r.stdout[-2000:])
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


# ---------- VirSorter2 ----------------------------------------------------

def run_virsorter2(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Run VirSorter2 (Guo et al. 2021, GPL-2). Biocontainer:
    ``quay.io/biocontainers/virsorter:2.2.4--pyhdfd78af1_2``.

    CLI::

        virsorter run --keep-original-seq -i contigs.fa -w out \\
            --include-groups dsDNAphage,ssDNA,NCLDV,RNA,lavidaviridae \\
            --min-length 1500 --min-score 0.5 -j N all

    Output parse: ``out/final-viral-boundary.tsv`` columns
    ``seqname, trim_bp_start, trim_bp_end, trim_pr, group, shape, partial``
    where ``trim_pr`` is already 0-1, ``partial==1`` ⇒ prophage,
    ``partial==0`` ⇒ whole-contig viral.
    """
    if not args.virsorter2_sif and not shutil.which("virsorter"):
        raise SystemExit("virsorter2 not on PATH and no --virsorter2-sif")
    if not args.db_path:
        raise SystemExit("virsorter2 requires --db-path (~11 GB DB)")

    for row in rows:
        out = output_path(row, "virsorter2")
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            in_dir = row.genome_fasta.parent.resolve()
            if args.virsorter2_sif:
                cmd = ["singularity", "exec",
                       "--bind", f"{in_dir}:/in",
                       "--bind", f"{tmp_path}:/out",
                       "--bind", f"{args.db_path}:/db",
                       args.virsorter2_sif,
                       "virsorter", "run",
                       "--keep-original-seq",
                       "-i", f"/in/{row.genome_fasta.name}",
                       "-w", "/out/vs2",
                       "-d", "/db",
                       "--min-length", "1500",
                       "--min-score", str(args.min_score),
                       "-j", str(args.threads),
                       "all"]
            else:
                cmd = ["virsorter", "run",
                       "--keep-original-seq",
                       "-i", str(row.genome_fasta),
                       "-w", str(tmp_path / "vs2"),
                       "-d", str(args.db_path),
                       "--min-length", "1500",
                       "--min-score", str(args.min_score),
                       "-j", str(args.threads),
                       "all"]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                LOG.warning("%s virsorter2 exit %d: %s",
                            row.tag, r.returncode, r.stderr[-2000:])
                continue

            boundary = tmp_path / "vs2" / "final-viral-boundary.tsv"
            fh, w = open_output(out)
            n = 0
            if boundary.exists():
                with boundary.open() as fin:
                    header = fin.readline().rstrip("\n").split("\t")
                    h = {c: i for i, c in enumerate(header)}
                    for line in fin:
                        f = line.rstrip("\n").split("\t")
                        if len(f) < len(header):
                            continue
                        # Strip VS2's "||full"/"||partial" seq-name suffix
                        contig = f[h.get("seqname", 0)].split("||")[0]
                        try:
                            start = int(f[h["trim_bp_start"]])
                            end   = int(f[h["trim_bp_end"]])
                            score = float(f[h["trim_pr"]])
                        except (ValueError, KeyError):
                            continue
                        if score < args.min_score:
                            continue
                        partial = f[h.get("partial", -1)] if "partial" in h else "0"
                        rtype = "prophage" if partial == "1" else "viral_contig"
                        attrs = fmt_attrs({
                            "tool":  "virsorter2",
                            "group": f[h.get("group", -1)] if "group" in h else "",
                            "shape": f[h.get("shape", -1)] if "shape" in h else "",
                        })
                        w.writerow([contig, start, end, rtype,
                                    f"{score:.4f}", attrs])
                        n += 1
            fh.close()
        LOG.info("%s virsorter2: %d genomic regions -> %s", row.tag, n, out)


# ---------- VIBRANT ------------------------------------------------------

# VIBRANT has no per-call probability; map its Quality string to a 0-1 score.
VIBRANT_QUALITY_TO_SCORE = {
    "complete circular":   1.00,
    "high quality draft":  0.90,
    "medium quality draft": 0.70,
    "low quality draft":   0.50,
}


def run_vibrant(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Run VIBRANT (Kieft et al. 2020, GPL-3). Biocontainer:
    ``quay.io/biocontainers/vibrant:1.2.1--hdfd78af_2``.

    CLI::

        VIBRANT_run.py -i contigs.fa -t N -folder out \\
            -d <db>/databases -m <db>/files -no_plot

    Output parse:
      - prophage coords:
        ``out/VIBRANT_<name>/VIBRANT_results_<name>/VIBRANT_integrated_prophage_coordinates_<name>.tsv``
        columns: ``scaffold, fragment, nucleotide start, nucleotide stop, ...``
      - whole-contig viral calls:
        ``out/VIBRANT_<name>/VIBRANT_results_<name>/VIBRANT_genome_quality_<name>.tsv``
        columns: ``scaffold, type, Quality``

    VIBRANT doesn't emit a per-call probability — we map Quality to a
    score via VIBRANT_QUALITY_TO_SCORE.
    """
    if not args.vibrant_sif and not shutil.which("VIBRANT_run.py"):
        raise SystemExit("VIBRANT not on PATH and no --vibrant-sif")
    if not args.db_path:
        raise SystemExit("vibrant requires --db-path (HMM databases dir)")

    for row in rows:
        out = output_path(row, "vibrant")
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            in_dir = row.genome_fasta.parent.resolve()
            if args.vibrant_sif:
                cmd = ["singularity", "exec",
                       "--bind", f"{in_dir}:/in",
                       "--bind", f"{tmp_path}:/out",
                       "--bind", f"{args.db_path}:/db",
                       args.vibrant_sif,
                       "VIBRANT_run.py",
                       "-i",      f"/in/{row.genome_fasta.name}",
                       "-folder", "/out",
                       "-d",      "/db/databases",
                       "-m",      "/db/files",
                       "-t",      str(args.threads),
                       "-no_plot"]
            else:
                cmd = ["VIBRANT_run.py",
                       "-i",      str(row.genome_fasta),
                       "-folder", str(tmp_path),
                       "-d",      str(args.db_path / "databases"),
                       "-m",      str(args.db_path / "files"),
                       "-t",      str(args.threads),
                       "-no_plot"]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                LOG.warning("%s vibrant exit %d: %s",
                            row.tag, r.returncode, r.stderr[-2000:])
                continue

            stem = row.genome_fasta.stem
            results_dir = (tmp_path / f"VIBRANT_{stem}"
                           / f"VIBRANT_results_{stem}")
            prophage_tsv = (results_dir
                / f"VIBRANT_integrated_prophage_coordinates_{stem}.tsv")
            quality_tsv  = (results_dir
                / f"VIBRANT_genome_quality_{stem}.tsv")

            fh, w = open_output(out)
            n = 0
            # Prophage rows (partial-contig integrated prophages)
            if prophage_tsv.exists():
                with prophage_tsv.open() as fin:
                    header = fin.readline().rstrip("\n").split("\t")
                    h = {c: i for i, c in enumerate(header)}
                    for line in fin:
                        f = line.rstrip("\n").split("\t")
                        if len(f) < len(header):
                            continue
                        contig = f[h.get("scaffold", 0)]
                        try:
                            start = int(f[h["nucleotide start"]])
                            end   = int(f[h["nucleotide stop"]])
                        except (ValueError, KeyError):
                            continue
                        # No native score; emit 0.9 (medium-high confidence)
                        score = 0.9
                        if score < args.min_score:
                            continue
                        attrs = fmt_attrs({
                            "tool":     "vibrant",
                            "fragment": f[h["fragment"]] if "fragment" in h else "",
                        })
                        w.writerow([contig, start, end, "prophage",
                                    f"{score:.4f}", attrs])
                        n += 1
            # Whole-contig viral calls (quality-graded)
            if quality_tsv.exists():
                with quality_tsv.open() as fin:
                    header = fin.readline().rstrip("\n").split("\t")
                    h = {c: i for i, c in enumerate(header)}
                    for line in fin:
                        f = line.rstrip("\n").split("\t")
                        if len(f) < len(header):
                            continue
                        contig = f[h.get("scaffold", 0)]
                        quality = f[h.get("Quality", -1)] if "Quality" in h else ""
                        score = VIBRANT_QUALITY_TO_SCORE.get(quality.lower(), 0.5)
                        if score < args.min_score:
                            continue
                        # VIBRANT quality TSV doesn't repeat coords; flag
                        # the whole contig and let downstream consumers
                        # join in a length lookup if needed.
                        attrs = fmt_attrs({
                            "tool":    "vibrant",
                            "quality": quality,
                            "type":    f[h["type"]] if "type" in h else "",
                        })
                        w.writerow([contig, 1, 0, "viral_contig",
                                    f"{score:.4f}", attrs])
                        n += 1
            fh.close()
        LOG.info("%s vibrant: %d genomic regions -> %s", row.tag, n, out)


# ---------- registry + CLI ------------------------------------------------

RUNNERS: dict[str, Callable] = {
    "genomad":    run_genomad,
    "checkv":     run_checkv,
    "phispy":     run_phispy,
    "virsorter2": run_virsorter2,
    "vibrant":    run_vibrant,
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
    ap.add_argument("--virsorter2-sif",
                    help="Singularity image for VirSorter2 (biocontainer recommended)")
    ap.add_argument("--vibrant-sif",
                    help="Singularity image for VIBRANT (biocontainer recommended)")
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
