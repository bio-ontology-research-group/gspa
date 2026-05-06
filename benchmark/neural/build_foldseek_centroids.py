#!/usr/bin/env python3
"""Build a FoldSeek function-centroid database from SwissProt + AlphaFoldDB.

For each GO term / EC number that has ≥ ``--min-class-size`` SwissProt
entries with an AlphaFold model, pick a *medoid* protein and stage its
PDB into a per-term file named::

    GO:0003824_medoid_<ACC>.pdb
    EC:1.1.1.1_medoid_<ACC>.pdb

Then run ``foldseek createdb`` over all staged files to produce a FoldSeek
database whose entry IDs match the
``(GO:\\d+|EC:[\\d.\\-]+)_medoid_<ACC>`` shape expected by the JVM-side
``FoldSeekPredictor`` running in ``CentroidMode``.

Medoid selection
----------------
*v1 (this script)*: deterministic representative selection — the
accession with the smallest string among the class members that has an
AlphaFold model available. This is a placeholder that's reproducible
without requiring intra-class all-vs-all structural alignment.

*v2 (TODO)*: real structural medoid via intra-class foldseek all-vs-all
alignment + pick the protein whose mean TM-score against every other
class member is largest. Strictly better; strictly slower (O(N²) per
class). Leaving the medoid selector behind a function hook so swap-in
is one-line.

Inputs
------
--swissprot-tsv      TSV with columns ``accession`` ``go_terms`` ``ec_numbers``
                     semicolon-separated lists. One row per SwissProt entry.
                     Build with ``make_swissprot_tsv.py`` or lift from the
                     goa / uniprot flat files with awk.
--alphafold-dir      Directory of AlphaFold models. Expected filename
                     pattern: ``AF-<ACC>-F1-model_v*.pdb`` (single-chain)
                     or ``AF-<ACC>-F1-model_v*.pdb.gz`` (gzipped).
--min-class-size     Minimum class cardinality to build a medoid for (default 10).
--outdir             Output directory.
--foldseek           Path to the foldseek binary (default: ``foldseek`` on PATH).
--skip-ec / --skip-go  Restrict to one annotation type (defaults to both).

Outputs
-------
``<outdir>/stage/`` — renamed PDBs (one per centroid).
``<outdir>/centroid_db/*`` — FoldSeek database files.
``<outdir>/centroid_metadata.tsv`` — per-centroid metadata: ``centroid_id``,
``term``, ``annotation_type``, ``n_class_members``, ``rep_accession``.
``<outdir>/build.log`` — log of the run.
"""
from __future__ import annotations

import argparse
import gzip
import logging
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterator, Optional

LOG = logging.getLogger("build_foldseek_centroids")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--swissprot-tsv", type=Path, required=True,
                    help="TSV: accession\\tgo_terms;...\\tec_numbers;...")
    ap.add_argument("--alphafold-dir", type=Path, required=True,
                    help="Directory with AF-<ACC>-F1-model_v*.pdb[.gz] files")
    ap.add_argument("--min-class-size", type=int, default=10,
                    help="Minimum class size; smaller classes are skipped (default 10)")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--foldseek", default="foldseek",
                    help="Path to the foldseek binary (default: 'foldseek' on PATH)")
    ap.add_argument("--skip-go", action="store_true",
                    help="Build only EC centroids, no GO")
    ap.add_argument("--skip-ec", action="store_true",
                    help="Build only GO centroids, no EC")
    ap.add_argument("--dry-run", action="store_true",
                    help="Stage PDBs and write metadata but do not call foldseek createdb")
    ap.add_argument("--threads", type=int, default=4)
    return ap.parse_args()


def configure_logging(log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler(sys.stderr),
        ],
    )


def iter_swissprot(path: Path) -> Iterator[tuple[str, list[str], list[str]]]:
    """Yield (accession, go_terms, ec_numbers) rows from the SwissProt TSV."""
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            acc_idx = header.index("accession")
            go_idx = header.index("go_terms")
            ec_idx = header.index("ec_numbers")
        except ValueError as exc:
            raise SystemExit(f"swissprot-tsv header missing required column: {exc}")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(acc_idx, go_idx, ec_idx):
                continue
            acc = fields[acc_idx].strip()
            gos = [t.strip() for t in fields[go_idx].split(";") if t.strip()]
            ecs = [t.strip() for t in fields[ec_idx].split(";") if t.strip()]
            if not acc:
                continue
            yield acc, gos, ecs


def find_alphafold_pdb(af_dir: Path, accession: str) -> Optional[Path]:
    """Return the AlphaFold PDB for the given accession, or None."""
    for suffix in (".pdb", ".pdb.gz"):
        # Match any model version (v1, v2, v3, v4, …)
        hits = sorted(af_dir.glob(f"AF-{accession}-F1-model_v*{suffix}"))
        if hits:
            # If multiple versions exist, take the highest one
            return hits[-1]
    return None


def group_by_term(
    rows: Iterator[tuple[str, list[str], list[str]]],
    do_go: bool,
    do_ec: bool,
) -> dict[str, list[str]]:
    """Build {term -> [accession...]} from the SwissProt rows."""
    groups: dict[str, set[str]] = defaultdict(set)
    for acc, gos, ecs in rows:
        if do_go:
            for go in gos:
                if go.startswith("GO:"):
                    groups[go].add(acc)
        if do_ec:
            for ec in ecs:
                # Skip EC numbers with hyphen placeholders (partial classifications)
                if not ec.startswith("EC:"):
                    ec = f"EC:{ec}"
                if "-" in ec.split(":")[1]:
                    continue
                groups[ec].add(acc)
    return {k: sorted(v) for k, v in groups.items()}


def pick_medoid(term: str, members: list[str], af_dir: Path) -> Optional[tuple[str, Path]]:
    """v1: pick the alphabetically-first member that has an AF model."""
    for acc in members:  # already sorted
        pdb = find_alphafold_pdb(af_dir, acc)
        if pdb is not None:
            return acc, pdb
    return None


def decompress_if_needed(src: Path, dst: Path) -> None:
    if src.suffix == ".gz":
        with gzip.open(src, "rb") as fin, dst.open("wb") as fout:
            shutil.copyfileobj(fin, fout)
    else:
        shutil.copy2(src, dst)


def stage_centroids(
    groups: dict[str, list[str]],
    af_dir: Path,
    min_class_size: int,
    stage_dir: Path,
) -> list[dict[str, object]]:
    """Stage one PDB per term that meets the min class size; return metadata rows."""
    stage_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for term in sorted(groups):
        members = groups[term]
        if len(members) < min_class_size:
            continue
        picked = pick_medoid(term, members, af_dir)
        if picked is None:
            LOG.warning("term %s has %d members but no AlphaFold model available", term, len(members))
            continue
        acc, src_pdb = picked
        centroid_id = f"{term}_medoid_{acc}"
        dst = stage_dir / f"{centroid_id}.pdb"
        try:
            decompress_if_needed(src_pdb, dst)
        except OSError as exc:
            LOG.warning("failed to stage %s (%s): %s", term, src_pdb, exc)
            continue
        ann_type = "GO" if term.startswith("GO:") else "EC"
        rows.append({
            "centroid_id": centroid_id,
            "term": term,
            "annotation_type": ann_type,
            "n_class_members": len(members),
            "rep_accession": acc,
            "src_pdb": str(src_pdb),
        })
    return rows


def write_metadata(rows: list[dict[str, object]], path: Path) -> None:
    headers = ["centroid_id", "term", "annotation_type", "n_class_members", "rep_accession", "src_pdb"]
    with path.open("w") as fh:
        fh.write("\t".join(headers) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[h]) for h in headers) + "\n")


def run_foldseek_createdb(foldseek: str, stage_dir: Path, db_path: Path, threads: int) -> None:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    pdbs = sorted(stage_dir.glob("*.pdb"))
    if not pdbs:
        raise SystemExit("no staged PDBs; nothing to build")
    cmd = [foldseek, "createdb", *[str(p) for p in pdbs], str(db_path), "--threads", str(threads)]
    LOG.info("running foldseek createdb (%d entries) → %s", len(pdbs), db_path)
    subprocess.run(cmd, check=True)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    configure_logging(args.outdir / "build.log")

    if args.skip_go and args.skip_ec:
        raise SystemExit("--skip-go and --skip-ec together leave nothing to build")

    LOG.info("reading SwissProt TSV: %s", args.swissprot_tsv)
    groups = group_by_term(
        iter_swissprot(args.swissprot_tsv),
        do_go=not args.skip_go,
        do_ec=not args.skip_ec,
    )
    LOG.info("grouped %d terms", len(groups))
    n_above_threshold = sum(1 for m in groups.values() if len(m) >= args.min_class_size)
    LOG.info("  %d terms meet min-class-size=%d", n_above_threshold, args.min_class_size)

    stage_dir = args.outdir / "stage"
    LOG.info("staging centroids into %s", stage_dir)
    rows = stage_centroids(groups, args.alphafold_dir, args.min_class_size, stage_dir)
    LOG.info("staged %d centroids", len(rows))

    meta_path = args.outdir / "centroid_metadata.tsv"
    write_metadata(rows, meta_path)
    LOG.info("wrote %s", meta_path)

    if args.dry_run:
        LOG.info("--dry-run: skipping foldseek createdb")
        return

    db_path = args.outdir / "centroid_db" / "centroid_db"
    run_foldseek_createdb(args.foldseek, stage_dir, db_path, args.threads)
    LOG.info("done: centroid DB at %s", db_path)


if __name__ == "__main__":
    main()
