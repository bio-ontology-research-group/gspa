#!/usr/bin/env python3
"""Sidecar for site-level FOSS protein predictors.

Output schema (5 columns)::

    protein_id<TAB>position<TAB>site_type<TAB>score<TAB>annotation_type

Position is 1-based residue index. ``annotation_type`` is one of
``PTM_SITE``, ``PPI_INTERFACE``.

Predictors
----------
- ``musitedeep`` — PTM sites (phospho-S/T/Y by default; configurable)
                  — duolinwang/MusiteDeep_web (MIT)
- ``scannet``   — PPI interface residues, structure-based
                  — jertubiana/ScanNet (Apache-2.0)

Usage
-----
::

    run_site_predictors.py --predictor musitedeep \\
        --manifest manifest.tsv --residue-types Phosphoserine_Phosphothreonine_Phosphotyrosine

    run_site_predictors.py --predictor scannet \\
        --manifest manifest.tsv --structure-dir /path/to/pdbs
"""
from __future__ import annotations

import argparse
import csv
import logging
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

LOG = logging.getLogger("run_site_predictors")

OUTPUT_HEADER = ["protein_id", "position", "site_type", "score", "annotation_type"]


@dataclass
class ManifestRow:
    tag: str
    fasta_path: Path
    output_dir: Path


def read_manifest(path: Path) -> list[ManifestRow]:
    rows: list[ManifestRow] = []
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {k: header.index(k) for k in ("tag", "fasta_path", "output_dir")}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            rows.append(ManifestRow(
                tag=f[idx["tag"]],
                fasta_path=Path(f[idx["fasta_path"]]),
                output_dir=Path(f[idx["output_dir"]]),
            ))
    return rows


def output_path(row: ManifestRow, predictor: str) -> Path:
    row.output_dir.mkdir(parents=True, exist_ok=True)
    return row.output_dir / f"{row.tag}.{predictor}.tsv"


def open_output(path: Path):
    fh = path.open("w", newline="")
    writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
    writer.writerow(OUTPUT_HEADER)
    return fh, writer


# ---------- MusiteDeep ----------------------------------------------------

def run_musitedeep(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Run MusiteDeep_web's predict_batch.py.

    Output format per residue type: a CSV with columns
    ID, Position, Residue, Score, Cutoff (depends on release).
    """
    if not args.model_dir:
        raise SystemExit("musitedeep requires --model-dir (path to MusiteDeep_web)")

    runner = Path(args.model_dir) / "MusiteDeep" / "predict_multi_batch.py"
    if not runner.exists():
        runner = Path(args.model_dir) / "predict_batch.py"
    if not runner.exists():
        raise SystemExit(f"musitedeep predict script not found in {args.model_dir}")

    for row in rows:
        out = output_path(row, "musitedeep")
        with tempfile.TemporaryDirectory() as tmp:
            cmd = ["python3", str(runner),
                   "-input", str(row.fasta_path),
                   "-output", tmp,
                   "-model-prefix", "MusiteDeep",
                   "-residue-types", args.residue_types]
            subprocess.run(cmd, check=True, cwd=args.model_dir)

            fh, w = open_output(out)
            n = 0
            for csv_file in Path(tmp).glob("*_results.txt"):
                # Format example: protein_id<TAB>pos<TAB>residue<TAB>score<TAB>cutoff<TAB>type
                with csv_file.open() as fin:
                    next(fin, None)  # header
                    for line in fin:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 5:
                            continue
                        try:
                            pid = parts[0]
                            pos = int(parts[1])
                            score = float(parts[3])
                        except ValueError:
                            continue
                        if score < args.min_score:
                            continue
                        site_type = parts[5] if len(parts) >= 6 else "phosphosite"
                        w.writerow([pid, pos, site_type, f"{score:.4f}", "PTM_SITE"])
                        n += 1
            fh.close()
        LOG.info("%s musitedeep: %d sites -> %s", row.tag, n, out)


# ---------- ScanNet -------------------------------------------------------

def run_scannet(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Run ScanNet on per-protein structures.

    For each tag, ScanNet expects a directory of PDB/CIF files under
    ``--structure-dir/<tag>/``. We invoke ScanNet's ``predict_features.py``
    once per protein file. Output: per-residue PPI-interface probability.
    """
    if not args.model_dir:
        raise SystemExit("scannet requires --model-dir (path to ScanNet repo)")
    if not args.structure_dir:
        raise SystemExit("scannet requires --structure-dir")

    runner = Path(args.model_dir) / "predict_features.py"
    if not runner.exists():
        raise SystemExit(f"ScanNet predict_features.py not found at {runner}")

    for row in rows:
        out = output_path(row, "scannet")
        struct_root = Path(args.structure_dir) / row.tag
        if not struct_root.exists():
            LOG.warning("no structure dir for %s at %s; skipping",
                        row.tag, struct_root)
            continue
        pdbs = sorted(list(struct_root.glob("*.pdb")) +
                       list(struct_root.glob("*.cif")))
        if not pdbs:
            LOG.warning("%s: no PDB/CIF in %s", row.tag, struct_root)
            continue

        fh, w = open_output(out)
        n = 0
        for pdb in pdbs:
            pid = pdb.stem
            with tempfile.TemporaryDirectory() as tmp:
                cmd = ["python3", str(runner),
                       "--pdb", str(pdb),
                       "--mode", "interface",
                       "--out", tmp]
                try:
                    subprocess.run(cmd, check=True, cwd=args.model_dir)
                except subprocess.CalledProcessError as exc:
                    LOG.warning("scannet failed on %s: %s", pid, exc)
                    continue
                # Output: a per-residue TSV (residue_index, score)
                for tsv in Path(tmp).glob("*.tsv"):
                    with tsv.open() as fin:
                        next(fin, None)
                        for line in fin:
                            parts = line.rstrip("\n").split("\t")
                            if len(parts) < 2:
                                continue
                            try:
                                pos = int(parts[0])
                                score = float(parts[1])
                            except ValueError:
                                continue
                            if score < args.min_score:
                                continue
                            w.writerow([pid, pos, "ppi_interface",
                                        f"{score:.4f}", "PPI_INTERFACE"])
                            n += 1
        fh.close()
        LOG.info("%s scannet: %d sites -> %s", row.tag, n, out)


# ---------- registry + CLI ------------------------------------------------

RUNNERS: dict[str, Callable] = {
    "musitedeep": run_musitedeep,
    "scannet":    run_scannet,
}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--predictor", required=True, choices=sorted(RUNNERS.keys()))
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--min-score", type=float, default=0.5)
    ap.add_argument("--model-dir", help="MusiteDeep_web or ScanNet repo path")
    ap.add_argument("--residue-types", default="Phosphoserine_Phosphothreonine",
                    help="MusiteDeep PTM types underscore-joined; default: Phosphoserine_Phosphothreonine")
    ap.add_argument("--structure-dir",
                    help="ScanNet: dir with <tag>/*.pdb subdirs")
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
