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
    """Run MusiteDeep_web's predict_multi_batch.py.

    The upstream image ``duolinwang/musitedeep_backend:2.0`` ships
    ``predict_multi_batch.py`` at ``/musite/MusiteDeep/`` with 13 PTM
    model directories under ``/musite/MusiteDeep/models/<PTM_type>/``.
    The CLI takes:
      - ``-input``  protein FASTA
      - ``-output`` PREFIX (not a directory) for ``<prefix>_results.txt``
      - ``-model-prefix`` path to one model dir, or ``;``-separated paths
        for multi-PTM joint prediction.

    There is no ``-residue-types`` flag (that was a training-script flag);
    --residue-types here selects which model directories to point at.

    Output format (single file ``<prefix>_results.txt``)::

        ID<TAB>Position<TAB>Residue<TAB>PTMscores<TAB>Cutoff=0.5
        >WP_xxx ...                          # interleaved FASTA-style headers
        WP_xxx<TAB>13<TAB>S<TAB>Phosphoserine:0.639<TAB>Phosphoserine:0.639
        WP_xxx<TAB>17<TAB>T<TAB>Phosphothreonine:0.099<TAB>None

    Score column is ``<PTM_type>:<probability>`` per residue.
    """
    if not args.musitedeep_sif and not args.model_dir:
        raise SystemExit(
            "musitedeep requires --musitedeep-sif "
            "(or --model-dir for native install)")

    # Map --residue-types to model-prefix paths
    rt_to_dir = {
        "Phosphoserine_Phosphothreonine": "Phosphoserine_Phosphothreonine",
        "Phosphotyrosine":                "Phosphotyrosine",
        "N6-acetyllysine":                "N6-acetyllysine",
        "Methyllysine":                   "Methyllysine",
        "Methylarginine":                 "Methylarginine",
        "Ubiquitination":                 "Ubiquitination",
        "SUMOylation":                    "SUMOylation",
        "Hydroxylysine":                  "Hydroxylysine",
        "Hydroxyproline":                 "Hydroxyproline",
        "N-linked_glycosylation":         "N-linked_glycosylation",
        "O-linked_glycosylation":         "O-linked_glycosylation",
        "Pyrrolidone_carboxylic_acid":    "Pyrrolidone_carboxylic_acid",
        "S-palmitoyl_cysteine":           "S-palmitoyl_cysteine",
    }
    if args.musitedeep_sif:
        models_root = "/musite/MusiteDeep/models"
        runner_in_container = "/musite/MusiteDeep/predict_multi_batch.py"
    else:
        models_root = str(Path(args.model_dir) / "MusiteDeep" / "models")
        runner_in_container = None  # not used in native mode

    requested = [rt.strip() for rt in args.residue_types.split(",") if rt.strip()]
    bad = [rt for rt in requested if rt not in rt_to_dir]
    if bad:
        raise SystemExit(f"unknown MusiteDeep PTM types: {bad}; "
                         f"available: {sorted(rt_to_dir)}")
    model_prefix = ";".join(f"{models_root}/{rt_to_dir[rt]}" for rt in requested)

    for row in rows:
        out = output_path(row, "musitedeep")
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            out_prefix = f"{tmp_path}/{row.tag}"
            if args.musitedeep_sif:
                in_dir = row.fasta_path.parent.resolve()
                # Inside container: bind in-dir → /in, tmp → /out
                container_prefix = f"/out/{row.tag}"
                cmd = ["singularity", "exec",
                       "--bind", f"{in_dir}:/in",
                       "--bind", f"{tmp_path}:/out",
                       args.musitedeep_sif,
                       "sh", "-c",
                       f"cd /musite/MusiteDeep && python3 {runner_in_container} "
                       f"-input /in/{row.fasta_path.name} "
                       f"-output {container_prefix} "
                       f"-model-prefix '{model_prefix}'"]
            else:
                runner = Path(args.model_dir) / "MusiteDeep" / "predict_multi_batch.py"
                cmd = ["python3", str(runner),
                       "-input", str(row.fasta_path),
                       "-output", out_prefix,
                       "-model-prefix", model_prefix]
            subprocess.run(cmd, check=True)

            fh, w = open_output(out)
            n = 0
            results_file = tmp_path / f"{row.tag}_results.txt"
            if not results_file.exists():
                LOG.warning("%s: no MusiteDeep results file at %s",
                            row.tag, results_file)
                fh.close()
                continue
            with results_file.open() as fin:
                next(fin, None)  # header line
                for line in fin:
                    if line.startswith(">") or not line.strip():
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 5:
                        continue
                    pid = parts[0]
                    try:
                        pos = int(parts[1])
                    except ValueError:
                        continue
                    cutoff_field = parts[4]
                    if cutoff_field == "None":
                        # below 0.5 — skip unless user wants raw scores
                        if args.min_score >= 0.5:
                            continue
                    # PTMscores column: '<type>:<prob>'
                    score_field = parts[3]
                    site_type = "phosphosite"
                    score = 0.0
                    if ":" in score_field:
                        ptm_type, _, prob_str = score_field.partition(":")
                        try:
                            score = float(prob_str)
                        except ValueError:
                            continue
                        site_type = ptm_type.lower()
                    if score < args.min_score:
                        continue
                    w.writerow([pid, pos, site_type,
                                f"{score:.4f}", "PTM_SITE"])
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
                    help="MusiteDeep PTM types comma-joined; default: "
                         "Phosphoserine_Phosphothreonine. Multiple e.g. "
                         "Phosphoserine_Phosphothreonine,Phosphotyrosine")
    ap.add_argument("--musitedeep-sif",
                    help="Path to duolinwang/musitedeep_backend Singularity image")
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
