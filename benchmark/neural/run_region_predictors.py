#!/usr/bin/env python3
"""Sidecar for region-level FOSS protein predictors.

Mirrors the manifest contract of ``run_neural_predictors.py`` but emits a
5-column region TSV::

    protein_id<TAB>region_start<TAB>region_end<TAB>region_type<TAB>score

(Coordinates are 1-based inclusive.)

Predictors
----------
- ``metapredict`` — disorder regions (idptools/metapredict, MIT)
- ``deepsig``     — Sec/Tat signal peptide (BolognaBiocomp/DeepSig, GPL-3)
- ``tmbed``       — TM helix segments via ProtT5 (BernhoferM/TMbed, Apache-2)
- ``tppred3``     — N-terminal targeting peptide (BolognaBiocomp/TPpred3, GPL-3)

Each runner shells out to the upstream tool with subprocess (mirrors the
clean / proteinfer pattern in run_neural_predictors.py) and parses its
output back into the 5-column TSV. No predictor changes the manifest or
output filename conventions.

Usage
-----
::

    run_region_predictors.py --predictor metapredict \\
        --manifest manifest.tsv --min-region-len 10 --min-score 0.5

    run_region_predictors.py --predictor deepsig \\
        --manifest manifest.tsv --kingdom gramn

    run_region_predictors.py --predictor tmbed \\
        --manifest manifest.tsv --prott5-model /path/to/prott5

    run_region_predictors.py --predictor tppred3 \\
        --manifest manifest.tsv --kingdom plant
"""
from __future__ import annotations

import argparse
import csv
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterator, Optional

LOG = logging.getLogger("run_region_predictors")

OUTPUT_HEADER = ["protein_id", "region_start", "region_end", "region_type", "score"]


@dataclass
class ManifestRow:
    tag: str
    fasta_path: Path
    output_dir: Path


def read_manifest(path: Path) -> list[ManifestRow]:
    rows: list[ManifestRow] = []
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        required = {"tag", "fasta_path", "output_dir"}
        missing = required - set(header)
        if missing:
            raise SystemExit(f"manifest missing columns: {sorted(missing)}")
        idx = {k: header.index(k) for k in required}
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


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    name: Optional[str] = None
    chunks: list[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            yield name, "".join(chunks)


def output_path(row: ManifestRow, predictor: str) -> Path:
    row.output_dir.mkdir(parents=True, exist_ok=True)
    return row.output_dir / f"{row.tag}.{predictor}.tsv"


def open_output(path: Path):
    fh = path.open("w", newline="")
    writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
    writer.writerow(OUTPUT_HEADER)
    return fh, writer


# ---------- Metapredict ---------------------------------------------------

def run_metapredict(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Use metapredict's Python API to compute IDR boundaries."""
    import metapredict as meta

    for row in rows:
        out = output_path(row, "metapredict")
        fh, w = open_output(out)
        n_regions = 0
        for pid, seq in iter_fasta(row.fasta_path):
            if not seq:
                continue
            try:
                # metapredict.predict_disorder_domains returns object with .disordered_domain_boundaries
                result = meta.predict_disorder_domains(seq)
                boundaries = result.disordered_domain_boundaries  # list of (start, end) 0-based half-open
                scores = result.disorder  # per-residue scores
            except Exception as exc:
                LOG.warning("metapredict failed on %s: %s", pid, exc)
                continue
            for s, e in boundaries:
                if (e - s) < args.min_region_len:
                    continue
                # Convert to 1-based inclusive
                start_1 = s + 1
                end_1 = e
                mean_score = float(sum(scores[s:e]) / max(1, e - s))
                if mean_score < args.min_score:
                    continue
                w.writerow([pid, start_1, end_1, "disorder", f"{mean_score:.4f}"])
                n_regions += 1
        fh.close()
        LOG.info("%s metapredict: %d regions -> %s", row.tag, n_regions, out)


# ---------- DeepSig -------------------------------------------------------

def run_deepsig(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Run DeepSig (BolognaBiocomp/DeepSig).

    Real CLI (deepsig-biocomp 0.9):
        deepsig -f INPUT -o OUTPUT -k {euk,gramp,gramn} [-m {json,gff3}]

    Output is GFF3 by default. Each predicted feature is one line:
        seqid  source  type  start  end  score  strand  phase  attributes

    The ``type`` column is human-readable with spaces ("Signal peptide",
    "Lipoprotein signal peptide", ...). We lowercase + underscore so the
    GSPA report's REGION_TYPE_TO_CLASS map can resolve it.

    DeepSig requires the ``DEEPSIG_ROOT`` env var; auto-set from the pip
    install location when missing.
    """
    if not shutil.which("deepsig"):
        raise SystemExit("deepsig binary not on PATH")

    env = os.environ.copy()
    if "DEEPSIG_ROOT" not in env:
        try:
            import deepsig as _ds
            env["DEEPSIG_ROOT"] = str(Path(_ds.__file__).parent)
            LOG.info("auto-set DEEPSIG_ROOT=%s", env["DEEPSIG_ROOT"])
        except ImportError:
            raise SystemExit("DEEPSIG_ROOT not set and deepsig package not importable")

    for row in rows:
        out = output_path(row, "deepsig")
        with tempfile.TemporaryDirectory() as tmp:
            tmp_out = Path(tmp) / f"{row.tag}.deepsig.gff3"
            cmd = ["deepsig", "-f", str(row.fasta_path),
                   "-k", args.kingdom, "-o", str(tmp_out),
                   "-m", "gff3"]
            subprocess.run(cmd, check=True, env=env)

            fh, w = open_output(out)
            n_regions = 0
            with tmp_out.open() as fin:
                for line in fin:
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 6:
                        continue
                    pid = parts[0]
                    feat_type = parts[2].lower().replace(" ", "_")
                    try:
                        start = int(parts[3])
                        end = int(parts[4])
                        score = float(parts[5])
                    except ValueError:
                        continue
                    if score < args.min_score:
                        continue
                    # Normalise common DeepSig type strings to GSPA region types
                    if "signal_peptide" in feat_type or feat_type.startswith("signal_"):
                        region_type = "signal_peptide"
                    elif "lipoprotein" in feat_type:
                        region_type = "lipo"
                    elif "tat" in feat_type or "twin-arginine" in feat_type:
                        region_type = "tat_signal_peptide"
                    elif feat_type == "transit_peptide":
                        region_type = "transit_peptide"
                    else:
                        # Unknown DeepSig feature type — keep verbatim (lower+underscore)
                        region_type = feat_type
                    w.writerow([pid, start, end, region_type, f"{score:.4f}"])
                    n_regions += 1
            fh.close()
        LOG.info("%s deepsig: %d regions -> %s", row.tag, n_regions, out)


# ---------- TMbed ---------------------------------------------------------

def run_tmbed(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Run TMbed via its ``tmbed`` CLI; parse 3-line topology output.

    TMbed CLI: ``tmbed predict -f input.fasta -o out.txt``
    Output is FASTA-like 3-line per protein: header, sequence, topology
    where topology chars are: H (TM helix), B (TM β-strand), S (signal
    peptide), I (inside), O (outside).
    """
    if not shutil.which("tmbed"):
        raise SystemExit("tmbed binary not on PATH")

    for row in rows:
        out = output_path(row, "tmbed")
        with tempfile.TemporaryDirectory() as tmp:
            tmp_out = Path(tmp) / f"{row.tag}.tmbed.3line"
            cmd = ["tmbed", "predict", "-f", str(row.fasta_path),
                   "-p", str(tmp_out), "--out-format", "0"]  # 3-line plain
            subprocess.run(cmd, check=True)

            fh, w = open_output(out)
            n_regions = 0
            with tmp_out.open() as fin:
                lines = [ln.rstrip("\n") for ln in fin if ln.strip()]
            i = 0
            while i + 2 < len(lines):
                header = lines[i]
                # seq    = lines[i+1]
                topo   = lines[i+2]
                i += 3
                if not header.startswith(">"):
                    continue
                pid = header[1:].split()[0]
                # Extract runs of H, B, S
                for code, label in (("H", "tm_helix"),
                                     ("B", "tm_beta"),
                                     ("S", "signal_peptide")):
                    for m in re.finditer(f"{code}+", topo):
                        s, e = m.start() + 1, m.end()  # 1-based inclusive
                        if (e - s + 1) < args.min_region_len and label != "signal_peptide":
                            continue
                        w.writerow([pid, s, e, label, "1.0000"])
                        n_regions += 1
            fh.close()
        LOG.info("%s tmbed: %d regions -> %s", row.tag, n_regions, out)


# ---------- TPpred 3 ------------------------------------------------------

def run_tppred3(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Run TPpred3 (BolognaBiocomp).

    CLI: ``tppred3.py -f input.fasta -k {plant,nonplant} -o out.tsv``
    Output has cleavage-site predictions for transit peptides.
    """
    if not shutil.which("tppred3.py"):
        raise SystemExit("tppred3.py not on PATH")

    for row in rows:
        out = output_path(row, "tppred3")
        with tempfile.TemporaryDirectory() as tmp:
            tmp_out = Path(tmp) / f"{row.tag}.tppred3.tsv"
            subprocess.run(
                ["tppred3.py", "-f", str(row.fasta_path),
                 "-k", args.kingdom, "-o", str(tmp_out)],
                check=True,
            )
            fh, w = open_output(out)
            n_regions = 0
            with tmp_out.open() as fin:
                for line in fin:
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 5:
                        continue
                    pid = parts[0]
                    # TPpred3 columns: id, prediction, cleavage_site, score, ...
                    pred = parts[1].lower()
                    try:
                        cs = int(parts[2])
                        score = float(parts[3])
                    except ValueError:
                        continue
                    if score < args.min_score:
                        continue
                    if pred in ("none", "no", "n", "0"):
                        continue
                    region_type = (
                        "mito_targeting"   if "mito"   in pred else
                        "chloro_targeting" if "chloro" in pred else
                        "targeting_peptide"
                    )
                    w.writerow([pid, 1, cs, region_type, f"{score:.4f}"])
                    n_regions += 1
            fh.close()
        LOG.info("%s tppred3: %d regions -> %s", row.tag, n_regions, out)


# ---------- registry + CLI ------------------------------------------------

RUNNERS: dict[str, Callable[[list[ManifestRow], argparse.Namespace], None]] = {
    "metapredict": run_metapredict,
    "deepsig":     run_deepsig,
    "tmbed":       run_tmbed,
    "tppred3":     run_tppred3,
}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--predictor", required=True, choices=sorted(RUNNERS.keys()))
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--min-score", type=float, default=0.5,
                    help="Drop predictions below this score (default 0.5)")
    ap.add_argument("--min-region-len", type=int, default=10,
                    help="Drop regions shorter than this many residues (default 10)")
    ap.add_argument("--kingdom", default="gramn",
                    help="DeepSig: gramp|gramn|euk; TPpred3: plant|nonplant")
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
