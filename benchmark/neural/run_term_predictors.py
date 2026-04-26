#!/usr/bin/env python3
"""Sidecar for term-level FOSS function predictors that emit the standard
4-column TSV (so they auto-join the v1.1.0 ensemble)::

    protein_id<TAB>term<TAB>score<TAB>annotation_type

Predictors
----------
- ``psortb``   — bacterial subcellular localization (PSORTb 3.0, GPL-3)
- ``deepfri``  — sequence-only GO predictions (BSD-3-Clause)
- ``deepec``   — EC-number CNN (AGPL-3.0)
- ``deeparg``  — antimicrobial-resistance gene calls (MIT)

Each runner shells out to the upstream tool and rewrites its native output
into the GSPA 4-column TSV.

Usage
-----
::

    run_term_predictors.py --predictor psortb \\
        --manifest manifest.tsv --gram gramn

    run_term_predictors.py --predictor deepfri \\
        --manifest manifest.tsv --model-dir /path/to/DeepFRI

    run_term_predictors.py --predictor deepec \\
        --manifest manifest.tsv --model-dir /path/to/deepec

    run_term_predictors.py --predictor deeparg \\
        --manifest manifest.tsv --model-dir /path/to/deeparg --type prot
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
from typing import Callable, Iterator, Optional

LOG = logging.getLogger("run_term_predictors")

OUTPUT_HEADER = ["protein_id", "term", "score", "annotation_type"]


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


# ---------- PSORTb 3.0 ----------------------------------------------------

# Predicted localization label -> GO cellular_component term
PSORTB_LOC_TO_GO = {
    "cytoplasmic":          "GO:0005737",
    "cytoplasmicmembrane":  "GO:0005886",
    "periplasmic":          "GO:0042597",
    "outermembrane":        "GO:0019867",
    "extracellular":        "GO:0005576",
    "cellwall":             "GO:0005618",
    "cytoplasm":            "GO:0005737",
}

def run_psortb(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Run PSORTb 3.0 (CLI: ``psortb``).

    Native output is plain text per-protein. We parse the
    "Final Prediction" stanzas. PSORTb runs as a Docker/Singularity
    image typically; this sidecar invokes the CLI directly.
    """
    if not shutil.which("psortb"):
        raise SystemExit("psortb binary not on PATH (install brinkmanlab/psortb_commandline)")

    for row in rows:
        out = output_path(row, "psortb")
        with tempfile.TemporaryDirectory() as tmp:
            tmp_out = Path(tmp) / f"{row.tag}.psortb.txt"
            cmd = ["psortb", f"--{args.gram}", "-i", str(row.fasta_path), "-r", tmp]
            subprocess.run(cmd, check=True)
            # Find the long-format output file
            outputs = list(Path(tmp).glob("*.txt"))
            if not outputs:
                LOG.warning("no psortb output for %s", row.tag)
                continue
            tmp_out = outputs[0]

            fh, w = open_output(out)
            n = 0
            current_pid = None
            with tmp_out.open() as fin:
                for line in fin:
                    s = line.rstrip("\n")
                    if s.startswith("SeqID:"):
                        current_pid = s.split(":", 1)[1].strip().split()[0]
                    elif s.strip().startswith("Final Prediction:") and current_pid:
                        rest = s.split("Final Prediction:", 1)[1].strip()
                        # rest is like "Cytoplasmic 9.97" or "Unknown" — split last token as score
                        bits = rest.split()
                        if len(bits) >= 2:
                            label = "".join(bits[:-1]).lower()
                            try:
                                score = float(bits[-1])
                            except ValueError:
                                continue
                            # PSORTb scores are 0-10 — normalize to [0, 1]
                            norm = score / 10.0
                            if norm < args.min_score:
                                continue
                            go_term = PSORTB_LOC_TO_GO.get(label)
                            if go_term:
                                w.writerow([current_pid, go_term, f"{norm:.4f}", "GO"])
                                n += 1
            fh.close()
        LOG.info("%s psortb: %d term rows -> %s", row.tag, n, out)


# ---------- DeepFRI -------------------------------------------------------

def run_deepfri(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Run DeepFRI in sequence-only mode (no structures required).

    DeepFRI ships ``predict.py`` accepting --seqres FASTA. We run it once
    per genome and rewrite its CSV into our 4-col TSV.
    """
    if not args.model_dir:
        raise SystemExit("deepfri requires --model-dir")

    predict = Path(args.model_dir) / "predict.py"
    if not predict.exists():
        raise SystemExit(f"deepfri predict.py not found at {predict}")

    for row in rows:
        out = output_path(row, "deepfri")
        with tempfile.TemporaryDirectory() as tmp:
            cmd = ["python3", str(predict),
                   "--seqres", str(row.fasta_path),
                   "--ont", "mf,bp,cc",
                   "--output_dir", tmp,
                   "--saliency", "False"]
            subprocess.run(cmd, check=True, cwd=args.model_dir)

            fh, w = open_output(out)
            n = 0
            for csv_file in Path(tmp).glob("*_predictions.csv"):
                with csv_file.open() as fin:
                    reader = csv.DictReader(fin)
                    for r in reader:
                        try:
                            score = float(r.get("Score", 0))
                        except ValueError:
                            continue
                        if score < args.min_score:
                            continue
                        pid = r.get("Protein", "").split()[0]
                        go = r.get("GO_term", "")
                        if pid and go.startswith("GO:"):
                            w.writerow([pid, go, f"{score:.4f}", "GO"])
                            n += 1
            fh.close()
        LOG.info("%s deepfri: %d term rows -> %s", row.tag, n, out)


# ---------- DeepEC --------------------------------------------------------

def run_deepec(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Run DeepEC (Ryu et al. 2019).

    CLI: ``python3 deepec.py -i input.fasta -o output_dir``
    Produces ``DeepEC_Result.txt`` with two columns: protein_id, EC_number.
    """
    if not args.model_dir:
        raise SystemExit("deepec requires --model-dir (path to deepec repo)")
    runner = Path(args.model_dir) / "deepec.py"
    if not runner.exists():
        raise SystemExit(f"deepec.py not found at {runner}")

    for row in rows:
        out = output_path(row, "deepec")
        with tempfile.TemporaryDirectory() as tmp:
            cmd = ["python3", str(runner),
                   "-i", str(row.fasta_path),
                   "-o", tmp]
            subprocess.run(cmd, check=True, cwd=args.model_dir)

            fh, w = open_output(out)
            n = 0
            res = Path(tmp) / "DeepEC_Result.txt"
            if res.exists():
                with res.open() as fin:
                    next(fin, None)  # header
                    for line in fin:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 2:
                            continue
                        pid = parts[0]
                        ec_raw = parts[1].strip()
                        if not ec_raw or ec_raw.lower() == "none":
                            continue
                        ec = ec_raw if ec_raw.startswith("EC:") else f"EC:{ec_raw}"
                        # DeepEC doesn't emit per-call scores in some configs;
                        # default to 1.0 unless a score column is present.
                        score = 1.0
                        if len(parts) >= 3:
                            try:
                                score = float(parts[2])
                            except ValueError:
                                pass
                        if score < args.min_score:
                            continue
                        w.writerow([pid, ec, f"{score:.4f}", "EC"])
                        n += 1
            fh.close()
        LOG.info("%s deepec: %d EC rows -> %s", row.tag, n, out)


# ---------- DeepARG -------------------------------------------------------

def run_deeparg(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Run DeepARG (Arango-Argoty et al. 2018).

    CLI: ``deeparg predict --model SS --type prot --input X --output OUT``
    produces ``OUT.mapping.ARG`` with ARG class predictions.
    """
    if not shutil.which("deeparg"):
        raise SystemExit("deeparg binary not on PATH")

    for row in rows:
        out = output_path(row, "deeparg")
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp) / row.tag
            cmd = ["deeparg", "predict",
                   "--model", "SS",
                   "--type", args.deeparg_type,
                   "--input", str(row.fasta_path),
                   "--output", str(base),
                   "-d", args.model_dir or os.environ.get("DEEPARG_DB", "")]
            subprocess.run(cmd, check=True)

            fh, w = open_output(out)
            n = 0
            mapping = base.with_suffix(".mapping.ARG")
            if mapping.exists():
                with mapping.open() as fin:
                    next(fin, None)  # header
                    for line in fin:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 5:
                            continue
                        # cols: ARG_NAME, QUERY, predicted_ARG-class, best-hit, identity, evalue
                        pid = parts[1]
                        arg_class = parts[2]
                        try:
                            identity = float(parts[4])
                        except ValueError:
                            continue
                        norm = identity / 100.0
                        if norm < args.min_score:
                            continue
                        # Use an "AMR:<class>" term (no canonical IRI for AMR families)
                        w.writerow([pid, f"AMR:{arg_class}", f"{norm:.4f}", "AMR"])
                        n += 1
            fh.close()
        LOG.info("%s deeparg: %d AMR rows -> %s", row.tag, n, out)


# ---------- registry + CLI ------------------------------------------------

RUNNERS: dict[str, Callable] = {
    "psortb":  run_psortb,
    "deepfri": run_deepfri,
    "deepec":  run_deepec,
    "deeparg": run_deeparg,
}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--predictor", required=True, choices=sorted(RUNNERS.keys()))
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--min-score", type=float, default=0.1)
    ap.add_argument("--model-dir", help="DeepFRI / DeepEC / DeepARG model directory")
    ap.add_argument("--gram", default="negative",
                    help="PSORTb: positive | negative | archaea")
    ap.add_argument("--deeparg-type", default="prot",
                    help="DeepARG input type: prot | nucl")
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
