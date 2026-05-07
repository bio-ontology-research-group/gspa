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
    """Run PSORTb 3.0 via the brinkmanlab/psortb_commandline Singularity
    image (GPL-3.0).

    The image's perl wrapper writes long-format output to
    ``/results/<timestamp>_psortb_<gram>.txt`` inside the container, so
    we bind a host scratch dir to ``/results`` and parse that file. The
    input FASTA is bound as ``/in/<basename>``.

    Args:
        --psortb-sif <path>   Singularity image (mandatory for psortb)
        --gram {positive,negative,archaea}
    """
    if not args.psortb_sif:
        raise SystemExit("psortb requires --psortb-sif <SIF path>")
    if not Path(args.psortb_sif).exists():
        raise SystemExit(f"psortb SIF not found: {args.psortb_sif}")
    if not shutil.which("singularity"):
        raise SystemExit("singularity binary not on PATH")

    gram_flag = f"--{args.gram}"  # --positive | --negative | --archaea

    for row in rows:
        out = output_path(row, "psortb")
        # The brinkmanlab perl wrapper hardcodes /tmp/results inside the
        # container. Singularity 3.9 will not auto-create bind targets, so
        # we bind a host dir to /tmp and pre-create /tmp/results on it.
        bind_root = Path(tempfile.mkdtemp(prefix=f"psortb-{row.tag}-"))
        results_dir = bind_root / "results"
        results_dir.mkdir(parents=True, exist_ok=True)
        try:
            in_path = row.fasta_path.resolve()
            in_dir = in_path.parent
            cmd = ["singularity", "exec",
                   "--bind", f"{bind_root}:/tmp",
                   "--bind", f"{in_dir}:/in",
                   args.psortb_sif,
                   "/usr/local/psortb/bin/psort",
                   gram_flag,
                   "-i", f"/in/{in_path.name}",
                   "-o", "long"]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                LOG.warning("%s: psort exit %d; stderr=%s", row.tag, r.returncode, r.stderr[:300])
                continue

            result_files = sorted(results_dir.glob("*.txt"))
            if not result_files:
                LOG.warning("%s: no PSORTb output in %s", row.tag, results_dir)
                continue

            fh, w = open_output(out)
            n = 0
            # PSORTb 3.0 ``-o long`` actually emits a single TSV with one row
            # per protein and per-sub-tool columns plus three summary columns:
            # ``Final_Localization``, ``Final_Localization_Details``, ``Final_Score``.
            with result_files[0].open() as fin:
                header_line = fin.readline().rstrip("\n").split("\t")
                try:
                    i_id    = header_line.index("SeqID")
                    i_loc   = header_line.index("Final_Localization")
                    i_score = header_line.index("Final_Score")
                except ValueError:
                    LOG.warning("%s: PSORTb header missing expected columns", row.tag)
                    continue
                for line in fin:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) <= max(i_id, i_loc, i_score):
                        continue
                    pid = fields[i_id].split()[0]   # take first whitespace token
                    label_raw = fields[i_loc].strip()
                    if not label_raw or label_raw.lower() == "unknown":
                        continue
                    try:
                        score = float(fields[i_score])
                    except ValueError:
                        continue
                    norm = score / 10.0
                    if norm < args.min_score:
                        continue
                    label = label_raw.replace(" ", "").lower()
                    go_term = PSORTB_LOC_TO_GO.get(label)
                    if go_term:
                        w.writerow([pid, go_term, f"{norm:.4f}", "GO"])
                        n += 1
            fh.close()
        finally:
            shutil.rmtree(bind_root, ignore_errors=True)
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
            # Upstream DeepFRI flags: --fasta_fn (not --seqres), -ont (single
            # dash) accepts multiple choices, -o is a file prefix.
            out_prefix = str(Path(tmp) / row.tag)
            cmd = ["python3", str(predict),
                   "--fasta_fn", str(row.fasta_path),
                   "-ont", "mf", "bp", "cc",
                   "-o", out_prefix]
            subprocess.run(cmd, check=True, cwd=args.model_dir)

            fh, w = open_output(out)
            n = 0
            # DeepFRI CSV: line 1 is "### Predictions made by DeepFRI." comment
            # line 2 is header: Protein,GO_term/EC_number,Score,...
            for csv_file in Path(tmp).glob("*_predictions.csv"):
                with csv_file.open() as fin:
                    first = fin.readline()
                    if not first.startswith("###"):
                        fin.seek(0)
                    reader = csv.DictReader(fin)
                    for r in reader:
                        try:
                            score = float(r.get("Score", 0))
                        except (ValueError, TypeError):
                            continue
                        if score < args.min_score:
                            continue
                        pid = (r.get("Protein") or "").split()[0]
                        go = r.get("GO_term/EC_number") or r.get("GO_term") or ""
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
            # DeepEC pickles use 'sklearn.preprocessing.label' (private path
            # before sklearn 0.22). Modern sklearn moved it to '_label'.
            # Inject a shim + add deepec/ to sys.path so its 'from deepec
            # import ...' works when invoked via exec().
            wrapper = Path(tmp) / "_deepec_wrapper.py"
            wrapper.write_text(
                "import sys, sklearn.preprocessing as _sp\n"
                "sys.modules['sklearn.preprocessing.label'] = _sp._label\n"
                f"sys.path.insert(0, {str(args.model_dir)!r})\n"
                f"sys.argv = [{str(runner)!r}, '-i', {str(row.fasta_path)!r}, '-o', {tmp!r}]\n"
                f"exec(open({str(runner)!r}).read())\n"
            )
            cmd = ["python3", str(wrapper)]
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
    """Run DeepARG (Arango-Argoty et al. 2018, MIT).

    CLI (pip or Singularity-bound):
        deeparg predict --model SS --type prot --input X --output OUT -d DB

    Output ``<OUT>.mapping.ARG`` is a TSV with ARG class predictions.

    With ``--deeparg-sif`` set, run via Singularity (binds the input dir
    and the model DB into the container). Else delegate to a pip install.
    """
    use_sif = bool(args.deeparg_sif)
    if not use_sif and not shutil.which("deeparg"):
        raise SystemExit("deeparg not on PATH (and no --deeparg-sif provided)")

    db_dir = args.model_dir or os.environ.get("DEEPARG_DB", "")
    if not db_dir:
        raise SystemExit("deeparg requires --model-dir <database dir>")

    for row in rows:
        out = output_path(row, "deeparg")
        with tempfile.TemporaryDirectory() as tmp:
            base_name = row.tag
            if use_sif:
                in_dir = row.fasta_path.parent.resolve()
                cmd = ["singularity", "exec",
                       "--bind", f"{in_dir}:/in",
                       "--bind", f"{tmp}:/out",
                       "--bind", f"{db_dir}:/db",
                       args.deeparg_sif,
                       "deeparg", "predict",
                       "--model", "SS",
                       "--type", args.deeparg_type,
                       "--input", f"/in/{row.fasta_path.name}",
                       "--output", f"/out/{base_name}",
                       "-d", "/db"]
                base = Path(tmp) / base_name
            else:
                base = Path(tmp) / base_name
                cmd = ["deeparg", "predict",
                       "--model", "SS",
                       "--type", args.deeparg_type,
                       "--input", str(row.fasta_path),
                       "--output", str(base),
                       "-d", db_dir]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                LOG.warning("%s: deeparg exit %d; stderr=%s",
                            row.tag, r.returncode, r.stderr[:300])
                continue

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


# ---------- mDeepFRI ------------------------------------------------------

_MDF_ASPECT = {
    "GO Biological Process": "GO",
    "GO Molecular Function": "GO",
    "GO Cellular Component": "GO",
    "Enzyme Commission": "EC",
}


def run_mdf(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Run metagenomic-deepFRI (Bezshapkin et al. 2026, BSD-3-Clause).

    CLI: ``mDeepFRI predict-function -i FASTA -o OUTDIR -w WEIGHTS [-d FOLDCOMP_DB] [--skip-pdb]``

    Native output is ``OUTDIR/results.tsv`` with columns
    ``protein, network_type, prediction_mode, go_term, score, ...``.
    Each row is a single (protein, term, aspect) prediction. We rewrite
    GO and EC rows into the canonical 4-column TSV; the prediction_mode
    column tells us BP/MF/CC vs. EC. Scores below ``--min-score`` are
    dropped.

    Required: ``--mdf-weights`` (path to mDeepFRI weights dir from
    ``mDeepFRI get-models --version 1.0``). Use ``--mdf-skip-pdb`` for
    sequence-only mode (no FoldComp DB needed). For the structure path,
    pass ``--mdf-foldcomp-db /path/to/db`` (FoldComp-format).
    """
    if not args.mdf_weights:
        raise SystemExit("mdf requires --mdf-weights (mDeepFRI weights dir)")

    mdf_bin = args.mdf_executable or "mDeepFRI"

    for row in rows:
        out = output_path(row, "mdf")
        with tempfile.TemporaryDirectory() as tmp:
            cmd = [mdf_bin, "predict-function",
                   "-i", str(row.fasta_path),
                   "-o", tmp,
                   "-w", args.mdf_weights,
                   "-t", str(args.mdf_threads)]
            if args.mdf_skip_pdb:
                cmd.append("--skip-pdb")
            if args.mdf_foldcomp_db:
                cmd.extend(["-d", args.mdf_foldcomp_db])
            subprocess.run(cmd, check=True)

            results_tsv = Path(tmp) / "results.tsv"
            if not results_tsv.exists():
                LOG.error("%s mdf: no results.tsv at %s", row.tag, results_tsv)
                continue

            fh, w = open_output(out)
            n = 0
            seen: dict[tuple[str, str], float] = {}
            with results_tsv.open() as fin:
                reader = csv.DictReader(fin, delimiter="\t")
                for r in reader:
                    pid = (r.get("protein") or "").split()[0]
                    term = r.get("go_term") or ""
                    mode = r.get("prediction_mode") or ""
                    ann_type = _MDF_ASPECT.get(mode)
                    if not (pid and term and ann_type):
                        continue
                    try:
                        score = float(r.get("score", 0))
                    except (ValueError, TypeError):
                        continue
                    if score < args.min_score:
                        continue
                    # Dedup (protein, term) keeping the highest score.
                    key = (pid, term)
                    prev = seen.get(key)
                    if prev is None or score > prev:
                        seen[key] = score
            for (pid, term), score in seen.items():
                # Annotation type follows the term namespace: GO:* → GO,
                # EC:* → EC. Recovers the type from the term itself.
                ann_type = "EC" if term.startswith("EC:") else "GO"
                w.writerow([pid, term, f"{score:.4f}", ann_type])
                n += 1
            fh.close()
        LOG.info("%s mdf: %d term rows -> %s", row.tag, n, out)


# ---------- registry + CLI ------------------------------------------------

RUNNERS: dict[str, Callable] = {
    "psortb":  run_psortb,
    "deepfri": run_deepfri,
    "deepec":  run_deepec,
    "deeparg": run_deeparg,
    "mdf":     run_mdf,
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
    ap.add_argument("--psortb-sif",
                    help="Path to brinkmanlab/psortb_commandline Singularity SIF")
    ap.add_argument("--deeparg-type", default="prot",
                    help="DeepARG input type: prot | nucl")
    ap.add_argument("--deeparg-sif",
                    help="Path to deeparg Singularity SIF (optional)")
    ap.add_argument("--mdf-weights",
                    help="mDeepFRI: directory with model weights (from "
                         "`mDeepFRI get-models --version 1.0 --output ...`)")
    ap.add_argument("--mdf-executable", default=None,
                    help="mDeepFRI: path to the mDeepFRI binary (default: "
                         "first on PATH; use this when running from a venv)")
    ap.add_argument("--mdf-skip-pdb", action="store_true",
                    help="mDeepFRI: skip PDB100 search (sequence-only path)")
    ap.add_argument("--mdf-foldcomp-db", default=None,
                    help="mDeepFRI: FoldComp structures database path "
                         "(structure-supplied path; AFDB / ESM-Atlas)")
    ap.add_argument("--mdf-threads", type=int, default=4,
                    help="mDeepFRI: thread count")
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
