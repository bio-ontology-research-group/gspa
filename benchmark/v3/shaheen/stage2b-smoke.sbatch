#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-stage2b-smoke
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/price-stage2b-rc1/slurm/stage2b-smoke-%j.log
# Stage 2B (flux-certified refiner) smoke test on the 22 Price genomes' raw claims on Shaheen III.
# When invoked outside Slurm (SLURM_JOB_ID unset and RUN_PAYLOAD!=1):
#   - default (no arguments): prints the exact sbatch submission command without submitting;
#   - --submit: submits this script via sbatch.
# Inside Slurm (or with RUN_PAYLOAD=1):
#   1. Builds --fba-data-dir from the pinned /project/k10018/sw/gapsmith-eq-ibex/data_merged.tar
#      using benchmark/v3/stage2b/build_fba_data.py --verify-pinned-sha256.
#   2. Streams each of the 22 Price genomes' raw claims.jsonl under
#      /project/k10018/gspa-v3/price-full-matrix-v2/parts/<predictor>/<assembly>/claims.jsonl
#      into the sidecar payload format (gene_ec_scores[proteinId][clean_ec] = max score).
#   3. Runs scripts/stoichiometric_locus_refiner.py (--gram neg --medium minimal --export-model)
#      across all 22 assemblies and records smoke-summary.json + receipt.sha256.
set -euo pipefail

OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/price-stage2b-rc1}
PARTS_ROOT=${PARTS_ROOT:-/project/k10018/gspa-v3/price-full-matrix-v2/parts}
GAPSMITH_TAR=${GAPSMITH_TAR:-/project/k10018/sw/gapsmith-eq-ibex/data_merged.tar}
REPO_ROOT=${REPO_ROOT:-/project/k10018/code/gspa-v3}
FBA_PYTHON=${FBA_PYTHON:-/project/k10018/sw/envs/gspa-v3/bin/python}
PREDICTORS=${PREDICTORS:-clean}
GRAM=${GRAM:-neg}
MEDIUM=${MEDIUM:-minimal}
MAX_PARALLEL=${MAX_PARALLEL:-8}

if [ -z "${SLURM_JOB_ID:-}" ] && [ "${RUN_PAYLOAD:-0}" != 1 ]; then
  submit=0
  case "${1:-}" in
    --submit) submit=1 ;;
    "") ;;
    *) echo "usage: $0 [--submit]" >&2; exit 2 ;;
  esac
  script_path=$(cd "$(dirname "$0")" && pwd)/$(basename "$0")
  cmd=(
    sbatch --parsable
    --account=k10018
    --partition=workq
    --nodes=1
    --time=02:00:00
    --no-requeue
    --exclude=nid00022,nid00497,nid01033
    --job-name=gspa-stage2b-smoke
    --output="$OUT_ROOT/slurm/stage2b-smoke-%j.log"
    "$script_path"
  )
  if [ "$submit" = 1 ]; then
    command -v sbatch >/dev/null || { echo "sbatch not found: run on a Shaheen login node" >&2; exit 2; }
    mkdir -p "$OUT_ROOT/slurm"
    jid=$("${cmd[@]}")
    echo "submitted stage2b-smoke job=$jid"
  else
    echo "# Stage 2B smoke test on 22 Price genomes (prints by default; pass --submit to submit)"
    printf 'mkdir -p %q\n' "$OUT_ROOT/slurm"
    printf '%q ' "${cmd[@]}"
    printf '\n'
  fi
  exit 0
fi

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
JOB=${SLURM_JOB_ID:-local$$}
SMOKE_OUT="$OUT_ROOT/smoke/$JOB"
FBA_DATA_DIR="$SMOKE_OUT/fba_data"
mkdir -p "$SMOKE_OUT" "$FBA_DATA_DIR"
date --iso-8601=seconds > "$SMOKE_OUT/started-at.txt"
hostname > "$SMOKE_OUT/hostname.txt"
cp "$0" "$SMOKE_OUT/executed.sbatch"

# Step 1: Build --fba-data-dir from pinned data_merged.tar
/usr/bin/time -v -o "$SMOKE_OUT/build_fba_data.resources.txt" \
  "$FBA_PYTHON" "$REPO_ROOT/benchmark/v3/stage2b/build_fba_data.py" \
  --gapsmith-source "$GAPSMITH_TAR" \
  --output-dir "$FBA_DATA_DIR" \
  --verify-pinned-sha256 > "$SMOKE_OUT/build_fba_data.log" 2>&1

# Step 2 & 3: Stream raw claims.jsonl for each of the 22 Price genomes and run the sidecar
export SMOKE_OUT FBA_DATA_DIR PARTS_ROOT REPO_ROOT FBA_PYTHON PREDICTORS GRAM MEDIUM MAX_PARALLEL
"$FBA_PYTHON" - <<'PY'
import concurrent.futures
import hashlib
import json
import os
import subprocess
import time
from pathlib import Path

smoke_out = Path(os.environ["SMOKE_OUT"])
fba_data_dir = Path(os.environ["FBA_DATA_DIR"])
parts_root = Path(os.environ["PARTS_ROOT"])
repo_root = Path(os.environ["REPO_ROOT"])
fba_python = os.environ["FBA_PYTHON"]
predictors = [p.strip() for p in os.environ["PREDICTORS"].replace(",", " ").split() if p.strip()]
gram = os.environ["GRAM"]
medium = os.environ["MEDIUM"]
max_parallel = int(os.environ.get("MAX_PARALLEL", "8"))
sidecar_script = repo_root / "scripts" / "stoichiometric_locus_refiner.py"


def extract_gene_ec_scores(claims_path: Path) -> tuple[dict[str, dict[str, float]], int]:
    scores: dict[str, dict[str, float]] = {}
    n_ec_claims = 0
    with open(claims_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if '"EC' not in line:
                continue
            rec = json.loads(line)
            # Raw predictor claims (price-full-matrix-v2 parts): function_type, function_id, protein_id, raw_score.
            if rec.get("function_type") != "EC":
                continue
            clean_ec = str(rec.get("function_id", "")).strip()
            if clean_ec.startswith("EC:"):
                clean_ec = clean_ec[3:]
            parts = clean_ec.split(".")
            if len(parts) != 4 or not all(p.isdigit() for p in parts):
                continue
            pid = str(rec.get("protein_id", "")).strip()
            if not pid:
                continue
            prob = float(rec.get("raw_score", 0.0))
            if prob <= 0.0:
                continue
            n_ec_claims += 1
            gmap = scores.setdefault(pid, {})
            if prob > gmap.get(clean_ec, 0.0):
                gmap[clean_ec] = prob
    return scores, n_ec_claims


def run_one(task: tuple[str, str, Path]) -> dict:
    predictor, assembly, claims_path = task
    cell_dir = smoke_out / "cells" / predictor / assembly
    cell_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    gene_ec_scores, n_ec_claims = extract_gene_ec_scores(claims_path)
    in_json = cell_dir / "fba_input.json"
    out_json = cell_dir / "fba_output.json"
    model_xml = cell_dir / "metabolic_model.xml"
    res_txt = cell_dir / "resources.txt"
    log_txt = cell_dir / "sidecar.log"

    in_json.write_text(
        json.dumps(
            {
                "gram": gram,
                "medium": medium,
                "tau_core": 0.50,
                "alpha_bio": 0.10,
                "gene_ec_scores": gene_ec_scores,
            }
        ),
        encoding="utf-8",
    )
    cmd = [
        "/usr/bin/time",
        "-v",
        "-o",
        str(res_txt),
        fba_python,
        str(sidecar_script),
        "--fba-data-dir",
        str(fba_data_dir),
        "--input-json",
        str(in_json),
        "--output-json",
        str(out_json),
        "--export-model",
        str(model_xml),
        "--gram",
        gram,
        "--medium",
        medium,
    ]
    with open(log_txt, "w", encoding="utf-8") as lf:
        subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, check=True)
    elapsed = round(time.time() - t0, 3)
    out = json.loads(out_json.read_text(encoding="utf-8"))
    return {
        "predictor": predictor,
        "assembly": assembly,
        "claims_path": str(claims_path),
        "n_ec_claims": n_ec_claims,
        "n_genes_with_ec": len(gene_ec_scores),
        "elapsed_sec": elapsed,
        "core_bio": out["core_bio"],
        "max_univ_bio": out["max_univ_bio"],
        "target_bio": out["target_bio"],
        "n_core_rxns": out["n_core_rxns"],
        "n_repaired_rxns": out["n_repaired_rxns"],
        "n_flux_certified_repaired": out["n_flux_certified_repaired"],
        "n_updated_genes": len(out.get("updated_gene_ec_scores", {})),
        "metabolic_model_xml_sha256": hashlib.sha256(model_xml.read_bytes()).hexdigest(),
        "reactions_tsv_sha256": hashlib.sha256((cell_dir / "reactions.tsv").read_bytes()).hexdigest(),
    }


tasks: list[tuple[str, str, Path]] = []
for pred in predictors:
    pred_dir = parts_root / pred
    if not pred_dir.is_dir():
        raise FileNotFoundError(f"Predictor parts directory not found: {pred_dir}")
    for asm_dir in sorted(pred_dir.iterdir()):
        claims_file = asm_dir / "claims.jsonl"
        if asm_dir.is_dir() and claims_file.is_file():
            tasks.append((pred, asm_dir.name, claims_file))

if len(tasks) != 22 * len(predictors):
    raise RuntimeError(f"Expected {22 * len(predictors)} Price genome claim files, found {len(tasks)}")

results = []
with concurrent.futures.ThreadPoolExecutor(max_workers=max_parallel) as pool:
    for rec in pool.map(run_one, tasks):
        results.append(rec)

summary = {
    "schema": "gspa.stage2b-smoke/1",
    "job": os.environ.get("SLURM_JOB_ID", "local"),
    "fba_data_receipt_sha256": hashlib.sha256((fba_data_dir / "receipt.json").read_bytes()).hexdigest(),
    "predictors": predictors,
    "n_cells": len(results),
    "gram": gram,
    "medium": medium,
    "cells": results,
}
summary_path = smoke_out / "smoke-summary.json"
summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY

sha256sum "$SMOKE_OUT/smoke-summary.json" "$FBA_DATA_DIR/receipt.json" > "$SMOKE_OUT/receipt.sha256"
date --iso-8601=seconds > "$SMOKE_OUT/finished-at.txt"
