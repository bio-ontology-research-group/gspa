#!/bin/bash
# P08 -- Price-149 (RB-TnSeq) EC benchmark: re-extraction of the six predictors' raw outputs
# on ibex + local evaluation (Table 6, Fig. 6b).
#
# Usage
#   On the ibex login node (submits jobs only; never compute there):
#     bash run.sh submit            # 6-task debug array (one predictor each) + dependent merge job
#   Inside Slurm (what `submit` sbatches):
#     run.sh extract <predictor>    # writes $OUT/parts/<predictor>.pkl.gz and sha256 TSV
#     run.sh merge                  # merges the parts -> $OUT/price149_all6_compact.pkl.gz
#   Locally (repo root = ../../..):
#     bash run.sh fetch             # rsync the merged compact pickle + sha256 TSVs from ibex
#     bash run.sh eval              # python scripts/eval_price149_real_all6.py on the fetched pickle
#
# Copy this folder to $PIPE (below) on ibex before `submit`.
#SBATCH -p debug
#SBATCH -A pi-hohndor
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 01:00:00

set -euo pipefail
ROOT=/ibex/scratch/projects/c2014/hohndor/price149_repro_2026-10-19   # only write location on ibex
PIPE=$ROOT/P08_price149
OUT=$ROOT/out
MANIFEST=/ibex/scratch/projects/c2014/hohndor/suite6_price149_manifest.tsv
PY=/home/hohndor/.conda/envs/mean/bin/python    # Python 3.8.13, pandas 1.5.3, numpy 1.23.5
PREDICTORS="clean dpz enzbert graphec mapred topec"

mode=${1:-help}
case "$mode" in
  submit)
    mkdir -p "$ROOT/logs" "$OUT/parts"
    jid=$(sbatch --parsable -J p08_extract --array=0-5 -o "$ROOT/logs/p08_extract_%A_%a.out" \
          "$PIPE/run.sh" extract-array)
    echo "extract array: $jid"
    sbatch -J p08_merge -d "afterok:$jid" -o "$ROOT/logs/p08_merge_%j.out" --mem=16G -t 00:30:00 \
          "$PIPE/run.sh" merge
    ;;
  extract-array)
    set -- $PREDICTORS; shift "$SLURM_ARRAY_TASK_ID"; exec bash "$PIPE/run.sh" extract "$1" ;;
  extract)
    b=$2
    echo "host=$(hostname) job=${SLURM_JOB_ID:-none} predictor=$b start=$(date -Is)"
    $PY -u "$PIPE/extract_price149_all6.py" --manifest "$MANIFEST" --predictors "$b" \
        --out "$OUT/parts/$b.pkl.gz" --sha256-tsv "$OUT/parts/$b.sha256.tsv"
    echo "end=$(date -Is)"
    ;;
  merge)
    parts=(); for b in $PREDICTORS; do parts+=("$OUT/parts/$b.pkl.gz"); done
    $PY -u "$PIPE/extract_price149_all6.py" --merge "${parts[@]}" --out "$OUT/price149_all6_compact.pkl.gz"
    { head -1 "$OUT/parts/clean.sha256.tsv"
      for b in $PREDICTORS; do tail -n +2 "$OUT/parts/$b.sha256.tsv" | grep -v '^manifest'; done
      grep '^manifest' "$OUT/parts/clean.sha256.tsv"; } > "$OUT/input_sha256.tsv"
    sha256sum "$OUT/price149_all6_compact.pkl.gz" > "$OUT/price149_all6_compact.pkl.gz.sha256"
    ;;
  fetch)
    REPO=$(cd "$(dirname "$0")/../../.." && pwd)
    DEST=$REPO/benchmark/manifests/reference_data/ibex_extracted
    [ -f "$DEST/price149_all6_compact.pkl.gz.old" ] || cp -p "$DEST/price149_all6_compact.pkl.gz" "$DEST/price149_all6_compact.pkl.gz.old"
    rsync -av "ibex:$OUT/price149_all6_compact.pkl.gz" "$DEST/"
    rsync -av "ibex:$OUT/input_sha256.tsv" "$(dirname "$0")/input_sha256.tsv"
    ;;
  eval)
    REPO=$(cd "$(dirname "$0")/../../.." && pwd)
    cd "$REPO"
    python3 scripts/eval_price149_real_all6.py \
        --compact benchmark/manifests/reference_data/ibex_extracted/price149_all6_compact.pkl.gz
    ;;
  *)
    sed -n '2,15p' "$0"; exit 1 ;;
esac
