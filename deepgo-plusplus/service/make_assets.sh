#!/bin/bash
# Build the runtime asset bundle for the DG++-Light webservice.
#
# Produces (into $OUT):
#   train_db.dmnd        DIAMOND DB of the pre-t0 train proteins (for diam + bridge)
#   train_net_index.tsv  precomputed STRING-neighbour vote per train node (the bridge)
#   train_terms.tsv      pre-t0 GO labels per train protein (for diam BLAST-KNN)
#   go-dag.tsv           child<TAB>ancestor closure (true-path propagation)
#   go.obo               GO ontology (term names; optional)
#
# train_net_index.tsv is the expensive, one-time precompute — it reads STRING once
# (~45 min). Reuse it across releases until STRING/the train set changes:
#   build_net_component.py --queries <train accessions with a string_id> \
#       --index text_string_index.tsv --train-terms train_terms.tsv \
#       --string-dir <STRING per-species dir> --out train_net_index.tsv
#
# Usage: make_assets.sh <OUT_DIR> <TRAIN_FASTA> <TRAIN_NET_INDEX> <TRAIN_TERMS> <GO_DAG> [GO_OBO]
set -euo pipefail
OUT="${1:?out dir}"; TRAIN_FASTA="${2:?train fasta}"; IDX="${3:?train_net_index.tsv}"
TT="${4:?train_terms.tsv}"; DAG="${5:?go-dag.tsv}"; OBO="${6:-}"
mkdir -p "$OUT"
echo "[make_assets] DIAMOND DB from $TRAIN_FASTA"
diamond makedb --in "$TRAIN_FASTA" -d "$OUT/train_db" --quiet
echo "[make_assets] copying index + labels + dag"
cp "$IDX" "$OUT/train_net_index.tsv"
cp "$TT"  "$OUT/train_terms.tsv"
cp "$DAG" "$OUT/go-dag.tsv"
[ -n "$OBO" ] && cp "$OBO" "$OUT/go.obo" || echo "[make_assets] (no go.obo — term names disabled)"
echo "[make_assets] done -> $OUT"
ls -la "$OUT"
