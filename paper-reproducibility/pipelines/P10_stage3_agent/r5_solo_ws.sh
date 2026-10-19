#!/usr/bin/env bash
# Round 5 (check): re-run 2 benchmark TP cases and 1 decoy case alone (batch size 1) to verify that
# packing non-interacting cases into one integrate run does not change the extractor record.
set -e
source "$(dirname "$0")/ws_env.sh"
bash $P10_PIPE/holdout_round.sh r5solotp tp PATHWAY_GAP $P10_ROOT/cases/solo_tp.txt 3 1
bash $P10_PIPE/holdout_round.sh r5solodec decoy PATHWAY_GAP $P10_ROOT/cases/solo_dec.txt 3 1
echo "[r5] done $(date -Is)"
