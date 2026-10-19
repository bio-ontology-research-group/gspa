#!/usr/bin/env bash
# Round 3 on ws: spec-literal tp variant (all claims for the target family deleted) on every
# PATHWAY_GAP case not yet run with tp. Reason: the tpw construction leaves the held-out protein
# with exactly one supporting claim, which separates TP from decoy records (leak check AUC 0.98).
set -e
source "$(dirname "$0")/ws_env.sh"
bash $P10_PIPE/holdout_round.sh r3tp tp PATHWAY_GAP $P10_ROOT/cases/r3_include.txt ${NPAR:-12} 8
echo "[r3] done $(date -Is)"
