#!/usr/bin/env bash
# Round 2b on ws: decoy variant for the PATHWAY_GAP decoy half fixed by `plan` (plan.json).
set -e
source "$(dirname "$0")/ws_env.sh"
bash $P10_PIPE/holdout_round.sh r2dec decoy PATHWAY_GAP $P10_ROOT/cases/decoy_include.txt ${NPAR:-12} 8
echo "[r2b] done $(date -Is)"
