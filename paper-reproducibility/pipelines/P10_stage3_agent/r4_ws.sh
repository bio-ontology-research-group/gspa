#!/usr/bin/env bash
# Round 4 on ws: decoy variant for the PATHWAY_GAP decoy half of the tp plan (plan.json).
set -e
source "$(dirname "$0")/ws_env.sh"
bash $P10_PIPE/holdout_round.sh r4dec decoy PATHWAY_GAP $P10_ROOT/cases/decoy_include_tp_pg.txt ${NPAR:-12} 8
echo "[r4] done $(date -Is)"
