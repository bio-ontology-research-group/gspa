#!/usr/bin/env bash
# Round 2a on ws: enlarge the PATHWAY_GAP hold-out population (tpw on the next 45 cases per genome).
set -e
source "$(dirname "$0")/ws_env.sh"
while pgrep -f r1b_ws.sh > /dev/null; do sleep 30; done
bash $P10_PIPE/holdout_round.sh r2tpw tpw PATHWAY_GAP $P10_ROOT/cases/r2_include.txt ${NPAR:-12} 8
