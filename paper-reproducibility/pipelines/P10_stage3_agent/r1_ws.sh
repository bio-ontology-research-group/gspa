#!/usr/bin/env bash
# Round 1 on ws: TP variants (tp for PATHWAY_GAP, tpw for all categories) + conflict injections.
set -e
source "$(dirname "$0")/ws_env.sh"
bash $P10_PIPE/holdout_round.sh r1tp tp PATHWAY_GAP $P10_ROOT/cases/r1_include.txt ${NPAR:-12} 8
bash $P10_PIPE/holdout_round.sh r1tpw tpw all $P10_ROOT/cases/r1_include.txt ${NPAR:-12} 8
bash $P10_PIPE/conflict_round.sh 7
