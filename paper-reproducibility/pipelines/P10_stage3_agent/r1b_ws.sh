#!/usr/bin/env bash
# Round 1b on ws: homology hold-out (tph) TP variant for INCOMPLETE_COMPLEX + BORDERLINE_ESSENTIAL
# (tp/tpw produced <=1 eligible case there), and the ONTOLOGY_CONFLICT injections with Stage 1 only.
set -e
source "$(dirname "$0")/ws_env.sh"
bash $P10_PIPE/holdout_round.sh r1tph tph INCOMPLETE_COMPLEX,BORDERLINE_ESSENTIAL - ${NPAR:-12} 8
bash $P10_PIPE/conflict_round.sh 7
