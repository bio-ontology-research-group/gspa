#!/usr/bin/env bash
# One hold-out round: write edited inputs for <variant> on the selected cases, then run them.
#   holdout_round.sh <round_name> <variant tp|tpw|decoy> <categories|all> <include_file|-> <n_par> [batch_size]
# Environment (set by node007_job.sbatch or run.sh): P10_PIPE, JAVA, JAR, REF, XMX; P10_ROOT = mirror root of
# /data/hohndor/gspa/proteomes/pillar4_repro_2026-10-19 on the executing node.
set -euo pipefail
round=$1; variant=$2; cats=$3; inc=$4; npar=$5; bs=${6:-8}
ROOT=${P10_ROOT:?}
B="$P10_PIPE/repro/benchmark/unified/pillar4/build_pillar4_dataset_v2.py"
[[ -f "$B" ]] || B="$P10_PIPE/repo/benchmark/unified/pillar4/build_pillar4_dataset_v2.py"
args=(variants --genomes "$P10_PIPE/genomes_v2_node.tsv" --go-obo "$REF/go.obo" --cases "$ROOT/cases/cases_all.jsonl"
      --variant "$variant" --round "$round" --baseline-root "$ROOT/baseline" --batch-size "$bs"
      --work-root "$ROOT/holdouts" --runlist "$ROOT/runlists/runlist_${round}.tsv")
[[ "$cats" != "all" ]] && args+=(--categories "$cats")
[[ "$inc" != "-" ]] && args+=(--include "$inc")
mkdir -p "$ROOT/runlists" "$ROOT/holdouts"
python3 "$B" "${args[@]}"
grep -v '^run_id' "$ROOT/runlists/runlist_${round}.tsv" | xargs -P "$npar" -L 1 "$P10_PIPE/integrate_one.sh"
echo "[holdout_round] $round done $(date -Is)"
