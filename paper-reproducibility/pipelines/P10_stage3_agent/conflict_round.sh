#!/usr/bin/env bash
# ONTOLOGY_CONFLICT injection runs (one run per genome, 40 injected forbidden terms each).
#   conflict_round.sh <n_par>     (env as for holdout_round.sh)
set -euo pipefail
npar=${1:-7}
ROOT=${P10_ROOT:?}
B="$P10_PIPE/repo/benchmark/unified/pillar4/build_pillar4_dataset_v2.py"
mkdir -p "$ROOT/runlists" "$ROOT/holdouts"
python3 "$B" conflicts --genomes "$P10_PIPE/genomes_v2_node.tsv" --go-obo "$REF/go.obo" \
  --tags ecoli bsubtilis paeruginosa mtb scerevisiae celegans athaliana --n-runs 1 --n-inject 40 \
  --work-root "$ROOT/holdouts" --runlist "$ROOT/runlists/runlist_conf.tsv"
# Stage 2 consistency enforcement relaxes every injected violation before the extractor runs
# (observed in a first pass with Stage 2 on), so the injection runs use Stage 1 only.
grep -v '^run_id' "$ROOT/runlists/runlist_conf.tsv" | P10_STAGE2=0 xargs -P "$npar" -L 1 "$P10_PIPE/integrate_one.sh"
echo "[conflict_round] done $(date -Is)"
