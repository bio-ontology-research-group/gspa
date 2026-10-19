#!/bin/bash
# P05 -- Six-eukaryote anatomical / developmental constraint stress test (Table 4, Fig. 4c-d).
# Host: unimatrix01 (SLURM, partition debug). Never run compute on the login node.
#
# Usage (from the unimatrix01 login node, after copying this folder + the entry point to $PIPE):
#   bash run.sh stage                 # (local machine) copy run.sh + entry point to $PIPE on unimatrix01
#   bash run.sh submit                # sbatch the full run (real + synthetic input modes, one job)
#   sbatch ... run.sh run [mode]      # what `submit` sbatches; mode = both (default) | real | synthetic
#   bash run.sh fetch                 # (local machine) copy the three result TSVs into results/
#
# Entry point: benchmark/unified/pillar3/run_pillar3_eukaryotic_anatomy.py (stdlib-only Python 3).
#SBATCH -p debug
#SBATCH --exclude=unimatrix,node007
#SBATCH -c 16
#SBATCH --mem=90G
#SBATCH -t 06:00:00

set -euo pipefail
PY=/usr/bin/python3
REPRO=/data/hohndor/gspa/proteomes/pillar3_unified/repro_2026-10-19
PIPE=$REPRO/P05_eukaryote_anatomy
WORK=$PIPE/work   # all P05 outputs: real/ synthetic/ jar_resources/ results/ run_info.json diagnostics.json
JAR=/data/hohndor/gspa/bin/gspa-unified.jar   # v2.1.0, Oct 17 12:48, sha256 e415d9c6...72a7
ENTRY=run_pillar3_eukaryotic_anatomy.py

mode=${1:-help}
case "$mode" in
  stage)   # run on the local checkout
    here=$(cd "$(dirname "$0")" && pwd); root=$(cd "$here/../../.." && pwd)
    ssh unimatrix01 "mkdir -p $PIPE $REPRO/slurm"
    scp "$here/run.sh" "$root/benchmark/unified/pillar3/$ENTRY" unimatrix01:$PIPE/
    ;;
  submit)
    mkdir -p "$REPRO/slurm"
    # node007 has no /data mount -> excluded. Resources can be overridden, e.g.
    #   SBATCH_EXTRA="-c 8 --mem=15G" WORKERS=1 XMX=11g bash run.sh submit
    sbatch -J p05_anat --exclude=unimatrix,node007 ${SBATCH_EXTRA:-} --export=ALL \
      -o "$REPRO/slurm/p05_anat_%j.out" -e "$REPRO/slurm/p05_anat_%j.err" \
      "$PIPE/run.sh" run "${2:-both}"
    ;;
  run)
    echo "host=$(hostname) job=${SLURM_JOB_ID:-none} start=$(date -Is)"
    sha256sum "$JAR"
    $PY "$PIPE/$ENTRY" --input-source "${2:-both}" --work-dir "$WORK" --results-dir "$WORK/results" \
        --jar "$JAR" --workers "${WORKERS:-6}" --xmx "${XMX:-12g}"
    echo "end=$(date -Is)"
    ;;
  fetch)   # run on the local checkout
    here=$(cd "$(dirname "$0")" && pwd); root=$(cd "$here/../../.." && pwd)
    for f in pillar3_eukaryote_anatomy.tsv pillar3_relax_vs_remove.tsv pillar3_eukaryote_endogenous_vs_injected.tsv; do
      scp unimatrix01:$WORK/results/$f "$root/results/$f"
    done
    ;;
  *)
    sed -n '2,11p' "$0"
    ;;
esac
