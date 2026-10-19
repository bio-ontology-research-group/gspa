#!/bin/bash
# P14 — pathway gap resolution across uncultured Empty Quarter MAGs (ibex, SLURM).
#
# From the laptop/repo root:
#   bash paper-reproducibility/pipelines/P14_mag_pathway_gaps/run.sh sync     # rsync code to ibex
# then on ibex (sbatch only; nothing heavy runs on the login node):
#   bash $OUT/code/run.sh select    # 00: references snapshot + checksums + MAG selection (debug)
#   bash $OUT/code/run.sh test      # 3 MAGs (first three of mag_list.tsv) on debug + aggregate
#   bash $OUT/code/run.sh all       # 100-MAG array (%CONC concurrent) -> aggregate (afterany)
#   bash $OUT/code/run.sh resume    # resubmit only MAGs without work/<mag>/DONE, then aggregate
#   bash $OUT/code/run.sh rescore   # re-score finished MAGs only (MAG_LIST/WORK_DIR/TAG for the test set)
#   bash $OUT/code/run.sh aggregate # aggregate whatever has finished
# back on the laptop:
#   bash paper-reproducibility/pipelines/P14_mag_pathway_gaps/run.sh fetch
set -euo pipefail
export OUT=${OUT:-/ibex/scratch/projects/c2014/hohndor/mags_repro_2026-10-19}
HERE=$(cd "$(dirname "$0")" && pwd)
J=$OUT/code/jobs
SB="sbatch --parsable --account=pi-hohndor --export=ALL"
CONC=${CONC:-10}
# PART=debug runs the per-MAG array on the debug partition (1 h per MAG; a MAG takes ~6 min),
# which was the only partition scheduling when the batch CPU quota was saturated.
POPT=${PART:+-p $PART -t 01:00:00}

case "${1:-}" in
  sync)
    REPO=$(cd "$HERE/../../.." && pwd)
    ssh ibex "mkdir -p $OUT/code/jobs"
    rsync -a "$HERE/run.sh" "$HERE/README.md" "$HERE/OUTPUTS.md" "$HERE/mag_list.tsv" ibex:$OUT/code/ 2>/dev/null || \
      rsync -a "$HERE/run.sh" ibex:$OUT/code/
    rsync -a "$HERE/jobs/" ibex:$OUT/code/jobs/
    rsync -a "$REPO/benchmark/unified/pillar3/mag_pathway_gaps.py" \
             "$REPO/benchmark/unified/pillar2_fba_lro_price149/bgc_benchmark.py" ibex:$OUT/code/
    ;;
  select)
    mkdir -p $OUT/logs
    $SB -o $OUT/logs/select-%j.out $J/00_select.sbatch
    ;;
  test)
    mkdir -p $OUT/logs $OUT/test
    L=$OUT/test/mag_list_test.tsv
    { head -1 $OUT/mag_list.tsv; awk -F'\t' 'NR>1 && ($7=="Pseudomonadota" || $7=="Actinomycetota" || $6=="Archaea") && !seen[($6=="Archaea") ? "A" : $7]++' $OUT/mag_list.tsv | head -3; } > $L
    a=$($SB -p debug -t 02:00:00 --array=0-2 -o $OUT/logs/test-mag-%A_%a.out \
        --export=ALL,MAG_LIST=$L,WORK_DIR=$OUT/test/work $J/01_mag.sbatch)
    g=$($SB -p debug -d afterany:$a -o $OUT/logs/test-agg-%j.out \
        --export=ALL,MAG_LIST=$L,WORK_DIR=$OUT/test/work,TAG=test $J/02_aggregate.sbatch)
    echo "test jobs: mag-array=$a aggregate=$g"
    ;;
  all)
    mkdir -p $OUT/logs
    n=$(( $(wc -l < $OUT/mag_list.tsv) - 2 ))
    a=$($SB $POPT --array=0-$n%$CONC -o $OUT/logs/mag-%A_%a.out $J/01_mag.sbatch)
    g=$($SB -d afterany:$a -o $OUT/logs/agg-%j.out $J/02_aggregate.sbatch)
    echo "jobs: mag-array=$a aggregate=$g"
    ;;
  resume)
    idx=$(awk -F'\t' -v w=$OUT/work 'NR>1{ if (system("test -s " w "/" $1 "/DONE")) printf "%s%d", (k++?",":""), NR-2 }' $OUT/mag_list.tsv)
    [[ -n "$idx" ]] || { echo "all MAGs done"; exit 0; }
    a=$($SB $POPT --array=$idx%$CONC -o $OUT/logs/mag-%A_%a.out $J/01_mag.sbatch)
    g=$($SB -d afterany:$a -o $OUT/logs/agg-%j.out $J/02_aggregate.sbatch)
    echo "resubmitted indices $idx: mag-array=$a aggregate=$g"
    ;;
  rescore)   # re-run only the score step on finished MAGs (after a scoring change), then aggregate
    L=${MAG_LIST:-$OUT/mag_list.tsv}; WD=${WORK_DIR:-$OUT/work}; TG=${TAG:-v2}
    n=$(( $(wc -l < $L) - 2 ))
    a=$($SB -p debug -t 00:30:00 --mem=8G -c 1 --array=0-$n%$CONC -o $OUT/logs/rescore-%A_%a.out \
        --export=ALL,RESCORE=1,MAG_LIST=$L,WORK_DIR=$WD $J/01_mag.sbatch)
    g=$($SB -p debug -d afterany:$a -o $OUT/logs/agg-%j.out --export=ALL,MAG_LIST=$L,WORK_DIR=$WD,TAG=$TG $J/02_aggregate.sbatch)
    echo "rescore=$a aggregate=$g"
    ;;
  aggregate)
    $SB -o $OUT/logs/agg-%j.out $J/02_aggregate.sbatch
    ;;
  fetch)
    REPO=$(cd "$HERE/../../.." && pwd)
    rsync -a ibex:$OUT/results/pillar3_dark_matter_mags_v2.tsv ibex:$OUT/results/pillar3_dark_matter_mags_v2_per_mag.tsv "$REPO/results/"
    mkdir -p "$HERE/outputs"
    rsync -a ibex:$OUT/results/{mag_gap_fills_v2.tsv.gz,timing_v2.tsv,input_sha256_per_mag_v2.tsv} \
             ibex:$OUT/ref/{sha256.txt,tool_versions.txt} ibex:$OUT/{select.log,eligible_mags.tsv} "$HERE/outputs/"
    rsync -a ibex:$OUT/mag_list.tsv "$HERE/mag_list.tsv"
    ;;
  *) echo "usage: $0 {sync|select|test|all|resume|rescore|aggregate|fetch}"; exit 1;;
esac
