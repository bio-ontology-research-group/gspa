#!/bin/bash
#SBATCH --job-name=eq-gapsmith
#SBATCH --partition=batch
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
# One genome per row (array tasks take blocks of rows, see below), after 01 of the same index (aftercorr; needs annot/<mag>.faa):
#   metabolic model for the gap prior, computed uniformly for every bin with gapsmith doall (find,
#   find-transport, draft, medium, fill; DIAMOND aligner; -t Archaea for archaeal bins) on the bin's proteome.
#   A gapseq 1.4.0 model of the HQ run (05_metabolism), when one exists, is linked under gapseq_ref/ for
#   comparison only (run.sh compare); it is not used as the prior, so that all bins share one model source.
# Output under $WORK/<mag>/gapseq/: <mag>.xml (gap-filled), <mag>-draft.xml, <mag>-all-{Pathways,Reactions}.tbl, ...
set -euo pipefail
# start from a clean module environment: the submitting shell may export LD_LIBRARY_PATH/PYTHONPATH of
# loaded modules (python/gcc), which breaks the conda tools and the static-linked gapsmith binary
module purge >/dev/null 2>&1 || true; unset LD_LIBRARY_PATH PYTHONPATH PYTHONHOME
source ${OUT:-/encrypted/e3008/gspa/eq_mags}/code/common.sh
# Array tasks work on blocks of GENOMES_PER_TASK consecutive rows of $BINS (ibex caps a user at 2,000 submitted
# jobs and every array task counts): the task re-runs this script once per row with ROW set.
if [[ -z "${ROW:-}" ]]; then
  K=${GENOMES_PER_TASK:-50}; fail=0; first=$(( SLURM_ARRAY_TASK_ID * K ))
  for (( r = first; r < first + K; r++ )); do ROW=$r bash "$0" || fail=$(( fail + 1 )); done
  echo "block ${SLURM_ARRAY_TASK_ID} (rows $first-$(( first + K - 1 ))): $fail failures $(date -Is)"; exit $(( fail > 0 ))
fi
bin_row $ROW || { [[ $ROW -ge $(( $(wc -l < $BINS) - 1 )) ]] && exit 0; echo "no genome for row $ROW"; exit 1; }
T=${SLURM_CPUS_PER_TASK:-4}
D=$WORK/$G; mkdir -p $D/gapseq $D/times; cd $D
# FORCE_COMPUTE=1 (validation): run gapsmith into gapseq_gapsmith/ even when a model can be linked,
# so that compare_models.py can set it against the gapseq model; the gapseq/ directory is left alone.
if [[ "${FORCE_COMPUTE:-0}" == "1" ]]; then
  # GM_TAG names the output set (gapseq_gapsmith$GM_TAG); TAX_OVERRIDE forces the reference taxonomy
  OD=gapseq_gapsmith${GM_TAG:-}
  [[ -s DONE_gapsmith_test${GM_TAG:-} ]] && { echo "[skip] $G gapsmith test done"; exit 0; }
  [[ -s annot/$G.faa ]] || { echo "$G: no proteome yet"; exit 1; }
  export PATH=$(dirname $DIAMOND_FOR_GAPSMITH):$PATH
  TAX=Bacteria; [[ "$DOMAIN" == "Archaea" ]] && TAX=Archaea; TAX=${TAX_OVERRIDE:-$TAX}
  mkdir -p $OD/run
  /usr/bin/time -v -o times/gapsmith_test${GM_TAG:-}.time $GAPSMITH doall annot/$G.faa -A diamond -K $T -t $TAX --data-dir $GAPSMITH_DATA -f $OD/run > $OD/gapsmith.log 2>&1 \
    || { echo "gapsmith failed"; tail -5 $OD/gapsmith.log; exit 1; }
  ls $OD/run; grep 'Elapsed (wall' times/gapsmith_test${GM_TAG:-}.time
  printf '%s\t%s\t%s\n' "$G" "$TAX" "$(date -Is)" > DONE_gapsmith_test${GM_TAG:-}; exit 0
fi
[[ -s DONE_gapseq ]] && { echo "[skip] $G model done"; exit 0; }
echo "== $G ($CLASS, $DOMAIN) model $(date -Is) host=$(hostname) cpus=$T"
if [[ -n "$GD" && -s $GD/$G.xml ]]; then
  mkdir -p gapseq_ref; for f in $GD/$G*; do ln -sfn "$f" "gapseq_ref/$(basename "$f")"; done
fi
[[ -s annot/$G.faa ]] || { echo "$G: no proteome yet (01 not finished)"; exit 1; }
export PATH=$(dirname $DIAMOND_FOR_GAPSMITH):$PATH
TAX=Bacteria; [[ "$DOMAIN" == "Archaea" ]] && TAX=Archaea
TV="/usr/bin/time -v"
mkdir -p gapseq/run
$TV -o times/gapsmith.time $GAPSMITH doall annot/$G.faa -A diamond -K $T -t $TAX --data-dir $GAPSMITH_DATA -f gapseq/run > gapseq/gapsmith.log 2>&1 \
  || { echo "gapsmith failed (see gapseq/gapsmith.log)"; tail -5 gapseq/gapsmith.log; exit 1; }
# gapsmith names its outputs after the FASTA stem ($G); expose the gap-filled model as <mag>.xml
for f in gapseq/run/$G*; do ln -sfn run/$(basename $f) gapseq/$(basename $f); done
ln -sfn run/$G-filled.xml gapseq/$G.xml
[[ -s gapseq/$G.xml && -s gapseq/$G-draft.xml && -s gapseq/$G-all-Reactions.tbl ]] || { echo "gapsmith outputs incomplete"; ls gapseq/run; exit 1; }
$GAPSMITH --version > tool_versions_gapsmith.txt 2>&1 || true
printf 'gapsmith\t%s\t%s\n' "$(grep 'Elapsed (wall' times/gapsmith.time | awk '{print $NF}')" "$(grep 'Maximum resident' times/gapsmith.time | awk '{print $NF}')" > timing_gapseq.tsv
printf '%s\tcomputed_gapsmith\t%s\t%s\n' "$G" "$D/gapseq/run" "$(date -Is)" > DONE_gapseq
echo "DONE_gapseq $G gapsmith $(date -Is)"
