#!/bin/bash
#SBATCH --job-name=eq-gapseq
#SBATCH --partition=batch
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
# One bin per array task: gapseq 1.4.0 draft + gap-filled model (find, find-transport, draft, medium,
# fill), unless the HQ run of 05_metabolism already holds a model for this bin (then it is linked).
# Output: $WORK/<mag>/gapseq/<mag>.xml, <mag>-draft.xml, <mag>-all-{Pathways,Reactions}.tbl, ...
set -euo pipefail
source ${OUT:-/encrypted/e3008/gspa/eq_mags}/code/common.sh
bin_row ${SLURM_ARRAY_TASK_ID} || { echo "no bin for index ${SLURM_ARRAY_TASK_ID}"; exit 0; }
T=${SLURM_CPUS_PER_TASK:-2}
D=$WORK/$G; mkdir -p $D/gapseq $D/times; cd $D
[[ -s DONE_gapseq ]] && { echo "[skip] $G gapseq done"; exit 0; }
echo "== $G ($CLASS, $DOMAIN) gapseq $(date -Is) host=$(hostname) cpus=$T"
if [[ -n "$GD" && -s $GD/$G.xml ]]; then
  for f in $GD/$G*; do ln -sfn "$f" "gapseq/$(basename "$f")"; done
  printf '%s\tlinked\t%s\t%s\n' "$G" "$GD" "$(date -Is)" > DONE_gapseq
  echo "linked existing gapseq model from $GD"; exit 0
fi
# genome FASTA (uncompressed, named <mag>.fna so gapseq names its outputs <mag>-*)
if [[ ! -s gapseq/$G.fna ]]; then
  if [[ -s annot/$G.fna ]]; then ln -sf ../annot/$G.fna gapseq/$G.fna; else zcat $FA > gapseq/$G.fna.part && mv gapseq/$G.fna.part gapseq/$G.fna; fi
fi
cd gapseq
export PATH=$(dirname $GAPSEQ):$PATH    # gapseq's draft/medium/fill call Rscript from its conda env
TAX=Bacteria; DRAFTB=$GRAM
[[ "$DOMAIN" == "Archaea" ]] && { TAX=Archaea; DRAFTB=archaea; }
TV="/usr/bin/time -v"; run() { local n=$1; shift; $TV -o ../times/$n.time "$@"; }
export OMP_NUM_THREADS=$T
[[ -s $G-all-Pathways.tbl && -s $G-all-Reactions.tbl ]] || run gapseq_find $GAPSEQ find -p all -t $TAX -K $T $G.fna > find.log 2>&1
[[ -s $G-Transporter.tbl ]] || run gapseq_transport $GAPSEQ find-transport -K $T $G.fna > transport.log 2>&1
[[ -s $G-draft.RDS ]] || run gapseq_draft $GAPSEQ draft -r $G-all-Reactions.tbl -t $G-Transporter.tbl -p $G-all-Pathways.tbl -c $G.fna -b $DRAFTB > draft.log 2>&1
[[ -s $G-medium.csv ]] || run gapseq_medium $GAPSEQ medium -m $G-draft.RDS -p $G-all-Pathways.tbl > medium.log 2>&1
[[ -s $G.xml ]] || run gapseq_fill $GAPSEQ fill -m $G-draft.RDS -n $G-medium.csv -c $G-rxnWeights.RDS -g $G-rxnXgenes.RDS > fill.log 2>&1
[[ -s $G.xml && -s $G-draft.xml ]] || { echo "gapseq did not produce $G.xml / $G-draft.xml"; exit 1; }
$GAPSEQ -v 2>&1 | head -3 > ../tool_versions_gapseq.txt || true
cd ..
for f in times/gapseq_*.time; do printf '%s\t%s\t%s\n' $(basename $f .time) \
  "$(grep 'Elapsed (wall' $f | awk '{print $NF}')" "$(grep 'Maximum resident' $f | awk '{print $NF}')"; done > timing_gapseq.tsv
printf '%s\tcomputed\t%s\t%s\n' "$G" "$D/gapseq" "$(date -Is)" > DONE_gapseq
echo "DONE_gapseq $G $(date -Is)"
