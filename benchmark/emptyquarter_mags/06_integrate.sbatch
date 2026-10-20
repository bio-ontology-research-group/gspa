#!/bin/bash
#SBATCH --job-name=eq-integrate
#SBATCH --partition=batch
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
# One genome per row (array tasks take blocks of rows, see below), after the scatter job (claims) and the gapsmith task of the same index (model):
#   claims of the bin from its sample file -> prep (layout, operons, gap reactions) -> GSPA Stage 0 ->
#   Stage 1 + 2A (six priors, relax enforcement, RLGC suggester) + Stage 2B (StoichiometricLocusRefiner)
#   -> P14 score (gap fills, pathway and essential completeness).
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
module load python/3.11.0 >/dev/null 2>&1
bin_row $ROW || { [[ $ROW -ge $(( $(wc -l < $BINS) - 1 )) ]] && exit 0; echo "no genome for row $ROW"; exit 1; }
T=${SLURM_CPUS_PER_TASK:-4}
D=$WORK/$G; mkdir -p $D/times; cd $D
[[ -s DONE ]] && { echo "[skip] $G integrated"; exit 0; }
[[ -s DONE_proteins ]] || { echo "$G: proteins not finished"; exit 1; }
[[ -s $CBS/DONE && -s $CBS/$SAMPLE.jsonl.gz ]] || { echo "$G: claims of sample $SAMPLE not scattered"; exit 1; }
GSQ=""; [[ -s DONE_gapseq && -s gapseq/$G.xml ]] && GSQ=$D/gapseq
echo "== $G ($CLASS, $TAXON, gram $GRAM, profile $PROF, model ${GSQ:-none}) $(date -Is) host=$(hostname) cpus=$T"
TV="/usr/bin/time -v"; run() { local n=$1; shift; $TV -o times/$n.time "$@"; }

# claims of this bin: lines whose protein_id starts with "<mag>|", prefix stripped
if [[ ! -s ${G}_claims.jsonl ]]; then
  zcat $CBS/$SAMPLE.jsonl.gz | grep -F "\"protein_id\":\"$G|" | sed "s/\"protein_id\":\"$G|/\"protein_id\":\"/" > ${G}_claims.jsonl.part
  mv ${G}_claims.jsonl.part ${G}_claims.jsonl
fi
NCLAIMS=$(wc -l < ${G}_claims.jsonl)
[[ $NCLAIMS -gt 0 ]] || echo "WARNING: $G has no claims"
run prep python3 $PREP --mag-id $G --annot-dir annot --mag-fasta $FA --ec2go $REF/ec2go.txt --out-dir $D \
    ${GSQ:+--gapseq-dir $GSQ --p14 $P14PY --ec-aliases $REF/seed_Enzyme_Class_Reactions_Aliases_unique.tsv} > prep.log 2>&1
cat prep.log
GAPS=(); [[ -s ${G}_gaps.jsonl ]] && GAPS=(--gaps ${G}_gaps.jsonl)
JOPT="-Xmx14g -XX:ActiveProcessorCount=$T"
[[ -s ${G}_stage0.tsv || -s ${G}_stage0.tsv.gz ]] || run stage0 $JAVA $JOPT -jar $JAR integrate --claims ${G}_claims.jsonl \
    --go-owl $REF/go.owl --lite --enable-priors none --out ${G}_stage0.tsv > stage0.log 2>&1
[[ -s ${G}_stage12.tsv || -s ${G}_stage12.tsv.gz ]] || run stage12 $JAVA $JOPT -jar $JAR integrate --claims ${G}_claims.jsonl \
    --go-owl $REF/go.owl --lite \
    --enable-priors essentiality,coherence,consistency,anatomical,gap_filling,genomic_context \
    --essential-profile $PROF --taxon $TAXON \
    --pathways $REF/kegg_pathways.tsv --modules $REF/kegg_modules.tsv --ec2go $REF/ec2go.txt \
    --operons ${G}_operons.tsv "${GAPS[@]}" \
    --rlc-suggester --reaction-graph $REF/seed_reactions_corrected.tsv --diffusion-mets $REF/diffusion_mets.tsv \
    --reaction-ec-aliases $REF/seed_Enzyme_Class_Reactions_Aliases_unique.tsv --genome-layout ${G}_layout.tsv \
    --suggestions-out ${G}_rlgc_suggestions.tsv \
    --enable-fba-refiner --fba-data-dir $REF/fba --fba-sidecar $FBA_SIDECAR --fba-python-bin $FBA_PY \
    --fba-gram $GRAM --fba-medium minimal \
    --enforce-stage2 --enforce-consistency-mode relax --enforcement-actions-out ${G}_stage2_actions.tsv \
    --out ${G}_stage12.tsv > stage12.log 2>&1
grep -h "RLGC\|Metabolic gaps\|Produced\|Suggestions\|Operon\|layout\|Stage 2B\|refiner\|Exception" stage12.log | head -20 || true
if [[ -n "$GSQ" ]]; then
  run score python3 $P14PY score --mag-id $G --prep-dir $D --stage0 ${G}_stage0.tsv --stage12 ${G}_stage12.tsv \
      --suggestions ${G}_rlgc_suggestions.tsv --ec2go $REF/ec2go.txt --go-dag $REF/go-dag.tsv \
      --essential $REF/essential/$PROF.tsv --out ${G}_score.tsv --fills-out ${G}_fills.tsv > score.log 2>&1 || echo "score failed (see score.log)"
fi
gzip -f ${G}_claims.jsonl ${G}_stage0.tsv ${G}_stage12.tsv 2>/dev/null || true
for f in times/*.time; do printf '%s\t%s\t%s\n' $(basename $f .time) \
  "$(grep 'Elapsed (wall' $f | awk '{print $NF}')" "$(grep 'Maximum resident' $f | awk '{print $NF}')"; done > timing.tsv
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$G" "$(grep -c '^>' annot/$G.faa)" "$NCLAIMS" "$TAXON" "${GSQ:+model}" "$(date -Is)" > DONE
echo "DONE $G claims=$NCLAIMS $(date -Is)"
