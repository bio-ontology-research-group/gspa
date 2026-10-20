#!/bin/bash
#SBATCH --job-name=eq-integrate
#SBATCH --partition=batch
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
# One bin per array task, after 01 (claims) and 02 (gapseq) of the same index (aftercorr):
#   prep (layout, operons, gapseq gaps) -> GSPA Stage 0 (noisy-OR only) -> Stage 1 + 2A (all six priors,
#   relax enforcement, RLGC suggester) + Stage 2B (StoichiometricLocusRefiner, flux-certified repair)
#   -> P14 score (gap fills, pathway and essential completeness) when a gapseq model exists.
set -euo pipefail
source ${OUT:-/encrypted/e3008/gspa/eq_mags}/code/common.sh
module load python/3.11.0 >/dev/null 2>&1
bin_row ${SLURM_ARRAY_TASK_ID} || { echo "no bin for index ${SLURM_ARRAY_TASK_ID}"; exit 0; }
T=${SLURM_CPUS_PER_TASK:-4}
D=$WORK/$G; cd $D
[[ -s DONE ]] && { echo "[skip] $G integrated"; exit 0; }
[[ -s DONE_channels ]] || { echo "$G: channels not finished"; exit 1; }
[[ -s ${G}_claims.jsonl || -s ${G}_claims.jsonl.gz ]] || { echo "$G: no claims"; exit 1; }
[[ -s ${G}_claims.jsonl ]] || gunzip -k ${G}_claims.jsonl.gz
GSQ=""; [[ -s DONE_gapseq && -s gapseq/$G.xml ]] && GSQ=$D/gapseq
echo "== $G ($CLASS, $TAXON, gram $GRAM, profile $PROF, gapseq ${GSQ:-none}) $(date -Is) host=$(hostname) cpus=$T"
TV="/usr/bin/time -v"; run() { local n=$1; shift; $TV -o times/$n.time "$@"; }

run prep python3 $PREP --mag-id $G --annot-dir annot --mag-fasta $FA --ec2go $REF/ec2go.txt --out-dir $D \
    ${GSQ:+--gapseq-dir $GSQ --p14 $P14PY} > prep.log 2>&1
cat prep.log
GAPS=(); [[ -s ${G}_gaps.jsonl ]] && GAPS=(--gaps ${G}_gaps.jsonl)
JOPT="-Xmx14g -XX:ActiveProcessorCount=$T"
# Stage 0: noisy-OR pooling of the three channels (ablation A1)
[[ -s ${G}_stage0.tsv ]] || run stage0 $JAVA $JOPT -jar $JAR integrate --claims ${G}_claims.jsonl \
    --go-owl $REF/go.owl --lite --enable-priors none --out ${G}_stage0.tsv > stage0.log 2>&1
# Stage 1 + 2A + 2B (configuration A9 with the stoichiometric refiner)
[[ -s ${G}_stage12.tsv ]] || run stage12 $JAVA $JOPT -jar $JAR integrate --claims ${G}_claims.jsonl \
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
grep -h "RLGC\|Metabolic gaps\|Produced\|Suggestions\|Operon\|layout\|Stage 2B\|refiner" stage12.log | head -20 || true
if [[ -n "$GSQ" ]]; then
  run score python3 $P14PY score --mag-id $G --prep-dir $D --stage0 ${G}_stage0.tsv --stage12 ${G}_stage12.tsv \
      --suggestions ${G}_rlgc_suggestions.tsv --ec2go $REF/ec2go.txt --go-dag $REF/go-dag.tsv \
      --essential $REF/essential/$PROF.tsv --out ${G}_score.tsv --fills-out ${G}_fills.tsv > score.log 2>&1 || echo "score failed (see score.log)"
fi
gzip -f ${G}_claims.jsonl ${G}_stage0.tsv ${G}_stage12.tsv 2>/dev/null || true
rm -f ${G}_claims.jsonl
for f in times/*.time; do printf '%s\t%s\t%s\n' $(basename $f .time) \
  "$(grep 'Elapsed (wall' $f | awk '{print $NF}')" "$(grep 'Maximum resident' $f | awk '{print $NF}')"; done > timing.tsv
printf '%s\t%s\t%s\t%s\t%s\n' "$G" "$(grep -c '^>' annot/$G.faa)" "$TAXON" "${GSQ:+gapseq}" "$(date -Is)" > DONE
echo "DONE $G $(date -Is)"
