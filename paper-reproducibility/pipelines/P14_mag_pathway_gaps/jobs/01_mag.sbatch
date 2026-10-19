#!/bin/bash
#SBATCH --job-name=p14-mag
#SBATCH --partition=batch
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
# One MAG per array task (index = data row of $LIST, 0-based):
#   prep (layout/operons/gaps) -> DIAMOND + HMMER -> claims -> GSPA Stage 0 and
#   Stage 1+2A(+RLGC) -> score.  Idempotent: finished steps are skipped.
#   RESCORE=1 re-runs only the score step on finished MAGs.
set -euo pipefail
OUT=${OUT:-/ibex/scratch/projects/c2014/hohndor/mags_repro_2026-10-19}
source $OUT/code/jobs/common.sh
module load python/3.11.0 diamond/2.1.16 hmmer/3.3.2 >/dev/null 2>&1
LIST=${MAG_LIST:-$OUT/mag_list.tsv}
W=${WORK_DIR:-$WORK}
i=${SLURM_ARRAY_TASK_ID}
row=$(awk -F'\t' -v i=$i 'NR>1 && NR-2==i' $LIST)
[[ -n "$row" ]] || { echo "no MAG for index $i"; exit 0; }
hdr=$(head -1 $LIST)
col() { paste <(tr '\t' '\n' <<< "$hdr") <(tr '\t' '\n' <<< "$row") | awk -F'\t' -v c="$1" '$1==c{print $2}'; }
G=$(col mag_id); BIN=$(col bin); FA=$(col mag_fasta); BD=$(col bakta_dir); GD=$(col gapseq_dir)
TAXON=$(col taxon); PROF=$(col essential_profile)
T=${SLURM_CPUS_PER_TASK:-4}
D=$W/$G; mkdir -p $D/times; cd $D
TV="/usr/bin/time -v"
run() { local n=$1; shift; $TV -o times/$n.time "$@"; }
score() {
  run score python3 $PY score --mag-id $G --prep-dir $D --stage0 ${G}_stage0.tsv --stage12 ${G}_stage12.tsv \
      --suggestions ${G}_rlgc_suggestions.tsv --ec2go $REF/ec2go.txt --go-dag $REF/go-dag.tsv \
      --essential $REF/essential/$PROF.tsv --out ${G}_score.tsv --fills-out ${G}_fills.tsv
}
if [[ -s $D/DONE ]]; then
  if [[ "${RESCORE:-0}" == "1" ]]; then echo "[rescore] $G"; score; exit 0; fi
  echo "[skip] $G done"; exit 0
fi
echo "== $G ($TAXON, $PROF) $(date -Is) host=$(hostname) cpus=$T"
FAA=$BD/$BIN.faa
sha256sum $FA $FAA $BD/$BIN.gff3 $BD/$BIN.fna $GD/$G.xml $GD/$G-draft.xml $GD/$G-all-Pathways.tbl $GD/$G-all-Reactions.tbl > input_sha256.txt
{ diamond version; hmmsearch -h | sed -n 2p; } > tool_versions.txt 2>&1 || true

run prep python3 $PY prep --mag-id $G --bin $BIN --mag-fasta $FA --bakta-dir $BD --gapseq-dir $GD \
    --ec2go $REF/ec2go.txt --out-dir $D > prep.log 2>&1
cat prep.log
[[ -s ${G}_diamond.tsv ]] || { run diamond diamond blastp --query $FAA --db $REF/sprot_exp.dmnd --threads $T \
    --outfmt 6 --evalue 1e-5 --max-target-seqs 5 --quiet --out ${G}_diamond.tsv.part \
    && mv ${G}_diamond.tsv.part ${G}_diamond.tsv; }
[[ -s ${G}_pfam.domtbl ]] || { run hmmsearch hmmsearch --domtblout ${G}_pfam.domtbl.part --noali --cpu $T --cut_ga \
    $PFAM_HMM $FAA > /dev/null && mv ${G}_pfam.domtbl.part ${G}_pfam.domtbl; }
run claims python3 $BGC build-claims --obo $REF/go.obo --diamond ${G}_diamond.tsv --sprot-go $REF/sprot_exp_go.tsv \
    --domtbl ${G}_pfam.domtbl --pfam2go $REF/pfam2go.txt --out ${G}_claims.jsonl
JOPT="-Xmx14g -XX:ActiveProcessorCount=$T"
# Stage 0: noisy-OR only (ablation A1)
[[ -s ${G}_stage0.tsv ]] || run stage0 $JAVA $JOPT -jar $JAR integrate --claims ${G}_claims.jsonl \
    --go-owl $REF/go.owl --lite --enable-priors none --out ${G}_stage0.tsv > stage0.log 2>&1
# Stage 1+2A (A9) + gapseq gaps + RLGC
[[ -s ${G}_stage12.tsv ]] || run stage12 $JAVA $JOPT -jar $JAR integrate --claims ${G}_claims.jsonl \
    --go-owl $REF/go.owl --lite \
    --enable-priors essentiality,coherence,consistency,anatomical,gap_filling,genomic_context \
    --essential-profile $PROF --taxon $TAXON \
    --pathways $REF/kegg_pathways.tsv --modules $REF/kegg_modules.tsv --ec2go $REF/ec2go.txt \
    --operons ${G}_operons.tsv --gaps ${G}_gaps.jsonl \
    --rlc-suggester --reaction-graph $REF/seed_reactions_corrected.tsv --diffusion-mets $REF/diffusion_mets.tsv \
    --reaction-ec-aliases $REF/seed_Enzyme_Class_Reactions_Aliases_unique.tsv --genome-layout ${G}_layout.tsv \
    --suggestions-out ${G}_rlgc_suggestions.tsv \
    --enforce-stage2 --enforce-consistency-mode relax --enforcement-actions-out ${G}_stage2_actions.tsv \
    --out ${G}_stage12.tsv > stage12.log 2>&1
grep -h "RLGC\|Metabolic gaps\|Produced\|Suggestions\|Operon\|layout" stage12.log || true
score
gzip -f ${G}_diamond.tsv ${G}_pfam.domtbl ${G}_claims.jsonl 2>/dev/null || true
for f in times/*.time; do printf '%s\t%s\t%s\n' $(basename $f .time) \
  "$(grep 'Elapsed (wall' $f | awk '{print $NF}')" "$(grep 'Maximum resident' $f | awk '{print $NF}')"; done > timing.tsv
printf '%s\t%s\t%s\n' "$G" "$(grep -c '^>' $FAA)" "$(date -Is)" > DONE
echo "DONE $G $(date -Is)"
