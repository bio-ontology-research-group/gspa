#!/bin/bash
# P12 step 02 (IBEX, one genome per Slurm array task): GSPA Stage 0 and
# Stage 1+2A on one RefSeq assembly.
#
#   Stage 0 channels : DIAMOND blastp vs reference_loo9 (Swiss-Prot minus the
#                      nine bench genomes; same DB and flags as pillar1_2),
#                      HMMER3 hmmsearch Pfam-A + pfam2go, InterProScan
#                      (re-used: annot_evaluation/preds_IPscan/<g>_IPscan.tsv)
#   Stage 0 output   : integrate --enable-priors none          (= ablation A1)
#   Stage 1+2A output: integrate, all six priors + --enforce-stage2
#                      --enforce-consistency-mode relax         (= ablation A9)
#                      operons from the RefSeq GFF3; NO gapseq gaps (not run).
#
# Every step runs under /usr/bin/time -v; times/<step>.time holds wall clock
# and max RSS.  usage: 02_gspa_genome.sh [GENOME]   (default: line
# SLURM_ARRAY_TASK_ID+OFFSET+1 of genomes.txt)
set -euo pipefail
R=${R:-/ibex/scratch/projects/c2014/hohndor/refseq1000_repro_2026-10-19}
C=/encrypted/e3008/annot_evaluation
JAVA=${JAVA:-/ibex/scratch/projects/c2014/rob/jdk-21/bin/java}
JAR=$R/bin/gspa-unified.jar
REF=$R/ref
P=$R/scripts
T=${SLURM_CPUS_PER_TASK:-4}
OFFSET=${OFFSET:-0}
type module >/dev/null 2>&1 || source /etc/profile.d/modules.sh
module load diamond/2.1.8 hmmer/3.3.2
command -v diamond >/dev/null && command -v hmmsearch >/dev/null

if [[ $# -ge 1 ]]; then G=$1; else
  G=$(sed -n "$(( SLURM_ARRAY_TASK_ID + OFFSET + 1 ))p" $R/genomes.txt); fi
[[ -n "$G" ]] || { echo "no genome"; exit 1; }
W=$R/work/$G
mkdir -p $W/preds $W/times
cd $W
[[ -s $W/DONE ]] && { echo "[skip] $G done"; exit 0; }
echo "[$(date +%T)] $G on $(hostname) cpus=$T"
FAA=$C/done_proteins/${G}_protein.faa
TAXON=$(awk -F'\t' -v g="$G" '$1==g{print $3}' $R/taxon_map.tsv); TAXON=${TAXON:-NCBITaxon_2}
TV="/usr/bin/time -v"

run() { local name=$1; shift; $TV -o times/$name.time "$@"; }

# --- Stage 0 channels -------------------------------------------------------
[[ -s preds/diamond.ok ]] || { run diamond diamond blastp --db $R/db/reference_loo9.dmnd \
    --query $FAA --out preds/diamond_results.tsv \
    --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore stitle \
    --evalue 1e-5 --max-target-seqs 50 --query-cover 50 --subject-cover 50 --id 30 \
    --threads $T > logs_diamond.txt 2>&1 && touch preds/diamond.ok; }
[[ -s preds/pfam.ok ]] || { run hmmsearch hmmsearch --domtblout preds/pfam_results.domtbl \
    --noali -E 1e-5 --cpu $T $R/db/Pfam-A.hmm $FAA > /dev/null && touch preds/pfam.ok; }
run parse_claims python3 $P/02b_parse_predictors_to_claims.py --results-dir preds \
    --goa $R/db/goa_reference_loo9_subset.gaf --pfam2go $REF/pfam2go.txt \
    --output claims_dp.jsonl > logs_parse.txt 2>&1
run parse_ipr python3 $P/parse_interproscan_claims.py --input $C/preds_IPscan/${G}_IPscan.tsv \
    --output claims_ipr.jsonl > logs_ipr.txt 2>&1
cat claims_dp.jsonl claims_ipr.jsonl > claims_merged.jsonl

OPS=()
[[ -s $R/operons/${G}_operons.tsv ]] && OPS=(--operons $R/operons/${G}_operons.tsv)

# --- Stage 0 (A1) ------------------------------------------------------------
run integrate_stage0 $JAVA -Xmx12g -XX:ActiveProcessorCount=$T -jar $JAR integrate \
    --claims claims_merged.jsonl --go-owl $REF/go.owl --out stage0_integrated.tsv --lite \
    --enable-priors none > logs_stage0.txt 2>&1

# --- Stage 1+2A (A9) ---------------------------------------------------------
run integrate_stage12 $JAVA -Xmx12g -XX:ActiveProcessorCount=$T -jar $JAR integrate \
    --claims claims_merged.jsonl --go-owl $REF/go.owl --out stage12_integrated.tsv --lite \
    --enable-priors essentiality,coherence,consistency,anatomical,gap_filling,genomic_context \
    --essential-profile bacteria --taxon $TAXON \
    --pathways $REF/kegg_pathways.tsv --modules $REF/kegg_modules.tsv --ec2go $REF/ec2go.txt \
    --enforce-stage2 --enforce-consistency-mode relax \
    --enforcement-actions-out stage12_actions.tsv "${OPS[@]}" > logs_stage12.txt 2>&1

gzip -f claims_dp.jsonl claims_ipr.jsonl; rm -f claims_merged.jsonl
gzip -f preds/diamond_results.tsv preds/pfam_results.domtbl
printf '%s\t%s\t%s\t%s\n' "$G" "$TAXON" "$(grep -c '^>' $FAA)" "$(date -Is)" > DONE
echo "[$(date +%T)] $G done"
