#!/bin/bash
# Smoke test: build a tiny synthetic gaps.jsonl and verify that
# GapFillingPrior and GenomicContextPrior actually fire end-to-end.
set -u
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
JAVA=/data/hohndor/envs/java21/bin/java
JAR=/data/hohndor/gspa/bin/gspa.jar
BENCH=/data/hohndor/gspa/proteomes/bench9
REF=/data/hohndor/gspa/reference
ROOT=/data/hohndor/gspa/proteomes

TAG=hpylori
OUT=${BENCH}/smoke
mkdir -p ${OUT}

# Fake gaps targeting 3 real GO terms from the hpylori truth set
cat > ${OUT}/${TAG}_fake_gaps.jsonl <<EOF
{"pathway_id":"TEST:GLYCOLYSIS","reaction_id":"RXN-HEXOKINASE","ec_number":"EC:2.7.1.1","go_term":"GO:0004396","gapseq_guessed":false}
{"pathway_id":"TEST:TCA","reaction_id":"RXN-ACONITASE","ec_number":"EC:4.2.1.3","go_term":"GO:0003994","gapseq_guessed":false}
{"pathway_id":"TEST:RESP","reaction_id":"RXN-NADH","ec_number":"EC:7.1.1.2","go_term":"GO:0008137","gapseq_guessed":true}
EOF

${JAVA} -jar ${JAR} integrate \
  --claims ${BENCH}/${TAG}_claims.jsonl \
  --out ${OUT}/${TAG}_smoke_integrated.tsv \
  --go-owl ${REF}/go-plus.owl --lite \
  --essential-profile bacteria \
  --pathways ${REF}/kegg_pathways.tsv \
  --ec2go ${REF}/ec2go.txt \
  --reasoner-cache ${REF}/reasoner-cache \
  --operons ${ROOT}/operons/${TAG}_operons.tsv \
  --gaps ${OUT}/${TAG}_fake_gaps.jsonl \
  --enable-priors essentiality,coherence,consistency,gap_filling,genomic_context \
  2>&1 | grep -E 'Pathways|Operons|gaps|taxon|priors_fired|SAT|Converged|Refining|Essential|Loaded' | tail -20

echo "--- rows with non-empty priors_fired ---"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if ($i=="priors_fired") col=i; next} $col!="" && $col!="priors_fired" {print $1"\t"$3"\t"$col}' ${OUT}/${TAG}_smoke_integrated.tsv | head -20
echo "--- total rows with priors fired ---"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if ($i=="priors_fired") col=i; next} $col!="" {n++} END {print n}' ${OUT}/${TAG}_smoke_integrated.tsv
