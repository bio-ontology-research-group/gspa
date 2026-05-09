#!/bin/bash
#SBATCH --job-name=gspa-full3
#SBATCH --partition=debug
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/data/hohndor/gspa-tutorial-MR59-6/scripts/full3-%j.log
set -eo pipefail
WORK=/data/hohndor/gspa-tutorial-MR59-6
JAVA=/data/hohndor/envs/jdk-21.0.10+7/bin/java
JAR=/data/hohndor/gspa/bin/gspa-1.5.0.jar
REF=/data/hohndor/gspa/reference
OUT=${WORK}/gspa_full3_out
mkdir -p ${OUT}

# 0. Reuse the merged claims.jsonl from the previous run — same predictor
#    outputs, same parser. The fix this round is in gspa.jar (ClaimExtractor
#    now keeps mdf/proteinfer/clean instead of silently dropping them).
cp -n ${WORK}/gspa_full2_out/claims.jsonl ${OUT}/claims.jsonl

# 1. Predict operons from the Prokka GFF (intergenic distance + same strand).
echo "[$(date +%T)] predicting operons"
python3 ${WORK}/scripts/predict_operons.py \
    --gff ${WORK}/prokka_out/MR59-6.gff \
    --out-dir ${OUT}

# 2. gspa integrate — new jar, claims unchanged, plus operon membership
#    so GenomicContextPrior can fire on co-operonic transfers.
echo "[$(date +%T)] gspa integrate (new jar, with operons)"
${JAVA} -Xmx48g -jar ${JAR} integrate \
    --claims ${OUT}/claims.jsonl \
    --go-owl ${REF}/go-plus.owl \
    --pathways ${REF}/kegg_pathways.tsv \
    --ec2go ${REF}/ec2go.txt \
    --essential-profile bacteria \
    --enable-priors essentiality,coherence,gap_filling,genomic_context \
    --reasoner-cache ${REF}/reasoner-cache \
    --operons ${OUT}/operons_for_integrate.tsv \
    --out ${OUT}/integrated.tsv \
    --provenance ${OUT}/provenance.json 2>&1 | tail -25

echo "[$(date +%T)] DONE"
echo "  integrated rows: $(wc -l < ${OUT}/integrated.tsv)"
echo "  unique sources in supporting_sources:"
python3 -c "
import json, collections
c = collections.Counter()
with open('${OUT}/provenance.json') as f:
    for prov in json.load(f):
        for s in prov.get('supporting_sources', []):
            c[s] += 1
for s, n in c.most_common(): print(f'    {s:18s} {n:>10,}')
"
