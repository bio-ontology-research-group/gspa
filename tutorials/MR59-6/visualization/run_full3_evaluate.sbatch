#!/bin/bash
#SBATCH --job-name=gspa-full3-eval
#SBATCH --partition=debug
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/data/hohndor/gspa-tutorial-MR59-6/scripts/full3-eval-%j.log
set -eo pipefail
WORK=/data/hohndor/gspa-tutorial-MR59-6
JAVA=/data/hohndor/envs/jdk-21.0.10+7/bin/java
JAR=/data/hohndor/gspa/bin/gspa-1.5.0.jar
REF=/data/hohndor/gspa/reference
OUT=${WORK}/gspa_full3_out

# 1. integrated.tsv -> GAF (filter at posterior >= 0.5).
echo "[$(date +%T)] integrated -> GAF"
python3 ${WORK}/scripts/integrated_to_gaf.py \
    --integrated ${OUT}/integrated.tsv \
    --threshold 0.5 \
    --db GSPA --taxon "taxon:0" \
    --output ${OUT}/MR59-6.gspa.gaf

# 2. gspa evaluate -> quality_gspa.json with the new GAEF detail (named
#    missing essentials, incoherent process pairs, incoherent pathways).
echo "[$(date +%T)] gspa evaluate"
${JAVA} -Xmx24g -jar ${JAR} evaluate \
    --input ${WORK}/prokka_out/MR59-6.gff \
    --annotations ${OUT}/MR59-6.gspa.gaf \
    --go-owl ${REF}/go-plus.owl \
    --pathways ${REF}/kegg_pathways.tsv \
    --ec2go ${REF}/ec2go.txt \
    --kingdom bacteria \
    --reasoner-cache ${REF}/reasoner-cache \
    --output ${OUT}/quality_gspa.json 2>&1 | tail -25

echo "[$(date +%T)] DONE"
echo "  GAF rows: $(wc -l < ${OUT}/MR59-6.gspa.gaf)"
python3 -c "
import json
q = json.load(open('${OUT}/quality_gspa.json'))
s = q['summary']; c = q['completeness']; co = q['coherence']
print(f'  Coverage: {s[\"annotated_proteins\"]}/{s[\"total_proteins\"]} ({100*s[\"annotation_coverage\"]:.1f}%)')
print(f'  Completeness: {100*c[\"score\"]:.1f}% — {c[\"present_count\"]}/{c[\"present_count\"]+c[\"missing_count\"]} essentials')
print(f'  Process coherence: {100*co[\"process_coherence\"]:.1f}%')
print(f'  Pathway coherence: {100*co[\"pathway_coherence\"]:.1f}%')
print(f'  Composite: {s[\"composite_score\"]:.3f}')
print(f'  Has process_unsatisfied_pairs key: {bool(co.get(\"process_unsatisfied_pairs\"))} ({len(co.get(\"process_unsatisfied_pairs\") or [])} items)')
print(f'  Has pathway_detail key: {bool(co.get(\"pathway_detail\"))} ({len(co.get(\"pathway_detail\") or [])} pathways)')
print(f'  Has missing_functions_named key: {bool(c.get(\"missing_functions_named\"))} ({len(c.get(\"missing_functions_named\") or [])} entries)')
"
