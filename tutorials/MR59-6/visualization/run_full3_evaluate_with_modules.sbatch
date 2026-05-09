#!/bin/bash
#SBATCH --job-name=gspa-full3-eval-mods
#SBATCH --partition=debug
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/data/hohndor/gspa-tutorial-MR59-6/scripts/full3-eval-mods-%j.log
set -eo pipefail
WORK=/data/hohndor/gspa-tutorial-MR59-6
JAVA=/data/hohndor/envs/jdk-21.0.10+7/bin/java
JAR=/data/hohndor/gspa/bin/gspa-1.5.0.jar
REF=/data/hohndor/gspa/reference
OUT=${WORK}/gspa_full3_out

# Re-run gspa evaluate with KEGG main pathways + KEGG Modules stacked.
# Modules add ~520 narrow units (5–15 enzymes each) on top of the
# 169 KEGG main pathways, so per-pathway enrichment for an operon
# actually has a chance of passing the strict bar (k>=2, coverage>=25%).
echo "[$(date +%T)] gspa evaluate (KEGG main + Modules)"
${JAVA} -Xmx24g -jar ${JAR} evaluate \
    --input ${WORK}/prokka_out/MR59-6.gff \
    --annotations ${OUT}/MR59-6.gspa.gaf \
    --go-owl ${REF}/go-plus.owl \
    --pathways ${REF}/kegg_pathways.tsv \
    --modules  ${REF}/kegg_modules.tsv \
    --ec2go ${REF}/ec2go.txt \
    --kingdom bacteria \
    --reasoner-cache ${REF}/reasoner-cache \
    --output ${OUT}/quality_gspa.json 2>&1 | tail -25

echo "[$(date +%T)] DONE"
python3 -c "
import json
q = json.load(open('${OUT}/quality_gspa.json'))
co = q['coherence']
detail = co.get('pathway_detail') or []
print(f'Pathways in quality JSON: {len(detail)}')
print(f'  KEGG main maps : {sum(1 for p in detail if p[\"id\"].startswith(\"KEGG:\") and not p[\"id\"].startswith(\"KEGG:M\"))}')
print(f'  KEGG Modules   : {sum(1 for p in detail if p[\"id\"].startswith(\"KEGG:M\"))}')
print(f'  Triggered (any enzyme present): {len(detail)}')
print(f'  Incoherent (<100% cov)         : {sum(1 for p in detail if p[\"completeness\"] < 1.0)}')
print(f'  100% covered                   : {sum(1 for p in detail if p[\"completeness\"] >= 1.0)}')
"
