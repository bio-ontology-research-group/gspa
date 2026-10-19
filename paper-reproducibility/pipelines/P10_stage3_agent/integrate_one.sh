#!/usr/bin/env bash
# Run one `gspa integrate` (P10 configuration) through the cap-lifting driver.
# Args: run_id taxon essential_profile claims operons gaps outdir
# operons / gaps may be '-'. A gaps file switches on the dark-matter suggester (unless
# P10_DM=0); its suggestions become extractor candidates.
# P10 configuration = the 88-proteome "iter3_relax" configuration
# (theta max_iter=3, Stage 2 consistency enforcement in relax mode, all priors,
# GO OWL lite, KEGG pathways + modules + ec2go, essential profile, asserted taxon).
set -uo pipefail
run_id=$1; taxon=$2; prof=$3; claims=$4; operons=$5; gaps=$6; out=$7
HERE=${P10_PIPE:-$(cd "$(dirname "$0")" && pwd)}
JAVA=${JAVA:-/data/hohndor/envs/jdk-21.0.10+7/bin/java}
JAR=${JAR:-/data/hohndor/gspa/bin/gspa-unified.jar}
REF=${REF:-/data/hohndor/gspa/reference}
XMX=${XMX:-10g}
mkdir -p "$out"
args=(--claims "$claims" --out "$out/integrated.tsv" --theta "$HERE/theta_iter3.json"
      --go-owl "$REF/go.owl" --lite --essential-profile "$prof"
      --pathways "$REF/kegg_pathways.tsv" --modules "$REF/kegg_modules.tsv" --ec2go "$REF/ec2go.txt"
      --taxon "$taxon"
      --provenance "$out/provenance.json" --anomalies-out "$out/anomalies.jsonl")
# P10_STAGE2=0: Stage 1 only (used for the ONTOLOGY_CONFLICT injections, which Stage 2 would
# otherwise relax before the extractor sees them)
[[ "${P10_STAGE2:-1}" == "1" ]] && args+=(--enforce-stage2 --enforce-consistency-mode relax)
[[ "$operons" != "-" ]] && args+=(--operons "$operons")
if [[ "$gaps" != "-" ]]; then
  args+=(--gaps "$gaps")
  # P10_DM=0 reproduces the batch (88-proteome) configuration, which has no dark-matter suggester
  [[ "${P10_DM:-1}" == "1" ]] && args+=(--dark-matter --suggestions-out "$out/suggestions.tsv")
fi
t0=$(date +%s)
P10_MAX_ANOMALIES=${P10_MAX_ANOMALIES:-100000} "$JAVA" -Xmx$XMX -cp "$JAR" groovy.ui.GroovyMain \
    "$HERE/p10_integrate.groovy" "${args[@]}" > "$out/integrate.log" 2>&1
rc=$?
sed -i '/DEBUG/d;/SLF4J/d' "$out/integrate.log"
# Edited claim copies (holdout runs) are regenerable from the case files; record their sha256 and
# optionally delete them to save disk (P10_CLEAN_CLAIMS=1).
if [[ "$claims" == "$(dirname "$out")/claims.jsonl" ]]; then
  sha256sum "$claims" > "$out/claims.sha256"
  [[ "${P10_CLEAN_CLAIMS:-0}" == "1" && $rc -eq 0 ]] && rm -f "$claims"
fi
echo -e "$run_id\t$rc\t$(( $(date +%s) - t0 ))\t$(hostname)" | tee -a "$out/run_status.tsv"
exit $rc
