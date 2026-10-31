#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --job-name=atbp-S0
#SBATCH --no-requeue
# S0: pilot selection (00_select_pilot.py) plus a preflight that every reference, tool and code path the
# later steps read exists, and the reference receipt (sha256 of every reference file; the JAR checked against
# JAR_SHA256). Fails before any node is spent when something is missing. Needs HERE_CODE and OUT (driver.sh).
set -euo pipefail
source "${HERE_CODE:?}/common.sh"
CORES_PER_NODE=1   # shared partition: charged per core, not per node (2 CPUs = 1 core)
make_local; install_finish S0
mkdir -p "$LISTS" "$SHARDS_DIR" "$RECEIPTS" "$LOGS"
phase preflight start
missing=0
for p in "$ATB_CATALOG" "$ATB_MANIFESTS" "$ATB_ARCHIVES" "$JAR" "$GO_OWL" "$GO_OBO" "$GO_DAG" "$EC2GO" "$PFAM2GO" \
         "$KEGG_PATHWAYS" "$KEGG_MODULES" "$ESSENTIAL_DIR/bacteria.tsv" "$TAXON_LABELS" "$SEED_REACTIONS" "$DIFFUSION_METS" \
         "$EC_ALIASES" "$FBA_DATA/universal_metadata.json" "$FBA_DATA/fba_matrices_v6RE.pkl" "$SPROT_DMND" "$SPROT_GO" \
         "$SPROT_EC" "$PFAM_HMM" "$GAPSMITH" "$GAPSMITH_DATA/seq/Bacteria" "$TOOLS/bin/pyrodigal" "$TOOLS/bin/mmseqs" \
         "$TOOLS/bin/diamond" "$TOOLS/bin/hmmsearch" "$IPS_HOME/interproscan.sh" "$PY_HOST" "$IMAGE" \
         "$PKG_ROOT/home/leechuck/tmp/gspa-v3-testenv/bin/python" "$PKG_ROOT/usr/lib/jvm/java-21-openjdk-amd64/bin/java" \
         "$CODE/benchmark/unified/pillar2_fba_lro_price149/bgc_benchmark.py" "$CODE/benchmark/unified/pillar3/mag_pathway_gaps.py" \
         "$CODE/benchmark/emptyquarter_mags/eq_prep.py" "$CODE/benchmark/emptyquarter_mags/add_ec_claims.py" \
         "$CODE/scripts/stoichiometric_locus_refiner.py" "$HERE_CODE/receipt.py" "$HERE_CODE/atb_scatter.py"; do
  [[ -e "$p" ]] || { echo "MISSING $p"; missing=$((missing + 1)); }
done
jar_sha=$(sha256sum "$JAR" 2>/dev/null | cut -d' ' -f1 || true)
[[ "$jar_sha" == "$JAR_SHA256" ]] || { echo "JAR sha256 $jar_sha != expected $JAR_SHA256"; missing=$((missing + 1)); }
phase preflight end
(( missing == 0 )) || { echo "preflight: $missing problems, not selecting"; exit 2; }
phase select start
timed select select pilot "$PY_HOST" "$HERE_CODE/00_select_pilot.py" --catalog "$ATB_CATALOG" --manifests "$ATB_MANIFESTS" \
  --labels "$TAXON_LABELS" --seed "$SEED" --shards "$N_SHARDS" --genomes-per-shard "$GENOMES_PER_SHARD" --out-dir "$LISTS"
phase select end
RARGS+=(--count "n_genomes=$(( $(wc -l < "$LISTS/pilot_genomes.tsv") - 1 ))" --count "n_shards=$(( $(wc -l < "$LISTS/shards.tsv") - 1 ))"
        --param "seed=$SEED" --param "genomes_per_shard=$GENOMES_PER_SHARD" --param "chunk_size=$CHUNK_SIZE" --param "jar=$JAR"
        --input-sha "$LISTS/pilot_genomes.tsv" --input-sha "$LISTS/shards.tsv" --input-sha "$LISTS/selection.json")
# reference receipt (hashed when the EXIT trap writes the receipt): every reference file the pipeline reads (the IPS data tree and the gapsmith reference are
# recorded through their manifests; hashing their ~50 k files per run is not worth it)
for p in "$JAR" "$GO_OWL" "$GO_OBO" "$GO_DAG" "$EC2GO" "$PFAM2GO" "$KEGG_PATHWAYS" "$KEGG_MODULES" "$ESSENTIAL_DIR"/bacteria.tsv \
         "$ESSENTIAL_DIR"/archaea.tsv "$TAXON_LABELS" "$SEED_REACTIONS" "$DIFFUSION_METS" "$EC_ALIASES" "$FBA_DATA"/*.pkl \
         "$FBA_DATA/universal_metadata.json" "$SPROT_DMND" "$SPROT_GO" "$SPROT_EC" "$PFAM_HMM" "$GAPSMITH" \
         "$IPS_HOME/interproscan.properties" /project/k10018/gspa-v3/gapsmith-reference-bacteria-v1/gapsmith-ref-bacteria.sha256; do
  RARGS+=(--input-sha "$p")
done
version jar_sha256 echo "$jar_sha"
echo "S0 done: $(( $(wc -l < "$LISTS/pilot_genomes.tsv") - 1 )) genomes"
