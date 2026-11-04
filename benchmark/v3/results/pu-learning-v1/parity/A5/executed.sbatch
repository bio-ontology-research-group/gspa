#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-pu-replay
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/pu-learning-v1/slurm/%x-%j.log
#
# Production replay of cohort genomes with a PU-trained learned model (plan tasks/v3-row4-pu-learning-2026-11-02.md,
# P1 one-genome parity and metric 9), and the JVM/Python table parity check.
#   MODE=replay   gspa integrate --batch-manifest on the genomes of $JOBDIR/manifest.tsv (built on the login node by
#                 build_pu_cohort_inputs.py) with --learned-model $MODEL; flags mirror the archived family-excluded
#                 verification command (cohort inference, CLADE hold-out, feature export), engine = release-candidate-3.
#   MODE=parity   check_context_mode_models.groovy (unchanged) over $TABLE with $MODEL -> $JOBDIR/parity.jsonl.
# REPO_ROOT must be passed explicitly. Nothing runs on login nodes. No label is read here.
set -euo pipefail
: "${REPO_ROOT:?pass REPO_ROOT explicitly}"
: "${JOBDIR:?job directory under /scratch/project/k10018/gspa-v3/out/pu-learning-v1}"
: "${MODEL:?learned model artifact}"
: "${EXPECT_MODEL_SHA256:?}"
MODE=${MODE:-replay}
COHORT=${COHORT:-/project/k10018/gspa-v3/training-cohort}
REFERENCE=${REFERENCE:-/project/k10018/gspa-v3/price-production-v1/reference}
GSREF=${GSREF:-/project/k10018/gspa-v3/gapsmith-reference-bacteria-v1/reference}
JAR=${JAR:-/project/k10018/gspa-v3/release-candidate-3/gspa-3.0.0.jar}
EXPECT_JAR_SHA256=${EXPECT_JAR_SHA256:-8391d73fc1f6826492918a6116d40de4ba8147607ab406a895c29c32e067e80d}
FSROOT=${FSROOT:-/scratch/project/k10018/gspa-v3/shared/arm-fs-d04c9e11ee64db47}
JAVA21=$FSROOT/usr/lib/jvm/java-21-openjdk-amd64/bin
THREADS=${THREADS:-8}
HEAP=${HEAP:-48g}
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS PYTHONPATH PYTHONHOME
mkdir -p "$JOBDIR/tmp"
export TMPDIR="$JOBDIR/tmp"
hostname > "$JOBDIR/hostname.txt"; date --iso-8601=seconds > "$JOBDIR/started-at.txt"
echo "$EXPECT_JAR_SHA256  $JAR" | sha256sum -c - | tee "$JOBDIR/jar-check.log"
echo "$EXPECT_MODEL_SHA256  $MODEL" | sha256sum -c - | tee -a "$JOBDIR/jar-check.log"
"$JAVA21/java" -version 2>&1 | tee "$JOBDIR/java21.version"
cp "$0" "$JOBDIR/executed.sbatch"

if [ "$MODE" = "parity" ]; then
  : "${TABLE:?training table jsonl}"
  sha256sum "$TABLE" "$MODEL" "$REPO_ROOT/benchmark/v3/check_context_mode_models.groovy" > "$JOBDIR/parity-inputs.sha256"
  /usr/bin/time -v "$JAVA21/java" -Xmx8g -cp "$JAR" groovy.ui.GroovyMain "$REPO_ROOT/benchmark/v3/check_context_mode_models.groovy" \
      "$MODEL" "$TABLE" "$JOBDIR/parity.jsonl" > "$JOBDIR/parity.log" 2> "$JOBDIR/parity.time" \
    && echo "parity exit 0" | tee "$JOBDIR/exit.txt" || { echo "parity exit $?" | tee "$JOBDIR/exit.txt"; tail -n 30 "$JOBDIR/parity.log"; exit 1; }
  wc -l "$JOBDIR/parity.jsonl" "$TABLE"; date --iso-8601=seconds > "$JOBDIR/finished-at.txt"; exit 0
fi

[ -f "$JOBDIR/manifest.tsv" ] || { echo "manifest missing: $JOBDIR/manifest.tsv" >&2; exit 2; }
sha256sum "$JOBDIR/manifest.tsv" "$MODEL" > "$JOBDIR/inputs.sha256"
cut -f2 "$JOBDIR/manifest.tsv" | tail -n +2 | xargs sha256sum >> "$JOBDIR/inputs.sha256"
cmd=("$JAVA21/java" "-Xmx$HEAP" -jar "$JAR" integrate
  --batch-manifest "$JOBDIR/manifest.tsv" --batch-threads "$THREADS"
  --go-owl "$REFERENCE/go.owl" --ec2go "$REFERENCE/ec2go.txt"
  --pathways "$REFERENCE/kegg_pathways.tsv" --modules "$REFERENCE/kegg_modules.tsv"
  --reasoner-cache "$COHORT/reasoner-cache"
  --reaction-graph "$GSREF/seed_reactions.tsv" --diffusion-mets "$GSREF/diffusion_mets.tsv"
  --reaction-ec-aliases "$GSREF/seed_Enzyme_Class_Reactions_Aliases_unique.tsv"
  --cohort --cohort-metadata "$COHORT/cohort_inputs/cohort_metadata.json"
  --cohort-orthogroups "$COHORT/cohort_inputs/cohort_orthogroups.tsv" --cohort-holdout CLADE
  --learned-model "$MODEL"
  --feature-export-dir "$JOBDIR/features" --batch-stats-out "$JOBDIR/inference/batch_stats.tsv")
printf '%q ' "${cmd[@]}" > "$JOBDIR/command.txt"; echo >> "$JOBDIR/command.txt"
echo "== integrate: $(date)"
/usr/bin/time -v "${cmd[@]}" > "$JOBDIR/integrate.log" 2> "$JOBDIR/integrate.time" \
  && echo "integrate exit 0" | tee "$JOBDIR/exit.txt" || { echo "integrate exit $?" | tee "$JOBDIR/exit.txt"; tail -n 40 "$JOBDIR/integrate.log"; exit 1; }
ls "$JOBDIR/inference" | head -50
date --iso-8601=seconds > "$JOBDIR/finished-at.txt"
rm -rf "$JOBDIR/tmp"
echo "== done: $(date)"
