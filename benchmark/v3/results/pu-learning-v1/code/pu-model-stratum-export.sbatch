#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-pu-stratum
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/pu-learning-v1/slurm/%x-%j.log
#
# Feature export for the model-genome evaluation stratum (plan tasks/v3-row4-pu-learning-2026-11-02.md §4.3):
# gspa integrate --feature-export-dir on E. coli K-12, B. subtilis 168, M. tuberculosis H37Rv with the family-excluded
# contract claims; no cohort inputs (these genomes are not cohort members), no learned model (scoring is done offline
# by every arm on the exported candidates). No label is read. Engine = release-candidate-3.
set -euo pipefail
JOBDIR=${JOBDIR:-/scratch/project/k10018/gspa-v3/out/pu-learning-v1/model-stratum/export}
REFERENCE=${REFERENCE:-/project/k10018/gspa-v3/price-production-v1/reference}
GSREF=${GSREF:-/project/k10018/gspa-v3/gapsmith-reference-bacteria-v1/reference}
COHORT=${COHORT:-/project/k10018/gspa-v3/training-cohort}
JAR=${JAR:-/project/k10018/gspa-v3/release-candidate-3/gspa-3.0.0.jar}
EXPECT_JAR_SHA256=${EXPECT_JAR_SHA256:-8391d73fc1f6826492918a6116d40de4ba8147607ab406a895c29c32e067e80d}
FSROOT=${FSROOT:-/scratch/project/k10018/gspa-v3/shared/arm-fs-d04c9e11ee64db47}
JAVA21=$FSROOT/usr/lib/jvm/java-21-openjdk-amd64/bin
THREADS=${THREADS:-3}
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS PYTHONPATH PYTHONHOME
mkdir -p "$JOBDIR/tmp" "$JOBDIR/inference" "$JOBDIR/features"
export TMPDIR="$JOBDIR/tmp"
hostname > "$JOBDIR/hostname.txt"; date --iso-8601=seconds > "$JOBDIR/started-at.txt"; cp "$0" "$JOBDIR/executed.sbatch"
echo "$EXPECT_JAR_SHA256  $JAR" | sha256sum -c - | tee "$JOBDIR/jar-check.log"
[ -f "$JOBDIR/manifest.tsv" ] || { echo "manifest missing" >&2; exit 2; }
sha256sum "$JOBDIR/manifest.tsv" > "$JOBDIR/inputs.sha256"; cut -f2,12 "$JOBDIR/manifest.tsv" | tail -n +2 | tr '\t' '\n' | xargs sha256sum >> "$JOBDIR/inputs.sha256"
cmd=("$JAVA21/java" -Xmx32g -jar "$JAR" integrate
  --batch-manifest "$JOBDIR/manifest.tsv" --batch-threads "$THREADS"
  --go-owl "$REFERENCE/go.owl" --ec2go "$REFERENCE/ec2go.txt"
  --pathways "$REFERENCE/kegg_pathways.tsv" --modules "$REFERENCE/kegg_modules.tsv"
  --reasoner-cache "$COHORT/reasoner-cache"
  --reaction-graph "$GSREF/seed_reactions.tsv" --diffusion-mets "$GSREF/diffusion_mets.tsv"
  --reaction-ec-aliases "$GSREF/seed_Enzyme_Class_Reactions_Aliases_unique.tsv"
  --feature-export-dir "$JOBDIR/features" --batch-stats-out "$JOBDIR/inference/batch_stats.tsv")
printf '%q ' "${cmd[@]}" > "$JOBDIR/command.txt"; echo >> "$JOBDIR/command.txt"
/usr/bin/time -v "${cmd[@]}" > "$JOBDIR/integrate.log" 2> "$JOBDIR/integrate.time" \
  && echo "integrate exit 0" | tee "$JOBDIR/exit.txt" || { echo "integrate exit $?" | tee "$JOBDIR/exit.txt"; tail -n 40 "$JOBDIR/integrate.log"; exit 1; }
ls -la "$JOBDIR/features"/*/ | head -20; date --iso-8601=seconds > "$JOBDIR/finished-at.txt"; rm -rf "$JOBDIR/tmp"
