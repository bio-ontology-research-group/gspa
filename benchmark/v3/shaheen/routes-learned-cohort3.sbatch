#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-routes-learned-cohort3
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/routes-learned-v1/slurm/cohort3-%j.log
#
# Three cohort routes on three Price genomes with contract-bearing claims (plan tasks/v3-routes-learned-cohort-2026-11-01.md,
# check A6): CLI batch, Nextflow cohort.nf and the web service's cohort job, same JAR, same inputs, one node.
# REPO_ROOT must be passed explicitly (Slurm spools this script).
set -euo pipefail
trap 'exit 143' TERM INT
: "${REPO_ROOT:?pass REPO_ROOT explicitly (the staged code slice)}"
ROOT=${ROOT:-/scratch/project/k10018/gspa-v3/out/routes-learned-v1}
export JOB=${JOB:-$ROOT/cohort3}
export INPUTS=${INPUTS:-$ROOT/inputs}
export LAYOUT_METADATA=${LAYOUT_METADATA:-$INPUTS/layout_metadata.json}
export GENOMES=${GENOMES:-"GCF_000236665.1=640081 GCF_000015245.1=326297 GCA_003443575.1=221822"}
export JAR=${JAR:-$ROOT/engine/gspa-3.0.0.jar}
export EXPECT_JAR_SHA256=${EXPECT_JAR_SHA256:-08f8b5a132f4f1e5bf458f45b7604b1aa3b94388c96f170682bd2c27f2559c5a}
FSROOT=${FSROOT:-/scratch/project/k10018/gspa-v3/shared/arm-fs-d04c9e11ee64db47}
export JAVA_BIN=${JAVA_BIN:-$FSROOT/usr/lib/jvm/java-21-openjdk-amd64/bin}
export PYTHON=${PYTHON:-/project/k10018/sw/envs/gspa-v3/bin/python}
export GO_OWL=${GO_OWL:-/project/k10018/gspa-v3/price-production-v1/reference/go.owl}
export EC2GO=${EC2GO:-/project/k10018/gspa-v3/price-production-v1/reference/ec2go.txt}
export REACTION_GRAPH=${REACTION_GRAPH:-/project/k10018/gspa-v3/gapsmith-reference-bacteria-v1/reference/seed_reactions.tsv}
export REACTION_ALIASES=${REACTION_ALIASES:-/project/k10018/gspa-v3/gapsmith-reference-bacteria-v1/reference/seed_Enzyme_Class_Reactions_Aliases_unique.tsv}
export MODEL=${MODEL:-$ROOT/model_logistic.json}
export GO_OBO=${GO_OBO:-/project/k10018/gspa-v3/training-cohort/go-basic.obo}
export THREADS=${THREADS:-16}
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS PYTHONPATH PYTHONHOME
module load nextflow/25.10.4
export NEXTFLOW=nextflow
mkdir -p "$JOB"
sha256sum "$GO_OWL" "$EC2GO" "$REACTION_GRAPH" "$REACTION_ALIASES" "$MODEL" "$GO_OBO" > "$JOB/references.sha256"
exec bash "$REPO_ROOT/benchmark/v3/routes_closure/run_cohort_routes.sh"
