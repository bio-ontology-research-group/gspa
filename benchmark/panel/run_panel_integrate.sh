#!/usr/bin/env bash
#SBATCH --job-name=panel-integrate
#SBATCH --partition=debug
#SBATCH --exclude=node003
#SBATCH --output=/data/hohndor/gspa/proteomes/bench_gtdb30/integrate-%A_%a.out
#SBATCH --error=/data/hohndor/gspa/proteomes/bench_gtdb30/integrate-%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
# Usage: sbatch --array=0-29 run_panel_integrate.sh

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics

ROOT=/data/hohndor/gspa/proteomes/bench_gtdb30
JAR=/data/hohndor/gspa/bin/gspa-phase12.jar
JAVA=/storage/miniforge3/envs/metagenomics/bin/java
GO_OWL=/data/hohndor/gspa/reference/go.owl
PATHWAYS=/data/hohndor/gspa/reference/kegg_pathways.tsv
EC2GO=/data/hohndor/gspa/reference/ec2go.txt

TAGS=($(tail -n +2 $ROOT/panel_manifest.tsv | cut -f1))
tag=${TAGS[$SLURM_ARRAY_TASK_ID]}

CLAIMS=$ROOT/claims/${tag}_dp_claims.jsonl
OUTDIR=$ROOT/integrated
OUT=$OUTDIR/${tag}_integrated.tsv
mkdir -p $OUTDIR

[[ -s $CLAIMS ]] || { echo "missing claims for $tag"; exit 0; }

echo "=== $tag: whole-genome integrate ==="
date
$JAVA -Xmx10g -jar $JAR integrate \
    --claims $CLAIMS --out $OUT \
    --go-owl $GO_OWL --lite \
    --essential-profile bacteria \
    --pathways $PATHWAYS --ec2go $EC2GO
wc -l $OUT
echo "DONE"
date
