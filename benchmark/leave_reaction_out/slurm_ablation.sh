#!/usr/bin/env bash
#SBATCH -A kaust
#SBATCH -p debug
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 24:00:00
#SBATCH -J lro_ablation
#SBATCH -o /data/hohndor/gspa/proteomes/bench_ecoli/leave_reaction_out/logs/lro_%A_%a.out
#SBATCH -e /data/hohndor/gspa/proteomes/bench_ecoli/leave_reaction_out/logs/lro_%A_%a.err
#
# Usage:
#   sbatch --array=0-4 slurm_ablation.sh <mode>    # mode=protein or function
#
# 5 shards across the test cases; each shard handles cases [shard*N/5, (shard+1)*N/5).

set -euo pipefail

MODE=${1:-protein}
ROOT=/data/hohndor/gspa/proteomes/bench_ecoli
LRO=$ROOT/leave_reaction_out
CASES=$LRO/cases_mg1655.tsv
TAG=mg1655

JAR=/data/hohndor/gspa/bin/gspa-phase11.jar
JAVA=/storage/miniforge3/envs/metagenomics/bin/java
GO_OWL=/data/hohndor/gspa/reference/go.owl
PATHWAYS=/data/hohndor/gspa/reference/kegg_pathways.tsv
EC2GO=/data/hohndor/gspa/reference/ec2go.txt

# Shard
N=$(tail -n +2 $CASES | wc -l)
SHARDS=${SLURM_ARRAY_TASK_COUNT:-5}
ID=${SLURM_ARRAY_TASK_ID:-0}
CHUNK=$(( (N + SHARDS - 1) / SHARDS ))
START=$(( ID * CHUNK ))
END=$(( START + CHUNK ))
if [[ $END -gt $N ]]; then END=$N; fi

OUT=$LRO/runs/${MODE}/shard_${ID}
mkdir -p $OUT $LRO/logs

echo "== mode=$MODE shard=$ID start=$START end=$END =="

python3 $LRO/run_ablation.py \
    --cases $CASES \
    --root $ROOT \
    --tag $TAG \
    --jar $JAR \
    --go-owl $GO_OWL \
    --pathways $PATHWAYS \
    --ec2go $EC2GO \
    --mode $MODE \
    --out-dir $OUT \
    --start $START --end $END \
    --java $JAVA
