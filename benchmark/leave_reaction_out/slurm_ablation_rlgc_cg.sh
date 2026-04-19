#!/usr/bin/env bash
#SBATCH -A kaust
#SBATCH -p debug
#SBATCH -c 4
#SBATCH --mem=20G
#SBATCH -t 24:00:00
#SBATCH -J lro_rlgc_cg
#SBATCH -o /data/hohndor/gspa/proteomes/bench_ecoli/leave_reaction_out/logs/lro_rlgc_cg_%A_%a.out
#SBATCH -e /data/hohndor/gspa/proteomes/bench_ecoli/leave_reaction_out/logs/lro_rlgc_cg_%A_%a.err

set -euo pipefail

MODE=${1:-protein}
ROOT=/data/hohndor/gspa/proteomes/bench_ecoli
LRO=$ROOT/leave_reaction_out
CASES=$LRO/cases_mg1655.tsv
TAG=mg1655

JAR=/data/hohndor/gspa/bin/gspa-phase12.jar
JAVA=/storage/miniforge3/envs/metagenomics/bin/java
GO_OWL=/data/hohndor/gspa/reference/go.owl
PATHWAYS=/data/hohndor/gspa/reference/kegg_pathways.tsv
EC2GO=/data/hohndor/gspa/reference/ec2go.txt
GAPS_DATA=/data/hohndor/gspa/bin/gapsmith/data_merged
LAYOUT=/data/hohndor/gspa/proteomes/bench_gtdb30/layout/mg1655_layout.tsv
CATALOG=/data/hohndor/gspa/proteomes/bench_gtdb30/catalog_panel_excl_mg1655.tsv
ORTHO=/data/hohndor/gspa/proteomes/bench_gtdb30/ortho/orthogroup_map_50.tsv

N=$(tail -n +2 $CASES | wc -l)
SHARDS=${SLURM_ARRAY_TASK_COUNT:-5}
ID=${SLURM_ARRAY_TASK_ID:-0}
CHUNK=$(( (N + SHARDS - 1) / SHARDS ))
START=$(( ID * CHUNK ))
END=$(( START + CHUNK ))
if [[ $END -gt $N ]]; then END=$N; fi

OUT=$LRO/runs_rlgc_cg/${MODE}/shard_${ID}
mkdir -p $OUT $LRO/logs

echo "== RLGC+CG mode=$MODE shard=$ID start=$START end=$END =="
date

python3 $LRO/run_ablation_rlgc_cg.py \
    --cases $CASES \
    --root $ROOT \
    --tag $TAG \
    --jar $JAR \
    --go-owl $GO_OWL \
    --pathways $PATHWAYS \
    --ec2go $EC2GO \
    --reaction-graph $GAPS_DATA/seed_reactions.tsv \
    --diffusion-mets $GAPS_DATA/diffusion_mets.tsv \
    --ec-aliases $GAPS_DATA/seed_Enzyme_Class_Reactions_Aliases_unique.tsv \
    --genome-layout $LAYOUT \
    --rxn-catalog $CATALOG \
    --cg-lambda 1.0 \
    --orthogroup-map $ORTHO \
    --mode $MODE \
    --out-dir $OUT \
    --start $START --end $END \
    --java $JAVA
date
