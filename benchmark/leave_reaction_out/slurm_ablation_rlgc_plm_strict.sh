#!/usr/bin/env bash
#SBATCH -A kaust
#SBATCH -p debug
#SBATCH --exclude=node003
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 24:00:00
#SBATCH -J lro_plm_Astrict
#SBATCH -o /data/hohndor/gspa/proteomes/bench_ecoli/leave_reaction_out/logs/lro_plm_Astrict_%A_%a.out
#SBATCH -e /data/hohndor/gspa/proteomes/bench_ecoli/leave_reaction_out/logs/lro_plm_Astrict_%A_%a.err

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
set -u

MODE=protein
ROOT=/data/hohndor/gspa/proteomes/bench_ecoli
LRO=$ROOT/leave_reaction_out
CASES=$LRO/cases_mg1655.tsv

JAR=/data/hohndor/gspa/bin/gspa-phase12.jar
JAVA=/storage/miniforge3/envs/metagenomics/bin/java
GO_OWL=/data/hohndor/gspa/reference/go.owl
PATHWAYS=/data/hohndor/gspa/reference/kegg_pathways.tsv
EC2GO=/data/hohndor/gspa/reference/ec2go.txt
GAPS_DATA=/data/hohndor/gspa/bin/gapsmith/data_merged
LAYOUT=/data/hohndor/gspa/proteomes/bench_gtdb30/layout/mg1655_layout.tsv
MODEL=/data/hohndor/gspa/proteomes/bench_gtdb30/ml/model_m4.txt
CENTROIDS=/data/hohndor/gspa/proteomes/bench_gtdb30/plm_centroids
PLM_NPY=/data/hohndor/gspa/proteomes/bench_gtdb30/plm/mg1655_esm2t30.npy
PLM_IDX=/data/hohndor/gspa/proteomes/bench_gtdb30/plm/mg1655_esm2t30.index.tsv
PANEL_NPY=/data/hohndor/gspa/proteomes/bench_gtdb30/plm_catalysts/ec_panel.npy
CATALYSTS=/data/hohndor/gspa/proteomes/bench_gtdb30/plm_catalysts/ec_catalysts.tsv
PANEL_IDX=/data/hohndor/gspa/proteomes/bench_gtdb30/plm_catalysts/ec_panel.index.tsv

N=$(tail -n +2 $CASES | wc -l)
SHARDS=${SLURM_ARRAY_TASK_COUNT:-5}
ID=${SLURM_ARRAY_TASK_ID:-0}
CHUNK=$(( (N + SHARDS - 1) / SHARDS ))
START=$(( ID * CHUNK ))
END=$(( START + CHUNK ))
if [[ $END -gt $N ]]; then END=$N; fi

OUT=$LRO/runs_rlgc_plm/track_A_strict/shard_${ID}
mkdir -p $OUT $LRO/logs

echo "== PLM strict Track A shard=$ID start=$START end=$END =="
date
python3 $LRO/run_ablation_rlgc_plm.py \
    --cases $CASES --root $ROOT --tag mg1655 \
    --jar $JAR --go-owl $GO_OWL \
    --pathways $PATHWAYS --ec2go $EC2GO \
    --reaction-graph $GAPS_DATA/seed_reactions.tsv \
    --diffusion-mets $GAPS_DATA/diffusion_mets.tsv \
    --ec-aliases $GAPS_DATA/seed_Enzyme_Class_Reactions_Aliases_unique.tsv \
    --genome-layout $LAYOUT \
    --model $MODEL \
    --centroids-dir $CENTROIDS \
    --target-plm-npy $PLM_NPY \
    --target-plm-index $PLM_IDX \
    --panel-npy $PANEL_NPY \
    --panel-index $PANEL_IDX \
    --catalysts-tsv $CATALYSTS \
    --track A_strict \
    --track-a-threshold 0.5 \
    --mode $MODE --out-dir $OUT \
    --start $START --end $END \
    --java $JAVA
echo "DONE"; date
