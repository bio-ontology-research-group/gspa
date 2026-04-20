#!/usr/bin/env bash
#SBATCH --job-name=panel-gtdbtk
#SBATCH --partition=debug
#SBATCH --exclude=node003
#SBATCH -c 16
#SBATCH --mem=100G
#SBATCH -t 24:00:00
#SBATCH -o /data/hohndor/gspa/proteomes/culture_panel/phase1/logs/gtdbtk-%j.out
#SBATCH -e /data/hohndor/gspa/proteomes/culture_panel/phase1/logs/gtdbtk-%j.err

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
set -u

ROOT=/data/hohndor/gspa/proteomes/culture_panel/phase1
STAGE=$ROOT/staged
OUT=$ROOT/gtdbtk
export GTDBTK_DATA_PATH=/storage/software/databases/gtdbtk

mkdir -p $OUT

echo "gtdbtk $(gtdbtk --version 2>&1 | head -1)"
echo "genomes:        $(ls $STAGE | wc -l)"
echo "GTDBTK_DATA:    $GTDBTK_DATA_PATH"
date

# classify_wf includes ANI screen + placement; --skip_ani_screen is
# faster but less accurate for genomes near GTDB reps.  Keep default.
gtdbtk classify_wf \
    --genome_dir $STAGE \
    --out_dir $OUT \
    --extension fna \
    --cpus 16 \
    --pplacer_cpus 2 \
    --skip_ani_screen \
    --force 2>&1 | tail -50

echo
echo "=== summary ==="
if [[ -f $OUT/gtdbtk.bac120.summary.tsv ]]; then
    wc -l $OUT/gtdbtk.bac120.summary.tsv
    echo "phyla:"
    awk -F'\t' 'NR>1 { split($2, a, ";"); for (i in a) if (a[i] ~ /^p__/) print a[i] }' \
        $OUT/gtdbtk.bac120.summary.tsv | sort | uniq -c | sort -rn
fi
if [[ -f $OUT/gtdbtk.ar53.summary.tsv ]]; then
    wc -l $OUT/gtdbtk.ar53.summary.tsv
fi

echo DONE; date
