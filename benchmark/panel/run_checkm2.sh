#!/usr/bin/env bash
#SBATCH --job-name=panel-checkm2
#SBATCH --partition=debug
#SBATCH --exclude=node003
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o /data/hohndor/gspa/proteomes/culture_panel/phase1/logs/checkm2-%j.out
#SBATCH -e /data/hohndor/gspa/proteomes/culture_panel/phase1/logs/checkm2-%j.err

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate checkm2-v3
set -u

ROOT=/data/hohndor/gspa/proteomes/culture_panel/phase1
STAGE=$ROOT/staged
OUT=$ROOT/checkm2
DB=/data/databases/checkm2/uniref100.KO.1.dmnd

mkdir -p $OUT

echo "checkm2 $(checkm2 --version 2>&1 | head -1)"
echo "genomes: $(ls $STAGE | wc -l)"
echo "db:      $DB"
date

checkm2 predict \
    --input $STAGE \
    --output-directory $OUT \
    --database_path $DB \
    --threads 16 \
    --extension .fna \
    --force 2>&1 | tail -40

echo
echo "=== quality_report.tsv summary ==="
wc -l $OUT/quality_report.tsv
head -1 $OUT/quality_report.tsv
awk -F'\t' 'NR>1 { if ($2 >= 90 && $3 <= 5) h++;
                  else if ($2 >= 70 && $3 <= 10) m++;
                  else if ($2 >= 50 && $3 <= 15) l++;
                  else e++
             }
    END { printf "high=%d medium=%d low=%d excluded=%d\n", h, m, l, e }
' $OUT/quality_report.tsv

echo DONE; date
