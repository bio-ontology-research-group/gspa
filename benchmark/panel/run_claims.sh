#!/bin/bash
#SBATCH --job-name=panel-claims
#SBATCH --partition=debug
#SBATCH --output=/data/hohndor/gspa/proteomes/bench_gtdb30/claims-%A_%a.out
#SBATCH --error=/data/hohndor/gspa/proteomes/bench_gtdb30/claims-%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-29

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics

ROOT=/data/hohndor/gspa/proteomes/bench_gtdb30
BP=/data/hohndor/gspa/benchmark-py
cd $ROOT
mkdir -p claims preds

TAGS=($(tail -n +2 panel_manifest.tsv | cut -f1))
tag=${TAGS[$SLURM_ARRAY_TASK_ID]}
faa=$ROOT/proteomes/${tag}.faa
[[ -s $faa ]] || { echo "missing $faa; skip"; exit 0; }

echo "=== $tag: claims ==="
date

workdir=preds/${tag}_work
mkdir -p $workdir

dm_out=$workdir/diamond_results.tsv
if [[ ! -s $dm_out ]]; then
    diamond blastp --query $faa --db $ROOT/reference_panel --threads 8 \
        --outfmt 6 --evalue 1e-5 --max-target-seqs 5 --quiet --out $dm_out
fi
echo "  diamond: $(wc -l < $dm_out) rows"

hmm_out=$workdir/pfam_results.domtbl
if [[ ! -s $hmm_out ]]; then
    hmmsearch --domtblout $hmm_out --noali --cpu 8 --cut_ga \
        /storage/software/databases/gtdbtk/markers/pfam/Pfam-A.hmm $faa > /dev/null
fi
echo "  pfam:    $(wc -l < $hmm_out) rows"

# Empty placeholders so the parser finds all expected files
touch $workdir/mmseqs2_results.tsv
touch $workdir/eggnog_results.emapper.annotations

python3 $BP/02b_parse_predictors_to_claims.py \
    --results-dir $workdir \
    --goa /data/hohndor/gspa/benchmark/goa_uniprot_all.gaf.gz \
    --pfam2go /data/hohndor/gspa/reference/pfam2go.txt \
    --output claims/${tag}_dp_claims.jsonl
echo "  claims:  $(wc -l < claims/${tag}_dp_claims.jsonl)"

date
echo DONE_$tag
