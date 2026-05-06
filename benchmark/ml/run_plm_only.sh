#!/usr/bin/env bash
#SBATCH --job-name=m4-plmonly
#SBATCH --partition=debug
#SBATCH --exclude=node003
#SBATCH -c 4
#SBATCH --mem=12G
#SBATCH -t 01:00:00
#SBATCH -o /data/hohndor/gspa/proteomes/bench_gtdb30/plmonly-%j.out
#SBATCH -e /data/hohndor/gspa/proteomes/bench_gtdb30/plmonly-%j.err

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
set -u

ROOT=/data/hohndor/gspa/proteomes/bench_gtdb30
cd $ROOT/ml

# Strip train/valid to only {qid, label, string-cols, plm_cos_centroid_EC, plm_has_emb}
# We keep the string-cols so the training script's header position indices match.
python3 - <<'PY'
IN_COLS_KEEP = {
    'qid', 'label', 'protein_id', 'reaction_id', 'pathway_id', 'ec', 'go_term',
    'plm_cos_centroid_EC', 'plm_has_emb',
}
for src, dst in [('train_m4.tsv', 'train_plmonly.tsv'),
                 ('valid_m4.tsv', 'valid_plmonly.tsv')]:
    with open(src) as fin, open(dst, 'w') as fout:
        header = fin.readline().rstrip('\n').split('\t')
        keep_idx = [i for i, c in enumerate(header) if c in IN_COLS_KEEP]
        fout.write('\t'.join([header[i] for i in keep_idx]) + '\n')
        for line in fin:
            parts = line.rstrip('\n').split('\t')
            fout.write('\t'.join([parts[i] for i in keep_idx]) + '\n')
PY

python3 /data/hohndor/gspa/benchmark-py/train_lambdamart.py \
    --train train_plmonly.tsv \
    --valid valid_plmonly.tsv \
    --out model_plmonly.txt \
    --iters 500 --early-stop 30 \
    --feature-importance importance_plmonly.json

echo "DONE"
date
