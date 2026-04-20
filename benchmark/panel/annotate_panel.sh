#!/usr/bin/env bash
#SBATCH --job-name=panel-anno
#SBATCH --partition=debug
#SBATCH --exclude=node003
#SBATCH --output=/data/hohndor/gspa/proteomes/culture_panel/phase2/logs/anno-%A_%a.out
#SBATCH --error=/data/hohndor/gspa/proteomes/culture_panel/phase2/logs/anno-%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
set -u

# Usage: sbatch --array=0-N%K annotate_panel.sh <manifest_tsv>
# manifest format: tag<TAB>fasta_path
#
# Adapted from benchmark/cultures/annotate_culture.sh to write into
# a per-genome subdir under culture_panel/phase2/<tag>/.

MANIFEST=${1:-/data/hohndor/gspa/proteomes/culture_panel/phase2_manifest.tsv}
ROOT=/data/hohndor/gspa/proteomes/culture_panel/phase2
TAGS=($(tail -n +2 $MANIFEST | cut -f1))
FASTAS=($(tail -n +2 $MANIFEST | cut -f2))
tag=${TAGS[$SLURM_ARRAY_TASK_ID]}
fasta=${FASTAS[$SLURM_ARRAY_TASK_ID]}

WD=$ROOT/$tag
mkdir -p $WD/{prodigal,preds,claims,layout,gapsmith,integrated}
cd $WD

DMND_DB=/data/hohndor/gspa/benchmark/benchmark_data/reference_db.dmnd
PFAM=/storage/software/databases/gtdbtk/markers/pfam/Pfam-A.hmm
GAPSMITH=/data/hohndor/gspa/bin/gapsmith/gapsmith-linux
GS_DATA=/data/hohndor/gspa/bin/gapsmith/data_merged
JAR=/data/hohndor/gspa/bin/gspa-phase12.jar
JAVA=/storage/miniforge3/envs/metagenomics/bin/java
BP=/data/hohndor/gspa/benchmark-py
GO_OWL=/data/hohndor/gspa/reference/go.owl
EC2GO=/data/hohndor/gspa/reference/ec2go.txt
PATHWAYS=/data/hohndor/gspa/reference/kegg_pathways.tsv
GOA=/data/hohndor/gspa/benchmark/goa_uniprot_all.gaf.gz

echo "=== [$tag] $(date) annotate start ==="
echo "fasta: $fasta"
echo "wd:    $WD"

# 1. prodigal gene calling (meta mode is safe for MAGs + isolates)
if [[ ! -s prodigal/$tag.faa ]]; then
    echo "[1/6] prodigal"
    prodigal -i $fasta -a prodigal/$tag.faa -f gff \
        -o prodigal/$tag.gff -p meta -q > prodigal/$tag.log 2>&1
fi
N_CDS=$(grep -c '^>' prodigal/$tag.faa)
echo "  $N_CDS proteins called"

# 2. layout.tsv from GFF
if [[ ! -s layout/$tag\_layout.tsv ]]; then
    echo "[2/6] layout"
    python3 - <<PY > layout/$tag\_layout.tsv
import re
print("protein_id\tcontig\tstart\tend\tstrand")
with open("prodigal/$tag.gff") as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.rstrip().split("\t")
        if len(p) < 9 or p[2] != "CDS": continue
        m = re.search(r"ID=([^;]+)", p[8])
        if not m: continue
        print(f"{m.group(1)}\t{p[0]}\t{p[3]}\t{p[4]}\t{p[6]}")
PY
fi

# 3. DIAMOND + Pfam
if [[ ! -s preds/diamond_results.tsv ]]; then
    echo "[3a/6] DIAMOND"
    diamond blastp -q prodigal/$tag.faa -d $DMND_DB \
        -o preds/diamond_results.tsv \
        --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
        --max-target-seqs 50 --evalue 1e-5 --threads 8 --quiet
fi
if [[ ! -s preds/pfam_results.domtbl ]]; then
    echo "[3b/6] Pfam hmmsearch"
    hmmsearch --cpu 8 --domtblout preds/pfam_results.domtbl \
        --cut_ga $PFAM prodigal/$tag.faa > preds/pfam.log 2>&1
fi
touch preds/mmseqs_results.tsv preds/eggnog_results.tsv

# 4. claims
if [[ ! -s claims/$tag\_dp_claims.jsonl ]]; then
    echo "[4/6] claims"
    python3 $BP/02b_parse_predictors_to_claims.py \
        --results-dir preds --goa $GOA \
        --output claims/$tag\_dp_claims.jsonl > claims/parse.log 2>&1
fi
wc -l claims/$tag\_dp_claims.jsonl

# 5. gapsmith find -p all
if [[ ! -s gapsmith/$tag-all-Reactions.tbl ]]; then
    echo "[5/6] gapsmith find"
    $GAPSMITH --data-dir $GS_DATA find -p all -t Bacteria \
        -o gapsmith -u all -A diamond -K 8 \
        prodigal/$tag.faa > gapsmith/$tag.log 2>&1 || true
fi

# 6. integrate (whole-genome posteriors)
if [[ ! -s integrated/$tag\_integrated.tsv ]]; then
    echo "[6/6] integrate"
    $JAVA -Xmx10g -jar $JAR integrate \
        --claims claims/$tag\_dp_claims.jsonl \
        --out integrated/$tag\_integrated.tsv \
        --go-owl $GO_OWL --lite \
        --essential-profile bacteria \
        --pathways $PATHWAYS --ec2go $EC2GO
fi

echo "=== [$tag] $(date) DONE ==="
