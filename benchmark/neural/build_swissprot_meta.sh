#!/usr/bin/env bash
#SBATCH --job-name=gspa-sprot-meta
#SBATCH --partition=debug
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/sprot-meta-%j.out
#SBATCH -e /data/hohndor/gspa-neural/logs/sprot-meta-%j.err

# Produce SwissProt-wide metadata artifacts:
#  - sprot_accs.txt           (UniProt ACs, one per line)
#  - sprot_refseq_map.tsv     RefSeq_id -> SwissProt_acc (from collab, filtered)
#  - sprot_go.tsv             UniProt_acc -> {aspect, GO:...} rows (GOA scan)
#  - sprot_ec.tsv             UniProt_acc -> EC:... rows (uniprot_sprot.dat scan)
#  - swissprot_go_ec.tsv      (accession, go_terms, ec_numbers) — centroid input

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix
export PYTHONUNBUFFERED=1

WORK=/data/hohndor/gspa-neural/sprot
REF=/data/hohndor/gspa/reference
BENCH=/data/hohndor/gspa/benchmark
NEURAL=/data/hohndor/gspa-neural/benchmark/neural
mkdir -p $WORK

# 1) SwissProt accession list (UniProt entries currently reviewed)
grep "^>" $BENCH/uniprot_sprot.fasta \
    | awk -F'|' '{print $2}' \
    | sort -u > $WORK/sprot_accs.txt
echo "sprot_accs: $(wc -l < $WORK/sprot_accs.txt)"

# 2) SwissProt-filtered RefSeq → UniProt map (for better panel truth)
zcat $REF/gene_refseq_uniprotkb_collab.gz \
    | awk -F'\t' 'NR==FNR{a[$0]=1; next} FNR==1 || $2 in a' \
        $WORK/sprot_accs.txt - \
    > $WORK/sprot_refseq_map.tsv
echo "sprot_refseq_map: $(wc -l < $WORK/sprot_refseq_map.tsv)"

# 3) GOA scan → sprot_go.tsv (UniProt→GO; aspect as MF/BP/CC)
python /data/hohndor/gspa-neural/benchmark/panel/extract_goa_dual.py \
    --goa $BENCH/goa_uniprot_all.gaf.gz \
    --accessions sprot:$WORK/sprot_accs.txt \
    --out-dir $WORK
mv $WORK/sprot_truth_all.tsv $WORK/sprot_go.tsv
echo "sprot_go: $(wc -l < $WORK/sprot_go.tsv)"

# 4) EC scan → sprot_ec.tsv
python $NEURAL/extract_ec_truth.py \
    --dat $REF/uniprot_sprot.dat.gz \
    --accessions $WORK/sprot_accs.txt \
    --out $WORK/sprot_ec.tsv
echo "sprot_ec: $(wc -l < $WORK/sprot_ec.tsv)"

# 5) Merge GO+EC per accession → swissprot_go_ec.tsv
python - <<'EOF'
from collections import defaultdict
WORK = "/data/hohndor/gspa-neural/sprot"
gos = defaultdict(set)
ecs = defaultdict(set)
with open(f"{WORK}/sprot_go.tsv") as fh:
    next(fh, None)
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3 and p[2].startswith("GO:"):
            gos[p[0]].add(p[2])
with open(f"{WORK}/sprot_ec.tsv") as fh:
    next(fh, None)
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3 and p[2].startswith("EC:"):
            ecs[p[0]].add(p[2][3:])  # strip "EC:" for centroid TSV format
accs = set(gos) | set(ecs)
with open(f"{WORK}/swissprot_go_ec.tsv", "w") as fh:
    fh.write("accession\tgo_terms\tec_numbers\n")
    for acc in sorted(accs):
        fh.write(f"{acc}\t{';'.join(sorted(gos.get(acc, ())))}\t{';'.join(sorted(ecs.get(acc, ())))}\n")
print(f"wrote swissprot_go_ec.tsv: {len(accs)} accessions")
EOF

echo DONE; date
