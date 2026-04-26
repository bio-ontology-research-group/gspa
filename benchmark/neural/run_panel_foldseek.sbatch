#!/usr/bin/env bash
#SBATCH --job-name=gspa-panel-foldseek
#SBATCH --partition=debug
#SBATCH --gres=gpu:1
#SBATCH --exclude=node001,node002,node003,node007
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 04:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/panel-foldseek-%A_%a.out
#SBATCH -e /data/hohndor/gspa-neural/logs/panel-foldseek-%A_%a.err

# FoldSeek homology-transfer predictor:
#   query FAA → ProstT5 → 3Di → easy-search vs afdb_swissprot → top-K hits
#   → UniProt hit AC → look up GO+EC from swissprot_go_ec.tsv → emit annotations
#
# Score used: bits rescaled to [0, 1] via min-max within the genome.
# The centroid-mode column in the GSPA predictor maps to this equivalent
# view when the DB is SwissProt-wide (every hit has its own label set,
# not a per-class medoid).

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix
export PYTHONUNBUFFERED=1

PANEL=/data/hohndor/gspa-neural/panel
NEURAL=/data/hohndor/gspa-neural/benchmark/neural
MANIFEST=$PANEL/panel_manifest.tsv
AFDB=/data/hohndor/gspa-neural/models/tmp_sp/latest/afdb_swissprot
PROSTT5=/data/hohndor/gspa-neural/models/ProstT5
META=/data/hohndor/gspa-neural/sprot/swissprot_go_ec.tsv

ROW=$(( SLURM_ARRAY_TASK_ID + 1 ))
TAG=$(awk -F'\t' -v r=$ROW 'NR==r{print $1}' $MANIFEST)
FASTA=$PANEL/fastas/${TAG}_protein.faa
OUT=$PANEL/preds/$TAG
TMP=$(mktemp -d -p /tmp)
mkdir -p $OUT

[[ -s $FASTA ]] || { echo "FASTA missing: $FASTA"; exit 2; }
[[ -s $META ]]  || { echo "swissprot meta missing: $META"; exit 2; }
echo "=== $TAG ==="

ALN=$OUT/${TAG}.foldseek_aln.tsv
foldseek easy-search \
    $FASTA \
    $AFDB \
    $ALN \
    $TMP \
    --prostt5-model $PROSTT5 \
    --format-mode 4 \
    --format-output "query,target,fident,alnlen,evalue,bits" \
    --max-seqs 5 \
    --threads 8

rm -rf $TMP

# Post-process: for each hit, emit GO+EC annotations from SwissProt meta.
python - <<'EOF'
import os
from pathlib import Path
from collections import defaultdict

OUT = os.environ["OUT"]
TAG = os.environ["TAG"]
META = os.environ["META"]
ALN = Path(f"{OUT}/{TAG}.foldseek_aln.tsv")

# Load UniProt → (go_terms, ec_numbers) once.
lookup = {}
with open(META) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 3: continue
        acc = p[0]
        gos = [g for g in p[1].split(";") if g.startswith("GO:")]
        ecs = ["EC:" + e for e in p[2].split(";") if e and "-" not in e]
        lookup[acc] = (gos, ecs)
print(f"loaded meta: {len(lookup)} accs", flush=True)

# Parse alignment, aggregate (query, term) → max bitscore
agg = defaultdict(lambda: defaultdict(float))  # pid → {term: score}
ann_type_of = {}
n_hits = 0
with ALN.open() as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    cols = {c: i for i, c in enumerate(hdr)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < len(hdr): continue
        q = f[cols["query"]]
        t_raw = f[cols["target"]]
        # target may be "AF-P12345-F1-model_v4" — extract the ACC
        acc = None
        if t_raw.startswith("AF-"):
            acc = t_raw.split("-")[1]
        else:
            # pipe-delimited fallback
            if "|" in t_raw:
                acc = t_raw.split("|")[1] if len(t_raw.split("|")) > 1 else t_raw
            else:
                acc = t_raw
        try:
            bits = float(f[cols["bits"]])
        except ValueError:
            continue
        gos, ecs = lookup.get(acc, ([], []))
        for term in gos:
            key = (q, term)
            if bits > agg[q].get(term, 0.0):
                agg[q][term] = bits
                ann_type_of[(q, term)] = "GO"
        for term in ecs:
            key = (q, term)
            if bits > agg[q].get(term, 0.0):
                agg[q][term] = bits
                ann_type_of[(q, term)] = "EC"
        n_hits += 1
print(f"hits: {n_hits}  proteins with hits: {len(agg)}", flush=True)

# Normalize bits to [0, 1] per genome (simple min-max)
all_scores = [s for m in agg.values() for s in m.values()]
if not all_scores:
    print("no hits; empty output", flush=True)
    with open(f"{OUT}/{TAG}.foldseek.tsv", "w") as fh:
        fh.write("protein_id\tterm\tscore\tannotation_type\n")
    exit(0)
mn, mx = min(all_scores), max(all_scores)
rng = mx - mn if mx > mn else 1.0

with open(f"{OUT}/{TAG}.foldseek.tsv", "w") as fh:
    fh.write("protein_id\tterm\tscore\tannotation_type\n")
    for pid, per in sorted(agg.items()):
        for term, sc in per.items():
            norm = (sc - mn) / rng
            fh.write(f"{pid}\t{term}\t{norm:.4f}\t{ann_type_of[(pid, term)]}\n")
print(f"wrote {OUT}/{TAG}.foldseek.tsv", flush=True)
EOF

wc -l $OUT/${TAG}.foldseek.tsv
echo DONE $TAG; date
