#!/usr/bin/env bash
#SBATCH --job-name=panel-skani
#SBATCH --partition=debug
#SBATCH --exclude=node003
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 06:00:00
#SBATCH -o /data/hohndor/gspa/proteomes/culture_panel/phase1/logs/skani-%j.out
#SBATCH -e /data/hohndor/gspa/proteomes/culture_panel/phase1/logs/skani-%j.err

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
set -u

ROOT=/data/hohndor/gspa/proteomes/culture_panel/phase1
STAGE=$ROOT/staged
OUT=$ROOT/skani

mkdir -p $OUT

echo "skani $(skani -V 2>&1)"
echo "genomes: $(ls $STAGE | wc -l)"
date

# All-vs-all ANI on staged FASTAs
ls $STAGE/*.fna > $OUT/input_list.txt

# skani triangle (sparse all-vs-all ANI)
skani triangle \
    -l $OUT/input_list.txt \
    -o $OUT/ani_matrix.sparse.tsv \
    --sparse -t 16 2>&1 | tail -20

wc -l $OUT/ani_matrix.sparse.tsv

# Cluster at 95% ANI (species-level). Using skani dereplicate-style
# greedy clustering directly here; skani has a built-in "cluster":
skani dist \
    --ql $OUT/input_list.txt \
    --rl $OUT/input_list.txt \
    -o $OUT/ani_pairs.tsv \
    -t 16 --min-af 50 2>&1 | tail -5

# Greedy ANI-95 dereplication with size bias toward large contigs
# (prefer higher-N50 members as cluster reps).
python3 - <<'PY'
import csv
from pathlib import Path
import subprocess

ani_file = Path("/data/hohndor/gspa/proteomes/culture_panel/phase1/skani/ani_pairs.tsv")
out = Path("/data/hohndor/gspa/proteomes/culture_panel/phase1/skani/clusters.tsv")
input_list = Path("/data/hohndor/gspa/proteomes/culture_panel/phase1/skani/input_list.txt")

genomes = [line.strip() for line in open(input_list) if line.strip()]

# ANI >= 95, AF >= 50
edges = {}  # g1 -> set of g2
with open(ani_file) as f:
    rdr = csv.DictReader(f, delimiter="\t")
    for r in rdr:
        try:
            ani = float(r["ANI"])
            af_r = float(r.get("Align_fraction_ref", "0"))
            af_q = float(r.get("Align_fraction_query", "0"))
        except (KeyError, ValueError):
            continue
        if ani < 95 or max(af_r, af_q) < 50:
            continue
        a = r["Ref_file"]; b = r["Query_file"]
        if a == b:
            continue
        edges.setdefault(a, set()).add(b)
        edges.setdefault(b, set()).add(a)

# Prefer genomes with fewer contigs as rep (better assembly)
def n_contigs(path):
    try:
        return sum(1 for L in open(path) if L.startswith(">"))
    except OSError:
        return 10**9


def size_or_zero(path):
    try:
        return Path(path).stat().st_size
    except (OSError, FileNotFoundError):
        return 0


ranks = sorted(genomes, key=lambda p: (n_contigs(p), -size_or_zero(p)))

seen = set()
cluster_of = {}
for g in ranks:
    if g in seen:
        continue
    # new cluster rep
    cluster_of[g] = g
    seen.add(g)
    # anyone connected (incl. multi-hop) gets this rep
    stack = list(edges.get(g, set()))
    while stack:
        h = stack.pop()
        if h in seen:
            continue
        cluster_of[h] = g
        seen.add(h)
        stack.extend(edges.get(h, set()) - seen)

with open(out, "w") as f:
    f.write("member\tcluster_rep\n")
    for m, r in cluster_of.items():
        f.write(f"{m}\t{r}\n")

reps = set(cluster_of.values())
print(f"[skani-derep] {len(genomes)} genomes -> {len(reps)} clusters "
      f"at ANI-95")
PY

echo DONE; date
