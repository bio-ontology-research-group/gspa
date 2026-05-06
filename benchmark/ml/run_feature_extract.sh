#!/bin/bash
#SBATCH -A kaust
#SBATCH -p debug
#SBATCH --exclude=node003
#SBATCH -c 4
#SBATCH --mem=12G
#SBATCH -t 24:00:00
#SBATCH -J m3_extract
#SBATCH -o /data/hohndor/gspa/proteomes/bench_gtdb30/m3_extract-%A_%a.out
#SBATCH -e /data/hohndor/gspa/proteomes/bench_gtdb30/m3_extract-%A_%a.err

# Usage: sbatch --array=0-28 run_feature_extract.sh  (29 panel genomes excluding mg1655)

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics

ROOT=/data/hohndor/gspa/proteomes/bench_gtdb30
JAR=/data/hohndor/gspa/bin/gspa-phase12.jar
JAVA=/storage/miniforge3/envs/metagenomics/bin/java
GO_OWL=/data/hohndor/gspa/reference/go.owl
PATHWAYS=/data/hohndor/gspa/reference/kegg_pathways.tsv
EC2GO=/data/hohndor/gspa/reference/ec2go.txt
GAPS_DATA=/data/hohndor/gspa/bin/gapsmith/data_merged

# Panel tags excluding mg1655 (held out for M3 evaluation)
TAGS=($(tail -n +2 $ROOT/panel_manifest.tsv | awk '$1 != "mg1655"' | cut -f1))
tag=${TAGS[$SLURM_ARRAY_TASK_ID]}

CLAIMS=$ROOT/claims/${tag}_dp_claims.jsonl
LAYOUT=$ROOT/layout/${tag}_layout.tsv
TRUTH=$ROOT/truth/${tag}_truth_all.tsv
REACTIONS_TBL=$ROOT/gapsmith/${tag}/${tag}-all-Reactions.tbl
OUTDIR=$ROOT/m3_features/${tag}
mkdir -p $OUTDIR

[[ -s $CLAIMS && -s $LAYOUT && -s $TRUTH && -s $REACTIONS_TBL ]] || { echo "missing data for $tag"; exit 0; }

echo "=== $tag: M3 feature extraction ==="
date

# Build per-genome test set: good_blast rows with EC + GO + truth match
python3 - <<PY
import csv, re
truth = set()
with open("${TRUTH}") as f:
    next(f)
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) >= 3: truth.add((p[0], p[2]))

ec2go = {}
with open("${EC2GO}") as f:
    for line in f:
        if line.startswith('!') or not line.strip(): continue
        m = re.match(r'^EC:(\S+)\s*>\s*GO:[^;]+;\s*(GO:\d+)', line)
        if m: ec2go[m.group(1)] = m.group(2)

cases = []
with open("${REACTIONS_TBL}") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row['status'] != 'good_blast': continue
        ec = row['ec'].strip()
        prot = row['stitle'].strip()
        if not ec or not prot: continue
        if '-' in ec: continue
        go = ec2go.get(ec)
        if not go: continue
        if (prot, go) not in truth: continue
        cases.append((row['pathway'], row['rxn'], ec, go, prot))

# Cap at 100 cases per genome
import random
random.seed(42)
if len(cases) > 100:
    cases = random.sample(cases, 100)

with open("${OUTDIR}/cases.tsv", 'w') as out:
    out.write('pathway\trxn\tec\tgo\tprotein\n')
    for c in cases:
        out.write('\t'.join(c) + '\n')
print(f"[${tag}] {len(cases)} cases")
PY

# For each case: ablate → integrate with --features-out → tag rows with qid+label
N=0
> $OUTDIR/training.tsv
WROTE_HEADER=0
while IFS=$'\t' read pathway rxn ec go prot; do
    [[ "$pathway" == "pathway" ]] && continue
    N=$((N+1))
    WD=$OUTDIR/case_$(printf %04d $N)
    mkdir -p $WD
    abl=$WD/claims.jsonl
    gap=$WD/gap.jsonl
    feats=$WD/features.tsv
    python3 -c "
import json
with open('${CLAIMS}') as fin, open('${abl}', 'w') as fout:
    for line in fin:
        try: obj = json.loads(line)
        except: fout.write(line); continue
        if obj.get('protein_id') == '${prot}': continue
        fout.write(line)
"
    python3 -c "
import json
print(json.dumps({'pathway_id':'${pathway}','reaction_id':'${rxn}','ec_number':'${ec}','go_term':'${go}','gapseq_guessed':False}), end='')
" > $gap
    echo "" >> $gap
    $JAVA -Xmx6g -jar $JAR integrate \
        --claims $abl --out $WD/integrated.tsv \
        --suggestions-out $WD/sugg.tsv --features-out $feats \
        --go-owl $GO_OWL --lite \
        --essential-profile bacteria \
        --pathways $PATHWAYS --ec2go $EC2GO \
        --gaps $gap \
        --reaction-graph $GAPS_DATA/seed_reactions.tsv \
        --diffusion-mets $GAPS_DATA/diffusion_mets.tsv \
        --reaction-ec-aliases $GAPS_DATA/seed_Enzyme_Class_Reactions_Aliases_unique.tsv \
        --genome-layout $LAYOUT \
        --rlc-suggester \
        --enable-priors essentiality,coherence,gap_filling,genomic_context > $WD/integrate.log 2>&1 || true
    if [[ -s $feats ]]; then
        if [[ $WROTE_HEADER == 0 ]]; then
            head -1 $feats | awk -v OFS='\t' '{print "qid","label",$0}' > $OUTDIR/training.tsv
            WROTE_HEADER=1
        fi
        awk -v qid="${tag}_$N" -v prot="$prot" 'NR>1 { lbl = ($1 == prot) ? 1 : 0; print qid"\t"lbl"\t"$0 }' $feats >> $OUTDIR/training.tsv
    fi
    rm -rf $WD
done < $OUTDIR/cases.tsv

echo "[${tag}] wrote $(tail -n +2 $OUTDIR/training.tsv | wc -l) rows"
date
