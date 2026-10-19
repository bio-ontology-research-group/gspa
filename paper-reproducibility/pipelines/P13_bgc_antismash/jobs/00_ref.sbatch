#!/bin/bash
#SBATCH --job-name=p13-ref
#SBATCH --partition=batch
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
# Build references: genome list, experimental-GO Swiss-Prot DIAMOND DB, class-term check, checksums.
set -euo pipefail
OUT=${OUT:-/ibex/scratch/projects/c2014/hohndor/bgc_repro_2026-10-19}
source $OUT/code/jobs/common.sh
module load python/3.11.0 diamond/2.1.16 >/dev/null 2>&1
mkdir -p $REF $WORK $INTEG $RES $LOGS
cd $REF
# genome list = antiSMASH output dirs joined with the METEOR taxonomy table
ls $ANTISMASH | sort > antismash_dirs.txt
awk -F'\t' 'NR==FNR{d[$1]=1;next} FNR>1 && ($1 in d){print $1"\t"$3"\t"$2}' antismash_dirs.txt $TAXTAB \
  | sort > genomes_body.tsv
( printf 'genome\ttaxid\torganism\n'; cat genomes_body.tsv ) > $OUT/genomes.tsv
rm genomes_body.tsv
echo "genomes: $(($(wc -l < $OUT/genomes.tsv)-1)) (antiSMASH dirs: $(wc -l < antismash_dirs.txt))"
cp $GO_OBO_SRC go.obo
for f in go.owl pfam2go.txt kegg_pathways.tsv kegg_modules.tsv ec2go.txt; do cp $ATB_REF/$f .; done
python3 $BGC verify-terms --obo go.obo | tee class_terms_verified.tsv
python3 $BGC build-ref --sprot-dat $SPROT_DAT --out-fasta sprot_exp.fasta --out-go sprot_exp_go.tsv | tee build_ref.log
diamond makedb --in sprot_exp.fasta -d sprot_exp --threads 8 --quiet
diamond version >> build_ref.log
python3 - <<'PY' >> build_ref.log
import json,glob
d=json.load(open(sorted(glob.glob("/ibex/scratch/projects/c2014/kexin/funcarve/paperA_2026/results/antismash_gc/*/*.json"))[0]))
print("antismash_version", d.get("version"), "schema", d.get("schema"), "taxon", d.get("taxon"))
PY
sha256sum $JAR $SPROT_DAT $PFAM_HMM go.obo go.owl pfam2go.txt kegg_pathways.tsv kegg_modules.tsv ec2go.txt \
  sprot_exp.fasta sprot_exp_go.tsv sprot_exp.dmnd > sha256.txt
grep -c '^NAME' $PFAM_HMM | sed 's/^/pfam_hmm_models /' >> build_ref.log
cat build_ref.log
