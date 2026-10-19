#!/bin/bash
#SBATCH --job-name=p14-select
#SBATCH --partition=debug
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
# Snapshot references into $REF, record checksums/versions, choose the MAG set.
set -euo pipefail
OUT=${OUT:-/ibex/scratch/projects/c2014/hohndor/mags_repro_2026-10-19}
source $OUT/code/jobs/common.sh
module load python/3.11.0 diamond/2.1.16 hmmer/3.3.2 >/dev/null 2>&1
mkdir -p $REF $WORK $RES $LOGS
cd $REF
for f in sprot_exp.dmnd sprot_exp_go.tsv go.obo go.owl; do [[ -s $f ]] || cp $P13_REF/$f .; done
for f in go-dag.tsv pfam2go.txt ec2go.txt kegg_pathways.tsv kegg_modules.tsv; do [[ -s $f ]] || cp $ATB_REF/$f .; done
for f in seed_reactions_corrected.tsv diffusion_mets.tsv seed_Enzyme_Class_Reactions_Aliases_unique.tsv; do [[ -s $f ]] || cp $GAPSEQ_DAT/$f .; done
mkdir -p essential
for k in bacteria archaea; do unzip -p $JAR essential-functions/$k.tsv > essential/$k.tsv; done
{
  echo "date $(date -Is) host $(hostname)"
  echo "gspa_jar $(basename $JAR) $(unzip -p $JAR META-INF/MANIFEST.MF 2>/dev/null | grep -i 'version' | tr -d '\r' | tr '\n' ' ')" || true
  $JAVA -version 2>&1 | sed -n 1p || true
  diamond version
  hmmsearch -h | sed -n 2p
  python3 --version
  grep -m1 'gapseq version' $GAPSEQ_OUT/10Dr2__SemiBin_12/10Dr2__SemiBin_12-all-Pathways.tbl | sed 's/^.*# //' || true
  echo "pfam_hmm_models $(grep -c '^NAME' $PFAM_HMM)"
  grep -m1 '^data-version' go.obo
  echo "sprot_exp_entries $(grep -c '^>' $P13_REF/sprot_exp.fasta)"
} > tool_versions.txt
cat tool_versions.txt
sha256sum $JAR $PFAM_HMM $HQ_LIST $CHECKM2 $GTDB_BAC $GTDB_AR $CODE/*.py \
  sprot_exp.dmnd sprot_exp_go.tsv go.obo go.owl go-dag.tsv pfam2go.txt ec2go.txt kegg_pathways.tsv kegg_modules.tsv \
  seed_reactions_corrected.tsv diffusion_mets.tsv seed_Enzyme_Class_Reactions_Aliases_unique.tsv essential/*.tsv > sha256.txt
cat sha256.txt
python3 $PY select --hq-list $HQ_LIST --checkm2 $CHECKM2 --gtdb $GTDB_BAC $GTDB_AR --gapseq-out $GAPSEQ_OUT \
  --n 100 --full-set-max 150 --min-phylum 10 --base 4 --seed 20261019 \
  --eligible-out $OUT/eligible_mags.tsv --out $OUT/mag_list.tsv | tee $OUT/select.log
