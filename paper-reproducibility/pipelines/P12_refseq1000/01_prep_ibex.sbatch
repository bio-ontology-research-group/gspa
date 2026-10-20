#!/bin/bash
#SBATCH -J p12_prep
#SBATCH -p batch
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -o %x_%j.log
# P12 step 01 (IBEX): build DIAMOND db, press Pfam-A, download RefSeq GFF3s,
# build operons and the organism-taxon map.  Idempotent.
set -euo pipefail
R=${R:-/ibex/scratch/projects/c2014/hohndor/refseq1000_repro_2026-10-19}
C=/encrypted/e3008/annot_evaluation
P=$R/scripts
type module >/dev/null 2>&1 || source /etc/profile.d/modules.sh
module load diamond/2.1.8 hmmer/3.3.2
cd $R
cp -n $C/1000files genomes.txt
{ diamond version; hmmsearch -h | sed -n 2p; } > tool_versions.txt
[[ -s db/reference_loo9.dmnd ]] || diamond makedb --in db/reference_loo9.fasta -d db/reference_loo9 --threads 4
[[ -s db/Pfam-A.hmm.h3m ]] || hmmpress db/Pfam-A.hmm
grep -c '^NAME' db/Pfam-A.hmm > db/pfam_nmodels.txt
# GFF3 for operons (same RefSeq assembly as done_proteins/*.faa)
mkdir -p gff operons
while read -r g; do
  [[ -s gff/${g}_genomic.gff.gz ]] && continue
  acc=${g%%_ASM*}; acc=$(echo "$g" | cut -d_ -f1,2); num=${acc#GCF_}; num=${num%%.*}
  url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${num:0:3}/${num:3:3}/${num:6:3}/${g}/${g}_genomic.gff.gz"
  curl -sSfL --retry 3 -o gff/${g}_genomic.gff.gz.tmp "$url" && mv gff/${g}_genomic.gff.gz.tmp gff/${g}_genomic.gff.gz || echo "GFF_FAIL $g" >> gff/failed.txt
done < genomes.txt
while read -r g; do
  [[ -s gff/${g}_genomic.gff.gz && ! -s operons/${g}_operons.tsv ]] || continue
  zcat gff/${g}_genomic.gff.gz > $TMPDIR/${g}.gff
  python3 $P/make_operons.py $TMPDIR/${g}.gff operons/${g}_operons.tsv 300 > /dev/null
  rm -f $TMPDIR/${g}.gff
done < genomes.txt
unzip -p bin/gspa-unified.jar taxon-constraints/ncbi-taxon-hierarchy.tsv > ref/ncbi-taxon-hierarchy.from-jar.tsv
python3 $P/01_taxon_map.py genomes.txt $C/taxa_constraints/propagated_taxon_ids.tsv ref/ncbi-taxon-hierarchy.from-jar.tsv taxon_map.tsv
sha256sum bin/gspa-unified.jar db/reference_loo9.fasta db/reference_loo9.dmnd db/Pfam-A.hmm ref/go.owl ref/pfam2go.txt ref/kegg_pathways.tsv ref/kegg_modules.tsv ref/ec2go.txt > input_sha256_ibex.txt
echo "gff ok: $(ls gff/*.gff.gz | wc -l); operons: $(ls operons | wc -l)"
