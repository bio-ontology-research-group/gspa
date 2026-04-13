#!/bin/bash
# Download 10 new genomes with PGAP GO annotations for benchmarking
# All genomes selected because their RefSeq GFF contains go_function/go_process/go_component
set -euo pipefail

ROOT=/data/hohndor/gspa/proteomes/bench10
mkdir -p ${ROOT}
cd ${ROOT}

# Genome definitions: tag | accession | assembly_name | FTP_path_suffix | kingdom
declare -A GENOME_FTP=(
  [vcholerae]='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/745/GCF_000006745.1_ASM674v1'
  [saureus]='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/645/GCF_000009645.1_ASM964v1'
  [spneumoniae]='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/885/GCF_000006885.1_ASM688v1'
  [ccrescentus]='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/905/GCF_000006905.1_ASM690v1'
  [rprowazekii]='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/735/GCF_000195735.1_ASM19573v1'
  [tpallidum]='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/605/GCF_000008605.1_ASM860v1'
  [tthermophilus]='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/091/545/GCF_000091545.1_ASM9154v1'
  [dradiodurans]='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/565/GCF_000008565.1_ASM856v1'
  [scoelicolor]='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/203/835/GCF_000203835.1_ASM20383v1'
  [pfuriosus]='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/305/GCF_000007305.1_ASM730v1'
)

declare -A GENOME_ACC=(
  [vcholerae]='GCF_000006745.1'
  [saureus]='GCF_000009645.1'
  [spneumoniae]='GCF_000006885.1'
  [ccrescentus]='GCF_000006905.1'
  [rprowazekii]='GCF_000195735.1'
  [tpallidum]='GCF_000008605.1'
  [tthermophilus]='GCF_000091545.1'
  [dradiodurans]='GCF_000008565.1'
  [scoelicolor]='GCF_000203835.1'
  [pfuriosus]='GCF_000007305.1'
)

declare -A GENOME_KINGDOM=(
  [vcholerae]='bacteria'
  [saureus]='bacteria'
  [spneumoniae]='bacteria'
  [ccrescentus]='bacteria'
  [rprowazekii]='bacteria'
  [tpallidum]='bacteria'
  [tthermophilus]='bacteria'
  [dradiodurans]='bacteria'
  [scoelicolor]='bacteria'
  [pfuriosus]='archaea'
)

TAGS=(vcholerae saureus spneumoniae ccrescentus rprowazekii tpallidum tthermophilus dradiodurans scoelicolor pfuriosus)

for tag in "${TAGS[@]}"; do
  echo "=== ${tag} (${GENOME_ACC[$tag]}) ==="
  base="${GENOME_FTP[$tag]}"
  aname="${base##*/}"

  # 1. Download genomic FASTA
  if [[ ! -s ${tag}_genomic.fna ]]; then
    echo "  downloading genomic FASTA..."
    curl -sSL "${base}/${aname}_genomic.fna.gz" -o ${tag}_genomic.fna.gz
    gunzip -f ${tag}_genomic.fna.gz
  fi
  echo "  FASTA: $(grep -c '^>' ${tag}_genomic.fna) contigs"

  # 2. Download GFF
  if [[ ! -s ${tag}_genomic.gff ]]; then
    echo "  downloading GFF..."
    curl -sSL "${base}/${aname}_genomic.gff.gz" -o ${tag}_genomic.gff.gz
    gunzip -f ${tag}_genomic.gff.gz
  fi
  echo "  GFF: $(wc -l <${tag}_genomic.gff) lines"
  echo "  GO annotations: $(grep -cE 'go_function|go_process|go_component' ${tag}_genomic.gff) lines"

  # 3. Download NCBI protein FASTA (RefSeq protein IDs)
  if [[ ! -s ${tag}_protein.faa ]]; then
    echo "  downloading protein FASTA..."
    curl -sSL "${base}/${aname}_protein.faa.gz" -o ${tag}_protein.faa.gz
    gunzip -f ${tag}_protein.faa.gz
  fi
  echo "  proteins: $(grep -c '^>' ${tag}_protein.faa)"

  # 4. Extract RefSeq protein IDs
  grep '^>' ${tag}_protein.faa | awk '{print substr($1,2)}' > ${tag}_refseq_ids.txt
  echo "  RefSeq IDs: $(wc -l <${tag}_refseq_ids.txt)"

  echo ""
done

echo "=== ALL DOWNLOADS COMPLETE ==="
wc -l *_refseq_ids.txt
