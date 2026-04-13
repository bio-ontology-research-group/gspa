#!/bin/bash
# bench10_pipeline.sh — Complete benchmark of GSPA vs PGAP on 10 new genomes
# Runs entirely locally on unimatrix01, no APIs.
set -euo pipefail

source /storage/miniforge3/etc/profile.d/conda.sh
conda activate bench10

JAVA=/data/hohndor/envs/jdk-21.0.10+7/bin/java
JAR=/data/hohndor/gspa/bin/gspa.jar
REF=/data/hohndor/gspa/reference
BP=/data/hohndor/gspa/benchmark-py
GOA=/data/hohndor/gspa/benchmark/goa_uniprot_all.gaf.gz
COLLAB=${REF}/gene_refseq_uniprotkb_collab.gz
SPROT=/data/hohndor/gspa/benchmark/uniprot_sprot.fasta.gz
PFAM_DB=/storage/software/databases/hmmer/Pfam-A.hmm
IPRSCAN=/data/hohndor/envs/interproscan-5.72-103.0/interproscan.sh

ROOT=/data/hohndor/gspa/proteomes/bench10
mkdir -p ${ROOT}
cd ${ROOT}

# ============================================================
# GENOME DEFINITIONS
# ============================================================
TAGS=(vcholerae saureus spneumoniae ccrescentus rprowazekii tpallidum tthermophilus dradiodurans scoelicolor pfuriosus)

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

# ============================================================
# STEP 1: Download genomes
# ============================================================
echo "========================================"
echo "STEP 1: Downloading 10 genomes"
echo "========================================"
for tag in "${TAGS[@]}"; do
  echo "--- ${tag} ---"
  base="${GENOME_FTP[$tag]}"
  aname="${base##*/}"
  [[ -s ${tag}_genomic.fna ]] || { curl -sSL "${base}/${aname}_genomic.fna.gz" | gunzip > ${tag}_genomic.fna; }
  [[ -s ${tag}_genomic.gff ]] || { curl -sSL "${base}/${aname}_genomic.gff.gz" | gunzip > ${tag}_genomic.gff; }
  [[ -s ${tag}_protein.faa ]] || { curl -sSL "${base}/${aname}_protein.faa.gz" | gunzip > ${tag}_protein.faa; }
  echo "  proteins: $(grep -c '^>' ${tag}_protein.faa), GO lines: $(grep -cE 'go_function|go_process|go_component' ${tag}_genomic.gff || echo 0)"
done

# ============================================================
# STEP 2: Build RefSeq→UniProt mappings from local NCBI collab file
# ============================================================
echo "========================================"
echo "STEP 2: RefSeq→UniProt mappings (local collab file)"
echo "========================================"
mkdir -p maps

# Extract RefSeq IDs from protein FASTAs
for tag in "${TAGS[@]}"; do
  grep '^>' ${tag}_protein.faa | awk '{print substr($1,2)}' > ${tag}_refseq_ids.txt
done

# Scan collab file once: extract all RefSeq→UniProt mappings for our genomes
if [[ ! -s maps/.done ]]; then
  python3 -c "
import gzip, sys, os
# Load all RefSeq IDs per genome
genome_refseq = {}
tags = '${TAGS[*]}'.split()
for tag in tags:
    ids = set()
    with open(f'{tag}_refseq_ids.txt') as f:
        for line in f:
            ids.add(line.strip())
    genome_refseq[tag] = ids
    print(f'  {tag}: {len(ids):,} RefSeq IDs', flush=True)
all_refseq = set()
for s in genome_refseq.values():
    all_refseq |= s
print(f'  union: {len(all_refseq):,} RefSeq IDs', flush=True)

genome_maps = {t: {} for t in tags}
n = 0; kept = 0
with gzip.open('${COLLAB}', 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        n += 1
        if n % 20_000_000 == 0:
            print(f'  {n:,} rows, {kept:,} kept', flush=True)
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 5: continue
        refseq, uniprot, method = fields[0], fields[1], fields[4]
        if refseq not in all_refseq: continue
        for tag in tags:
            if refseq in genome_refseq[tag]:
                existing = genome_maps[tag].get(refseq)
                if existing is None or (existing[1] != 'identical' and method == 'identical'):
                    genome_maps[tag][refseq] = (uniprot, method)
                    kept += 1

print(f'Total: {n:,} rows, {kept:,} kept', flush=True)
os.makedirs('maps', exist_ok=True)
for tag in tags:
    out_path = f'maps/{tag}.refseq_to_uniprot.tsv'
    with open(out_path, 'w') as out:
        for refseq, (uniprot, method) in sorted(genome_maps[tag].items()):
            out.write(f'{refseq}\t{uniprot}\n')
    print(f'  {tag}: {len(genome_maps[tag]):,} mappings -> {out_path}')
"
  touch maps/.done
fi

for tag in "${TAGS[@]}"; do
  echo "  ${tag}: $(wc -l < maps/${tag}.refseq_to_uniprot.tsv 2>/dev/null || echo 0) mappings"
done

# Build UniProt accession lists from the mappings
for tag in "${TAGS[@]}"; do
  awk '{print $2}' maps/${tag}.refseq_to_uniprot.tsv | sort -u > ${tag}_uniprot_accs.txt
  echo "  ${tag} UniProt accs: $(wc -l < ${tag}_uniprot_accs.txt)"
done

# ============================================================
# STEP 3: Build UniProt protein FASTA for DIAMOND
# ============================================================
echo "========================================"
echo "STEP 3: Preparing protein FASTA (UniProt IDs) for DIAMOND"
echo "========================================"
# Use the NCBI protein FASTA but remap headers to UniProt accessions
for tag in "${TAGS[@]}"; do
  if [[ ! -s ${tag}.faa ]]; then
    python3 -c "
import sys
mapping = {}
with open('maps/${tag}.refseq_to_uniprot.tsv') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            mapping[parts[0]] = parts[1]
seen = set()
with open('${tag}_protein.faa') as fin, open('${tag}.faa', 'w') as fout:
    current_acc = None
    for line in fin:
        if line.startswith('>'):
            refseq = line.split()[0][1:]
            uni = mapping.get(refseq)
            if uni and uni not in seen:
                current_acc = uni
                seen.add(uni)
                fout.write(f'>{uni} {line[1:].strip()}\n')
            else:
                current_acc = None
        elif current_acc:
            fout.write(line)
print(f'  ${tag}: {len(seen)} UniProt-mapped proteins')
"
  else
    echo "  ${tag}.faa already exists ($(grep -c '^>' ${tag}.faa) proteins)"
  fi
done

# ============================================================
# STEP 4: Build leave-19-out DIAMOND reference DB
# ============================================================
echo "========================================"
echo "STEP 4: Building leave-19-out DIAMOND DB"
echo "========================================"
# Combine exclusion lists from original 9 + new 10 genomes
cat /data/hohndor/gspa/proteomes/all9_excluded.txt > all19_excluded.txt
for tag in "${TAGS[@]}"; do
  cat ${tag}_uniprot_accs.txt >> all19_excluded.txt
done
sort -u all19_excluded.txt -o all19_excluded.txt
echo "  total excluded: $(wc -l < all19_excluded.txt) accessions"

if [[ ! -s reference_loo19.dmnd ]]; then
  echo "  filtering Swiss-Prot..."
  python3 ${BP}/filter_fasta_by_exclude.py \
    --fasta ${SPROT} \
    --exclude all19_excluded.txt \
    --output reference_loo19.fasta
  echo "  building DIAMOND DB..."
  diamond makedb --in reference_loo19.fasta --db reference_loo19 --threads 16
  echo "  done: $(grep -c '^>' reference_loo19.fasta) proteins in reference"
else
  echo "  reference_loo19.dmnd already exists"
fi

# ============================================================
# STEP 5: Run predictors (DIAMOND + HMMER + InterProScan)
# ============================================================
echo "========================================"
echo "STEP 5: Running predictors"
echo "========================================"
mkdir -p preds

for tag in "${TAGS[@]}"; do
  echo "--- ${tag} ---"
  faa=${ROOT}/${tag}.faa

  # DIAMOND blastp
  if [[ ! -s preds/${tag}_diamond.tsv ]]; then
    echo "  DIAMOND..."
    diamond blastp \
      --query ${faa} \
      --db ${ROOT}/reference_loo19 \
      --out preds/${tag}_diamond.tsv \
      --outfmt 6 qseqid sseqid pident length mismatch gapopen evalue bitscore stitle \
      --max-target-seqs 10 \
      --evalue 1e-5 \
      --threads 16
    echo "    $(wc -l < preds/${tag}_diamond.tsv) hits"
  else
    echo "  DIAMOND already done ($(wc -l < preds/${tag}_diamond.tsv) hits)"
  fi

  # HMMER/Pfam
  if [[ ! -s preds/${tag}_pfam.domtbl ]]; then
    echo "  HMMER/Pfam..."
    hmmsearch \
      --domtblout preds/${tag}_pfam.domtbl \
      --noali -E 1e-5 --cpu 16 \
      ${PFAM_DB} ${faa} > /dev/null
    echo "    $(grep -c '^[^#]' preds/${tag}_pfam.domtbl) domain hits"
  else
    echo "  HMMER already done"
  fi

  # InterProScan
  if [[ ! -s preds/${tag}_interproscan.tsv ]]; then
    echo "  InterProScan..."
    conda deactivate
    conda activate interproscan
    ${IPRSCAN} -i ${faa} -o preds/${tag}_interproscan.tsv -f TSV --goterms --cpu 8 -dp 2>/dev/null
    conda deactivate
    conda activate bench10
    echo "    $(wc -l < preds/${tag}_interproscan.tsv) lines"
  else
    echo "  InterProScan already done ($(wc -l < preds/${tag}_interproscan.tsv) lines)"
  fi
done

# ============================================================
# STEP 6: Parse claims (DIAMOND + Pfam + InterProScan → JSONL)
# ============================================================
echo "========================================"
echo "STEP 6: Parsing claims"
echo "========================================"
mkdir -p claims

# Build GOA subset for DIAMOND target lookup
echo "  building diamond_targets list..."
for tag in "${TAGS[@]}"; do
  awk -F'\t' '{print $2}' preds/${tag}_diamond.tsv | sort -u
done | sort -u > diamond_targets.txt
echo "  $(wc -l < diamond_targets.txt) unique DIAMOND targets"

# Build GOA subset for the new 10 genomes (scan full GOA once, filter by diamond targets)
if [[ ! -s goa_subset.gaf ]]; then
  echo "  scanning GOA for DIAMOND target annotations..."
  python3 -c "
import gzip, sys
targets = set()
with open('diamond_targets.txt') as f:
    for line in f:
        t = line.strip()
        # Extract bare accession from sp|ACC|NAME format
        if '|' in t:
            t = t.split('|')[1]
        targets.add(t)
print(f'  {len(targets)} unique targets', flush=True)
n = 0; kept = 0
with gzip.open('${GOA}', 'rt') as fin, open('goa_subset.gaf', 'w') as fout:
    for line in fin:
        n += 1
        if n % 50_000_000 == 0:
            print(f'  {n:,} lines, kept {kept:,}', flush=True)
        if line.startswith('!'):
            continue
        fields = line.split('\t', 5)
        if len(fields) < 5:
            continue
        if fields[1] in targets:
            fout.write(line)
            kept += 1
print(f'  Total: {n:,} lines, kept {kept:,}', flush=True)
"
fi
echo "  GOA subset: $(wc -l < goa_subset.gaf) lines"

for tag in "${TAGS[@]}"; do
  echo "--- ${tag} ---"

  # Organize predictor outputs into per-genome results dir for 02b parser
  resdir=preds/${tag}_results
  mkdir -p ${resdir}
  ln -sf ${ROOT}/preds/${tag}_diamond.tsv ${resdir}/diamond_results.tsv 2>/dev/null || true
  ln -sf ${ROOT}/preds/${tag}_pfam.domtbl ${resdir}/pfam_results.domtbl 2>/dev/null || true

  # Parse DIAMOND + Pfam claims via standard parser
  if [[ ! -s claims/${tag}_claims.jsonl ]]; then
    python3 ${BP}/02b_parse_predictors_to_claims.py \
      --results-dir ${resdir} \
      --goa goa_subset.gaf \
      --pfam2go ${REF}/pfam2go.txt \
      --test-accs ${ROOT}/${tag}_uniprot_accs.txt \
      --output claims/${tag}_claims.jsonl
  fi
  echo "  D+P claims: $(wc -l < claims/${tag}_claims.jsonl)"

  # Parse InterProScan claims
  if [[ -s preds/${tag}_interproscan.tsv ]] && [[ ! -s claims/${tag}_iprscan_claims.jsonl ]]; then
    python3 ${BP}/parse_interproscan_claims.py \
      --input preds/${tag}_interproscan.tsv \
      --output claims/${tag}_iprscan_claims.jsonl
  fi
  [[ -s claims/${tag}_iprscan_claims.jsonl ]] && echo "  IPR claims: $(wc -l < claims/${tag}_iprscan_claims.jsonl)"

  # Merge all claims
  cat claims/${tag}_claims.jsonl > claims/${tag}_merged.jsonl
  [[ -s claims/${tag}_iprscan_claims.jsonl ]] && cat claims/${tag}_iprscan_claims.jsonl >> claims/${tag}_merged.jsonl
  echo "  merged: $(wc -l < claims/${tag}_merged.jsonl) total claims"
done

# ============================================================
# STEP 7: Extract PGAP GO annotations from GFF
# ============================================================
echo "========================================"
echo "STEP 7: Extracting PGAP GO from GFF"
echo "========================================"
mkdir -p pgap

for tag in "${TAGS[@]}"; do
  if [[ ! -s pgap/${tag}_pgap.tsv ]]; then
    python3 ${BP}/extract_pgap_annotations.py \
      --gff ${ROOT}/${tag}_genomic.gff \
      --refseq-to-uniprot ${ROOT}/maps/${tag}.refseq_to_uniprot.tsv \
      --go-owl-aspect-map ${REF}/go_aspect_map.tsv \
      --out pgap/${tag}_pgap.tsv
  fi
  echo "  ${tag}: $(wc -l < pgap/${tag}_pgap.tsv) PGAP annotations"
done

# ============================================================
# STEP 8: Extract GOA ground truth (experimental + all)
# ============================================================
echo "========================================"
echo "STEP 8: Extracting GOA ground truth"
echo "========================================"
mkdir -p truth

python3 ${BP}/extract_goa_dual.py \
  --goa ${GOA} \
  --accessions $(for t in "${TAGS[@]}"; do echo -n "${t}:${ROOT}/${t}_uniprot_accs.txt "; done) \
  --out-dir ${ROOT}/truth/

for tag in "${TAGS[@]}"; do
  echo "  ${tag} exp: $(wc -l < truth/${tag}_truth_exp.tsv) all: $(wc -l < truth/${tag}_truth_all.tsv)"
done

# ============================================================
# STEP 9: Run GSPA integration (Noisy-OR + priors)
# ============================================================
echo "========================================"
echo "STEP 9: GSPA integration"
echo "========================================"
mkdir -p integrated

for tag in "${TAGS[@]}"; do
  echo "--- ${tag} ---"
  k=${GENOME_KINGDOM[$tag]}
  claims=claims/${tag}_merged.jsonl
  [[ -s ${claims} ]] || { echo "  no claims, skip"; continue; }

  # Baseline (no priors)
  if [[ ! -s integrated/${tag}_integrated.tsv ]]; then
    ${JAVA} -jar ${JAR} integrate \
      --claims ${claims} --out integrated/${tag}_integrated.tsv \
      --go-owl ${REF}/go.owl --lite \
      > integrated/${tag}_base.log 2>&1
  fi

  # With priors (essentiality + coherence)
  if [[ ! -s integrated/${tag}_integrated_priors.tsv ]]; then
    ${JAVA} -jar ${JAR} integrate \
      --claims ${claims} --out integrated/${tag}_integrated_priors.tsv \
      --go-owl ${REF}/go.owl --lite \
      --essential-profile ${k} \
      --pathways ${REF}/kegg_pathways.tsv --ec2go ${REF}/ec2go.txt \
      --enable-priors essentiality,coherence \
      > integrated/${tag}_priors.log 2>&1
  fi

  echo "  base=$(wc -l < integrated/${tag}_integrated.tsv) priors=$(wc -l < integrated/${tag}_integrated_priors.tsv)"
done

# ============================================================
# STEP 10: Compute F-max with bootstrap CIs
# ============================================================
echo "========================================"
echo "STEP 10: F-max comparison (GSPA vs PGAP)"
echo "========================================"
mkdir -p results

for tag in "${TAGS[@]}"; do
  truth_exp=truth/${tag}_truth_exp.tsv
  truth_all=truth/${tag}_truth_all.tsv
  priors=integrated/${tag}_integrated_priors.tsv
  pgap_file=pgap/${tag}_pgap.tsv

  [[ -s ${priors} ]] || { echo "  ${tag}: no integrated file, skip"; continue; }
  [[ -s ${truth_all} ]] || { echo "  ${tag}: no truth, skip"; continue; }

  pgap_arg=''
  if [[ -s ${pgap_file} ]] && [[ $(wc -l < ${pgap_file}) -gt 2 ]]; then
    pgap_arg="--pgap ${pgap_file}"
  fi

  truth_args="--truth all:${truth_all}"
  [[ -s ${truth_exp} ]] && [[ $(wc -l < ${truth_exp}) -gt 2 ]] && truth_args="${truth_args} --truth exp:${truth_exp}"

  python3 ${BP}/benchmark_pgap_v2.py \
    ${truth_args} \
    --gspa integrated/${tag}_integrated.tsv \
    --gspa-priors ${priors} \
    ${pgap_arg} \
    --tag ${tag} --n-bootstrap 200 \
    > results/${tag}_fmax.json 2>/dev/null

  echo "  ${tag}: done"
done

# ============================================================
# STEP 11: Print results
# ============================================================
echo "========================================"
echo "FINAL RESULTS"
echo "========================================"
python3 ${BP}/print_fmax.py results/*_fmax.json

echo ""
echo "ALL DONE. Results in ${ROOT}/results/"
