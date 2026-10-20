#!/bin/bash
#SBATCH --job-name=eq-channels
#SBATCH --partition=batch
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
# One bin per array task (index = data row of $BINS, 0-based): proteins (Bakta annotation when
# present, pyrodigal otherwise) -> DIAMOND vs experimental-GO Swiss-Prot, HMMER Pfam-A, InterProScan
# -> Stage 0 claims (DIAMOND + Pfam + InterProScan). Idempotent: finished steps are skipped.
set -euo pipefail
source ${OUT:-/encrypted/e3008/gspa/eq_mags}/code/common.sh
module load python/3.11.0 diamond/2.1.16 hmmer/3.3.2 >/dev/null 2>&1
bin_row ${SLURM_ARRAY_TASK_ID} || { echo "no bin for index ${SLURM_ARRAY_TASK_ID}"; exit 0; }
T=${SLURM_CPUS_PER_TASK:-8}
D=$WORK/$G; mkdir -p $D/annot $D/times; cd $D
[[ -s DONE_channels ]] && { echo "[skip] $G channels done"; exit 0; }
echo "== $G ($CLASS, $DOMAIN, $TAXON, gram $GRAM) $(date -Is) host=$(hostname) cpus=$T"
TV="/usr/bin/time -v"; run() { local n=$1; shift; $TV -o times/$n.time "$@"; }

# --- proteins: annot/$G.{faa,gff3,fna} ------------------------------------------------------
if [[ ! -s annot/$G.faa || ! -s annot/$G.gff3 || ! -s annot/$G.fna ]]; then
  if [[ -n "$BD" && -s $BD/$BIN.faa ]]; then
    ln -sf $BD/$BIN.faa annot/$G.faa
    if [[ -s $BD/$BIN.gff3 ]]; then ln -sf $BD/$BIN.gff3 annot/$G.gff3; else zcat $BD/$BIN.gff3.gz > annot/$G.gff3.part && mv annot/$G.gff3.part annot/$G.gff3; fi
    if [[ -s $BD/$BIN.fna ]]; then ln -sf $BD/$BIN.fna annot/$G.fna; else zcat $BD/$BIN.fna.gz > annot/$G.fna.part && mv annot/$G.fna.part annot/$G.fna; fi
    echo bakta > annot/SOURCE
  else
    zcat $FA > annot/$G.fna.part && mv annot/$G.fna.part annot/$G.fna
    # the conda pyrodigal must not see the python module's environment (SRE module mismatch otherwise)
    run pyrodigal env -u PYTHONPATH -u PYTHONHOME -u PYTHONSTARTUP -u PYTHON_LIBS $PYRODIGAL -i annot/$G.fna -a annot/$G.faa.part -f gff -o annot/$G.gff3.part -p single -j $T > pyrodigal.log 2>&1
    mv annot/$G.faa.part annot/$G.faa; mv annot/$G.gff3.part annot/$G.gff3
    echo pyrodigal > annot/SOURCE
  fi
fi
FAA=annot/$G.faa
NPROT=$(grep -c '^>' $FAA)
sha256sum $FA $FAA annot/$G.gff3 annot/$G.fna > input_sha256.txt
{ diamond version; hmmsearch -h | sed -n 2p; $PYRODIGAL --version; } > tool_versions.txt 2>&1 || true

# --- channels -------------------------------------------------------------------------------
[[ -s ${G}_diamond.tsv ]] || { run diamond diamond blastp --query $FAA --db $REF/sprot_exp.dmnd --threads $T \
    --outfmt 6 --evalue 1e-5 --max-target-seqs 5 --quiet --out ${G}_diamond.tsv.part \
    && mv ${G}_diamond.tsv.part ${G}_diamond.tsv; }
[[ -s ${G}_pfam.domtbl ]] || { run hmmsearch hmmsearch --domtblout ${G}_pfam.domtbl.part --noali --cpu $T --cut_ga \
    $PFAM_HMM $FAA > /dev/null && mv ${G}_pfam.domtbl.part ${G}_pfam.domtbl; }
if [[ ! -s ${G}_ips.tsv ]]; then
  module load $IPS_MODULE >/dev/null 2>&1
  command -v interproscan.sh >/dev/null || { echo "interproscan.sh not on PATH after module load $IPS_MODULE"; exit 1; }
  TMPD=${TMPDIR:-/tmp}/ips_${G}_$$; mkdir -p $TMPD
  sed '/^>/!s/\*//g' $FAA > $TMPD/in.faa
  run interproscan interproscan.sh -i $TMPD/in.faa -f tsv -o ${G}_ips.tsv.part -appl $IPS_APPL -goterms -iprlookup -dp \
      -cpu $T -T $TMPD > ${G}_ips.log 2>&1
  mv ${G}_ips.tsv.part ${G}_ips.tsv; rm -rf $TMPD
  interproscan.sh --version 2>/dev/null | head -2 >> tool_versions.txt || true
fi
run claims python3 $BGC build-claims --obo $REF/go.obo --diamond ${G}_diamond.tsv --sprot-go $REF/sprot_exp_go.tsv \
    --domtbl ${G}_pfam.domtbl --pfam2go $REF/pfam2go.txt --ips ${G}_ips.tsv --out ${G}_claims.jsonl
NCLAIMS=$(wc -l < ${G}_claims.jsonl)
gzip -f ${G}_diamond.tsv ${G}_pfam.domtbl ${G}_ips.tsv 2>/dev/null || true
for f in times/*.time; do printf '%s\t%s\t%s\n' $(basename $f .time) \
  "$(grep 'Elapsed (wall' $f | awk '{print $NF}')" "$(grep 'Maximum resident' $f | awk '{print $NF}')"; done > timing_channels.tsv
printf '%s\t%s\t%s\t%s\t%s\n' "$G" "$(cat annot/SOURCE)" "$NPROT" "$NCLAIMS" "$(date -Is)" > DONE_channels
echo "DONE_channels $G proteins=$NPROT claims=$NCLAIMS $(date -Is)"
