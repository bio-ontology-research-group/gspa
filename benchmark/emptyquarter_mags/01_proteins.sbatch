#!/bin/bash
#SBATCH --job-name=eq-proteins
#SBATCH --partition=batch
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
# One genome per row (array tasks take blocks of rows, see below): annot/<mag>.{faa,gff3,fna} from the Bakta annotation of the metagenomics tree
# when present (symlinks), otherwise pyrodigal 3.7.1 on the bin FASTA. Idempotent.
set -euo pipefail
# start from a clean module environment: the submitting shell may export LD_LIBRARY_PATH/PYTHONPATH of
# loaded modules (python/gcc), which breaks the conda tools and the static-linked gapsmith binary
module purge >/dev/null 2>&1 || true; unset LD_LIBRARY_PATH PYTHONPATH PYTHONHOME
source ${OUT:-/encrypted/e3008/gspa/eq_mags}/code/common.sh
# Array tasks work on blocks of GENOMES_PER_TASK consecutive rows of $BINS (ibex caps a user at 2,000 submitted
# jobs and every array task counts): the task re-runs this script once per row with ROW set.
if [[ -z "${ROW:-}" ]]; then
  K=${GENOMES_PER_TASK:-50}; fail=0; first=$(( SLURM_ARRAY_TASK_ID * K ))
  for (( r = first; r < first + K; r++ )); do ROW=$r bash "$0" || fail=$(( fail + 1 )); done
  echo "block ${SLURM_ARRAY_TASK_ID} (rows $first-$(( first + K - 1 ))): $fail failures $(date -Is)"; exit $(( fail > 0 ))
fi
bin_row $ROW || { [[ $ROW -ge $(( $(wc -l < $BINS) - 1 )) ]] && exit 0; echo "no genome for row $ROW"; exit 1; }
T=${SLURM_CPUS_PER_TASK:-2}
D=$WORK/$G; mkdir -p $D/annot $D/times; cd $D
[[ -s DONE_proteins ]] && { echo "[skip] $G proteins done"; exit 0; }
if [[ ! -s annot/$G.faa || ! -s annot/$G.gff3 || ! -s annot/$G.fna ]]; then
  if [[ -n "$BD" && -s $BD/$BIN.faa ]]; then
    ln -sfn $BD/$BIN.faa annot/$G.faa
    if [[ -s $BD/$BIN.gff3 ]]; then ln -sfn $BD/$BIN.gff3 annot/$G.gff3; else zcat $BD/$BIN.gff3.gz > annot/$G.gff3.part && mv annot/$G.gff3.part annot/$G.gff3; fi
    if [[ -s $BD/$BIN.fna ]]; then ln -sfn $BD/$BIN.fna annot/$G.fna; else zcat $BD/$BIN.fna.gz > annot/$G.fna.part && mv annot/$G.fna.part annot/$G.fna; fi
    echo bakta > annot/SOURCE
  else
    gzip -cdf $FA > annot/$G.fna.part && mv annot/$G.fna.part annot/$G.fna   # bins are .fa.gz, culture assemblies plain .fa
    # the conda pyrodigal must not see a loaded python module's environment (SRE module mismatch otherwise)
    /usr/bin/time -v -o times/pyrodigal.time env -u PYTHONPATH -u PYTHONHOME -u PYTHONSTARTUP -u PYTHON_LIBS \
      $PYRODIGAL -i annot/$G.fna -a annot/$G.faa.part -f gff -o annot/$G.gff3.part -p single -j $T > pyrodigal.log 2>&1
    mv annot/$G.faa.part annot/$G.faa; mv annot/$G.gff3.part annot/$G.gff3
    echo pyrodigal > annot/SOURCE
  fi
fi
NPROT=$(grep -c '^>' annot/$G.faa)
[[ $NPROT -gt 0 ]] || { echo "$G: empty proteome"; exit 1; }
sha256sum $FA annot/$G.faa annot/$G.gff3 annot/$G.fna > input_sha256.txt
printf '%s\t%s\t%s\t%s\n' "$G" "$(cat annot/SOURCE)" "$NPROT" "$(date -Is)" > DONE_proteins
echo "DONE_proteins $G $(cat annot/SOURCE) proteins=$NPROT"
