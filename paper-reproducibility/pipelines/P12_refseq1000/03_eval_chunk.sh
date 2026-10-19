#!/bin/bash
# P12 step 03 (IBEX, one chunk of genomes per Slurm array task): score every
# method with (a) the annot_evaluation structural suite re-implemented in
# benchmark/unified/pillar1_5_6/p12_structural_suite.py and (b) GSPA's own
# evaluator (`gspa-unified.jar evaluate --lite`, the GAEF path of
# benchmark/unified/pillar1_2/03_run_gaef_eval.py) on the leaf GO sets.
#
# usage: METHODS=tools|gspa|all CHUNK=25 [GAEF_METHODS=regex] sbatch --array=0-39 03_eval_chunk.sh
set -euo pipefail
R=${R:-/ibex/scratch/projects/c2014/hohndor/refseq1000_repro_2026-10-19}
C=/ibex/scratch/projects/c2014/annot_evaluation
JAVA=${JAVA:-/ibex/scratch/projects/c2014/rob/jdk-21/bin/java}
JAR=$R/bin/gspa-unified.jar
REF=$R/ref
P=$R/scripts
METHODS=${METHODS:-tools}
CHUNK=${CHUNK:-25}
K=${SLURM_ARRAY_TASK_ID:-0}
T=${SLURM_CPUS_PER_TASK:-4}
OUT=$R/eval/$METHODS
mkdir -p $OUT/suite $OUT/leaf $OUT/gaef $OUT/manifests
mapfile -t GS < <(sed -n "$(( K*CHUNK + 1 )),$(( (K+1)*CHUNK ))p" $R/genomes.txt)
[[ ${#GS[@]} -gt 0 ]] || exit 0
if [[ $METHODS == gspa ]]; then   # only genomes whose GSPA run finished
  OK=(); for g in "${GS[@]}"; do [[ -s $R/work/$g/DONE ]] && OK+=("$g"); done; GS=("${OK[@]}")
  [[ ${#GS[@]} -gt 0 ]] || { echo "no finished GSPA genomes in chunk $K"; exit 0; }
fi
J=$OUT/suite/chunk_${K}.jsonl
if [[ ${GAEF_METHODS:-.} == . || ! -s $J ]]; then   # a GAEF-only re-run keeps the suite output
  rm -f $J.tmp
  /usr/bin/time -v -o $OUT/suite/chunk_${K}.time python3 $P/p12_structural_suite.py --genome "${GS[@]}" \
    --methods $METHODS --annot-eval $C --gspa-work $R/work --out-jsonl $J.tmp --leaf-tsv-dir $OUT/leaf
  mv $J.tmp $J
fi
MAN=$OUT/manifests/chunk_${K}.tsv
printf 'gff_path\tgaf_path\tout_json_path\tkingdom\ttaxon\n' > $MAN
for g in "${GS[@]}"; do
  TAXON=$(awk -F'\t' -v g="$g" '$1==g{print $3}' $R/taxon_map.tsv); TAXON=${TAXON:-NCBITaxon_2}
  for lf in $OUT/leaf/${g}.*.leaf.tsv; do
    m=$(basename $lf .leaf.tsv); m=${m#${g}.}
    [[ $m =~ ${GAEF_METHODS:-.} ]] || continue
    # DeepFRI predictions are keyed by UniProt accession (structure models), not by
    # the RefSeq WP_ ids of the FASTA: the synthetic GFF lists the FASTA ids plus
    # every accession of the leaf set, so no prediction is dropped.
    IDS=$OUT/gaef/${g}.${m}.ids.faa
    { grep '^>' $C/done_proteins/${g}_protein.faa | cut -d' ' -f1; \
      tail -n +2 $lf | cut -f1 | sort -u | sed 's/^/>/'; } | awk '!s[$0]++' > $IDS
    python3 $P/make_synth_gff_gaf.py --fasta $IDS --tsv $lf \
      --schema pgap --gff-out $OUT/gaef/${g}.${m}.gff --gaf-out $OUT/gaef/${g}.${m}.gaf > /dev/null
    printf '%s\t%s\t%s\tbacteria\t%s\n' $OUT/gaef/${g}.${m}.gff $OUT/gaef/${g}.${m}.gaf \
      $OUT/gaef/${g}.${m}_quality.json $TAXON >> $MAN
  done
done
/usr/bin/time -v -o $OUT/gaef/chunk_${K}.time $JAVA -Xmx8g -XX:ActiveProcessorCount=$T -jar $JAR evaluate --lite \
  --go-owl $REF/go.owl --batch-manifest $MAN -k bacteria \
  --pathways $REF/kegg_pathways.tsv --modules $REF/kegg_modules.tsv --ec2go $REF/ec2go.txt \
  -o $OUT/gaef/chunk_${K}_unused.json > $OUT/gaef/chunk_${K}.log 2>&1
for g in "${GS[@]}"; do rm -f $OUT/gaef/${g}.*.gff $OUT/gaef/${g}.*.ids.faa; done   # synthetic, regenerated on demand
echo "chunk $K: ${#GS[@]} genomes done"
