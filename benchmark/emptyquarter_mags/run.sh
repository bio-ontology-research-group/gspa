#!/bin/bash
# Empty Quarter genome annotation with the full GSPA pipeline (ibex, Slurm; never the login node):
# the HQ/MQ soil bins and the culture genomes of the RAK-PR catalogue (enrichments, isolates, thermophiles).
# Clustered design: proteins per bin -> one linclust (90 % id, 90 % cov) over all proteomes -> DIAMOND, Pfam,
# InterProScan once per representative chunk -> claims scattered to bins -> gapsmith per bin -> GSPA per bin.
#
#   bash run.sh stage                 # local: copy code, jar and reference files to $OUT on ibex
#   bash run.sh select                # ibex job: build lists/bins.tsv (16,333 HQ+MQ soil bins + 722 RAK-PR culture genomes)
#   bash run.sh test [N]              # ibex, small test set (default 5 genomes): the whole chain under $OUT/test
#   bash run.sh launch [AFTER_JOBID]  # ibex: full chain; the proteins array starts after AFTER_JOBID ends
#                                     # (arrays take GENOMES_PER_TASK=50 genomes per task: ibex allows 2,000 submitted jobs)
#   bash run.sh rebalance THR ARRAY   # ibex: change an array's throttle (chunks / gapsmith / integrate)
#   bash run.sh status                # ibex: job states and DONE counts
#   bash run.sh compare [WORKDIR]     # ibex: gapsmith vs gapseq reaction sets for bins that have both
#   bash run.sh aggregate             # ibex job: results/eq_mags_{per_bin,summary}.tsv
#   bash run.sh fetch                 # local: copy results/ to benchmark/emptyquarter_mags/results/
#   python3 07_report.py              # local: results/eq_mags_{by_dataset,phyla}.tsv from the per-genome table
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd); root=$(cd "$here/../.." && pwd)
OUT=${OUT:-/encrypted/e3008/gspa/eq_mags}
source "$here/common.sh"
mode=${1:-help}
case "$mode" in
  stage)
    ssh ibex "mkdir -p $CODE $REF/fba $REF/essential $WORK $LISTS $LOGS $RES $OUT/bin $OUT/test"
    scp -q "$here"/*.sbatch "$here"/*.py "$here"/common.sh "$here"/run.sh \
      "$root/benchmark/unified/pillar3/mag_pathway_gaps.py" "$root/benchmark/unified/pillar2_fba_lro_price149/bgc_benchmark.py" \
      "$root/scripts/stoichiometric_locus_refiner.py" ibex:$CODE/
    scp -q "$root/gspa-unified.jar" ibex:$OUT/bin/gspa-unified.jar
    scp -q "$root/gspa-core/src/main/resources/taxon-constraints/ncbi-taxon-labels.tsv" ibex:$REF/
    ssh ibex "cp -n $P14_REF/{go.owl,go.obo,go-dag.tsv,ec2go.txt,pfam2go.txt,kegg_pathways.tsv,kegg_modules.tsv,sprot_exp.dmnd,sprot_exp_go.tsv,seed_reactions_corrected.tsv,diffusion_mets.tsv,seed_Enzyme_Class_Reactions_Aliases_unique.tsv} $REF/ 2>/dev/null; \
      cp -n $P14_REF/essential/* $REF/essential/ 2>/dev/null; cp -n $FBA_DATA_SRC/* $REF/fba/ 2>/dev/null; \
      sha256sum $OUT/bin/gspa-unified.jar $REF/*.* $REF/fba/* $REF/essential/* > $REF/sha256.txt; ls $REF $REF/fba | head -40" ;;
  select)
    ssh ibex "sbatch -p debug -t 00:30:00 -c 1 --mem=4G -J eq-select -o $LOGS/select_%j.out --wrap='module load python/3.11.0; \
      python3 $CODE/00_select.py --checkm2 $CHECKM2 --gtdb-bac $GTDB_BAC --gtdb-ar $GTDB_AR --labels $REF/ncbi-taxon-labels.tsv \
      --gapseq-out $GAPSEQ_OUT_HQ --catalogue $RAKPR_CATALOGUE --inputs $RAKPR_INPUTS --culture-bakta $RAKPR_BAKTA --out $LISTS/bins.tsv'" ;;
  test)
    n=${2:-5}; TO=$OUT/test
    ssh ibex "module load python/3.11.0; mkdir -p $TO/lists $TO/work $TO/logs $TO/cluster $TO/chunks $TO/claims_by_sample; python3 - <<'PY'
import csv
rows=list(csv.DictReader(open('$LISTS/bins.tsv'),delimiter='\t'))
pick=[]
for want in (lambda r: r['class']=='HQ' and r['gapseq_dir'] and r['bakta_dir'],
             lambda r: r['class']=='MQ' and r['bakta_dir'] and not r['gapseq_dir'],
             lambda r: r['dataset']=='soil_bin' and not r['bakta_dir'],
             lambda r: r['dataset']=='culture_isolate' and r['bakta_dir'],
             lambda r: r['dataset']=='culture_thermophilic'):
    for r in rows:
        if want(r) and r not in pick: pick.append(r); break
pick=(pick+[r for r in rows if r not in pick])[:$n]
with open('$TO/lists/bins.tsv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys()),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(pick)
print('\n'.join(r['mag_id']+' '+r['dataset']+' '+r['class']+' bakta='+str(bool(r['bakta_dir']))+' gapseq='+str(bool(r['gapseq_dir'])) for r in pick))
PY
      K=${GENOMES_PER_TASK_TEST:-2}; nt=\$(( ($n + K - 1) / K )); \
      E=ALL,OUT=$OUT,BINS=$TO/lists/bins.tsv,WORK=$TO/work,LISTS=$TO/lists,CL=$TO/cluster,CHUNKS=$TO/chunks,CBS=$TO/claims_by_sample,LOGS=$TO/logs,JOBS_TSV=$TO/jobs.tsv,CHUNK_SIZE=${CHUNK_SIZE_TEST:-3000},GENOMES_PER_TASK=\$K; \
      j1=\$(sbatch --parsable -p debug -t 01:00:00 --array=0-\$(( nt - 1 )) --export=\$E -o $TO/logs/pr_%A_%a.out $CODE/01_proteins.sbatch); \
      j5=\$(sbatch --parsable -p debug -t 01:00:00 --array=0-\$(( nt - 1 )) --dependency=aftercorr:\$j1 --export=\$E -o $TO/logs/gm_%A_%a.out $CODE/05_gapsmith.sbatch); \
      j2=\$(sbatch --parsable -p debug -t 02:00:00 -c 8 --mem=32G --dependency=afterok:\$j1 --export=\$E,GAPSMITH_JOB=\$j5,THROTTLE_CH=4,THROTTLE_IN=4 -o $TO/logs/cluster_%j.out $CODE/02_cluster.sbatch); \
      printf 'proteins\t%s\ngapsmith\t%s\ncluster\t%s\n' \$j1 \$j5 \$j2 > $TO/jobs.tsv; cat $TO/jobs.tsv" ;;
  launch)
    after=${2:-}; dep=${after:+--dependency=afterany:$after}
    ssh ibex "n=\$(( \$(wc -l < $LISTS/bins.tsv) - 1 )); nt=\$(( (n + $GENOMES_PER_TASK - 1) / $GENOMES_PER_TASK )); echo \"\$n genomes in \$nt tasks of $GENOMES_PER_TASK\"; \
      j1=\$(sbatch --parsable $dep --array=0-\$(( nt - 1 ))%${THROTTLE_PR:-100} --export=ALL,OUT=$OUT -o $LOGS/pr_%A_%a.out $CODE/01_proteins.sbatch); \
      j5=\$(sbatch --parsable --array=0-\$(( nt - 1 ))%${THROTTLE_GM:-40} -t 06:00:00 --dependency=aftercorr:\$j1 --export=ALL,OUT=$OUT -o $LOGS/gm_%A_%a.out $CODE/05_gapsmith.sbatch); \
      j2=\$(sbatch --parsable --dependency=afterok:\$j1 --export=ALL,OUT=$OUT,GAPSMITH_JOB=\$j5,THROTTLE_CH=${THROTTLE_CH:-140},THROTTLE_IN=${THROTTLE_IN:-60} -o $LOGS/cluster_%j.out $CODE/02_cluster.sbatch); \
      printf 'proteins\t%s\ngapsmith\t%s\ncluster\t%s\nafter\t%s\n' \$j1 \$j5 \$j2 '$after' > $OUT/jobs.tsv; cat $OUT/jobs.tsv" ;;
  rebalance)
    thr=${2:-600}; which=${3:-gapsmith}
    ssh ibex "j=\$(awk -F'\t' -v w=$which '\$1==w{print \$2}' $OUT/jobs.tsv); scontrol update JobId=\$j ArrayTaskThrottle=$thr && echo \"$which array \$j throttle -> $thr\"" ;;
  compare)
    ssh ibex "module load python/3.11.0; python3 $CODE/compare_models.py --work ${2:-$OUT/test/work} --gapseq-out $GAPSEQ_OUT_HQ" ;;
  status)
    ssh ibex "for j in \$(cut -f2 $OUT/jobs.tsv | grep -E '^[0-9]+$'); do echo \"== \$j \$(sacct -j \$j -X -n -o JobName%12 | head -1)\"; sacct -j \$j -X -P -n --format=State | sort | uniq -c | tr '\n' ' '; echo; done; \
      echo \"DONE_channels \$(ls $WORK/*/DONE_channels 2>/dev/null | wc -l)  DONE_gapseq \$(ls $WORK/*/DONE_gapseq 2>/dev/null | wc -l)  DONE \$(ls $WORK/*/DONE 2>/dev/null | wc -l)  of \$(( \$(wc -l < $LISTS/bins.tsv) - 1 ))\"" ;;
  aggregate)
    ssh ibex "sbatch -p debug -t 02:00:00 -c 2 --mem=8G -J eq-aggregate -o $LOGS/aggregate_%j.out --wrap='module load python/3.11.0; \
      python3 $CODE/04_aggregate.py --bins $LISTS/bins.tsv --work $WORK --out-dir $RES'" ;;
  fetch)
    mkdir -p "$here/results"; scp -q "ibex:$RES/eq_mags_*.tsv" "$here/results/"; ls -la "$here/results" ;;
  *) sed -n '2,12p' "$0" ;;
esac
