#!/bin/bash
# Row-3 certification of the gspa-nf predictor containers on Shaheen: run main.nf's default
# predictor processes (PYRODIGAL, DIAMOND_BLASTP, HMMSEARCH, BARRNAP, MINCED, MERGE_ANNOTATIONS)
# for one real genome (E. coli BW25113) with Singularity images pulled from quay.io by the head
# process on the login node, tasks submitted to Slurm (account k10018, nid00022/nid00497/nid01033 excluded).
# Run detached on a login node: nohup bash nf-shaheen-run.sh > logs/nf-shaheen.log 2>&1 &
# The certification run is a fresh run (every process must appear as COMPLETED in trace.tsv; a
# resumed run lists cached tasks as CACHED). Set NF_RESUME=-resume to resume instead.
set -u
ROOT=/scratch/hohndor/gspa-cert
export JAVA_HOME=/project/k10018/sw/jdk/current
export PATH=$JAVA_HOME/bin:/project/k10018/sw/envs/nf-core/bin:$PATH
export NXF_HOME=$ROOT/nf-home NXF_SINGULARITY_CACHEDIR=$ROOT/nf-cache TMPDIR=$ROOT/tmp
# keep every Singularity layer/build cache off $HOME (its quota broke the first SIF build)
export SINGULARITY_CACHEDIR=$ROOT/sing-cache SINGULARITY_TMPDIR=$ROOT/tmp APPTAINER_CACHEDIR=$ROOT/sing-cache APPTAINER_TMPDIR=$ROOT/tmp
mkdir -p $ROOT/sing-cache
module load singularity
# (exported again after the module: it must reach both the pull below and the Nextflow head)
export SINGULARITY_CACHEDIR=$ROOT/sing-cache SINGULARITY_TMPDIR=$ROOT/tmp APPTAINER_CACHEDIR=$ROOT/sing-cache APPTAINER_TMPDIR=$ROOT/tmp
mkdir -p $ROOT/nf-cache $ROOT/nf-work $ROOT/nf-out $ROOT/logs $ROOT/db/nf
cd $ROOT
# DIAMOND database for the pipeline, built with the pipeline's own diamond container
if [ ! -s $ROOT/db/nf/swissprot_nade_ecn.dmnd ]; then
  DIMG=$ROOT/nf-cache/quay.io-biocontainers-diamond-2.1.9--h43eeafb_0.img
  [ -s $DIMG ] || singularity pull --name $DIMG docker://quay.io/biocontainers/diamond:2.1.9--h43eeafb_0 > $ROOT/logs/nf-diamond-pull.log 2>&1
  singularity exec --bind $ROOT:$ROOT $DIMG \
    diamond makedb --in $ROOT/reference/swissprot_nade_ecn.faa -d $ROOT/db/nf/swissprot_nade_ecn > $ROOT/logs/nf-diamond-makedb.log 2>&1
  echo "diamond makedb rc=$?"
fi
nextflow -version | grep version > $ROOT/nf-out/nextflow_version.txt
singularity --version >> $ROOT/nf-out/nextflow_version.txt
t0=$(date +%s)
nextflow run $ROOT/gspa-nf/main.nf \
  --input $ROOT/bw25113/GCF_000750555.1.fna \
  --outdir $ROOT/nf-out/results \
  --diamond_db $ROOT/db/nf/swissprot_nade_ecn.dmnd \
  --pfam_db /project/k10018/sw/db/pfam-ibex/pfam/Pfam-A.hmm \
  ${NF_RESUME:-} -c $ROOT/nf-shaheen.config -with-trace $ROOT/nf-out/trace.tsv -with-report $ROOT/nf-out/report.html
rc=$?
echo "nextflow rc=$rc seconds=$(( $(date +%s) - t0 ))"
ls -la $ROOT/nf-cache
echo NF_SHAHEEN_DONE
