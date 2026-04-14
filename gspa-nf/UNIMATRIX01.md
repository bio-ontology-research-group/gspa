# Running the GSPA Nextflow pipeline on unimatrix01 (SLURM + Singularity)

This is a worked example of running the pipeline end-to-end on the
unimatrix01 cluster, using the SLURM executor and Singularity
containers (no Docker daemon needed on workers).

## Cluster constraints we hit

- `/data/` is GlusterFS — only mounted on `node005` and `node006`
  (workers `node002`–`node004` and `node007` lack the GlusterFS
  client). The `clusterOptions = '--nodelist=node005,node006'` line in
  `slurm_singularity.config` constrains every job to those two nodes.
- `/home/<user>/` is local SSD on the login node, NOT shared with
  workers — never put the work directory or the Singularity cache
  there. Use `/data/<user>/...`.
- `/storage/` is shared but read-only for normal users.
- Singularity 4.2.2 is available on login node and on `node005` /
  `node006`.
- `python:3.12-slim` lacks `procps` (`ps`), which Nextflow needs for
  task-metric collection. The slurm overlay sets `MERGE_ANNOTATIONS`
  to run on the host (no container) so the worker's `/usr/bin/ps` is
  used.
- `gluster` occasionally surfaces a 0-byte "linkto" stub alongside the
  real file in `ls`. Singularity tolerates it, but if a pull is
  interrupted the *real* file may end up zero-bytes — pre-pull all
  images cleanly before launching Nextflow (`prepull_singularity.sh`).

## One-time setup

```bash
# Java 21 for Nextflow (login node)
conda activate java21
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/

# Pre-pull all Singularity images on the login node
bash gspa-nf/prepull_singularity.sh
```

## Running the M. genitalium reference test

```bash
cd /data/<user>/gspa
bash gspa-nf/run_unimatrix01.sh
```

The launcher script (`run_unimatrix01.sh`) runs:

```bash
nextflow run gspa-nf/main.nf \
  -c gspa-nf/nextflow.config \
  -c gspa-nf/slurm_singularity.config \
  --input /data/<user>/gspa/proteomes/mgenitalium_genomic.fna \
  --diamond_db /data/<user>/gspa/proteomes/reference_loo9.dmnd \
  --pfam_db /storage/software/databases/hmmer/Pfam-A.hmm \
  --outdir results \
  --kingdom bacteria \
  -profile singularity
```

## Verified (2026-04-14)

| Stage              | Container                                            | SLURM jobid | Time |
|--------------------|------------------------------------------------------|-------------|------|
| PYRODIGAL          | quay.io/biocontainers/pyrodigal:3.7.1--py312h247cb63 | 1083        | 2 s  |
| BARRNAP            | quay.io/biocontainers/barrnap:0.8                    | 1084        | 5 s  |
| DIAMOND_BLASTP     | quay.io/biocontainers/diamond:2.1.9                  | 1085        | 36 s |
| MINCED             | quay.io/biocontainers/minced:0.3.0                   | 1086        | 1 s  |
| HMMSEARCH (Pfam)   | quay.io/biocontainers/hmmer:3.4                      | 1087        | 1m54s |
| MERGE_ANNOTATIONS  | (host)                                               | 1088        | 2 s  |

End-to-end wall time: **2m 22s**. Output:
`results/mgenitalium_genomic/`:
- `gene_calling/{*.faa, *.gff, *.fna}` — 995 proteins called
- `diamond/mgenitalium_genomic_diamond.tsv` — DIAMOND hits
- `hmmer/mgenitalium_genomic_pfam.domtbl` — Pfam domain hits
- `barrnap/`, `minced/` — rRNA + CRISPR
- `mgenitalium_genomic_annotations.tsv` — merged predictions
  (3,289 annotations on 657 / 995 proteins; 627 with Pfam)
- `mgenitalium_genomic_annotated.gff3` — annotated GFF
- `mgenitalium_genomic_summary.txt` — counts
