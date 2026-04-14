# Running the GSPA Nextflow pipeline on unimatrix01 (SLURM + Singularity)

This is a worked example of running the pipeline end-to-end on the
unimatrix01 cluster, using the SLURM executor and Singularity
containers (no Docker daemon needed on workers).

## Cluster constraints we hit

- `/data/` is GlusterFS, mounted on all worker nodes (`node002`-`node007`).
  Singularity 4.2.2 is available on all workers. No nodelist restriction
  is needed.
- `/home/<user>/` is local SSD on the login node, NOT shared with
  workers — never put the work directory or the Singularity cache
  there. Use `/data/<user>/...`.
- `/storage/` is shared but read-only for normal users.
- `python:3.12-slim` lacks `procps` (`ps`), which Nextflow needs for
  task-metric collection. Use full `python:3.12` instead — its biggest
  downside is the 367MB image size vs 41MB slim. `MERGE_ANNOTATIONS`
  runs in this full python container so every pipeline stage is in
  Singularity.
- GlusterFS sometimes surfaces duplicate directory entries (one real
  file + one 0-byte "linkto" stub) after writes; Singularity resolves
  this fine in practice, but *interrupted* pulls can leave a 0-byte
  real file. Pull images fresh into `/tmp` and then `cp` to the
  shared cache — the `stage_images.sh` script automates that. Run
  once before launching Nextflow, and re-run if any image fails to
  `singularity inspect`.
- GlusterFS read inconsistency: a worker occasionally sees a stale /
  corrupt version of a shared-cache file even though the login node
  sees a valid one. The `beforeScript` in `slurm_singularity.config`
  does a full `cat "$img" > /dev/null` pre-read of every cached image
  to force a fresh FUSE fetch before Singularity opens it. Combined
  with `errorStrategy = 'retry'` (maxRetries = 2), the occasional
  stale read is absorbed without human intervention.
- `node007` on the current cluster is stuck in SLURM `COMPLETING`
  state (slurmd epilog stall). The config sets
  `clusterOptions = '--exclude=node007'` to route jobs to healthy
  nodes. Remove this when node007 recovers.

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

All 6 stages containerized via Singularity, distributed across SLURM
workers. `node007` excluded due to stuck COMPLETING state.

| Stage              | Container                                            | Time  |
|--------------------|------------------------------------------------------|-------|
| PYRODIGAL          | quay.io/biocontainers/pyrodigal:3.7.1--py312h247cb63 | 22 s  |
| BARRNAP            | quay.io/biocontainers/barrnap:0.8                    | ~10 s |
| DIAMOND_BLASTP     | quay.io/biocontainers/diamond:2.1.9                  | ~40 s |
| MINCED             | quay.io/biocontainers/minced:0.3.0                   | ~5 s  |
| HMMSEARCH (Pfam)   | quay.io/biocontainers/hmmer:3.4                      | ~3m   |
| MERGE_ANNOTATIONS  | python:3.12 (with procps)                            | ~5 s  |

End-to-end wall time: **3m 53s** (vs. 2m 22s when MERGE ran on host —
the extra ~90s is the per-container pre-warm + python:3.12 image load).

Output:
`results/mgenitalium_genomic/`:
- `gene_calling/{*.faa, *.gff, *.fna}` — 995 proteins called
- `diamond/mgenitalium_genomic_diamond.tsv` — DIAMOND hits
- `hmmer/mgenitalium_genomic_pfam.domtbl` — Pfam domain hits
- `barrnap/`, `minced/` — rRNA + CRISPR
- `mgenitalium_genomic_annotations.tsv` — merged predictions
  (3,289 annotations on 657 / 995 proteins; 627 with Pfam)
- `mgenitalium_genomic_annotated.gff3` — annotated GFF
- `mgenitalium_genomic_summary.txt` — counts
