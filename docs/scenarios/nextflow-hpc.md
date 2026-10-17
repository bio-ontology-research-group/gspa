# Nextflow and HPC

Use this workflow to run GSPA at scale on High Performance Computing (HPC)
clusters or cloud environments using Nextflow DSL2 (`gspa-nf`).

While the standalone JVM CLI (`gspa annotate` / `gspa integrate`) is ideal for
single genomes and local developer workflows, processing large cohorts across a
compute cluster benefits from workflow management: containerized tool dependencies,
automatic job submission to Slurm or PBS, and checkpoint/resume capabilities.

`gspa-nf` orchestrates the complete annotation stack across modular containerized
tasks:

- **Gene calling**: Pyrodigal.
- **Homology and structure**: DIAMOND, MMseqs2, FoldSeek.
- **Protein domains**: HMMER (Pfam-A), InterProScan, eggNOG-mapper, dbCAN.
- **Specialized features**: AMRFinderPlus, antiSMASH, CheckM2, GTDB-Tk, MinCED.
- **Neural predictors**: ESM2-DeepGOPlus, CLEAN, ProteInfer.
- **Evidence integration**: automated claims extraction, operon prediction,
  Bayesian integration, and self-contained HTML browser generation.

## When to use Nextflow vs. JVM CLI

| Scenario | Recommended interface |
|---|---|
| Single isolate or MAG on local workstation | JVM CLI (`gspa annotate` / `gspa integrate`) |
| Batch cohort on a Slurm/PBS cluster | `nextflow run gspa-nf/main.nf` |
| Execution without pre-installing 10+ bioinformatics tools | Nextflow with Docker or Singularity |
| Resuming failed pipeline runs after hardware interruptions | `nextflow run ... -resume` |

## What you need

| Dependency | Minimum version | Purpose |
|---|---|---|
| Nextflow | $\ge$ 22.10.0 | Pipeline workflow manager |
| Singularity or Apptainer | $\ge$ 3.8.0 | HPC container execution (or Docker on local machines) |
| Reference databases | variable | Downloaded via `gspa-nf/fetch_databases.sh` |
| Input data | nucleotide FASTA | Single `.fna` file or a two-column `samplesheet.csv` |

## Step 1: fetch reference databases

`gspa-nf` tracks all external reference files in a versioned manifest
(`gspa-nf/database_manifest.tsv`), specifying exact URLs, release versions,
and SHA-256 checksums.

Use the `fetch_databases.sh` helper to download the required databases:

```bash
# Fetch core databases: DIAMOND UniProt, Pfam-A, and ProteInfer
bash gspa-nf/fetch_databases.sh \
  --outdir /data/reference/gspa_dbs \
  --diamond \
  --pfam \
  --proteinfer
```

To fetch all databases including neural model weights and specialized tools:

```bash
bash gspa-nf/fetch_databases.sh \
  --outdir /data/reference/gspa_dbs \
  --all
```

## Step 2: run locally with Docker

To annotate a single genome on a multi-core desktop or server using Docker:

```bash
nextflow run gspa-nf/main.nf \
  --input genome.fna \
  --diamond_db /data/reference/gspa_dbs/uniprot_sprot.dmnd \
  --pfam_db /data/reference/gspa_dbs/Pfam-A.hmm \
  --outdir results/ \
  -profile docker
```

### Batch execution with a samplesheet

For multiple genomes, create a CSV samplesheet:

```text
sample_id,fasta
isolate_A,/data/genomes/isolate_A.fna
isolate_B,/data/genomes/isolate_B.fna
isolate_C,/data/genomes/isolate_C.fna
```

Run the pipeline:

```bash
nextflow run gspa-nf/main.nf \
  --input samplesheet.csv \
  -c gspa-nf/databases.config \
  --outdir batch_results/ \
  -profile docker
```

## Step 3: run on an HPC cluster with Slurm and Singularity

### 1. Pre-stage container images

On shared cluster filesystems, pulling Singularity images concurrently from
hundreds of compute nodes can saturate network bandwidth. Pre-pull all images
into a shared directory using `prepull_singularity.sh`:

```bash
export NXF_SINGULARITY_CACHEDIR=/data/containers/singularity_cache
bash gspa-nf/prepull_singularity.sh
```

### 2. Configure cluster resources

`gspa-nf/slurm_singularity.config` defines resource allocations and Slurm
partition directives:

```groovy
process {
    executor = 'slurm'
    queue    = 'batch'
    time     = '4h'
    memory   = '16.GB'
    cpus     = 8

    withName: 'DIAMOND_BLASTP' {
        cpus   = 16
        memory = '32.GB'
    }
    withName: 'INTERPROSCAN' {
        cpus   = 16
        memory = '48.GB'
        time   = '8h'
    }
    withName: 'INTEGRATE' {
        memory = '16.GB'
        cpus   = 4
    }
}

singularity {
    enabled    = true
    autoMounts = true
    cacheDir   = '/data/containers/singularity_cache'
}
```

### 3. Launch the Slurm pipeline

Launch the workflow within a `tmux` session or submit it as a lightweight batch job:

```bash
nextflow run gspa-nf/main.nf \
  -c gspa-nf/slurm_singularity.config \
  -c gspa-nf/databases.config \
  --input samplesheet.csv \
  --outdir /scratch/user/gspa_output \
  --run_diamond \
  --run_hmmer \
  --run_interproscan \
  --run_integrate \
  -profile singularity_slurm \
  -resume
```

The `-resume` flag ensures that if any node fails or exceeds its wall-time limit,
re-launching the exact same command skips completed tasks and continues from
the last valid checkpoint.

## Pipeline module organization

`gspa-nf/main.nf` coordinates 22 processes structured into modular sub-workflows:

| Sub-workflow | Modules | Behavior |
|---|---|---|
| `gene_calling.nf` | `PYRODIGAL` | Always executed: calls CDS features and writes translated proteins and GFF3. |
| `similarity.nf` | `DIAMOND_BLASTP`, `MMSEQS2_SEARCH`, `FOLDSEEK` | Controlled by `--run_diamond`, `--run_mmseqs2`, `--run_foldseek`. |
| `domains.nf` | `HMMSEARCH`, `INTERPROSCAN`, `EGGNOG_MAPPER`, `DBCAN` | Controlled by `--run_hmmer`, `--run_interproscan`, `--run_eggnog`, `--run_dbcan`. |
| `specialized.nf` | `AMRFINDER`, `ANTISMASH`, `CHECKM2`, `GTDBTK` | Specialized AMR, secondary metabolite, and quality tools. |
| `neural.nf` | `ESM2_DEEPGOPLUS`, `CLEAN`, `PROTEINFER` | Deep learning models running on CPU or GPU worker nodes. |
| `integrate.nf` | `BUILD_CLAIMS`, `OPERONS`, `INTEGRATE`, `VISUALIZE` | Converts raw tool outputs to `claims.jsonl`, runs Bayesian integration, and compiles `browser.html`. |

When `--run_integrate` is enabled, the pipeline automatically feeds tool outputs
into `gspa-cli integrate`, executing Stage 1 Bayesian integration, Stage 2
symbolic enforcement, and emitting the self-contained HTML browser.

## Output structure

Upon completion, Nextflow writes organized outputs to `--outdir`:

```text
results/
  sample_A/
    proteins.faa                       Pyrodigal protein translations
    annotation.gff3                    Pyrodigal structural coordinates
    diamond_results.tsv                Raw DIAMOND alignments
    pfam_results.domtbl                HMMER Pfam domain table
    interproscan.tsv                   InterProScan annotations
    claims.jsonl                       Consolidated GSPA evidence claims
    operons.tsv                        Predicted operons
    integrated.tsv                     Calibrated posterior probabilities
    provenance.json                    Per-call evidence audit trail
    quality_report.json                GAEF quality metrics
    sample_A_browser.html              Self-contained HTML5 interactive browser
  pipeline_info/
    execution_timeline.html            Wall-clock timing per process
    execution_report.html              CPU and memory usage per task
```

## Pitfalls

- **Shared filesystem lockups**: SQLite-based tools (such as eggNOG-mapper or GTDB-Tk)
  can crash when running concurrently on distributed NFS or GlusterFS filesystems.
  Configure `beforeScript` in `slurm_singularity.config` to copy temporary databases
  to `/tmp` or local node storage.
- **Unset container cache**: if `NXF_SINGULARITY_CACHEDIR` is unset, Nextflow
  downloads Singularity images into the current working directory, rapidly
  exhausting user home directory quotas.
- **Memory limits on InterProScan**: InterProScan requires at least 32 GB to 48 GB
  of RAM when scanning full proteomes. Ensure cluster resource limits match
  tool requirements.
