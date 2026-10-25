# Nextflow parameters

The `gspa-nf` directory contains an enterprise Nextflow DSL2 pipeline implementing
the full GSPA workflow across high-performance computing clusters and cloud environments.
This page documents all pipeline parameters (`params.*`), execution profiles,
database manifest specifications, and helper staging scripts.

## Pipeline execution overview

The pipeline is invoked with `nextflow run`:

```bash
# Run with Docker on a local workstation
nextflow run gspa-nf/main.nf \
  -c gspa-nf/databases.config \
  --input genome.fna \
  -profile docker

# Run with Singularity and SLURM on an HPC cluster
nextflow run gspa-nf/main.nf \
  -c gspa-nf/databases.config \
  -c gspa-nf/slurm_singularity.config \
  --input samplesheet.csv \
  -profile singularity,gpu
```

The primary input `--input` accepts either a single FASTA file (`.fna`, `.fasta`, `.faa`)
or a two-column CSV samplesheet (`sample_id,fasta`).

## Pipeline parameters (`params.*`)

### Input and organism options

| Parameter | Type | Default | Description |
|---|---|---|---|
| `input` | string | `null` | Path to query nucleotide FASTA file or samplesheet CSV. |
| `kingdom` | string | `bacteria` | Organism taxonomic kingdom: `bacteria`, `archaea`, `eukaryote`, `virus`, or `auto`. |
| `mag` | boolean | `false` | Enable metagenome-assembled genome mode (relaxes completeness criteria). |
| `outdir` | string | `results` | Output directory where result tables and annotations are saved. |

### Core alignment and domain databases

| Parameter | Type | Default | Description |
|---|---|---|---|
| `diamond_db` | file | `null` | DIAMOND database (`.dmnd`) compiled from Swiss-Prot or UniRef90. |
| `pfam_db` | file | `null` | Pfam-A HMM database (`Pfam-A.hmm`), required to be formatted with `hmmpress`. |
| `interproscan_dir` | directory | `null` | InterProScan installation root containing `interproscan.sh` and `data/`. |
| `foldseek_db` | directory | `null` | FoldSeek tertiary structure database (e.g. AlphaFold DB Swiss-Prot). |
| `prostt5_model` | directory | `null` | ProstT5 model for sequence-to-structure search without precomputed structures. |
| `mmseqs2_db` | directory | `null` | MMseqs2 target sequence database (alternative to DIAMOND). |
| `eggnog_db` | directory | `null` | eggNOG-mapper diamond and taxon database directory. |
| `amrfinder_db` | directory | `null` | NCBI AMRFinderPlus database directory. |
| `antismash_db` | directory | `null` | antiSMASH cluster detection database directory. |
| `dbcan_db` | directory | `null` | dbCAN CAZyme HMM and diamond database directory. |
| `gtdbtk_db` | directory | `null` | GTDB-Tk reference taxonomy database (`GTDBTK_DATA_PATH`). |
| `checkm2_db` | file | `null` | CheckM2 DIAMOND database (`uniref100.KO.1.dmnd`). |
| `signalp_path` | file | `null` | Path to licensed SignalP v6 binary executable. |

### Quality assessment and ground truth evaluation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `go_owl` | file | `null` | Gene Ontology OWL or OBO file for completeness, coherence, and consistency. |
| `go_obo` | file | `null` | Gene Ontology OBO file used for namespace aspect lookup in evaluation. |
| `truth_dir` | directory | `null` | Directory containing `<sample_id>_truth.tsv` curated reference files. |
| `eval_truth_name` | string | `truth` | Label identifier stamped on ground truth comparisons. |
| `eval_n_bootstrap` | integer | `0` | Number of bootstrap resamples for confidence intervals (0 = point estimates). |
| `report_min_score` | float | `0.0` | Minimum score threshold for rows included in the multi-format report. |

### Neural functional predictors

| Parameter | Type | Default | Description |
|---|---|---|---|
| `esm2_dgp_ckpt` | file | `null` | Trained ESM-2 DeepGO-Plus classification head weights (`.pt`). |
| `esm2_dgp_terms` | file | `null` | GO term vocabulary mapping matching the classifier head output layer. |
| `esm2_dgp_model` | string | `esm2_t33_650M_UR50D` | Pretrained ESM-2 backbone transformer architecture. |
| `esm2_dgp_batch_size` | integer | `16` | Batch size for ESM-2 DeepGO-Plus sequence inference. |
| `esm2_centroid_db` | file | `null` | NPZ archive with precomputed functional centroids. |
| `esm2_centroid_model` | string | `esm2_t12_35M_UR50D` | Lightweight ESM-2 model used for centroid query embedding. |
| `esm2_centroid_top_k` | integer | `5` | Nearest functional centroids retained per query protein. |
| `esm2_centroid_batch_size` | integer | `32` | Batch size for centroid embedding extraction. |
| `clean_model_dir` | directory | `null` | CLEAN checkpoint directory containing model weights and distance maps. |
| `clean_batch_size` | integer | `16` | Batch size for CLEAN contrastive inference. |
| `proteinfer_model_dir` | directory | `null` | ProteInfer TensorFlow SavedModel directory. |
| `proteinfer_batch_size` | integer | `16` | Batch size for ProteInfer CNN prediction. |
| `neural_min_score` | float | `0.1` | Minimum prediction score emitted by neural sidecar runners. |

### Ensemble fusion

| Parameter | Type | Default | Description |
|---|---|---|---|
| `ensemble_mode` | string | `mean` | Score aggregation across per-protein predictors: `max`, `mean`, or `rank`. |
| `ensemble_min_score` | float | `0.0` | Minimum consensus score retained in ensemble outputs. |
| `genomic_ensemble_mode` | string | `max` | Region aggregation mode: `max`, `mean`, or `wcov` (weighted coverage). |
| `genomic_ensemble_min_overlap` | float | `0.5` | Minimum reciprocal overlap fraction required to cluster genomic calls. |
| `genomic_ensemble_min_score` | float | `0.0` | Minimum score cutoff for genomic ensemble calls. |
| `run_ensemble_genomic` | boolean | `false` | Enable multi-caller fusion across viral and prophage predictors. |

### Open-source auxiliary predictors

| Parameter | Type | Default | Description |
|---|---|---|---|
| `metapredict_min_score` | float | `0.5` | Minimum disorder score threshold for Metapredict. |
| `metapredict_min_region_len` | integer | `10` | Minimum amino acid length for intrinsically disordered segments. |
| `deepsig_min_score` | float | `0.5` | Minimum confidence score for DeepSig signal peptide calls. |
| `deepsig_kingdom` | string | `gramn` | DeepSig model kingdom: `gramp`, `gramn`, or `euk`. |
| `tmbed_min_score` | float | `0.5` | Minimum confidence score for TMbed transmembrane predictions. |
| `tmbed_prott5_dir` | directory | `null` | Path to ProtT5 weights (~5 GB). TMbed downloads if unset. |
| `tppred3_min_score` | float | `0.5` | Minimum confidence score for TPpred3 targeting peptide calls. |
| `tppred3_kingdom` | string | `nonplant` | TPpred3 model kingdom: `plant` or `nonplant`. |
| `psortb_min_score` | float | `0.7` | Minimum score threshold for PSORTb localization predictions. |
| `psortb_gram` | string | `negative` | PSORTb target organism type: `positive`, `negative`, or `archaea`. |
| `deepfri_model_dir` | directory | `null` | DeepFRI model repository containing pretrained weights. |
| `deepfri_min_score` | float | `0.5` | Minimum score threshold for DeepFRI predictions. |
| `deepec_model_dir` | directory | `null` | DeepEC repository directory. |
| `deepec_min_score` | float | `0.5` | Minimum score threshold for DeepEC predictions. |
| `deeparg_model_dir` | directory | `null` | DeepARG model and database directory. |
| `deeparg_min_score` | float | `0.5` | Minimum score threshold for DeepARG predictions. |
| `deeparg_type` | string | `prot` | Query sequence type for DeepARG: `prot` or `nucl`. |
| `mdf_weights` | directory | `null` | mDeepFRI weights directory. |
| `mdf_skip_pdb` | boolean | `true` | When true, runs mDeepFRI sequence-only path without FoldComp DB. |
| `mdf_foldcomp_db` | directory | `null` | FoldComp database directory for mDeepFRI structural path. |
| `mdf_min_score` | float | `0.1` | Minimum score threshold for mDeepFRI predictions. |
| `musitedeep_model_dir` | directory | `null` | MusiteDeep model repository directory. |
| `musitedeep_residue_types` | string | `Phosphoserine_Phosphothreonine` | Target post-translational modification residue types. |
| `musitedeep_min_score` | float | `0.5` | Minimum score threshold for MusiteDeep PTM site calls. |
| `scannet_model_dir` | directory | `null` | ScanNet model weights directory. |
| `scannet_min_score` | float | `0.5` | Minimum score threshold for ScanNet binding site predictions. |
| `structures_from` | string | `none` | Source of tertiary protein structures: `none`, `esmfold`, or `afdb`. |
| `genomad_db` | directory | `null` | geNomad viral and plasmid model database path. |
| `genomad_min_score` | float | `0.7` | Minimum virus/plasmid score threshold for geNomad. |
| `checkv_db` | directory | `null` | CheckV provirus reference database directory. |
| `checkv_min_score` | float | `0.5` | Minimum completeness score for CheckV provirus calls. |
| `phispy_input` | file | `null` | Optional GenBank file with CDS features for PhiSpy prophage calling. |
| `phispy_min_score` | float | `0.5` | Minimum score threshold for PhiSpy prophage calls. |
| `virsorter2_db` | directory | `null` | VirSorter2 database directory (~11 GB). |
| `virsorter2_min_score` | float | `0.5` | Minimum score threshold for VirSorter2 viral calls. |
| `vibrant_db` | directory | `null` | VIBRANT HMM model directory. |
| `vibrant_min_score` | float | `0.5` | Minimum score threshold for VIBRANT virus annotations. |

### Process execution switches (`run_*`)

Individual tool processes are toggled using boolean flags. By default, core tools
are enabled while specialized, neural, and viral predictors are opt-in:

- Core tools: `run_diamond` (true), `run_hmmer` (true), `run_barrnap` (true), `run_minced` (true).
- Domain and structure: `run_interproscan` (false), `run_foldseek` (false), `run_mmseqs2` (false), `run_eggnog` (false).
- Specialized: `run_amrfinder` (false), `run_antismash` (false), `run_dbcan` (false), `run_gtdbtk` (false), `run_checkm2` (false), `run_signalp` (false).
- Neural models: `run_esm2_deepgoplus` (false), `run_esm2_centroid` (false), `run_clean` (false), `run_proteinfer` (false), `run_ensemble` (false), `run_eval` (false), `run_report` (false).
- Region and site tools: `run_metapredict` (false), `run_deepsig` (false), `run_tmbed` (false), `run_tppred3` (false), `run_psortb` (false), `run_deepfri` (false), `run_deepec` (false), `run_deeparg` (false), `run_mdf` (false), `run_musitedeep` (false), `run_scannet` (false).
- Viral predictors: `run_genomad` (false), `run_checkv` (false), `run_phispy` (false), `run_flatten_prophage_cds` (false), `run_virsorter2` (false), `run_vibrant` (false).

### Bayesian evidence integration parameters (`run_integrate = true`)

When `run_integrate` is enabled, the pipeline executes `BUILD_CLAIMS` and `INTEGRATE`
to combine all evidence into unified posterior annotations:

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gspa_jar` | file | `null` | Path to compiled `gspa-2.2.2.jar`. |
| `goa` | file | `null` | Reference Swiss-Prot GOA annotation file (`.gaf` or `.gaf.gz`). |
| `pfam2go` | file | `null` | External Pfam-to-GO mapping file. |
| `ec2go` | file | `null` | External EC-to-GO mapping file. |
| `pathways` | file | `null` | KEGG metabolic pathway TSV mapping reactions to pathways. |
| `theta_file` | file | `null` | Custom hyperparameters JSON file (`theta.json`). |
| `essential_profile` | string | `bacteria` | Completeness essentiality profile: `bacteria`, `archaea`, `eukaryote`. |
| `enable_priors` | string | `essentiality,coherence,gap_filling,genomic_context` | Comma-delimited list of active genome-scale priors. |

### Computational resource constraints

| Parameter | Default | Description |
|---|---|---|
| `max_cpus` | Host CPU count | Maximum CPUs available to any single task. |
| `max_memory` | `32.GB` | Maximum RAM allocated to any single task. |
| `max_time` | `24.h` | Maximum execution time before a task times out. |

## Execution profiles

Profiles are selected on the command line using `-profile <name>[,<name>...]`:

### `standard` and `docker`
Enables Docker container execution. Mounts the current user ID and group ID
(`-u $(id -u):$(id -g)`) to ensure output files retain proper host filesystem
permissions.

### `singularity`
Enables Singularity container execution with automatic filesystem mounting
(`singularity.autoMounts = true`). Suitable for shared multi-tenant clusters where
Docker daemon access is restricted.

### `local`
Runs all tools directly on the host machine without container isolation. Requires
all binaries (DIAMOND, HMMER, Pyrodigal) to be installed on system `PATH`.

### `gpu`
Can be combined with container profiles (e.g. `-profile docker,gpu` or
`-profile singularity,gpu`). Passes GPU devices into containerized deep learning
processes (`--gpus all` for Docker, `--nv` for Singularity). Applies to
`ESM2_DEEPGOPLUS`, `ESM2_CENTROID`, `CLEAN`, and `PROTEINFER`.

### `test`
Self-contained integration test profile. Configures inputs and reference databases
to annotate the bundled *Mycoplasma genitalium* test genome (`m_genitalium.fna`)
with DIAMOND and Pfam.

### `smoke_proteinfer`
Lightweight CPU-only integration test. Disables DIAMOND and Pfam to run ProteInfer,
ensemble fusion, evaluation, and reporting on a synthetic test fixture.

### `slurm_singularity.config`
Overlay configuration for SLURM workload managers:
- Sets `process.executor = 'slurm'`.
- Configures default partitions (`debug`), memory, and wall-clock times.
- Defines per-process resource reservations (e.g. 8 CPUs and 8 GB for `DIAMOND_BLASTP`,
  8 CPUs and 12 GB for `HMMSEARCH`).
- Implements transparent task retry on exit codes 1 to 255 to recover from transient
  network filesystem stale-cache errors.
- Sets a Singularity image cache directory and executes a pre-warm read loop
  before task execution.

## Database manifest (`database_manifest.tsv`)

The versioned manifest `gspa-nf/database_manifest.tsv` tracks all neural model
checkpoints, weights, and truth sets required by optional pipeline modules:

| Predictor | Artefact | Version | Default URL base | Size | Notes |
|---|---|---|---|---|---|
| `esm2-deepgoplus` | `head_ckpt` | v1.0 | `${GSPA_DB_BASE}/esm2_dgp/v1.0/head.pt` | ~5.0 GB | Frozen ESM2-t33 + FC head. GPU recommended. |
| `esm2-deepgoplus` | `go_terms` | v1.0 | `${GSPA_DB_BASE}/esm2_dgp/v1.0/terms.txt` | ~200 KB | One GO term per line, matching head columns. |
| `esm2-centroid` | `centroid_db` | v1.0 | `${GSPA_DB_BASE}/esm2_centroid/v1.0/db.npz` | ~75 MB | NPZ with centroids, terms, and types. |
| `proteinfer` | `saved_model` | v1.0 | `${GSPA_DB_BASE}/proteinfer/v1.0/model.tar.gz` | ~200 MB | TensorFlow SavedModel format. |
| `clean` | `weights` | upstream | Upstream GitHub archive | ~1.0 GB | Upstream CLEAN checkpoint repository. |
| `truth` | `swissprot_go_ec` | 2026.04 | `${GSPA_DB_BASE}/truth/2026.04/swissprot_go_ec.tsv.gz` | ~50 MB | Accession to GO and EC mapping for evaluation. |
| `truth` | `smoke_fixture` | 2026.04 | `${GSPA_DB_BASE}/truth/2026.04/smoke_fixture.tsv` | ~5 KB | Ground truth for smoke testing. |
| `go_aspect_map` | `go_obo` | 2026.04 | `${GSPA_DB_BASE}/go/2026.04/go.obo` | ~150 MB | Gene Ontology definition file for GAEF eval. |
| `tmbed` | `prott5` | v1.0 | `${GSPA_DB_BASE}/tmbed/v1.0/prott5.tar.gz` | ~5.5 GB | ProtT5 transformer weights. |
| `deepfri` | `model` | v1.0 | `${GSPA_DB_BASE}/deepfri/v1.0/DeepFRI.tar.gz` | ~200 MB | DeepFRI GCN model weights. |
| `deepec` | `model` | v1.0 | `${GSPA_DB_BASE}/deepec/v1.0/deepec.tar.gz` | ~50 MB | DeepEC CNN enzyme classifier. |
| `deeparg` | `model` | v1.0 | Upstream Virginia Tech mirror | ~100 MB | DeepARG antimicrobial resistance database. |
| `musitedeep` | `model` | v1.0 | `${GSPA_DB_BASE}/musitedeep/v1.0/MusiteDeep_web.tar.gz` | ~50 MB | Phosphorylation prediction models. |
| `tppred3` | `model` | v1.0 | `${GSPA_DB_BASE}/tppred3/v1.0/tppred3.tar.gz` | ~10 MB | Targeting peptide prediction models. |
| `scannet` | `model` | v1.0 | `${GSPA_DB_BASE}/scannet/v1.0/ScanNet_models.tar.gz` | ~50 MB | Geometric deep learning models for PDB structures. |
| `psortb` | `bundled` | v3.0.6 | Container image | N/A | Pre-bundled inside Docker container image. |

The default base URL is `https://gspa.bio2vec.net/db`, customizable via the
`GSPA_DB_BASE` environment variable.

## Helper scripts

### `fetch_databases.sh`
Downloads and unpacks model weights based on `database_manifest.tsv`:

```bash
# Download and extract ProteInfer and ground truth evaluation data
bash gspa-nf/fetch_databases.sh \
  --predictor proteinfer \
  --predictor truth \
  --dest ~/gspa-db
```

Features:
- Resumable file transfers using `curl -C -`.
- Automatic SHA-256 integrity verification.
- Automated extraction of `.tar.gz` and `.gz` archives.
- Outputs corresponding `databases.config` parameter export snippets to stdout.

### `prepull_singularity.sh`
Pre-populates a local or shared Singularity cache directory with container images
specified in `nextflow.config`:

```bash
export SINGULARITY_CACHEDIR=/shared/singularity_cache
bash gspa-nf/prepull_singularity.sh
```

Pre-pulls biocontainers for Pyrodigal, DIAMOND, HMMER, Barrnap, MinCED, and
Python base environments, naming them according to Nextflow's internal cache
conventions.

### `stage_images.sh`
Stages Singularity container images cleanly through local storage:

```bash
bash gspa-nf/stage_images.sh
```

Pulls images to a local `/tmp` scratch location before copying to shared storage,
verifying image readability using `singularity inspect`. This prevents file
truncation and corruption issues commonly encountered when multiple cluster
worker nodes pull images simultaneously across network file systems.
