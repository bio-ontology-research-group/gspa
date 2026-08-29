# DeepGO experimental LAFA containers

This directory contains LAFA batch entrypoints for two submitted methods:

- `deepgo-experimental-light-lafa`: a conservative beta predictor using DIAMOND
  homology transfer plus the lightweight DeepGO experimental integrator. All
  runtime assets are derived from LAFA inputs. CPU only.
- `deepgo-experimental-full-lafa`: an experimental CPU/GPU batch predictor using
  the same LAFA-built homology assets plus an ESM2-35M kNN store built from
  `--train_sequences` and cached under `--cache_dir`, and the shipped
  hierarchy-aware 1D-CNN (`models/weights/cnn_mcm.pt`).

Both containers implement the standard LAFA command-line interface and write the
required headerless 3-column output:

```text
Query_ID<TAB>GO_Term<TAB>Score
```

Published images (linux/amd64):

| image | tag | manifest digest | compressed / on disk |
| --- | --- | --- | --- |
| `docker.io/leechuck/deepgo-experimental-light-lafa` | `v1` | `sha256:e1f13f614f26a47aca0279cb03100c66021f843c2fb7cf9222fa7fc033193256` | 85 MB / 240 MB |
| `docker.io/leechuck/deepgo-experimental-full-lafa` | `v1` | `sha256:9b2037491e4c754eca12de5d2b50d54667f02b3a641bb65383447cdab2086199` | 5.3 GB / 8.9 GB |

## Build

Run from the `deepgo-plusplus/` directory:

```bash
docker build -f lafa/Dockerfile.light -t deepgo-experimental-light-lafa:v1 .
docker build -f lafa/Dockerfile.full  -t deepgo-experimental-full-lafa:v1 .
```

The light image is `python:3.12-slim` plus a static DIAMOND 2.1.9 binary; the
light path uses only the Python standard library. The full image is
`nvidia/cuda:12.4.1-cudnn-runtime-ubuntu22.04` (CUDA 12.4.1, cuDNN 9.1.0) with
Ubuntu's Python 3.10 and unpinned `numpy`, `fair-esm` and `torch` from the
`cu124` wheel index; the ESM2-35M weights (`esm2_t12_35M_UR50D`) are downloaded
at build time into `TORCH_HOME=/app/.cache/torch`, so no network access is
needed at run time. The published `v1` images contain
Python 3.10.12, torch 2.6.0+cu124 (compiled for sm_50 to sm_90), fair-esm 2.0.0,
numpy 2.2.6 and triton 3.2.0 (read from the pulled image on 2026-08-29).

## Run

```bash
docker run --rm \
  -v /path/to/lafa_data:/app/data:ro \
  -v /path/to/output:/app/output:rw \
  deepgo-experimental-light-lafa:v1 \
  --query_file /app/data/test_sequences.fasta \
  --train_sequences /app/data/train_sequences.fasta \
  --annot_file /app/data/train_terms.tsv \
  --graph /app/data/go-basic.obo \
  --output_file /app/output/deepgopp_light.tsv.gz \
  --num_threads 8
```

Full model:

```bash
docker run --rm --gpus all \
  -v /path/to/lafa_data:/app/data:ro \
  -v /path/to/output:/app/output:rw \
  deepgo-experimental-full-lafa:v1 \
  --query_file /app/data/test_sequences.fasta \
  --train_sequences /app/data/train_sequences.fasta \
  --annot_file /app/data/train_terms.tsv \
  --graph /app/data/go-basic.obo \
  --output_file /app/output/deepgopp_full.tsv.gz \
  --cache_dir /app/output/dgpp_cache \
  --num_threads 8
```

Without `--gpus all` (or with `--device cpu`) the full model runs entirely on
the CPU; the output is the same up to floating-point differences in the ESM2
embeddings.

## Command-line arguments

`lafa_main.py` is the image ENTRYPOINT, so every argument below is appended
directly to `docker run ... <image>`.

Required:

| argument | meaning |
| --- | --- |
| `--query_file`, `-q` | query FASTA (plain text, not gzipped). The query ID is the first whitespace-delimited token after `>`. The whole file is read into memory and passed to DIAMOND as one job. |
| `--train_sequences` | training FASTA. Used to build the DIAMOND database and, in full mode, the ESM2 embedding store. |
| `--annot_file`, `-a` | training labels: either a CAFA-style `train_terms.tsv` (`EntryID<TAB>term<TAB>aspect`, header optional and auto-detected, plain or `.gz`) or a GAF/GOA file (auto-detected when a data line has at least 9 columns with a `GO:` id in column 5; `NOT` qualifiers are skipped; column 2 must match the FASTA IDs of `--train_sequences`). Propagated files are accepted but unnecessary: the wrapper propagates labels to ancestors itself. |
| `--graph` | GO ontology in OBO format (`go-basic.obo`). `is_a` and `part_of` edges and the namespaces are used; obsolete terms are skipped. |
| `--output_file`, `-o` | output TSV; gzipped when the name ends in `.gz`. The parent directory is created if missing. |

Optional:

| argument | default | meaning |
| --- | --- | --- |
| `--mode {light,full}` | env `DGPP_LAFA_MODE` (light image: `light`, full image: `full`) | which predictor to run. Selecting `full` inside the light image fails because torch is not installed there. |
| `--cache_dir` | env `DGPP_LAFA_CACHE` (`/app/cache`) | where the derived assets are kept: `assets-<hash>/` (DIAMOND database, `train_terms.tsv`, `go-dag.tsv`, `go.obo`, empty `train_net_index.tsv`) and, in full mode, `train_esm2_<hash>.npz`. The hash covers the absolute path, size and mtime of the training inputs (plus the ESM2 model name and layer for the store), so mount a persistent directory and mount the inputs at the same paths to reuse it across runs. Must be writable. |
| `--num_threads` | env `DGPP_THREADS` (`8`) | DIAMOND threads (`-p`) and, in full mode, the torch CPU thread count. |
| `--min_score` | `0.1` | rows with an integrated score below this are not written. |
| `--topk` | `5` | number of best DIAMOND hits per query whose labels vote. |
| `--diamond` | env `DGPP_DIAMOND` (`diamond`) | path to the DIAMOND binary. |
| `--device` | env `DGPP_DEVICE` (`auto`) | full mode only: `auto` (CUDA if `torch.cuda.is_available()`, else CPU), `cuda`, or `cpu`. |
| `--batch_size` | `16` | full mode only: batch size (sequences) for the ESM2 forward passes when building the training store. |
| `--esm2_model` | `esm2_t12_35M_UR50D` | full mode only: `esm.pretrained` model name. Only the default is baked into the image; any other name would try to download weights. |
| `--esm2_layer` | `12` | full mode only: representation layer that is mean-pooled into the embedding. |

Fixed internals (not exposed): DIAMOND runs `blastp --very-sensitive -k 25
--evalue 1e-3`; the ESM2 kNN uses the 10 nearest training proteins by cosine
similarity; sequences are truncated to 1022 residues for ESM2 and 1000 for the
CNN; the three aspect root terms are never written; output rows are sorted by
query ID and then by descending score (the order among equal scores is not
stable between runs).

All wrapper messages go to stderr with the prefix `[lafa-dgpp]`; nothing is
written to stdout. The container runs as root, so files it writes to mounted
directories are owned by root; `--user $(id -u):$(id -g)` works as long as
`--cache_dir` points to a directory that user can write.

## Tested commands and measured resources

Measured 2026-08-29 on a laptop with an Intel Core Ultra 9 285H (16 logical
CPUs), 62 GB RAM, an NVIDIA RTX PRO 2000 Blackwell GPU (8 GB), Debian 13,
Docker 26.1.5, NVIDIA driver 610.43. Input bundle: `train_sequences.fasta` with
3,000 SwissProt proteins and their CAFA6 `train_terms.tsv` labels (21,437
rows), `go-basic.obo` (release 2025-06-01), and held-out SwissProt queries
(40 for the small test, 2,000 for the large one). Wall time is the whole
`docker run`, including container start and the asset build; peak memory is
from `docker stats`.

Light image (exact command as in "Run" above, with the `small` and `large`
query files; the large run added `--cache_dir /app/output/dgpp_cache`):

| queries | threads | peak RAM | cache dir | wall time | rows written |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 40 | 8 | 49 MiB | 44 MB | 4.7 s | 776 (9 queries with a DIAMOND hit) |
| 2,000 | 8 | 125 MiB | 44 MB | 4.7 s | 58,386 (409 queries with a hit) |

Full image:

| queries | device | ESM2 store | threads | peak RAM | cache dir | wall time | rows written |
| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 40 | `--gpus all`, `--device auto` | building | 8 | 376 MiB | 44 MB | 6.9 s, exit 1 | none: `RuntimeError: CUDA error: no kernel image is available for execution on the device` (Blackwell sm_120 is not in the cu124 build's `sm_50` to `sm_90` list; see below) |
| 40 | `--device cpu` | built in this run (20.5 min for 3,000 proteins) | 8 | 4.19 GiB | 49 MB (store 5.3 MB) | 20 min 45 s | 114,068 (all 40 queries) |
| 2,000 | `--device cpu` | cached | 8 | 4.76 GiB | 49 MB | 15 min 55 s | 5,731,940 (all 2,000 queries; 28 MB gzipped) |

For reference, the same `predict_full` code path in the GSPA service image with a
precomputed 82,404-protein store scored 256 proteins in 16.2 s on a GPU and
148.4 s on CPU (commit 4ccf4f3).

Guidance derived from these runs:

- Light: pass the whole query set in one call (no internal batching); 2 GB RAM
  and 1 GB scratch are ample; minutes of wall time even for tens of thousands
  of queries. The cache directory is dominated by the copied OBO and the GO
  closure table (about 42 MB); the DIAMOND database adds roughly 0.6 bytes per
  training residue.
- Full: the ESM2 store of the training set is the dominant one-off cost, about
  0.4 s per training protein on 8 CPU threads (roughly 9 hours for an 82,000
  protein training set) versus minutes on a supported GPU; build it once into a
  persistent `--cache_dir` (ideally with a GPU) and reuse it. With the store
  cached the query side costs about 0.5 s per query on CPU, so keep CPU batches
  at a few thousand queries; 8 GB RAM and 2 GB scratch are sufficient. Every
  query receives rows (ESM2 kNN and CNN score any sequence), about 2,900 rows
  per query at `--min_score 0.1`; raise `--min_score` to reduce output volume.
- GPU support: `torch.cuda.get_arch_list()` in the image is `sm_50, sm_60,
  sm_70, sm_75, sm_80, sm_86, sm_90` (Maxwell to Hopper, e.g. A100, H100).
  Blackwell GPUs (sm_120) are detected as available but have no kernels, so
  `--device auto` crashes; use `--device cpu` there or rebuild the image with a
  cu128 torch wheel. The host driver must support CUDA 12.4 (driver 550 or
  newer).

## Provenance and LAFA input policy

The wrapper follows the strict LAFA-input build mode:

- builds `train_db.dmnd` from `--train_sequences`;
- converts `--annot_file` from either train-terms TSV or GAF/GOA to
  `train_terms.tsv`;
- builds the GO ancestor closure from `--graph`;
- creates an empty `train_net_index.tsv`, because LAFA's standard data bundle does
  not provide STRING. This makes the `net_union` component deterministically zero
  instead of using outside-release evidence.

The full image ships the learned DG++ weights and builds the ESM2 kNN reference
embeddings from the LAFA training FASTA on first run. Mount a writable cache/output
directory so subsequent LAFA runs for the same release reuse those embeddings.

Frozen parts and the releases they were trained on (see `../VERSIONS.md` and
`../TRAINING.md`):

| artefact | used by | training inputs |
| --- | --- | --- |
| `models/deepgo_plusplus_light_fast.json` (integrator over `diam`, `net_union`) | light | CAFA6 official `train_terms.tsv` (Kaggle release of 2025-12-10, 82,404 proteins, cutoff t0 = 2026-02-02); SwissProt sequences from UniProt release 2025_03 (`uniprot_sprot.dat.gz`, entry versions up to 18-JUN-2025); STRING v12.0 (2023) for the `net_union` training feature; fitting targets = experimental annotations gained after t0 in `goa_uniprot_all.gaf.gz` generated 2026-06-17 (GAF 2.2, `!go-version` 2026-06-15); GO `go-basic.obo` release 2025-10-10 |
| `models/deepgo_plusplus_integrator_cpu_lean_mcm.json` (integrator over `diam`, `interpro`, `cnn`, `net_union`, `esm2_knn`, `proteinfer`, `esm2_head`; fit 2026-06-24) | full | same inputs as above |
| `models/weights/cnn_mcm.pt` (1D-CNN, C-HMCNN loss, 5,265 terms) | full | pre-t0 SwissProt sequences (UniProt 2025_03) of the CAFA6 training proteins with `train_terms.tsv` labels propagated over `go-basic.obo` 2025-10-10 |
| `esm2_t12_35M_UR50D` (fair-esm) | full | Meta AI pretraining on UniRef50, unchanged |

In the LAFA containers only `diam` (light) and `diam`, `esm2_knn`, `cnn` (full)
carry live evidence; `net_union` is zero by construction and the `interpro`,
`proteinfer` and `esm2_head` features of the full integrator are not computed,
so the absolute calibration differs from the offline benchmark numbers in
`../VERSIONS.md`.
