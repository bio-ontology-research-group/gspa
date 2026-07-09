# DeepGO-PlusPlus LAFA containers

This directory contains LAFA batch entrypoints for two submitted methods:

- `deepgopp-light-lafa`: DIAMOND homology transfer plus the DG++Light learned
  integrator. All runtime assets are derived from LAFA inputs.
- `deepgopp-full-lafa`: DG++ full CPU/GPU batch predictor using the same LAFA-built
  homology assets plus an ESM2-35M kNN store built from `--train_sequences` and
  cached under `--cache_dir`.

Both containers implement the standard LAFA command-line interface and write the
required headerless 3-column output:

```text
Query_ID<TAB>GO_Term<TAB>Score
```

## Build

Run from the `deepgo-plusplus/` directory:

```bash
docker build -f lafa/Dockerfile.light -t deepgopp-light-lafa:v1 .
docker build -f lafa/Dockerfile.full  -t deepgopp-full-lafa:v1 .
```

## Run

```bash
docker run --rm \
  -v /path/to/lafa_data:/app/data:ro \
  -v /path/to/output:/app/output:rw \
  deepgopp-light-lafa:v1 \
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
  deepgopp-full-lafa:v1 \
  --query_file /app/data/test_sequences.fasta \
  --train_sequences /app/data/train_sequences.fasta \
  --annot_file /app/data/train_terms.tsv \
  --graph /app/data/go-basic.obo \
  --output_file /app/output/deepgopp_full.tsv.gz \
  --cache_dir /app/output/dgpp_cache \
  --num_threads 8
```

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
