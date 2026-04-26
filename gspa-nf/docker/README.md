# GSPA neural-predictor Docker images

Three images cover the four neural predictors plus ensemble + eval:

| Image | Wraps | Base | Size (rough) |
|---|---|---|---|
| `gspa-esm2-stack` | `esm2-deepgoplus`, `esm2-centroid`, `clean` | `pytorch/pytorch:2.4.0-cuda12.1-cudnn9-runtime` | ~6 GB |
| `gspa-proteinfer-stack` | `proteinfer` | `tensorflow/tensorflow:2.15.0` | ~2 GB |
| `gspa-eval-stack` | `ensemble`, `eval` | `python:3.12-slim` | ~150 MB |

## Build

From the repo root (so `benchmark/neural/*` is in the build context):

```bash
docker build -f gspa-nf/docker/Dockerfile.esm2-stack       -t leechuck/gspa-esm2-stack:0.1       .
docker build -f gspa-nf/docker/Dockerfile.proteinfer-stack -t leechuck/gspa-proteinfer-stack:0.1 .
docker build -f gspa-nf/docker/Dockerfile.eval-stack       -t leechuck/gspa-eval-stack:0.1       .
```

## Push

```bash
docker push leechuck/gspa-esm2-stack:0.1
docker push leechuck/gspa-proteinfer-stack:0.1
docker push leechuck/gspa-eval-stack:0.1
```

## Wire into Nextflow

Edit `gspa-nf/nextflow.config`. Replace the placeholder image names in the
`process { withName: ... }` block:

```groovy
withName: 'ESM2_DEEPGOPLUS|ESM2_CENTROID|CLEAN' { container = 'leechuck/gspa-esm2-stack:0.1' }
withName: 'PROTEINFER'                          { container = 'leechuck/gspa-proteinfer-stack:0.1' }
withName: 'ENSEMBLE_PREDS|EVAL_PGAP'            { container = 'leechuck/gspa-eval-stack:0.1' }
```

## GPU runtime

Both stacks work CPU-only. For GPU, run with `-profile docker,gpu` (sets
`containerOptions = '--gpus all'`) or `-profile singularity,gpu` (sets
`runOptions = '--nv'`). See `nextflow.config`.

## Tag conventions

Use semver minor for the wrapped sidecar version:

- `:0.1` first published cut tracking the validated 21-genome benchmark.
- Bump `:0.2` when `run_neural_predictors.py` changes its CLI.
- Bump major (`:1.0`) on schema-breaking changes (column rename, etc.).

## Image size note

The ESM2 stack ships the **runtime** PyTorch/CUDA; weights are NOT baked in
— they're downloaded once via `gspa-nf/fetch_databases.sh` and bind-mounted
by Nextflow. This keeps the image rebuildable and lets you swap weights
without a rebuild.
