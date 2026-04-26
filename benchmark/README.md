# `benchmark/` — evaluation harness + neural sidecar

Python-side tooling that sits outside the Gradle build. Two audiences:

1. **Scorers** — `benchmark_pgap_v2.py` and friends produce F-max / Smin /
   EC-per-level numbers used by `benchmark/RESULTS.md`.
2. **Neural sidecars** — `benchmark/neural/run_neural_predictors.py` is the
   Python entry point that the JVM-side `DeepGoPlusEsm2Predictor` /
   `ProteInferPredictor` / `CleanPredictor` / `Esm2CentroidPredictor` wrap
   over subprocess.

## Directory layout

| Path | What lives there |
|---|---|
| `benchmark/*.py` | Shared scorers and ad-hoc scripts: `benchmark_pgap_v2.py` (F-max / Smin / EC), `extract_goa_dual.py` (GOA → truth TSV), `01_..` through `04_..` (older pipeline steps). These are the **top-level utilities**; anything in a subdirectory can assume they live here. |
| `benchmark/neural/` | Neural predictor sidecar + DB builders: `run_neural_predictors.py` (the Python entry invoked by JVM predictors), `build_esm2_centroids.py`, `build_foldseek_centroids.py`, `propagate_truth.py`, `evaluate_panel.py`, plus SLURM `*.sbatch` wrappers for panel runs. |
| `benchmark/panel/` | 21-genome reference panel: `panel_manifest.tsv`, `download_panel.sh`, `build_refseq_uniprot_map.py`, per-genome annotation / integrate drivers. |
| `benchmark/cultures/` | MAG / culture-genome-specific tooling (Phase 12 panel expansion). |
| `benchmark/cross_genome/` | Cross-genome LRO / prior-learning scripts. |
| `benchmark/ml/` | Prior-weight learning harness. |
| `benchmark/leave_reaction_out/` | LRO harness (gapsmith). |

## Path conventions for scripts

- Shared utilities (e.g. `extract_goa_dual.py`, `benchmark_pgap_v2.py`) stay
  at `benchmark/` root and are called by absolute path. **Do not use
  `$SCRIPTS/../extract_goa_dual.py` style traversal from subdirectory
  scripts** — that form has silently bitten SLURM jobs twice because the
  deploy location on the cluster has a different relative layout.
- New SLURM sbatch wrappers should reference tooling by absolute path
  (e.g. `/data/hohndor/gspa-neural/benchmark/neural/run_neural_predictors.py`),
  not relative-traversal. A pinned `BENCH=/data/...` variable at the top of
  the sbatch is idiomatic.

## The neural sidecar protocol

`run_neural_predictors.py --predictor {esm2-deepgoplus|proteinfer|clean|esm2-centroid}`
reads a manifest TSV with columns `tag`, `fasta_path`, `output_dir` and
emits one TSV per genome:

```
protein_id<TAB>term<TAB>score<TAB>annotation_type
```

This is the contract the JVM predictors parse. Adding a new neural
predictor is a matter of adding a `run_<name>(...)` function + a CLI case
— no new sidecar script.

## F-max / Smin definitions

See [`RESULTS.md`](RESULTS.md) for the exact CAFA III/IV protocols and
why both a per-genome micro F-max and a protein-centric CAFA F-max are
reported side by side.
