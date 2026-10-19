# P08 outputs

## On ibex (`/ibex/scratch/projects/c2014/hohndor/price149_repro_2026-10-19/`)

| File | Content |
| --- | --- |
| `out/parts/<pred>.pkl.gz`, `out/parts/<pred>.sha256.tsv` | Per-predictor extraction (array task) and input checksums |
| `out/price149_all6_compact.pkl.gz` | Merged compact pickle (schema in README.md) |
| `out/input_sha256.tsv` | sha256 / size / mtime of the 132 raw pickles + manifest (copied here as `input_sha256.tsv`) |
| `logs/p08_extract_<array>_<task>.out`, `logs/p08_merge_<job>.out` | Slurm logs |

## In the repository

| File | Written by | Content |
| --- | --- | --- |
| `benchmark/manifests/reference_data/ibex_extracted/price149_all6_compact.pkl.gz` | `run.sh fetch` | v2 compact pickle (`_meta.extraction_version = v2-2026-10-19`) |
| `benchmark/manifests/reference_data/ibex_extracted/price149_all6_compact.pkl.gz.old` | `run.sh fetch` (first time only) | v1 pickle (2026-10-15, `scripts/extract_ibex_compact_matrices.py`), kept for comparison |
| `results/pillar2_price149_all6_real.tsv` | `scripts/eval_price149_real_all6.py` | One row per (baseline, arm); 6 × 13 rows |
| `results/pillar2_price149_by_stratum.tsv` | same | Same arms split by metabolic stratum (A core biomass, B catabolic assay-specific, C secondary/peripheral) |
| `results/pillar2_price149_per_protein_real.tsv` | same | One row per (baseline, arm, target protein) |

### `pillar2_price149_all6_real.tsv` columns

`extraction_version` (v2-2026-10-19; v1-legacy for the old pickle), `baseline`, `method`
(arm), `pooled_predictors` (predictors whose scores enter the arm; for the Multi-Pred arms
the baseline first, then the noisy-OR pool; `pool_without_this_baseline` lists only the
other predictors), `n_proteins` (148), `top1_count`, `top1_acc`, `top1_ci95_lo/hi`
(bootstrap, 2,000 resamples, seed 42), `top5_count`, `top5_acc`, `fmax` (mean per-protein
F-max over each protein's own thresholds), `fmax_dataset_macro` (protein-centric CAFA
F-max at one global τ ∈ {0.01..0.99}), `fmax_dataset_micro` (micro-averaged over all
(protein, EC) pairs at one global τ), `best_tau` (τ of the macro F-max), `f1_at_05`, `mrr`,
`n_improved_vs_base` (McNemar b: baseline wrong → arm right at Top-1),
`n_degraded_vs_base` (c), `mcnemar_p` (exact two-sided binomial test on b vs c, 4
significant digits).

Arms per baseline: Baseline; + METEOR v8 (Iter 1); + METEOR v8 (Iter 2); + GSPA SLR
(Minimal Glucose); + GSPA SLR + RLGC (Stoichiometric Local Context); + GSPA Single-Pred
(Minimal Glucose); + Full GSPA (Iter 1/2/3); + Stage 1 Noisy-OR Ensemble (Sequence Only);
+ Full GSPA Multi-Pred (Minimal Glucose); + Full GSPA Multi-Pred (Iter 3);
pool_without_this_baseline.

`pillar2_price149_by_stratum.tsv` and `pillar2_price149_per_protein_real.tsv` also carry
`extraction_version` as their first column.
