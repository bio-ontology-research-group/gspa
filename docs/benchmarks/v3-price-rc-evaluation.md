# Price enzyme benchmark under the release-candidate JAR

Status: complete for the six single-predictor panels (2026-10-30). Collection: Shaheen
`/scratch/project/k10018/gspa-v3/out/price-production-rc1` (132 cells = 6 predictors × 22 genomes, two arms each,
1,526,919,081 claims; registration 90a960f8…, release-candidate JAR 3fe6ba0a…; every cell's publication verified by
the unchanged finalizer, job 15708110). Scoring: the unchanged evaluator `benchmark/v3/evaluate_price_reference.py`
(15/15 tests) in Shaheen job 15713513 via `benchmark/v3/shaheen/price-rc-evaluate.sbatch`; the independent arithmetic
check passes for every predictor. Results: `benchmark/v3/results/price-reference-rc1/` (`summary_exact.tsv` and each
predictor's `report.json`).

## What is compared

Truth: the 148 Price et al. (2018) proteins with 151 experimentally assigned enzyme numbers in 22 bacteria (one
replicon row excluded as documented). For each predictor: its raw exported scores, GSPA integration without priors
(stage 0) and GSPA integration with the default priors (stage 1). Top-1 counts a protein as correct when its highest-
scoring enzyme matches the reference exactly. "Calls" are integrated scores of at least 0.5, a fixed diagnostic
threshold, not a calibrated operating point. A call on a reference protein that names a different enzyme is not
counted as correct but is also not proven wrong (the protein may have a second activity).

## Results (exact enzyme numbers; 148 proteins, every one has candidates)

| Predictor | Top-1 correct | Top-1 share, 95 % interval (genome bootstrap) | Top-3 | Top-5 | Mean reciprocal rank | Calls matching the reference | Other calls on reference proteins |
|---|---:|---|---:|---:|---:|---:|---:|
| CLEAN | 60 | 0.405 [0.330, 0.485] | 87 | 95 | 0.516 | 40 | 50 |
| DeepProZyme | 46 | 0.311 [0.248, 0.371] | 90 | 101 | 0.479 | 14 | 26 |
| EnzBERT | 53 | 0.358 [0.280, 0.429] | 76 | 84 | 0.443 | 41 | 34 |
| GraphEC | 18 | 0.122 [0.047, 0.193] | 37 | 52 | 0.217 | 8 | 20 |
| MAPred | 50 | 0.338 [0.265, 0.410] | 95 | 99 | 0.494 | 61 | 136 |
| TopEC | 5 | 0.034 [0.006, 0.070] | 5 | 5 | 0.034 | 2 | 83 |

The table is the same for the raw export, stage 0 and stage 1 of every predictor (all three columns sets are in
`summary_exact.tsv`): stage 1 is byte-identical to stage 0 in all 132 cells, because the default priors act on GO
claims only, these cells hold enzyme claims only, and no gap or operon inputs were given.

## What this shows

GSPA's integration without the metabolic refiner neither helps nor harms any of the six enzyme predictors on the Price
proteins; the comparison measures the predictors. CLEAN ranks the true enzyme first for 60 of 148 proteins, as in the
earlier development run. The test of GSPA's metabolic claim is the separate Stage 2B collection
(`price-stage2b-rc1`, flux-certified refinement switched on), reported in `docs/benchmarks/v3-price-stage2b.md` when
complete.
