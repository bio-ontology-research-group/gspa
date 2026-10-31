# Stage 2B flux-certified refinement on the Price enzymes

Status: complete (2026-10-31). Collection: Shaheen `/scratch/project/k10018/gspa-v3/out/price-stage2b-rc1`
(132 cells = 6 predictors × 22 genomes; arms stage 0 = integration without priors, stage 2B = default priors plus
`--enable-fba-refiner`; verified by the collection finalizer, job 15712378). Stage 2B data directory built
reproducibly from the pinned gapsmith reference (`benchmark/v3/stage2b/build_fba_data.py`; data_merged sha256 25bad8d2…;
Gram-negative biomass, minimal glucose medium; all 22 genomes are Gram-negative Pseudomonadota). Scoring:
`benchmark/v3/stage2b/evaluate_price_stage2b.py` with the same 148-protein reference, enzyme-identity audit (5302494b…)
and metrics as the release-candidate Price report; jobs 15724205 and 15729831. Results:
`benchmark/v3/results/price-stage2b-rc1/` (`summary_exact.tsv`, one `report.json` per predictor).

## Why this test

The README's v2.1.0 table credits "GSPA Full" (Stage 2B refinement) with a gain of about +0.10 F-max over each enzyme
predictor on Price-149. That run's data directory came from a colleague's partial snapshot and could not be reproduced.
This is the first reproducible test of the claim.

## Results (exact enzyme numbers; 148 proteins)

| Predictor | Top-1, stage 0 → 2B | Top-3, stage 0 → 2B | Mean reciprocal rank, stage 0 → 2B | Top call changed | Corrected | Harmed | CAFA F-max change, 95 % interval |
|---|---|---|---|---:|---:|---:|---|
| CLEAN | 60 → 59 | 87 → 80 | 0.516 → 0.489 | 10 | 0 | 1 | 0.000 [0.000, 0.000] |
| DeepProZyme | 46 → 46 | 90 → 90 | 0.479 → 0.479 | 0 | 0 | 0 | 0.000 [0.000, 0.000] |
| EnzBERT | 53 → 53 | 76 → 75 | 0.443 → 0.441 | 0 | 0 | 0 | 0.000 [0.000, 0.000] |
| GraphEC | 18 → 18 | 37 → 37 | 0.217 → 0.217 | 0 | 0 | 0 | 0.000 [0.000, 0.000] |
| MAPred | 50 → 50 | 95 → 95 | 0.494 → 0.494 | 0 | 0 | 0 | 0.000 [0.000, 0.000] |
| TopEC | 5 → 5 | 5 → 5 | 0.034 → 0.034 | 0 | 0 | 0 | 0.000 [0.000, 0.000] |

Micro F-max is also unchanged for every predictor (`summary_exact.tsv`).

## What this shows

Stage 2B does not improve enzyme prediction on the Price proteins for any predictor; for CLEAN it changes ten top
calls, none for the better and one for the worse, and pushes correct enzymes out of the top three. The refiner runs
in every cell (models, repairs and flux certificates are written), but it only adds the few reactions a model needs to
grow on the default medium and raises the scores of genes for those reactions; in a typical cell one score of about
fifty thousand changes. Almost none of the Price enzymes lies on that growth path (2 of 148 in the off-growth
experiment), so biomass-driven refinement cannot reach them. The README's +0.10 F-max is not reproduced.

## Limits

* One medium (minimal glucose) and one biomass template; the Price assays used other substrates.
* The universal model and bounds are rebuilt from gapsmith rather than the original snapshot, so this tests Stage 2B
  as it can be reproduced, not the unreproducible original run.
