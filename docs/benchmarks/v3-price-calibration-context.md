# Calibration-only thresholds and reaction-context modes on the Price benchmark (rows 1 and 5)

Status: rows 1 and 5 complete (2026-11-02).
Pre-registered plan: `tasks/v3-row1-row5-price-calibration-context-2026-11-01.md` (SHA-256
`02caeeb9…`, recorded in `benchmark/v3/results/price-calibration-context-v1/PLAN.sha256` before any check ran;
deviations are appended there, dated). Evidence: `benchmark/v3/results/price-calibration-context-v1/` (`row1/`,
`row5/`, each with `SHA256SUMS`). Code: `benchmark/v3/price_calibration_context/` (`row1_calibration.py`,
`build_context_inputs.py`, `run_context_cell.py`, `evaluate_context.py`, `archive_results.py`), tests in
`benchmark/v3/tests/test_price_calibration_context.py` (13 tests), Shaheen wrapper
`benchmark/v3/shaheen/price-context-cells.sbatch`. Every number below is read from files in the evidence folder.

The labels are the 148 reannotated Price et al. (2018) proteins (151 enzyme numbers) in the 22 benchmark genomes,
scored with the unchanged Price evaluator conventions (`price_reference_metrics.py`, exact-ID policy). A call on a
reference protein that names another enzyme number is "unconfirmed", never a negative. No AllTheBacteria data is
used; GSPA defaults are unchanged.

## Row 1: calibration-only thresholds and calibration metrics

### What was done

No inference was rerun. The per-protein score files that the Price evaluator retains for the 148 reference proteins
(every candidate enzyme number with its score) were read for four arms of the existing release-candidate
collections: the predictor's own exported scores (`raw_export`), integration without priors (`stage0`; `stage1`
with default priors is byte-identical to it in all 132 cells and is not shown separately), the taxon/EC consistency
fix (`stage1_taxon`, release candidate 2) and the flux-certified refiner (`stage2b`). Six predictors: CLEAN,
DeepProZyme, EnzBERT, GraphEC, MAPred, TopEC. Check C6: the top-1 counts recomputed from these files equal the
published summaries for every predictor and arm (`row1/report.json: integrity_all_equal = true`).

**Split (C1).** Genus-blocked, with Stutzerimonas merged into the Pseudomonas block (12 blocks). The pre-registered
hash-order greedy rule failed its own minimum (36 calibration proteins; `row1/split_a2_rule_failed.json`), so the
dated replacement rule, decided before any threshold or metric, walks the blocks from the largest down with the
same greedy balancing. Result (`row1/split.json`): calibration = the seven Pseudomonas/Stutzerimonas genomes
(78 reference proteins), test = the other 15 genomes from 11 genera (70 proteins, 73 reference pairs). Thresholds
are thus chosen on one genus and tested on eleven others. The swapped assignment is reported as a sensitivity
analysis.

**Threshold rule (C2).** The existing fixed-precision rule (`evaluate_wp5_recovery.py`: at least 10 calls with at
least 90 % demonstrated support, then most supported calls, fewest unsupported, highest threshold), applied to the
candidate (protein, enzyme number) pairs of the calibration reference proteins only; thresholds were written and
hashed (`row1/policy.json`, `eeba6b02…`) before any test metric.

### Result: no predictor and no arm obtains a calibrated threshold

Under the rule, **every one of the 24 predictor × arm combinations abstains** (`no_supported_threshold_abstain`),
in the primary and in the swapped split. The diagnostic recorded with each abstention shows how far the calls are
from the bar: at the first threshold reaching 10 calls, the best support is 6/10 (CLEAN, DeepProZyme), and the
highest support any threshold with at least 10 calls achieves on the calibration genomes is 12/18 (0.667, CLEAN) in
the primary split and 8/11 (0.727, EnzBERT) in the swapped split; TopEC reaches 0/10. The pairs that compete for
the top of the ranking are mostly alternative enzyme numbers of the same reference proteins (CLEAN: 408,876 candidate
pairs on 78 calibration proteins, 78 of them reference pairs).

Consequently the calibrated-call metrics of the plan (calls, supported share, reference recall at the calibrated
threshold; corrections and harms at that threshold) are empty by construction for every arm, and only the declared
fixed 0.5 diagnostic and the threshold-free top-k and calibration metrics carry information. This is the row-1
answer: on the Price benchmark the fixed-precision policy that GSPA uses to turn scores into calls cannot be
demonstrated for any of the six enzyme predictors, in any processing stage, from a calibration set of this size.

### Test-genome metrics (primary split; 15 genomes, 70 proteins, 73 reference pairs)

Top-1 is the tie-aware share of reference proteins whose top-ranked enzyme number is a reference number
(genome-bootstrap 95 % interval, 2,000 resamples, seed 13). Fixed 0.5 counts (protein, enzyme) calls with score
≥ 0.5: M = reference matches, U = unconfirmed. Brier (positive pairs) is the mean squared distance of the reference
enzyme number's score from 1 (a missing candidate scores 0; label-safe). Top-call Brier and 10-bin ECE score the
top-ranked call of each reference protein against "equals a reference number" (reference support, not error rate).

| Predictor | Arm | Top-1 share [95 % CI] | Top-3 | 0.5: M/U | Brier (positive pairs) | Top-call Brier | ECE |
|---|---|---|---:|---|---:|---:|---:|
| CLEAN | raw_export | 0.343 [0.253, 0.475] | 37 | 15/28 | 0.765 | 0.377 | 0.392 |
| CLEAN | stage0 | 0.343 [0.253, 0.475] | 37 | 13/18 | 0.708 | 0.255 | 0.198 |
| CLEAN | stage1_taxon | 0.343 [0.253, 0.475] | 37 | 13/18 | 0.708 | 0.255 | 0.198 |
| CLEAN | stage2b | 0.343 [0.253, 0.475] | 33 | 13/18 | 0.708 | 0.255 | 0.198 |
| DeepProZyme | raw_export | 0.257 [0.204, 0.294] | 41 | 18/46 | 0.586 | 0.374 | 0.406 |
| DeepProZyme | stage0 / taxon / 2B | 0.257 [0.204, 0.294] | 41 | 5/15 | 0.658 | 0.222 | 0.174 |
| EnzBERT | raw_export | 0.271 [0.185, 0.386] | 31 | 17/27 | 0.729 | 0.323 | 0.355 |
| EnzBERT | stage0 / taxon | 0.271 [0.185, 0.386] | 31 | 16/20 | 0.738 | 0.173 | 0.119 |
| EnzBERT | stage2b | 0.271 [0.185, 0.386] | 30 | 16/20 | 0.738 | 0.173 | 0.119 |
| GraphEC | raw_export | 0.029 [0.000, 0.067] | 6 | 2/15 | 0.895 | 0.133 | 0.302 |
| GraphEC | stage0 / taxon / 2B | 0.029 [0.000, 0.067] | 6 | 1/4 | 0.783 | 0.082 | 0.231 |
| MAPred | raw_export | 0.329 [0.227, 0.477] | 43 | 39/180 | 0.410 | 0.518 | 0.550 |
| MAPred | stage0 / taxon / 2B | 0.329 [0.227, 0.477] | 43 | 30/59 | 0.528 | 0.258 | 0.210 |
| TopEC | raw_export | 0.043 [0.000, 0.128] | 3 | 3/47 | 0.963 | 0.519 | 0.643 |
| TopEC | stage0 / taxon / 2B | 0.043 [0.000, 0.128] | 3 | 2/39 | 0.969 | 0.263 | 0.445 |

Rows merged with "/" are identical in every column. Full tables: `row1/thresholds_primary.tsv`,
`row1/calibration_primary.tsv`, `row1/reliability_primary.tsv` (per-bin counts for the top calls and for all
candidate pairs); the swapped split in the `_swapped` files.

What the table shows, and what it does not. Integration (Stage 0) leaves the ranking of every predictor unchanged on
the test genomes (identical top-1 and, except for Stage 2B, top-3), lowers the number of unconfirmed calls at 0.5
for every predictor while keeping between 1 and 9 fewer reference matches, and lowers the top-call Brier and ECE for
every predictor. The reliability tables explain the ECE change: after integration the top-call scores of all six
predictors lie in [0, 0.6), with the bulk of the confident calls in the single bin [0.5, 0.6) (CLEAN 24 of 70 top
calls, 11 supported; EnzBERT 36, 16 supported; MAPred 55, 20 supported), so the integrated score is a compressed
rescaling whose nominal value is closer to the observed support share than the predictors' own near-1 scores. The
positive-pair Brier moves in both directions (lower for CLEAN and GraphEC, higher for the other four), because
compressing the score of the reference enzyme number away from 1 costs under this metric. None of this is a
calibration guarantee: the test set has 70 proteins, the intervals are wide, and "supported" means reference
agreement, not correctness.

The taxon/EC fix changes no score on the 148 proteins (as in the taxon report), and Stage 2B changes the ranking of
four CLEAN proteins (top-3 37 → 33) and one EnzBERT protein without any top-1 change; the contrasts
(`row1/contrasts_primary.tsv`) record 0 corrections and 0 harms at 0.5 for both, and for `raw_export → stage0`
between 1 and 13 reference pairs that are called at 0.5 only before integration (DeepProZyme 13, MAPred 9), none
only after.

Swapped split (thresholds on the 15 non-Pseudomonas genomes, test on the 7 Pseudomonas genomes, 78 proteins): the
same abstention for all 24 combinations; Stage 0 top-1 shares 0.462 (CLEAN), 0.436 (EnzBERT), 0.359
(DeepProZyme), 0.346 (MAPred), 0.205 (GraphEC), 0.026 (TopEC).

### Row-1 verdict against the plan

C1 split: pass under the dated deviation (the pre-registered rule failed and is recorded). C2 thresholds: pass as a
procedure, every combination abstains. C3 test metrics: produced; the calibrated-call part is empty by construction.
C4 calibration metrics: pass. C5 contrasts: pass. C6 integrity: pass. The acceptance element "calibration-only
thresholds on the Price benchmark" is therefore shown as *measured and not obtainable* for these predictors under
GSPA's own call rule, not as a threshold that could be reported.

## Row 5: reaction-context modes with actual gaps

Results: `benchmark/v3/results/price-calibration-context-v1/row5/` (`SHA256SUMS`; plan with deviations
`a1dddbfa…`, `PLAN-with-deviations.sha256`). Cells: 132 (6 predictors × 22 genomes), 462 arms, run on Shaheen.
Raw cell outputs stay on Shaheen (`/scratch/project/k10018/gspa-v3/out/price-context-v1`) and on ws
(`~/tmp/price-row5-out`, `~/tmp/price-context-cells`).

What ran (see the B4 deviations): `stage0` and `gaps` for all six predictors; the three context modes
(`universal`, `pathway`, `organism`) for TopEC, MAPred and EnzBERT only. For CLEAN, DeepProZyme and GraphEC a
five-arm cell would take 9 to 42 hours, beyond the pre-registered 12-hour rule, so the mode contrast on those three
predictors was not measured.

### Gaps and power

| Quantity | Value |
|---|---|
| Genomes with gapsmith observations converted | 22 / 22 |
| Gap records (all genomes) | 1,714 |
| Reference proteins (Price loci) | 148 |
| Gap-linked Price loci (a reference EC is an EC of a gap record of its genome) | 8 |
| Gap-linked loci with context available, universal / pathway / organism | 8 / 8 / 8 |
| Power gate (pre-registered) | 20 gap-linked loci |

With 8 gap-linked loci the mode ablation is **unpowered** under the pre-registered gate: most Price enzymes are not
gap fillers of their genome's gapsmith model.

### Checks

| Check | Result |
|---|---|
| D1 gaps | pass: 22 genomes converted with table hashes; gap context `available` in all 330 non-`stage0` arms |
| D2 cells | pass: 462 / 462 arms exit 0; prediction rows equal input claims and the candidate universe is preserved in every arm; `stage0` equals rc1's retained scores for all six predictors; gap records identical across predictors for every genome |
| D3 prior contrast | pass: produced for all six predictors (table below) |
| D4 mode contrast | pass for TopEC, MAPred, EnzBERT (table below); not measured for CLEAN, DeepProZyme, GraphEC (scale rule) |
| D5 power | pass: produced; gate not met (8 < 20) |

### Price metrics per arm (148 reference proteins, 22 genomes)

Every arm of a predictor gives the same numbers, so one row per predictor is shown; arms listed are those that ran.

| Predictor | Arms run | Top-1 match | Top-1 share [95 % CI] | Top-3 match | Matched at 0.5 | Unconfirmed at 0.5 |
|---|---|---|---|---|---|---|
| CLEAN | stage0, gaps | 60 | 0.405 [0.330, 0.485] | 87 | 40 | 50 |
| DeepProZyme | stage0, gaps | 46 | 0.311 [0.248, 0.371] | 90 | 14 | 26 |
| EnzBERT | all five | 53 | 0.358 [0.280, 0.429] | 76 | 41 | 34 |
| GraphEC | stage0, gaps | 18 | 0.122 [0.047, 0.193] | 37 | 8 | 20 |
| MAPred | all five | 50 | 0.338 [0.265, 0.410] | 95 | 61 | 136 |
| TopEC | all five | 5 | 0.034 [0.006, 0.070] | 5 | 2 | 83 |

No predictor has a row-1 calibrated threshold (all abstain, see row 1), so the calibrated-threshold contrasts are
empty by construction.

### Gap prior: `stage0` → `gaps`

| Predictor | Output rows | Scores changed (all proteins) | Price loci with a changed score | Top-1 gained / lost | Corrected / harmed at 0.5 |
|---|---|---|---|---|---|
| CLEAN | 605,634,470 | 75,852 | 1 | 0 / 0 | 0 / 0 |
| DeepProZyme | 323,731,872 | 50,568 | 1 | 0 / 0 | 0 / 0 |
| EnzBERT | 11,553,500 | 2,019 | 0 | 0 / 0 | 0 / 0 |
| GraphEC | 583,869,869 | 75,762 | 1 | 0 / 0 | 0 / 0 |
| MAPred | 2,014,189 | 234 | 0 | 0 / 0 | 0 / 0 |
| TopEC | 115,181 | 29 | 0 | 0 / 0 | 0 / 0 |

The gap prior acts on actual gaps, but on the Price loci it changes one score per dense predictor, and no call or
top-1 answer.

### Context modes (TopEC, MAPred, EnzBERT)

| Predictor | Output rows per arm | Differing scores, each of the six mode pairs (gaps vs each mode, and mode vs mode) | Learned stage: candidates scored / abstained, each mode |
|---|---|---|---|
| EnzBERT | 11,553,500 | 0 | 0 / 11,553,500 |
| MAPred | 2,014,189 | 0 | 0 / 2,014,189 |
| TopEC | 115,181 | 0 | 0 / 115,181 |

Abstention reason in all 198 mode arms: "no candidates supported by the model input-coverage manifest; baseline
scores retained". The modes feed only the learned stage, and the learned models carry no input-coverage contract for
these predictors, so the three modes cannot change any score. This is the mechanism the plan predicted (B1).

### Compute

| Predictor | Arm | Total time (22 cells) | Slowest cell |
|---|---|---|---|
| CLEAN | stage0 / gaps | 3.2 h / 17.9 h | 13 / 79 min |
| DeepProZyme | stage0 / gaps | 2.0 h / 10.3 h | 8 / 45 min |
| GraphEC | stage0 / gaps | 3.2 h / 17.5 h | 13 / 80 min |
| EnzBERT | stage0 / gaps / each mode | 0.5 h / 0.6 h / 10.4–11.0 h | 2 / 3 / 47–51 min |
| MAPred | stage0 / gaps / each mode | 0.2 h / 0.5 h / 1.9–2.3 h | 1 / 2 / 9–11 min |
| TopEC | stage0 / gaps / each mode | 0.1 h / 0.3 h / 0.4–0.5 h | 1 / 2 / 2 min |

### Row-5 verdict against the plan

All five checks pass as procedures. The result is negative on two counts. First, it is unpowered: 8 gap-linked Price
loci against a gate of 20. Second, it is structurally null: the context modes change 0 scores because the learned
stage abstains on every candidate, and the gap prior changes no Price call. The acceptance bar (pathway-only versus
universal versus organism-specific ablation on held-out loci with actual gaps) is therefore **not met**. A fair test
needs two things: a predictor whose claims the learned models accept, for example native CLEAN with the pinned cohort
adapter and the predictor contract, and a locus set with at least 20 gap-linked enzymes. The dense-predictor mode
arms also need roughly 9 to 42 hours per cell, or a faster cohort pass.
