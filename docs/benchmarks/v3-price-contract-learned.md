# The learned stage on the Price benchmark with claims from its own training predictor (rows 1, 4, 5, 11)

Date: 2026-11-02. Frozen plan: `tasks/v3-price-contract-learned-2026-11-02.md` (SHA-256 of the text up to the
deviations heading `0b1096bf…`, recorded in `benchmark/v3/results/price-contract-learned-v1/PLAN.sha256` before any
check ran; dated deviations appended). Evidence: `benchmark/v3/results/price-contract-learned-v1/` (`SHA256SUMS`,
223 files; claim files and integrated outputs over 5 MB are listed in `LARGE_FILES.tsv` and kept on ws under
`/mnt/data2/leechuck/gspa-v3-data/price-contract-learned-v1/` and on Shaheen under
`/scratch/project/k10018/gspa-v3/out/price-contract-learned-v1/`). Code: `benchmark/v3/price_contract_learned/`
(`leakage_check.py`, `build_stage.py`, `run_cell.py`, `run_anomalies.py`, `evaluate.py`, `row11_truth.py`), tests
`benchmark/v3/tests/test_price_contract_learned.py` (4 tests), Shaheen wrapper
`benchmark/v3/shaheen/price-contract-learned.sbatch`. Every number below is read from files in the evidence folder.

In the earlier Price analyses the learned stage abstained because none of the six Price enzyme predictors carries a
predictor contract. Here the Price genomes get claims from the one predictor the contract-3 models were trained on,
the genus-fenced Swiss-Prot DIAMOND transfer, with its contract. Labels and scoring are unchanged: the 148
reannotated Price et al. (2018) proteins (151 enzyme numbers), exact-ID policy of `price_reference_metrics.py`,
genome bootstrap (2,000 resamples, seed 13), the frozen row-1 genus-blocked split
(`price-calibration-context-v1/row1/split.json`) and the WP5 call rule (≥ 10 calls at ≥ 90 % reference support).
A call naming another enzyme on a reference protein is "unconfirmed", never a negative. Engine: the pinned
release-candidate JAR `08f8b5a1…` (not rebuilt); models: contract-3 logistic `c9df93d6…` and GBDT `d319630b…`.

## Verdict per check

| Check | Result |
|---|---|
| L leakage | Pass (measured). No Price-148 protein has an exact or ≥ 95 % identical reference entry; 25 reference entries belong to a Price organism itself (below) |
| C claims | Pass. 22 / 22 genomes, every contract equals the model's (`ae823453…`); the three genomes run before are byte-identical (claims and hits) |
| I cells | Pass. 22 genomes × 14 arms, all exit 0 (Shaheen 15753627, 15753628); the learned stage executed in 6 arm types on 22 / 22 genomes, with 0 candidates lacking a contract in every learned arm |
| E rows 1 and 4 | Pass as a procedure. Policy written and hashed before test metrics (`evaluation/policy.json`, `f75fb98b…`). The call rule finds no threshold for any arm in either split |
| F row 5 | Pass as a procedure; the ablation is unpowered (8 gap-linked loci, gate 20) and the three modes give identical scores |
| G row 11 | G1, G2 pass: 484 real cases, 0 candidate pairs with known truth. G3 not applicable under the pre-registered condition (no known truth in either split); the agent endpoint was reachable but not used |

## Leakage (step L)

| Quantity | Value |
|---|---|
| Price-148 proteins searched (layout proteome) | 148 |
| with an exact reference sequence | 0 |
| with a reference entry ≥ 95 % identity, ≥ 80 % query and subject coverage | 0 |
| best identity of any Price-148 protein to the reference | 94.1 % (one protein, GCA_000006965.1 CAC47870.1, Q06191) |
| Price-148 proteins with best identity ≥ 80 % / ≥ 50 % / no hit | 3 / 32 / 13 |
| reference entries whose Swiss-Prot taxid is a Price organism | 25, all *Sinorhizobium meliloti* 1021 (genome GCA_000006965.1) |
| Swiss-Prot experimental Price labels whose accession is a reference entry | 15 (all *S. meliloti*) |

The genus fence of the training reference compares the first word of the Swiss-Prot organism name; *S. meliloti* is
named "Rhizobium meliloti (…) (Ensifer meliloti) (Sinorhizobium meliloti)" there, and only "ensifer" is fenced, so
its entries passed. None of them is a Price-148 protein, so the pre-registered exclusion set is empty; the dated
deviation adds a sensitivity analysis without that genome (below), and row-11 truth from those 15 labels counts as
circular.

## Claims (step C, ws, DIAMOND 2.1.11 `d6bea1f3…`)

| Quantity | Value |
|---|---|
| Genomes | 22 |
| Proteins / with at least one hit | 115,542 / 71,537 |
| Claims (GO / EC) | 1,234,318 (1,150,060 / 84,258) |
| Price-148 proteins with at least one EC candidate (any hit: 135) | 126 |

Per-genome counts and hashes: `claims/claims_report.json`.

## Rows 1 and 4: the learned stage on Price (step E)

### Learned stage coverage (all 22 genomes)

| Arm | Learned stage executed (genomes) | Candidates scored / unsupported | Price-148 EC candidates scored / unsupported | Main unsupported reasons |
|---|---|---|---|---|
| `learned` (logistic) | 22 | 919,227 / 315,091 | 460 / 2 | GO profile lacks both classes in training and calibration 296,245; unseen availability profile 16,518; score outside the calibration range 2,328 |
| `learned_gbdt` | 22 | 919,227 / 315,091 | 460 / 2 | same |
| `rlc_learned` (suggester, no gap observations) | 22 | 919,227 / 315,091 | 460 / 2 | same; scores identical to `learned` |
| `obs_learned`, `gaps_learned` (gap observations) | 0 | 0 / 1,234,318 | 0 / 462 | unseen availability profile, every candidate |
| `universal`, `pathway`, `organism` (gaps + cohort) | 0 | 0 / 1,234,318 | 0 / 462 | unseen availability profile / unseen reaction-context mode |
| `cohort_universal`, `cohort_pathway`, `cohort_organism` (cohort, no gaps) | 22 | 863,298 / 371,020 | 149 / 313 | unseen availability profile (universal) or unseen reaction-context mode (pathway, organism) |

Candidates without a predictor contract: 0 in every learned arm. The gap observations alone make every candidate's
feature-availability profile one the model never saw in training, so the model declines all of them; the
reaction-local-context suggester alone does not.

### Calibration-only thresholds (the call rule on the calibration genomes)

| Split | Arm | Calibration candidate pairs (reference pairs) | Highest support at ≥ 10 calls | Rule result |
|---|---|---|---|---|
| primary (7 Pseudomonas/Stutzerimonas genomes) | raw / stage0 / stage1 | 243 (38) | 14 / 30 (0.467) | abstain |
| primary | learned | 243 (38) | 19 / 55 (0.345) | abstain |
| primary | learned_gbdt | 243 (38) | 10 / 25 (0.400) | abstain |
| swapped (15 other genomes) | raw / stage0 / stage1 | 219 (29) | 5 / 12 (0.417) | abstain |
| swapped | learned | 219 (29) | 10 / 21 (0.476) | abstain |
| swapped | learned_gbdt | 219 (29) | 4 / 10 (0.400) | abstain |

No arm reaches 90 % support at any threshold, so calibrated calls, and corrections and harms at a calibrated
threshold, are empty by construction for every arm (as for the six enzyme predictors).

### Test-genome metrics

S/U = supported / unconfirmed calls. "Model threshold" is each model's frozen display threshold, chosen on its own
training calibration split (logistic 0.321, GBDT 0.223); S/(S+U) there is the supported precision (a lower bound on
precision, since unconfirmed calls are not shown to be wrong). Brier on positive pairs is label-safe; top-call Brier
and 10-bin ECE score the top call against reference support.

| Split (test proteins) | Arm | Top-1 share [95 % CI] | Top-3 | At 0.5: S/U | At model threshold: S/U | Brier (positive pairs) | Top-call Brier | ECE |
|---|---|---|---:|---|---|---:|---:|---:|
| primary (70) | raw | 0.271 [0.198, 0.388] | 28 | 7 / 11 | | 0.761 | 0.170 | 0.095 |
| primary (70) | stage0 / stage1 | 0.271 [0.198, 0.388] | 28 | 1 / 3 | | 0.813 | 0.161 | 0.061 |
| primary (70) | learned | 0.271 [0.198, 0.388] | 28 | 5 / 8 | 11 / 19 | 0.806 | 0.152 | 0.067 |
| primary (70) | learned_gbdt | 0.251 [0.176, 0.367] | 26 | 1 / 1 | 10 / 21 | 0.879 | 0.147 | 0.116 |
| swapped (78) | raw | 0.333 [0.230, 0.430] | 37 | 4 / 14 | | 0.684 | 0.196 | 0.194 |
| swapped (78) | stage0 / stage1 | 0.333 [0.230, 0.430] | 37 | 0 / 3 | | 0.758 | 0.193 | 0.177 |
| swapped (78) | learned | 0.321 [0.209, 0.424] | 37 | 2 / 11 | 13 / 25 | 0.765 | 0.204 | 0.181 |
| swapped (78) | learned_gbdt | 0.325 [0.231, 0.414] | 32 | 0 / 0 | 12 / 27 | 0.844 | 0.200 | 0.185 |

Restricted to the candidates the learned model scored, the learned rows are unchanged except one call at 0.5
(`evaluation/test_metrics.tsv`, columns `ls_*`). Reliability tables: `evaluation/reliability.tsv`.

### Corrections and harms (test genomes, fixed 0.5; no calibrated threshold exists)

Corrections are reference pairs called only after, harms only before.

| Split | Contrast | Corrected | Harmed | Net [95 % CI] | Top-1 gained / lost |
|---|---|---:|---:|---|---|
| primary | raw → stage1 | 0 | 6 | −6 [−11, −2] | 0 / 0 |
| primary | stage0 → stage1 | 0 | 0 | 0 | 0 / 0 |
| primary | stage1 → learned | 4 | 0 | 4 [1, 8] | 0 / 0 |
| primary | stage1 → learned_gbdt | 1 | 1 | 0 [−3, 3] | 0 / 1 |
| swapped | raw → stage1 | 0 | 4 | −4 [−8, −1] | 0 / 0 |
| swapped | stage1 → learned | 2 | 0 | 2 [0, 4] | 0 / 2 |
| swapped | stage1 → learned_gbdt | 0 | 0 | 0 | 0 / 1 |

### All 22 genomes (148 proteins)

| Arm | Top-1 | Top-1 share [95 % CI] | Top-3 | Matched / unconfirmed at 0.5 |
|---|---:|---|---:|---|
| raw | 49 | 0.304 [0.236, 0.378] | 65 | 11 / 25 |
| stage0, stage1, gaps and the five gap-based learned arms | 49 | 0.304 [0.236, 0.378] | 65 | 1 / 6 |
| learned, rlc_learned | 47 | 0.297 [0.226, 0.373] | 65 | 7 / 19 |
| learned_gbdt | 47 | 0.290 [0.224, 0.361] | 58 | 1 / 1 |
| cohort_universal, cohort_pathway, cohort_organism | 49 | 0.311 [0.249, 0.379] | 64 | 1 / 6 |

For comparison, the six enzyme predictors' stage-0 top-1 on the same 148 proteins is 5 to 60 (CLEAN 60,
`v3-price-calibration-context.md`).

### Sensitivity without GCA_000006965.1 (primary split, 52 test proteins; the swapped test set does not contain it)

| Arm | Top-1 share [95 % CI] | At 0.5: S/U | At model threshold: S/U | Brier (positive pairs) | ECE | stage1 → arm: corrected / harmed |
|---|---|---|---|---:|---:|---|
| raw | 0.292 [0.189, 0.426] | 5 / 7 | | 0.742 | 0.132 | |
| stage1 | 0.292 [0.189, 0.426] | 0 / 0 | | 0.799 | 0.079 | |
| learned | 0.292 [0.189, 0.426] | 4 / 4 | 8 / 12 | 0.797 | 0.085 | 4 / 0 |
| learned_gbdt | 0.264 [0.159, 0.402] | 1 / 0 | 7 / 15 | 0.871 | 0.124 | 1 / 0 |

The direction of every contrast is the same without the genome whose own proteins are in the reference.

What this shows. The learned stage now runs on Price and scores almost every enzyme candidate of the reference
proteins (460 of 462). Its probabilities are not well enough separated for GSPA's call rule on either split; at the
model's own frozen threshold about one call in three is supported (0.34–0.37). Compared with integration without
the model, the logistic model adds 6 reference pairs at 0.5 with no harm and loses 2 top-1 answers, all in the
swapped split; the tree model lowers top-3 and gains nothing. Integration itself (raw → stage 1) compresses the
DIAMOND scores and drops 10 reference pairs below 0.5. ECE moves in both directions and the intervals are wide
(70 and 78 test proteins).

## Row 5: context modes with the learned stage (step F)

### Score differences between modes (whole outputs, 1,234,318 rows; and the 148 loci)

| Comparison | Differing scores (all rows) | Price loci with a changed score |
|---|---:|---:|
| stage1 → gaps (gap prior) | 1,133 | 0 |
| gaps_learned vs universal / pathway / organism, and the three mode pairs | 0 | 0 |
| cohort_universal vs cohort_pathway / cohort_organism, and pathway vs organism | 0 | 0 |
| learned → cohort_universal (cohort context added) | 55,929 | 119 |
| learned → gaps_learned (gap observations added) | 916,737 | 126 |

With gaps the learned stage declines every candidate, so the mode arms equal the gap arm. Without gaps the learned
stage scores, but the three modes still give identical scores: the reference EC candidates whose features depend on
the mode are declined (pathway and organism: unseen reaction-context mode; universal: unseen availability profile)
and the others are scored the same in all three. Price metrics: the cohort arms give top-1 49, the same as stage 1,
with 3 top-1 gained and 1 lost against `learned` and 6 reference pairs no longer called at 0.5.

### Power gate (unchanged gap-linked loci; gap records equal the row-5 run for all 22 genomes)

| Quantity | Value |
|---|---|
| Gap-linked Price loci | 8 |
| of which these claims propose the reference EC | 3 (two glucokinases 2.7.1.2, top-1 in every arm; one 1.3.8.7, not top-1) |
| Pre-registered gate | 20 |
| Powered | no |

Row 5 stays unpowered and structurally null on this predictor as on the six enzyme predictors.

## Row 11: anomaly extraction on GO-bearing outputs (step G)

| Quantity | Value |
|---|---|
| Arms run through the production extractor (Shaheen 15753686) | stage1, learned_logistic, learned_gbdt; 22 / 22 genomes each, receipts verified by `prepare` |
| Cases (calibration / test) | 484 (239 / 245); 0 genomes at the 150-case cap |
| Case types | pathway gap 464, borderline essential 9, borderline conflict 6, incomplete complex 5; all targets GO |
| Candidate pairs | 1,167 |
| Pairs on proteins with an experimental label (Price-148 or Swiss-Prot) | 3 |
| Known truth (positive via exact EC or ec2go; no GO negative derived) | 0 |
| Shared production wall time (three arms) | 800 s |

The three labelled pairs are Price proteins proposed for a different function (two L-arabinonate dehydratases,
EC 4.2.1.25, proposed for GO:0004456; one EC 2.6.1.57 protein proposed for GO:0004069): unconfirmed, not negative.
The unchanged EC wrapper `prepare_price_agent_truth.py` reports 0 pairs on this bundle because it reads the
candidates at the wrong nesting level; the plan's own join reads all 1,167 (both outputs kept). With no known truth
in either split the six methods and the agent were not run (pre-registered condition); correctness and harm remain
unidentifiable. The agent endpoint answered non-interactively (HTTP 200), so a future run needs truth, not access.

## Compute

| Step | Where | Time |
|---|---|---|
| DIAMOND claims, 22 genomes | ws, 4 threads | 13 min |
| Leakage search | ws, 16 threads | under 5 min |
| 22 cells × 14 arms | Shaheen, 2 workq nodes, 16 GB heap per cell | 31 min wall; slowest arm 3.1 min |
| Anomaly extractor, 3 arms × 22 genomes | Shaheen, 1 node | 13 min |

## Limits

* The predictor is a homology transfer from experimental Swiss-Prot labels of other genera; its EC top-1 on Price
  (49 / 148) is below CLEAN's. The result measures GSPA's learned stage on its training predictor, not an enzyme
  predictor comparison.
* The reference holds 25 entries of *S. meliloti* 1021 (fence gap); no Price-148 protein is among them and the
  sensitivity analysis without that genome leaves every direction unchanged. The fence defect itself is in
  `homology_baseline.py` (first-word genus) and affects the training reference too.
* No taxon constraint is applied in these single-genome commands (strain taxa not in the bundled hierarchy;
  enforcement not enabled), as in the row-5 cells.
* With actual gap observations the learned models abstain on every candidate; a model trained with gap context
  available would be needed for a learned reaction-context ablation.
* The row-1 split puts one genus in calibration; intervals are wide.
