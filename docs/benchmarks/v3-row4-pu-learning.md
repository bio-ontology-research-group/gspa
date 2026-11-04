# Positive-unlabelled training of the learned stage and of the no-candidate outcome (row 4)

Plan: `tasks/v3-row4-pu-learning-2026-11-02.md`, frozen 2026-11-02 15:01 (sha256 `9ef5576a…`, recorded in
`benchmark/v3/results/pu-learning-v1/PLAN_SHA256` before any fit); deviations appended, dated, in the same file.
Results folder: `benchmark/v3/results/pu-learning-v1/` (SHA256SUMS; files over 5 MB listed in `LARGE_FILES.tsv`
and kept in `/mnt/data2/leechuck/gspa-v3-data/pu-learning-v1/`). Code: `benchmark/v3/pu_learning/`, tests
`benchmark/v3/tests/test_pu_*.py` (written by agy-auto from the supervisor's spec, reviewed and tested by the
supervisor; the production trainer `scripts/train_context_model.py` and the JVM tree are unchanged).

## 1. Question and method

The complete-locus plan of the same day showed that closed-world completeness cannot supply labels for the
no-candidate ("none of these is right") outcome: 0 all-negative competitions and 6 contested usable ones in
941 complete competitions, because nearly every competition has unlabelled members. This work changes the
learning problem instead. A candidate row is positive when Swiss-Prot carries an experimental annotation of
that function on that protein, and **unlabelled** otherwise; it is never a negative. The candidate layer is
fitted with the unbiased positive-unlabelled risk of du Plessis et al. (2015) with the logistic loss, which is
convex for the linear model, and the non-negative correction of Kiryo et al. (2017), the estimator PU-GO uses
(Zhapa-Camacho, Tang, Kulmanov, Hoehndorf, Bioinformatics 2024, doi:10.1093/bioinformatics/btae237), is applied
only when the fitted negative-class risk goes below zero. PU-GO's ranking form (positives must outrank
unlabelled rows of the same genome and function type) is a second arm, and the same risk drives a tree arm.
The class prior π comes from the train split only (Elkan–Noto: labelled rate divided by the mean out-of-fold
score of labelled rows, c), and is swept. Probabilities are checked with the identity P(labelled | x) = c · P(true | x)
("PU reliability"); no Platt refit against closed-world zeros is made, only one temperature chosen on the
calibration split by that identity. The null layer keeps the production softmax form and is fitted by the PU
marginal likelihood: a competition with one labelled positive keeps its target, a competition with no label
contributes P(null) + (1 − c) Σ_j P(j). Experimental negatives (GOA NOT rows and accepted Swiss-Prot negations
from the complete-locus registry) enter evaluation only, as harm counts.

Everything else is the frozen recipe of the family-excluded training: contract 3 (160 features), the 11/3/4
genome split disjoint on genome, group, genus and family, family masking, group weighting, five grouped folds,
seed 42, and the frozen PN artifacts as the comparison arms. Both claim regimes are run (family-excluded,
primary; homology, which is the regime the Price cells use).

## 2. Data

| quantity | family-excluded | homology |
|---|---:|---:|
| rows (train / calibration / test) | 1,683 / 345 / 416 | 3,038 / 616 / 608 |
| positives | 183 / 29 / 34 | 394 / 72 / 75 |
| unlabelled (closed-world unconfirmed) | 1,499 / 316 / 382 | 2,643 / 544 / 533 |
| explicit negatives (kept out of every fit) | 1 / 0 / 0 | 1 / 0 / 0 |
| competitions with ≥ 2 members | 283 / 38 / 74 | 482 / 84 / 117 |
| of which no labelled member | 240 / 34 / 67 | 398 / 75 / 100 |
| of which exactly one labelled positive | 30 / 1 / 4 | 52 / 3 / 7 |
| of which two or more labelled positives (excluded from the null layer) | 13 / 3 / 3 | 32 / 6 / 10 |
| Elkan–Noto estimate on train OOF: labelled rate, c, π̂ | 0.125, 0.199, 0.626 | 0.146, 0.256, 0.568 |

The plan's §3.2 counted competitions with at least one positive; the single-positive counts above are the ones
the null layer can use, and both calibration splits fall short of the trainer's gate of five (deviation of
15:50). Experimental negatives available for evaluation: 7 accepted (1 cohort, 6 model-genome), of which 6 join
a scored candidate; P9WKI3 / GO:0003938 has no family-excluded claim.

## 3. Arms

A1/A2 are the frozen PN logistic and tree artifacts of each regime (never refitted). A3 is the PU logistic
model at the Elkan–Noto prior π̂, A4 the PU-GO ranking form, A5 the PU tree; S1 repeats A3 at π equal to the
labelled rate, 0.20 and 0.80; S2 adds up to ten unlabelled export candidates per table competition (4,605
train and 917 calibration rows on proteins without any experimental annotation) and re-estimates the prior;
S3 refits the PN and PU logistic models leaving one genus out at a time (18 folds). The non-negative correction
fired for A3 in both regimes (the unbiased fit's negative-class risk was −0.054 on family-excluded train OOF)
and for the family-excluded tree. All artifacts were hashed (`SHA256SUMS.models`) before any test row was read.

## 4. Family-excluded regime, cohort test split (4 genomes, 416 rows, 299 supported by every arm)

Supported precision counts labelled positives among calls at each arm's own frozen display threshold; the
PU-adjusted expected precision is the mean probability among the calls, the model's own claim about hidden
positives, and is not evidence. Negatives are the six accepted experimental negatives that join a scored
candidate (one cohort train-genome row, five model-genome rows; section 7), scored with every arm.

| arm | π | threshold | calls | labelled positives among calls (supported precision) | PU-adjusted expected precision | recall (34 test positives) | AUROC positive vs unlabelled | PU-ECE | PN-ECE (reference) | positives ranked first in their competition | negatives called / 6 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A1 PN logistic (frozen) | — | 0.298 | 4 | 1 (0.250) | 0.428 | 0.029 | 0.743 | 0.094 | 0.050 | 1.000 | 0 |
| A2 PN tree (frozen) | — | 0.216 | 42 | 12 (0.286) | 0.274 | 0.353 | 0.752 | 0.091 | 0.031 | 1.000 | 4 |
| A3 PU logistic, π̂ | 0.626 | 0.498 | 39 | 7 (0.179) | 0.715 | 0.206 | 0.712 | 0.081 | 0.209 | 1.000 | 4 |
| A4 PU-rank logistic, π̂ | 0.626 | 0.999 | 100 | 20 (0.200) | 0.999 | 0.588 | 0.738 | 0.092 | 0.892 | 1.000 | 6 |
| A5 PU tree, π̂ | 0.626 | 0.633 | 128 | 23 (0.180) | 0.857 | 0.676 | 0.682 | 0.029 | 0.527 | 0.250 | 6 |
| S1 PU logistic, π = labelled rate | 0.125 | 0.448 | 9 | 4 (0.444) | 0.484 | 0.118 | 0.739 | 0.063 | 0.244 | 1.000 | 0 |
| S1 PU logistic, π = 0.20 | 0.200 | 0.508 | 11 | 4 (0.364) | 0.587 | 0.118 | 0.740 | 0.063 | 0.200 | 1.000 | 0 |
| S1 PU logistic, π = 0.80 | 0.800 | 0.506 | 35 | 6 (0.171) | 0.694 | 0.176 | 0.713 | 0.057 | 0.336 | 1.000 | 4 |
| S2 PU logistic + unlabelled export pool | 0.381 | 0.988 | 13 | 5 (0.385) | 1.000 | 0.147 | 0.710 | 0.091 | 0.167 | 1.000 | 0 |

Pairwise contrasts on the 299 rows supported by both arms (corrections = positives called only by the first
arm; harms = experimental negatives called only by the first arm; unconfirmed rows called only by one arm are
listed, never counted as harms):

| contrast | corrections | harms (negatives) | unconfirmed called only by the first arm | both call |
|---|---:|---:|---:|---:|
| A3 (π̂) vs A1, test rows | 6 | 0 among test rows; 4 among the six negatives | 29 | 1 positive, 3 unconfirmed |
| A5 (π̂) vs A2, test rows | 11 | 0 among test rows; 2 more negatives than A2 | 75 | — |
| S1 π = 0.20 vs A1, test rows | 3 | 0; 0 negatives | 4 | — |

Reading. At the Elkan–Noto prior the PU model calls ten times as many rows as the PN logistic model and recovers
six more labelled positives, but its labelled precision falls below the PN model's and it calls four of the six
experimentally denied functions (probabilities 0.51–0.82; the two not called sit at 0.38 and 0.46). The model's own PU-adjusted precision of 0.72 says
that most of the 29 extra unconfirmed calls are hidden positives; the negatives say the prior is too high, as
§11 of the plan anticipated (the "selected completely at random" assumption fails in the direction of
over-estimating π for peripheral functions, and the estimator is biased upward when the classifier is weak).
At π = 0.20 or at the labelled rate the PU model is better than the PN model on every cohort criterion that
can be measured: higher labelled precision at four times the recall, lower PU-ECE, no experimental negative
called, three corrections and no harm. The PU-GO ranking form (A4) and the PU tree (A5) are not probability
models at their thresholds: they call a third to a half of all supported rows and every negative. The
leave-one-genus-out spread (18 folds, PU at π̂ minus PN) has a median PU-ECE difference of −0.038 (15 of 18
folds lower), a median AUROC difference of −0.006 (11 folds lower, 5 higher), a median supported-precision
difference of −0.009 and a median recall difference of +0.021; the PU model is better calibrated by the PU
identity and otherwise indistinguishable from the PN model across genera.

## 5. Null layer and deletion test (family-excluded)

At the trainer's gates the null layer does not train in either regime (one and three single-positive
calibration competitions against a gate of five), so P2 fails as written. The gate-deviation fits (calibration
gate lowered to the available count; deviation of 15:50) and S2 train a null layer; the table reports them as
secondary results.

| fit | calibration single-positive groups | calibration marginal NLL vs uniform | test marginal NLL vs uniform (≥ 2 members, fully supported groups) | mean P(null), no-label test groups | mean P(null), single-positive test groups | deletion test: competitions, AUROC of P(null), AUROC of naive null 1 − max p̂ |
|---|---:|---:|---:|---:|---:|---|
| A3 at the gates | 1 | not trained | — | — | — | 4, —, 0.813 (A1 naive 0.813) |
| A3 gate deviation | 1 | 0.152 vs 0.179 | 0.235 vs 0.241 (39 groups) | 0.499 | 0.496 | 4, 1.000, 0.813 |
| S2 unlabelled pool | 1 | 0.175 vs 0.175 | 0.143 vs 0.161 (39 groups) | 3 × 10⁻⁶ | 3 × 10⁻⁶ | 4, 1.000, 0.813 |

The gate-deviation null layer puts about one half on "none of these" for every competition and does not
separate no-label from single-positive groups; the S2 null layer collapsed to zero null probability (its
deletion AUROC of 1.0 rests on differences of 10⁻⁶). With four deletion competitions neither figure is
evidence; P3 is not met by a null layer that trains at the gates, and the honest statement is that the
family-excluded cohort does not have enough contested single-positive competitions in calibration or test to
train or to judge the null outcome, with or without PU learning.

## 6. Homology regime and model-genome stratum

### 6.1 Homology regime, cohort test split (4 genomes, 608 rows, 462 supported; 75 test positives)

| arm | π | threshold | calls | labelled positives among calls (supported precision) | PU-adjusted expected precision | recall | AUROC positive vs unlabelled | PU-ECE | PN-ECE (reference) | positives ranked first in their competition | negatives called / 1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A1 PN logistic (frozen) | — | 0.321 | 31 | 14 (0.452) | 0.540 | 0.187 | 0.756 | 0.132 | 0.075 | 0.667 | 0 |
| A2 PN tree (frozen) | — | 0.223 | 38 | 14 (0.368) | 0.329 | 0.187 | 0.720 | 0.134 | 0.074 | 0.500 | 0 |
| A3 PU logistic, π̂ | 0.568 | 0.957 | 40 | 15 (0.375) | 0.993 | 0.200 | 0.735 | 0.062 | 0.234 | 0.611 | 0 |
| A4 PU-rank logistic, π̂ | 0.568 | 0.999 | 122 | 38 (0.311) | 0.999 | 0.507 | 0.742 | 0.100 | 0.843 | 0.667 | 0 |
| A5 PU tree, π̂ | 0.568 | 0.995 | 48 | 17 (0.354) | 0.995 | 0.227 | 0.741 | 0.041 | 0.418 | 0.722 | 0 |
| S1 PU logistic, π = labelled rate | 0.146 | 0.512 | 30 | 13 (0.433) | 0.638 | 0.173 | 0.749 | 0.097 | 0.093 | 0.611 | 0 |
| S1 PU logistic, π = 0.20 | 0.200 | 0.674 | 27 | 12 (0.444) | 0.808 | 0.160 | 0.749 | 0.091 | 0.126 | 0.611 | 0 |
| S1 PU logistic, π = 0.80 | 0.800 | 0.645 | 33 | 12 (0.364) | 0.763 | 0.160 | 0.731 | 0.092 | 0.456 | 0.611 | 0 |

Contrasts on the 462 rows supported by
both arms: A3 against A1 corrects 3 and loses 2 labelled positives with no harm (14 and 6 unconfirmed rows
called by only one side); A5 against A2 corrects 10 and loses 7; S1 at π = 0.20 against A1 corrects 1 and loses
3. The only joined negative (the cohort ThrRS editing subunit) is below every arm's threshold. The
leave-one-genus-out spread (PU at π̂ minus PN, 18 folds) has a median supported-precision difference of −0.155
(13 folds lower), a median recall difference of 0.000, a median AUROC difference of −0.024 (13 lower) and a
median PU-ECE difference of −0.059 (16 lower). In this regime the PU model at π̂ is better calibrated by the PU
identity and worse on labelled precision; the lower priors recover the PN model's precision at slightly lower
recall. Nothing here distinguishes the arms beyond the intervals a 75-positive test split allows.

Null layer, homology: not trained at the gates (three single-positive calibration competitions). The
gate-deviation fit (calibration gate 3) trains: calibration marginal NLL 0.130 against 0.233 for the uniform
model (49 competitions of ≥ 2 members, 46 without a label), test marginal NLL 0.287 against 0.294 (70
competitions); mean P(null) 0.76 on no-label test groups and 0.77 on single-positive groups; deletion test on 7
competitions: AUROC of P(null) 0.796 against 0.673 for the naive null of the PN model (mean P(null) 0.767
intact, 0.857 after deleting the positive). The PU tree's null layer gives 0.735 against 0.561. These are the
first non-degenerate null-outcome figures in the project, but they rest on three calibration and seven test
competitions and on a gate that the plan lowered after the fact; they are reported, not counted.

### 6.2 Model-genome stratum (E. coli K-12, B. subtilis 168, M. tuberculosis H37Rv; held out, development-examined, never trained on)

Family-excluded contract claims (3 proteomes, contracts matching the family-excluded model), production feature
export on Shaheen without cohort inputs, 69,789 candidates, 47,781 (68 %) supported by every arm; labels from the
complete-locus registry's model set: 508 positive rows after is_a/part_of propagation (450 exact), 5 accepted
negatives joined (the sixth has no claim), no conflicts. Labelled positives are 0.7 % of rows here, so labelled
precision is small for every arm and only the relative figures matter.

| arm | calls | labelled positives among calls (supported precision) | PU-adjusted expected precision | recall (508 positives) | AUROC positive vs unlabelled | PU-ECE | negatives called / 5 |
|---|---:|---:|---:|---:|---:|---:|---:|
| A1 PN logistic | 1,340 | 52 (0.039) | 0.348 | 0.102 | 0.756 | 0.085 | 0 |
| A2 PN tree | 5,483 | 124 (0.023) | 0.266 | 0.244 | 0.724 | 0.083 | 4 |
| A3 PU logistic, π̂ | 10,209 | 243 (0.024) | 0.658 | 0.478 | 0.735 | 0.063 | 3 |
| A4 PU-rank | 17,202 | 326 (0.019) | 0.999 | 0.642 | 0.737 | 0.189 | 5 |
| A5 PU tree | 22,855 | 371 (0.016) | 0.878 | 0.730 | 0.712 | 0.122 | 5 |
| S1 PU logistic, π = 0.20 | 2,649 | 97 (0.037) | 0.534 | 0.191 | 0.754 | 0.058 | 0 |

The five negatives: the PN logistic model scores all of them at 0.12–0.24, below its threshold; the PU model at
π̂ scores them at 0.38–0.82 and calls three (the lipoate ligase activity denied for the octanoyltransferase LipM
at 0.82; the phosphoglycerate mutase activity denied for O07617 and LipM's GO term, both at 0.51); at π = 0.20
all five fall below the threshold (0.40–0.49). In the three competitions where the negative is the only
candidate, the gate-deviation null layer puts 0.64–0.65 on "none of these"; in the two where a labelled positive
competes with the negative (B. subtilis LipM against the true lipoate ligase), the positive ranks first in every
arm and P(null) is 0.45–0.48. Across the stratum's 6,233 competitions with ≥ 2 supported members the
gate-deviation null layer's marginal NLL is 0.145 against 0.174 for the uniform model.

## 7. JVM parity and production replay (P1)

The PU artifacts keep the artifact schema, so the unchanged JVM reader loads them. Two checks were run on
Shaheen with the release-candidate-3 engine: the table check (`check_context_mode_models.groovy`, unchanged,
every row of the family-excluded table) and a full production replay of the 18 cohort genomes
(`gspa integrate --batch-manifest` with cohort inference, CLADE hold-out and feature export, on the
contract-bearing family-excluded claims, `--learned-model` set to the PU artifact), compared with the Python
reference predictor on the archived export rows.

| check | artifact | rows compared | support status agreement | rows supported by both | max abs. probability difference | null probabilities compared (max abs. difference) |
|---|---|---:|---:|---:|---:|---|
| table rows | A3 PU logistic | 2,444 | 2,444 | 1,779 | 3.3 × 10⁻¹⁶ | — (null layer not trained) |
| table rows | A5 PU tree | 2,444 | 2,444 | 1,779 | 1.1 × 10⁻¹⁶ | — |
| table rows | A3 gate-deviation (null trained) | 2,444 | 2,444 | 1,779 | 3.3 × 10⁻¹⁶ | — (table rows are scored singly) |
| production, 18 genomes | A3 PU logistic | 253,743 | 253,743 | 179,623 | 3.3 × 10⁻¹⁶ | — |
| production, 18 genomes | A3 gate-deviation | 253,743 | 253,743 | 179,623 | 3.3 × 10⁻¹⁶ | 169,562 (2.2 × 10⁻¹⁶) |

The production learned stage reports the same support structure as the archived PN verification: 272,366
candidates per run, 196,621 scored (193,888 rescored, 2,733 scored without projection), 75,745 abstaining for
the same three reasons and counts (70,536 profile lacking both classes, 3,538 unseen profile, 1,671 outside the
frozen score range). With the gate-deviation artifact the JVM null layer ran on 185,063 candidates (median
P(null) 0.34, 26 % above 0.5). The 18,623 exported candidates that did not join the production report by key are exactly the EC candidates: the
engine writes `EC:` in front of the EC number in its candidate key and the archived export does not; after
removing the prefix the two key sets are identical for every genome (checked on GCA_000007005.1: 15,664 keys on
each side, 1,220 EC keys differing only by the prefix). The EC rows are covered numerically by the table check,
which holds EC rows and agrees on all of them. P1 is met: the unchanged reader loads the PU artifacts, the
probabilities agree to floating-point precision on every compared candidate, and the production support
structure is the one archived for the PN models.

## 8. Price benchmark (22 genomes, homology-regime artifacts, frozen row-1 genus split)

The other session's contract-bearing DIAMOND claims and cells (`price-contract-learned-v1`, read only) supply
the `raw`, `stage0`, `stage1` and PN learned arms; the PU arms were run with the unchanged cell runner and
wrapper on Shaheen (jobs 15768034, 15768035; 22 cells each, no failure) with the homology-regime PU artifacts.
The evaluator is the row-1 code (`row1_calibration`), the split is `split.json` (sha256 `95020999…`), the policy
was written and hashed (`cddf53f6…`) before any test metric. The WP5 calibrated-threshold rule (≥ 10 calls at
≥ 90 % reference support on the calibration genomes) finds no supported threshold for any arm, PN or PU, in
either split direction, exactly as for the PN arms in the other session's evaluation; the table therefore
reports top-1, the fixed 0.5 cut and each learned arm's own display threshold. Learned-stage coverage is
identical for all five learned arms (1,234,318 candidates, 919,227 scored, 315,091 abstained with the same
reasons), as the input-support gates are inherited from the training table.

| split direction | arm | reference proteins | top-1 matches (share, 95 % CI) | top-3 matches | display threshold | display calls (supported / unconfirmed) | positive-pair Brier | top-call ECE |
|---|---|---:|---|---:|---:|---|---:|---:|
| primary (15 test genomes) | stage1 | 70 | 22 (0.271, 0.198–0.388) | 28 | — | — | 0.813 | 0.061 |
| primary | learned PN logistic | 70 | 22 (0.271, 0.198–0.388) | 28 | 0.321 | 30 (11 / 19) | 0.806 | 0.067 |
| primary | learned PN tree | 70 | 21 (0.251, 0.176–0.367) | 26 | 0.223 | 31 (10 / 21) | 0.879 | 0.116 |
| primary | learned PU logistic (π̂ 0.568) | 70 | 22 (0.271, 0.198–0.388) | 28 | 0.957 | 18 (7 / 11) | 0.617 | 0.456 |
| primary | learned PU tree | 70 | 21 (0.238, 0.164–0.357) | 28 | 0.995 | 26 (10 / 16) | 0.637 | 0.461 |
| primary | learned PU-rank | 70 | 22 (0.271, 0.198–0.388) | 28 | 0.999 | 19 (9 / 10) | 0.603 | 0.599 |
| swapped (7 Pseudomonas genomes as test) | stage1 | 78 | 27 (0.333, 0.230–0.430) | 37 | — | — | 0.758 | 0.177 |
| swapped | learned PN logistic | 78 | 25 (0.321, 0.209–0.424) | 37 | 0.321 | 38 (13 / 25) | 0.765 | 0.181 |
| swapped | learned PU logistic | 78 | 25 (0.321, 0.209–0.424) | 37 | 0.957 | 18 (4 / 14) | 0.522 | 0.426 |
| swapped | learned PU tree | 78 | 15 (0.208, 0.128–0.287) | 37 | 0.995 | 19 (6 / 13) | 0.524 | 0.552 |

Contrasts at the fixed 0.5 cut: PU logistic corrects 24 (primary) and 36 (swapped) reference pairs relative
to the PN logistic model and harms none, because its probabilities are higher throughout; at the arms' own
display thresholds it calls fewer reference pairs than the PN model (4 and 9 supported calls lost, none gained).
Top-1 ranking is unchanged for the logistic arms and worse for the PU tree in the swapped direction (12 top-1
lost, 1 gained). The PU probabilities are not calibrated against labelled reference pairs (top-call ECE 0.43–0.60
against 0.06–0.18): they estimate the probability that the function is true, which under the PU model is
several times the probability that it is labelled, and the Price evaluator measures against labels. Without
GCA_000006965.1 (the Sinorhizobium reference leak noted by the Price plan) every figure moves within its interval.
Price adds no evidence for or against the null outcome (no competition there has a label for every member),
and it confirms the cohort finding that the Elkan–Noto prior sharpens probabilities beyond what labels support.

## 9. Verdict against the frozen criteria

| id | criterion (short) | family-excluded | homology | status |
|---|---|---|---|---|
| P0 | plan and artifacts hashed before any test read | plan sha256 recorded 15:01; `SHA256SUMS.models` before evaluate | same | PASS |
| P1 | unchanged JVM reader; table and production parity | 2,444 rows and 253,743 production candidates agree to ≤ 3.3 × 10⁻¹⁶; support structure identical to the PN archive | not run (same reader, same artifact schema) | PASS |
| P2 | null layer trains at the default gates, calibration NLL below uniform | not trained (1 single-positive calibration competition against a gate of 5) | not trained (3 against 5) | FAIL |
| P3 | deletion AUROC ≥ 0.70 and above the PN naive null, for a gated null layer | no gated null layer; 4 competitions | no gated null layer; gate deviation 0.796 vs 0.673 on 7 competitions, reported only | FAIL |
| P4 | supported precision of A3 ≥ A1 at ≥ recall | 0.179 vs 0.250 (recall 0.206 vs 0.029) | 0.375 vs 0.452 (0.200 vs 0.187) | FAIL |
| P5 | A3 calls no more experimental negatives than A1 | 4 of 6 vs 0 of 6 | 0 of 1 vs 0 of 1 (stratum: 3 of 5 vs 0 of 5) | FAIL |
| P6 | PU-ECE of A3 ≤ A1 on test | 0.081 vs 0.094 | 0.062 vs 0.132 | PASS |
| P7 | P3–P5 hold at π̂ and at π = 0.20; LOGO spread reported | at π = 0.20: precision 0.364 vs 0.250, recall 0.118 vs 0.029, 0 negatives called; at π̂ the opposite | at π = 0.20: 0.444 vs 0.452, recall 0.160 vs 0.187 | FAIL (the direction depends on the prior) |
| P8 | Price scored under the frozen split | homology artifacts in 22 cells, all arms scored, no criterion attached | — | PASS (reported) |

What the measurements say. (1) PU learning is implementable inside the existing artifact and engine with no
JVM change, and its probabilities survive the PU reliability identity better than the PN models' in both
regimes (P6). (2) The Elkan–Noto prior (0.63 and 0.57) is too high: at that prior the PU logistic model calls
experimentally denied functions that the PN model does not (four of six cohort-joined negatives, three of five
in the model genomes), and its labelled precision falls below the PN model's. The plan named this risk; the
negatives, which never entered any fit, are what showed it. (3) At a low prior (0.20, or the labelled rate) the
PU logistic model is better than the PN logistic model on every family-excluded criterion that can be measured
(precision 0.36–0.44 against 0.25 at four times the recall, lower PU-ECE, no negative called, three corrections
and no harm) and equal to it in the homology regime; but π = 0.20 was one of four pre-registered sweep values,
and choosing it from these results would select on the test negatives. A prior fixed from the negatives would
need a negative set that is not also the harm count; with seven accepted negatives in all, that split does not
exist yet. (4) The ranking form and the tree are not probability models at their thresholds here. (5) The
no-candidate outcome still cannot be trained or judged at the gates the plan kept: one and three contested
single-positive calibration competitions. With the gate lowered, the PU marginal likelihood trains a null
layer that beats the uniform model on calibration and test and separates deleted-positive competitions from
intact ones (0.796 on seven competitions), and in the model genomes puts 0.64–0.65 on "none of these" where
the only candidate is an experimentally denied one. That is the first non-degenerate null-outcome signal in the
project, and it is too small to count.

Row 4 stays Not met. What this round adds to the row: an honest, reproducible way to train and to measure the
null outcome without complete locus labels; a measured answer that the learned stage's precision is bounded by
the prior, not by the loss; and a harm count against experimental negatives for every arm. What would move it:
more accepted experimental negatives (the complete-locus review queue holds 22; GOA and Swiss-Prot hold
hundreds for the three model genomes) so that the prior can be chosen on one negative set and the harm counted
on another; and the contested single-positive competitions in calibration that only labels on competing
candidates can give.

## 10. Reproduction

Stages, in order (`benchmark/v3/pu_learning/run_pu_learning.py`; work folder `--work`): `prepare`, `train
--folds 5 --seed 42` (family-excluded with `--export-dir` for the unlabelled pool), `evaluate` (with the same
export for the negatives), `report`; the gate-deviation fits use `--min-cal-pos-groups 1` (family-excluded) and
`3` (homology) with `--skip-s3`. Model-genome claims: `results/pu-learning-v1/code/predict_fe_contract_claims.py`;
feature export and production replays: `code/pu-model-stratum-export.sbatch`, `code/pu-cohort-replay.sbatch`,
`code/build_pu_cohort_inputs.py` (Shaheen, k10018); stratum scoring `pu_stratum_evaluate.py`; Price
`pu_price_evaluate.py` on the other session's PN cells and on the PU cells run with its unchanged wrapper. Tests:
24 in `benchmark/v3/tests/test_pu_*.py`, and the 41 trainer tests unchanged. The results folder holds 983 hashed
files (`SHA256SUMS`), 72 files over 5 MB listed in `LARGE_FILES.tsv` and kept in the work folder.

