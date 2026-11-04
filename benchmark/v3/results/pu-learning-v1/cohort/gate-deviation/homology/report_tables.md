# PU Learning Evaluation Tables — Regime `homology`

## 1. Candidate-Layer Test Metrics (Main Arms)

| Arm | Kind | Threshold | Supported Precision | PU-Adjusted Precision | Recall | AUROC (Pos vs Unl) | PU-ECE | PN-ECE (Ref) | Calls | Abstained |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A1 | logistic | 0.3212 | 0.4516 | 0.5401 | 0.1867 | 0.7556 | 0.1316 | 0.0748 | 31 | 146 |
| A2 | gbdt | 0.2232 | 0.3684 | 0.3293 | 0.1867 | 0.7195 | 0.1337 | 0.0744 | 38 | 146 |
| A3 | logistic | 0.9566 | 0.3750 | 0.9935 | 0.2000 | 0.7346 | 0.0618 | 0.2342 | 40 | 146 |
| A4 | logistic | 0.9988 | 0.3115 | 0.9990 | 0.5067 | 0.7416 | 0.1000 | 0.8428 | 122 | 146 |
| A5 | gbdt | 0.9951 | 0.3542 | 0.9951 | 0.2267 | 0.7408 | 0.0406 | 0.4177 | 48 | 146 |

## 2. Pairwise Contrasts on Test Split

| Contrast (A vs B) | Both Supported | Corrections (A only) | Corrections (B only) | Harms (A only) | Harms (B only) | Unconfirmed (A only) | Unconfirmed (B only) |
|---|---:|---:|---:|---:|---:|---:|---:|
| A3_vs_A1_test | 462 | 3 | 2 | 0 | 0 | 14 | 6 |
| A3_vs_A2_test | 462 | 6 | 5 | 0 | 0 | 16 | 15 |
| A4_vs_A1_test | 462 | 28 | 4 | 0 | 0 | 68 | 1 |
| A4_vs_A2_test | 462 | 25 | 1 | 0 | 0 | 62 | 2 |
| A5_vs_A1_test | 462 | 9 | 6 | 0 | 0 | 21 | 7 |
| A5_vs_A2_test | 462 | 10 | 7 | 0 | 0 | 22 | 15 |
| S1_pi_0.1456_vs_A1_test | 462 | 1 | 2 | 0 | 0 | 4 | 4 |
| S1_pi_0.1456_vs_A2_test | 462 | 4 | 5 | 0 | 0 | 6 | 13 |
| S1_pi_0.20_vs_A1_test | 462 | 1 | 3 | 0 | 0 | 3 | 5 |
| S1_pi_0.20_vs_A2_test | 462 | 4 | 6 | 0 | 0 | 5 | 14 |
| S1_pi_0.5681_vs_A1_test | 462 | 3 | 2 | 0 | 0 | 14 | 6 |
| S1_pi_0.5681_vs_A2_test | 462 | 6 | 5 | 0 | 0 | 16 | 15 |
| S1_pi_0.8000_vs_A1_test | 462 | 3 | 5 | 0 | 0 | 14 | 10 |
| S1_pi_0.8000_vs_A2_test | 462 | 6 | 8 | 0 | 0 | 16 | 19 |

## 3. Experimental Negatives Specificity

- Accepted negatives in registry + table: 7
- Joined and scored rows: 1
- Unjoined negatives: 6

| Arm | Threshold | Supported Negatives | Below Threshold | Called | Specificity | Probabilities |
|---|---:|---:|---:|---:|---:|---|
| A1 | 0.3212 | 1 | 1 | 0 | 1.0000 | 0.0860 |
| A2 | 0.2232 | 1 | 1 | 0 | 1.0000 | 0.0702 |
| A3 | 0.9566 | 1 | 1 | 0 | 1.0000 | 0.4163 |
| A4 | 0.9988 | 1 | 1 | 0 | 1.0000 | 0.9988 |
| A5 | 0.9951 | 1 | 1 | 0 | 1.0000 | 0.9126 |

## 4. Null Layer (Gates and Marginal NLL vs Baselines)

| Arm | Split | Trained | Groups (>=2) | Single-Pos | No-Label | Marginal NLL | Uniform Marginal NLL | Mean log(n+1) | Mean P(null) Single-Pos | Mean P(null) No-Label | Spearman P(null) vs log1p(n) |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A1 | calibration | false | 49 | 3 | 46 | — | 0.2477 | 1.1501 | — | — | — |
| A1 | test | false | 70 | 7 | 63 | — | 0.2937 | 1.1962 | — | — | — |
| A2 | calibration | false | 49 | 3 | 46 | — | 0.2477 | 1.1501 | — | — | — |
| A2 | test | false | 70 | 7 | 63 | — | 0.2937 | 1.1962 | — | — | — |
| A3 | calibration | true | 49 | 3 | 46 | 0.1764 | 0.2477 | 1.1501 | 0.6447 | 0.7597 | -0.2538 |
| A3 | test | true | 70 | 7 | 63 | 0.2871 | 0.2937 | 1.1962 | 0.7675 | 0.7644 | -0.5416 |
| A4 | calibration | true | 49 | 3 | 46 | 0.3161 | 0.2477 | 1.1501 | 0.0000 | 0.0000 | -0.2487 |
| A4 | test | true | 70 | 7 | 63 | 0.3416 | 0.2937 | 1.1962 | 0.0000 | 0.0000 | -0.4560 |
| A5 | calibration | true | 49 | 3 | 46 | 0.1670 | 0.2477 | 1.1501 | 0.4769 | 0.6864 | -0.2003 |
| A5 | test | true | 70 | 7 | 63 | 0.2395 | 0.2937 | 1.1962 | 0.4718 | 0.6218 | -0.3667 |

## 5. Deletion Test (Single-Positive Test Competitions, >=2 Members)

| Arm | Locus Trained | Competitions | Mean P(null) Intact | Mean P(null) Deleted | AUROC P(null) | Mean Naive Null Intact | Mean Naive Null Deleted | AUROC Naive Null |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| A1 | false | 7 | — | — | — | 0.8503 | 0.8892 | 0.6735 |
| A2 | false | 7 | — | — | — | 0.8263 | 0.8364 | 0.5408 |
| A3 | true | 7 | 0.7675 | 0.8572 | 0.7959 | 0.5288 | 0.6465 | 0.6735 |
| A4 | true | 7 | 0.0000 | 0.0000 | 0.7959 | 0.0012 | 0.0012 | 0.6735 |
| A5 | true | 7 | 0.4718 | 0.6079 | 0.7347 | 0.1440 | 0.1790 | 0.5612 |

## 6. S1 Prior Sweep and S2 Unlabelled Pool

| Arm | Role | pi | Supported Precision | PU-Adjusted Precision | Recall | AUROC (Pos vs Unl) | PU-ECE | Deletion AUROC P(null) | Negatives Called |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| S1_pi_0.1456 | labelled_rate | 0.1456 | 0.4333 | 0.6378 | 0.1733 | 0.7489 | 0.0971 | 0.7959 | 0 |
| S1_pi_0.20 | 0.20 | 0.2000 | 0.4444 | 0.8084 | 0.1600 | 0.7487 | 0.0908 | 0.7959 | 0 |
| S1_pi_0.5681 | pi_hat | 0.5681 | 0.3750 | 0.9935 | 0.2000 | 0.7346 | 0.0618 | 0.7959 | 0 |
| S1_pi_0.8000 | min_2pi_hat_0.8 | 0.8000 | 0.3636 | 0.7628 | 0.1600 | 0.7310 | 0.0925 | 0.7959 | 0 |

## 7. S3 Leave-One-Genus-Out (LOGO) Spread (PU minus PN)

| Metric | Folds | Median (PU - PN) | Min (PU - PN) | Max (PU - PN) | Folds > 0 | Folds = 0 | Folds < 0 |
|---|---:|---:|---:|---:|---:|---:|---:|
| auroc_positive_vs_unlabelled | 0 | — | — | — | 0 | 0 | 0 |
| pu_ece | 0 | — | — | — | 0 | 0 | 0 |
| pu_expected_precision | 0 | — | — | — | 0 | 0 | 0 |
| recall | 0 | — | — | — | 0 | 0 | 0 |
| supported_precision | 0 | — | — | — | 0 | 0 | 0 |

## 8. Pass Criteria Checklist (P0–P8)

| Criterion | Measured Value | Status |
|---|---|---|
| P0 (Plan & artifact hashes frozen before test) | plan_sha256=True, SHA256SUMS.models=present | PASS |
| P1 (JVM reader & production genome replay parity) | Python artifact schema/calibration verified; JVM replay external | PENDING_JVM_REPLAY |
| P2 (Null layer trained & cal marginal NLL < uniform) | trained=True, cal_nll=0.1764 vs uniform=0.2477 | PASS |
| P3 (Deletion test AUROC >= 0.70 and > A1 naive null) | A3_p_null=0.7959, A3_naive=0.6735, A1_naive=0.6735 | PASS |
| P4 (Supported precision A3 >= A1 at >= recall) | A3_prec=0.3750, A1_prec=0.4516, A3_rec=0.2000, A1_rec=0.1867 | FAIL |
| P5 (Harms on negatives A3 <= A1) | A3_called=0, A1_called=0 | PASS |
| P6 (Test PU-ECE A3 <= A1) | A3_pu_ece=0.0618, A1_pu_ece=0.1316 | PASS |
| P7 (Robustness at pi_hat and pi=0.20 + LOGO reported) | pi_hat(P3=True,P4=False,P5=True), pi_0.20(P3=True,P4=False,P5=True) | FAIL |
| P8 (Price-148 evaluation under frozen Price split) | External to cohort table (awaits frozen Price Step C claims) | PENDING_PRICE |
