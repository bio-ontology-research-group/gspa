# PU Learning Evaluation Tables — Regime `family_excluded`

## 1. Candidate-Layer Test Metrics (Main Arms)

| Arm | Kind | Threshold | Supported Precision | PU-Adjusted Precision | Recall | AUROC (Pos vs Unl) | PU-ECE | PN-ECE (Ref) | Calls | Abstained |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A1 | logistic | 0.2982 | 0.2500 | 0.4280 | 0.0294 | 0.7429 | 0.0944 | 0.0501 | 4 | 117 |
| A2 | gbdt | 0.2158 | 0.2857 | 0.2738 | 0.3529 | 0.7517 | 0.0914 | 0.0314 | 42 | 117 |
| A3 | logistic | 0.4981 | 0.1795 | 0.7150 | 0.2059 | 0.7119 | 0.0812 | 0.2087 | 39 | 117 |
| A4 | logistic | 0.9988 | 0.2000 | 0.9990 | 0.5882 | 0.7379 | 0.0919 | 0.8916 | 100 | 117 |
| A5 | gbdt | 0.6326 | 0.1797 | 0.8568 | 0.6765 | 0.6824 | 0.0292 | 0.5275 | 128 | 117 |

## 2. Pairwise Contrasts on Test Split

| Contrast (A vs B) | Both Supported | Corrections (A only) | Corrections (B only) | Harms (A only) | Harms (B only) | Unconfirmed (A only) | Unconfirmed (B only) |
|---|---:|---:|---:|---:|---:|---:|---:|
| A3_vs_A1_test | 299 | 6 | 0 | 0 | 0 | 29 | 0 |
| A3_vs_A2_test | 299 | 1 | 6 | 0 | 0 | 17 | 15 |
| A4_vs_A1_test | 299 | 19 | 0 | 0 | 0 | 78 | 1 |
| A4_vs_A2_test | 299 | 8 | 0 | 0 | 0 | 50 | 0 |
| A5_vs_A1_test | 299 | 22 | 0 | 0 | 0 | 102 | 0 |
| A5_vs_A2_test | 299 | 11 | 0 | 0 | 0 | 75 | 0 |
| S1_pi_0.1246_vs_A1_test | 299 | 3 | 0 | 0 | 0 | 2 | 0 |
| S1_pi_0.1246_vs_A2_test | 299 | 0 | 8 | 0 | 0 | 2 | 27 |
| S1_pi_0.20_vs_A1_test | 299 | 3 | 0 | 0 | 0 | 4 | 0 |
| S1_pi_0.20_vs_A2_test | 299 | 0 | 8 | 0 | 0 | 4 | 27 |
| S1_pi_0.6258_vs_A1_test | 299 | 6 | 0 | 0 | 0 | 29 | 0 |
| S1_pi_0.6258_vs_A2_test | 299 | 1 | 6 | 0 | 0 | 17 | 15 |
| S1_pi_0.8000_vs_A1_test | 299 | 5 | 0 | 0 | 0 | 26 | 0 |
| S1_pi_0.8000_vs_A2_test | 299 | 1 | 7 | 0 | 0 | 17 | 18 |

## 3. Experimental Negatives Specificity

- Accepted negatives in registry + table: 7
- Joined and scored rows: 1
- Unjoined negatives: 6

| Arm | Threshold | Supported Negatives | Below Threshold | Called | Specificity | Probabilities |
|---|---:|---:|---:|---:|---:|---|
| A1 | 0.2982 | 1 | 1 | 0 | 1.0000 | 0.2166 |
| A2 | 0.2158 | 1 | 1 | 0 | 1.0000 | 0.1805 |
| A3 | 0.4981 | 1 | 0 | 1 | 0.0000 | 0.5099 |
| A4 | 0.9988 | 1 | 0 | 1 | 0.0000 | 0.9991 |
| A5 | 0.6326 | 1 | 0 | 1 | 0.0000 | 0.9405 |

## 4. Null Layer (Gates and Marginal NLL vs Baselines)

| Arm | Split | Trained | Groups (>=2) | Single-Pos | No-Label | Marginal NLL | Uniform Marginal NLL | Mean log(n+1) | Mean P(null) Single-Pos | Mean P(null) No-Label | Spearman P(null) vs log1p(n) |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A1 | calibration | false | 16 | 1 | 15 | — | 0.2034 | 1.1166 | — | — | — |
| A1 | test | false | 39 | 4 | 35 | — | 0.2413 | 1.1117 | — | — | — |
| A2 | calibration | false | 16 | 1 | 15 | — | 0.2034 | 1.1166 | — | — | — |
| A2 | test | false | 39 | 4 | 35 | — | 0.2413 | 1.1117 | — | — | — |
| A3 | calibration | true | 16 | 1 | 15 | 0.1872 | 0.2034 | 1.1166 | 0.5079 | 0.4797 | -0.3643 |
| A3 | test | true | 39 | 4 | 35 | 0.2349 | 0.2413 | 1.1117 | 0.4965 | 0.4990 | -0.2747 |
| A4 | calibration | true | 16 | 1 | 15 | 0.2229 | 0.2034 | 1.1166 | 0.0000 | 0.0000 | -0.4204 |
| A4 | test | true | 39 | 4 | 35 | 0.2641 | 0.2413 | 1.1117 | 0.0000 | 0.0000 | -0.1735 |
| A5 | calibration | true | 16 | 1 | 15 | 0.2214 | 0.2034 | 1.1166 | 0.7205 | 0.5066 | -0.4313 |
| A5 | test | true | 39 | 4 | 35 | 0.2165 | 0.2413 | 1.1117 | 0.5545 | 0.6213 | -0.1490 |

## 5. Deletion Test (Single-Positive Test Competitions, >=2 Members)

| Arm | Locus Trained | Competitions | Mean P(null) Intact | Mean P(null) Deleted | AUROC P(null) | Mean Naive Null Intact | Mean Naive Null Deleted | AUROC Naive Null |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| A1 | false | 4 | — | — | — | 0.9269 | 0.9439 | 0.8125 |
| A2 | false | 4 | — | — | — | 0.9344 | 0.9666 | 0.8125 |
| A3 | true | 4 | 0.4965 | 0.6621 | 1.0000 | 0.6906 | 0.7392 | 0.8125 |
| A4 | true | 4 | 0.0000 | 0.0000 | 1.0000 | 0.0012 | 0.0013 | 0.8125 |
| A5 | true | 4 | 0.5545 | 0.7774 | 1.0000 | 0.2904 | 0.3653 | 0.5312 |

## 6. S1 Prior Sweep and S2 Unlabelled Pool

| Arm | Role | pi | Supported Precision | PU-Adjusted Precision | Recall | AUROC (Pos vs Unl) | PU-ECE | Deletion AUROC P(null) | Negatives Called |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| S1_pi_0.1246 | labelled_rate | 0.1246 | 0.4444 | 0.4842 | 0.1176 | 0.7389 | 0.0633 | 1.0000 | 0 |
| S1_pi_0.20 | 0.20 | 0.2000 | 0.3636 | 0.5867 | 0.1176 | 0.7402 | 0.0626 | 1.0000 | 0 |
| S1_pi_0.6258 | pi_hat | 0.6258 | 0.1795 | 0.7150 | 0.2059 | 0.7119 | 0.0812 | 1.0000 | 1 |
| S1_pi_0.8000 | min_2pi_hat_0.8 | 0.8000 | 0.1714 | 0.6940 | 0.1765 | 0.7130 | 0.0569 | 1.0000 | 1 |

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
| P2 (Null layer trained & cal marginal NLL < uniform) | trained=True, cal_nll=0.1872 vs uniform=0.2034 | PASS |
| P3 (Deletion test AUROC >= 0.70 and > A1 naive null) | A3_p_null=1.0000, A3_naive=0.8125, A1_naive=0.8125 | PASS |
| P4 (Supported precision A3 >= A1 at >= recall) | A3_prec=0.1795, A1_prec=0.2500, A3_rec=0.2059, A1_rec=0.0294 | FAIL |
| P5 (Harms on negatives A3 <= A1) | A3_called=1, A1_called=0 | FAIL |
| P6 (Test PU-ECE A3 <= A1) | A3_pu_ece=0.0812, A1_pu_ece=0.0944 | PASS |
| P7 (Robustness at pi_hat and pi=0.20 + LOGO reported) | pi_hat(P3=True,P4=False,P5=False), pi_0.20(P3=True,P4=True,P5=True) | FAIL |
| P8 (Price-148 evaluation under frozen Price split) | External to cohort table (awaits frozen Price Step C claims) | PENDING_PRICE |
