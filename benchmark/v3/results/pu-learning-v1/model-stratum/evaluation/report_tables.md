# Model-Genome Stratum Evaluation Tables (Row 4, §4.3, P5)

- Stratum classification: **held-out, development-examined; never trained on**
- Genomes: GCF_000005845.2, GCF_000009045.1, GCF_000195955.2
- Total export rows: 69789
- Kept rows: 69789 (positives: 508 [exact: 450, propagated ancestors: 58], negatives: 5, unlabelled: 69276)
- Excluded positive/negative conflicts: 0
- Unjoined accepted model negatives: 1

## 1. Pooled Candidate-Layer & Negative-Specificity Metrics

| Arm | Kind | c (Note) | Threshold | Supported | Abstained | Calls | Pos Called | Unconf Called | Neg Called | Supported Prec | PU Expected Prec | Recall | AUROC (Pos vs Unl) | PU-ECE | Neg Specificity |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A1 | logistic | 1.0000 (PU) | 0.2982 | 47781 | 22008 | 1340 | 52 | 1288 | 0 | 0.0388 | 0.3476 | 0.1024 | 0.7556 | 0.0846 | 1.0000 (5/5) |
| A2 | gbdt | 1.0000 (PU) | 0.2158 | 47781 | 22008 | 5483 | 124 | 5355 | 4 | 0.0226 | 0.2657 | 0.2441 | 0.7242 | 0.0832 | 0.2000 (1/5) |
| A3 | logistic | 0.1992 (PU) | 0.4981 | 47781 | 22008 | 10209 | 243 | 9963 | 3 | 0.0238 | 0.6584 | 0.4783 | 0.7350 | 0.0633 | 0.4000 (2/5) |
| A4 | logistic | 0.1992 (PU) | 0.9988 | 47781 | 22008 | 17202 | 326 | 16871 | 5 | 0.0190 | 0.9990 | 0.6417 | 0.7373 | 0.1893 | 0.0000 (0/5) |
| A5 | gbdt | 0.1992 (PU) | 0.6326 | 47781 | 22008 | 22855 | 371 | 22479 | 5 | 0.0162 | 0.8781 | 0.7303 | 0.7122 | 0.1218 | 0.0000 (0/5) |
| A3_gate | logistic | 0.1992 (PU) | 0.4981 | 47781 | 22008 | 10209 | 243 | 9963 | 3 | 0.0238 | 0.6584 | 0.4783 | 0.7350 | 0.0633 | 0.4000 (2/5) |
| S1_pi_0.20 | logistic | 0.1992 (PU) | 0.5083 | 47781 | 22008 | 2649 | 97 | 2552 | 0 | 0.0366 | 0.5336 | 0.1909 | 0.7537 | 0.0580 | 1.0000 (5/5) |

## 2. Per-Genome Candidate-Layer & Negative-Specificity Metrics

| Genome | Arm | Supported | Abstained | Positives (Sup/All) | Negatives (Sup/All) | Calls | Pos Called | Neg Called | Supported Prec | Recall | AUROC | Neg Specificity |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GCF_000005845.2 | A1 | 17413 | 6210 | 337/368 | 1/1 | 440 | 33 | 0 | 0.0750 | 0.0897 | 0.7423 | 1.0000 |
| GCF_000005845.2 | A2 | 17413 | 6210 | 337/368 | 1/1 | 1982 | 78 | 1 | 0.0394 | 0.2120 | 0.7142 | 0.0000 |
| GCF_000005845.2 | A3 | 17413 | 6210 | 337/368 | 1/1 | 3623 | 167 | 0 | 0.0461 | 0.4538 | 0.7186 | 1.0000 |
| GCF_000005845.2 | A4 | 17413 | 6210 | 337/368 | 1/1 | 6394 | 240 | 1 | 0.0375 | 0.6522 | 0.7310 | 0.0000 |
| GCF_000005845.2 | A5 | 17413 | 6210 | 337/368 | 1/1 | 8390 | 262 | 1 | 0.0312 | 0.7120 | 0.6982 | 0.0000 |
| GCF_000005845.2 | A3_gate | 17413 | 6210 | 337/368 | 1/1 | 3623 | 167 | 0 | 0.0461 | 0.4538 | 0.7186 | 1.0000 |
| GCF_000005845.2 | S1_pi_0.20 | 17413 | 6210 | 337/368 | 1/1 | 865 | 66 | 0 | 0.0763 | 0.1793 | 0.7394 | 1.0000 |
| GCF_000009045.1 | A1 | 13957 | 8633 | 24/30 | 3/3 | 407 | 6 | 0 | 0.0147 | 0.2000 | 0.8311 | 1.0000 |
| GCF_000009045.1 | A2 | 13957 | 8633 | 24/30 | 3/3 | 1629 | 7 | 2 | 0.0043 | 0.2333 | 0.7140 | 0.3333 |
| GCF_000009045.1 | A3 | 13957 | 8633 | 24/30 | 3/3 | 2873 | 18 | 3 | 0.0063 | 0.6000 | 0.8564 | 0.0000 |
| GCF_000009045.1 | A4 | 13957 | 8633 | 24/30 | 3/3 | 5047 | 15 | 3 | 0.0030 | 0.5000 | 0.6827 | 0.0000 |
| GCF_000009045.1 | A5 | 13957 | 8633 | 24/30 | 3/3 | 6517 | 21 | 3 | 0.0032 | 0.7000 | 0.7506 | 0.0000 |
| GCF_000009045.1 | A3_gate | 13957 | 8633 | 24/30 | 3/3 | 2873 | 18 | 3 | 0.0063 | 0.6000 | 0.8564 | 0.0000 |
| GCF_000009045.1 | S1_pi_0.20 | 13957 | 8633 | 24/30 | 3/3 | 767 | 9 | 0 | 0.0117 | 0.3000 | 0.8303 | 1.0000 |
| GCF_000195955.2 | A1 | 16411 | 7165 | 98/110 | 1/1 | 493 | 13 | 0 | 0.0264 | 0.1182 | 0.7941 | 1.0000 |
| GCF_000195955.2 | A2 | 16411 | 7165 | 98/110 | 1/1 | 1872 | 39 | 1 | 0.0208 | 0.3545 | 0.7679 | 0.0000 |
| GCF_000195955.2 | A3 | 16411 | 7165 | 98/110 | 1/1 | 3713 | 58 | 0 | 0.0156 | 0.5273 | 0.7679 | 1.0000 |
| GCF_000195955.2 | A4 | 16411 | 7165 | 98/110 | 1/1 | 5761 | 71 | 1 | 0.0123 | 0.6455 | 0.7812 | 0.0000 |
| GCF_000195955.2 | A5 | 16411 | 7165 | 98/110 | 1/1 | 7948 | 88 | 1 | 0.0111 | 0.8000 | 0.7541 | 0.0000 |
| GCF_000195955.2 | A3_gate | 16411 | 7165 | 98/110 | 1/1 | 3713 | 58 | 0 | 0.0156 | 0.5273 | 0.7679 | 1.0000 |
| GCF_000195955.2 | S1_pi_0.20 | 16411 | 7165 | 98/110 | 1/1 | 1017 | 22 | 0 | 0.0216 | 0.2000 | 0.7915 | 1.0000 |

## 3. Pooled Pairwise Contrasts (PU vs PN)

| Contrast (A vs B) | Both Supported | Corrections (A only) | Corrections (B only) | Harms on Negatives (A only) | Harms on Negatives (B only) | Unconfirmed (A only) | Unconfirmed (B only) |
|---|---:|---:|---:|---:|---:|---:|---:|
| A3_vs_A1 | 47781 | 191 | 0 | 3 | 0 | 8675 | 0 |
| A3_vs_A2 | 47781 | 160 | 41 | 1 | 2 | 7235 | 2627 |
| A4_vs_A1 | 47781 | 274 | 0 | 5 | 0 | 15584 | 1 |
| A4_vs_A2 | 47781 | 202 | 0 | 1 | 0 | 11591 | 75 |
| A5_vs_A1 | 47781 | 319 | 0 | 5 | 0 | 21191 | 0 |
| A5_vs_A2 | 47781 | 247 | 0 | 1 | 0 | 17124 | 0 |
| A3_gate_vs_A1 | 47781 | 191 | 0 | 3 | 0 | 8675 | 0 |
| A3_gate_vs_A2 | 47781 | 160 | 41 | 1 | 2 | 7235 | 2627 |
| S1_pi_0.20_vs_A1 | 47781 | 45 | 0 | 0 | 0 | 1264 | 0 |
| S1_pi_0.20_vs_A2 | 47781 | 63 | 90 | 0 | 4 | 1916 | 4719 |

## 4. Negative-Bearing Competitions & Null-Layer Summary (Pooled)

| Arm | Locus Trained | Neg-Bearing Comps | All-Neg Comps | Sup Neg Members | Mean Neg Rank | Mean Neg Rel Rank | Neg Ranked First | Mean P(null) on Neg Comps | Mean Naive Null on Neg Comps | All Groups (>=2) Marginal NLL |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A1 | false | 5 | 3 | 5 | 1.4000 | 1.0000 | 0.6000 | — | 0.7990 | — |
| A2 | false | 5 | 3 | 5 | 1.4000 | 1.0000 | 0.6000 | — | 0.7505 | — |
| A3 | false | 5 | 3 | 5 | 1.4000 | 1.0000 | 0.6000 | — | 0.4564 | — |
| A4 | false | 5 | 3 | 5 | 1.4000 | 1.0000 | 0.6000 | — | 0.0010 | — |
| A5 | false | 5 | 3 | 5 | 1.2000 | 0.9000 | 0.6000 | — | 0.0595 | — |
| A3_gate | true | 5 | 3 | 5 | 1.4000 | 1.0000 | 0.6000 | 0.5753 | 0.4564 | 0.1445 |
| S1_pi_0.20 | false | 5 | 3 | 5 | 1.4000 | 1.0000 | 0.6000 | — | 0.5388 | — |

## 5. Abstention Reasons by Arm (Pooled)

| Arm | Total Rows | Supported | Abstained | Reasons |
|---|---:|---:|---:|---|
| A1 | 69789 | 47781 | 22008 | profile_lacks_both_classes_in_train_and_calibration:15890, sequence_score_outside_training_calibration_range:381, unseen_source_function_role_or_availability_profile:5737 |
| A2 | 69789 | 47781 | 22008 | profile_lacks_both_classes_in_train_and_calibration:15890, sequence_score_outside_training_calibration_range:381, unseen_source_function_role_or_availability_profile:5737 |
| A3 | 69789 | 47781 | 22008 | profile_lacks_both_classes_in_train_and_calibration:15890, sequence_score_outside_training_calibration_range:381, unseen_source_function_role_or_availability_profile:5737 |
| A4 | 69789 | 47781 | 22008 | profile_lacks_both_classes_in_train_and_calibration:15890, sequence_score_outside_training_calibration_range:381, unseen_source_function_role_or_availability_profile:5737 |
| A5 | 69789 | 47781 | 22008 | profile_lacks_both_classes_in_train_and_calibration:15890, sequence_score_outside_training_calibration_range:381, unseen_source_function_role_or_availability_profile:5737 |
| A3_gate | 69789 | 47781 | 22008 | profile_lacks_both_classes_in_train_and_calibration:15890, sequence_score_outside_training_calibration_range:381, unseen_source_function_role_or_availability_profile:5737 |
| S1_pi_0.20 | 69789 | 47781 | 22008 | profile_lacks_both_classes_in_train_and_calibration:15890, sequence_score_outside_training_calibration_range:381, unseen_source_function_role_or_availability_profile:5737 |
