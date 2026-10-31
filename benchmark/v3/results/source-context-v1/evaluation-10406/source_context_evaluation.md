# Source x context evaluation: evaluation.json

Cells: metric (supported/total test rows). Precision/recall cells also give the number of supported rows predicted positive at the frozen display threshold. Metric scope as recorded: "supported rows only; input_coverage reports abstentions over the full test set".

## AUROC: logistic

| source | SEQONLY | UNIVERSAL | PATHWAY | ORGANISM |
|---|---|---|---|---|
| diamond_homology | 0.723 (463/608) | 0.761 (462/608) | 0.759 (461/608) | 0.766 (460/608) |
| diamond_family_excluded | 0.682 (300/416) | 0.751 (299/416) | 0.737 (297/416) | 0.752 (299/416) |
| clean | 0.768 (38/38) | 0.042 (36/38) | 0.022 (35/38) | 0.032 (30/38) |
| clean_diamond_homology | 0.778 (483/631) | 0.776 (459/631) | 0.768 (461/631) | 0.759 (456/631) |
| clean_diamond_family_excluded | 0.808 (333/449) | 0.809 (317/449) | 0.801 (312/449) | 0.797 (311/449) |

## AUROC: gbdt

| source | SEQONLY | UNIVERSAL | PATHWAY | ORGANISM |
|---|---|---|---|---|
| diamond_homology | 0.688 (463/608) | 0.719 (462/608) | 0.710 (461/608) | 0.715 (460/608) |
| diamond_family_excluded | 0.713 (300/416) | 0.745 (299/416) | 0.737 (297/416) | 0.756 (299/416) |
| clean | 0.525 (38/38) | 0.628 (36/38) | 0.536 (35/38) | 0.582 (30/38) |
| clean_diamond_homology | 0.776 (483/631) | 0.741 (459/631) | 0.723 (461/631) | 0.713 (456/631) |
| clean_diamond_family_excluded | 0.830 (333/449) | 0.814 (317/449) | 0.814 (312/449) | 0.818 (311/449) |

## Log-loss: logistic

| source | SEQONLY | UNIVERSAL | PATHWAY | ORGANISM |
|---|---|---|---|---|
| diamond_homology | 0.432 (463/608) | 0.402 (462/608) | 0.413 (461/608) | 0.400 (460/608) |
| diamond_family_excluded | 0.344 (300/416) | 0.309 (299/416) | 0.300 (297/416) | 0.307 (299/416) |
| clean | 0.755 (38/38) | 0.707 (36/38) | 0.722 (35/38) | 0.656 (30/38) |
| clean_diamond_homology | 0.385 (483/631) | 0.366 (459/631) | 0.366 (461/631) | 0.369 (456/631) |
| clean_diamond_family_excluded | 0.329 (333/449) | 0.305 (317/449) | 0.310 (312/449) | 0.309 (311/449) |

## Log-loss: gbdt

| source | SEQONLY | UNIVERSAL | PATHWAY | ORGANISM |
|---|---|---|---|---|
| diamond_homology | 0.462 (463/608) | 0.424 (462/608) | 0.432 (461/608) | 0.426 (460/608) |
| diamond_family_excluded | 0.329 (300/416) | 0.311 (299/416) | 0.302 (297/416) | 0.305 (299/416) |
| clean | 0.744 (38/38) | 0.727 (36/38) | 0.866 (35/38) | 0.725 (30/38) |
| clean_diamond_homology | 0.392 (483/631) | 0.416 (459/631) | 0.413 (461/631) | 0.413 (456/631) |
| clean_diamond_family_excluded | 0.330 (333/449) | 0.350 (317/449) | 0.339 (312/449) | 0.331 (311/449) |

## Precision / recall at the display threshold: logistic

| source | SEQONLY | UNIVERSAL | PATHWAY | ORGANISM |
|---|---|---|---|---|
| diamond_homology | 0.333 / 0.236 (51 pred.+; 463/608) | 0.424 / 0.194 (33 pred.+; 462/608) | 0.424 / 0.194 (33 pred.+; 461/608) | 0.500 / 0.113 (16 pred.+; 460/608) |
| diamond_family_excluded | 0.136 / 0.188 (44 pred.+; 300/416) | 0.188 / 0.188 (32 pred.+; 299/416) | 0.148 / 0.133 (27 pred.+; 297/416) | 0.200 / 0.156 (25 pred.+; 299/416) |
| clean | 0.692 / 0.750 (26 pred.+; 38/38) | 0.667 / 1.000 (36 pred.+; 36/38) | 0.636 / 0.913 (33 pred.+; 35/38) | 0.690 / 0.952 (29 pred.+; 30/38) |
| clean_diamond_homology | 0.677 / 0.266 (31 pred.+; 483/631) | 0.431 / 0.314 (51 pred.+; 459/631) | 1.000 / 0.147 (10 pred.+; 461/631) | 0.362 / 0.262 (47 pred.+; 456/631) |
| clean_diamond_family_excluded | 0.700 / 0.412 (30 pred.+; 333/449) | 0.944 / 0.370 (18 pred.+; 317/449) | 0.941 / 0.364 (17 pred.+; 312/449) | 0.682 / 0.349 (22 pred.+; 311/449) |

## Precision / recall at the display threshold: gbdt

| source | SEQONLY | UNIVERSAL | PATHWAY | ORGANISM |
|---|---|---|---|---|
| diamond_homology | 0.303 / 0.139 (33 pred.+; 463/608) | 0.400 / 0.278 (50 pred.+; 462/608) | 0.343 / 0.319 (67 pred.+; 461/608) | 0.370 / 0.282 (54 pred.+; 460/608) |
| diamond_family_excluded | 0.143 / 0.031 (7 pred.+; 300/416) | 0.279 / 0.375 (43 pred.+; 299/416) | 0.267 / 0.267 (30 pred.+; 297/416) | 0.275 / 0.344 (40 pred.+; 299/416) |
| clean | 0.708 / 0.708 (24 pred.+; 38/38) | 0.731 / 0.792 (26 pred.+; 36/38) | 0.720 / 0.783 (25 pred.+; 35/38) | 0.750 / 0.857 (24 pred.+; 30/38) |
| clean_diamond_homology | 0.690 / 0.253 (29 pred.+; 483/631) | 0.529 / 0.129 (17 pred.+; 459/631) | 0.556 / 0.147 (18 pred.+; 461/631) | 0.467 / 0.108 (15 pred.+; 456/631) |
| clean_diamond_family_excluded | 0.677 / 0.412 (31 pred.+; 333/449) | 0.586 / 0.370 (29 pred.+; 317/449) | 0.640 / 0.364 (25 pred.+; 312/449) | 0.625 / 0.349 (24 pred.+; 311/449) |

## Unsupported test rows by reason

| source | context | kind | total | supported | `profile_lacks_both_classes_in_train_and_calibration` | `sequence_score_outside_training_calibration_range` |
|---|---|---|---|---|---|---|
| diamond_homology | SEQONLY | logistic | 608 | 463 | 145 | 0 |
| diamond_homology | SEQONLY | gbdt | 608 | 463 | 145 | 0 |
| diamond_homology | UNIVERSAL | logistic | 608 | 462 | 146 | 0 |
| diamond_homology | UNIVERSAL | gbdt | 608 | 462 | 146 | 0 |
| diamond_homology | PATHWAY | logistic | 608 | 461 | 145 | 2 |
| diamond_homology | PATHWAY | gbdt | 608 | 461 | 145 | 2 |
| diamond_homology | ORGANISM | logistic | 608 | 460 | 147 | 1 |
| diamond_homology | ORGANISM | gbdt | 608 | 460 | 147 | 1 |
| diamond_family_excluded | SEQONLY | logistic | 416 | 300 | 114 | 2 |
| diamond_family_excluded | SEQONLY | gbdt | 416 | 300 | 114 | 2 |
| diamond_family_excluded | UNIVERSAL | logistic | 416 | 299 | 115 | 2 |
| diamond_family_excluded | UNIVERSAL | gbdt | 416 | 299 | 115 | 2 |
| diamond_family_excluded | PATHWAY | logistic | 416 | 297 | 117 | 2 |
| diamond_family_excluded | PATHWAY | gbdt | 416 | 297 | 117 | 2 |
| diamond_family_excluded | ORGANISM | logistic | 416 | 299 | 115 | 2 |
| diamond_family_excluded | ORGANISM | gbdt | 416 | 299 | 115 | 2 |
| clean | SEQONLY | logistic | 38 | 38 | 0 | 0 |
| clean | SEQONLY | gbdt | 38 | 38 | 0 | 0 |
| clean | UNIVERSAL | logistic | 38 | 36 | 0 | 2 |
| clean | UNIVERSAL | gbdt | 38 | 36 | 0 | 2 |
| clean | PATHWAY | logistic | 38 | 35 | 0 | 3 |
| clean | PATHWAY | gbdt | 38 | 35 | 0 | 3 |
| clean | ORGANISM | logistic | 38 | 30 | 5 | 3 |
| clean | ORGANISM | gbdt | 38 | 30 | 5 | 3 |
| clean_diamond_homology | SEQONLY | logistic | 631 | 483 | 145 | 3 |
| clean_diamond_homology | SEQONLY | gbdt | 631 | 483 | 145 | 3 |
| clean_diamond_homology | UNIVERSAL | logistic | 631 | 459 | 167 | 5 |
| clean_diamond_homology | UNIVERSAL | gbdt | 631 | 459 | 167 | 5 |
| clean_diamond_homology | PATHWAY | logistic | 631 | 461 | 165 | 5 |
| clean_diamond_homology | PATHWAY | gbdt | 631 | 461 | 165 | 5 |
| clean_diamond_homology | ORGANISM | logistic | 631 | 456 | 169 | 6 |
| clean_diamond_homology | ORGANISM | gbdt | 631 | 456 | 169 | 6 |
| clean_diamond_family_excluded | SEQONLY | logistic | 449 | 333 | 114 | 2 |
| clean_diamond_family_excluded | SEQONLY | gbdt | 449 | 333 | 114 | 2 |
| clean_diamond_family_excluded | UNIVERSAL | logistic | 449 | 317 | 127 | 5 |
| clean_diamond_family_excluded | UNIVERSAL | gbdt | 449 | 317 | 127 | 5 |
| clean_diamond_family_excluded | PATHWAY | logistic | 449 | 312 | 134 | 3 |
| clean_diamond_family_excluded | PATHWAY | gbdt | 449 | 312 | 134 | 3 |
| clean_diamond_family_excluded | ORGANISM | logistic | 449 | 311 | 134 | 4 |
| clean_diamond_family_excluded | ORGANISM | gbdt | 449 | 311 | 134 | 4 |

## Limitations

Quoted from the evaluation JSON.

- scope: "Conditional candidate-table metrics on a reused reference-overlapping cohort; not the primary endpoint."
- limitations: "Previously examined cohort that overlaps the CLEAN reference (exact known sequences, not unseen families); development reuse limits apply. Candidate-table metrics are conditional on the exported candidates; no biological accuracy, calibration or promotion claim. SEQONLY retains the declared organism taxon and essential-profile metadata and function definitions, but disables priors and context stages."
- pending_stages: `full_reference_ec_endpoint`: not_run, `production_replay`: not_run
