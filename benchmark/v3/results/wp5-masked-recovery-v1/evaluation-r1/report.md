# WP5 masked recovery: evaluation tables

Function recovery, phenotype agreement and model-conditional capacity are separate. Simulation agreement is not biological confirmation.

## Function recovery at fixed precision (test split)

Denominators include every eligible test case, including candidate-generation failures. Random seeds are averaged per case.

| Arm | Scope | Policy | Threshold | N | Recovery | 95% CI | Harm | Abstention | False null | Correct negative | Candidate recall | Conditional recovery |
|---|---|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|
| evidence:k3 | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.769 | 0.000 | 0.000 | 0.692 | 0.000 |
| evidence:k3 | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.909 | 0.000 | 0.000 | 0.818 | 0.000 |
| evidence:k6 | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.769 | 0.000 | 0.000 | 0.692 | 0.000 |
| evidence:k6 | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.909 | 0.000 | 0.000 | 0.818 | 0.000 |
| exhaustive | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.769 | 0.000 | 0.000 | 0.692 | 0.000 |
| exhaustive | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.909 | 0.000 | 0.000 | 0.818 | 0.000 |
| fixed_panel:k3 | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.769 | 0.000 | 0.000 | 0.692 | 0.000 |
| fixed_panel:k3 | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.909 | 0.000 | 0.000 | 0.818 | 0.000 |
| fixed_panel:k6 | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.769 | 0.000 | 0.000 | 0.692 | 0.000 |
| fixed_panel:k6 | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.909 | 0.000 | 0.000 | 0.818 | 0.000 |
| minimal_medium:k3 | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.769 | 0.000 | 0.000 | 0.692 | 0.000 |
| minimal_medium:k3 | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.909 | 0.000 | 0.000 | 0.818 | 0.000 |
| minimal_medium:k6 | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.769 | 0.000 | 0.000 | 0.692 | 0.000 |
| minimal_medium:k6 | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.909 | 0.000 | 0.000 | 0.818 | 0.000 |
| qwen:k3 | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.538 | 0.000 | 0.000 | 0.692 | 0.000 |
| qwen:k3 | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.455 | 0.000 | 0.000 | 0.818 | 0.000 |
| qwen:k6 | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.615 | 0.000 | 0.000 | 0.692 | 0.000 |
| qwen:k6 | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.000 | 0.000 | 0.000 | 0.818 | 0.000 |
| random:k3 | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.769 | 0.000 | 0.000 | 0.692 | 0.000 |
| random:k3 | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.909 | 0.000 | 0.000 | 0.818 | 0.000 |
| random:k6 | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.769 | 0.000 | 0.000 | 0.692 | 0.000 |
| random:k6 | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.900 | 0.009 | 0.000 | 0.818 | 0.000 |
| sequence_only | primary | selected_on_calibration | 0.864 | 13 | 0.231 | 0.082–0.503 | 0.000 | 0.538 | 0.000 | 0.000 | 0.692 | 0.333 |
| sequence_only | sensitivity | selected_on_calibration | 0.864 | 11 | 0.273 | 0.097–0.566 | 0.000 | 0.636 | 0.000 | 0.000 | 0.818 | 0.333 |

## Descriptive test precision–recovery curves (primary scope)

These curves do not select or replace the calibration-frozen operating point.

| Arm | Threshold | Decisions | Precision | Recovery | Harm |
|---|---:|---:|---:|---:|---:|
| evidence:k3 | 0.430 | 3.000 | 1.000 | 0.231 | 0.000 |
| evidence:k3 | 0.430 | 2.000 | 1.000 | 0.154 | 0.000 |
| evidence:k3 | 0.762 | 1.000 | 1.000 | 0.077 | 0.000 |
| evidence:k6 | 0.430 | 3.000 | 1.000 | 0.231 | 0.000 |
| evidence:k6 | 0.430 | 2.000 | 1.000 | 0.154 | 0.000 |
| evidence:k6 | 0.762 | 1.000 | 1.000 | 0.077 | 0.000 |
| exhaustive | 0.762 | 1.000 | 1.000 | 0.077 | 0.000 |
| fixed_panel:k3 | 0.430 | 3.000 | 1.000 | 0.231 | 0.000 |
| fixed_panel:k3 | 0.430 | 2.000 | 1.000 | 0.154 | 0.000 |
| fixed_panel:k3 | 0.762 | 1.000 | 1.000 | 0.077 | 0.000 |
| fixed_panel:k6 | 0.430 | 3.000 | 1.000 | 0.231 | 0.000 |
| fixed_panel:k6 | 0.430 | 2.000 | 1.000 | 0.154 | 0.000 |
| fixed_panel:k6 | 0.762 | 1.000 | 1.000 | 0.077 | 0.000 |
| minimal_medium:k3 | 0.430 | 3.000 | 1.000 | 0.231 | 0.000 |
| minimal_medium:k3 | 0.430 | 2.000 | 1.000 | 0.154 | 0.000 |
| minimal_medium:k3 | 0.762 | 1.000 | 1.000 | 0.077 | 0.000 |
| minimal_medium:k6 | 0.430 | 3.000 | 1.000 | 0.231 | 0.000 |
| minimal_medium:k6 | 0.430 | 2.000 | 1.000 | 0.154 | 0.000 |
| minimal_medium:k6 | 0.762 | 1.000 | 1.000 | 0.077 | 0.000 |
| qwen:k6 | 0.430 | 1.000 | 1.000 | 0.077 | 0.000 |
| sequence_only | 0.430 | 7.000 | 1.000 | 0.538 | 0.000 |
| sequence_only | 0.430 | 6.000 | 1.000 | 0.462 | 0.000 |
| sequence_only | 0.594 | 5.000 | 1.000 | 0.385 | 0.000 |
| sequence_only | 0.762 | 4.000 | 1.000 | 0.308 | 0.000 |
| sequence_only | 0.937 | 3.000 | 1.000 | 0.231 | 0.000 |
| sequence_only | 0.951 | 2.000 | 1.000 | 0.154 | 0.000 |
| sequence_only | 0.960 | 1.000 | 1.000 | 0.077 | 0.000 |

## Paired comparisons (primary scope)

| Arm | Comparator | McNemar only-arm | only-comparator | p | Bootstrap 95% recovery difference |
|---|---|---:|---:|---:|---|
| evidence:k3 | sequence_only | 0 | 3 | 0.250 | -0.462 to 0.000 |
| evidence:k3 | minimal_medium:k3 | 0 | 0 | 1.000 | 0.000 to 0.000 |
| evidence:k6 | sequence_only | 0 | 3 | 0.250 | -0.462 to 0.000 |
| evidence:k6 | minimal_medium:k6 | 0 | 0 | 1.000 | 0.000 to 0.000 |
| exhaustive | sequence_only | 0 | 3 | 0.250 | -0.462 to 0.000 |
| fixed_panel:k3 | sequence_only | 0 | 3 | 0.250 | -0.462 to 0.000 |
| fixed_panel:k3 | minimal_medium:k3 | 0 | 0 | 1.000 | 0.000 to 0.000 |
| fixed_panel:k6 | sequence_only | 0 | 3 | 0.250 | -0.462 to 0.000 |
| fixed_panel:k6 | minimal_medium:k6 | 0 | 0 | 1.000 | 0.000 to 0.000 |
| minimal_medium:k3 | sequence_only | 0 | 3 | 0.250 | -0.462 to 0.000 |
| minimal_medium:k6 | sequence_only | 0 | 3 | 0.250 | -0.462 to 0.000 |
| qwen:k3 | sequence_only | 0 | 3 | 0.250 | -0.462 to 0.000 |
| qwen:k3 | minimal_medium:k3 | 0 | 0 | 1.000 | 0.000 to 0.000 |
| qwen:k3 | fixed_panel:k3 | 0 | 0 | 1.000 | 0.000 to 0.000 |
| qwen:k3 | random:k3 | – | – | – | 0.000 to 0.000 |
| qwen:k3 | evidence:k3 | 0 | 0 | 1.000 | 0.000 to 0.000 |
| qwen:k6 | sequence_only | 0 | 3 | 0.250 | -0.462 to 0.000 |
| qwen:k6 | minimal_medium:k6 | 0 | 0 | 1.000 | 0.000 to 0.000 |
| qwen:k6 | fixed_panel:k6 | 0 | 0 | 1.000 | 0.000 to 0.000 |
| qwen:k6 | random:k6 | – | – | – | 0.000 to 0.000 |
| qwen:k6 | evidence:k6 | 0 | 0 | 1.000 | 0.000 to 0.000 |
| random:k3 | sequence_only | – | – | – | -0.462 to 0.000 |
| random:k3 | minimal_medium:k3 | – | – | – | 0.000 to 0.000 |
| random:k6 | sequence_only | – | – | – | -0.462 to 0.000 |
| random:k6 | minimal_medium:k6 | – | – | – | 0.000 to 0.000 |

## Phenotype agreement

Called-agreement uses binary-evaluable pairs with a called prediction; coverage shows the called-pair fraction of all binary pairs. Partial growth is excluded and reported separately.

| Arm | Scope | Set | Binary pairs | Called pairs | Called agree | Called agreement | Coverage | Constant growth | Constant no-growth | Partial pairs |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| evidence:k3 | primary | visited | 22 | 0 | 0 | – | 0.000 | 0.182 | 0.818 | 0 |
| evidence:k3 | primary | common exhaustive | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| evidence:k3 | sensitivity | visited | 22 | 0 | 0 | – | 0.000 | 0.182 | 0.818 | 0 |
| evidence:k3 | sensitivity | common exhaustive | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| evidence:k6 | primary | visited | 28 | 0 | 0 | – | 0.000 | 0.357 | 0.643 | 0 |
| evidence:k6 | primary | common exhaustive | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| evidence:k6 | sensitivity | visited | 28 | 0 | 0 | – | 0.000 | 0.357 | 0.643 | 0 |
| evidence:k6 | sensitivity | common exhaustive | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| exhaustive | primary | visited | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| exhaustive | primary | common exhaustive | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| exhaustive | sensitivity | visited | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| exhaustive | sensitivity | common exhaustive | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| fixed_panel:k3 | primary | visited | 21 | 0 | 0 | – | 0.000 | 0.143 | 0.857 | 0 |
| fixed_panel:k3 | primary | common exhaustive | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| fixed_panel:k3 | sensitivity | visited | 21 | 0 | 0 | – | 0.000 | 0.143 | 0.857 | 0 |
| fixed_panel:k3 | sensitivity | common exhaustive | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| fixed_panel:k6 | primary | visited | 27 | 0 | 0 | – | 0.000 | 0.333 | 0.667 | 0 |
| fixed_panel:k6 | primary | common exhaustive | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| fixed_panel:k6 | sensitivity | visited | 27 | 0 | 0 | – | 0.000 | 0.333 | 0.667 | 0 |
| fixed_panel:k6 | sensitivity | common exhaustive | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| minimal_medium:k3 | primary | visited | 10 | 0 | 0 | – | 0.000 | 0.100 | 0.900 | 0 |
| minimal_medium:k3 | primary | common exhaustive | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| minimal_medium:k3 | sensitivity | visited | 10 | 0 | 0 | – | 0.000 | 0.100 | 0.900 | 0 |
| minimal_medium:k3 | sensitivity | common exhaustive | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| minimal_medium:k6 | primary | visited | 10 | 0 | 0 | – | 0.000 | 0.100 | 0.900 | 0 |
| minimal_medium:k6 | primary | common exhaustive | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| minimal_medium:k6 | sensitivity | visited | 10 | 0 | 0 | – | 0.000 | 0.100 | 0.900 | 0 |
| minimal_medium:k6 | sensitivity | common exhaustive | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| qwen:k3 | primary | visited | 13 | 0 | 0 | – | 0.000 | 0.077 | 0.923 | 0 |
| qwen:k3 | primary | common exhaustive | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| qwen:k3 | sensitivity | visited | 8 | 0 | 0 | – | 0.000 | 0.125 | 0.875 | 0 |
| qwen:k3 | sensitivity | common exhaustive | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| qwen:k6 | primary | visited | 31 | 0 | 0 | – | 0.000 | 0.161 | 0.839 | 0 |
| qwen:k6 | primary | common exhaustive | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| qwen:k6 | sensitivity | visited | 0 | 0 | 0 | – | – | – | – | 0 |
| qwen:k6 | sensitivity | common exhaustive | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| random:k3 | primary | visited | 79 | 0 | 0 | – | 0.000 | 0.582 | 0.418 | 0 |
| random:k3 | primary | common exhaustive | 7600 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| random:k3 | sensitivity | visited | 79 | 0 | 0 | – | 0.000 | 0.582 | 0.418 | 0 |
| random:k3 | sensitivity | common exhaustive | 7600 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| random:k6 | primary | visited | 150 | 0 | 0 | – | 0.000 | 0.600 | 0.400 | 0 |
| random:k6 | primary | common exhaustive | 7600 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| random:k6 | sensitivity | visited | 150 | 2 | 2 | 1.000 | 0.013 | 0.600 | 0.400 | 0 |
| random:k6 | sensitivity | common exhaustive | 7600 | 80 | 44 | 0.550 | 0.011 | 0.533 | 0.467 | 0 |
| sequence_only | primary | visited | 0 | 0 | 0 | – | – | – | – | 0 |
| sequence_only | primary | common exhaustive | 760 | 0 | 0 | – | 0.000 | 0.533 | 0.467 | 0 |
| sequence_only | sensitivity | visited | 0 | 0 | 0 | – | – | – | – | 0 |
| sequence_only | sensitivity | common exhaustive | 760 | 240 | 235 | 0.979 | 0.316 | 0.533 | 0.467 | 0 |

The unmasked-original comparator is unavailable: the public-only §8 runner has no such stage. No private original model is silently introduced into evaluation.

## Model-conditional capacity and cost

| Arm | Scope | C1 cases | C1 discriminating | C1 fraction | LPs used | LP cap | Unused LPs | Runs | Failed runs | Qwen requests | Qwen failures |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| evidence:k3 | primary | 10 | 4.000 | 0.400 | 182 | 270 | 88 | 14 | 0 | 0 | 0 |
| evidence:k3 | sensitivity | 10 | 8.000 | 0.800 | 270 | 270 | 0 | 12 | 0 | 0 | 0 |
| evidence:k6 | primary | 10 | 4.000 | 0.400 | 364 | 540 | 176 | 14 | 0 | 0 | 0 |
| evidence:k6 | sensitivity | 10 | 8.000 | 0.800 | 540 | 540 | 0 | 12 | 0 | 0 | 0 |
| exhaustive | primary | 10 | 4.000 | 0.400 | 17766 | 26656 | 8890 | 14 | 0 | 0 | 0 |
| exhaustive | sensitivity | 10 | 8.000 | 0.800 | 25996 | 26656 | 660 | 12 | 0 | 0 | 0 |
| fixed_panel:k3 | primary | 10 | 4.000 | 0.400 | 196 | 270 | 74 | 14 | 0 | 0 | 0 |
| fixed_panel:k3 | sensitivity | 10 | 8.000 | 0.800 | 270 | 270 | 0 | 12 | 0 | 0 | 0 |
| fixed_panel:k6 | primary | 10 | 4.000 | 0.400 | 378 | 540 | 162 | 14 | 0 | 0 | 0 |
| fixed_panel:k6 | sensitivity | 10 | 8.000 | 0.800 | 540 | 540 | 0 | 12 | 0 | 0 | 0 |
| minimal_medium:k3 | primary | 10 | 4.000 | 0.400 | 70 | 270 | 200 | 14 | 0 | 0 | 0 |
| minimal_medium:k3 | sensitivity | 10 | 8.000 | 0.800 | 90 | 270 | 180 | 12 | 0 | 0 | 0 |
| minimal_medium:k6 | primary | 10 | 4.000 | 0.400 | 70 | 540 | 470 | 14 | 0 | 0 | 0 |
| minimal_medium:k6 | sensitivity | 10 | 8.000 | 0.800 | 90 | 540 | 450 | 12 | 0 | 0 | 0 |
| qwen:k3 | primary | 10 | 0.000 | 0.000 | 148 | 270 | 122 | 14 | 3 | 31 | 1 |
| qwen:k3 | sensitivity | 10 | 1.000 | 0.100 | 154 | 270 | 116 | 12 | 6 | 31 | 1 |
| qwen:k6 | primary | 10 | 2.000 | 0.200 | 310 | 540 | 230 | 14 | 2 | 60 | 0 |
| qwen:k6 | sensitivity | 10 | 0.000 | 0.000 | 132 | 540 | 408 | 12 | 11 | 31 | 3 |
| random:k3 | primary | 10 | 1.900 | 0.190 | 1806 | 2700 | 894 | 140 | 0 | 0 | 0 |
| random:k3 | sensitivity | 10 | 8.000 | 0.800 | 2666 | 2700 | 34 | 120 | 0 | 0 | 0 |
| random:k6 | primary | 10 | 2.300 | 0.230 | 3584 | 5400 | 1816 | 140 | 0 | 0 | 0 |
| random:k6 | sensitivity | 10 | 8.000 | 0.800 | 5292 | 5400 | 108 | 120 | 0 | 0 | 0 |
| sequence_only | primary | 10 | 0.000 | 0.000 | 0 | 0 | 0 | 14 | 0 | 0 | 0 |
| sequence_only | sensitivity | 10 | 0.000 | 0.000 | 0 | 0 | 0 | 12 | 0 | 0 | 0 |

## Per-stratum and per-species tables (primary scope)

| Arm | Group | N | Recovery | Harm | False null | Candidate recall | Conditional recovery |
|---|---|---:|---:|---:|---:|---:|---:|
| evidence:k3 | species:Bs | 6 | 0.000 | 0.000 | 0.000 | 0.833 | 0.000 |
| evidence:k3 | species:Ec | 7 | 0.000 | 0.000 | 0.000 | 0.571 | 0.000 |
| evidence:k3 | stratum:catalytic | 11 | 0.000 | 0.000 | 0.000 | 0.818 | 0.000 |
| evidence:k3 | stratum:transport_only | 2 | 0.000 | 0.000 | 0.000 | 0.000 | – |
| evidence:k6 | species:Bs | 6 | 0.000 | 0.000 | 0.000 | 0.833 | 0.000 |
| evidence:k6 | species:Ec | 7 | 0.000 | 0.000 | 0.000 | 0.571 | 0.000 |
| evidence:k6 | stratum:catalytic | 11 | 0.000 | 0.000 | 0.000 | 0.818 | 0.000 |
| evidence:k6 | stratum:transport_only | 2 | 0.000 | 0.000 | 0.000 | 0.000 | – |
| exhaustive | species:Bs | 6 | 0.000 | 0.000 | 0.000 | 0.833 | 0.000 |
| exhaustive | species:Ec | 7 | 0.000 | 0.000 | 0.000 | 0.571 | 0.000 |
| exhaustive | stratum:catalytic | 11 | 0.000 | 0.000 | 0.000 | 0.818 | 0.000 |
| exhaustive | stratum:transport_only | 2 | 0.000 | 0.000 | 0.000 | 0.000 | – |
| fixed_panel:k3 | species:Bs | 6 | 0.000 | 0.000 | 0.000 | 0.833 | 0.000 |
| fixed_panel:k3 | species:Ec | 7 | 0.000 | 0.000 | 0.000 | 0.571 | 0.000 |
| fixed_panel:k3 | stratum:catalytic | 11 | 0.000 | 0.000 | 0.000 | 0.818 | 0.000 |
| fixed_panel:k3 | stratum:transport_only | 2 | 0.000 | 0.000 | 0.000 | 0.000 | – |
| fixed_panel:k6 | species:Bs | 6 | 0.000 | 0.000 | 0.000 | 0.833 | 0.000 |
| fixed_panel:k6 | species:Ec | 7 | 0.000 | 0.000 | 0.000 | 0.571 | 0.000 |
| fixed_panel:k6 | stratum:catalytic | 11 | 0.000 | 0.000 | 0.000 | 0.818 | 0.000 |
| fixed_panel:k6 | stratum:transport_only | 2 | 0.000 | 0.000 | 0.000 | 0.000 | – |
| minimal_medium:k3 | species:Bs | 6 | 0.000 | 0.000 | 0.000 | 0.833 | 0.000 |
| minimal_medium:k3 | species:Ec | 7 | 0.000 | 0.000 | 0.000 | 0.571 | 0.000 |
| minimal_medium:k3 | stratum:catalytic | 11 | 0.000 | 0.000 | 0.000 | 0.818 | 0.000 |
| minimal_medium:k3 | stratum:transport_only | 2 | 0.000 | 0.000 | 0.000 | 0.000 | – |
| minimal_medium:k6 | species:Bs | 6 | 0.000 | 0.000 | 0.000 | 0.833 | 0.000 |
| minimal_medium:k6 | species:Ec | 7 | 0.000 | 0.000 | 0.000 | 0.571 | 0.000 |
| minimal_medium:k6 | stratum:catalytic | 11 | 0.000 | 0.000 | 0.000 | 0.818 | 0.000 |
| minimal_medium:k6 | stratum:transport_only | 2 | 0.000 | 0.000 | 0.000 | 0.000 | – |
| qwen:k3 | species:Bs | 6 | 0.000 | 0.000 | 0.000 | 0.833 | 0.000 |
| qwen:k3 | species:Ec | 7 | 0.000 | 0.000 | 0.000 | 0.571 | 0.000 |
| qwen:k3 | stratum:catalytic | 11 | 0.000 | 0.000 | 0.000 | 0.818 | 0.000 |
| qwen:k3 | stratum:transport_only | 2 | 0.000 | 0.000 | 0.000 | 0.000 | – |
| qwen:k6 | species:Bs | 6 | 0.000 | 0.000 | 0.000 | 0.833 | 0.000 |
| qwen:k6 | species:Ec | 7 | 0.000 | 0.000 | 0.000 | 0.571 | 0.000 |
| qwen:k6 | stratum:catalytic | 11 | 0.000 | 0.000 | 0.000 | 0.818 | 0.000 |
| qwen:k6 | stratum:transport_only | 2 | 0.000 | 0.000 | 0.000 | 0.000 | – |
| random:k3 | species:Bs | 6 | 0.000 | 0.000 | 0.000 | 0.833 | 0.000 |
| random:k3 | species:Ec | 7 | 0.000 | 0.000 | 0.000 | 0.571 | 0.000 |
| random:k3 | stratum:catalytic | 11 | 0.000 | 0.000 | 0.000 | 0.818 | 0.000 |
| random:k3 | stratum:transport_only | 2 | 0.000 | 0.000 | 0.000 | 0.000 | – |
| random:k6 | species:Bs | 6 | 0.000 | 0.000 | 0.000 | 0.833 | 0.000 |
| random:k6 | species:Ec | 7 | 0.000 | 0.000 | 0.000 | 0.571 | 0.000 |
| random:k6 | stratum:catalytic | 11 | 0.000 | 0.000 | 0.000 | 0.818 | 0.000 |
| random:k6 | stratum:transport_only | 2 | 0.000 | 0.000 | 0.000 | 0.000 | – |
| sequence_only | species:Bs | 6 | 0.000 | 0.000 | 0.000 | 0.833 | 0.000 |
| sequence_only | species:Ec | 7 | 0.429 | 0.000 | 0.000 | 0.571 | 0.750 |
| sequence_only | stratum:catalytic | 11 | 0.273 | 0.000 | 0.000 | 0.818 | 0.333 |
| sequence_only | stratum:transport_only | 2 | 0.000 | 0.000 | 0.000 | 0.000 | – |

## Frozen protocol and limitations

- Keio (calibration) measured pairs were chosen by an iML1515 design that used the target GPR; the menu offered the full substrate grid plus glucose.
- Nutrient assays use oxygen uptake 20 mmol gDW^-1 h^-1; printed is primary and additional_inorganic is sensitivity.
- For supplemented nutrient conditions, the wild type is the corresponding unsupplemented wild type.
- Sulfonate sulfur uptake ceiling is 1.0; tauA and ssuA share a target-free background menu.
- No oxygen sweep, product objective, reaction-level condition knockout, learned predictor or Tier B case is admitted.
- Small denominators: every proportion carries its Wilson interval and denominator.
- Random arm: per-case mean over seeds; it enters bootstrap differences, not exact McNemar.
- Candidate-generation failures and transport controls remain in the all-eligible denominator.
- The unmasked-original phenotype comparator has no §8 runner stage and is marked unavailable, never invented.

## Calibration policy

> Keio (calibration) measured pairs were chosen by an iML1515 design that used the target GPR; the menu offered the full substrate grid plus glucose.

| Arm | Status | Threshold | Calibration counts |
|---|---|---:|---|
| evidence:k3 | no_supported_threshold_abstain | – | `{}` |
| evidence:k6 | no_supported_threshold_abstain | – | `{}` |
| exhaustive | no_supported_threshold_abstain | – | `{}` |
| fixed_panel:k3 | no_supported_threshold_abstain | – | `{}` |
| fixed_panel:k6 | no_supported_threshold_abstain | – | `{}` |
| minimal_medium:k3 | no_supported_threshold_abstain | – | `{}` |
| minimal_medium:k6 | no_supported_threshold_abstain | – | `{}` |
| qwen:k3 | no_supported_threshold_abstain | – | `{}` |
| qwen:k6 | no_supported_threshold_abstain | – | `{}` |
| random:k3 | no_supported_threshold_abstain | – | `{}` |
| random:k6 | no_supported_threshold_abstain | – | `{}` |
| sequence_only | selected_on_calibration | 0.864 | `{"abstain_policy": 6.0, "abstain_tied": 2.0, "correct_call": 9.0, "no_admissible_experiment": 23.0, "wrong_call": 1.0}` |
