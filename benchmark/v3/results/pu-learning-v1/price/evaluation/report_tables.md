# Price-148 PU Evaluation Tables (Row 4, §4.2, P8)

- Convention: unconfirmed calls are not negatives; support share counts them as unsupported (reference support)
- Policy sha256: `cddf53f637e20c33c23b434e0bbc8c43ff4c0d221e725526ca523aab05be5d77`
- Split sha256: `950209995fa673648f499121f92a48a9cee9288b9026a1f1ff808c7c504d5c55`

## 1. Arms × Test Metrics (Calibrated WP5 Policy, Fixed 0.5, and Artifact Display Threshold)

| Analysis | Arm | Status | Proteins | Cal Threshold | Cal Calls | Cal Supported | Cal Unconfirmed | Support Share [95% CI] | Ref Recall [95% CI] | Top-1 Match | Top-1 Share [95% CI] | Top-3 Match | Display Thr | Display Calls (Sup / Unconf) | Pos-Pair Brier | Top-Call ECE | LS Top-Call ECE |
|---|---|---|---:|---:|---:|---:|---:|---|---|---:|---|---:|---:|---|---:|---:|---:|
| primary/all | raw | no_supported_threshold_abstain | 70 | — | — | — | — | — | — | 22 | 0.271 [0.198, 0.388] | 28 | — | — | 0.7608 | 0.0954 | — |
| primary/all | stage0 | no_supported_threshold_abstain | 70 | — | — | — | — | — | — | 22 | 0.271 [0.198, 0.388] | 28 | — | — | 0.8134 | 0.0612 | — |
| primary/all | stage1 | no_supported_threshold_abstain | 70 | — | — | — | — | — | — | 22 | 0.271 [0.198, 0.388] | 28 | — | — | 0.8134 | 0.0612 | — |
| primary/all | learned_pn | no_supported_threshold_abstain | 70 | — | — | — | — | — | — | 22 | 0.271 [0.198, 0.388] | 28 | 0.3212 | 30 (11 / 19) | 0.8059 | 0.0671 | 0.0657 |
| primary/all | learned_pn_gbdt | no_supported_threshold_abstain | 70 | — | — | — | — | — | — | 21 | 0.251 [0.176, 0.367] | 26 | 0.2232 | 31 (10 / 21) | 0.8787 | 0.1156 | 0.1182 |
| primary/all | learned_pu | no_supported_threshold_abstain | 70 | — | — | — | — | — | — | 22 | 0.271 [0.198, 0.388] | 28 | 0.9566 | 18 (7 / 11) | 0.6170 | 0.4562 | 0.4562 |
| primary/all | learned_pu_gbdt | no_supported_threshold_abstain | 70 | — | — | — | — | — | — | 21 | 0.238 [0.164, 0.357] | 28 | 0.9951 | 26 (10 / 16) | 0.6368 | 0.4613 | 0.4613 |
| primary/all | learned_pu_rank | no_supported_threshold_abstain | 70 | — | — | — | — | — | — | 22 | 0.271 [0.198, 0.388] | 28 | 0.9988 | 19 (9 / 10) | 0.6027 | 0.5988 | 0.5988 |
| swapped/all | raw | no_supported_threshold_abstain | 78 | — | — | — | — | — | — | 27 | 0.333 [0.230, 0.430] | 37 | — | — | 0.6840 | 0.1938 | — |
| swapped/all | stage0 | no_supported_threshold_abstain | 78 | — | — | — | — | — | — | 27 | 0.333 [0.230, 0.430] | 37 | — | — | 0.7582 | 0.1768 | — |
| swapped/all | stage1 | no_supported_threshold_abstain | 78 | — | — | — | — | — | — | 27 | 0.333 [0.230, 0.430] | 37 | — | — | 0.7582 | 0.1768 | — |
| swapped/all | learned_pn | no_supported_threshold_abstain | 78 | — | — | — | — | — | — | 25 | 0.321 [0.209, 0.424] | 37 | 0.3212 | 38 (13 / 25) | 0.7653 | 0.1805 | 0.1805 |
| swapped/all | learned_pn_gbdt | no_supported_threshold_abstain | 78 | — | — | — | — | — | — | 26 | 0.325 [0.231, 0.414] | 32 | 0.2232 | 39 (12 / 27) | 0.8439 | 0.1854 | 0.1854 |
| swapped/all | learned_pu | no_supported_threshold_abstain | 78 | — | — | — | — | — | — | 25 | 0.321 [0.209, 0.424] | 37 | 0.9566 | 18 (4 / 14) | 0.5224 | 0.4261 | 0.4261 |
| swapped/all | learned_pu_gbdt | no_supported_threshold_abstain | 78 | — | — | — | — | — | — | 15 | 0.208 [0.128, 0.287] | 37 | 0.9951 | 19 (6 / 13) | 0.5235 | 0.5522 | 0.5522 |
| swapped/all | learned_pu_rank | no_supported_threshold_abstain | 78 | — | — | — | — | — | — | 25 | 0.321 [0.209, 0.424] | 37 | 0.9988 | 26 (8 / 18) | 0.5128 | 0.5118 | 0.5118 |
| primary/without_GCA_000006965.1 | raw | no_supported_threshold_abstain | 52 | — | — | — | — | — | — | 18 | 0.292 [0.189, 0.426] | 23 | — | — | 0.7420 | 0.1319 | — |
| primary/without_GCA_000006965.1 | stage0 | no_supported_threshold_abstain | 52 | — | — | — | — | — | — | 18 | 0.292 [0.189, 0.426] | 23 | — | — | 0.7994 | 0.0793 | — |
| primary/without_GCA_000006965.1 | stage1 | no_supported_threshold_abstain | 52 | — | — | — | — | — | — | 18 | 0.292 [0.189, 0.426] | 23 | — | — | 0.7994 | 0.0793 | — |
| primary/without_GCA_000006965.1 | learned_pn | no_supported_threshold_abstain | 52 | — | — | — | — | — | — | 18 | 0.292 [0.189, 0.426] | 23 | 0.3212 | 20 (8 / 12) | 0.7969 | 0.0850 | 0.0850 |
| primary/without_GCA_000006965.1 | learned_pn_gbdt | no_supported_threshold_abstain | 52 | — | — | — | — | — | — | 17 | 0.264 [0.159, 0.402] | 21 | 0.2232 | 22 (7 / 15) | 0.8709 | 0.1237 | 0.1237 |
| primary/without_GCA_000006965.1 | learned_pu | no_supported_threshold_abstain | 52 | — | — | — | — | — | — | 18 | 0.292 [0.189, 0.426] | 23 | 0.9566 | 13 (5 / 8) | 0.5789 | 0.4714 | 0.4714 |
| primary/without_GCA_000006965.1 | learned_pu_gbdt | no_supported_threshold_abstain | 52 | — | — | — | — | — | — | 17 | 0.253 [0.146, 0.397] | 23 | 0.9951 | 19 (7 / 12) | 0.6007 | 0.4921 | 0.4921 |
| primary/without_GCA_000006965.1 | learned_pu_rank | no_supported_threshold_abstain | 52 | — | — | — | — | — | — | 18 | 0.292 [0.189, 0.426] | 23 | 0.9988 | 13 (6 / 7) | 0.5636 | 0.6110 | 0.6110 |
| swapped/without_GCA_000006965.1 | raw | no_supported_threshold_abstain | 78 | — | — | — | — | — | — | 27 | 0.333 [0.230, 0.430] | 37 | — | — | 0.6840 | 0.1938 | — |
| swapped/without_GCA_000006965.1 | stage0 | no_supported_threshold_abstain | 78 | — | — | — | — | — | — | 27 | 0.333 [0.230, 0.430] | 37 | — | — | 0.7582 | 0.1768 | — |
| swapped/without_GCA_000006965.1 | stage1 | no_supported_threshold_abstain | 78 | — | — | — | — | — | — | 27 | 0.333 [0.230, 0.430] | 37 | — | — | 0.7582 | 0.1768 | — |
| swapped/without_GCA_000006965.1 | learned_pn | no_supported_threshold_abstain | 78 | — | — | — | — | — | — | 25 | 0.321 [0.209, 0.424] | 37 | 0.3212 | 38 (13 / 25) | 0.7653 | 0.1805 | 0.1805 |
| swapped/without_GCA_000006965.1 | learned_pn_gbdt | no_supported_threshold_abstain | 78 | — | — | — | — | — | — | 26 | 0.325 [0.231, 0.414] | 32 | 0.2232 | 39 (12 / 27) | 0.8439 | 0.1854 | 0.1854 |
| swapped/without_GCA_000006965.1 | learned_pu | no_supported_threshold_abstain | 78 | — | — | — | — | — | — | 25 | 0.321 [0.209, 0.424] | 37 | 0.9566 | 18 (4 / 14) | 0.5224 | 0.4261 | 0.4261 |
| swapped/without_GCA_000006965.1 | learned_pu_gbdt | no_supported_threshold_abstain | 78 | — | — | — | — | — | — | 15 | 0.208 [0.128, 0.287] | 37 | 0.9951 | 19 (6 / 13) | 0.5235 | 0.5522 | 0.5522 |
| swapped/without_GCA_000006965.1 | learned_pu_rank | no_supported_threshold_abstain | 78 | — | — | — | — | — | — | 25 | 0.321 [0.209, 0.424] | 37 | 0.9988 | 26 (8 / 18) | 0.5128 | 0.5118 | 0.5118 |

## 2. Pairwise Contrasts on Test Genomes

| Analysis | Contrast | Thr Before | Thr After | Cal Corrected | Cal Harmed | Fixed 0.5 Corrected | Fixed 0.5 Harmed | Fixed 0.5 Net [95% CI] | Display Corrected | Display Harmed | Top-1 Gained | Top-1 Lost | Top-1 Net [95% CI] | Changed Top-1 |
|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---|---:|
| primary/all | stage1_to_learned_pu | — | — | 0 | 0 | 28 | 0 | 28 [15.00, 43.00] | 6 | 0 | 0 | 0 | 0 [0.00, 0.00] | 1 |
| primary/all | learned_pn_to_learned_pu | — | — | 0 | 0 | 24 | 0 | 24 [13.00, 37.00] | 0 | 4 | 0 | 0 | 0 [0.00, 0.00] | 2 |
| primary/all | learned_pn_gbdt_to_learned_pu_gbdt | — | — | 0 | 0 | 20 | 0 | 20 [10.00, 33.00] | 3 | 3 | 0 | 0 | 0 [0.00, 0.00] | 6 |
| primary/all | learned_pu_to_learned_pu_rank | — | — | 0 | 0 | 0 | 0 | 0 [0.00, 0.00] | 2 | 0 | 0 | 0 | 0 [0.00, 0.00] | 1 |
| swapped/all | stage1_to_learned_pu | — | — | 0 | 0 | 38 | 0 | 38 [24.00, 54.00] | 4 | 0 | 0 | 2 | -2 [-4.00, 0.00] | 3 |
| swapped/all | learned_pn_to_learned_pu | — | — | 0 | 0 | 36 | 0 | 36 [23.00, 50.00] | 0 | 9 | 0 | 0 | 0 [0.00, 0.00] | 0 |
| swapped/all | learned_pn_gbdt_to_learned_pu_gbdt | — | — | 0 | 0 | 36 | 0 | 36 [24.00, 49.00] | 3 | 9 | 1 | 12 | -11 [-18.00, -5.00] | 16 |
| swapped/all | learned_pu_to_learned_pu_rank | — | — | 0 | 0 | 0 | 0 | 0 [0.00, 0.00] | 4 | 0 | 0 | 0 | 0 [0.00, 0.00] | 1 |
| primary/without_GCA_000006965.1 | stage1_to_learned_pu | — | — | 0 | 0 | 24 | 0 | 24 [12.00, 38.00] | 5 | 0 | 0 | 0 | 0 [0.00, 0.00] | 0 |
| primary/without_GCA_000006965.1 | learned_pn_to_learned_pu | — | — | 0 | 0 | 20 | 0 | 20 [10.00, 32.00] | 0 | 3 | 0 | 0 | 0 [0.00, 0.00] | 1 |
| primary/without_GCA_000006965.1 | learned_pn_gbdt_to_learned_pu_gbdt | — | — | 0 | 0 | 16 | 0 | 16 [7.00, 27.00] | 3 | 3 | 0 | 0 | 0 [0.00, 0.00] | 5 |
| primary/without_GCA_000006965.1 | learned_pu_to_learned_pu_rank | — | — | 0 | 0 | 0 | 0 | 0 [0.00, 0.00] | 1 | 0 | 0 | 0 | 0 [0.00, 0.00] | 1 |
| swapped/without_GCA_000006965.1 | stage1_to_learned_pu | — | — | 0 | 0 | 38 | 0 | 38 [24.00, 54.00] | 4 | 0 | 0 | 2 | -2 [-4.00, 0.00] | 3 |
| swapped/without_GCA_000006965.1 | learned_pn_to_learned_pu | — | — | 0 | 0 | 36 | 0 | 36 [23.00, 50.00] | 0 | 9 | 0 | 0 | 0 [0.00, 0.00] | 0 |
| swapped/without_GCA_000006965.1 | learned_pn_gbdt_to_learned_pu_gbdt | — | — | 0 | 0 | 36 | 0 | 36 [24.00, 49.00] | 3 | 9 | 1 | 12 | -11 [-18.00, -5.00] | 16 |
| swapped/without_GCA_000006965.1 | learned_pu_to_learned_pu_rank | — | — | 0 | 0 | 0 | 0 | 0 [0.00, 0.00] | 4 | 0 | 0 | 0 | 0 [0.00, 0.00] | 1 |

## 3. Learned Arm Candidate & Abstention Counts (All Genomes and Without GCA_000006965.1)

| Scope | Arm | Cells | Executed | Candidates | Scored | Unsupported | Rescored | Contract Unavail | Ref EC Scored | Ref EC Unsupported | Unsupported Reasons |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| all | learned_pn | 22 | 22 | 1234318 | 919227 | 315091 | 916705 | 0 | 460 | 2 | EC/sequence_score_outside_training_calibration_range:822, EC/unseen_source_function_role_or_availability_profile:46, GO/profile_lacks_both_classes_in_train_and_calibration:296245, GO/sequence_score_outside_training_calibration_range:1506, GO/unseen_source_function_role_or_availability_profile:16472 |
| all | learned_pn_gbdt | 22 | 22 | 1234318 | 919227 | 315091 | 916705 | 0 | 460 | 2 | EC/sequence_score_outside_training_calibration_range:822, EC/unseen_source_function_role_or_availability_profile:46, GO/profile_lacks_both_classes_in_train_and_calibration:296245, GO/sequence_score_outside_training_calibration_range:1506, GO/unseen_source_function_role_or_availability_profile:16472 |
| all | learned_pu | 22 | 22 | 1234318 | 919227 | 315091 | 916705 | 0 | 460 | 2 | EC/sequence_score_outside_training_calibration_range:822, EC/unseen_source_function_role_or_availability_profile:46, GO/profile_lacks_both_classes_in_train_and_calibration:296245, GO/sequence_score_outside_training_calibration_range:1506, GO/unseen_source_function_role_or_availability_profile:16472 |
| all | learned_pu_gbdt | 22 | 22 | 1234318 | 919227 | 315091 | 916705 | 0 | 460 | 2 | EC/sequence_score_outside_training_calibration_range:822, EC/unseen_source_function_role_or_availability_profile:46, GO/profile_lacks_both_classes_in_train_and_calibration:296245, GO/sequence_score_outside_training_calibration_range:1506, GO/unseen_source_function_role_or_availability_profile:16472 |
| all | learned_pu_rank | 22 | 22 | 1234318 | 919227 | 315091 | 916705 | 0 | 460 | 2 | EC/sequence_score_outside_training_calibration_range:822, EC/unseen_source_function_role_or_availability_profile:46, GO/profile_lacks_both_classes_in_train_and_calibration:296245, GO/sequence_score_outside_training_calibration_range:1506, GO/unseen_source_function_role_or_availability_profile:16472 |
| without_GCA_000006965.1 | learned_pn | 21 | 21 | 1167496 | 868832 | 298664 | 866310 | 0 | 409 | 1 | EC/sequence_score_outside_training_calibration_range:758, EC/unseen_source_function_role_or_availability_profile:46, GO/profile_lacks_both_classes_in_train_and_calibration:280747, GO/sequence_score_outside_training_calibration_range:1423, GO/unseen_source_function_role_or_availability_profile:15690 |
| without_GCA_000006965.1 | learned_pn_gbdt | 21 | 21 | 1167496 | 868832 | 298664 | 866310 | 0 | 409 | 1 | EC/sequence_score_outside_training_calibration_range:758, EC/unseen_source_function_role_or_availability_profile:46, GO/profile_lacks_both_classes_in_train_and_calibration:280747, GO/sequence_score_outside_training_calibration_range:1423, GO/unseen_source_function_role_or_availability_profile:15690 |
| without_GCA_000006965.1 | learned_pu | 21 | 21 | 1167496 | 868832 | 298664 | 866310 | 0 | 409 | 1 | EC/sequence_score_outside_training_calibration_range:758, EC/unseen_source_function_role_or_availability_profile:46, GO/profile_lacks_both_classes_in_train_and_calibration:280747, GO/sequence_score_outside_training_calibration_range:1423, GO/unseen_source_function_role_or_availability_profile:15690 |
| without_GCA_000006965.1 | learned_pu_gbdt | 21 | 21 | 1167496 | 868832 | 298664 | 866310 | 0 | 409 | 1 | EC/sequence_score_outside_training_calibration_range:758, EC/unseen_source_function_role_or_availability_profile:46, GO/profile_lacks_both_classes_in_train_and_calibration:280747, GO/sequence_score_outside_training_calibration_range:1423, GO/unseen_source_function_role_or_availability_profile:15690 |
| without_GCA_000006965.1 | learned_pu_rank | 21 | 21 | 1167496 | 868832 | 298664 | 866310 | 0 | 409 | 1 | EC/sequence_score_outside_training_calibration_range:758, EC/unseen_source_function_role_or_availability_profile:46, GO/profile_lacks_both_classes_in_train_and_calibration:280747, GO/sequence_score_outside_training_calibration_range:1423, GO/unseen_source_function_role_or_availability_profile:15690 |
