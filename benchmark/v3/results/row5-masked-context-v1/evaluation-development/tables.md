# Row-5 Masked-Context Evaluation (`scope=development`)

Chosen context window (bp): `10000`

## 1. Contrasts (Overall and Per Claim Set)

| Contrast | Claim set | Cells | Improved | Harmed | Unchanged | Discordant | Critical k | p (improve) | Harm elsewhere | Verdict |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `C1_stage0_to_gaps` | **all** | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `diamond` | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `C1b_gaps_to_gapsfix` | **all** | 58 | 22 | 0 | 36 | 22 | 16 | 2.384e-07 | 0 | systematic improvement |
| `C1b_gaps_to_gapsfix` | `diamond` | 58 | 22 | 0 | 36 | 22 | 16 | 2.384e-07 | 0 | systematic improvement |
| `C1c_stage0_to_gapsfix` | **all** | 58 | 22 | 0 | 36 | 22 | 16 | 2.384e-07 | 0 | systematic improvement |
| `C1c_stage0_to_gapsfix` | `diamond` | 58 | 22 | 0 | 36 | 22 | 16 | 2.384e-07 | 0 | systematic improvement |
| `C5_universal_w5000` | **all** | 58 | 0 | 7 | 51 | 7 | 7 | 1 | 2 | harm |
| `C5_universal_w5000` | `diamond` | 58 | 0 | 7 | 51 | 7 | 7 | 1 | 2 | harm |
| `C5_pathway_w5000` | **all** | 58 | 0 | 6 | 52 | 6 | 6 | 1 | 1 | harm |
| `C5_pathway_w5000` | `diamond` | 58 | 0 | 6 | 52 | 6 | 6 | 1 | 1 | harm |
| `C5_organism_w5000` | **all** | 58 | 0 | 6 | 52 | 6 | 6 | 1 | 1 | harm |
| `C5_organism_w5000` | `diamond` | 58 | 0 | 6 | 52 | 6 | 6 | 1 | 1 | harm |
| `C5total_universal_w5000` | **all** | 58 | 16 | 0 | 42 | 16 | 12 | 1.526e-05 | 1 | systematic improvement |
| `C5total_universal_w5000` | `diamond` | 58 | 16 | 0 | 42 | 16 | 12 | 1.526e-05 | 1 | systematic improvement |
| `C5_pathway_vs_universal_w5000` | **all** | 58 | 1 | 0 | 57 | 1 | - | 0.5 | 1 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | `diamond` | 58 | 1 | 0 | 57 | 1 | - | 0.5 | 1 | no detectable effect |
| `C5_organism_vs_universal_w5000` | **all** | 58 | 1 | 0 | 57 | 1 | - | 0.5 | 1 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `diamond` | 58 | 1 | 0 | 57 | 1 | - | 0.5 | 1 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | **all** | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `diamond` | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w10000` | **all** | 58 | 0 | 7 | 51 | 7 | 7 | 1 | 2 | harm |
| `C5_universal_w10000` | `diamond` | 58 | 0 | 7 | 51 | 7 | 7 | 1 | 2 | harm |
| `C5_pathway_w10000` | **all** | 58 | 0 | 6 | 52 | 6 | 6 | 1 | 1 | harm |
| `C5_pathway_w10000` | `diamond` | 58 | 0 | 6 | 52 | 6 | 6 | 1 | 1 | harm |
| `C5_organism_w10000` | **all** | 58 | 0 | 6 | 52 | 6 | 6 | 1 | 1 | harm |
| `C5_organism_w10000` | `diamond` | 58 | 0 | 6 | 52 | 6 | 6 | 1 | 1 | harm |
| `C5total_universal_w10000` | **all** | 58 | 16 | 0 | 42 | 16 | 12 | 1.526e-05 | 1 | systematic improvement |
| `C5total_universal_w10000` | `diamond` | 58 | 16 | 0 | 42 | 16 | 12 | 1.526e-05 | 1 | systematic improvement |
| `C5_pathway_vs_universal_w10000` | **all** | 58 | 1 | 0 | 57 | 1 | - | 0.5 | 1 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | `diamond` | 58 | 1 | 0 | 57 | 1 | - | 0.5 | 1 | no detectable effect |
| `C5_organism_vs_universal_w10000` | **all** | 58 | 1 | 0 | 57 | 1 | - | 0.5 | 1 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `diamond` | 58 | 1 | 0 | 57 | 1 | - | 0.5 | 1 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | **all** | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `diamond` | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w20000` | **all** | 58 | 0 | 7 | 51 | 7 | 7 | 1 | 2 | harm |
| `C5_universal_w20000` | `diamond` | 58 | 0 | 7 | 51 | 7 | 7 | 1 | 2 | harm |
| `C5_pathway_w20000` | **all** | 58 | 0 | 6 | 52 | 6 | 6 | 1 | 1 | harm |
| `C5_pathway_w20000` | `diamond` | 58 | 0 | 6 | 52 | 6 | 6 | 1 | 1 | harm |
| `C5_organism_w20000` | **all** | 58 | 0 | 6 | 52 | 6 | 6 | 1 | 1 | harm |
| `C5_organism_w20000` | `diamond` | 58 | 0 | 6 | 52 | 6 | 6 | 1 | 1 | harm |
| `C5total_universal_w20000` | **all** | 58 | 16 | 0 | 42 | 16 | 12 | 1.526e-05 | 1 | systematic improvement |
| `C5total_universal_w20000` | `diamond` | 58 | 16 | 0 | 42 | 16 | 12 | 1.526e-05 | 1 | systematic improvement |
| `C5_pathway_vs_universal_w20000` | **all** | 58 | 1 | 0 | 57 | 1 | - | 0.5 | 1 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `diamond` | 58 | 1 | 0 | 57 | 1 | - | 0.5 | 1 | no detectable effect |
| `C5_organism_vs_universal_w20000` | **all** | 58 | 1 | 0 | 57 | 1 | - | 0.5 | 1 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `diamond` | 58 | 1 | 0 | 57 | 1 | - | 0.5 | 1 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | **all** | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `diamond` | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal` | **all** | 58 | 0 | 7 | 51 | 7 | 7 | 1 | 2 | harm |
| `C5_universal` | `diamond` | 58 | 0 | 7 | 51 | 7 | 7 | 1 | 2 | harm |
| `C5_pathway` | **all** | 58 | 0 | 6 | 52 | 6 | 6 | 1 | 1 | harm |
| `C5_pathway` | `diamond` | 58 | 0 | 6 | 52 | 6 | 6 | 1 | 1 | harm |
| `C5_organism` | **all** | 58 | 0 | 6 | 52 | 6 | 6 | 1 | 1 | harm |
| `C5_organism` | `diamond` | 58 | 0 | 6 | 52 | 6 | 6 | 1 | 1 | harm |
| `C5_pathway_vs_universal` | **all** | 58 | 1 | 0 | 57 | 1 | - | 0.5 | 1 | no detectable effect |
| `C5_pathway_vs_universal` | `diamond` | 58 | 1 | 0 | 57 | 1 | - | 0.5 | 1 | no detectable effect |
| `C5_organism_vs_universal` | **all** | 58 | 1 | 0 | 57 | 1 | - | 0.5 | 1 | no detectable effect |
| `C5_organism_vs_universal` | `diamond` | 58 | 1 | 0 | 57 | 1 | - | 0.5 | 1 | no detectable effect |
| `C5_organism_vs_pathway` | **all** | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway` | `diamond` | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | **all** | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `diamond` | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_gapsfix` | **all** | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_gapsfix` | `diamond` | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | **all** | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `diamond` | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | **all** | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | `diamond` | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | **all** | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `diamond` | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | **all** | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | `diamond` | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | **all** | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `diamond` | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | **all** | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | `diamond` | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | **all** | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `diamond` | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | **all** | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `diamond` | 58 | 0 | 0 | 58 | 0 | - | 1 | 0 | no detectable effect |

## 2. Harm Elsewhere and Cost

| Contrast | Improved | Harmed | Harm elsewhere | Unconfirmed calls created | Delta calls >= 0.5 | Changed scores | Bootstrap 95% CI | Wilcoxon p (1-sided) |
|---|---:|---:|---:|---:|---:|---:|---|---:|
| `C1_stage0_to_gaps` | 0 | 0 | 0 | 0 | 1856 | 115099 | [0.0, 0.0] | 1 |
| `C1b_gaps_to_gapsfix` | 22 | 0 | 0 | 8 | 68 | 506 | [10.0, 40.0] | 0.03394 |
| `C1c_stage0_to_gapsfix` | 22 | 0 | 0 | 8 | 1924 | 115167 | [10.0, 40.0] | 0.03394 |
| `C5_universal_w5000` | 0 | 7 | 2 | 0 | -346 | 428 | [-17.0, 0.0] | 1 |
| `C5_pathway_w5000` | 0 | 6 | 1 | 0 | -324 | 407 | [-16.0, 0.0] | 1 |
| `C5_organism_w5000` | 0 | 6 | 1 | 0 | -172 | 229 | [-16.0, 0.0] | 1 |
| `C5total_universal_w5000` | 16 | 0 | 1 | 8 | -278 | 752 | [9.0, 25.0] | 0.03394 |
| `C5_pathway_vs_universal_w5000` | 1 | 0 | 1 | 0 | 22 | 163 | [0.0, 3.0] | 1 |
| `C5_organism_vs_universal_w5000` | 1 | 0 | 1 | 0 | 174 | 258 | [0.0, 3.0] | 1 |
| `C5_organism_vs_pathway_w5000` | 0 | 0 | 0 | 0 | 152 | 283 | [0.0, 0.0] | 1 |
| `C5_universal_w10000` | 0 | 7 | 2 | 0 | -321 | 420 | [-17.0, 0.0] | 1 |
| `C5_pathway_w10000` | 0 | 6 | 1 | 0 | -302 | 408 | [-16.0, 0.0] | 1 |
| `C5_organism_w10000` | 0 | 6 | 1 | 0 | -166 | 223 | [-16.0, 0.0] | 1 |
| `C5total_universal_w10000` | 16 | 0 | 1 | 8 | -253 | 746 | [9.0, 25.0] | 0.03394 |
| `C5_pathway_vs_universal_w10000` | 1 | 0 | 1 | 0 | 19 | 190 | [0.0, 3.0] | 1 |
| `C5_organism_vs_universal_w10000` | 1 | 0 | 1 | 0 | 155 | 274 | [0.0, 3.0] | 1 |
| `C5_organism_vs_pathway_w10000` | 0 | 0 | 0 | 0 | 136 | 306 | [0.0, 0.0] | 1 |
| `C5_universal_w20000` | 0 | 7 | 2 | 0 | -329 | 436 | [-17.0, 0.0] | 1 |
| `C5_pathway_w20000` | 0 | 6 | 1 | 0 | -304 | 414 | [-16.0, 0.0] | 1 |
| `C5_organism_w20000` | 0 | 6 | 1 | 0 | -166 | 250 | [-16.0, 0.0] | 1 |
| `C5total_universal_w20000` | 16 | 0 | 1 | 8 | -261 | 767 | [9.0, 25.0] | 0.03394 |
| `C5_pathway_vs_universal_w20000` | 1 | 0 | 1 | 0 | 25 | 302 | [0.0, 3.0] | 1 |
| `C5_organism_vs_universal_w20000` | 1 | 0 | 1 | 0 | 163 | 313 | [0.0, 3.0] | 1 |
| `C5_organism_vs_pathway_w20000` | 0 | 0 | 0 | 0 | 138 | 367 | [0.0, 0.0] | 1 |
| `C5_universal` | 0 | 7 | 2 | 0 | -321 | 420 | [-17.0, 0.0] | 1 |
| `C5_pathway` | 0 | 6 | 1 | 0 | -302 | 408 | [-16.0, 0.0] | 1 |
| `C5_organism` | 0 | 6 | 1 | 0 | -166 | 223 | [-16.0, 0.0] | 1 |
| `C5_pathway_vs_universal` | 1 | 0 | 1 | 0 | 19 | 190 | [0.0, 3.0] | 1 |
| `C5_organism_vs_universal` | 1 | 0 | 1 | 0 | 155 | 274 | [0.0, 3.0] | 1 |
| `C5_organism_vs_pathway` | 0 | 0 | 0 | 0 | 136 | 306 | [0.0, 0.0] | 1 |
| `S1_stage0_to_sham_gaps` | 0 | 0 | 0 | 0 | 1837 | 115080 | [0.0, 0.0] | 1 |
| `S1_sham_gaps_to_sham_gapsfix` | 0 | 0 | 0 | 0 | 55 | 477 | [0.0, 0.0] | 1 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | 0 | 0 | 0 | 0 | -337 | 417 | [0.0, 0.0] | 1 |
| `S1_stage0_to_sham_ctx_universal_w5000` | 0 | 0 | 0 | 0 | 1555 | 115135 | [0.0, 0.0] | 1 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | 0 | 0 | 0 | 0 | -312 | 409 | [0.0, 0.0] | 1 |
| `S1_stage0_to_sham_ctx_universal_w10000` | 0 | 0 | 0 | 0 | 1580 | 115135 | [0.0, 0.0] | 1 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | 0 | 0 | 0 | 0 | -320 | 425 | [0.0, 0.0] | 1 |
| `S1_stage0_to_sham_ctx_universal_w20000` | 0 | 0 | 0 | 0 | 1572 | 115135 | [0.0, 0.0] | 1 |
| `S1_sham_gaps_to_sham_ctx_universal` | 0 | 0 | 0 | 0 | -312 | 409 | [0.0, 0.0] | 1 |
| `U1_unmasked_gaps_vs_gaps` | 0 | 0 | 0 | 0 | 32 | 32 | [0.0, 0.0] | 1 |

## 3. Strata Breakdown

| Contrast | Stratum axis | Stratum | Cells | Improved | Harmed | Unchanged | Discordant |
|---|---|---|---:|---:|---:|---:|---:|
| `C1_stage0_to_gaps` | `domain` | `archaea` | 54 | 0 | 0 | 54 | 0 |
| `C1_stage0_to_gaps` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C1_stage0_to_gaps` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C1_stage0_to_gaps` | `gap_kind` | `maskable_gap` | 56 | 0 | 0 | 56 | 0 |
| `C1_stage0_to_gaps` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C1_stage0_to_gaps` | `pathway_robust` | `robust` | 38 | 0 | 0 | 38 | 0 |
| `C1b_gaps_to_gapsfix` | `domain` | `archaea` | 54 | 20 | 0 | 34 | 20 |
| `C1b_gaps_to_gapsfix` | `domain` | `bacteria` | 4 | 2 | 0 | 2 | 2 |
| `C1b_gaps_to_gapsfix` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C1b_gaps_to_gapsfix` | `gap_kind` | `maskable_gap` | 56 | 22 | 0 | 34 | 22 |
| `C1b_gaps_to_gapsfix` | `pathway_robust` | `marginal` | 20 | 5 | 0 | 15 | 5 |
| `C1b_gaps_to_gapsfix` | `pathway_robust` | `robust` | 38 | 17 | 0 | 21 | 17 |
| `C1c_stage0_to_gapsfix` | `domain` | `archaea` | 54 | 20 | 0 | 34 | 20 |
| `C1c_stage0_to_gapsfix` | `domain` | `bacteria` | 4 | 2 | 0 | 2 | 2 |
| `C1c_stage0_to_gapsfix` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C1c_stage0_to_gapsfix` | `gap_kind` | `maskable_gap` | 56 | 22 | 0 | 34 | 22 |
| `C1c_stage0_to_gapsfix` | `pathway_robust` | `marginal` | 20 | 5 | 0 | 15 | 5 |
| `C1c_stage0_to_gapsfix` | `pathway_robust` | `robust` | 38 | 17 | 0 | 21 | 17 |
| `C5_universal_w5000` | `domain` | `archaea` | 54 | 0 | 7 | 47 | 7 |
| `C5_universal_w5000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_universal_w5000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_universal_w5000` | `gap_kind` | `maskable_gap` | 56 | 0 | 7 | 49 | 7 |
| `C5_universal_w5000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_universal_w5000` | `pathway_robust` | `robust` | 38 | 0 | 7 | 31 | 7 |
| `C5_pathway_w5000` | `domain` | `archaea` | 54 | 0 | 6 | 48 | 6 |
| `C5_pathway_w5000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_pathway_w5000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_pathway_w5000` | `gap_kind` | `maskable_gap` | 56 | 0 | 6 | 50 | 6 |
| `C5_pathway_w5000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_pathway_w5000` | `pathway_robust` | `robust` | 38 | 0 | 6 | 32 | 6 |
| `C5_organism_w5000` | `domain` | `archaea` | 54 | 0 | 6 | 48 | 6 |
| `C5_organism_w5000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_organism_w5000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_organism_w5000` | `gap_kind` | `maskable_gap` | 56 | 0 | 6 | 50 | 6 |
| `C5_organism_w5000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_organism_w5000` | `pathway_robust` | `robust` | 38 | 0 | 6 | 32 | 6 |
| `C5total_universal_w5000` | `domain` | `archaea` | 54 | 14 | 0 | 40 | 14 |
| `C5total_universal_w5000` | `domain` | `bacteria` | 4 | 2 | 0 | 2 | 2 |
| `C5total_universal_w5000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5total_universal_w5000` | `gap_kind` | `maskable_gap` | 56 | 16 | 0 | 40 | 16 |
| `C5total_universal_w5000` | `pathway_robust` | `marginal` | 20 | 5 | 0 | 15 | 5 |
| `C5total_universal_w5000` | `pathway_robust` | `robust` | 38 | 11 | 0 | 27 | 11 |
| `C5_pathway_vs_universal_w5000` | `domain` | `archaea` | 54 | 1 | 0 | 53 | 1 |
| `C5_pathway_vs_universal_w5000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_pathway_vs_universal_w5000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_pathway_vs_universal_w5000` | `gap_kind` | `maskable_gap` | 56 | 1 | 0 | 55 | 1 |
| `C5_pathway_vs_universal_w5000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_pathway_vs_universal_w5000` | `pathway_robust` | `robust` | 38 | 1 | 0 | 37 | 1 |
| `C5_organism_vs_universal_w5000` | `domain` | `archaea` | 54 | 1 | 0 | 53 | 1 |
| `C5_organism_vs_universal_w5000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_organism_vs_universal_w5000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_organism_vs_universal_w5000` | `gap_kind` | `maskable_gap` | 56 | 1 | 0 | 55 | 1 |
| `C5_organism_vs_universal_w5000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_organism_vs_universal_w5000` | `pathway_robust` | `robust` | 38 | 1 | 0 | 37 | 1 |
| `C5_organism_vs_pathway_w5000` | `domain` | `archaea` | 54 | 0 | 0 | 54 | 0 |
| `C5_organism_vs_pathway_w5000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_organism_vs_pathway_w5000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_organism_vs_pathway_w5000` | `gap_kind` | `maskable_gap` | 56 | 0 | 0 | 56 | 0 |
| `C5_organism_vs_pathway_w5000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_organism_vs_pathway_w5000` | `pathway_robust` | `robust` | 38 | 0 | 0 | 38 | 0 |
| `C5_universal_w10000` | `domain` | `archaea` | 54 | 0 | 7 | 47 | 7 |
| `C5_universal_w10000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_universal_w10000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_universal_w10000` | `gap_kind` | `maskable_gap` | 56 | 0 | 7 | 49 | 7 |
| `C5_universal_w10000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_universal_w10000` | `pathway_robust` | `robust` | 38 | 0 | 7 | 31 | 7 |
| `C5_pathway_w10000` | `domain` | `archaea` | 54 | 0 | 6 | 48 | 6 |
| `C5_pathway_w10000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_pathway_w10000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_pathway_w10000` | `gap_kind` | `maskable_gap` | 56 | 0 | 6 | 50 | 6 |
| `C5_pathway_w10000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_pathway_w10000` | `pathway_robust` | `robust` | 38 | 0 | 6 | 32 | 6 |
| `C5_organism_w10000` | `domain` | `archaea` | 54 | 0 | 6 | 48 | 6 |
| `C5_organism_w10000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_organism_w10000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_organism_w10000` | `gap_kind` | `maskable_gap` | 56 | 0 | 6 | 50 | 6 |
| `C5_organism_w10000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_organism_w10000` | `pathway_robust` | `robust` | 38 | 0 | 6 | 32 | 6 |
| `C5total_universal_w10000` | `domain` | `archaea` | 54 | 14 | 0 | 40 | 14 |
| `C5total_universal_w10000` | `domain` | `bacteria` | 4 | 2 | 0 | 2 | 2 |
| `C5total_universal_w10000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5total_universal_w10000` | `gap_kind` | `maskable_gap` | 56 | 16 | 0 | 40 | 16 |
| `C5total_universal_w10000` | `pathway_robust` | `marginal` | 20 | 5 | 0 | 15 | 5 |
| `C5total_universal_w10000` | `pathway_robust` | `robust` | 38 | 11 | 0 | 27 | 11 |
| `C5_pathway_vs_universal_w10000` | `domain` | `archaea` | 54 | 1 | 0 | 53 | 1 |
| `C5_pathway_vs_universal_w10000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_pathway_vs_universal_w10000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_pathway_vs_universal_w10000` | `gap_kind` | `maskable_gap` | 56 | 1 | 0 | 55 | 1 |
| `C5_pathway_vs_universal_w10000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_pathway_vs_universal_w10000` | `pathway_robust` | `robust` | 38 | 1 | 0 | 37 | 1 |
| `C5_organism_vs_universal_w10000` | `domain` | `archaea` | 54 | 1 | 0 | 53 | 1 |
| `C5_organism_vs_universal_w10000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_organism_vs_universal_w10000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_organism_vs_universal_w10000` | `gap_kind` | `maskable_gap` | 56 | 1 | 0 | 55 | 1 |
| `C5_organism_vs_universal_w10000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_organism_vs_universal_w10000` | `pathway_robust` | `robust` | 38 | 1 | 0 | 37 | 1 |
| `C5_organism_vs_pathway_w10000` | `domain` | `archaea` | 54 | 0 | 0 | 54 | 0 |
| `C5_organism_vs_pathway_w10000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_organism_vs_pathway_w10000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_organism_vs_pathway_w10000` | `gap_kind` | `maskable_gap` | 56 | 0 | 0 | 56 | 0 |
| `C5_organism_vs_pathway_w10000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_organism_vs_pathway_w10000` | `pathway_robust` | `robust` | 38 | 0 | 0 | 38 | 0 |
| `C5_universal_w20000` | `domain` | `archaea` | 54 | 0 | 7 | 47 | 7 |
| `C5_universal_w20000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_universal_w20000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_universal_w20000` | `gap_kind` | `maskable_gap` | 56 | 0 | 7 | 49 | 7 |
| `C5_universal_w20000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_universal_w20000` | `pathway_robust` | `robust` | 38 | 0 | 7 | 31 | 7 |
| `C5_pathway_w20000` | `domain` | `archaea` | 54 | 0 | 6 | 48 | 6 |
| `C5_pathway_w20000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_pathway_w20000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_pathway_w20000` | `gap_kind` | `maskable_gap` | 56 | 0 | 6 | 50 | 6 |
| `C5_pathway_w20000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_pathway_w20000` | `pathway_robust` | `robust` | 38 | 0 | 6 | 32 | 6 |
| `C5_organism_w20000` | `domain` | `archaea` | 54 | 0 | 6 | 48 | 6 |
| `C5_organism_w20000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_organism_w20000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_organism_w20000` | `gap_kind` | `maskable_gap` | 56 | 0 | 6 | 50 | 6 |
| `C5_organism_w20000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_organism_w20000` | `pathway_robust` | `robust` | 38 | 0 | 6 | 32 | 6 |
| `C5total_universal_w20000` | `domain` | `archaea` | 54 | 14 | 0 | 40 | 14 |
| `C5total_universal_w20000` | `domain` | `bacteria` | 4 | 2 | 0 | 2 | 2 |
| `C5total_universal_w20000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5total_universal_w20000` | `gap_kind` | `maskable_gap` | 56 | 16 | 0 | 40 | 16 |
| `C5total_universal_w20000` | `pathway_robust` | `marginal` | 20 | 5 | 0 | 15 | 5 |
| `C5total_universal_w20000` | `pathway_robust` | `robust` | 38 | 11 | 0 | 27 | 11 |
| `C5_pathway_vs_universal_w20000` | `domain` | `archaea` | 54 | 1 | 0 | 53 | 1 |
| `C5_pathway_vs_universal_w20000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_pathway_vs_universal_w20000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_pathway_vs_universal_w20000` | `gap_kind` | `maskable_gap` | 56 | 1 | 0 | 55 | 1 |
| `C5_pathway_vs_universal_w20000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_pathway_vs_universal_w20000` | `pathway_robust` | `robust` | 38 | 1 | 0 | 37 | 1 |
| `C5_organism_vs_universal_w20000` | `domain` | `archaea` | 54 | 1 | 0 | 53 | 1 |
| `C5_organism_vs_universal_w20000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_organism_vs_universal_w20000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_organism_vs_universal_w20000` | `gap_kind` | `maskable_gap` | 56 | 1 | 0 | 55 | 1 |
| `C5_organism_vs_universal_w20000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_organism_vs_universal_w20000` | `pathway_robust` | `robust` | 38 | 1 | 0 | 37 | 1 |
| `C5_organism_vs_pathway_w20000` | `domain` | `archaea` | 54 | 0 | 0 | 54 | 0 |
| `C5_organism_vs_pathway_w20000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_organism_vs_pathway_w20000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_organism_vs_pathway_w20000` | `gap_kind` | `maskable_gap` | 56 | 0 | 0 | 56 | 0 |
| `C5_organism_vs_pathway_w20000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_organism_vs_pathway_w20000` | `pathway_robust` | `robust` | 38 | 0 | 0 | 38 | 0 |
| `C5_universal` | `domain` | `archaea` | 54 | 0 | 7 | 47 | 7 |
| `C5_universal` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_universal` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_universal` | `gap_kind` | `maskable_gap` | 56 | 0 | 7 | 49 | 7 |
| `C5_universal` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_universal` | `pathway_robust` | `robust` | 38 | 0 | 7 | 31 | 7 |
| `C5_pathway` | `domain` | `archaea` | 54 | 0 | 6 | 48 | 6 |
| `C5_pathway` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_pathway` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_pathway` | `gap_kind` | `maskable_gap` | 56 | 0 | 6 | 50 | 6 |
| `C5_pathway` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_pathway` | `pathway_robust` | `robust` | 38 | 0 | 6 | 32 | 6 |
| `C5_organism` | `domain` | `archaea` | 54 | 0 | 6 | 48 | 6 |
| `C5_organism` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_organism` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_organism` | `gap_kind` | `maskable_gap` | 56 | 0 | 6 | 50 | 6 |
| `C5_organism` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_organism` | `pathway_robust` | `robust` | 38 | 0 | 6 | 32 | 6 |
| `C5_pathway_vs_universal` | `domain` | `archaea` | 54 | 1 | 0 | 53 | 1 |
| `C5_pathway_vs_universal` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_pathway_vs_universal` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_pathway_vs_universal` | `gap_kind` | `maskable_gap` | 56 | 1 | 0 | 55 | 1 |
| `C5_pathway_vs_universal` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_pathway_vs_universal` | `pathway_robust` | `robust` | 38 | 1 | 0 | 37 | 1 |
| `C5_organism_vs_universal` | `domain` | `archaea` | 54 | 1 | 0 | 53 | 1 |
| `C5_organism_vs_universal` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_organism_vs_universal` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_organism_vs_universal` | `gap_kind` | `maskable_gap` | 56 | 1 | 0 | 55 | 1 |
| `C5_organism_vs_universal` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_organism_vs_universal` | `pathway_robust` | `robust` | 38 | 1 | 0 | 37 | 1 |
| `C5_organism_vs_pathway` | `domain` | `archaea` | 54 | 0 | 0 | 54 | 0 |
| `C5_organism_vs_pathway` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `C5_organism_vs_pathway` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `C5_organism_vs_pathway` | `gap_kind` | `maskable_gap` | 56 | 0 | 0 | 56 | 0 |
| `C5_organism_vs_pathway` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `C5_organism_vs_pathway` | `pathway_robust` | `robust` | 38 | 0 | 0 | 38 | 0 |
| `S1_stage0_to_sham_gaps` | `domain` | `archaea` | 54 | 0 | 0 | 54 | 0 |
| `S1_stage0_to_sham_gaps` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `S1_stage0_to_sham_gaps` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `S1_stage0_to_sham_gaps` | `gap_kind` | `maskable_gap` | 56 | 0 | 0 | 56 | 0 |
| `S1_stage0_to_sham_gaps` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `S1_stage0_to_sham_gaps` | `pathway_robust` | `robust` | 38 | 0 | 0 | 38 | 0 |
| `S1_sham_gaps_to_sham_gapsfix` | `domain` | `archaea` | 54 | 0 | 0 | 54 | 0 |
| `S1_sham_gaps_to_sham_gapsfix` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `S1_sham_gaps_to_sham_gapsfix` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `S1_sham_gaps_to_sham_gapsfix` | `gap_kind` | `maskable_gap` | 56 | 0 | 0 | 56 | 0 |
| `S1_sham_gaps_to_sham_gapsfix` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `S1_sham_gaps_to_sham_gapsfix` | `pathway_robust` | `robust` | 38 | 0 | 0 | 38 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `domain` | `archaea` | 54 | 0 | 0 | 54 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `gap_kind` | `maskable_gap` | 56 | 0 | 0 | 56 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `pathway_robust` | `robust` | 38 | 0 | 0 | 38 | 0 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `domain` | `archaea` | 54 | 0 | 0 | 54 | 0 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `gap_kind` | `maskable_gap` | 56 | 0 | 0 | 56 | 0 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `pathway_robust` | `robust` | 38 | 0 | 0 | 38 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `domain` | `archaea` | 54 | 0 | 0 | 54 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `gap_kind` | `maskable_gap` | 56 | 0 | 0 | 56 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `pathway_robust` | `robust` | 38 | 0 | 0 | 38 | 0 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `domain` | `archaea` | 54 | 0 | 0 | 54 | 0 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `gap_kind` | `maskable_gap` | 56 | 0 | 0 | 56 | 0 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `pathway_robust` | `robust` | 38 | 0 | 0 | 38 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `domain` | `archaea` | 54 | 0 | 0 | 54 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `gap_kind` | `maskable_gap` | 56 | 0 | 0 | 56 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `pathway_robust` | `robust` | 38 | 0 | 0 | 38 | 0 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `domain` | `archaea` | 54 | 0 | 0 | 54 | 0 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `gap_kind` | `maskable_gap` | 56 | 0 | 0 | 56 | 0 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `pathway_robust` | `robust` | 38 | 0 | 0 | 38 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal` | `domain` | `archaea` | 54 | 0 | 0 | 54 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal` | `gap_kind` | `maskable_gap` | 56 | 0 | 0 | 56 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal` | `pathway_robust` | `robust` | 38 | 0 | 0 | 38 | 0 |
| `U1_unmasked_gaps_vs_gaps` | `domain` | `archaea` | 54 | 0 | 0 | 54 | 0 |
| `U1_unmasked_gaps_vs_gaps` | `domain` | `bacteria` | 4 | 0 | 0 | 4 | 0 |
| `U1_unmasked_gaps_vs_gaps` | `gap_kind` | `already_gap` | 2 | 0 | 0 | 2 | 0 |
| `U1_unmasked_gaps_vs_gaps` | `gap_kind` | `maskable_gap` | 56 | 0 | 0 | 56 | 0 |
| `U1_unmasked_gaps_vs_gaps` | `pathway_robust` | `marginal` | 20 | 0 | 0 | 20 | 0 |
| `U1_unmasked_gaps_vs_gaps` | `pathway_robust` | `robust` | 38 | 0 | 0 | 38 | 0 |

## 4. Learned-Stage Candidate Counts

| Arm | Candidates | Rescored | Unsupported | Unchanged unsupported | Retracted | Scored only |
|---|---:|---:|---:|---:|---:|---:|
| - | 0 | 0 | 0 | 0 | 0 | 0 |

## 5. Masking & Build Audit Summary

| Check | Count |
|---|---:|
| `total_cells` | 58 |
| `complete_cells` | 58 |
| `candidate_universe_preserved_cells` | 58 |
| `masked_gap_present_cells` | 58 |
| `gaps_ctxjar_identical_cells` | 58 |
| `already_gap_unmasked_total` | 2 |
| `already_gap_unmasked_matches` | 2 |
