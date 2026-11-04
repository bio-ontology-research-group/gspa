# Off-Growth Enzyme Discovery Report

This report summarizes the pre-registered experiment testing whether targeted gap-filling beyond single-medium biomass growth recovers experimentally verified bacterial enzymes after per-case masking.

## 1. Enzyme Recovery by Stratum

Each truth case is classified before evaluation into the on-growth stratum (deleting the truth reaction from the unmasked model drops default-medium biomass growth) or the off-growth stratum (default-medium growth persists without the truth reaction). Confidence intervals come from the genome-block bootstrap.

| Arm | Total Cases | Total Recovered | Total Recovery [95% CI] | On-Growth Cases | On-Growth Recovered | On-Growth Recovery [95% CI] | Off-Growth Cases | Off-Growth Recovered | Off-Growth Recovery [95% CI] |
| --- | ---: | ---: | --- | ---: | ---: | --- | ---: | ---: | --- |
| `A0` | 520 | 24 | 0.046 [0.018, 0.079] | 54 | 22 | 0.407 [0.229, 0.581] | 466 | 2 | 0.004 [0.000, 0.011] |
| `A1` | 520 | 56 | 0.108 [0.058, 0.160] | 54 | 32 | 0.593 [0.371, 0.742] | 466 | 24 | 0.052 [0.032, 0.071] |
| `A2` | 520 | 46 | 0.088 [0.042, 0.137] | 54 | 30 | 0.556 [0.303, 0.722] | 466 | 16 | 0.034 [0.020, 0.048] |
| `A3` | 520 | 36 | 0.069 [0.032, 0.109] | 54 | 21 | 0.389 [0.148, 0.543] | 466 | 15 | 0.032 [0.018, 0.046] |
| `A1_A2_A3` | 520 | 66 | 0.127 [0.074, 0.180] | 54 | 32 | 0.593 [0.371, 0.742] | 466 | 34 | 0.073 [0.049, 0.097] |
| `A4@10%` | 520 | 49 | 0.094 [0.045, 0.148] | 54 | 32 | 0.593 [0.371, 0.742] | 466 | 17 | 0.036 [0.019, 0.057] |
| `A4@25%` | 520 | 52 | 0.100 [0.050, 0.154] | 54 | 32 | 0.593 [0.371, 0.742] | 466 | 20 | 0.043 [0.024, 0.063] |

### Recovery by Truth Case Source Stratum

Stratification of recovery across the Price-148 truth set and the Swiss-Prot experimental enzyme truth set.

| Arm | Case Source | Cases | Recovered | Recovery [95% CI] | Off-Growth Cases | Off-Growth Recovered | Off-Growth Recovery [95% CI] |
| --- | --- | ---: | ---: | --- | ---: | ---: | --- |
| `A0` | `price148` | 148 | 0 | 0.000 [0.000, 0.000] | 146 | 0 | 0.000 [0.000, 0.000] |
| `A0` | `swissprot` | 372 | 24 | 0.065 [0.028, 0.107] | 320 | 2 | 0.006 [0.000, 0.015] |
| `A1` | `price148` | 148 | 2 | 0.014 [0.000, 0.038] | 146 | 2 | 0.014 [0.000, 0.038] |
| `A1` | `swissprot` | 372 | 54 | 0.145 [0.082, 0.211] | 320 | 22 | 0.069 [0.043, 0.095] |
| `A2` | `price148` | 148 | 3 | 0.020 [0.000, 0.042] | 146 | 3 | 0.021 [0.000, 0.043] |
| `A2` | `swissprot` | 372 | 43 | 0.116 [0.054, 0.180] | 320 | 13 | 0.041 [0.022, 0.059] |
| `A3` | `price148` | 148 | 2 | 0.014 [0.000, 0.032] | 146 | 2 | 0.014 [0.000, 0.032] |
| `A3` | `swissprot` | 372 | 34 | 0.091 [0.041, 0.143] | 320 | 13 | 0.041 [0.023, 0.058] |
| `A1_A2_A3` | `price148` | 148 | 4 | 0.027 [0.006, 0.055] | 146 | 4 | 0.027 [0.006, 0.055] |
| `A1_A2_A3` | `swissprot` | 372 | 62 | 0.167 [0.099, 0.234] | 320 | 30 | 0.094 [0.061, 0.125] |
| `A4@10%` | `price148` | 148 | 1 | 0.007 [0.000, 0.025] | 146 | 1 | 0.007 [0.000, 0.026] |
| `A4@10%` | `swissprot` | 372 | 48 | 0.129 [0.065, 0.196] | 320 | 16 | 0.050 [0.026, 0.077] |
| `A4@25%` | `price148` | 148 | 2 | 0.014 [0.000, 0.038] | 146 | 2 | 0.014 [0.000, 0.038] |
| `A4@25%` | `swissprot` | 372 | 50 | 0.134 [0.070, 0.201] | 320 | 18 | 0.056 [0.032, 0.083] |

## 2. Call Quality, Harm, and Unconfirmed Calls

Added reactions are assigned to the top-scoring masked candidate gene for a matching enzyme class. Calls on proteins with known experimental labels are scored as supported when the reaction matches the experimental enzyme class and harmful when it disagrees. Calls on proteins without experimental labels are listed as unconfirmed and are never counted as negatives.

| Arm | Gene Calls | Reaction-Only Calls | Supported Calls | Harmful Calls | Unconfirmed Calls | Harm Rate [95% CI] | Exceeds Baseline Harm Gate |
| --- | ---: | ---: | ---: | ---: | ---: | --- | --- |
| `A0` | 25961 | 44224 | 116 | 273 | 25572 | 0.702 [0.538, 0.884] | no |
| `A1` | 93172 | 289501 | 383 | 1565 | 91224 | 0.803 [0.736, 0.848] | yes |
| `A2` | 95093 | 266600 | 267 | 1287 | 93539 | 0.828 [0.685, 0.895] | yes |
| `A3` | 64511 | 162106 | 253 | 941 | 63317 | 0.788 [0.566, 0.897] | no |
| `A1_A2_A3` | 143357 | 513443 | 509 | 2301 | 140547 | 0.819 [0.721, 0.863] | yes |
| `A4@10%` | 72373 | 173408 | 287 | 1149 | 70937 | 0.800 [0.752, 0.865] | no |
| `A4@25%` | 87620 | 221247 | 310 | 1313 | 85997 | 0.809 [0.754, 0.873] | yes |

## 3. Simulation and Linear Programming Cost

Every linear program and weighted gap-fill solve is counted together with processor time, wall-clock time, and peak resident memory.

| Arm | Tests Considered | Tests Executed | Total LPs | Total Gap-Fills | CPU Seconds | Wall Seconds | Peak RSS (KiB) | LPs per Recovered Enzyme |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `A0` | 520 | 520 | 2600 | 520 | 3413.050 | 3413.939 | 3198604 | 108.3 |
| `A1` | 152360 | 152360 | 451298 | 126614 | 373871.652 | 379814.586 | 3198604 | 8058.9 |
| `A2` | 395166 | 395166 | 1564808 | 260349 | 635413.434 | 636153.174 | 3198604 | 34017.6 |
| `A3` | 164546 | 164546 | 267952 | 93644 | 197076.096 | 197099.730 | 3198604 | 7443.1 |
| `A1_A2_A3` | 712072 | 712072 | 2284058 | 480607 | 1206361.182 | 1213067.490 | 3198604 | 34606.9 |
| `A4@10%` | 712072 | 71447 | 202909 | 63349 | 191976.985 | 194721.772 | 3198604 | 4141.0 |
| `A4@25%` | 712072 | 178183 | 417364 | 131690 | 360771.569 | 363547.447 | 3198604 | 8026.2 |

## 4. Pre-Registered Statistical Tests

The paired comparisons test whether the targeted gap-filling arms recover more off-growth enzymes than default-medium biomass gap-filling alone, and whether the deterministic rule ordering retains most of the off-growth recovery within a reduced test budget.

| Hypothesis | Treatment Arm | Baseline Arm | Stratum | Paired Cases | Both Recovered | Treatment Only | Baseline Only | Neither | Raw p-value | Holm p-value | Significant |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `H1` | `A1_A2_A3` | `A0` | `off_growth` | 466 | 2 | 32 | 0 | 432 | 2.328e-10 | 2.328e-10 | yes |
| `H2_A1` | `A1` | `A0` | `off_growth` | 466 | 2 | 22 | 0 | 442 | 2.384e-07 | 7.153e-07 | yes |
| `H2_A2` | `A2` | `A0` | `off_growth` | 466 | 2 | 14 | 0 | 450 | 6.104e-05 | 0.0001221 | yes |
| `H2_A3` | `A3` | `A0` | `off_growth` | 466 | 2 | 13 | 0 | 451 | 0.0001221 | 0.0001221 | yes |

| Budget Arm | Full Arm | Off-Growth Recovered (Budget) | Off-Growth Recovered (Full) | Retention Ratio [95% CI] | Tests Executed / Considered | Budget Fraction | Meets Retention Target |
| --- | --- | ---: | ---: | --- | ---: | ---: | --- |
| `A4@25%` | `A1_A2_A3` | 20 | 34 | 0.588 [0.424, 0.750] | 178183 / 712072 | 0.250 | no |

## 5. Controls and Solver Replay

The shuffled-candidate control permutes candidate protein identifiers within each genome across fixed random seeds and reassigns calls on the added reactions. The solver replay control re-solves a deterministic sample of weighted gap-fill linear programs with the secondary solver.

| Shuffled Seed | Arm | Total Cases | Total Recovered | Total Recovery Rate | Off-Growth Cases | Off-Growth Recovered | Off-Growth Recovery Rate |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | `A1_A2_A3` | 520 | 0 | 0.000 | 466 | 0 | 0.000 |
| 2 | `A1_A2_A3` | 520 | 0 | 0.000 | 466 | 0 | 0.000 |
| 3 | `A1_A2_A3` | 520 | 0 | 0.000 | 466 | 0 | 0.000 |
| 4 | `A1_A2_A3` | 520 | 0 | 0.000 | 466 | 0 | 0.000 |
| 5 | `A1_A2_A3` | 520 | 0 | 0.000 | 466 | 0 | 0.000 |

| Primary Solver | Replay Solver | Highspy Installed | Replayed Gap-Fills | Agreed Gap-Fills | Agreement Rate | Max Relative Objective Diff |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| `glpk` | `scipy` | no | 27027 | 10999 | 0.407 | 8.915e-02 |

### Unmasked Reachability Control

Running every arm on the unmasked draft model with unmasked candidates and scoring the same way separates reaction reachability in the universal database from masking loss.

| Arm | Total Cases | Unmasked Recovered | Unmasked Recovery Rate | Off-Growth Cases | Unmasked Off-Growth Recovered | Unmasked Off-Growth Rate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `A0` | 520 | 4 | 0.008 | 466 | 3 | 0.006 |
| `A1` | 520 | 14 | 0.027 | 466 | 13 | 0.028 |
| `A2` | 520 | 11 | 0.021 | 466 | 11 | 0.024 |
| `A3` | 520 | 10 | 0.019 | 466 | 9 | 0.019 |
| `A1_A2_A3` | 520 | 19 | 0.037 | 466 | 18 | 0.039 |
| `A4@10%` | 520 | 8 | 0.015 | 466 | 7 | 0.015 |
| `A4@25%` | 520 | 10 | 0.019 | 466 | 9 | 0.019 |

### Masked-Candidate Attribution Control

Scoring the masked arm runs against the masked candidate table (where the target protein's truth enzyme claim has been deleted) yields zero recovery by construction.

| Arm | Total Cases | Masked-Candidate Recovered | Masked-Candidate Recovery Rate | Off-Growth Cases | Masked-Candidate Off-Growth Recovered | Masked-Candidate Off-Growth Rate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `A0` | 520 | 0 | 0.000 | 466 | 0 | 0.000 |
| `A1` | 520 | 0 | 0.000 | 466 | 0 | 0.000 |
| `A2` | 520 | 0 | 0.000 | 466 | 0 | 0.000 |
| `A3` | 520 | 0 | 0.000 | 466 | 0 | 0.000 |
| `A1_A2_A3` | 520 | 0 | 0.000 | 466 | 0 | 0.000 |
| `A4@10%` | 520 | 0 | 0.000 | 466 | 0 | 0.000 |
| `A4@25%` | 520 | 0 | 0.000 | 466 | 0 | 0.000 |

## 6. Unmapped or Unrouted Source Panel Compounds

Every compound from the input source panel that could not be mapped to a ModelSEED identifier with a verified elemental formula and universal transport route is listed below rather than dropped silently.

| Plate | Well | Element | Compound Name | SEED ID | Reason | Detail |
| --- | --- | --- | --- | --- | --- | --- |
| `PM1` | `A01` | `C` | Negative Control | `NA` | `negative_control` | Negative control well excluded from nutrient source panel |
| `PM1` | `B08` | `C` | D-Xylose | `cpd01422` | `no_transport_route` | No usable transport reaction for cpd01422 in universal SEED DB |
| `PM1` | `C02` | `C` | D-Galactonic Acid-γ-Lactone | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM1` | `D05` | `C` | Tween 40 | `cpd24417` | `no_transport_route` | No usable transport reaction for cpd24417 in universal SEED DB |
| `PM1` | `D08` | `C` | α-Methyl-D-Galactoside | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM1` | `E03` | `C` | D-Glucose-1-Phosphate | `cpd00501` | `no_transport_route` | No usable transport reaction for cpd00501 in universal SEED DB |
| `PM1` | `E06` | `C` | α-Hydroxy Glutaric Acid-γ-Lactone | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM1` | `F01` | `C` | Glycyl-L-Aspartic Acid | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM1` | `F06` | `C` | Bromo Succinic Acid | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM1` | `G01` | `C` | Glycyl-L-Glutamic Acid | `cpd24418` | `no_transport_route` | No usable transport reaction for cpd24418 in universal SEED DB |
| `PM1` | `G08` | `C` | N-Acetyl-ß-D-Mannosamine | `cpd27563` | `no_transport_route` | No usable transport reaction for cpd27563 in universal SEED DB |
| `PM1` | `G09` | `C` | Mono Methyl Succinate | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM1` | `G10` | `C` | Methyl Pyruvate | `cpd24420` | `no_transport_route` | No usable transport reaction for cpd24420 in universal SEED DB |
| `PM1` | `H01` | `C` | Glycyl-L-Proline | `cpd22832` | `no_transport_route` | No usable transport reaction for cpd22832 in universal SEED DB |
| `PM1` | `H03` | `C` | m-Hydroxy Phenyl Acetic Acid | `cpd03320` | `no_transport_route` | No usable transport reaction for cpd03320 in universal SEED DB |
| `PM1` | `H07` | `C` | Glucuronamide | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM1` | `H09` | `C` | L-Galactonic Acid-γ-Lactone | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM2A` | `A01` | `C` | Negative Control | `NA` | `negative_control` | Negative control well excluded from nutrient source panel |
| `PM2A` | `A02` | `C` | Chondroitin Sulfate C | `cpd26720` | `formula_lacks_element` | Formula 'null' does not contain element 'C' |
| `PM2A` | `A03` | `C` | α-Cyclodextrin | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM2A` | `A04` | `C` | ß-Cyclodextrin | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM2A` | `A05` | `C` | γ-Cyclodextrin | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM2A` | `A07` | `C` | Gelatin | `cpd11879` | `formula_lacks_element` | Formula 'null' does not contain element 'C' |
| `PM2A` | `B12` | `C` | 3-0-ß-D-Galactopyranosyl- D-Arabinose | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM2A` | `C02` | `C` | L-Glucose | `cpd24421` | `no_transport_route` | No usable transport reaction for cpd24421 in universal SEED DB |
| `PM2A` | `C05` | `C` | Maltitol | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM2A` | `C08` | `C` | 3-Methyl Glucose | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM2A` | `C09` | `C` | ß-Methyl-D-Glucuronic Acid | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM2A` | `C10` | `C` | α-Methyl-D-Mannoside | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM2A` | `C11` | `C` | ß-Methyl-D-Xyloside | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM2A` | `D03` | `C` | Sedoheptulosan | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM2A` | `D07` | `C` | Turanose | `cpd20885` | `no_transport_route` | No usable transport reaction for cpd20885 in universal SEED DB |
| `PM2A` | `D09` | `C` | N-Acetyl-D-Glucosaminitol | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM2A` | `D11` | `C` | δ-Amino Valeric Acid | `cpd17076` | `no_transport_route` | No usable transport reaction for cpd17076 in universal SEED DB |
| `PM2A` | `E04` | `C` | Citramalic Acid | `cpd00607` | `no_transport_route` | No usable transport reaction for cpd00607 in universal SEED DB |
| `PM2A` | `F01` | `C` | D-Lactic Acid Methyl Ester | `cpd24428` | `no_transport_route` | No usable transport reaction for cpd24428 in universal SEED DB |
| `PM2A` | `F03` | `C` | Melibionic Acid | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM2A` | `F05` | `C` | Oxalomalic Acid | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM2A` | `F07` | `C` | D-Ribono-1,4-Lactone | `cpd23861` | `no_transport_route` | No usable transport reaction for cpd23861 in universal SEED DB |
| `PM2A` | `F10` | `C` | Succinamic Acid | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM2A` | `F11` | `C` | D-Tartaric Acid | `cpd19018` | `no_transport_route` | No usable transport reaction for cpd19018 in universal SEED DB |
| `PM2A` | `G02` | `C` | L-Alaninamide | `cpd23030` | `no_transport_route` | No usable transport reaction for cpd23030 in universal SEED DB |
| `PM2A` | `H06` | `C` | Sec-Butylamine | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM3B` | `A01` | `N` | Negative Control | `NA` | `negative_control` | Negative control well excluded from nutrient source panel |
| `PM3B` | `D02` | `N` | N-Phthaloyl-L-Glutamic Acid | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM3B` | `D06` | `N` | N-Amylamine | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM3B` | `D07` | `N` | N-Butylamine | `cpd19969` | `no_transport_route` | No usable transport reaction for cpd19969 in universal SEED DB |
| `PM3B` | `D10` | `N` | Ethylenediamine | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM3B` | `E06` | `N` | Glucuronamide | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM3B` | `E07` | `N` | D,L-Lactamide | `cpd23860` | `no_transport_route` | No usable transport reaction for cpd23860 in universal SEED DB |
| `PM3B` | `G04` | `N` | Alloxan | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM3B` | `G06` | `N` | Parabanic Acid | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM3B` | `G09` | `N` | ε-Amino-N-Caproic Acid | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM3B` | `G10` | `N` | D,L-α-Amino-Caprylic Acid | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM3B` | `G12` | `N` | α-Amino-N- Valeric Acid | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM3B` | `H01` | `N` | Ala-Asp | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM3B` | `H03` | `N` | Ala-Glu | `cpd15388` | `no_transport_route` | No usable transport reaction for cpd15388 in universal SEED DB |
| `PM3B` | `H07` | `N` | Ala-Thr | `cpd23857` | `no_transport_route` | No usable transport reaction for cpd23857 in universal SEED DB |
| `PM3B` | `H08` | `N` | Gly-Asn | `cpd23855` | `no_transport_route` | No usable transport reaction for cpd23855 in universal SEED DB |
| `PM3B` | `H10` | `N` | Gly-Glu | `cpd24418` | `no_transport_route` | No usable transport reaction for cpd24418 in universal SEED DB |
| `PM3B` | `H12` | `N` | Met-Ala | `cpd23854` | `no_transport_route` | No usable transport reaction for cpd23854 in universal SEED DB |
| `PM4A` | `A01` | `P` | Negative Control | `NA` | `negative_control` | Negative control well excluded from nutrient source panel |
| `PM4A` | `A06` | `P` | Triethyl Phosphate | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM4A` | `A07` | `P` | Hypophosphite | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM4A` | `A08` | `P` | Adenosine-2’-monophosphate | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM4A` | `A12` | `P` | Adenosine-3’,5’-cyclic monophosphate | `cpd00446` | `no_transport_route` | No usable transport reaction for cpd00446 in universal SEED DB |
| `PM4A` | `B01` | `P` | Thiophosphate | `cpd30723` | `formula_lacks_element` | Formula 'null' does not contain element 'P' |
| `PM4A` | `B02` | `P` | Dithiophosphate | `cpd24445` | `no_transport_route` | No usable transport reaction for cpd24445 in universal SEED DB |
| `PM4A` | `B08` | `P` | Guanosine-2’-monophosphate | `cpd23708` | `no_transport_route` | No usable transport reaction for cpd23708 in universal SEED DB |
| `PM4A` | `C03` | `P` | D-Glucose-1-Phosphate | `cpd00501` | `no_transport_route` | No usable transport reaction for cpd00501 in universal SEED DB |
| `PM4A` | `C06` | `P` | D-Glucosamine-6-Phosphate | `cpd23898` | `no_transport_route` | No usable transport reaction for cpd23898 in universal SEED DB |
| `PM4A` | `C08` | `P` | Cytidine-2-monophosphate | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM4A` | `D03` | `P` | Cysteamine-S-Phosphate | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM4A` | `D05` | `P` | O-Phospho-D-Serine | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM4A` | `D12` | `P` | Uridine-3’,5’- cyclic monophosphate | `cpd24439` | `no_transport_route` | No usable transport reaction for cpd24439 in universal SEED DB |
| `PM4A` | `E01` | `P` | O-Phospho-D-Tyrosine | `cpd24440` | `no_transport_route` | No usable transport reaction for cpd24440 in universal SEED DB |
| `PM4A` | `E02` | `P` | O-Phospho-L-Tyrosine | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM4A` | `E08` | `P` | Methylene Diphosphonic Acid | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM4A` | `E09` | `P` | Thymidine-3’-monophosphate | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM4A` | `E11` | `P` | Inositol Hexaphosphate | `cpd00885` | `no_transport_route` | No usable transport reaction for cpd00885 in universal SEED DB |
| `PM4A` | `E12` | `P` | Thymidine 3’,5’-cyclic monophosphate | `cpd24443` | `no_transport_route` | No usable transport reaction for cpd24443 in universal SEED DB |
| `PM4A` | `F01` | `S` | Negative Control | `NA` | `negative_control` | Negative control well excluded from nutrient source panel |
| `PM4A` | `F03` | `S` | Sodium thiophosphate | `cpd30723` | `formula_lacks_element` | Formula 'null' does not contain element 'S' |
| `PM4A` | `F05` | `S` | Thiophosphate | `cpd30723` | `formula_lacks_element` | Formula 'null' does not contain element 'S' |
| `PM4A` | `F06` | `S` | Dithiophosphate | `cpd24445` | `no_transport_route` | No usable transport reaction for cpd24445 in universal SEED DB |
| `PM4A` | `G02` | `S` | S-Methyl-L-Cysteine | `cpd28135` | `no_transport_route` | No usable transport reaction for cpd28135 in universal SEED DB |
| `PM4A` | `G12` | `S` | L-Methionine Sulfone | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM4A` | `H01` | `S` | L-Djenkolic Acid | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM4A` | `H03` | `S` | 1-Thio-ß-D- Glucose | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
| `PM4A` | `H08` | `S` | m-Amino benzene sulfonic acid | `cpd24118` | `no_transport_route` | No usable transport reaction for cpd24118 in universal SEED DB |
| `PM4A` | `H12` | `S` | Tetramethylene Sulfone | `NA` | `unmapped_compound` | No matching ModelSEED compound ID found in reference tables |
