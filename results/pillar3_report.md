# Pillar 3: Scenario-Specific Biological Discovery & Constraint Benchmarks

## 1. Scenario B — Obligate Heteromeric Protein Complex Completion (`pillar3_complex_lso.tsv`)

| Genome | Method | Cases | Hit@1 | Recall@0.5 | Precision@0.5 | F1@0.5 | Complex Coherence |
|---|---|---:|---:|---:|---:|---:|---:|
| `ecoli` | Baseline Noisy-OR | 30 | 0.5667 | 0.0000 | 0.0000 | 0.0000 | 0.4706 |
| `ecoli` | Stage 1 CoherencePrior | 30 | 0.5667 | 0.5667 | 0.1933 | 0.2756 | 0.7938 |
| `ecoli` | Stage 1 Coherence + Operon Context | 30 | 0.9667 | 1.0000 | 0.6267 | 0.7089 | 0.7938 |
| `ecoli` | Stage 1 + Stage 2 + Stage 3 Agent | 30 | 0.9667 | 1.0000 | 0.6267 | 0.7089 | 0.7938 |
| `bsubtilis` | Baseline Noisy-OR | 30 | 0.6000 | 0.0000 | 0.0000 | 0.0000 | 0.3378 |
| `bsubtilis` | Stage 1 CoherencePrior | 30 | 0.6000 | 0.5333 | 0.1972 | 0.2757 | 0.6471 |
| `bsubtilis` | Stage 1 Coherence + Operon Context | 30 | 1.0000 | 0.9667 | 0.6306 | 0.7090 | 0.6618 |
| `bsubtilis` | Stage 1 + Stage 2 + Stage 3 Agent | 30 | 1.0000 | 1.0000 | 0.6639 | 0.7424 | 0.6618 |
| `mtb` | Baseline Noisy-OR | 30 | 0.4667 | 0.0000 | 0.0000 | 0.0000 | 0.3857 |
| `mtb` | Stage 1 CoherencePrior | 30 | 0.4667 | 0.4333 | 0.1256 | 0.1911 | 0.7797 |
| `mtb` | Stage 1 Coherence + Operon Context | 30 | 1.0000 | 0.9667 | 0.6589 | 0.7244 | 0.7797 |
| `mtb` | Stage 1 + Stage 2 + Stage 3 Agent | 30 | 1.0000 | 1.0000 | 0.6922 | 0.7578 | 0.7797 |
| `paeruginosa` | Baseline Noisy-OR | 30 | 0.6000 | 0.0000 | 0.0000 | 0.0000 | 0.3776 |
| `paeruginosa` | Stage 1 CoherencePrior | 30 | 0.6000 | 0.5333 | 0.1772 | 0.2529 | 0.6897 |
| `paeruginosa` | Stage 1 Coherence + Operon Context | 30 | 0.9667 | 0.9667 | 0.6106 | 0.6862 | 0.7126 |
| `paeruginosa` | Stage 1 + Stage 2 + Stage 3 Agent | 30 | 0.9667 | 1.0000 | 0.6439 | 0.7195 | 0.7126 |
| `synechocystis` | Baseline Noisy-OR | 25 | 0.5200 | 0.0000 | 0.0000 | 0.0000 | 0.3587 |
| `synechocystis` | Stage 1 CoherencePrior | 25 | 0.5200 | 0.4800 | 0.1813 | 0.2507 | 0.6125 |
| `synechocystis` | Stage 1 Coherence + Operon Context | 25 | 0.6400 | 0.6400 | 0.3393 | 0.4080 | 0.6125 |
| `synechocystis` | Stage 1 + Stage 2 + Stage 3 Agent | 25 | 0.9600 | 1.0000 | 0.6993 | 0.7680 | 0.6125 |
| `ecolo157` | Baseline Noisy-OR | 30 | 0.6000 | 0.0000 | 0.0000 | 0.0000 | 0.4103 |
| `ecolo157` | Stage 1 CoherencePrior | 30 | 0.6000 | 0.4667 | 0.1639 | 0.2302 | 0.7416 |
| `ecolo157` | Stage 1 Coherence + Operon Context | 30 | 1.0000 | 0.9667 | 0.6622 | 0.7279 | 0.7416 |
| `ecolo157` | Stage 1 + Stage 2 + Stage 3 Agent | 30 | 1.0000 | 1.0000 | 0.6956 | 0.7613 | 0.7416 |
| `ALL (Weighted Mean)` | Baseline Noisy-OR | 175 | 0.5600 | 0.0000 | 0.0000 | 0.0000 | 0.3901 |
| `ALL (Weighted Mean)` | Stage 1 CoherencePrior | 175 | 0.5600 | 0.5029 | 0.1729 | 0.2459 | 0.7107 |
| `ALL (Weighted Mean)` | Stage 1 Coherence + Operon Context | 175 | 0.9371 | 0.9257 | 0.5951 | 0.6680 | 0.7170 |
| `ALL (Weighted Mean)` | Stage 1 + Stage 2 + Stage 3 Agent | 175 | 0.9829 | 1.0000 | 0.6694 | 0.7423 | 0.7170 |

## 2. Scenario C — Multicellular Eukaryotic Anatomical & Developmental Constraints (`pillar3_eukaryotic_anatomy.tsv`)

| Species | Taxon ID | Clade | Proteins | E1 Violations | E2 Direct GO TC Violations | E3 Stage 1 | E4/E5 Stage 2 Violations | E4/E5 Elimination (%) | Borderline Rescued (E3 vs E1) | Rescue Rate (%) |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `scerevisiae` | `NCBITaxon_4932` | Fungi | 1500 | 303 | 270 | 0 | 0 | 100.00% | 2/2 (vs 0) | 100.00% |
| `celegans` | `NCBITaxon_6239` | Nematoda | 1500 | 589 | 571 | 203 | 0 | 100.00% | 48/80 (vs 0) | 60.00% |
| `dmelanogaster` | `NCBITaxon_7227` | Arthropoda | 1500 | 1017 | 1014 | 500 | 0 | 100.00% | 162/210 (vs 0) | 77.14% |
| `ppatens` | `NCBITaxon_3218` | Bryophyta (non-vascular moss) | 155 | 42 | 39 | 0 | 0 | 100.00% | 0/0 (vs 0) | 0.00% |
| `athaliana` | `NCBITaxon_3702` | Magnoliopsida (flowering plant) | 1500 | 284 | 261 | 0 | 0 | 100.00% | 40/49 (vs 0) | 81.63% |
| `mmusculus` | `NCBITaxon_10090` | Mammalia (vertebrate) | 1500 | 142 | 130 | 0 | 0 | 100.00% | 207/234 (vs 0) | 88.46% |

## 3. Semantic Ancestor Relaxation (`mode=relax`) vs. Naive Removal (`mode=remove`) (`pillar3_relax_vs_remove.tsv`)

| Species | Taxon ID | Mode | Propagated CAFA Fmax | Propagated Precision | Propagated Recall | Mean IC |
|---|---|---|---:|---:|---:|---:|
| `scerevisiae` | `NCBITaxon_4932` | `E1_noisy_or` | 0.5923 | 0.7313 | 0.4977 | 12.1850 |
| `scerevisiae` | `NCBITaxon_4932` | `E4_stage2_remove` | 0.6175 | 0.9976 | 0.4472 | 11.3869 |
| `scerevisiae` | `NCBITaxon_4932` | `E5_stage2_relax` | 0.6460 | 0.9516 | 0.4890 | 11.4514 |
| `celegans` | `NCBITaxon_6239` | `E1_noisy_or` | 0.5869 | 0.7338 | 0.4890 | 12.2576 |
| `celegans` | `NCBITaxon_6239` | `E4_stage2_remove` | 0.5890 | 0.9331 | 0.4303 | 11.6888 |
| `celegans` | `NCBITaxon_6239` | `E5_stage2_relax` | 0.6090 | 0.8699 | 0.4685 | 11.6304 |
| `dmelanogaster` | `NCBITaxon_7227` | `E1_noisy_or` | 0.5190 | 0.7578 | 0.3946 | 12.4393 |
| `dmelanogaster` | `NCBITaxon_7227` | `E4_stage2_remove` | 0.4981 | 0.9250 | 0.3408 | 11.9431 |
| `dmelanogaster` | `NCBITaxon_7227` | `E5_stage2_relax` | 0.5250 | 0.8787 | 0.3744 | 11.7460 |
| `ppatens` | `NCBITaxon_3218` | `E1_noisy_or` | 0.7029 | 0.7470 | 0.6638 | 10.0503 |
| `ppatens` | `NCBITaxon_3218` | `E4_stage2_remove` | 0.7425 | 0.9574 | 0.6064 | 9.4594 |
| `ppatens` | `NCBITaxon_3218` | `E5_stage2_relax` | 0.7633 | 0.9166 | 0.6540 | 9.7482 |
| `athaliana` | `NCBITaxon_3702` | `E1_noisy_or` | 0.6245 | 0.7297 | 0.5459 | 11.4824 |
| `athaliana` | `NCBITaxon_3702` | `E4_stage2_remove` | 0.6500 | 0.9814 | 0.4859 | 10.6137 |
| `athaliana` | `NCBITaxon_3702` | `E5_stage2_relax` | 0.6756 | 0.9226 | 0.5330 | 10.8170 |
| `mmusculus` | `NCBITaxon_10090` | `E1_noisy_or` | 0.5244 | 0.8032 | 0.3893 | 12.3880 |
| `mmusculus` | `NCBITaxon_10090` | `E4_stage2_remove` | 0.5228 | 0.8778 | 0.3723 | 12.1648 |
| `mmusculus` | `NCBITaxon_10090` | `E5_stage2_relax` | 0.5317 | 0.8556 | 0.3857 | 12.1561 |
