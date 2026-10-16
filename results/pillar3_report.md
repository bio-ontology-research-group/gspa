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
| `bsubtilis` | Stage 1 + Stage 2 + Stage 3 Agent | 30 | 1.0000 | 0.9667 | 0.6306 | 0.7090 | 0.6618 |
| `mtb` | Baseline Noisy-OR | 30 | 0.4667 | 0.0000 | 0.0000 | 0.0000 | 0.3857 |
| `mtb` | Stage 1 CoherencePrior | 30 | 0.4667 | 0.4333 | 0.1256 | 0.1911 | 0.7797 |
| `mtb` | Stage 1 Coherence + Operon Context | 30 | 1.0000 | 0.9667 | 0.6589 | 0.7244 | 0.7797 |
| `mtb` | Stage 1 + Stage 2 + Stage 3 Agent | 30 | 1.0000 | 0.9667 | 0.6589 | 0.7244 | 0.7797 |
| `paeruginosa` | Baseline Noisy-OR | 30 | 0.6000 | 0.0000 | 0.0000 | 0.0000 | 0.3776 |
| `paeruginosa` | Stage 1 CoherencePrior | 30 | 0.6000 | 0.5333 | 0.1772 | 0.2529 | 0.6897 |
| `paeruginosa` | Stage 1 Coherence + Operon Context | 30 | 0.9667 | 0.9667 | 0.6106 | 0.6862 | 0.7126 |
| `paeruginosa` | Stage 1 + Stage 2 + Stage 3 Agent | 30 | 0.9667 | 0.9667 | 0.6106 | 0.6862 | 0.7126 |
| `synechocystis` | Baseline Noisy-OR | 25 | 0.5200 | 0.0000 | 0.0000 | 0.0000 | 0.3587 |
| `synechocystis` | Stage 1 CoherencePrior | 25 | 0.5200 | 0.4800 | 0.1813 | 0.2507 | 0.6125 |
| `synechocystis` | Stage 1 Coherence + Operon Context | 25 | 0.6400 | 0.6400 | 0.3393 | 0.4080 | 0.6125 |
| `synechocystis` | Stage 1 + Stage 2 + Stage 3 Agent | 25 | 0.6800 | 0.6800 | 0.3793 | 0.4480 | 0.6125 |
| `ecolo157` | Baseline Noisy-OR | 30 | 0.6000 | 0.0000 | 0.0000 | 0.0000 | 0.4103 |
| `ecolo157` | Stage 1 CoherencePrior | 30 | 0.6000 | 0.4667 | 0.1639 | 0.2302 | 0.7416 |
| `ecolo157` | Stage 1 Coherence + Operon Context | 30 | 1.0000 | 0.9667 | 0.6622 | 0.7279 | 0.7416 |
| `ecolo157` | Stage 1 + Stage 2 + Stage 3 Agent | 30 | 1.0000 | 0.9667 | 0.6622 | 0.7279 | 0.7416 |
| `ALL (Weighted Mean)` | Baseline Noisy-OR | 175 | 0.5600 | 0.0000 | 0.0000 | 0.0000 | 0.3901 |
| `ALL (Weighted Mean)` | Stage 1 CoherencePrior | 175 | 0.5600 | 0.5029 | 0.1729 | 0.2459 | 0.7107 |
| `ALL (Weighted Mean)` | Stage 1 Coherence + Operon Context | 175 | 0.9371 | 0.9257 | 0.5951 | 0.6680 | 0.7170 |
| `ALL (Weighted Mean)` | Stage 1 + Stage 2 + Stage 3 Agent | 175 | 0.9429 | 0.9314 | 0.6009 | 0.6737 | 0.7170 |

## 2. Scenario C — Multicellular Eukaryotic Anatomical & Developmental Constraints (`pillar3_eukaryote_anatomy.tsv`)

| Species | Taxon ID | Clade | Proteins | E1 Violations | E2 Direct GO TC Violations | E3 Stage 1 | E4/E5 Stage 2 Violations | E4/E5 Elimination (%) | Borderline Rescued (E3 vs E1) | Rescue Rate (%) |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `scerevisiae` | `NCBITaxon_4932` | Fungi | 1500 | 66 | 0 | 0 | 0 | 100.00% | 2/3 (vs 0) | 66.67% |
| `celegans` | `NCBITaxon_6239` | Nematoda | 1500 | 389 | 190 | 201 | 0 | 100.00% | 44/66 (vs 0) | 66.67% |
| `dmelanogaster` | `NCBITaxon_7227` | Arthropoda | 1500 | 894 | 426 | 559 | 0 | 100.00% | 140/198 (vs 0) | 70.71% |
| `ppatens` | `NCBITaxon_3218` | Bryophyta (non-vascular moss) | 155 | 3 | 0 | 0 | 0 | 100.00% | 0/0 (vs 0) | 0.00% |
| `athaliana` | `NCBITaxon_3702` | Magnoliopsida (flowering plant) | 1500 | 59 | 26 | 1 | 0 | 100.00% | 62/69 (vs 0) | 89.86% |
| `mmusculus` | `NCBITaxon_10090` | Mammalia (vertebrate) | 1500 | 41 | 0 | 0 | 0 | 100.00% | 254/282 (vs 0) | 90.07% |

## 3. Semantic Ancestor Relaxation (`mode=relax`) vs. Naive Removal (`mode=remove`) (`pillar3_relax_vs_remove.tsv`)

| Species | Taxon ID | Mode | Propagated CAFA Fmax | Propagated Precision | Propagated Recall | Mean IC |
|---|---|---|---:|---:|---:|---:|
| `scerevisiae` | `NCBITaxon_4932` | `E1_noisy_or` | 0.9033 | 0.8259 | 0.9966 | 12.4792 |
| `scerevisiae` | `NCBITaxon_4932` | `E4_stage2_remove` | 0.9604 | 0.9989 | 0.9249 | 12.0152 |
| `scerevisiae` | `NCBITaxon_4932` | `E5_stage2_relax` | 0.9592 | 0.9848 | 0.9349 | 12.0734 |
| `celegans` | `NCBITaxon_6239` | `E1_noisy_or` | 0.9015 | 0.8406 | 0.9720 | 12.5216 |
| `celegans` | `NCBITaxon_6239` | `E4_stage2_remove` | 0.9567 | 0.9887 | 0.9268 | 12.1475 |
| `celegans` | `NCBITaxon_6239` | `E5_stage2_relax` | 0.9628 | 0.9786 | 0.9476 | 12.0943 |
| `dmelanogaster` | `NCBITaxon_7227` | `E1_noisy_or` | 0.9076 | 0.8793 | 0.9378 | 12.8646 |
| `dmelanogaster` | `NCBITaxon_7227` | `E4_stage2_remove` | 0.9238 | 0.9968 | 0.8608 | 12.5510 |
| `dmelanogaster` | `NCBITaxon_7227` | `E5_stage2_relax` | 0.9395 | 0.9933 | 0.8913 | 12.3485 |
| `ppatens` | `NCBITaxon_3218` | `E1_noisy_or` | 0.8949 | 0.8107 | 0.9986 | 10.1007 |
| `ppatens` | `NCBITaxon_3218` | `E4_stage2_remove` | 0.9828 | 0.9963 | 0.9697 | 9.6394 |
| `ppatens` | `NCBITaxon_3218` | `E5_stage2_relax` | 0.9792 | 0.9787 | 0.9796 | 9.7905 |
| `athaliana` | `NCBITaxon_3702` | `E1_noisy_or` | 0.8884 | 0.8196 | 0.9699 | 11.7586 |
| `athaliana` | `NCBITaxon_3702` | `E4_stage2_remove` | 0.9497 | 0.9986 | 0.9054 | 11.1288 |
| `athaliana` | `NCBITaxon_3702` | `E5_stage2_relax` | 0.9547 | 0.9884 | 0.9232 | 11.1768 |
| `mmusculus` | `NCBITaxon_10090` | `E1_noisy_or` | 0.9271 | 0.9022 | 0.9534 | 13.0142 |
| `mmusculus` | `NCBITaxon_10090` | `E4_stage2_remove` | 0.9672 | 0.9874 | 0.9478 | 12.8367 |
| `mmusculus` | `NCBITaxon_10090` | `E5_stage2_relax` | 0.9678 | 0.9814 | 0.9545 | 12.8466 |
