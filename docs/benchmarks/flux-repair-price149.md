# Flux-certified repair and fitness-validated enzymes

## Question

Does positive-cost lexicographic parsimonious flux balance analysis (pFBA) eliminate the zero-flux reaction padding produced by negative-cost gap filling, does competitive locus attribution improve recovery of fitness-validated bacterial enzymes in Price-149, and how does reconstruction accuracy compare against curated BiGG models?

Paper: Section 2.6, Table 5 (in latex `tab:price149`), and Figure 6a to 6c (`paper/main.tex`).

## Datasets

### 108-genome stoichiometric invariant audit

108 genome-scale metabolic models (GEMs) reconstructed across diverse bacterial phyla, evaluated for steady-state flux feasibility, mass and charge balance, and reaction padding. Source data: `results/pillar2_invariant_audit_summary.tsv`.

### Curated 6 BiGG reference GEMs

Six gold-standard, manually curated genome-scale metabolic reconstructions from the BiGG database:

- *Escherichia coli* K-12 MG1655 (iML1515)
- *Salmonella enterica* LT2 (STM_v1_0)
- *Klebsiella pneumoniae* MGH 78578 (iYL1228)
- *Pseudomonas putida* KT2440 (iJN1463)
- *Staphylococcus aureus* USA300 (iYS854)
- *Bacillus subtilis* 168 (iYO844)

Source data: `results/pillar2_curated6_gems_real.tsv`.

### Price-149 fitness-validated enzyme benchmark

148 divergent bacterial enzymes across 22 host genomes whose physiological and biochemical functions were established experimentally through random-barcoded transposon mutagenesis (RB-TnSeq, Price et al. 2018). These enzymes represent difficult functional assignments where sequence homology to characterized enzymes is distant or ambiguous; 55 of the 148 enzymes participate in alternative catabolic pathways inactive on minimal glucose. Source data: `results/pillar2_price149_all6_real.tsv`.

## Protocol

### Stoichiometric audit

Four gap-filling formulations are compared:

1. **Raw Sequence Thresholding:** Enzymes with DeepProZyme confidence \(\ge 0.50\).
2. **Negative-Cost Reaction-Indicator MILP (METEOR v8):** Assigns negative linear objective coefficients to candidate reactions above an addition threshold.
3. **Grow-Support Positive-Cost LP:** Positive-cost linear programming without lexicographic minimization.
4. **GSPA StoichiometricLocusRefiner (SLR):** Positive-cost lexicographic pFBA requiring that every added reaction carry steady-state flux in an elementally balanced network.

Metrics: mean active reactions, flux-carrying active reactions, zero-flux active reactions, added reactions, flux-certified added reactions, uncertified padding rate (%), elementally imbalanced reactions, dead-end metabolites, and LP biomass feasibility (%).

### Curated BiGG evaluation

Evaluates CarveMe (2018), Reconstructor (2023), and METEOR v8 (2026) against curated BiGG models on reaction precision, recall, and \(F_1\); EC precision, recall, and \(F_1\); and steady-state growth feasibility.

### Price-149 evaluation

Evaluates six standalone EC predictors: CLEAN (2023), DeepProZyme (2023), EnzBERT (2023), GraphEC (2024), MAPred (2025), and TopEC (2025), alongside ECrecer (2023) and DCGM / NeSy-Adapter (2024).

For each predictor, the evaluation compares:

- Standalone baseline
- Single-predictor GSPA under minimal glucose pFBA
- Single-predictor GSPA with cognate-substrate stoichiometric probing
- Multi-predictor GSPA combining Stage 0 pooling and competitive locus attribution

Metrics: Top-1 count and accuracy, Top-5 count and accuracy, micro and macro \(F_{\max}\), mean reciprocal rank (MRR), net improved versus degraded enzymes, and two-sided exact McNemar test \(p\)-values.

## Baselines

Negative-cost MILP (METEOR v8), CarveMe, Reconstructor, six standalone neural EC predictors, ECrecer, and DCGM.

## Reproduce

The complete Price-149 evaluation is executed by `scripts/eval_price149_real_all6.py`:

```bash
uv run --with numpy --with scipy --with pandas \
    python scripts/eval_price149_real_all6.py
```

This script imports `scripts/stoichiometric_locus_refiner.py`, computes metrics for all six predictors and ablations, and writes `results/pillar2_price149_all6_real.tsv`.

The stoichiometric invariant audit and BiGG comparisons are executed by `benchmark/unified/pillar2_fba_lro_price149/run_pillar2_all.py`:

```bash
uv run --with numpy --with scipy \
    python benchmark/unified/pillar2_fba_lro_price149/run_pillar2_all.py \
        --run-audit --run-bigg --output-dir results/
```

## Result tables

- `results/pillar2_invariant_audit_summary.tsv`: columns `method`, `n_genomes`, `mean_active_rxns`, `mean_flux_carrying_active`, `mean_zero_flux_active`, `mean_added_rxns`, `mean_flux_certified_added`, `mean_uncertified_padding`, `uncertified_padding_pct`, `mean_imbalanced_added_or_active`, `mean_deadend_active`, `lp_biomass_feasible_pct`. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar2_invariant_audit_summary.tsv))
- `results/pillar2_curated6_gems_real.tsv`: columns `organism`, `curated_gem`, `method`, `rxn_precision`, `rxn_recall`, `rxn_f1`, `ec_precision`, `ec_recall`, `ec_f1`, `fba_grows`. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar2_curated6_gems_real.tsv))
- `results/pillar2_price149_all6_real.tsv`: columns `baseline`, `method`, `n_proteins`, `top1_count`, `top1_acc`, `top1_ci95_lo`, `top1_ci95_hi`, `top5_count`, `top5_acc`, `fmax`, `fmax_dataset_macro`, `fmax_dataset_micro`, `best_tau`, `f1_at_05`, `mrr`, `n_improved_vs_base`, `n_degraded_vs_base`, `mcnemar_p`. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar2_price149_all6_real.tsv))

### 108-genome stoichiometric invariant audit

Values transcribed directly from `results/pillar2_invariant_audit_summary.tsv`:

| Method | Genomes | Active Rxns | Flux-Carrying | Zero-Flux | Added Rxns | Flux-Certified Added | Uncertified Padding Rate | Imbalanced Rxns | Dead-End Metabolites | Biomass Feasible |
|---|---|---|---|---|---|---|---|---|---|---|
| Raw Sequence (DPZ \(\ge 0.5\)) | 108 | 2,809.5 | 0.0 | 2,809.5 | 0.0 | 0.0 | 0.0% | 2,621.7 | 9.75 | 0.0% |
| METEOR v8 (Negative-Cost MILP) | 108 | 3,121.0 | 446.4 | 2,674.6 | 313.35 | 51.01 | **83.72%** | 2,597.8 | 2.36 | 100.0% |
| Grow-Support LP (Positive-Cost) | 108 | 2,584.8 | 446.4 | 2,138.4 | 81.31 | 81.31 | 0.0% | 0.0 | 2.23 | 100.0% |
| **GSPA StoichiometricLocusRefiner** | **108** | **2,815.2** | **300.6** | **2,514.6** | **5.7** | **5.7** | **0.0%** | **0.0** | **0.0** | **100.0%** |

### Curated 6 BiGG models evaluation (6-GEM means)

Values transcribed from `results/pillar2_curated6_gems_real.tsv`:

| Method | Rxn Precision | Rxn Recall | Rxn \(F_1\) | EC Precision | EC Recall | EC \(F_1\) | Growth Feasible |
|---|---|---|---|---|---|---|---|
| CarveMe (2018) | **0.6693** | **0.7305** | **0.6955** | **0.8217** | **0.8927** | **0.8536** | **6/6** |
| Reconstructor (2023) | ModelSEED | Namespace | (672.3 rxns) | 0.4387 | 0.8787 | 0.5844 | 4/6 |
| METEOR v8 (2026) | 0.1655 | 0.2900 | 0.2091 | 0.4383 | 0.8842 | 0.5851 | **6/6** |

### Price-149 fitness-validated enzyme recovery

Values transcribed from `results/pillar2_price149_all6_real.tsv` across the 148 proteins:

| Baseline | Method Configuration | Top-1 Count | Top-1 Acc. [95% CI] | Top-5 Count | Top-5 Acc. | \(F_{\max}\) | MRR | Improved / Degraded | McNemar \(p\) |
|---|---|---|---|---|---|---|---|---|---|
| CLEAN (2023) | Baseline | 60 | 0.4054 [0.324, 0.487] | 90 | 0.6081 | 0.5295 | 0.4999 | 0 / 0 | 1.0 |
| CLEAN (2023) | + METEOR v8 (Iter 2) | 61 | 0.4122 [0.331, 0.493] | 91 | 0.6149 | 0.5343 | 0.5056 | 1 / 0 | 1.0 |
| CLEAN (2023) | + GSPA Single-Pred (Cognate Substrate) | 63 | 0.4257 [0.345, 0.507] | 88 | 0.5946 | 0.5214 | 0.5028 | 3 / 0 | 0.25 |
| CLEAN (2023) | **+ Full GSPA Multi-Pred (Iter 3)** | **66** | **0.4459 [0.365, 0.527]** | **102** | **0.6892** | **0.6084** | **0.5629** | **7 / 1** | **0.070** |
| DeepProZyme (2023) | Baseline | 49 | 0.3311 [0.250, 0.405] | 101 | 0.6824 | 0.5574 | 0.4904 | 0 / 0 | 1.0 |
| DeepProZyme (2023) | + METEOR v8 (Iter 2) | 49 | 0.3311 [0.257, 0.405] | 99 | 0.6689 | 0.5427 | 0.4798 | 1 / 1 | 1.0 |
| DeepProZyme (2023) | + GSPA Single-Pred (Cognate Substrate) | 51 | 0.3446 [0.270, 0.426] | 107 | 0.7230 | 0.5772 | 0.5078 | 5 / 3 | 0.73 |
| DeepProZyme (2023) | **+ Full GSPA Multi-Pred (Iter 3)** | **59** | **0.3986 [0.317, 0.473]** | **110** | **0.7432** | **0.6136** | **0.5526** | **13 / 3** | **0.021** |
| EnzBERT (2023) | Baseline | 53 | 0.3581 [0.277, 0.439] | 84 | 0.5676 | 0.4786 | 0.4422 | 0 / 0 | 1.0 |
| EnzBERT (2023) | + METEOR v8 (Iter 2) | 52 | 0.3514 [0.277, 0.432] | 74 | 0.5000 | 0.4451 | 0.4173 | 0 / 1 | 1.0 |
| EnzBERT (2023) | **+ Stage 1 Noisy-OR Ensemble** | **62** | **0.4189 [0.338, 0.500]** | **95** | **0.6419** | **0.5535** | **0.5176** | **9 / 0** | **0.0039** |
| GraphEC (2024) | Baseline | 18 | 0.1216 [0.068, 0.176] | 52 | 0.3514 | 0.2634 | 0.2161 | 0 / 0 | 1.0 |
| GraphEC (2024) | + METEOR v8 (Iter 2) | 18 | 0.1216 [0.068, 0.176] | 40 | 0.2703 | 0.2250 | 0.1919 | 0 / 0 | 1.0 |
| GraphEC (2024) | **+ Stage 1 Noisy-OR Ensemble** | **48** | **0.3243 [0.243, 0.399]** | **83** | **0.5608** | **0.4756** | **0.4324** | **30 / 0** | **< 10^-6** |
| MAPred (2025) | Baseline | 50 | 0.3378 [0.264, 0.412] | 99 | 0.6689 | 0.5534 | 0.4938 | 0 / 0 | 1.0 |
| MAPred (2025) | + METEOR v8 (Iter 2) | 52 | 0.3514 [0.277, 0.426] | 98 | 0.6622 | 0.5588 | 0.5010 | 2 / 0 | 0.50 |
| MAPred (2025) | + GSPA Single-Pred (Cognate Substrate) | 56 | 0.3784 [0.297, 0.460] | 103 | 0.6959 | 0.5746 | 0.5201 | 8 / 2 | 0.109 |
| MAPred (2025) | **+ Full GSPA Multi-Pred (Iter 3)** | **73** | **0.4932 [0.412, 0.574]** | **100** | **0.6757** | **0.6178** | **0.5822** | **31 / 8** | **0.00029** |
| TopEC (2025) | Baseline | 5 | 0.0338 [0.007, 0.068] | 5 | 0.0338 | 0.0338 | 0.0338 | 0 / 0 | 1.0 |
| TopEC (2025) | + METEOR v8 (Iter 2) | 5 | 0.0338 [0.007, 0.068] | 5 | 0.0338 | 0.0353 | 0.0338 | 0 / 0 | 1.0 |
| TopEC (2025) | **+ Full GSPA Multi-Pred (Iter 3)** | **39** | **0.2635 [0.196, 0.338]** | **76** | **0.5135** | **0.4210** | **0.3785** | **34 / 0** | **< 10^-6** |
| DCGM / NeSy-Adapter (2024) | Baseline (Standalone) | 60 | 0.4054 [0.324, 0.486] | 92 | 0.6216 | 0.5380 | 0.5042 | 0 / 0 | 1.0 |
| ECrecer (2023) | Baseline (Standalone) | 38 | 0.2568 [0.189, 0.331] | 69 | 0.4662 | 0.3890 | 0.3480 | 0 / 0 | 1.0 |

## Main findings

1. **Elimination of zero-flux reaction padding:** Negative-cost MILP gap filling (METEOR v8) adds an average of 313.35 reactions per genome, of which 262.34 (83.72%) carry zero steady-state flux, accumulating disconnected dead-end reactions and 2,597.8 elementally imbalanced reactions. GSPA StoichiometricLocusRefiner adds an average of 5.7 reactions per genome, of which 100% (5.7) carry steady-state flux in an elementally balanced network (0.0% uncertified padding, 0.0 imbalanced reactions, and 100.0% biomass feasibility).
2. **BiGG reconstruction fidelity:** CarveMe maintains high precision on curated BiGG models (reaction \(F_1 = 0.6955\), EC \(F_1 = 0.8536\)). Negative-cost MILP degrades reaction \(F_1\) to 0.2091 because unconstrained negative addition rewards incentivize accumulating hundreds of unverified reactions.
3. **Fitness-validated enzyme discovery:** On 148 divergent bacterial enzymes from Price-149, single-medium MILP gap filling fails or degrades accuracy because 55 enzymes catalyze alternative catabolic reactions inactive on minimal glucose. Stage 2B cognate-substrate probing and competitive locus attribution improve Top-1 enzyme discovery across all top predictors:
    - **MAPred:** Rises from 50 (33.8%) standalone to 56 (37.8%) with single-predictor GSPA, and to 73 (49.3%, \(p = 0.00029\)) with multi-predictor GSPA.
    - **CLEAN:** Rises from 60 (40.5%) standalone to 63 (42.6%) single-predictor, and to 66 (44.6%) multi-predictor.
    - **DeepProZyme:** Rises from 49 (33.1%) standalone to 51 (34.5%) single-predictor, and to 59 (39.9%, \(p = 0.021\)) multi-predictor.
    - **EnzBERT and GraphEC:** Stage 0 multi-predictor pooling provides massive initial gains on lower-recall models (e.g. GraphEC rises from 18 to 48 correct; \(p < 10^{-6}\)).
