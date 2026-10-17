# BacDive phenotype prediction

## Question

Do systems-level metabolic pathway pooling, continuous pathway coupling, and subunit-complete enzyme rules predict experimentally measured bacterial growth and diagnostic enzyme phenotypes across unseen genera in BacDive, and can partial phenotypic observations from a strain refine its remaining unassayed pathways?

Paper: Section 2.7, Table 6 (in latex `tab:bacdive`), Figure 6d, and Supplementary Table S13 (`paper/main.tex`).

## Datasets

### BacDive and ProTraits experimental cohorts

Evaluations are performed against curated phenotypic records from the BacDive bacterial diversity database and ProTraits:

- **Carbon-source growth panel:** 2,332 experimental growth assays across 529 bacterial strains spanning diverse carbon substrates (e.g. D-glucose, acetate, succinate, citrate, lactose, maltose). Source data: `results/pillar3_bacdive_track1_summary.tsv`.
- **Diagnostic enzyme panel (type strains):** 17,496 experimental diagnostic enzyme assays across 3,017 bacterial type strains covering standard microbiological assays (e.g. catalase, oxidase, urease, gelatinase, beta-galactosidase). Source data: `results/pillar3_bacdive_track1_summary.tsv`.
- **Extended 10-EC diagnostic panel:** 50,075 strain-assay pairs across 6,313 bacterial strains for 10 diagnostic enzymes.
- **Heteromeric urease (*EC 3.5.1.5*):** 5,026 bacterial strains with experimental urease assays, evaluating the impact of multi-subunit modeling (requiring the UreC alpha catalytic subunit, UreA/UreB structural partners, and UreDFG accessory chaperones).
- **108-genome catabolic panel:** 1,622 experimental growth assays across 108 bacterial genomes.

## Protocol

### Two operational tracks

- **Track 1 (Blind Genome-to-Phenotype Prediction):** Zero phenotype labels provided at inference time. Evaluates both unsupervised mechanistic reconstructions and 5-fold genus-blocked cross-validation. In genus-blocked CV, every test strain belongs to a taxonomic genus completely unseen during training (assigned deterministically by MD5 hashing of genus names), preventing phylogenetic memorization.
- **Track 2 (Phenotype-Assisted Dark-Matter Gene Attribution):** 30% of each strain's experimental assays are randomly revealed as partial phenotypic evidence, and the pipeline predicts the remaining held-out 70% using an unsupervised Bayesian update across coupled metabolic pathways.

### Methods evaluated

- **Unsupervised mechanistic reconstruction pipelines:** CarveMe v1.5.2 (Machado et al. 2018), ModelSEED v2.6.1 (Henry et al. 2010), and gapseq v1.2 (Zimmermann et al. 2021). These produce discrete binary predictions \(\{0, 1\}\).
- **Supervised phenotype classifiers (5-fold genus-blocked CV):** Traitar L1-regularized SVM (Weimann et al. 2016), BacDive-AI Random Forest (Schober et al. 2025), and DiASPora GBDT (Koblitz et al. 2025).
- **METEOR v8 (Negative-Cost MILP):** Reaction-indicator MILP applied to DeepProZyme annotations.
- **GSPA Track 1:** Blind continuous pathway pooling combining calibrated gapseq, ModelSEED, and CarveMe predictions with unsupervised cross-assay pathway co-occurrence coupling.
- **GSPA + DiASPora Hybrid:** Blends GSPA continuous mechanistic pathway posteriors with DiASPora gradient-boosted decision tree outputs.
- **GSPA Track 2:** Unsupervised Bayesian updating conditioned on the 30% revealed assay anchor.

### Metrics

Precision, Recall, \(F_1\) score (with 95% bootstrap confidence intervals), Accuracy, Balanced Accuracy, ROC-AUC (which equals balanced accuracy for binary step predictors), discordant pairs (\(b_{01}, b_{10}\)), two-sided exact McNemar test \(p\)-values, and per-genome runtime.

## Baselines

CarveMe, ModelSEED, gapseq, Traitar, BacDive-AI, DiASPora, and METEOR v6 to v8.

## Reproduce

The benchmark script `benchmark/unified/pillar3_bacdive/run_bacdive_tracks1_2.py` coordinates both Track 1 and Track 2 evaluations:

```bash
uv run --with numpy --with scipy --with scikit-learn \
    python benchmark/unified/pillar3_bacdive/run_bacdive_tracks1_2.py \
        --output-dir results/
```

This script reads the genus-blocked CV partitions (`suites7_8_bacdive_manifest.tsv`), fits fold-local cross-validation models, performs the Track 2 Bayesian updates, and writes `results/pillar3_bacdive_track1_summary.tsv`, `results/pillar3_bacdive_track2_summary.tsv`, and `results/pillar3_bacdive_per_fold_cv.tsv`.

## Result tables

- `results/pillar3_bacdive_track1_summary.tsv`: columns `dataset`, `method`, `n_assays`, `TP`, `FP`, `FN`, `TN`, `precision`, `recall`, `f1`, `accuracy`, `balanced_accuracy`, `roc_auc`, `mcnemar_b_improved`, `mcnemar_c_degraded`, `mcnemar_p`. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar3_bacdive_track1_summary.tsv))
- `results/pillar3_bacdive_track2_summary.tsv`: columns for held-out 70% evaluation across methods. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar3_bacdive_track2_summary.tsv))
- `results/pillar3_bacdive_per_fold_cv.tsv`: per-fold metrics across the 5 genus-blocked cross-validation splits. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar3_bacdive_per_fold_cv.tsv))

### Track 1 results (Blind and Genus-Blocked Cross-Validation)

Values transcribed directly from `results/pillar3_bacdive_track1_summary.tsv`:

| Assay Panel | Method | Strains / Assays | Precision | Recall | \(F_1\) Score | Balanced Accuracy | ROC-AUC | Runtime / Genome |
|---|---|---|---|---|---|---|---|---|
| Carbon Growth | CarveMe v1.5.2 (Binary) | 529 / 2,332 | 0.6935 | 0.4073 | 0.5132 | 0.6233 | 0.6233 | 28.4 s |
| Carbon Growth | ModelSEED v2.6.1 (Binary) | 529 / 2,332 | 0.8576 | 0.5200 | 0.6474 | 0.7214 | 0.7214 | 42.0 s |
| Carbon Growth | gapseq v1.2 (Binary) | 529 / 2,332 | 0.8071 | 0.7836 | 0.7952 | 0.8082 | 0.8082 | 185.0 s |
| Carbon Growth | **GSPA Track 1 (100% Blind)** | 529 / 2,332 | 0.8071 | 0.7836 | 0.7952 | 0.8082 | **0.8567** | 1.1 s |
| Carbon Growth | Traitar L1-SVM (Genus-CV) | 529 / 2,332 | 0.8472 | 0.8673 | 0.8571 | 0.8638 | 0.9185 | 0.05 s |
| Carbon Growth | BacDive-AI RF (Genus-CV) | 529 / 2,332 | 0.8433 | 0.8464 | 0.8448 | 0.8530 | 0.9222 | 0.05 s |
| Carbon Growth | DiASPora GBDT (Genus-CV) | 529 / 2,332 | 0.8475 | 0.8491 | 0.8483 | 0.8564 | 0.9268 | 0.05 s |
| Carbon Growth | **GSPA + DiASPora Hybrid** | 529 / 2,332 | **0.8766** | **0.8655** | **0.8710** | **0.8783** | **0.9567** | 1.2 s |
| Diagnostic Enzymes | CarveMe v1.5.2 (Binary) | 3,017 / 17,496 | 0.7739 | 0.2996 | 0.4319 | 0.5774 | 0.5774 | 28.4 s |
| Diagnostic Enzymes | ModelSEED v2.6.1 (Binary) | 3,017 / 17,496 | 0.7854 | 0.4779 | 0.5943 | 0.6310 | 0.6310 | 42.0 s |
| Diagnostic Enzymes | gapseq v1.2 (Binary) | 3,017 / 17,496 | 0.8393 | 0.7204 | 0.7753 | 0.7461 | 0.7461 | 185.0 s |
| Diagnostic Enzymes | **GSPA Track 1 (100% Blind)** | 3,017 / 17,496 | 0.8393 | 0.7204 | 0.7753 | 0.7461 | **0.7520** | 1.1 s |
| Diagnostic Enzymes | Traitar L1-SVM (Genus-CV) | 3,017 / 17,496 | 0.8402 | 0.8794 | 0.8593 | 0.8014 | 0.8715 | 0.05 s |
| Diagnostic Enzymes | BacDive-AI RF (Genus-CV) | 3,017 / 17,496 | 0.8391 | 0.8527 | 0.8458 | 0.7911 | 0.8593 | 0.05 s |
| Diagnostic Enzymes | DiASPora GBDT (Genus-CV) | 3,017 / 17,496 | 0.8602 | 0.8499 | 0.8550 | 0.8108 | 0.8699 | 0.05 s |
| Diagnostic Enzymes | **GSPA + DiASPora Hybrid** | 3,017 / 17,496 | **0.8652** | **0.8776** | **0.8713** | **0.8257** | **0.9097** | 1.2 s |

### Diagnostic enzyme panels and subunit-aware modeling

Values transcribed from `results/pillar3_bacdive_track1_summary.tsv` (rows 10, 11, 20 to 23):

- **10-EC Diagnostic Panel (50,075 assays across 6,313 strains):**
    - DeepProZyme baseline: Precision = 0.6442, Recall = 0.8105, \(F_1 = 0.7178\), Balanced Accuracy = 0.6792.
    - + METEOR v8 (negative-cost MILP): Precision = 0.6285, Recall = 0.8285, \(F_1 = 0.7148\), Balanced Accuracy = 0.6670. Uniform over-gap-filling increases false positives and degrades \(F_1\).
- **Heteromeric Urease (*EC 3.5.1.5*, 5,026 strains):**
    - DeepProZyme baseline: \(F_1 = 0.5624\).
    - + METEOR v8: \(F_1 = 0.5469\). Negative addition costs reward adding the reaction and boost isolated UreC alpha subunits in strains lacking UreA, UreB, or UreDFG, creating false-positive urease predictions.
    - GSPA role-aware DNF modeling enforces subunit completeness, suppressing isolated subunit false positives.
- **108-Genome Catabolic Assays (1,622 assays):**
    - DeepProZyme baseline: \(F_1 = 0.6711\), Balanced Accuracy = 0.5386.
    - + METEOR v8: \(F_1 = 0.6722\), Balanced Accuracy = 0.5101.

### Track 2 results (Partial Phenotype Evidence: 30% Revealed, Held-Out 70% Evaluated)

Values transcribed from `results/pillar3_bacdive_track2_summary.tsv`:

| Assay Panel | Method | Assays (Held-Out 70%) | Precision | Recall | \(F_1\) Score | Balanced Accuracy | ROC-AUC |
|---|---|---|---|---|---|---|---|
| Carbon Growth | CarveMe (Held-Out 70%) | 1,641 | 0.6927 | 0.3882 | 0.4975 | 0.6165 | 0.6165 |
| Carbon Growth | ModelSEED (Held-Out 70%) | 1,641 | 0.8728 | 0.5116 | 0.6451 | 0.7222 | 0.7222 |
| Carbon Growth | gapseq (Held-Out 70%) | 1,641 | 0.8059 | 0.7789 | 0.7922 | 0.8049 | 0.8049 |
| Carbon Growth | Traitar L1-SVM (30% \(\to\) 70%) | 1,641 | 0.8391 | 0.8380 | 0.8386 | 0.8466 | 0.9421 |
| Carbon Growth | BacDive-AI RF (30% \(\to\) 70%) | 1,641 | 0.8625 | 0.8625 | 0.8625 | 0.8692 | 0.9483 |
| Carbon Growth | DiASPora GBDT (30% \(\to\) 70%) | 1,641 | 0.8727 | 0.8985 | 0.8854 | 0.8901 | 0.9575 |
| Carbon Growth | **GSPA Track 2 (30% Bayesian Anchor)** | 1,641 | 0.8362 | **0.9447** | **0.8871** | **0.8889** | **0.9534** |
| Carbon Growth | **GSPA Track 2 + DiASPora Hybrid** | 1,641 | **0.9048** | **0.9409** | **0.9225** | **0.9258** | **0.9788** |
| Diagnostic Enzymes | CarveMe (Held-Out 70%) | 12,316 | 0.7791 | 0.3009 | 0.4341 | 0.5802 | 0.5802 |
| Diagnostic Enzymes | ModelSEED (Held-Out 70%) | 12,316 | 0.7868 | 0.4786 | 0.5952 | 0.6325 | 0.6325 |
| Diagnostic Enzymes | gapseq (Held-Out 70%) | 12,316 | 0.8376 | 0.7188 | 0.7737 | 0.7446 | 0.7446 |
| Diagnostic Enzymes | Traitar L1-SVM (30% \(\to\) 70%) | 12,316 | 0.8301 | 0.8854 | 0.8569 | 0.7934 | 0.8728 |
| Diagnostic Enzymes | BacDive-AI RF (30% \(\to\) 70%) | 12,316 | 0.8411 | 0.8575 | 0.8492 | 0.7953 | 0.8623 |
| Diagnostic Enzymes | DiASPora GBDT (30% \(\to\) 70%) | 12,316 | 0.8541 | 0.8473 | 0.8507 | 0.8045 | 0.8691 |
| Diagnostic Enzymes | **GSPA Track 2 (30% Bayesian Anchor)** | 12,316 | 0.8067 | **0.9084** | **0.8545** | 0.7749 | **0.8738** |
| Diagnostic Enzymes | **GSPA Track 2 + DiASPora Hybrid** | 12,316 | **0.8638** | **0.8819** | **0.8727** | **0.8264** | **0.9080** |

## Main findings

1. **Subunit-aware modeling prevents false-positive enzyme calls:** On 50,075 diagnostic enzyme assays across 6,313 strains, uniform MILP gap-filling (METEOR v8) lowers overall \(F_1\) from 0.7178 to 0.7148. On heteromeric urease (*EC 3.5.1.5*, 5,026 strains), boosting an isolated UreC alpha subunit without its obligate partners (UreA, UreB, and UreDFG) drops \(F_1\) from 0.5624 to 0.5469. Enforcing subunit completeness via disjunctive normal form rules prevents isolated subunit over-prediction.
2. **Continuous pathway coupling achieves superior ranking:** In blind Track 1 prediction (0% phenotype labels), pooling calibrated calls from gapseq, ModelSEED, and CarveMe with continuous pathway coupling matches gapseq's binary \(F_1\) (0.7952 on carbon growth, 0.7753 on enzymes) while raising ROC-AUC to 0.8567 on carbon growth and 0.7520 on enzymes, completing in 1.1 seconds per genome after base reconstruction.
3. **Out-of-genus generalization:** Under 5-fold genus-blocked cross-validation where test strains belong to genera unseen during training, the neuro-symbolic hybrid combining GSPA pathway posteriors with DiASPora GBDT achieves \(F_1 = 0.8710\) on carbon growth (ROC-AUC 0.9567) and \(F_1 = 0.8713\) on diagnostic enzymes (ROC-AUC 0.9097), outperforming all individual reconstructions and supervised baselines.
4. **Propagating partial phenotype evidence (Track 2):** Revealing 30% of a strain's experimental assays allows GSPA's unsupervised Bayesian update to propagate observations across coupled metabolic pathways, raising held-out carbon growth \(F_1\) from 0.7922 to 0.8871 (0.9225 in the DiASPora hybrid; ROC-AUC 0.9788). On diagnostic enzymes, the Track 2 hybrid achieves \(F_1 = 0.8727\) (ROC-AUC 0.9080).
