# Enzyme EC holdouts and homology masking

## Question

When enzyme predictors are combined in Stage 0 and refined with metabolic context, does EC prediction improve on enzymes characterised after the training cutoff? And when metabolic gap filling proposes changes to the top-ranked EC number of a protein, which proposed changes should be accepted, especially as close homologues are removed from the reference databases?

Paper: Supplementary Tables S5 and S6 (Suites 1C and 1D in Supplementary Table S1).

## Datasets

### Post-cutoff Swiss-Prot holdout (Holdout-125)

125 Swiss-Prot enzymes from 74 species whose EC annotations were added after the training cutoff of the predictors, with 293 ground-truth protein-EC pairs. Only complete four-level EC numbers are used. The sequences and ground truth are archived under `benchmark/manifests/reference_data/paperA_2026/benchmark/swissprot_holdout/` (`holdout227.fasta`, `swissprot_holdout_sample.tsv`); the 125 evaluated proteins are those of the larger archived sample that carry a complete EC number. The archived predictor outputs are `benchmark/manifests/reference_data/ibex_extracted/holdout227_*.pkl`.

### 95 bacterial genomes under donor masking

All 69,666 proteins of 95 bacterial RefSeq genomes, annotated by three enzyme predictors (CLEAN, DeepProZyme, EnzBERT) under four donor-masking tiers: unfiltered, \( \le 70\% \), \( \le 50\% \), and \( \le 30\% \) identity to any reference donor (see [Donor-homology masking](index.md#donor-homology-masking)). For CLEAN in the unfiltered tier, 73 genomes with 52,814 proteins are available. The archived per-genome evaluation files are `benchmark/manifests/reference_data/ibex_extracted/full_genome_eval_gbff_fmax_{clean,dpz,enzbert}_{vanilla,filt70,filt50,filt30}.pkl`.

## Protocol

- **Holdout-125.** Each predictor is scored alone and, for DeepProZyme, after single-predictor refinement with METEOR v8 or with GSPA (locus refinement and coherence). The three predictors are then pooled in Stage 0 (Noisy-OR) and refined for one to three iterations. Metrics: Top-1, Top-3, and Top-5 accuracy (the true EC among the \(k\) highest-scoring ECs), mean reciprocal rank, and micro \( F_{\max} \) over protein-EC pairs. EnzBERT, CLEAN, and DeepProZyme are also evaluated with masked lookup databases at the three identity tiers. The targets are isolated sequences without genome coordinates, so genome-scale priors act only through EC-level pathway coherence and the [Stage 2B](../concepts/stage2b-fba.md) machinery.
- **95 genomes.** For each predictor and masking tier, METEOR v8 proposes changes to the top-ranked EC of individual proteins through reaction-level gap filling. Three acceptance policies are compared against the unrefined predictor: accept every METEOR change; accept a change only if GSPA's single-predictor stoichiometric and pathway filter supports it; accept a change only if GSPA's three-channel consensus (agreement of the other two predictors) supports it. For each policy the table reports Top-1 accuracy, \( F_{\max} \), the number of corrected (\( b_{01} \)) and newly broken (\( b_{10} \)) proteins, their difference (net), and a two-sided exact McNemar test on the discordant pairs.

## Baselines

CLEAN, DeepProZyme, and EnzBERT as published (with their own homology-split lookup protocols and frozen weights), and METEOR v8 negative-cost gap filling applied to each.

## Reproduce

```bash
uv run --with numpy --with scipy --with pandas \
    python scripts/eval_scenario1_holdouts_real.py
```

The same script also writes the [prospective GO holdout](prospective-go-holdouts.md) tables. For the two enzyme tables it recomputes all metrics from the archived prediction matrices and evaluation files listed above (it imports `scripts/stoichiometric_locus_refiner.py` and reads the reference metabolic data in `benchmark/manifests/reference_data/psb_revision/code_snapshot/data`). The masked CLEAN and DeepProZyme rows of Holdout-125 are read from the archived `benchmark/manifests/reference_data/paperA_2026/results/eval_consolidated_full.tsv`. Rerunning the predictors themselves requires their released models and is not part of the repository. The script sets `ROOT` to the authors' checkout path and writes into `results/`; see [How to reproduce](index.md#how-to-reproduce).

## Result tables

- `results/scenario1_swissprot_ec_holdout125_real.tsv`: one row per method and masking tier; columns `n_test`, `top1_hits`, `top1_acc`, `top3_hits`, `top3_acc`, `top5_hits`, `top5_acc`, `mrr`, `fmax`. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/scenario1_swissprot_ec_holdout125_real.tsv))
- `results/scenario1_95genomes_ec_homology_real.tsv`: one row per predictor and tier (`vanilla`, `filt70`, `filt50`, `filt30`); for the base predictor, METEOR, GSPA single-predictor, and GSPA consensus policies the columns `*_top1`, `*_top1_acc`, `*_fmax`, `*_corr`, `*_reg`, `*_net`, `*_mcnemar_p`. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/scenario1_95genomes_ec_homology_real.tsv))

## Main finding

On Holdout-125, Stage 0 pooling of the three predictors gives the highest micro \( F_{\max} \) (`fmax` 0.4381, row `Noisy-OR (3-Method)`), and refinement trades a little Top-1 accuracy for more correct ECs in the top three and top five.

On the 95 genomes, accepting every METEOR change lowers Top-1 accuracy for EnzBERT in every tier and for DeepProZyme in the unfiltered and \( \le 70\% \) tiers, because reaction-level boosts promote non-cognate paralogs. Accepting only changes supported by the three-channel consensus gives a positive net gain for all three predictors in all four tiers, largest for CLEAN (`gspa_net` 155 in row `CLEAN (2023)`, `vanilla`). Agreement among complementary predictors identifies which metabolic updates to trust.
