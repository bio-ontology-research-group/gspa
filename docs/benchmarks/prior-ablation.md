# Genome-scale prior ablation and perturbation

## Question

What does each Stage 1 whole-genome prior contribute to organism-level structural coherence, does iterative Bayesian refinement preserve the protein-level accuracy baseline, and do taxonomic and anatomical constraints eliminate cross-clade annotation errors?

Paper: Section 2.2, Table 1, and Figure 2a to 2c (`paper/main.tex`).

## Datasets

### Reference genome ablation panel

Nine bacterial and archaeal reference genomes representing diverse phylogenetic clades and metabolic styles, totaling 81 genome-scale inference runs across nine configurations:

- *Escherichia coli* K-12 MG1655 (calibration genome)
- *Escherichia coli* O157:H7 EDL933 (enteric pathogen)
- *Bacillus subtilis* 168 (Gram-positive model)
- *Pseudomonas aeruginosa* PAO1 (opportunistic pathogen)
- *Mycobacterium tuberculosis* H37Rv (slow-growing actinobacterium)
- *Helicobacter pylori* 26695 (microaerophilic pathogen)
- *Mycoplasma genitalium* G37 (minimal genome bacterium)
- *Synechocystis* sp. PCC 6803 (photosynthetic cyanobacterium)
- *Methanocaldococcus jannaschii* DSM 2661 (methanogenic archaeon)

Hyperparameters were calibrated on *E. coli* K-12 and held fixed for the other eight genomes (reported as the out-of-sample 8-genome mean).

### Cross-clade perturbation cohort

A controlled perturbation experiment injecting 100 cross-clade false-positive claims into each of the nine genomes (900 injected claims in total). Injected terms represent biological processes structurally impossible in the host clades, such as animal organogenesis (*heart development*, *neural crest cell migration*) and plant-specific cell wall synthesis.

## Protocol

### Nine ablation configurations

1. `A1_noisy_or`: Unconstrained Stage 0 Noisy-OR baseline pooling DIAMOND, Pfam, and InterProScan.
2. `A2_essentiality`: Stage 0 plus the domain essentiality prior (\(B_1\)).
3. `A3_coherence`: Stage 0 plus the pathway and heteromeric complex coherence prior (\(B_2\)).
4. `A4_gap_filling`: Stage 0 plus the metabolic gap-filling prior (\(B_5\)).
5. `A5_genomic_context`: Stage 0 plus the co-directional operon context prior (\(B_6\)).
6. `A6_consistency_anatomical`: Stage 0 plus taxonomic and anatomical consistency priors (\(B_3, B_4\)).
7. `A7_stage1_all`: Full Stage 1 integration combining all six priors (\(B_1\) to \(B_6\)).
8. `A8_stage2_remove`: Full Stage 1 plus Stage 2A deterministic violation deletion.
9. `A9_stage2_relax`: Full Stage 1 plus Stage 2A semantic ancestor relaxation.

### Metrics

- **CAFA \(F_{\max}\) and Micro \(F_{\max}\):** Evaluated against all GO annotations and experimental-only ground truth.
- **Essentiality completeness (\(C_{\text{ess}}\)):** Fraction of core domain-essential functional roles present.
- **Pathway coherence:** Mean completeness of active metabolic pathways.
- **Complex coherence:** Fraction of heteromeric complexes with all obligate subunits identified.
- **Consistency:** \(C_{\text{cons}} = \max(0, 1 - 0.01 V)\) for \(V\) formal ontology violations.
- **Composite structural audit score (GAEF):** \(0.30 C_{\text{ess}} + 0.40 C_{\text{func}} + 0.30 C_{\text{cons}}\).
- **Perturbation suppression rate:** Percentage of surviving cross-clade false claims suppressed.

## Baselines

The baseline is the unconstrained Stage 0 Noisy-OR ensemble. Each configuration is evaluated against this baseline and compared across individual prior additions and Stage 2A enforcement variants.

## Reproduce

The complete ablation pipeline is orchestrated by `benchmark/unified/pillar1_2/run_all_pillar1_2.sh`.

```bash
cd benchmark/unified/pillar1_2
bash run_all_pillar1_2.sh
```

This workflow executes five modular scripts:

1. `01_run_ablation_grid.py`: Runs `gspa-unified.jar integrate` across all 9 genomes and 9 configurations.
2. `02_run_cafa_eval.py`: Computes ancestor-propagated CAFA and micro \(F_{\max}\).
3. `03_run_gaef_eval.py`: Evaluates structural invariants (\(C_{\text{ess}}\), coherence, consistency).
4. `04_run_perturbation.py`: Injects 900 cross-clade claims and measures suppression rates.
5. `05_compile_results.py`: Compiles the summary tables into `results/`.

External cluster dependencies are required for step 1. If running offline with precomputed outputs, steps 2 to 5 re-tabulate the metrics directly.

## Result tables

- `results/pillar1_2_ablation_summary.tsv`: 9 ablation configurations with 9-genome averages and out-of-sample 8-genome rows; columns `config`, `label`, `fmax_micro_exp`, `delta_micro_exp`, `fmax_cafa_exp`, `delta_cafa_exp`, `fmax_micro_all`, `fmax_cafa_all`, `coverage_exp`, `ic_recall_exp`, `completeness`, `pathway_coherence`, `complex_coherence`, `total_violations`, `composite_score`. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar1_2_ablation_summary.tsv))
- `results/pillar1_2_ablation_per_genome.tsv`: per-genome metrics across all 81 runs. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar1_2_ablation_per_genome.tsv))
- `results/pillar1_2_perturbation_summary.tsv`: injected false claims, surviving violations, and suppression rates per genome and configuration. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar1_2_perturbation_summary.tsv))

### Empirical results

Values transcribed from `results/pillar1_2_ablation_summary.tsv`:

| Configuration | CAFA \(F_{\max}\) (All / Exp) | Micro \(F_{\max}\) (All / Exp) | Essential \(C_{\text{ess}}\) | Pathway Coherence | Complex Coherence | Violations (\(V\)) | Audit Score |
|---|---|---|---|---|---|---|---|
| Stage 0 Noisy-OR Base | 0.8041 (0.2361) | 0.7672 (0.2836) | 0.9400 | 0.5041 | 0.8636 | 25 | 0.8559 |
| + Essentiality Prior (\(B_1\)) | 0.8040 (0.2361) | 0.7683 (0.2838) | **0.9826** | 0.5041 | 0.8636 | 25 | 0.8686 |
| + Pathway/Module Coherence (\(B_2\)) | 0.8028 (0.2357) | 0.7671 (0.2815) | 0.9400 | 0.5616 | 0.9222 | 25 | 0.8812 |
| + Metabolic Gap Filling (\(B_5\)) | 0.8039 (0.2361) | 0.7685 (0.2838) | 0.9400 | 0.5100 | 0.8636 | 25 | 0.8573 |
| + Operon Genomic Context (\(B_6\)) | 0.8037 (0.2356) | 0.7675 (0.2815) | 0.9400 | 0.5336 | 0.8636 | 25 | 0.8616 |
| + Consistency / Anatomical (\(B_3, B_4\)) | **0.8060** (0.2371) | 0.7711 (0.2851) | 0.9400 | 0.5038 | 0.8636 | 9 | 0.8635 |
| All Stage 1 Priors (\(B_1\) to \(B_6\)) | 0.8043 (0.2368) | 0.7703 (0.2828) | **0.9826** | 0.5621 | 0.9222 | 9 | 0.9005 |
| Stage 1 + Stage 2A Removal | 0.8053 (0.2378) | 0.7710 (0.2852) | 0.9782 | **0.5622** | **0.9722** | **0** | 0.9090 |
| **Stage 1 + Stage 2A Relaxation** | 0.8056 (**0.2382**) | **0.7733** (**0.2863**) | **0.9826** | **0.5622** | **0.9722** | **0** | **0.9106** |
| *Out-of-Sample 8-Genome Mean (Base \(\to\) Full)* | *0.8143 \(\to\) **0.8160*** | *0.7764 \(\to\) **0.7830*** | *0.9325 \(\to\) **0.9805*** | *0.4904 \(\to\) **0.5514*** | *0.8743 \(\to\) **0.9844*** | *21 \(\to\) **0*** | *0.8531 \(\to\) **0.9096*** |

### Cross-clade perturbation results

Values transcribed from `results/pillar1_2_perturbation_summary.tsv` across the 900 injected claims (100 per genome):

| Genome | Injected Claims | Surviving in Stage 0 | Surviving in Stage 2A | Stage 2A Suppression Rate |
|---|---|---|---|---|
| *H. pylori* 26695 | 100 | 10 | 0 | 100.0% |
| *M. genitalium* G37 | 100 | 14 | 0 | 100.0% |
| *M. jannaschii* DSM 2661 | 100 | 11 | 2 | 98.0% |
| *E. coli* O157:H7 EDL933 | 100 | 15 | 1 | 99.0% |
| *E. coli* K-12 MG1655 | 100 | 14 | 0 | 100.0% |
| *B. subtilis* 168 | 100 | 14 | 0 | 100.0% |
| *Synechocystis* sp. PCC 6803 | 100 | 15 | 0 | 100.0% |
| *M. tuberculosis* H37Rv | 100 | 12 | 0 | 100.0% |
| *P. aeruginosa* PAO1 | 100 | 20 | 0 | 100.0% |
| **Total / Overall** | **900** | **125** | **3** | **97.6% (99.7% of all 900)** |

## Main findings

1. **Accuracy improvement and preservation:** Whole-genome inference with Stage 2A semantic ancestor relaxation improves both protein-centric CAFA \(F_{\max}\) (0.8041 to 0.8056 on all annotations; 0.2361 to 0.2382 on experimental-only annotations) and pair-level micro \(F_{\max}\) (0.7672 to 0.7733 on all annotations; 0.2836 to 0.2863 on experimental-only annotations; out-of-sample 8-genome micro \(F_{\max}\) 0.7764 to 0.7830).
2. **Selective prior impacts:** Each prior selectively improves its corresponding genome property:
    - Essentiality prior (\(B_1\)) raises domain-core completeness \(C_{\text{ess}}\) from 0.9400 to 0.9826 (out-of-sample 0.9325 to 0.9805).
    - Pathway and complex coherence prior (\(B_2\)) increases pathway coherence from 0.5041 to 0.5616 and heteromeric complex completeness from 0.8636 to 0.9222.
    - Metabolic gap-filling (\(B_5\)) and operon context (\(B_6\)) priors increase pathway coherence to 0.5100 and 0.5336, respectively.
    - Soft taxonomic and anatomical consistency priors (\(B_3, B_4\)) suppress 16 of the 25 endogenous taxonomic violations in Stage 1 (`25` to `9`), while Stage 2A deterministic enforcement resolves the remaining 9 high-confidence violations (`9` to `0`) and prunes isolated obligate-complex singletons to reach 0.9722 complex coherence.
3. **Semantic ancestor relaxation outperforms hard deletion:** Within Stage 2A, relaxing clade-incompatible terms to their nearest valid GO ancestors (`A9_stage2_relax`) outperforms hard term removal (`A8_stage2_remove`) across CAFA \(F_{\max}\) (0.8056 versus 0.8053), micro \(F_{\max}\) (0.7733 versus 0.7710), essential completeness (0.9826 versus 0.9782), experimental coverage (0.8057 versus 0.8000), and composite audit score (0.9106 versus 0.9090).
4. **Structural audit improvement:** The composite structural audit score increases from 0.8559 to 0.9106 across all nine genomes (\(p = 0.0039\), two-sided Wilcoxon signed-rank test; out-of-sample 8-genome mean 0.8531 to 0.9096).
5. **Cross-clade error elimination:** Of 900 injected cross-clade false-positive claims, 125 survive unconstrained Stage 0 calibration (86.1% baseline suppression). Stage 2A deterministic enforcement eliminates 122 of these 125 surviving claims (97.6% post-calibration suppression; 99.7% total suppression), leaving zero violations in seven of the nine species.
