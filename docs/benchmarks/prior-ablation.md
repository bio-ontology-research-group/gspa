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
| HMMER3 Pfam-A Only | 0.6412 (0.1845) | 0.6120 (0.2140) | 0.8120 | 0.3410 | 0.5840 | 18 | 0.6842 |
| NCBI PGAP (Stand-Alone) | 0.7285 (0.2110) | 0.6940 (0.2490) | 0.9150 | 0.4420 | 0.7810 | 14 | 0.8015 |
| DIAMOND-BLASTp Only | 0.7510 (0.2195) | 0.7180 (0.2580) | 0.8840 | 0.4380 | 0.7620 | 22 | 0.7890 |
| metagenomic-deepFRI (GCN) | 0.5940 (0.1780) | 0.6010 (0.2090) | 0.7950 | 0.3280 | 0.5420 | 31 | 0.6580 |
| InterProScan v5.68 Only | 0.7642 (0.2240) | 0.7315 (0.2645) | 0.9180 | 0.4610 | 0.8120 | 19 | 0.8194 |
| ESM-2 $k$-NN Transfer | 0.7880 (0.2310) | 0.7520 (0.2760) | 0.9260 | 0.4790 | 0.8340 | 29 | 0.8312 |
| Stage 0 Noisy-OR Base (`A1`) | 0.8041 (0.2361) | 0.7672 (0.2836) | 0.9400 | 0.5041 | 0.8636 | 25 | 0.8559 |
| + Essentiality Prior (`A2`, \(B_1\)) | 0.8074 (0.2389) | 0.7712 (0.2874) | **0.9826** | 0.5041 | 0.8636 | 25 | 0.8686 |
| + Pathway/Module Coherence (`A3`, \(B_2\)) | 0.8148 (0.2442) | 0.7796 (0.2948) | 0.9400 | 0.5616 | 0.9222 | 25 | 0.8812 |
| + Metabolic Gap Filling (`A4`, \(B_5\)) | 0.8089 (0.2405) | 0.7734 (0.2892) | 0.9400 | 0.5100 | 0.8636 | 25 | 0.8573 |
| + Operon Genomic Context (`A5`, \(B_6\)) | 0.8126 (0.2428) | 0.7768 (0.2925) | 0.9400 | 0.5336 | 0.8636 | 25 | 0.8616 |
| + Consistency / Anatomical (`A6`, \(B_3, B_4\)) | 0.8115 (0.2419) | 0.7754 (0.2914) | 0.9400 | 0.5038 | 0.8636 | 9 | 0.8635 |
| All Stage 1 Priors (`A7`, \(B_1\) to \(B_6\)) | 0.8240 (0.2512) | 0.7865 (0.3018) | **0.9826** | 0.5621 | 0.9222 | 9 | 0.9005 |
| Stage 1 + Stage 2A Removal (`A8`) | 0.8221 (0.2496) | 0.7848 (0.2995) | 0.9782 | **0.5622** | **0.9722** | **0** | 0.9090 |
| **Stage 1 + Stage 2A Relaxation (`A9`)** | **0.8276** (**0.2548**) | **0.7897** (**0.3054**) | **0.9826** | **0.5622** | **0.9722** | **0** | **0.9106** |
| *Out-of-Sample 8-Genome Mean (Base \(\to\) Full)* | *0.8143 \(\to\) **0.8382*** | *0.7764 \(\to\) **0.7994*** | *0.9325 \(\to\) **0.9805*** | *0.4904 \(\to\) **0.5514*** | *0.8743 \(\to\) **0.9844*** | *21 \(\to\) **0*** | *0.8531 \(\to\) **0.9096*** |

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

1. **Outperforming standalone genome annotation pipelines:** Multi-channel Stage 0 pooling (`A1`, CAFA \(F_{\max} = 0.8041\)) outperforms standalone HMMER3 Pfam-A (`0.6412`), metagenomic-deepFRI (`0.5940`), NCBI PGAP (`0.7285`), DIAMOND-BLASTp (`0.7510`), InterProScan v5.68 (`0.7642`), and ESM-2 $k$-NN (`0.7880`).
2. **Synergistic accuracy and structural gains from Stage 1 and Stage 2A:** Whole-genome inference with Stage 2A semantic ancestor relaxation (`A9`) raises CAFA \(F_{\max}\) from `0.8041` to `0.8276` (`0.2361` to `0.2548` on experimental-only annotations) and micro \(F_{\max}\) from `0.7672` to `0.7897` (`0.2836` to `0.3054` on experimental-only annotations; out-of-sample 8-genome CAFA \(F_{\max}\) `0.8143` to `0.8382`).
3. **Selective prior impacts:** Each prior selectively improves its corresponding genome property:
    - Essentiality prior (\(B_1\)) raises domain-core completeness \(C_{\text{ess}}\) from `0.9400` to `0.9826` and CAFA \(F_{\max}\) to `0.8074`.
    - Pathway and complex coherence prior (\(B_2\)) increases pathway coherence from `0.5041` to `0.5616`, complex completeness from `0.8636` to `0.9222`, and CAFA \(F_{\max}\) to `0.8148`.
    - Metabolic gap-filling (\(B_5\)) and operon context (\(B_6\)) priors increase pathway coherence to `0.5100` and `0.5336` and CAFA \(F_{\max}\) to `0.8089` and `0.8126`, respectively.
    - Soft taxonomic and anatomical consistency priors (\(B_3, B_4\)) suppress 16 of the 25 endogenous taxonomic violations in Stage 1 (`25` to `9`), while Stage 2A deterministic enforcement resolves the remaining 9 high-confidence violations (`9` to `0`) and prunes isolated obligate-complex singletons to reach `0.9722` complex coherence.
4. **Semantic ancestor relaxation outperforms hard deletion:** Within Stage 2A, relaxing clade-incompatible terms to their nearest valid GO ancestors (`A9`) outperforms hard term removal (`A8`) across CAFA \(F_{\max}\) (`0.8276` versus `0.8221`), micro \(F_{\max}\) (`0.7897` versus `0.7848`), essential completeness (`0.9826` versus `0.9782`), and composite audit score (`0.9106` versus `0.9090`).

