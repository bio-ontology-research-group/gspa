# Secondary metabolism and dark-matter MAGs

## Question

Does reaction-local genomic context (RLGC) and operon synteny extend beyond primary core metabolism to secondary metabolism (biosynthetic gene clusters) and uncultivated metagenome-assembled genomes (MAGs) from dark-matter phyla and desert extremophiles?

Paper: Section 2.5, Supplementary Table S1, and Supplementary Table S12, S16 (`paper/supplement.tex`).

## Datasets

### 106 MIBiG / antiSMASH biosynthetic gene clusters

106 curated biosynthetic gene clusters (BGCs) across 106 bacterial genomes identified through antiSMASH and MIBiG, spanning five major structural classes:

- **Non-ribosomal peptide synthetases (NRPS):** 38 BGCs, 114 core biosynthetic loci, 266 tailoring and transport loci.
- **Type I / II polyketide synthases (PKS):** 29 BGCs, 87 core biosynthetic loci, 203 tailoring and transport loci.
- **Terpenes and carotenoids:** 18 BGCs, 36 core biosynthetic loci, 90 tailoring and transport loci.
- **Siderophores (NRPS-independent / aerobactin-like):** 12 BGCs, 36 core biosynthetic loci, 72 tailoring and transport loci.
- **RiPPs / bacteriocins / lantipeptides:** 9 BGCs, 18 core biosynthetic loci, 54 tailoring and transport loci.
- **Total:** 106 BGCs, 291 core loci, 685 tailoring and transport loci. Source data: `results/pillar2_3_antismash_106bgc_specificity.tsv`.

### 97 uncultivated gut MAGs

97 high-quality metagenome-assembled genomes (MAGs) recovered from uncultured human and environmental gut microbiomes, spanning 16 bacterial and archaeal phyla including candidate phyla radiation (CPR) lineages:

- **Top 8 phyla:** Pseudomonadota (19 MAGs), Actinomycetota (12), Bacteroidota (11), Bacillota (9), Cyanobacteriota (8), Chloroflexota (7), Halobacteriota (7), Thermoproteota (6).
- **Other 8 phyla:** Patescibacteria / CPR (5), Verrucomicrobiota (3), Planctomycetota (2), Acidobacteriota (2), Spirochaetota (2), Campylobacterota (2), Asgardarchaeota (1), Thermodesulfobacteriota (1). Source data: `results/pillar3_dark_matter_97mags.tsv`.

### Empty Quarter (*Rub' al Khali*) desert extremophile vignette

Five bacterial isolates recovered from hyper-arid sands in Saudi Arabia with complete long-read assemblies:

- `EQ_01_Salinicola_sp` (*Salinicola* sp., Halomonadaceae, 3,842 CDS)
- `EQ_02_Halomonas_sp` (*Halomonas* sp., Halomonadaceae, 4,118 CDS)
- `EQ_03_Pseudomonas_sp` (*Pseudomonas* sp., Pseudomonadaceae, 5,490 CDS)
- `EQ_04_Bacillus_endophyticus_EQ` (*Bacillus* sp., Bacillaceae, 4,625 CDS)
- `EQ_05_Streptomyces_arid_EQ` (*Streptomyces* sp., Streptomycetaceae, 7,184 CDS)

Mean proteome size: 5,051.8 CDS. Source data: `results/pillar6_empty_quarter_vignette.tsv`.

## Protocol

- **BGC evaluation:** Discriminates core multi-modular assembly enzymes (NRPS/PKS condensation, adenylation, and ketosynthase domains) from auxiliary tailoring enzymes (glycosyltransferases, methyltransferases, tailoring oxidoreductases) and transporters within BGC boundaries versus non-BGC background genes. Evaluated metrics: core recall, tailoring recall, non-BGC false promotion rate, non-BGC specificity, and overall \(F_1\).
- **Metabolic dark matter in 97 MAGs:** Measures recovery of missing pathway reactions in uncultivated lineages lacking close cultured relatives. Evaluates candidate gap fills at posterior threshold \(q \ge 0.30\) and high-confidence fills at \(q \ge 0.50\), and reports pathway completeness before and after refinement.
- **Extremophile vignette:** Compares baseline sequence homology against Stage 0 pooling and full GSPA (Stages 1 to 3). Evaluates GAEF completeness, pathway coherence, pFBA-repaired reactions, Stage 3 promoted loci, and manual curation of complete extremophile stress-response operons.

## Baselines

Unconstrained Stage 0 pooling and standalone sequence homology search.

## Reproduce

As documented in `results/PROVENANCE_MAP.md`, these tables are Level 3 transcribed archives of cluster runs on `unimatrix01`.

To run BGC and operon-context annotation on an annotated assembly:

```bash
gspa integrate \
    --faa isolate.faa \
    --gff isolate.gff \
    --taxid 562 \
    --enable-fba-refiner \
    --output-dir results/bgc_annotation/
```

To run dark-matter gap-filling across a directory of MAG assemblies:

```bash
java -jar build/libs/gspa-unified.jar integrate \
    --batch-manifest mags_manifest.tsv \
    --output-dir results/mags_out/
```

## Result tables

- `results/pillar2_3_antismash_106bgc_specificity.tsv`: columns `bgc_class`, `n_genomes`, `n_bgcs`, `n_core_biosynthetic_loci`, `n_tailoring_transport_loci`, `stage0_core_recall`, `gspa_stage1_2_core_recall`, `stage0_tailoring_recall`, `gspa_stage1_2_tailoring_recall`, `non_bgc_false_promotion_rate`, `non_bgc_specificity`, `stage0_overall_f1`, `gspa_overall_f1`. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar2_3_antismash_106bgc_specificity.tsv))
- `results/pillar3_dark_matter_97mags.tsv`: columns `phylum_or_group`, `domain`, `n_mags`, `plot_group`, `candidate_fills_q030`, `hq_fills_q050`, `mean_fills_per_mag_q030`, `mean_fills_per_mag_q050`, `pathway_completeness_before`, `pathway_completeness_after`. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar3_dark_matter_97mags.tsv))
- `results/pillar6_empty_quarter_vignette.tsv`: columns `isolate_id`, `taxonomy`, `n_cds`, `baseline_annotated_pct`, `gspa_annotated_pct`, `gaef_completeness_raw_homology`, `gaef_completeness_stage0`, `gaef_completeness_gspa`, `gaef_pathway_coherence_base`, `gaef_pathway_coherence_gspa`, `fba_repaired_core_rxns`, `stage3_promoted_loci`, `completed_extremophile_operons`. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar6_empty_quarter_vignette.tsv))

### Biosynthetic gene cluster results (106 BGCs)

Values transcribed directly from `results/pillar2_3_antismash_106bgc_specificity.tsv`:

| BGC Class | Genomes | BGCs | Core Recall (Stage 0 \(\to\) GSPA) | Tailoring Recall (Stage 0 \(\to\) GSPA) | Non-BGC Specificity | Overall \(F_1\) (Stage 0 \(\to\) GSPA) |
|---|---|---|---|---|---|---|
| NRPS | 106 | 38 | 0.7807 \(\to\) **0.9211** | 0.5414 \(\to\) **0.8383** | **0.9958** | 0.7420 \(\to\) **0.8962** |
| Type I / II PKS | 106 | 29 | 0.7586 \(\to\) **0.9080** | 0.5222 \(\to\) **0.8227** | **0.9962** | 0.7227 \(\to\) **0.8865** |
| Terpene / Carotenoid | 106 | 18 | 0.8056 \(\to\) **0.9444** | 0.5667 \(\to\) **0.8556** | **0.9969** | 0.7619 \(\to\) **0.9136** |
| Siderophore | 106 | 12 | 0.8333 \(\to\) **0.9722** | 0.5972 \(\to\) **0.8889** | **0.9972** | 0.7935 \(\to\) **0.9384** |
| RiPP / Bacteriocin | 106 | 9 | 0.7222 \(\to\) **0.8889** | 0.4815 \(\to\) **0.7963** | **0.9955** | 0.6842 \(\to\) **0.8613** |
| **All 106 BGCs** | **106** | **106** | **0.7801 \(\to\) 0.9244** | **0.5387 \(\to\) 0.8380** | **0.9962** | **0.7395 \(\to\) 0.8959** |

### Microbial dark matter discovery in 97 uncultured MAGs

Values transcribed from `results/pillar3_dark_matter_97mags.tsv`:

| Phylum / Group | Domain | MAGs | Fills (\(q \ge 0.30\)) | HQ Fills (\(q \ge 0.50\)) | Mean Fills / MAG (\(q \ge 0.50\)) | Pathway Comp. Before | Pathway Comp. After |
|---|---|---|---|---|---|---|---|
| Pseudomonadota | Bacteria | 19 | 502 | 137 | 7.2 | 0.542 | **0.684** |
| Actinomycetota | Bacteria | 12 | 298 | 82 | 6.8 | 0.518 | **0.661** |
| Bacteroidota | Bacteria | 11 | 237 | 61 | 5.5 | 0.491 | **0.642** |
| Bacillota | Bacteria | 9 | 208 | 58 | 6.4 | 0.528 | **0.669** |
| Cyanobacteriota | Bacteria | 8 | 158 | 39 | 4.9 | 0.534 | **0.675** |
| Chloroflexota | Bacteria | 7 | 143 | 36 | 5.1 | 0.462 | **0.618** |
| Halobacteriota | Archaea | 7 | 132 | 32 | 4.6 | 0.489 | **0.631** |
| Thermoproteota | Archaea | 6 | 105 | 25 | 4.2 | 0.495 | **0.638** |
| Patescibacteria (CPR) | Bacteria | 5 | 104 | 26 | 5.2 | 0.384 | **0.529** |
| Verrucomicrobiota | Bacteria | 3 | 68 | 17 | 5.7 | 0.502 | **0.648** |
| Planctomycetota | Bacteria | 2 | 45 | 11 | 5.5 | 0.478 | **0.624** |
| Acidobacteriota | Bacteria | 2 | 44 | 11 | 5.5 | 0.486 | **0.632** |
| Spirochaetota | Bacteria | 2 | 43 | 11 | 5.5 | 0.512 | **0.654** |
| Campylobacterota | Bacteria | 2 | 45 | 11 | 5.5 | 0.548 | **0.682** |
| Asgardarchaeota | Archaea | 1 | 21 | 5 | 5.0 | 0.441 | **0.592** |
| Thermodesulfobacteriota | Bacteria | 1 | 22 | 6 | 6.0 | 0.498 | **0.644** |
| **All 97 MAGs** | **All** | **97** | **2,175** | **568** | **5.9** | **0.511** | **0.656** |

### Empty Quarter desert extremophile isolates

Values transcribed from `results/pillar6_empty_quarter_vignette.tsv`:

| Isolate ID | Taxonomy | CDS | Annotated (Base \(\to\) GSPA) | GAEF Comp. (Base \(\to\) GSPA) | Pathway Coh. (Base \(\to\) GSPA) | Repaired Core Rxns | Promoted Loci | Key Completed Operons |
|---|---|---|---|---|---|---|---|---|
| `EQ_01_Salinicola_sp` | *Salinicola* sp. | 3,842 | 74.2% \(\to\) **83.6%** | 0.881 \(\to\) **0.968** | 0.482 \(\to\) **0.614** | 6 | 5 | *ectABC-ask_ectD* (ectoine); *otsAB* (trehalose); *betIBA* (glycine betaine) |
| `EQ_02_Halomonas_sp` | *Halomonas* sp. | 4,118 | 75.8% \(\to\) **84.9%** | 0.892 \(\to\) **0.974** | 0.495 \(\to\) **0.628** | 5 | 6 | *ectABC* (ectoine); *crtEBIY* (carotenoids); *uvrABC-recA* (UV repair) |
| `EQ_03_Pseudomonas_sp` | *Pseudomonas* sp. | 5,490 | 78.4% \(\to\) **86.7%** | 0.914 \(\to\) **0.982** | 0.519 \(\to\) **0.641** | 4 | 4 | *treYZ-treS* (trehalose); *pvd* (pyoverdine siderophore); *katE-sodB* |
| `EQ_04_Bacillus_endophyticus` | *Bacillus* sp. | 4,625 | 76.9% \(\to\) **85.8%** | 0.901 \(\to\) **0.976** | 0.504 \(\to\) **0.632** | 5 | 7 | *opuAA-AB-AC* (osmoprotectant ABC); *spo0A-sigFGEK* (endospore); *dhbACEBF* |
| `EQ_05_Streptomyces_arid` | *Streptomyces* sp. | 7,184 | 68.5% \(\to\) **79.4%** | 0.854 \(\to\) **0.959** | 0.446 \(\to\) **0.598** | 8 | 9 | *ectABCD* (hydroxyectoine); *desABCD* (desferrioxamine); terpene BGCs |
| **Mean (5 Isolates)** | **All 5** | **5,051.8** | **74.76% \(\to\) 84.08%** | **0.8884 \(\to\) 0.9718** | **0.4892 \(\to\) 0.6226** | **5.6** | **6.2** | Complete osmoprotection, UV repair, and siderophore operons |

## Main findings

1. **Secondary metabolism specificity:** In 106 biosynthetic gene clusters, sequence-only Stage 0 identifies core assembly enzymes (78.01% recall) but misses auxiliary tailoring enzymes and transporters (53.87% recall). GSPA Stage 1+2 increases tailoring recall to 83.80% (+29.93 percentage points) and core recall to 92.44% while maintaining 99.62% non-BGC specificity (false promotion rate 0.0038), lifting overall BGC \(F_1\) from 0.7395 to 0.8959.
2. **Dark-matter gap filling across 16 phyla:** In 97 uncultivated gut MAGs, GSPA identifies 2,175 candidate gap fills at posterior \(q \ge 0.30\) (mean 22.4 per MAG) and 568 high-confidence fills at \(q \ge 0.50\) (mean 5.9 per MAG), increasing mean pathway completeness from 0.511 to 0.656 (+14.5 percentage points). In deep-branching CPR lineages (Patescibacteria), pathway completeness rises from 0.384 to 0.529.
3. **Completing extremophile operons:** In five hyper-arid desert isolates from the Empty Quarter, whole-genome inference increases annotated CDS from 74.76% to 84.08% and pathway coherence from 0.4892 to 0.6226. pFBA repairs a mean of 5.6 core reactions per isolate, and Stage 3 agent reasoning resolves 6.2 difficult loci per genome, completing interrupted operons for ectoine biosynthesis (*ectABC-ask_ectD*), trehalose synthesis (*otsAB*, *treYZ*), osmoprotectant transport (*opuAA-AC*), and siderophore synthesis (*desABCD*, *pvd*).
