# Reference panels and production pipelines

## Question

Does calibrating and pooling complementary per-protein predictors in Stage 0 establish a stronger protein-level baseline than individual predictors, production annotation pipelines, and genomic language models when evaluated under a uniform ancestor-propagated protocol?

Paper: Section 2.1, Figure 2d, Supplementary Table S2 (`paper/supplement.tex`).

## Datasets

Evaluations span reference proteomes with high-confidence experimental and curated annotations:

- **13-genome panel:** 11 prokaryotes and 2 eukaryotes (`Panel13_All13`), and the 11 prokaryotic genomes alone (`Panel13_Prok11`) to allow direct comparison with NCBI PGAP, which is restricted to prokaryotic genomes.
- **21-genome panel:** 21 reference proteomes (`Panel21_All21`) covering bacteria, archaea, and unicellular and multicellular eukaryotes, drawn from `benchmark/manifests/suite2_proteomes25_manifest.tsv`.
- **9-genome subset:** 9 reference bacterial and archaeal genomes (`Panel9_Sub9`) evaluated against genomic language models that require chromosomal contiguity.

All Stage 1 and Stage 2 pipeline hyperparameters were calibrated on *Escherichia coli* K-12 MG1655 and held frozen across all other genomes (see [Hyperparameters](../reference/theta.md)).

## Protocol

- **Annotation propagation:** Gene Ontology (GO) terms across Molecular Function (MF), Biological Process (BP), and Cellular Component (CC) are propagated to the ontology root according to the true-path rule prior to threshold sweeping over \(\theta \in [0.01, 0.99]\).
- **Enzyme Commission (EC) evaluation:** 4-digit EC numbers are evaluated identically across all eligible enzymatic loci.
- **Metrics:**
    - **Micro \(F_{\max}\):** Global pair-level \(F_{\max}\) pooled across all proteins and terms.
    - **CAFA \(F_{\max}\):** Protein-centric \(F_{\max}\) averaged over annotated proteins.
    - **Semantic distance (\(S_{\min}\)):** Minimum semantic distance between prediction and truth; lower values indicate higher semantic precision.
    - **Efficiency:** Per-genome wall-clock runtime in seconds and peak resident memory (RAM) in megabytes.

## Baselines

- **Production annotation pipelines:** NCBI PGAP (Li et al. 2021), InterProScan 5.72 (Jones et al. 2014), and metagenomic-deepFRI (Bezshapkin et al. 2026).
- **Genomic language models:** gLM and gLM2 (Hwang et al. 2024; Cornman et al. 2024).
- **Individual Stage 0 predictor channels:** DIAMOND-BLASTp sequence alignment, ProtT5 \(k\)-NN transfer (Elnaggar et al. 2021), ESM-2 \(k\)-NN transfer (Lin et al. 2023), ProteInfer (Sanderson et al. 2023), ESM2-DeepGOPlus (2024), FoldSeek + ProstT5 3Di structural transfer (2024), and CLEAN contrastive EC prediction (Yu et al. 2023).

## Reproduce

The frozen metrics are stored in `results/pillar1_suite2_panel13_21_sota.tsv`. As noted in [How to reproduce](index.md#how-to-reproduce), this table is a Level 3 transcribed archive from cluster runs because running external baselines such as NCBI PGAP requires cluster-side reference databases.

To run GSPA's evaluation engine on predictions and truth files:

```bash
gspa evaluate \
    --predictions predictions.tsv \
    --truth ground_truth.tsv \
    --ontology go.obo \
    --propagate
```

For batch execution across the panel proteomes:

```bash
java -jar build/libs/gspa-unified.jar integrate \
    --batch-manifest benchmark/manifests/suite2_proteomes25_manifest.tsv \
    --output-dir results/panel_evaluation/
```

## Result tables

- `results/pillar1_suite2_panel13_21_sota.tsv`: columns `panel`, `method`, `n_genomes`, `go_fmax_micro_all`, `go_fmax_cafa_all`, `go_fmax_strict_holdout`, `go_fmax_exp_prop`, `go_fmax_strict_leaf`, `go_smin_all`, `go_smin_exp_prop`, `ec_fmax_micro`, `ec_fmax_cafa`, `runtime_sec_per_genome`, `peak_ram_mb`. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar1_suite2_panel13_21_sota.tsv))

### Empirical results

Values transcribed directly from `results/pillar1_suite2_panel13_21_sota.tsv`:

| Panel | Method | Genomes | GO Micro \(F_{\max}\) | GO CAFA \(F_{\max}\) | GO \(S_{\min}\) | EC Micro \(F_{\max}\) | EC CAFA \(F_{\max}\) | Runtime / Genome | Peak RAM |
|---|---|---|---|---|---|---|---|---|---|
| Panel13_Prok11 | metagenomic-deepFRI | 11 | 0.6010 | 0.5940 | 95.20 | n/a | n/a | 38.0 s | 4,800 MB |
| Panel13_Prok11 | InterProScan 5.72 | 11 | 0.7410 | 0.6780 | 58.40 | n/a | n/a | 680.0 s | 8,200 MB |
| Panel13_Prok11 | NCBI PGAP | 11 | 0.7820 | 0.7680 | 54.20 | n/a | n/a | 1,800.0 s | 16,000 MB |
| Panel13_Prok11 | **GSPA Stage 1+2 (3-Channel Full)** | 11 | **0.8421** | **0.8605** | **21.40** | **0.8830** | **0.8870** | **5.4 s** | **4,812 MB** |
| Panel13_All13 | metagenomic-deepFRI | 13 | 0.5850 | 0.5790 | 98.40 | n/a | n/a | 42.0 s | 4,800 MB |
| Panel13_All13 | InterProScan 5.72 | 13 | 0.7220 | 0.6530 | 62.10 | n/a | n/a | 720.0 s | 8,200 MB |
| Panel13_All13 | **GSPA Stage 1+2 (3-Channel Full)** | 13 | **0.8180** | **0.8060** | **24.80** | **0.8830** | **0.8870** | **5.8 s** | **4,812 MB** |
| Panel21_All21 | FoldSeek + ProstT5 3Di | 21 | 0.6840 | 0.6310 | 83.24 | n/a | n/a | 10,800.0 s | 14,000 MB |
| Panel21_All21 | InterProScan 5.72 | 21 | 0.7180 | 0.6490 | 95.93 | 0.3120 | n/a | 720.0 s | 8,200 MB |
| Panel21_All21 | ESM2-DeepGOPlus | 21 | 0.7420 | 0.6950 | 105.29 | n/a | n/a | 420.0 s | 6,400 MB |
| Panel21_All21 | ProteInfer | 21 | 0.7950 | 0.7540 | 38.12 | 0.3880 | n/a | 120.0 s | 5,800 MB |
| Panel21_All21 | DIAMOND-BLASTp | 21 | 0.8400 | 0.7890 | 93.81 | 0.8550 | n/a | 17.0 s | 2,100 MB |
| Panel21_All21 | ProtT5 \(k\)-NN | 21 | 0.8530 | 0.8010 | 39.40 | n/a | n/a | 28.0 s | 6,200 MB |
| Panel21_All21 | ESM-2 \(k\)-NN | 21 | 0.8580 | 0.8080 | 37.85 | n/a | n/a | 24.0 s | 5,800 MB |
| Panel21_All21 | CLEAN Contrastive EC | 21 | n/a | n/a | n/a | 0.8530 | 0.8590 | 300.0 s | 4,200 MB |
| Panel21_All21 | **GSPA Stage 1+2 (3-Channel Core)** | 21 | **0.8690** | **0.8220** | **36.14** | **0.8830** | **0.8870** | **6.1 s** | **4,812 MB** |
| Panel21_All21 | **GSPA Stage 1+2 (7-Channel All-Pooled)** | 21 | **0.8705** | **0.8407** | **35.20** | **0.8860** | **0.8900** | **6.8 s** | **4,812 MB** |
| Panel9_Sub9 | gLM / gLM2 Genomic LM | 9 | 0.7280 | 0.7610 | 37.40 | n/a | n/a | 95.0 s | 7,200 MB |

## Main findings

1. **Multi-predictor pooling outperforms standalone pipelines:** Across the 11 prokaryotic genomes, GSPA using three channels reaches micro \(F_{\max} = 0.8421\) and CAFA \(F_{\max} = 0.8605\), improving over NCBI PGAP (0.7820 micro, 0.7680 CAFA) and InterProScan 5.72 (0.7410 micro, 0.6780 CAFA).
2. **Complementary channels beat single predictors:** On the 21-genome panel, 7-channel GSPA achieves micro \(F_{\max} = 0.8705\) and CAFA \(F_{\max} = 0.8407\), exceeding the strongest single channel, ESM-2 \(k\)-NN (0.8580 micro, 0.8080 CAFA), while reducing semantic error \(S_{\min}\) from 37.85 to 35.20.
3. **EC prediction accuracy:** In enzyme annotation, GSPA achieves EC micro \(F_{\max} = 0.8860\) and CAFA \(F_{\max} = 0.8900\), surpassing CLEAN (0.8530 micro, 0.8590 CAFA) and DIAMOND-BLASTp (0.8550 micro).
4. **Efficiency:** Stages 1 and 2 complete in 5.4 to 6.8 seconds per genome with 4,812 MB peak heap memory, compared to 1,800 seconds and 16,000 MB for NCBI PGAP, and 10,800 seconds and 14,000 MB for FoldSeek + ProstT5 3Di.
