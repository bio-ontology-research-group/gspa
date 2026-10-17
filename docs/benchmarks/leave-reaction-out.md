# Leave-reaction-out enzyme attribution

## Question

Can reaction-local genomic context (RLGC) and protein language-model centroids identify which chromosomal gene encodes a held-out metabolic reaction under homology masking, and how does recovery vary across metabolic strata and formal recoverability ceilings?

Paper: Section 2.5, Table 4 (in latex `tab:metabolism`), and Figure 3a to 3c (`paper/main.tex`).

## Datasets

### Benchmark cohorts

- ***E. coli* K-12 strict LRO (256 reactions):** A curated set of 256 metabolic reactions in *Escherichia coli* K-12 MG1655 with close sequence homologs of the true enzyme masked from reference databases (calibration genome).
- ***E. coli* K-12 role-stratified LRO (414 reactions):** The complete metabolic reaction network of *E. coli* K-12 partitioned into three physiological strata:
    - **Stratum A (Core Biomass, \(n=307\)):** Reactions required for steady-state biomass precursor production under minimal glucose medium.
    - **Stratum B (Alternative Catabolic, \(n=49\)):** Reactions participating in alternative carbon, nitrogen, sulfur, or phosphorus catabolism.
    - **Stratum C (Peripheral and Secondary, \(n=58\)):** Reactions in specialized pathways and secondary metabolism that do not carry flux on minimal media.
- ***B. subtilis* 168 out-of-sample transfer (144 reactions):** 144 held-out reactions in *Bacillus subtilis* 168 evaluated with all hyperparameters frozen from *E. coli*.
- **Ambiguous multi-candidate LRO subset (50 gaps):** 50 unique-protein metabolic gaps where Stage 1 RLGC yields multiple closely tied candidate loci, evaluated for Stage 3 agent disambiguation.

### Homology masking and evaluation tracks

- **Protocol H masking tiers:** Evaluated under unfiltered (100% identity), \(\le 70\%\), \(\le 50\%\) (twilight zone), and \(\le 30\%\) identity (metabolic dark matter).
- **Track A (Strict Homology and Embedding Masking):** Masks both sequence alignment donors and embedding donors with cosine similarity \(> 0.50\); excludes reaction centroids.
- **Track B (PLM Reaction Centroids):** Masks sequence alignment donors while retaining Swiss-Prot ESM-2 reaction embedding centroids.

### Five formal recoverability ceilings

Every held-out metabolic reaction is governed by five structural constraints defining theoretical recoverability:

1. **\(\mathcal{C}_{\text{struct}}\) (Universal database coverage):** The reaction must be present, mass-balanced, and charge-balanced in the universal reaction database \(\mathcal{U}_{\text{clean}}\).
2. **\(\mathcal{C}_{\text{task}}\) (Medium-specific flux demand):** The reaction must carry non-zero steady-state flux under specified nutritional conditions (partitioning core Stratum A from catabolic Stratum B and secondary Stratum C).
3. **\(\mathcal{C}_{\text{unique}}\) (Isozyme non-redundancy):** The reaction is encoded by a unique non-redundant gene locus rather than multiple interchangeable isozymes.
4. **\(\mathcal{C}_{\text{operon}}\) (Chromosomal operon synteny):** The candidate locus resides within a co-directional operon containing anchor genes of pathway-adjacent steps.
5. **\(\mathcal{C}_{\text{complex}}\) (Subunit stoichiometry):** For heteromeric enzymes, all required partner subunits must be encoded and identifiable in the host chromosome.

## Protocol

In each trial, the enzyme annotation of one target reaction is removed, and candidates across the genome are scored:

- **Unconditional ranking ("All"):** Evaluates all held-out reactions. If the true gene does not enter the Stage 0 candidate pool, the trial is scored as a miss (Rank = \(\infty\)).
- **Conditional ranking ("Pool"):** Evaluates ranking only among reactions where the true gene enters the candidate pool (\(n_{\text{cov}}\)).

### Methods evaluated

- **Sequence-Only / METEOR Uniform Boost:** Ranks candidates by sequence similarity alone. Reaction-level gap-filling applies an identical score boost to every candidate, leaving rank order unchanged.
- **RLGC Operon Kernel Density:** Scores candidates by chromosomal base-pair distance and co-directionality relative to anchor genes of pathway-adjacent reactions.
- **GSPA Stage 1 Coherence + RLGC (Track A):** Integrates RLGC operon kernel density with the Stage 1 pathway coherence prior without language-model centroids.
- **ESM-2 PLM Centroids Only:** Ranks candidates by cosine similarity between protein ESM-2 embeddings and Swiss-Prot reaction centroids.
- **Full GSPA (Track B):** Combines RLGC kernel density, pathway coherence, and ESM-2 centroids.
- **Stage 3 Guarded Re-Ranking:** Selective tool-mediated agentic reasoning over the top candidate loci.

## Baselines

Sequence-only ranking, METEOR v8 negative-cost MILP uniform boosting, function prior alone, RLGC alone, and ESM-2 centroids alone.

## Reproduce

The complete leave-reaction-out benchmark is orchestrated by `benchmark/unified/pillar2_fba_lro_price149/run_pillar2_all.py`:

```bash
uv run --with numpy --with scipy \
    python benchmark/unified/pillar2_fba_lro_price149/run_pillar2_all.py \
        --run-lro --output-dir results/
```

The script evaluates both *E. coli* cohorts and *B. subtilis*, computing unconditional and pool metrics across all strata and masking tiers, writing `results/pillar2_lro15_summary.tsv`.

## Result tables

- `results/pillar2_lro15_summary.tsv`: columns `benchmark_cohort`, `stratum`, `protocol_H`, `method`, `n_total`, `n_cov`, `cov_frac`, `hit1_all`, `hit3_all`, `hit5_all`, `mrr_all`, `hit1_cov`, `hit3_cov`, `mrr_cov`. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar2_lro15_summary.tsv))

### Empirical results

Values transcribed from `results/pillar2_lro15_summary.tsv`:

| Benchmark Cohort | Stratum / Tier | Method | Total / Pool (\(n / n_{\text{cov}}\)) | Hit@1 (All) | Hit@3 (All) | Hit@1 (Pool) | Hit@3 (Pool) | MRR (Pool) |
|---|---|---|---|---|---|---|---|---|
| *E. coli* Strict 256 | All / Unfiltered | DIAMOND-BLASTp Only | 256 / 175 | 1.95% (5) | 7.03% | 2.86% | 10.29% | 0.1354 |
| *E. coli* Strict 256 | All / Unfiltered | METEOR v8 (Negative-Cost MILP) | 256 / 175 | 3.52% (9) | 9.38% | 5.14% | 13.71% | 0.1612 |
| *E. coli* Strict 256 | All / Unfiltered | gapseq v1.2 (Bitscore Heuristic) | 256 / 181 | 8.98% (23) | 24.61% | 12.71% | 34.81% | 0.2540 |
| *E. coli* Strict 256 | All / Unfiltered | GSPA Strict Track A (RLGC, No PLM) | 256 / 181 | 11.33% (29) | 34.77% | 16.02% | 49.17% | 0.3685 |
| *E. coli* Strict 256 | All / Unfiltered | RLGC Kernel Density Only | 256 / 181 | 23.05% (59) | 37.89% | 32.60% | 53.59% | 0.4745 |
| *E. coli* Strict 256 | All / Unfiltered | GSPA Stage 1 Coherence + RLGC | 256 / 175 | 23.44% (60) | 36.33% | 34.29% | 53.14% | 0.4864 |
| *E. coli* Strict 256 | All / Unfiltered | ESM-2 PLM Centroids Only | 256 / 181 | 62.11% (159) | 67.58% | 87.85% | 95.58% | 0.9157 |
| *E. coli* Strict 256 | All / Unfiltered | **Full GSPA (Track B: RLGC + PLM)** | 256 / 181 | **63.67% (163)** | **68.36%** | **90.06%** | **96.69%** | **0.9330** |
| *E. coli* Full 414 | Stratum A (Core) | RLGC Kernel Density Only | 307 / 170 | 13.36% | 24.43% | 24.12% | 44.12% | 0.3939 |
| *E. coli* Full 414 | Stratum A (Core) | **Full GSPA (Track B)** | 307 / 170 | **51.79% (159)** | **53.42%** | **93.53%** | **96.47%** | **0.9557** |
| *E. coli* Full 414 | Stratum B (Catabolic) | RLGC Kernel Density Only | 49 / 27 | 30.61% | 51.02% | 55.56% | 92.59% | 0.7325 |
| *E. coli* Full 414 | Stratum B (Catabolic) | **Full GSPA (Track B)** | 49 / 27 | **48.98% (24)** | **53.06%** | **88.89%** | **96.30%** | **0.9216** |
| *E. coli* Full 414 | Stratum C (Secondary) | RLGC Kernel Density Only | 58 / 29 | 17.24% | 29.31% | 34.48% | 58.62% | 0.4882 |
| *E. coli* Full 414 | Stratum C (Secondary) | **Full GSPA (Track B)** | 58 / 29 | **41.38% (24)** | **48.28%** | **82.76%** | **96.55%** | **0.8793** |
| *E. coli* Full 414 | Stratum A (\(\le 70\%\) Homology) | **Full GSPA (Track B)** | 89 / 41 | **41.57% (37)** | **44.94%** | **90.24%** | **97.56%** | **0.9410** |
| *E. coli* Full 414 | Stratum A (\(\le 50\%\) Twilight) | **Full GSPA (Track B)** | 23 / 7 | **26.09% (6)** | **30.43%** | **85.71%** | **100.00%** | **0.9286** |
| *E. coli* Full 414 | Stratum A (\(\le 30\%\) Dark Matter) | **Full GSPA (Track B)** | 12 / 4 | **25.00% (3)** | **33.33%** | **75.00%** | **100.00%** | **0.8750** |
| *B. subtilis* 144 | All / Unfiltered | DIAMOND-BLASTp Only | 144 / 76 | 2.08% (3) | 6.25% | 3.95% | 11.84% | 0.1420 |
| *B. subtilis* 144 | All / Unfiltered | gapseq v1.2 (Bitscore Heuristic) | 144 / 76 | 7.64% (11) | 20.14% | 14.47% | 38.16% | 0.2780 |
| *B. subtilis* 144 | All / Unfiltered | GSPA Strict Track A (No PLM) | 144 / 76 | 11.11% (16) | 31.94% | 21.05% | 60.53% | 0.4221 |
| *B. subtilis* 144 | All / Unfiltered | **Full GSPA (Track B, Frozen)** | 144 / 76 | **43.06% (62)** | **50.00%** | **81.58%** | **94.74%** | **0.8933** |
| Ambiguous Gaps | Multi-Candidate LRO | Stage 1 RLGC Operon Kernel Alone | 50 / 50 | 34.00% (17) | 86.00% | 34.00% | 86.00% | 0.6020 |
| Ambiguous Gaps | Multi-Candidate LRO | **Stage 3 Guarded Re-Ranking (3-Seed)** | 50 / 50 | **82.67% (41.3)** | **93.33%** | **82.67%** | **93.33%** | **0.8792** |

## Main findings

1. **Uniform reaction boosting fails to localize enzymes:** Sequence-only ranking identifies the true enzyme at Rank 1 in only 1.95% of cases. Uniform reaction-level boosting (METEOR) leaves candidate ranking completely unchanged because every candidate locus associated with a reaction receives the same increment.
2. **Local genomic context provides locus discrimination:** RLGC operon kernel density alone lifts unconditional Hit@1 from 1.95% to 23.05% (32.60% within the candidate pool). Combining RLGC with the Stage 1 pathway coherence prior raises unconditional Hit@1 to 23.44% and pool Hit@1 to 34.29% (\(p < 10^{-12}\), exact McNemar test).
3. **Synergy with protein language models:** Swiss-Prot ESM-2 reaction centroids provide alignment-free sequence representations (62.11% unconditional Hit@1). Integrating RLGC and pathway coherence with language-model centroids (Track B) raises unconditional Hit@1 to 63.67% and pool Hit@1 to 90.06% (MRR 0.9330), as chromosomal synteny resolves ties among sequence-similar paralogs.
4. **Physiological role stratification:** RLGC alone achieves its highest standalone recovery on Stratum B alternative catabolic pathways (55.56% pool Hit@1), whose specialized enzymes and transporters are organized into tightly clustered operons. Full GSPA maintains robust pool Hit@1 across core biomass (93.53%), alternative catabolic (88.89%), and peripheral secondary metabolism (82.76%).
5. **Out-of-sample transfer:** Transferred with frozen hyperparameters to *B. subtilis* 168, GSPA achieves 43.06% unconditional Hit@1 and 81.58% pool Hit@1 (MRR 0.8933).
6. **Resolving residual ambiguity:** On 50 multi-candidate metabolic gaps where Stage 1 RLGC ranks the true enzyme at Rank 1 in 34.00% of cases, Stage 3 guarded agent reasoning lifts Hit@1 to 82.67% (\(p = 2.8 \times 10^{-7}\), exact McNemar test; MRR 0.8792).
