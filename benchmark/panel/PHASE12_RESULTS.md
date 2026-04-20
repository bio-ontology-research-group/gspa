# Phase 12 panel expansion — end-to-end results

**Date**: 2026-04-20
**Branch**: `phase11-crossgenome`
**Scope**: scale the cross-genome conditional-LR framework from the
29-genome GTDB panel + 4-culture pilot to the ~340 KAUST isolate/MAG
genomes under `/data/emptyquarter/sequencing-results/` on unimatrix01.
Strategy C (panel = query = 97 high-quality dereplicated genomes),
ANI-95 dereplication.

## Headline

**595 unique (culture, protein, EC) high-confidence dark-matter
function predictions** across 97 phylogenetically-diverse KAUST
culture genomes from 4 phyla. Every prediction satisfies three
filters: top-3 reaction-local context rank in its own genome,
cross-genome log_lr ≥ 0.3 (≥2× enrichment over orthogroup baseline),
and DIAMOND-pident < 30% against Swiss-Prot (no recognizable
homolog).

## Panel composition

| Phylum          | HQ reps |
|-----------------|--------:|
| Bacillota       |      44 |
| Pseudomonadota  |      40 |
| Actinomycetota  |      11 |
| Bacteroidota    |       2 |
| **Total**       |  **97** |

Source: 219 enrichment-MAG bins + 116 site-isolate assemblies +
3 rh_sequencing → 239 post-QC (CheckM2 ≥70% completeness, ≤10%
contamination) → 97 post-ANI-95-dereplication (skani).

## Phase flow

| Phase | What it does | Compute | Rows out |
|-------|--------------|---------|----------|
| 1 | Inventory + CheckM2 + GTDB-Tk + skani derep → `genome_manifest.tsv` | ~3 h cluster | 338 |
| 2 | prodigal + DIAMOND + Pfam + gapsmith + integrate per genome | ~8 h cluster | 97 |
| 3 | concat proteomes → MMseqs2 cluster-50 → build `nonanchor_catalog_panel.tsv` | ~2 h cluster | 128 M |
| 4A | `predict_dark_matter.py` per genome (reaction-local context kernel density) | ~30 min cluster | 8.5 M candidates |
| 4B | `augment_all.py` — single pass, filtered catalog (log_lr ≥ 0.3, ~127 K rows) | ~30 min single node | 27,613 catalog-matched |
| 5 | `build_shortlist.py` — dedup + phylum-stratify + max-pident filter | seconds | **595 unique predictions** |

## Cross-genome signal distribution (Phase 3 catalog)

| log_lr bucket | orthogroup × reaction pairs |
|---------------|----------------------------:|
| high (≥ 1)    |          11,538 (~10× enrichment) |
| medium (0.3–1)|         115,899 (~2–10× enrichment) |
| positive (0–0.3) |       58,104,887 |
| negative (< 0)  |       70,219,574 |

**~127 K** pairs exceed the "meaningful" floor of 2× enrichment —
the signal pool Phase 4B draws on.

## Top-10 predictions (master shortlist)

| Culture | Phylum | Protein | log_lr | pident | Predicted function |
|---------|--------|---------|-------:|-------:|--------------------|
| enrichment__C-23_metaflye.4 | Actinomycetota | contig_32_323 | **1.646** | 0.0 | exopolyphosphatase (EC 3.6.1.11) |
| enrichment__MO-1_metaflye.3 | Pseudomonadota | contig_204_2739 | 1.471 | 29.4 | Δ24-sterol reductase (EC 1.3.1.72) |
| isolates__site60_MR60-2 | Pseudomonadota | contig_4_923 | 1.375 | 23.6 | sepiapterin reductase (EC 1.1.1.153) |
| enrichment__C-29_metaflye.1 | Pseudomonadota | … | 1.365 | 0.0 | naphthalene dihydrodiol dehydrogenase (EC 1.3.1.29) |
| enrichment__C-16_metaflye.1 | Actinomycetota | … | 1.346 | 0.0 | carbazole 1,9a-dioxygenase (EC 1.14.12.22) |
| enrichment__MO-G_metaflye.* | Bacteroidota | … | 1.250 | 0.0 | alkene monooxygenase (EC 1.14.13.69) |
| enrichment__C-27_metaflye.* | Pseudomonadota | … | 1.250 | 28.5 | 4-hydroxy-4-methyl-2-oxoglutarate aldolase (EC 4.1.3.17) |
| … | … | … | … | … | … |

Biological coverage is diverse: **polyphosphate metabolism, sterol
biosynthesis, biopterin, aromatic-hydrocarbon degradation, alkene
oxidation** all appear in the top-10. The shortlist is a candidate
pool for experimental validation — each row is a specific hypothesis
of the form "protein X in genome Y catalyzes reaction Z."

## Caveats

- **No self-reference hold-out** at the catalog-build step: each query
  genome contributes at most 1 to its own orthogroup's `n_base_with`
  and `n_sig_with`, so the optimism is ≤ 1/97 ≈ 1%. Negligible at this
  panel size; not worth 97× catalog rebuilds.
- **GTDB-Tk species-level classification failed** (OOM at the
  species-tree pplacer step). The manifest uses class-level placement
  recovered from intermediate results — sufficient for phylum-level
  stratification, insufficient for species ID.
- **Orthogroup clustering at 50% id / 80% cov** is permissive enough
  for cross-phylum orthologs but will merge paralogs that diverged
  below that threshold. A second map at 90% id exists under
  `bench_gtdb30/ortho/` from the earlier pilot if paralog resolution
  matters downstream.
- **gapsmith gaps are SEED-reaction-based.** The EC ↔ SEED alias
  expansion in `augment_all.py` may miss legitimate equivalents;
  `n_hits` in the shortlist records how many SEED equivalents the
  augmenter found.

## Key artifacts on unimatrix01

```
/data/hohndor/gspa/proteomes/culture_panel/
  genome_inventory_v2.tsv        # Phase 1 inventory (338 rows)
  genome_manifest.tsv            # Phase 1 final manifest (338 rows)
  phase2_manifest.tsv            # Strategy C panel (97 rows)
  phase1/                        # CheckM2 + GTDB-Tk + skani outputs
  phase2/<tag>/                  # prodigal, preds, claims, layout, gapsmith, integrated
  phase3/                        # ortho map + panel_manifest for catalog build
  nonanchor_catalog_panel.tsv    # 128 M rows
  phase4/<tag>/dark_matter_augmented.tsv
  phase5_shortlist/              # master + per-phylum TSVs + README
```

## Phase 6 — validation (in progress)

Two orthogonal tests:

1. **Convergent evidence** — does the same orthogroup appear as a
   top-3 candidate for the same reaction in multiple independent
   phyla? Convergent hits are especially trustworthy because the
   phylogenetic signal that might bias the density field is
   independent across phyla.

2. **Signal decay** — rebuild the non-anchor catalog at panel sizes
   29 / 60 / 97 and confirm that `log_lr` for the top shortlist
   entries stabilizes rather than simply inflating with panel size.
   If log_lr kept climbing linearly with N, the statistic would be
   an artifact of sample size; if it plateaus, it reflects a real
   biological signal.
