# Phase 1 — gLM operon caller results

Generated 2026-05-05 from gLM v8473041 + ESM2-650M on ORIX (1× H200).
Operon predictor: shipped sklearn LogisticRegression on 19 layers × 10 heads
attention contacts (190-dim feature per adjacent-pair). Mode = real;
inference time ~2 min per genome.

## Per-genome operon stats

| Genome | Genes | Operons | Mean size | Mean confidence |
|---|---:|---:|---:|---:|
| bsubtilis     | 4240 | 1010 | 2.83 | 0.846 |
| ecoli         | 4303 | 1094 | 2.76 | 0.856 |
| ecolo157      | 5156 | 1282 | 2.90 | 0.845 |
| hpylori       | 1450 |  381 | 2.86 | 0.810 |
| mgenitalium   |  504 |  106 | 3.13 | 0.808 |
| mjannaschii   | 1811 |  484 | 2.83 | 0.842 |
| mtb           | 3906 | 1028 | 2.94 | 0.877 |
| paeruginosa   | 5573 | 1453 | 2.90 | 0.858 |
| rprowazekii   |  834 |  226 | 2.94 | 0.836 |
| saureus       | 2767 |  671 | 2.90 | 0.841 |
| synechocystis | 3635 |  993 | 2.85 | 0.860 |
| tpallidum     | 1009 |  256 | 2.90 | 0.831 |
| vcholerae     | 3591 |  903 | 2.84 | 0.863 |

Total: 9,887 operons across 13 genomes.

## Outputs (per genome)

- `operons.tsv`             — tab-sep RefSeq protein IDs per operon, drop-in for `make_operons.py`
- `operons_confidence.tsv`  — operon_idx, size, confidence (mean of 1-P(boundary) over internal pairs)
- `operons_centroids.npz`   *(on ORIX, not committed — 4.8 MB / genome)* — gLM contextualized centroid per operon
- `protein_embeddings.npz`  *(on ORIX, not committed — 41 MB / genome)* — ESM2 + gLM contextualized per protein

`*.npz` files live under `/mnt/data/u/hohndor/gspa-glm/phase1/preds/<tag>/` on ORIX
(consumed by phase 2 / phase 3 — not needed for the 1.G F-max delta).

## Validation

B. subtilis 168 sanity check (task 1.D):
- Mean operon size 2.83 — matches typical bacterial operon distribution
- Largest operons: 20, 16, 14, 9 genes — all on plausible polycistronic loci
- Mean confidence 0.846 — well above the 0.5 default emission threshold
- Distribution: 511 size-2, 297 size-3, 135 size-4, 78 size-≥5

## ID convention

All `operons.tsv` files use **NCBI RefSeq protein IDs** (e.g. `NP_387882.1`,
`YP_003097713.1`) since the gLM caller was driven by RefSeq FAA + GFF.
For the 1.E integrate step the existing `build_refseq_uniprot_map.py`
remap is required before feeding `--operons` (matches `make_operons.py` /
`run_integrate_full_priors.sh` workflow).
