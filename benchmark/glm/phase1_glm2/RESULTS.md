# Phase 1 — gLM2 operon caller F-max ablation

Run date: 2026-05-05.
gLM2: `tattabio/gLM2_650M` (HF), commit `ee8cea27` (TattaBio/gLM2 repo).
Reference paper: Cornman et al., bioRxiv 2024.08.14.607850 (ICLR 2025).

This is a redo of the gLM phase-1 ablation with **gLM2** as the
foundation model instead of gLM (Hwang & Ovchinnikov 2024). Same 8
benchmark genomes, same heuristic baseline, same integrate pipeline,
same F-max protocol. Only the FM that produces the operons differs.

## Methodology

- **Model:** gLM2 650M (`tattabio/gLM2_650M`), `bfloat16` on H200,
  4096-token context, mixed-modality input
  (`<+>AA<+>nucleotide_igs<->AA<->...`).
- **Operon segmentation:** sliding window of ≤12 genes per contig
  (4096 tokens at AA truncation 1000 + IGS truncation 200), overlap 1
  gene per window. For each adjacent same-contig same-strand pair,
  cosine similarity between the gLM2 contextualized per-protein
  embeddings (mean-pool over each protein's tokens). Cosine → P(break)
  via sigmoid centered at 0.85, sharpness 12.
- The center / sharpness defaults were chosen empirically so the
  operon size distribution lands in the same neighbourhood as gLM /
  the heuristic baseline (~3 genes / operon mean). A formal logreg
  calibration on E. coli operon ground truth (analogous to gLM's
  shipped predictor) is feasible but skipped here — the F-max signal
  on the comparison ended up well below noise either way.
- All other ablation knobs identical to gLM phase 1: same 8 of 9
  bench9 genomes (synechocystis dropped, same reason), same RefSeq →
  UniProt remap, same `--operon-caller heuristic` baseline run, same
  `gspa integrate` invocation with the full prior stack, same
  `benchmark_pgap_v2.py` against `truth_dual/<tag>_truth_{exp,all}.tsv`.

## Per-genome operon stats (post-remap to UniProt)

| Genome | gLM2 operons | mean size | gLM operons | heur operons |
|---|---:|---:|---:|---:|
| bsubtilis     |  711 | 3.78 | 1010 |  748 |
| ecoli         |  508 | 2.64 | 1094 | 1764 |
| ecolo157      |  606 | 3.08 | 1282 | 1775 |
| hpylori       |  269 | 3.66 |  381 |  ~?  |
| mgenitalium   |   86 | 4.73 |  106 |  ~?  |
| mjannaschii   |  269 | 3.44 |  484 |  ~?  |
| mtb           |  822 | 3.20 | 1028 |  ~?  |
| paeruginosa   | 1117 | 3.20 | 1453 |  ~?  |

gLM2 operon counts are between the heuristic baseline and gLM phase 1.

## F-max delta vs heuristic baseline

### truth = exp (Swiss-Prot experimentally-validated)

| Genome | μ heur | μ glm2 | Δμ | CAFA heur | CAFA glm2 | Δ CAFA |
|---|---:|---:|---:|---:|---:|---:|
| bsubtilis | 0.271 | 0.271 | +0.000 | 0.228 | 0.228 | +0.000 |
| ecoli | 0.408 | 0.408 | −0.000 | 0.390 | 0.390 | +0.000 |
| ecolo157 | 0.269 | 0.269 | +0.000 | 0.193 | 0.193 | +0.000 |
| hpylori | 0.147 | 0.147 | +0.000 | 0.087 | 0.087 | +0.000 |
| mgenitalium | 0.242 | 0.242 | +0.000 | 0.334 | 0.334 | +0.000 |
| mjannaschii | 0.352 | 0.352 | +0.000 | 0.320 | 0.320 | +0.000 |
| mtb | 0.336 | 0.336 | −0.000 | 0.243 | 0.243 | −0.000 |
| paeruginosa | 0.222 | 0.222 | −0.000 | 0.221 | 0.221 | +0.000 |
| **mean** | — | — | **−0.0000** | — | — | **−0.0000** |

Worst-genome Δμ: mtb (−0.0002).

### truth = all (Swiss-Prot all annotations, including IEA)

Mean Δμ = **−0.0000**, mean Δ CAFA = **+0.0000**, worst-genome Δμ =
mgenitalium (−0.0002). Per-genome table is identical-to-noise across
all 8 genomes for both metrics.

## Verdict

> **Bar:** mean micro F-max delta ≥ +0.005 across the genomes
> **and** no genome regresses by > 0.01.

- Mean Δ μ F-max (truth=exp) = **−0.0000** — fails the +0.005 threshold.
- Worst-genome Δ μ F-max = **−0.0002** — passes the −0.01 floor.

→ **Verdict: NO-GO.** Same outcome as gLM phase 1 (mean Δμ = +0.0002).

## Sanity / wiring

- 711 gLM2 operons loaded for bsubtilis; GenomicContextPrior boosted
  416 claims (vs 459 heuristic, 397 gLM).
- bsubtilis integrated TSV differs from heuristic on 854 of 28,075
  rows (3.0%) — slightly more divergence than gLM (798 rows, 2.8%).
- The wiring is correct, the null result is real for both gLMs.

## Comparison: gLM vs gLM2 vs heuristic (bsubtilis as representative)

| Metric | heuristic | gLM | gLM2 |
|---|---:|---:|---:|
| operons loaded into integrate | 748 | 883 | 711 |
| GenomicContextPrior claims boosted | 459 | 397 | 416 |
| rows differing from heuristic | (baseline) | 798 | 854 |
| μ F-max (truth=exp) | 0.271 | 0.271 | 0.271 |
| Δμ vs heuristic | — | +0.000 | +0.000 |

Both FMs produce *different* operons that *change* the integrated
posteriors, but neither changes them *toward truth* on these
Swiss-Prot-rich Phase-7 benchmark genomes.

## Why it didn't move (post-mortem, mostly identical to gLM phase 1)

The phase-1 post-mortem applies unchanged to gLM2:

1. **GenomicContextPrior is dominated by the homology stack.** With
   DIAMOND + Pfam + InterProScan already running, swapping the operon
   set is at best a soft re-weighting of an already-saturated
   posterior.
2. **Both FMs are *more selective* than the heuristic.** gLM had
   ~−42 ctx-fired claims/genome, gLM2 has ~−30. Fewer propagation
   events → fewer chances to either help or hurt.
3. **Both FMs were applied with downstream BP-transfer unchanged.**
   The transfer logic in `OperonPredictor.transferAnnotations` is the
   bottleneck the FMs cannot reach — the SPEC's design choice to
   keep that fixed for step-1 ablation purity is what makes the swap
   look like a no-op when the bottleneck is the transfer, not the
   detection.
4. **Calibration of the cos→P(break) sigmoid is approximate.** A
   formal logreg trained on E. coli operon ground truth (gLM's
   approach) might tighten the gLM2 operon set further but would not
   plausibly move μ F-max from `≈ 0` to `≥ +0.005` given the
   homology-stack saturation in (1).

## Decision

**NO-GO on phase 2 / phase 3 — same conclusion as gLM.** Two
independent foundation models (gLM, gLM2) produce different but
F-max-equivalent operon sets at the integrate level. The line is
closed.

If the gLM2 path is ever re-opened, the highest-leverage changes
would be (in priority order):
1. Run on the **dark-matter slice** (proteins with no DIAMOND/Pfam/IPS
   hit) — that's where the homology stack does *not* dominate, so the
   operon channel could matter.
2. Replace `transferAnnotations`'s coarse "any BP term spreads to
   operon members" with a **per-claim weighting** by gLM2 operon
   confidence.
3. Train a proper **logreg on E. coli operon ground truth** so the
   gLM2 operon set is calibrated rather than thresholded by a
   hand-tuned sigmoid center.
4. Bring **GENOMIC_LANGUAGE_MODEL** evidence (gLM2-derived per-protein
   function predictions, the original phase-2 plan) into Phase 7 as
   an *additional* claim type rather than as a re-weighting of an
   existing one.

## Artifacts

- `preds/<tag>/operons.tsv`            — gLM2 operons (RefSeq IDs)
- `preds/<tag>/operons_confidence.tsv` — per-operon confidence
- `results/<tag>_glm2_fmax.json`       — F-max with bootstrap CIs
- `/data/hohndor/gspa-glm/full_priors_glm2/<tag>_integrated.tsv` (unimatrix01)
- `/mnt/data/u/hohndor/gspa-glm/phase1_glm2/preds/<tag>/operons_centroids.npz` (ORIX, ~5 MB / genome)
- `/mnt/data/u/hohndor/gspa-glm/phase1_glm2/preds/<tag>/protein_embeddings.npz` (ORIX, ~40 MB / genome)
