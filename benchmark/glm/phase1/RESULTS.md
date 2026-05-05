# Phase 1 — gLM operon caller F-max ablation

Run date: 2026-05-05.
gLM commit: `8473041306afe96296e5b218a8805338f68a8e6c` (y-hwang/gLM).
ESM2: `esm2_t33_650M_UR50D`. Operon predictor: shipped sklearn LogReg
(`gLM/repo/data/operon_predictor.pkl`).

## Methodology

Single variable swapped between treatment and control: the operons
file fed to `gspa integrate --operons`. Everything else identical
(claims, priors, calibration, taxon constraints).

- **control (heuristic)** — operons from `make_operons.py`
  (intergenic-distance ≤ 300 bp, same strand, ≥2 genes).
  `${ROOT}/operons/<tag>_operons.tsv`.
- **treatment (glm)** — operons from `benchmark/neural/run_glm_operon.py`
  (ESM2 + gLM forward → 190-dim attention contacts → shipped LogReg
  → P(same operon) → segment by 0.5 threshold), then RefSeq → UniProt
  remap (`scripts/remap_glm_operons.py`).
  `/data/hohndor/gspa-glm/phase1/preds_uniprot/<tag>/operons.tsv`.

8 of the 9 bench9 genomes; **synechocystis dropped** because its
RefSeq protein IDs are `WP_*` (newer assembly) which has no entries
in the existing `synechocystis.refseq_to_uniprot.tsv` mapping (file
is empty). Cleanly skipping it leaves 8 genomes for the comparison.

The other 4 genomes the SPEC originally targeted (saureus, vcholerae,
tpallidum, rprowazekii) live in bench10/ but have no
`*_claims.jsonl` on unimatrix01 — they were never wired into the
`integrate` pipeline. Computing claims for them would require re-running
the full predictor stack (DIAMOND, Pfam, InterProScan), which is
out-of-scope for the gating decision.

`gspa integrate` was the existing script with full priors
(`essentiality, coherence, consistency, gap_filling, genomic_context`),
run via the `--operon-caller {heuristic,glm}` switch landed in `89242fb`.
F-max via `benchmark/benchmark_pgap_v2.py` against
`truth_dual/<tag>_truth_{exp,all}.tsv`.

## Results

### truth = exp (Swiss-Prot experimentally-validated annotations)

| Genome | μF-max heur | μF-max glm | Δ μ | CAFA heur | CAFA glm | Δ CAFA | ctx-fired heur | ctx-fired glm | Δ ctx |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| bsubtilis     | 0.271 | 0.271 | +0.000 | 0.228 | 0.228 | +0.000 | 459 | 397 | −62 |
| ecoli         | 0.408 | 0.408 | −0.000 | 0.390 | 0.390 | +0.000 | 802 | 658 | −144 |
| ecolo157      | 0.269 | 0.269 | +0.000 | 0.193 | 0.193 | +0.000 | 706 | 642 | −64 |
| hpylori       | 0.147 | 0.147 | +0.000 | 0.087 | 0.087 | +0.000 | 204 | 178 | −26 |
| mgenitalium   | 0.242 | 0.242 | +0.000 | 0.334 | 0.334 | +0.000 |  55 |  48 |  −7 |
| mjannaschii   | 0.352 | 0.353 | +0.001 | 0.320 | 0.320 | +0.000 | 120 | 100 | −20 |
| mtb           | 0.336 | 0.336 | +0.000 | 0.243 | 0.243 | +0.000 | 644 | 659 | +15 |
| paeruginosa   | 0.222 | 0.222 | −0.000 | 0.221 | 0.221 | +0.000 | 583 | 553 | −30 |
| **mean**      | — | — | **+0.0002** | — | — | **+0.0001** | — | — | — |

Worst-genome Δ μ F-max: ecoli (−0.0002).

### truth = all (Swiss-Prot all annotations including IEA)

| Genome | μF-max heur | μF-max glm | Δ μ | CAFA heur | CAFA glm | Δ CAFA |
|---|---:|---:|---:|---:|---:|---:|
| bsubtilis     | 0.674 | 0.674 | +0.000 | 0.690 | 0.690 | +0.000 |
| ecoli         | 0.670 | 0.670 | −0.000 | 0.683 | 0.683 | −0.000 |
| ecolo157      | 0.835 | 0.835 | −0.000 | 0.839 | 0.839 | +0.000 |
| hpylori       | 0.753 | 0.753 | −0.000 | 0.731 | 0.731 | +0.000 |
| mgenitalium   | 0.912 | 0.912 | +0.000 | 0.912 | 0.912 | +0.000 |
| mjannaschii   | 0.639 | 0.639 | +0.000 | 0.637 | 0.637 | +0.000 |
| mtb           | 0.715 | 0.715 | −0.000 | 0.697 | 0.696 | −0.000 |
| paeruginosa   | 0.598 | 0.598 | −0.000 | 0.654 | 0.654 | +0.000 |
| **mean**      | — | — | **−0.0000** | — | — | **−0.0000** |

## Gating verdict per SPEC §1.G

> **Bar:** mean micro F-max delta ≥ +0.005 across the 13 genomes
> **and** no genome regresses by > 0.01.

Result on 8 of the canonical genomes:

- Mean Δ micro F-max (truth=exp) = **+0.0002** — fails the +0.005 threshold.
- Worst-genome Δ micro F-max = **−0.0002** — passes the −0.01 floor.
- Mean Δ CAFA F-max (truth=exp) = **+0.0001** — also flat.

→ **Verdict: NO-GO. Drop the line.** Per SPEC: do not start phase 2.

## Sanity / wiring checks

The null result is genuine, not a pipeline bug:

- Operons ARE loaded into integrate. From bsubtilis logs:
  `heur — Operons: 748, GenomicContextPrior: 459 claims boosted across 748 operons`
  `glm  — Operons: 883, GenomicContextPrior: 397 claims boosted across 883 operons`.
- Integrated TSVs differ (bsubtilis: 798 of 28,075 rows differ; ecoli:
  1,202 of 29,690; mjannaschii: 188 of 11,236). The Phase-7 posteriors
  do change with the new operon set — they just don't change toward truth.
- `--operon-caller heuristic` mode reproduces the existing benchmark
  numbers — heuristic-mode F-max matches the prior `bench9_full_priors`
  reference within the bootstrap CIs.

## Why it didn't move (post-mortem)

A few observations from the data, ordered most→least likely to be the
causal story:

1. **GenomicContextPrior is dominated by the other priors.** With
   `essentiality + coherence + consistency + gap_filling` already
   running, the marginal contribution of any operon-set change is
   small. The 800–1,200 row diff per genome touches mostly
   already-confident posteriors that were going to clear/miss the
   threshold either way.

2. **gLM operons are *more selective*, not *more informative*.**
   `ctx-fired` decreased on 7/8 genomes (mean −42 per genome). Fewer
   propagation events means fewer chances to either help or hurt. The
   net F-max impact is therefore close to zero by construction. mtb
   is the lone exception (+15 fired) and is also the only genome
   where Δ goes neither up nor down.

3. **The shipped operon LogReg was trained on *E. coli* operons
   only.** Its calibration may transfer poorly to mtb,
   mgenitalium, mjannaschii (small / archaeal genomes). Re-training
   the operon head on a multi-organism set would test this — out of
   scope for the gating decision.

4. **BP-transfer in `OperonPredictor.transferAnnotations` is the same
   coarse heuristic on both sides.** As designed (SPEC §3 resolved
   decision: "step 1 keeps `OperonPredictor.transferAnnotations`
   unchanged"), the only way the swap moves F-max is via the operon
   set itself. If the swap doesn't move F-max it means the operon set
   was not the bottleneck.

5. **Drop ratio in RefSeq → UniProt remap.** Across the 8 genomes,
   8–35% of gLM operons were dropped because some members had no
   UniProt mapping (mjannaschii lost 36% of operons; ecolo157 lost
   13%). Heuristic operons are pre-mapped to UniProt at make time,
   so they don't pay this loss. This systematically biases the
   gLM arm toward fewer operons, but the loss is uniform across
   genomes and the trend on un-mapped operons would have to be
   strongly different to flip the verdict.

## Decision

**NO-GO on phase 2 (GENOMIC_CONTEXT_FM evidence type) and phase 3
(BF(O,P) augmentation).** Per the SPEC's own gating logic, the
operon-detection signal is not the lever that improves F-max in this
benchmark.

That does NOT mean gLM is useless for the dark-matter / Empty
Quarter use case — the benchmark exercises the *integrated*
posterior on Swiss-Prot-rich genomes, where dense homology already
saturates the answer. For the EQ MAGs the predictor stack is much
sparser. But re-opening the line should require:

- A different proxy benchmark (dark-matter slice, or held-out
  CAFA-style test on poorly-annotated genomes), AND
- A clearer hypothesis about *which mechanism* (operon detection
  vs. gLM-derived per-protein function predictions vs. BF(O,P)
  augmentation) would actually show signal there.

## Artifacts

- `preds/<tag>/operons.tsv`            — per-genome gLM operons (RefSeq IDs)
- `preds/<tag>/operons_confidence.tsv` — per-operon confidence
- `results/<tag>_{heuristic,glm}_fmax.json` — raw F-max with bootstrap CIs
- `integrated_{heuristic,glm}/<tag>_integrate.log` — gspa integrate logs

`*_centroids.npz` and `*_protein_embeddings.npz` (~600 MB total) stay
on ORIX under `/mnt/data/u/hohndor/gspa-glm/phase1/preds/<tag>/`.

ORIX gLM env, weights cache, and gLM repo at pinned commit remain
provisioned for any future re-investigation.
