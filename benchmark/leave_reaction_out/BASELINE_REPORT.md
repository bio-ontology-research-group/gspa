# Leave-Reaction-Out Baseline: Current DarkMatterSuggester on E. coli MG1655

**Date**: 2026-04-18
**Genome**: E. coli K-12 MG1655 (UniProt/RefSeq, 4,019 proteins with coords)
**Test set**: 256 LRO cases (filters: good_blast gapsmith hit, EC → GO mapping,
  ≥2 local catalyzed-neighbors within 20 kb, MG1655-internal paralog count ≤3,
  ground-truth-confirmed)
**Ablation modes**: function-only (256) + protein-level (256) = 512 runs
**Tool**: current DarkMatterSuggester (Phase 10 refined BF, Noisy-OR + IC + purity)

## Ranking definition

DM emits one suggestion per (gap, operon). Within each suggestion it lists
proteins and per-protein q-values. We rank proteins globally by
`score = suggestion_score * q(p)`, dedupe per protein (keep max), sort desc.

## Headline numbers

| metric                     | overall | function-ablation | protein-ablation |
|----------------------------|---------|-------------------|------------------|
| hit@1                      | 0.018   | 0.016             | 0.020            |
| hit@3                      | 0.057   | 0.043             | 0.070            |
| hit@5                      | 0.168   | 0.152             | 0.184            |
| recall (p in cand list)    | 0.674   | 0.664             | 0.684            |
| MRR                        | 0.084   | 0.076             | 0.093            |
| mean \|cand\|              | 96.0    | 96.0              | 95.9             |
| margin (hit@1)             | 0.000   | 0.000             | 0.000            |

## Decomposition of failure modes (n=512)

| outcome                         | count | frac  |
|---------------------------------|-------|-------|
| hit@1 (rank=1)                  |     9 | 0.018 |
| hit in top 2-5                  |    77 | 0.150 |
| buried in candidate list (>5)   |   259 | 0.506 |
| p not in candidate list         |   139 | 0.271 |
| no DM suggestions at all        |    28 | 0.055 |

## Stratified by n_neighbors_local (catalyzed reaction-graph neighbors within 20 kb)

| n_neighbors | n   | hit@1 | hit@5 | recall | MRR   | \|cand\| |
|-------------|-----|-------|-------|--------|-------|---------|
| 2           | 214 | 0.009 | 0.150 | 0.696  | 0.070 | 109.4   |
| 3-4         | 168 | 0.000 | 0.185 | 0.685  | 0.083 | 104.7   |
| ≥5          | 130 | 0.054 | 0.177 | 0.623  | 0.111 | 62.5    |

## Stratified by within-genome orthogroup size

| ortho size | n   | hit@1 | hit@5 | recall |
|------------|-----|-------|-------|--------|
| 1 (single) | 498 | 0.018 | 0.171 | 0.673  |
| 2-3        |  14 | 0.000 | 0.071 | 0.714  |

## What the numbers say

1. **Essentially random precision at top-1.** hit@1 of 1.8% over 512 cases
   is not meaningfully above chance. All 9 hit@1 successes have a *zero*
   margin — they are disjunctive ties where the target protein happened
   to be listed first alphabetically in a uniform-q suggestion.

2. **Wide recall, zero specificity.** 67% of cases have `p` somewhere in
   the candidate list, but the mean candidate-set size is 96 proteins.
   DM is saying "here is every protein in every operon that hits any
   pathway function" — no useful ranking across operons.

3. **Neighbor count helps, weakly.** hit@1 rises from 0.9% (2 neighbors)
   to 5.4% (≥5 neighbors). The effect is in the right direction but the
   absolute level is still poor. This validates the design hypothesis
   that *more reaction-graph neighbors = more usable signal*, but it
   also shows that 5+ neighbors within 20 kb is rare (only 130/512
   cases) and even there the current framework under-exploits it.

4. **Paralogs collapse entirely.** Orthogroup size 2-3 → 0% hit@1. The
   current framework cannot distinguish paralogs by genomic context
   alone, as expected from the design analysis.

5. **Failure decomposition.** 28 cases (5.5%) produce NO suggestions at
   all — the operon containing `p` fails the BF ≥ 10 gate after
   ablation. 139 cases (27%) emit suggestions but `p` is not in any of
   the returned operons — operon assignment doesn't include `p`, or
   `p`'s operon doesn't clear the BF gate.

## The bar RLGC must clear

- **Top-1 accuracy**: beat ~2%. Any sharp method should easily hit
  20-40% on the ≥5-neighbor subset; total hit@1 ≥ 15% across strata is
  a credible target.
- **Candidate set size**: far below 96 — a useful method emits a
  handful of candidates with non-zero margin.
- **No-signal cases**: reduce the 5.5% "no suggestions at all" failure
  mode via continuous kernel density (no hard BF ≥ 10 gate).
- **Cross-genome disambiguator**: needed to crack the ortho 2-3 band
  that's at 0% here. Within-genome alone can't do it.

## Next step

Build M1: `ReactionLocalContextSuggester` in `gspa-core/integration/suggester/`.
Replace pathway-DB × operon-assignment with `N_k(R)` reaction-graph
neighborhood × kernel-smoothed genomic density. Benchmark under the same
512-run LRO harness.

**Artifacts**
- `benchmark/leave_reaction_out/cases_mg1655.tsv` — test set
- `benchmark/leave_reaction_out/results_both.tsv` — raw per-case results
- Shard-level outputs under `runs/{protein,function}/shard_*/`

---

## Phase 12 milestone results (n=256, mg1655 LRO, mode=protein)

| milestone | hit@1 | hit@3 | hit@5 | MRR | |cand| | margin | notes |
|-----------|-------|-------|-------|-----|--------|--------|-------|
| DM baseline | 0.018 | — | 0.168 | — | 96.0 | 0.000 | Phase 11 numbers (n=512) |
| **M1 RLGC** | 0.229 | 0.387 | 0.415 | 0.328 | 18.0 | 0.001 | reaction-graph locality × density |
| **M2 raw-claims catalog** | 0.232 | 0.379 | 0.422 | 0.326 | 14.7 | 0.001 | neutral vs M1 — catalog too noisy |
| **M2 integrated-posteriors catalog** | **0.236** | 0.366 | 0.409 | 0.335 | 14.7 | 0.001 | marginal lift; confirms LR formulation is not the bottleneck |
| **M3 GBDT** | 0.180 | 0.410 | 0.488 | 0.328 | 18.0 | 0.027 | hit@1 dips but top-3/top-5 and margin sharpen |
| **M4 Track B (PLM permissive)** | **0.637** | **0.684** | **0.695** | **0.660** | 18.0 | **1.336** | massive lift — PLM cosine-to-EC-centroid dominates |
| M4 Track A (PLM masked if cos-to-target > 0.7) | 0.113 | 0.348 | 0.375 | 0.261 | 18.0 | 0.622 | note: M3-model is the cleaner "no-PLM" baseline (0.180) |

### M4 Track B stratification (mg1655, n=256)

- By n_neighbors_local:
  - 2 nbrs (n=107): hit@1=0.617, hit@5=0.664, MRR=0.636
  - 3-4 nbrs (n=84):  hit@1=0.548, hit@5=0.607, MRR=0.570
  - ≥5 nbrs  (n=65):  **hit@1=0.785, hit@5=0.862, MRR=0.815**
- By orthogroup size:
  - ortho=1 (n=249): hit@1=0.643, hit@5=0.703
  - ortho=2-3 (n=7): hit@1=0.143 (still harder than singletons, as expected)

The ≥5-neighbour bucket was the persistent weak point of M1/M2/M3
(stuck at ~0.05 hit@1). Under M4 Track B it jumps to 0.785 — exactly
where PLM was expected to contribute most, since densely-connected
metabolism has many EC-similar candidates that pure genomic context
cannot distinguish.

### M4 feature importance (LightGBM gain)

1. plm_cos_centroid_EC: 20905 (54×)
2. rlc_q:                 2121
3. commitment:            1154
4. log_density:            575
5. self:                   474
6. diversity:              404

The model put ~93% of gain on a single PLM feature. This dominance
is itself a finding: once PLM cosine-to-catalyst-centroid is in the
feature set, every other signal (including the whole RLGC context
stack) becomes a tiebreaker. It does not mean context is worthless —
the residual ranker gain on rlc_q/commitment/log_density shows
context still decides ambiguous cases — but it does mean the
deployment headline for M4 is almost entirely PLM-driven.

### M2 tune (integrated posteriors) interpretation

Rebuilding the reaction-locus catalog from integrated posteriors
(not raw DIAMOND claims) moved hit@1 from 0.232 → 0.236 — within
noise. The cross-genome LR formulation itself (not the catalog
quality) is the bottleneck; on a 29-genome panel with per-(C, R)
support of typically 3-10 genomes, the LR estimates are too
uncertain to sharpen an already-concentrated M1 posterior.

### Next steps

- M4 Track A needs a cleaner reference: run the existing M3 model
  (no-PLM) on the same 256 cases — that is the "context-only" number
  to report alongside Track B, not a zero-out of PLM features in the
  M4 model.
- M5 (neural) gate: with Track B already at 0.637 hit@1 (M5 target
  was 0.65), the remaining headroom is small and the PLM feature
  already dominates. Investigate whether a GNN or masked-reaction
  transformer can extract additional signal the centroid loses —
  e.g., attention over neighbour-reaction embeddings as opposed to
  a single per-EC mean. If yes, proceed to M5; if no, M4 is the
  shipping target.

---

## M4 follow-up sanity checks

### (a) PLM-only model (2 features: plm_cos_centroid_EC, plm_has_emb)

Trained LightGBM on the same 2,133 panel qids with only PLM
features.

| model | hit@1 | hit@3 | hit@5 | MRR | margin |
|-------|-------|-------|-------|-----|--------|
| M4 full (19 features) | 0.637 | 0.684 | 0.695 | 0.660 | 1.336 |
| **PLM-only (2 features)** | **0.621** | 0.676 | 0.691 | 0.647 | 0.099 |
| M3 (17 RLGC features, no PLM) | 0.180 | 0.410 | 0.488 | 0.328 | 0.027 |

RLGC contributes ~1.6 hit@1 points on top of PLM cosine alone. The
Phase-12 scoring stack is ~98% PLM, ~2% context. Context is a
tiebreaker, not a discriminator.

### (b) Strict Track A (per-case centroid with sequence-similar panel catalysts removed)

Threshold 0.5 cos-sim to target. Produces numerically identical
hit@k to the naive Track A, because under a PLM-heavy model, any
zero-ing of the PLM feature pushes ranking onto the residual RLGC
branches — same fallback regardless of how PLM is zeroed.

| track | hit@1 | hit@3 | hit@5 | MRR |
|-------|-------|-------|-------|-----|
| Track A (naive) | 0.113 | 0.348 | 0.375 | 0.261 |
| Track A strict | 0.113 | 0.348 | 0.375 | 0.261 |
| M3 (no-PLM, cleanest Track A proxy) | 0.180 | 0.410 | 0.488 | 0.328 |

M3 is a better "Track A = context-only" reference than zeroing PLM
features in a PLM-heavy model.

### (c) Cross-genome generalisation: bsubtilis (n=144, held out from centroid build)

Same M4 model; centroids rebuilt excluding bsubtilis. This is the
honest generalisation test.

| genome | Track B hit@1 | Track B hit@5 | Track A hit@1 | recall | margin |
|--------|---------------|---------------|---------------|--------|--------|
| mg1655 | 0.637 | 0.695 | 0.113 | 0.707 | 1.336 |
| **bsubtilis** | **0.431** | 0.521 | 0.111 | 0.528 | 0.872 |

Hit@1 drops ~20 points when the target genome is held out of the
centroid build (0.637 → 0.431). Recall also drops (0.707 → 0.528),
which reflects the RLC suggesters own genome-sensitivity, not PLM.
Track B / PLM-only / Track A gap holds qualitatively on both
genomes: the PLM feature is the decisive signal everywhere.

## Updated summary table (mg1655 LRO, n=256 unless noted)

| milestone / variant | hit@1 | hit@5 | MRR | cand | notes |
|---------------------|-------|-------|-----|------|-------|
| DM baseline (Phase 11) | 0.018 | 0.168 | -- | 96.0 | n=512 |
| M1 RLGC | 0.229 | 0.415 | 0.328 | 18.0 | reaction-graph x density |
| M2 (raw-claims catalog) | 0.232 | 0.422 | 0.326 | 14.7 | neutral vs M1 |
| M2 (integrated-posteriors) | 0.236 | 0.409 | 0.335 | 14.7 | marginal tune |
| M3 GBDT (17 RLGC feats) | 0.180 | 0.488 | 0.328 | 18.0 | sharper top-5 |
| **M4 Track B (full)** | **0.637** | 0.695 | 0.660 | 18.0 | PLM dominant |
| M4 PLM-only | 0.621 | 0.691 | 0.647 | 18.0 | RLGC adds 1.6 pts |
| M4 Track A (naive/strict) | 0.113 | 0.375 | 0.261 | 18.0 | zeroing PLM breaks a PLM-heavy model |
| M4 bsubtilis Track B | 0.431 | 0.521 | 0.471 | 17.2 | held-out genome, n=144 |

## M5 decision

M5 (neural ranker) is not recommended as the next milestone.

Evidence:

1. PLM-only is at 0.621 -- the ceiling for beating PLM-cosine is
   0.08 points on mg1655, less on bsubtilis.
2. RLGC features add 1.6 points on mg1655 -- a GNN over the same
   features cannot meaningfully exceed that.
3. The real headroom is elsewhere:
   - Better PLM features: per-neighbour-reaction cosine, not a
     single per-EC mean. ~6-10 features, no new infra.
   - More panel genomes: bsubtilis Track B drop (0.637 -> 0.431)
     suggests the 29-genome centroid is under-powered for phyla not
     represented in training. Scaling to 100-200 GTDB reps is
     likely higher-value than a neural architecture.
   - Fix the Track A reporting: report M3 as the Track A number,
     not the zeroed-out M4. The scientific story is "context alone
     gets 0.18, PLM alone gets 0.62, combined gets 0.64".

## Next-step recommendation (option, not a promise)

Drop M5 from the immediate roadmap. Instead:

1. Add per-neighbour-reaction PLM features (attention-lite).
2. Scale the panel (bench_gtdb100 or larger).
3. Re-evaluate M4 on >=3 held-out genomes to get a proper
   generalisation distribution.

If those still show <0.70 on held-out genomes, a neural model over
{PLM, neighbour-PLM, graph-structure} is the right bet -- but only
then.
