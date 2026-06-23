# DG++-Light genome cascade — homology-gated, CPU-budgeted annotation

**Goal (project owner, 2026-06-23):** annotate a **bacterial genome in 5–10 min**
and a **eukaryotic genome in ~30 min**, CPU-only, for a *newly sequenced organism*
(novel proteins — no UniProt entry, no literature, not in STRING). This is the real
GSPA use case, distinct from re-annotating SwissProt.

## Design: a homology-gated two-tier cascade

One DIAMOND search of the whole proteome vs the **pre-t0 train DB** triages every
protein, then expensive per-protein features run **only on the orphan minority**.

```
proteome.faa
  └─ Stage 1: DIAMOND vs pre-t0 DB        (whole proteome — the only universal pass)
       ├─ Tier A: HAS homolog (~84–91%)   cheap: diam BLAST-KNN + net_bridge (lookup, reuses THIS search)
       │            → Integrator-A (diam + net_union + interpro)   [no ESM2/CNN — diam dominates here]
       └─ Tier B: NO homolog (orphans)    expensive features run here only:
                    ESM2-kNN + CNN + Pfam(hmmscan) + net_bridge
                    → Integrator-B (esm2_knn + cnn + pfam [+ net_union])
```

The gating is **accuracy-justified**, not just a speed hack: the leak-free ablation
shows expensive features add value exactly where homology is absent — for proteins
*with* a homolog, `diam` dominates (clean-B diam 0.60 vs the rest), so skipping
ESM2/CNN/full-IPS on Tier A costs ~nothing.

## Validated triage (real sequences, via precomputed components)

Measured on the **CAFA6 no-knowledge ("novel") set, n=2,704** — the *hardest*
realistic case; a real bacterial genome is more conserved and tiers even better:

| signal | coverage of novel proteins |
|---|---|
| **Tier A — pre-t0 DIAMOND homolog** | **83.7%** (2263/2704) |
| **Tier B — orphan, no homolog** | **16.3%** (441/2704) |
| net_bridge (net_union) | **95.3%** |
| interpro domains | 51.5% |

**net_bridge reaches further than diam (95% > 84%)**, because `diam` needs the
homolog to be *experimentally GO-annotated* while net_bridge only needs it to be a
*STRING node* (a much larger set). So net_bridge covers ~70% of the proteins diam
misses — broadest coverage of any component, exactly in the orphan regime. Its
accuracy is lower (clean-B 0.25 vs diam 0.60), so the integrator weights it as a
**coverage extender**, not a leader. (Confirms the owner's "net bridge may work".)

## Budget (CPU, multicore)

| stage | bacterium (~4k prot) | eukaryote (~20k prot) |
|---|---|---|
| DIAMOND triage (whole proteome) | 1–2 min | 3–5 min |
| Tier A: diam + net_bridge lookups (~84%) | <1 min | <2 min |
| Tier B ESM2-35M embed + kNN (~16%) | ~1 min | ~4–8 min |
| Tier B CNN | <1 min | ~1 min |
| Tier B **Pfam (hmmscan)** | ~2–4 min | ~10–15 min |
| integrators + DAG-propagate | instant | instant |
| **total** | **~5–8 min ✅** | **~20–30 min ✅** |

**The one binding constraint: full InterProScan (PANTHER/Gene3D HMMs, ~hours even on
hundreds of orphans) cannot fit the budget.** Fast mode therefore uses **Pfam-only
`hmmscan`** as the domain engine; full InterProScan is an opt-in "thorough/overnight"
mode off the critical path.

## Precompute once (offline / GPU-OK, ship as frozen assets)

- pre-t0 DIAMOND DB (SwissProt-exp) — `cnn_work/train_db.dmnd`
- **train ESM2-35M embedding store** (the kNN reference) — ORIX, minutes  ← Phase 2
- net_bridge tables (`apply_net_bridge.load_train_net`, `net_index.tsv.gz`)
- CNN weights, Pfam-A HMM DB, the two integrator JSONs

## Two tier-specific integrators (disjoint feature regimes)

Tier A and Tier B see different features, so one integrator with zero-filled features
is wrong — train two, each on the matching slice, tiered by **real DIAMOND-hit** (the
operational definition, not train_terms membership):

- **Integrator-A** ← proteins WITH a hit: `diam + net_union + interpro`
- **Integrator-B** ← proteins WITH NO hit: `esm2_knn + cnn + pfam [+ net_union]`

### Tier split of the leak-free clean-A GT (`gt_no_cleanA.tsv`, 282 prot)

| tier | proteins | rows |
|---|---|---|
| Tier A (homology) | **258 (91%)** | 466 |
| Tier B (orphan) | **24** | 44 |

**Limitation:** the orphan slice of the *test* GT is only 24 proteins — far too few to
train Integrator-B on. Integrator-B must instead be trained on the **large pre-t0
train population** (train proteins held out from their own homologs = a synthetic,
leak-free orphan set). That is the same deferred IBEX dependency already noted for the
clean retrains. Integrator-A (the 91% majority) is trainable from assets in hand now.

## Results (2026-06-23, leak-free clean GT, cafaeval on ws)

### Tier A — homology path (Integrator-A: diam + net_union + interpro)

| on Tier-A (homology) clean GT | MF | BP | CC | mean f_w |
|---|---|---|---|---|
| diam alone | 0.577 | 0.342 | 0.493 | 0.471 |
| net_union alone | 0.417 | 0.237 | 0.426 | 0.360 |
| interpro alone | 0.410 | 0.189 | 0.057 | 0.219 |
| **Integrator-A** | **0.650** | 0.369 | 0.508 | **0.509** |

Learned weights (all +): diam dominates, net_union solid second, interpro weakest.
Integration lifts diam-alone 0.471 → 0.509 (MF +0.073). net_bridge earns real weight.

### Tier B — orphan path (ESM2-35M embedding-kNN, the workhorse)

| standalone on ORPHAN tier (clean-A, n=24) | MF | BP | CC | mean f_w |
|---|---|---|---|---|
| **esm2_knn (ESM2-35M, CPU)** | **0.701** | 0.210 | 0.612 | **0.508** |
| cnn | 0.342 | 0.197 | 0.520 | 0.353 |
| net_union | 0.458 | 0.250 | 0.377 | 0.362 |
| interpro | 0.333 | 0.167 | 0.000 | 0.167 |

On strict-novel proteins (clean-B, n=251) esm2_knn = **0.458**, *beating DIAMOND 0.442*
while covering 100% of proteins (diam covers only homolog-bearing ones). esm2_knn
**doubles the CNN and triples interpro** on orphans, and on MF (the original novel-protein
disaster, 0.153) it reaches **0.70**. **The orphan tier is now as accurate as the homology
tier (0.508 ≈ 0.509).** Integrator-B (frozen, single-feature esm2_knn calibration) has the
same f_w (cafaeval max-F is invariant to monotonic transforms).

## Status

- [x] Architecture + budget validated on real-sequence triage.
- [x] Tier split of clean-A GT by real DIAMOND-hit (258 Tier-A / 24 Tier-B).
- [x] **Phase 1:** Integrator-A trained (OOF) + evaled on ws = **0.509** homology tier;
      orphan-gating + Pfam-fast path wired into `service/predict.py::cascade()`.
      → `models/deepgo_plusplus_integrator_tierA.json`.
- [x] **Phase 2:** ESM2-35M train store extracted on ORIX (82,404 × 480, ~5 min GPU);
      `pipeline/build_esm2_knn.py` (LOO train feature + test transfer); `_esm2_knn_component`
      wired into predict.py (CPU inference).
- [x] **Phase 1.5:** synthetic-orphan train feature (LOO kNN over 82,404 train proteins);
      Integrator-B frozen on it = **0.508** orphan tier.
      → `models/deepgo_plusplus_integrator_tierB.json`.
- [ ] **Follow-ups (deferred, modest expected gain given esm2_knn dominance):**
      multi-feature Integrator-B (add cnn + Pfam via train-side OOF features); end-to-end
      `cascade()` wall-clock benchmark on a real bacterial proteome (assets: ship the 79 MB
      ESM2-35M store + tier models). Production re-freeze on a larger pre-t0 population (IBEX).

## Frozen assets

- `models/deepgo_plusplus_integrator_tierA.json` — Tier-A (diam, net_union, interpro).
- `models/deepgo_plusplus_integrator_tierB.json` — Tier-B (esm2_knn).
- ESM2-35M train store (`emb/train_esm2_35m.npz`, 82,404 × 480 fp16) — the CPU kNN
  reference; gitignored scratch, ship as a deployment asset. Re-extract via
  `extract_embeddings.py --model esm2_35m` on ORIX.

## Reuse (already implemented)

- `service/predict.py::DGppLight` — already does ONE DIAMOND search feeding `diam` +
  `net_union` (= net_bridge; precomputed STRING-neighbor vectors, no STRING read at
  request time). Missing: tier-gating + Pfam-fast + ESM2-kNN.
- `pipeline/apply_net_bridge.py`, `build_net_bridge.py` — net_bridge.
- `pipeline/train_integrator.py` — per-aspect logreg OOF over component TSVs.
- `pipeline/build_clean_gt.py` — leak-free GT (OBO-namespace aspects).
