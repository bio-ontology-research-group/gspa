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

## Auxiliary components (eggNOG, ProteInfer, localization, DeepFRI)

Four GO-emitting GSPA predictors that DG++ did not originally stack over, added as
components for complementary evidence. All are computable for a **novel protein from
sequence alone**, and split cleanly by inference hardware:

| component | GSPA wrapper | modality | inference | novel-protein | full DG++ | DG++-Light |
|---|---|---|---|---|---|---|
| `eggnog` | `EggNogMapperPredictor` | orthology (eggNOG OGs) | **CPU** | ✅ OG homolog from FASTA | ✅ | ✅ |
| `proteinfer` | `ProteInferPredictor` | shallow seq-CNN (~50 ms/prot) | **CPU** | ✅ sequence-only | ✅ | ✅ |
| `psortb` | `PSORTbPredictor` | subcellular localization → GO:CC | **CPU** | ✅ sequence-only (bacterial) | ✅ | ✅ |
| `deepfri` | `DeepFriPredictor` | learned-structure GCN | **GPU**¹ | ✅ if structure predicted | ✅ | ✗ |

¹ DeepFRI's GCN runs on CPU *if a structure is supplied*; a novel protein needs a
predicted structure first (ESMFold/AlphaFold = GPU) — the same gate as `foldseek`. Its
`seq` mode is CPU but collapses to a sequence-LM that overlaps `cnn`/`mlp`, so only the
structural (GPU) path is routed in, to the full model.

**Why complementary** (not already in DG++): `eggnog` is *orthology* (curated OGs), distinct
from `diam`'s raw BLAST-KNN — it catches function where similarity is ambiguous, and on
orphans it fires where DIAMOND found no hit. `proteinfer` is a different ab-initio
architecture than `cnn`. `psortb` serves the **CC aspect**, which the sequence/homology
components underserve. `deepfri` is *learned* structure (fires on novel folds), vs
`foldseek`'s structure-*homology* transfer.

### How they slot in (computable for novel proteins)

- **Offline component TSVs** (full DG++ + retraining): `pipeline/build_aux_components.py`
  runs each tool on a query FASTA and emits `<component>.tsv` (`protein⇥GO⇥score`).
  ProteInfer reuses the existing `run_neural_predictors.py --predictor proteinfer` runner.
- **DG++-Light cascade** (CPU three only): `service/predict.py` gains
  `_eggnog_component` / `_proteinfer_component` / `_psortb_component`, wired into
  `cascade()`'s **orphan (Tier B)** block, each gated on its tool being configured
  (constructor params `emapper`/`eggnog_data`, `proteinfer_dir`, `psortb`/`psortb_gram`).
  They are additive: an integrator that doesn't list a component simply ignores it.

### Retrain to weight them (the integrator is component-list-generic)

No integrator code change is needed — only recompute components on the train/test
populations and include their names in `--component-list`:

```bash
# 1. build the new component TSVs on a host with the tools + DBs (eggNOG DB, ProteInfer
#    model, PSORTb, DeepFRI). For the leak-free clean test set and the pre-t0 train set:
build_aux_components.py eggnog  --fasta train.faa --emapper emapper.py --data-dir <db> --out components/eggnog.tsv
build_aux_components.py psortb  --fasta train.faa --gram neg --out components/psortb.tsv
build_aux_components.py deepfri --fasta train.faa --deepfri-dir <repo> --structures <cmaps> --out components/deepfri.tsv
run_neural_predictors.py --predictor proteinfer --model-dir <proteinfer> ... # -> components/proteinfer.tsv

# 2. fold into the tier integrators (Light: + the CPU three; full: + all four incl deepfri)
train_integrator.py --components components --gt gt/gt_no_cleanA_tierB.tsv \
  --component-list esm2_knn,cnn,interpro,eggnog,proteinfer,psortb \
  --save-model models/deepgo_plusplus_integrator_tierB.json ...
```

### Measured standalone ablation (2026-06-23, leak-free clean GT, cafaeval)

Computed on the no-knowledge test set (eggNOG on IBEX KSL DB; DeepFRI seq-mode + PSORTb
via the IBEX FOSS harness; ProteInfer in a TF1.15 container). IA-weighted f_w:

| component (new) | clean-A (MF/BP/CC) | orphan-tier | strict-novel | coverage | note |
|---|---|---|---|---|---|
| **proteinfer** | **0.402** (0.44/0.28/0.48) | **0.368** | 0.393 | 2402/2704 | strongest new signal; beats the whole net family + cnn; CPU |
| **eggnog** | 0.349 (0.43/0.24/0.37) | — | 0.342 | 1816/2704 | net-class; orthology overlaps `diam` (misses diam-orphans) |
| **deepfri** (seq) | 0.287 (0.27/0.21/0.38) | 0.253 | 0.288 | 2704/2704 | modest, full coverage incl. orphans; CPU |
| **psortb** | 0.061 (0.00/0.00/0.18) | 0.008 | 0.062 | 1697/2704 | weak on taxon-mixed set (bacterial-only) → kingdom-gate |

Reference (same GT): PLM heads ~0.47, `lit` 0.459, `diam` 0.451, net family ~0.35,
`cnn` 0.268; orphan tier `esm2_knn` 0.508. **Verdicts:** proteinfer + eggnog are worth
folding into the integrator (proteinfer is the orphan-tier CPU runner-up to esm2_knn at
0.368); deepfri-seq is a marginal full-coverage CPU add (its GPU struct-mode untested);
psortb earns its place only on bacterial genomes. **DeepFRI seq-mode is CPU + full
coverage**, so it is a legitimate Light candidate, not GPU-only as first routed.

### Integrator-inclusion ablation (marginal lift on clean-A, 2026-06-23)

| integrator | mean f_w | Δ vs its base |
|---|---|---|
| FULL base (no aux) | 0.540 | — |
| FULL +aux (proteinfer+eggnog+deepfri) | 0.541 | **+0.001** |
| CPU base (diam,interpro,cnn,net_union,esm2_knn) | 0.517 | — |
| CPU +eggnog+deepfri (no proteinfer) | 0.519 | +0.002 |
| CPU +aux (all 3) | 0.530 | **+0.013** |

**Only `proteinfer` meaningfully helps, and only in the CPU model** (≈ the entire +0.013;
eggnog+deepfri add +0.002 together). `eggnog` is redundant with `diam`/`esm2_knn`
homology; `deepfri` here is **seq-mode (DeepCNN, CPU, no structure)** — another sequence
CNN, redundant with `cnn`/`proteinfer`. In the FULL model the aux add +0.001 because the
GPU PLM heads already capture the sequence signal. **DeepFRI's distinctive structural
signal (GraphConv GCN) needs a predicted structure (GPU) and was NOT benchmarked** — the
seq-mode number above is not the structural DeepFRI. Recommendation: fold `proteinfer`
into the CPU integrator; keep eggnog/deepfri optional; psortb kingdom-gated.

### Struct-mode DeepFRI (GraphConv GCN on ESMFold structures, 2026-06-23)

Tested the *real* structural DeepFRI: ESMFold (ORIX GPU) folded 251 clean-A proteins →
DeepFRI GraphConv GCN (IBEX) → ablate.

| config | mean f_w | vs |
|---|---|---|
| struct-DeepFRI standalone | **0.333** (MF .358 CC .409) | seq-mode DeepCNN 0.287 → **+0.046** |
| CPU base + struct | 0.526 | base 0.517 → **+0.009** |
| FULL base + struct | 0.541 | base 0.540 → **+0.001** |

**Structure is a real signal** (GCN 0.333 > seq-CNN 0.287), but it has no home: it's
**GPU-gated** (novel proteins need a predicted structure → can't run in the CPU/genome
pipeline where it'd add +0.009), and in the GPU FULL model where it *can* run it adds only
+0.001 because **`prostt5` (structure-aware PLM) already encodes the structural signal**.
Efficient structural signal = one ProstT5 forward pass, not ESMFold-per-protein + GCN.

**Status:** engine + builders wired; standalone + integrator-inclusion + struct-DeepFRI
ablations done (clean-A, n=272/251 interim). Lean CPU integrator frozen
(`deepgo_plusplus_integrator_cpu_lean.json`, 0.521 = base+proteinfer). GSPA CLI/sidecar
exposure of the cascade aux components is the remaining plumbing.

## Hierarchy-aware components (C-HMCNN, is_a ∪ part_of) — 2026-06-24

Every **trainable** head is retrained with a hierarchy-aware loss instead of flat BCE:
the **C-HMCNN Max-Constraint Module** (Giunchiglia & Lukasiewicz, NeurIPS 2020) over the
GO **is_a ∪ part_of** DAG. (The kNN/homology and external-tool components — diam, foldseek,
esm2_knn, net, eggnog, proteinfer, psortb, deepfri — are not trainable classifiers; their
hierarchy-awareness is the input max-propagation.) The MCM is a differentiable output layer
that enforces parent ≥ child and delegates each positive to its most-confident true
descendant; implemented with `scatter_reduce(amax)` over the per-aspect descendant edge list
so it scales to GO's thousands of terms (the reference (B,n,n) impl OOMs).

- `pipeline/train_head_hmcnn.py` — PLM heads (`--loss bce|mcm|softreg`, `--save-model` /
  `--load-model` apply). `pipeline/build_cnn_component.py --loss mcm` — the CPU 1D-CNN.
- **Standalone f_w (clean-A):** prostt5 0.474→0.480, esm2_3b 0.475→0.483, ESM2-650M
  0.466→0.474, **cnn 0.257→0.294 (+0.037)** — the weakest head gains most.
- **Integrated (clean-A):** FULL 0.545→**0.550** (prostt5 swap) / all-heads **0.548**;
  CPU `cpu_lean` 0.519→**0.524** (cnn swap). Hard max-constraint > soft penalty; gains in BP/CC.
- **Deployed:** `models/deepgo_plusplus_integrator_{cpu_lean,full_aux}_mcm.json` (committed);
  weights `models/weights/{cnn_mcm,head_prostt5_mcm,head_650m_mcm}.pt` (gitignored, rebuild at
  release). cnn served by `service/predict.py`; PLM heads by the `dgpp-head` sidecar runner.
- **Pending the existing IBEX dep:** tierA/tierB re-freeze on the pre-t0 orphan pop with MCM.

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
