# `deepgo-plusplus` — GSPA's CAFA6 learned-stacker predictor

A CAFA6-competitive GO predictor built **without** new model architecture: it
replaces naive max-merge of GSPA's existing component predictors with a learned
per-aspect logistic-regression stacker. On a faithful CAFA6 reconstruction
(GOA snapshot, t0=2026-02-02, official `cafaeval` + `IA.tsv`) learned integration
reaches **no-knowledge IA-weighted f_w ≈ 0.52** — **level with the 0.524
first-place entry (GOAlpha), not ahead of it.** See
`../benchmark/neural/cafa6_recon/` for the reconstruction harness.

> ## ⚠️ Correction (2026-06-23) — read this before trusting any number below
>
> An audit found **two reporting errors** in earlier versions of this file. The
> corrected numbers are in [Corrected results](#corrected-results-2026-06-23);
> older tables further down were computed on the **contaminated** GT and are kept
> only for the record (marked at the section break).
>
> **Error 1 — metric mismatch (the big one).** The CAFA6 leaderboard headline is
> the **no-knowledge** f_w. Proof: re-scoring our *actual* submission, the
> no-knowledge reconstruction (0.359) matches our real LB (0.37749, Δ0.019); the
> **3-class mean (0.575) does not** (Δ0.198). Earlier text compared our *3-class
> mean* (e.g. 0.647 with `net`) against GOAlpha's *no-knowledge* 0.524 — an
> apples-to-oranges comparison that made integration look ~0.12 ahead of the
> winner. **It is not.** Never compare the 3-class mean to GOAlpha; the 3-class
> mean is inflated by trivially-easy partial-knowledge proteins (0.774).
>
> **Error 2 — leaked no-knowledge MF ground truth (the GAF-date bug).** The GOA
> `date` column is *last-modified*, not *created*, so `build_groundtruth.py`
> re-dates pre-t0 annotations past t0 and counts them as novel. **95% of the
> "no-knowledge MF" targets (2658 → 143) already had that protein's MF term in the
> pre-t0 `train_terms.tsv`.** PPI/homology methods then *retrieve pre-known
> labels*: `net` MF **0.803 → 0.441** once the contamination is removed (its mean
> 0.475 → 0.347). This corroborates GOAlpha's own result that Net-KNN *drops* on
> the genuinely-novel private set (0.29→0.26). Fix + leak-free GT builder:
> `pipeline/build_clean_gt.py` (per-aspect filter on `train_terms`).
>
> **What survives the audit.** The *integrator* is robust to the leak — it blends
> `net` with non-leaking PLM/homology, so on the clean GT it falls only
> **0.532 → 0.521** (and its clean MF 0.651 still beats every single component, a
> genuine integration lift). So the headline "learned integration ≈ GOAlpha-level
> on novel proteins" holds; what does **not** hold is "beats GOAlpha" and the
> standalone `net` MF 0.803. The inflation was ~90% metric-mismatch, ~10%
> net-component leak.

## Corrected results (2026-06-23)

Official `cafaeval`, IA-weighted f_w, **no-knowledge** class, run on `ws`
(driver in `../benchmark/neural/cafa6_recon/`, copied to `~/cafa6_cleaneval/run.sh`).
`dirty` = the GAF-date-contaminated `gt_no.tsv`; `clean-A` = per-aspect
`train_terms` filter (`build_clean_gt.py`); `clean-B` = proteins entirely absent
from `train_terms`.

| predictor | GT | MF | BP | CC | **mean** |
|---|---|---|---|---|---|
| our 2025 submission | dirty | 0.153 | 0.369 | 0.554 | **0.359** ← matches real LB 0.377 |
| our 2025 submission | clean-A | 0.483 | 0.346 | 0.535 | 0.455 |
| `net` (component) | dirty | **0.803** | 0.251 | 0.372 | **0.475** |
| `net` (component) | clean-A | **0.441** | 0.238 | 0.361 | **0.347** ← leak removed |
| `mlp` (component) | dirty | 0.466 | 0.318 | 0.562 | 0.449 |
| `mlp` (component) | clean-A | 0.550 | 0.305 | 0.562 | 0.472 |
| **integrator (6+net)** | dirty | 0.667 | 0.357 | 0.572 | **0.532** |
| **integrator (6+net)** | clean-A | 0.651 | 0.341 | 0.572 | **0.521** |
| **integrator (6+net)** | clean-B | 0.655 | 0.331 | 0.573 | 0.520 |
| GOAlpha (1st place) | CAFA6 no-knowledge LB | — | — | — | **0.524** |

The clean set *helps* non-leaking methods (submission/mlp MF rise) and *hurts*
the one that leaks (`net` MF collapses) — the diagnostic signature of the leak.
**Caveat:** the integrator's clean numbers score a GroupKFold model that was still
*trained* on the contaminated labels; full rigor would retrain on the clean GT,
but scoring-on-clean already shows the headline is not leak-driven. The clean set
is small (282 proteins / 153 MF targets) so exact values are noisy (±0.02); the
directional findings are robust.

### Corrected DeepGO-PlusPlus-Light (no-GPU) — the "beats GPU" claim **reverses**

The Light panels are `net`/`net_union`-centric, so they are the **most**
leak-inflated. On the leak-free GT the old headline ("no-GPU panels beat the full
GPU model on novel proteins") **flips**:

**Shipped Light panels** (OOF integrator, no-knowledge f_w):

| model | components | dirty mean | **clean mean** | dirty MF | clean MF |
|---|---|---|---|---|---|
| full GPU integrator (6+net) | — | 0.532 | **0.521** | 0.667 | 0.651 |
| `light_afdb` = `deepgo_plusplus_light.json` | diam,foldseek,interpro,net | 0.550 | **0.488** | 0.785 | 0.637 |
| `light_cnn` = `..._light_cnn.json` | cnn,diam,foldseek,interpro,net | 0.516 | **0.470** | 0.686 | 0.585 |
| **`light_cpu`** = `..._light_cpu.json` *(was "recommended")* | diam,interpro,net_union | **0.564** | **0.464** | 0.838 | 0.577 |
| `light_fast` = `..._light_fast.json` *(webservice default)* | diam,net_union | 0.562 | **0.451** | 0.845 | 0.550 |

**Full standalone-component ablation** (every component scored directly with
`cafaeval`; no-knowledge f_w, dirty → clean). This is the table the audit was
missing — it includes the **CNN** and every other CPU component:

| component | type | dirty | **clean** | Δ (clean−dirty) |
|---|---|---|---|---|
| `lit` | CPU (text) | 0.445 | **0.459** | +0.014 |
| `diam` | CPU (homology) | 0.388 | **0.451** | +0.063 |
| `foldseek` | **GPU-gated** (structure)† | 0.363 | **0.428** | +0.065 |
| `interpro_lr` | CPU (domain LR) | 0.233 | 0.325 | +0.092 |
| `cnn` | CPU (1D-CNN) | 0.209 | **0.268** | +0.059 |
| `interpro` | CPU (domain raw) | 0.164 | 0.215 | +0.051 |
| `clean` | CPU (EC→GO) | 0.007 | 0.116 | +0.109 |
| `net_bridge` | CPU (DIAMOND→STRING bridge; any seq) | 0.504 | **0.358** | **−0.146** |
| `net_union` | CPU (direct + bridge; any seq) | 0.506 | **0.353** | **−0.153** |
| `net` | CPU‡ — **query must be in STRING** | 0.475 | **0.347** | **−0.128** |
| `mlp` | GPU PLM | 0.449 | 0.472 | +0.023 |
| `prostt5` | GPU PLM | 0.435 | 0.472 | +0.037 |
| `esm2_3b` | GPU PLM | 0.430 | 0.471 | +0.041 |

† **`foldseek` is *not* unconditionally CPU.** FoldSeek *search* is CPU, but each
query needs a **3D structure or a ProstT5-derived 3Di sequence** as input. That is a
CPU lookup only for proteins already in AlphaFold DB (or when structures/3Di are
supplied); for a genuinely novel sequence you must fold/encode it first (ESMFold or
ProstT5 3Di) — **GPU**. So `foldseek` counts as GPU-free *only given structures/3Di*,
not for arbitrary novel sequences.

#### Leak-free net-bridge accuracy (the actual number)

`clean-A` can still let the bridge **self-match** (a `clean-A` protein may remain in
`train_terms`/the DIAMOND DB via another aspect → its query self-matches → it votes
its *own* precomputed STRING-neighbour vector, i.e. direct net, not a real homology
hop). The fully leak-free set is **`clean-B`: proteins entirely absent from
`train_terms` → absent from the bridge's DIAMOND DB → self-match impossible**
(verified: 0/251 clean-B proteins are in `train_terms`). Split by STRING membership:

| component | all clean-B (251) | STRING-member (214) | **non-STRING (37)** |
|---|---|---|---|
| `net` (direct) | 0.339 | 0.380 | **0.000** |
| `net_union` (direct+bridge) | 0.351 | 0.379 | **0.252** |
| `net_bridge` (pure homology hop) | 0.353 | 0.377 | **0.252** |
| `diam` (reference) | 0.442 | 0.414 | **0.598** |

Self-match makes ~no difference (clean-A ≈ clean-B for all net variants), so **~0.35
is the genuine leak-free accuracy** of the bridge. Two findings:

- On non-STRING proteins direct `net` is **0.000** (confirmed N/A), and the bridge
  recovers a real **0.252** — so the bridge is *not* leak/zero; it genuinely adds a
  signal where direct net is blind.
- **But on those same non-STRING proteins plain `diam` scores 0.598 — it crushes the
  bridge (0.252).** The bridge DIAMOND-searches to a homolog and then votes that
  homolog's *PPI-neighbours'* labels; just transferring the homolog's *own* labels
  (`diam`) is far better. So the net/PPI family is **redundant with `diam` and
  dominated by it even on its one supposed unique use case** (non-STRING novel
  proteins). (n=37 non-STRING is small, ±noisy, but the gap is large and the
  direct-net=0 result is exact.)

Corrected conclusions for DG++-Light:

1. **The `net` family is the *only* group that drops on the clean GT** (net −0.13,
   net_union −0.15, net_bridge −0.15); **every other component rises**. That is the
   leak, isolated at component level: the PPI components were retrieving pre-known MF
   labels. Leak-free the whole family sits at ~0.35 — and (table above) it is
   **redundant with / dominated by `diam`**, including on non-STRING proteins. So the
   bridge is real but adds little `diam` doesn't already cover better.
2. **On genuinely-novel proteins the strongest signals are `lit` (0.459), `diam`
   (0.451) and `foldseek` (0.428)† — *not* `net` (0.347).** DG++-Light is built
   around `net_union`, i.e. the most leak-inflated and (on clean) one of the *weakest*
   CPU signals. The **strictly-no-GPU, any-sequence** core should centre on
   `diam`+`lit` (both pure CPU); **`foldseek` adds the next-best signal but is
   GPU-gated** (needs a structure/3Di — CPU only for AFDB-covered proteins, †above),
   so it belongs in a "no-GPU *given structures*" tier, not the strictly-CPU core.
3. **The CNN is the weakest sequence signal** (clean 0.268) but it does **not** leak
   (rises +0.06 on clean); its role stays *coverage* of true orphans, not accuracy.
4. **The full GPU model wins on genuinely-novel proteins** (0.521 vs the Light panels'
   0.45–0.49) — the PLM heads are *not* redundant; that earlier finding was the leak.
   The "recommended" `light_cpu` (0.464) and the webservice default `light_fast` (0.451)
   are the **weakest** panels on clean; `light_afdb` (0.488) is best.

DG++-Light remains a legitimate strictly-CPU predictor, but it is an honest
**no-GPU / fast-deploy option that costs ~0.03–0.07 f_w**, not a win. At the next
clean re-freeze it should be **redesigned around `diam` + `lit`** as the strictly-CPU,
any-novel-sequence core, with `cnn` for orphan coverage and `foldseek` only in a
"no-GPU *given structures/3Di*" tier. **The net/PPI family (`net`/`net_union`/
`net_bridge`) is a candidate to drop entirely**: leak-free it is ~0.35 and `diam`
dominates it everywhere — including on non-STRING proteins (diam 0.598 vs bridge
0.252), the one case the bridge was added for. **Two query-input gates to keep in mind
(both fail exactly on novel proteins):** `foldseek` needs a structure/3Di (GPU to
fold/encode a new sequence — †above), and **plain `net` needs the query to be a
STRING node** (novel proteins are *not* in STRING, so plain `net` is N/A for them, and
the CPU DIAMOND bridge that does fire is itself beaten by `diam` — ‡below). The
CLI/webservice wiring is unaffected (component list lives in the JSON); only the *reported numbers,
the default model, and the recommended panel* change. (Caveat: the clean set is
small, n≈282 / 153 MF targets, so panel values are ±0.02 noisy; the
net-drops/everything-else-rises direction is robust.)

‡ **`net` vs the bridge — what works for a novel protein.** Plain `net` (STRING
Net-KNN) votes the *query's own* STRING neighbours, so it only fires when the query
is itself a STRING node. The CAFA6 no-knowledge set is ~87 % STRING members (well-
studied organisms), so plain `net` still scores there — but a genuinely novel
sequence (e.g. from a freshly sequenced genome) is **not in STRING at all**, and
plain `net` returns nothing. The deployable net signal for such proteins is the
**homology bridge** (`net_union`/`net_bridge`): a CPU DIAMOND search maps the query
to a pre-t0 STRING-member homolog and votes *its* neighbours' labels. So for novel
proteins the net signal is always a DIAMOND-homology step (CPU), never the direct
STRING lookup — and even then its accuracy here was leak-inflated (clean ≈ 0.35).

### Retrained DG++-Light (net-free, leak-free) — 2026-06-23

Acting on the above: re-trained the Light integrator **without the net family**,
**on the leak-free clean GT** (GroupKFold OOF on `gt_no_cleanA`; both train and eval
leak-free). No-knowledge IA-weighted f_w:

| panel (no net) | components | MF | BP | CC | **mean** |
|---|---|---|---|---|---|
| `fold` → `deepgo_plusplus_light_clean.json` | diam,foldseek,interpro,lit | 0.616 | 0.362 | 0.547 | **0.508** |
| `cpu` → `deepgo_plusplus_light_cpu_clean.json` | diam,interpro,lit,cnn | 0.612 | 0.345 | 0.542 | **0.500** |
| (diam,lit,interpro) | — | 0.598 | 0.346 | 0.538 | 0.494 |
| (diam,lit) | — | 0.603 | 0.330 | 0.535 | 0.489 |
| *old net-based* `light_afdb` | diam,foldseek,interpro,net | — | — | — | 0.488 |
| *old net-based* `light_cpu` | diam,interpro,net_union | — | — | — | 0.464 |
| full GPU integrator (ref) | — | — | — | — | 0.521 |

**Dropping net + retraining leak-free *raised* DG++-Light from 0.464–0.488 to
0.50–0.51** — now within ~0.01 of the full GPU model. The CNN even *helps* marginally
here (no longer redundant with net). Shipped two re-frozen models:
`deepgo_plusplus_light_clean.json` (best, needs structures for foldseek) and
`deepgo_plusplus_light_cpu_clean.json` (strictly-CPU, any sequence). **Caveat:** these
are full-data fits on the small leak-free clean-A (n≈282) — interim v2 integrators;
the production fix is to retrain on a *large* pre-t0 leak-free population (deferred,
IBEX). The fair train-protocol comparison to the 0.521 GPU number also needs that.

### Better CPU predictor — literature-driven roadmap (2026-06-23)

The weak link is the **sequence model for orphans / low homology** (DeepGOCNN-style
1D-CNN, clean 0.268 — independent benchmarks confirm DeepGOCNN is among the weakest
learned methods, Fmax MFO ≈ 0.37, vs PLM-embedding heads ≈ 0.45–0.62). A literature
sweep (ProteInfer eLife 2023; goPredSim Sci Rep 2021; Light-Attention 2021; ESM2
small models; VESM distillation Nat Methods 2026; DL-AFP benchmark Brief Bioinform
2024) → ranked plan for a **CPU-inference** upgrade (training may use GPU):

1. **Embedding-kNN third channel (goPredSim-style) — best bet, low effort, no training.**
   Precompute **ESM2-35M** (CPU-deployable; or ProtT5, cached) embeddings for the pre-t0
   train set; at query time embed (ESM2-35M on CPU, sub-second) + cosine-kNN GO transfer.
   Fills exactly the BLAST-blind <20 %-identity gap where `diam` collapses to ~0.10 and
   PLM-kNN still predicts. **GSPA already ships `Esm2CentroidPredictor` (embedding-kNN)
   and a `ProteInferPredictor`** — we have the wrappers.
2. **Replace the 1D-CNN with a small-PLM embedding + light-attention head.** ESM2-8M/35M
   frozen embeddings → light-attention pooling (Stärk 2021) → DAG-aware sigmoid head.
   PLM-fed heads beat DeepGOCNN by ~0.08 MFO in independent benchmarks; light-attention
   is a near-free pooling upgrade over global-max-pool. Train on GPU (ORIX), infer CPU.
3. **Swap DeepGOCNN convs → ProteInfer dilated-residual CNN** (5 residual dilated blocks;
   5–7 MB model, proven <1.5 s CPU/browser inference; clustered-split Fmax ≫ DeepGOCNN).
4. **Distill ESM2-650M/3B or ProtT5 (teacher) → ESM2-35M / ProteInfer student** for GO
   (VESM-validated principle): highest ceiling, GPU-train / CPU-infer. (We have ESM2-3B +
   ProstT5 on ORIX already from Phase B.)
5. **Loss/head upgrades** (stack on any of the above): IA-weighted **asymmetric loss**,
   **hierarchical-violation penalty** (child ≤ parent, TALE/DeepGOZero), label smoothing
   + temperature calibration — cheap, target IA-weighted F-max + novel-protein calibration.

**Single best next experiment:** ESM2-35M embedding-kNN channel (#1) — CPU-only at
inference, no training, attacks the orphan regime; extract train embeddings on ORIX
(minutes, GPU-OK), add as a component to the leak-free integrator, re-eval. Sources in
`benchmark/neural/cafa6_recon/` notes; see the ranked table for effort/citations.

---

> **Everything below this line predates the 2026-06-23 audit.** The no-knowledge
> numbers (esp. MF, `net`, and any "3-class mean vs GOAlpha" comparison) are on the
> **contaminated** GT — read them as *contaminated-GT* values, superseded by the
> corrected table above.

## Components (produced upstream)

Each is a TSV `protein\tterm\tscore` named `<component>.tsv[.gz]` in a
components directory. The shipped model uses six:

| component | source | role |
|---|---|---|
| `diam`     | `DiamondPredictor`      | BLAST-KNN homology transfer |
| `foldseek` | `FoldSeekPredictor`     | FoldSeek-KNN structure transfer |
| `clean`    | `CleanPredictor`        | EC → GO |
| `interpro` | `InterProScanPredictor` | domain → GO |
| `mlp`      | ESM2 MLP head           | sequence-PLM |
| `prostt5`  | ProstT5 MLP head        | **structure-aware PLM** (the complementary lever) |
| `net`      | Net-KNN over STRING PPI | **guilt-by-association** (the network signal — *recommended*) |
| `lit`      | BM25 text-kNN           | literature/text (optional; see caveat) |

Per-component no-knowledge f_w: net 0.475 (**best single-component MF 0.803**),
mlp 0.447, lit 0.445 (MF 0.539), prostt5 0.435 (MF 0.466), diam 0.389,
foldseek 0.364, interpro 0.165, clean ~0. ProstT5 is the only PLM that
*complements* the MLP rather than cannibalising it (ESM2-3B regresses —
correlated-PLM collinearity).

### Network + literature components (added; see §3 for the build)

Two extra signals reproduce the rest of GOAlpha's heterogeneous panel:

* **`net` — Net-KNN over STRING.** For each query protein, vote its STRING v12
  neighbours' pre-t0 GO labels, weighted by the STRING combined score.
  **Leak-free**: STRING v12 (2023) edges and `train_terms` labels are both
  before t0 = 2026-02-02, and the query's own annotations are never used. The
  very high no-knowledge MF (0.803) is genuine guilt-by-association — novel
  proteins that are members of already-characterised complexes — and is
  concentrated in well-studied organisms with dense PPI (the test set is ~half
  human).
* **`lit` — BM25 text-kNN.** Transfer GO from textually-similar training
  proteins (GORetriever in spirit). **Caveat:** the query text is restricted to
  identification fields (protein name + gene name), never the post-t0
  `CC FUNCTION`, but a novel protein's *name* in the current SwissProt dump may
  still have been updated post-t0 — so its no-knowledge MF gain (0.539) carries
  a residual leakage risk. Treat `lit` as an optional no-knowledge booster, not
  a default.

### Results — adding net and lit (IA-weighted f_w; baseline = the 6-comp above)

| model | no-knowledge | limited | partial | **official 3-class mean** |
|---|---|---|---|---|
| 6-comp (baseline)        | 0.489 | 0.630 | 0.768 | 0.629 |
| **6-comp + net** *(ship)*| 0.538 | 0.654 | 0.748 | **0.647** (+0.018) |
| 6-comp + lit + net       | 0.553 | 0.646 | 0.715 | 0.638 |
| 6-comp + lit             | 0.507 | 0.631 | 0.720 | 0.620 |

**`net` is the clean win**: +0.018 official, gains on novel/limited proteins,
tiny partial dip, no leakage — shipped as
`models/deepgo_plusplus_integrator_net.json`. Adding `lit` pushes no-knowledge
higher (the real-LB proxy) but *drags* the 3-class mean by hurting
partial-knowledge proteins (the integrator is frozen on no-knowledge weights),
so the 8-comp model is offered separately as
`models/deepgo_plusplus_integrator_lit_net.json` for no-knowledge-focused use.
The partial-knowledge dip is the frozen-on-no-knowledge weighting artifact;
per-knowledge-class / pre-t0-population training (deferred) would remove it.

## Full ablation — which component / aggregator actually helps

This is the CAFA6-winner-style ablation: quantify the *marginal* contribution of
every component and compare the **logistic** aggregator against **gradient-
boosted trees** (GOAlpha's Learning-to-Rank choice). All rows are scored
identically — per-aspect candidate rows → GroupKFold-by-protein **out-of-fold**
scores → official `cafaeval` IA-weighted f_w on the **no-knowledge** GT (2,694
novel proteins; the metric the project optimises, where CAFA6 is decided).

Reproduce the whole table in one pass (loads each component once, scores all
configs with a single `cafaeval` call):

```bash
python pipeline/ablation.py --components <components-dir> --gt gt/gt_no.tsv \
    --dag ~/Public/software/cafa6/go-dag.tsv \
    --ia  ~/Public/software/cafa6/kaggle-official/IA.tsv \
    --obo ~/Public/software/cafa6/go.obo \
    --train-terms ~/Public/software/cafa6/kaggle-official/train_terms.tsv \
    --out ablation_no
```

Raw numbers committed to `ablation_no_results.tsv`. (CPU-heavy — run on a many-
core box, not a laptop: ~30 min here.)

### 1. Single components (each alone, no-knowledge f_w)

| component | MF | BP | CC | **mean** |
|---|---|---|---|---|
| `net`        | **0.803** | 0.251 | 0.370 | **0.475** |
| `mlp`        | 0.466 | 0.323 | 0.557 | **0.449** |
| `lit`        | 0.540 | 0.287 | 0.507 | **0.445** |
| `prostt5`    | 0.466 | 0.278 | 0.562 | **0.435** |
| `esm2_3b`    | 0.422 | 0.318 | 0.547 | **0.429** |
| `diam`       | 0.323 | 0.348 | 0.496 | **0.389** |
| `foldseek`   | 0.339 | 0.314 | 0.442 | **0.365** |
| `interpro_lr`| 0.099 | 0.287 | 0.335 | **0.240** |
| `interpro`   | 0.243 | 0.206 | 0.061 | **0.170** |
| `clean`      | 0.020 | 0.047 | 0.000 | **0.022** |

`net`'s MF 0.803 is genuine guilt-by-association (novel members of characterised
complexes); `clean` (EC→GO) is ~noise on novel proteins.

### 2. Marginal contribution — leave-one-out from the 7-comp panel

Panel = `diam,foldseek,clean,interpro,mlp,prostt5,net`; OOF mean **0.532**. Δ is
*panel mean − (panel minus that component)*, i.e. what the component is worth **on
top of all the others**.

| removed | MF | BP | CC | mean | **marginal Δ** |
|---|---|---|---|---|---|
| *(full 7-comp panel)* | 0.667 | 0.357 | 0.572 | 0.532 | — |
| − `net`      | 0.518 | 0.351 | 0.581 | 0.483 | **+0.049** |
| − `interpro` | 0.663 | 0.350 | 0.569 | 0.527 | **+0.005** |
| − `diam`     | 0.664 | 0.351 | 0.577 | 0.531 | +0.001 |
| − `foldseek` | 0.669 | 0.358 | 0.572 | 0.533 | −0.001 |
| − `mlp`      | 0.678 | 0.349 | 0.577 | 0.535 | −0.003 |
| − `clean`    | 0.666 | 0.366 | 0.572 | 0.535 | **−0.003** |
| − `prostt5`  | 0.711 | 0.387 | 0.536 | 0.545 | **−0.013** |

**`net` is the only component with a substantial positive marginal contribution
(+0.049).** Everything else lands within ±0.013 — and `clean`, `foldseek`, `mlp`
are at or *below* zero: removing them *helps or is neutral* once `net` is present.
`prostt5` is the strongest standalone-PLM lever but its no-knowledge **mean** goes
*negative* (−0.013) inside the full panel — `net` subsumes most of its MF signal,
and it only still earns its place on **CC** (0.572 → 0.536 when removed) and on the
mixed-knowledge (partial) eval. This is diminishing-returns / redundancy: six PLM
+ homology + domain channels carry largely the same information for novel proteins.

### 3. Cumulative build-up (best-first)

| running panel | mean | Δ |
|---|---|---|
| `mlp`                                   | 0.449 | — |
| `+prostt5`                              | 0.470 | +0.021 |
| `+diam`                                 | 0.480 | +0.010 |
| `+foldseek`                            | 0.483 | +0.003 |
| `+interpro`                            | 0.484 | +0.002 |
| `+clean`                               | 0.484 | **−0.001** |
| `+net`                                  | 0.532 | **+0.048** |

Each PLM/homology/domain channel adds a shrinking sliver; **`clean` is the one
component that *decreases* the running score**, and `net` is the single biggest
jump. The honest leak-free models all cluster in **0.48–0.53** — that tight spread
is itself the diagnosis: with these features there is no more juice to squeeze
without leaking (next row).

### 4. Aggregator — logistic regression vs gradient-boosted trees

Same 6-comp panel (`diam,foldseek,clean,interpro,mlp,prostt5`), same OOF protocol;
only the integrator model / feature set changes.

| aggregator | features | MF | BP | CC | **mean** |
|---|---|---|---|---|---|
| **logreg** *(shipped)* | 6 component scores | 0.518 | 0.351 | 0.581 | **0.483** |
| **xgboost** | 6 component scores | 0.600 | 0.328 | 0.505 | **0.478** |
| logreg | scores + IA + log-freq | 0.575 | 0.351 | 0.581 | **0.502** |
| xgboost | scores + IA + log-freq | **0.963** | 0.339 | 0.582 | **0.628** |

**Leak-free gradient boosting does *not* beat the linear stacker** (0.478 vs
0.483): with only six calibrated score features, XGBoost's extra capacity overfits
slightly — it lifts MF but loses more on BP/CC. The winning entry's XGBoost-LTR
edge came from a *rich heterogeneous* feature set + a test-matched validation
protocol, **not** tree capacity per se; on our score-only features the linear model
is the right, leak-resistant choice.

**Adding term-identity features (IA + log-freq) is where trees go off the rails.**
For the *linear* model they act as mild honest per-term priors (+0.019, MF 0.575).
For **XGBoost they leak**: MF rockets to **0.963** and the mean to 0.628. The OOF
folds split by *protein*, but **GO terms are shared across folds**, so a tree can
memorise *which terms are frequently positive among no-knowledge proteins* — that
is label leakage, not generalisation, and 0.963 MF is non-credible (the winner's
*entire* novel-protein score is ~0.5). This is exactly why DeepGO-PlusPlus ships
the **linear, scores-only** integrator and why an earlier XGBoost+IA/freq variant
was discarded. The takeaway mirrors the CAFA6 post-mortem: **the bottleneck is
generalisation and signal diversity (one real new channel, `net`), not aggregator
sophistication** — the only way to "beat" 0.53 on this feature set is to leak.

## DeepGO-PlusPlus-Light — a no-GPU predictor

The three strongest standalone components are GPU PLM heads (`mlp`, `prostt5`,
`esm2_3b`). DG++-Light targets deployment **without a GPU**, so it drops all of
them and relies on the CPU components — DIAMOND BLAST-KNN (`diam`), FoldSeek
structure-KNN (`foldseek`, over precomputed AlphaFold-DB structures), InterProScan
domains raw (`interpro`) and as a logistic model (`interpro_lr`), STRING Net-KNN
(`net`), optional BM25 `lit`, and a CPU **1D-CNN over sequence** (`cnn`, DeepGOCNN
style — `pipeline/build_cnn_component.py`) that *replaces* the PLM heads. Reproduce
the panel comparison with `pipeline/eval_light.py`.

### Dropping the PLM heads (no CNN yet) — no-knowledge f_w

| no-GPU panel | MF | BP | CC | **mean** |
|---|---|---|---|---|
| `diam+foldseek+interpro+net` | 0.785 | 0.366 | 0.498 | **0.550** |
| `diam+net` | 0.787 | 0.349 | 0.490 | **0.542** |
| `diam+foldseek+interpro_lr+net+lit` | 0.700 | 0.368 | 0.519 | **0.529** |
| `diam+foldseek+interpro+interpro_lr+net+lit` | 0.702 | 0.362 | 0.512 | **0.525** |
| `diam+interpro_lr+net` | 0.723 | 0.360 | 0.486 | **0.523** |
| `diam+foldseek+interpro+interpro_lr+net` | 0.726 | 0.362 | 0.476 | **0.521** |
| `diam+foldseek+interpro_lr+net` | 0.728 | 0.355 | 0.472 | **0.518** |

**The headline is counter-intuitive but follows directly from the ablation: on
novel proteins the best no-GPU panel (`diam+foldseek+interpro+net`, 0.550) *beats*
the full GPU 7-component model (0.532).** The PLM heads are redundant with `net` on
no-knowledge MF (the leave-one-out showed `−prostt5` *+0.013* and `−mlp` *+0.003*),
so for novel proteins removing them and trusting `net` + homology + structure +
domains is a net win. (The PLM heads still earn their place on the easier
limited/partial knowledge classes and on CC; DG++-Light is a *novel-protein-first*,
GPU-free model — for proteins with prior annotations the full model is still
preferred.) `interpro` (raw) beats `interpro_lr` inside the panel, and `lit` adds
CC but costs MF, so it is left optional.

### Adding the CPU 1D-CNN (`cnn`)

The 1D-CNN (embedding → parallel Conv1d k∈{8,16,24,32} → global-max-pool → dense →
5,265-term head; frequency-weighted BCE; trained on the 80,750 pre-t0 SwissProt
proteins, CPU-only, ~7 min/epoch on 56 cores, early-stopped at epoch 7) is the
sequence model meant to *replace* the PLM heads. Result — no-knowledge f_w:

| panel | MF | BP | CC | **mean** |
|---|---|---|---|---|
| `cnn` (standalone) | 0.088 | 0.186 | 0.345 | **0.206** |
| `diam+foldseek+interpro_lr+net` | 0.728 | 0.355 | 0.472 | **0.518** |
| ` + cnn` | 0.672 | 0.355 | 0.491 | **0.506** *(−0.012)* |
| `diam+foldseek+interpro+interpro_lr+net` | 0.726 | 0.362 | 0.476 | **0.521** |
| ` + cnn` | 0.672 | 0.365 | 0.488 | **0.508** *(−0.013)* |
| `cnn+diam+foldseek+interpro+net` | 0.686 | 0.369 | 0.493 | **0.516** |

**The CNN does *not* help on novel proteins — it consistently *lowers* the panel
by ~0.012.** Two reasons, and they are exactly the CAFA6 lesson: (1) novel
proteins have no close homologs, but a CNN trained on pre-t0 sequences *also*
fails to generalise to genuinely new functions (standalone 0.206, MF 0.088 — the
generalisation wall that homology/structure/network climb via guilt-by-
association). This is the same redundancy the main ablation found for the *GPU*
PLM heads — for novel proteins, learned sequence models add little over `net`.
(2) Its dense output (~2.5k terms/protein above 0.01) pollutes the integrator's
candidate union, costing precision.

**Where the CNN *does* earn its place: coverage.** It predicts for **all 2,694**
no-knowledge proteins, vs `diam` 2,253 / `net` 2,211 / `interpro` 1,390 — so for
true orphans (no homolog, no structure, no PPI, no domain) it is the *only* signal.

### Shipped DG++-Light models

| model | components | no-knowledge mean f_w | use |
|---|---|---|---|
| **`deepgo_plusplus_light_cpu.json`** | `diam,interpro,net_union` | **0.564** | **recommended** — strictly no-GPU (no structures needed) *and* works for proteins not in STRING (DIAMOND-bridged `net`) |
| `deepgo_plusplus_light.json` | `diam,foldseek,interpro,net` | 0.550 | no-GPU *given AFDB structures* for foldseek |
| `deepgo_plusplus_light_cnn.json` | `cnn,diam,foldseek,interpro,net` | 0.516 | coverage-first: adds the CPU 1D-CNN for orphan proteins (lower mean, full coverage) |

All three beat the full GPU model (0.532) on novel proteins.

Both run through the unchanged `DeepGoPlusPlusPredictor` (component list read from
the JSON). The `cnn` component is rebuilt at each release with
`pipeline/build_cnn_component.py` (train, or `--save-model`/`--load-model` to apply
a saved checkpoint without retraining); the other CPU components come from GSPA's
DIAMOND/FoldSeek/InterProScan wrappers and `make net`.

### How "no-GPU" really is each component (FoldSeek + STRING caveats)

Two components carry a hidden dependency worth stating plainly:

* **`foldseek` is CPU only if you already have query structures.** FoldSeek
  *search* is CPU, but you must give it a structure per query. With a precomputed
  structure (AlphaFold DB now covers ~all of UniProt, incl. the SwissProt CAFA6
  test set) it is a CPU lookup — no GPU. For a genuinely novel sequence **not in
  AFDB**, you must fold it first (ESMFold or FoldSeek's ProstT5 3Di) — that needs
  a **GPU**. So `foldseek` is GPU-free *given structures*, not unconditionally.
* **`net` only fires for proteins that are STRING nodes** (see the bridge below).

So there are two GPU-free tiers (no-knowledge mean f_w):

| tier | panel | f_w | GPU-free for… |
|---|---|---|---|
| with structures | `diam+foldseek+interpro+net` | **0.550** | any UniProt/AFDB-covered protein |
| **structure-free** | `diam+interpro+net` | **0.544** | **any sequence** (no folding at all) |

Dropping `foldseek` costs only **0.006** — the **structure-free** panel still
**beats the full GPU model (0.532)** on novel proteins. So a strictly-no-GPU,
any-input DG++-Light is `diam+interpro+net` (+ the CPU `cnn`/bridge below).

### Homology-bridged `net` — guilt-by-association for proteins not in STRING

Plain `net` is 0 for a query that is not itself a STRING node — exactly the novel
case. **356 of the 2,694** no-knowledge proteins have no STRING id. The bridge
(`pipeline/build_net_bridge.py`): DIAMOND the query against the **pre-t0 train**
proteins → take its STRING-member homolog(s) `h` → vote `h`'s STRING neighbours'
**pre-t0** GO labels, weighted by homology strength × STRING confidence.

**Leak-safe by construction** (the train/test split that matters here): targets
are pre-t0 train proteins (STRING v12 = 2023; labels = pre-t0 `train_terms`);
no-knowledge queries have *no* pre-t0 labels so they are **absent from the train
DB → DIAMOND can never self-match**; and `h ≠ q`, neighbour `≠ q`, q's own
(post-t0) truth is never read. The query is purely on the test side; everything
voted is training/pre-t0 side.

Result on the **356 no-STRING proteins** (isolated; this is the realistic
"novel, not in STRING" test):

| component on the 356 | MF | BP | CC | mean |
|---|---|---|---|---|
| plain `net` | 0.000 | 0.000 | 0.000 | **0.000** (no STRING node → nothing) |
| **`net_bridge`** | 0.684 | 0.176 | 0.405 | **0.422** |

The DIAMOND→STRING hop recovers strong guilt-by-association (MF 0.68) for proteins
plain STRING `net` cannot touch at all — covering **259 of the 356**.

**Controlled hold-out benchmark (the rigorous, larger-scale test).** Set A (356
genuinely-no-STRING proteins) is small and biased toward orphans. The clean test
is to take proteins that *are* in STRING, **delete their STRING node**, and see how
much of the direct-`net` signal the DIAMOND bridge recovers. `net_bridge` already
routes every query `q → homolog h → N(h)` and never uses `q`'s own edges; we
verified **0 of 2,014** no-knowledge set-B queries have a homolog mapping to `q`'s
own STRING id (0.0% contamination), so it *is* a clean node-removal simulation.

**Swept across all three knowledge classes** (every in-STRING test protein, node
held out; mean IA-weighted f_w; `cafaeval` run on the office box):

| knowledge class (in-STRING, node held out) | n | `net` direct (**oracle**) | `net_bridge` (node **removed**) | recovery |
|---|---|---|---|---|
| no-knowledge | 2,348  | 0.5213 | 0.5210 | **~100 %** |
| limited      | 7,319  | 0.5227 | 0.5217 | **~100 %** |
| partial      | 16,083 | 0.4490 | 0.4320 | **~96 %** |

**The DIAMOND bridge recovers essentially all of the direct-STRING signal with the
protein's own node deleted** — ~100 % for novel and limited proteins, ~96 % for
partial-knowledge (where the real STRING neighbourhood is more informative, so the
homology hop costs a little more: 0.449 → 0.432). Per-aspect, the hop trades a
small MF/CC loss for a BP gain. So a protein absent from STRING loses almost
nothing provided it has a present homolog — across the whole test set the bridge is
a near-perfect *substitute* for STRING membership, not just a fallback. (Caveat: it
still needs *some* STRING-member homolog; true no-homolog orphans remain the
`cnn`/coverage case.)

As a drop-in `net` replacement (`net_union` = direct where the protein is a STRING node, bridge
otherwise), in the full structure-free panel (no-knowledge mean f_w):

| structure-free panel | MF | BP | CC | **mean** |
|---|---|---|---|---|
| `diam+interpro+net`        | 0.784 | 0.352 | 0.495 | **0.544** |
| `diam+interpro+net_bridge` (pure homolog-hop) | 0.766 | 0.331 | 0.506 | **0.534** |
| **`diam+interpro+net_union`** (direct + bridge) | 0.838 | 0.358 | 0.496 | **0.564** |

**The bridge lifts the whole structure-free panel +0.020 (0.544 → 0.564)** — by
giving the 356 no-STRING proteins net's strong MF signal (panel MF 0.784 → 0.838),
it beats *both* the foldseek panel (0.550) and the full GPU model (0.532). Shipped
as **`deepgo_plusplus_light_cpu.json`** (`diam,interpro,net_union`) — the strictly
no-GPU, any-sequence model. `net_union` is rebuilt with `build_net_bridge.py`
(DIAMOND vs the pre-t0 train DB) + plain `net`. (`net_bridge` alone is below direct
`net` because the homology hop adds noise for proteins that *are* STRING nodes —
hence the union, not the pure bridge.)

## 1. Train + freeze the integrator (once)

```bash
python pipeline/train_integrator.py \
  --components ../benchmark/neural/cafa6_recon/components --gt ../benchmark/neural/cafa6_recon/gt/gt_no.tsv \
  --dag   ~/Public/software/cafa6/go-dag.tsv \
  --ia    ~/Public/software/cafa6/kaggle-official/IA.tsv \
  --train-terms ~/Public/software/cafa6/kaggle-official/train_terms.tsv \
  --taxon ~/Public/software/cafa6/kaggle-official/testsuperset-taxon-list.tsv \
  --model logreg --features scores \
  --component-list diam,foldseek,clean,interpro,mlp,prostt5 \
  --save-model models/deepgo_plusplus_integrator.json --out /tmp/ltr_run
```

`--model logreg --features scores` is **leak-free** — the linear model over
component scores cannot memorise term identity (an XGBoost variant with IA/freq
features leaked to MF f_w 0.963 in-distribution and was discarded). The OOF
GroupKFold-by-protein score in `/tmp/ltr_run/preds/ltr.tsv` is the honest
generalisation estimate (0.483). The frozen JSON is the full-data fit
(`models/deepgo_plusplus_integrator.json`, 2 KB) shipped for inference.

## 2. Apply — directly

```bash
python run_neural_predictors.py --predictor deepgo-plusplus \
  --manifest manifest.tsv --min-score 0.1 \
  --integrator models/deepgo_plusplus_integrator.json \
  --components-dir ../benchmark/neural/cafa6_recon/components \
  --dag ~/Public/software/cafa6/go-dag.tsv
```

Emits `<tag>.deepgo-plusplus.tsv` (`protein\tterm\tscore\tGO`). Applying the
frozen model to the no-knowledge set scores **0.489** (in-sample full-fit;
consistent with the 0.483 OOF), confirming the artifact is faithful.

## 2′. Apply — via GSPA

```bash
./gradlew :gspa-cli:run --args="annotate --input proteome.faa --output out \
  --neural-sidecar $PWD/benchmark/neural/run_neural_predictors.py \
  --deepgo-plusplus \
  --deepgo-plusplus-integrator $PWD/models/deepgo_plusplus_integrator.json \
  --deepgo-plusplus-components-dir <dir-of-component-scores> \
  --deepgo-plusplus-dag ~/Public/software/cafa6/go-dag.tsv"
```

(or set `predictors.neural.deepGoPlusPlus.{enabled,integrator,componentsDir,dag}`
in YAML). The sidecar propagates each component to GO ancestors (max), forms
per-(protein, term) candidates per aspect, and emits `sigmoid(w·x + b)`.

## 3. Building the `net` and `lit` components

Both derive from one streaming pass over the SwissProt flat file, which yields
per-accession identification text, the `DR STRING` xref, and the `OX` taxon:

```bash
# one parse -> text_string_index.tsv (accession\ttaxon\tstring_id\tname\tchar_text)
python pipeline/build_text_string_index.py uniprot_sprot.dat.gz text_string_index.tsv
```

**Literature** (`lit`) — BM25 text-kNN, CPU-only, `--shard i/N` for parallelism:
```bash
python pipeline/build_lit_component.py --index text_string_index.tsv \
  --train-terms train_terms.tsv --queries test_proteins.txt \
  --out lit.tsv --topk 30          # corpus uses full text; query uses NAME only
```

**Net-KNN** (`net`) — needs STRING per-species link files. The test set spans
only ~68 species (97 % of proteins have a STRING id), so download per species,
not the full dump. `../benchmark/neural/run_net_ws.sh` does the download (with
gzip integrity + retry, corrupt files dropped) and the build in one step:
```bash
# slim 3-column index is enough for net (no text):
cut -f1-3 text_string_index.tsv > net_index.tsv
python pipeline/build_net_component.py --index net_index.tsv \
  --train-terms train_terms.tsv --queries test_proteins.txt \
  --string-dir <dir of {taxid}.protein.links.v12.0.txt.gz> \
  --out net.tsv --min-conf 400 --topk 50
```

Gzip the outputs into the components dir (`components/{net,lit}.tsv.gz`) and add
`net` (and optionally `lit`) to `--component-list`. The Groovy
`DeepGoPlusPlusPredictor` needs **no change** — it reads the component list from
the integrator JSON, so swapping in `deepgo_plusplus_integrator_net.json` is
enough.

## Temporal integrity — blind pre-t0-validation check

The integrator's frozen weights are normally fit (GroupKFold) on the
no-knowledge **test** labels, which you would not have at submission time. To
quantify that optimism we trained a `{prostt5, esm2_3b, net}` integrator two
ways and applied both to the test set:

| integrator training population | no-know | limited | partial | **official 3-class** |
|---|---|---|---|---|
| **blind** — pre-t0 validation set (25k train proteins, never sees test) | 0.508 | 0.584 | 0.591 | **0.561** |
| test-trained — on no-knowledge test labels | 0.527 | 0.585 | 0.585 | 0.566 |

**Training the integrator blind costs only 0.005 on the official metric** (and
is +0.006 *better* on partial-knowledge). And `net`'s learned weight *increases*
under blind training (MF/BP/CC 0.81/0.78/0.92 blind vs 0.58/0.34/0.37
test-trained) — confirming it is a genuine pre-t0 signal, not a tuning artifact.
So the DeepGO-PlusPlus result is essentially CAFA-faithful; the only non-faithful
element is `lit` (query-name leakage), which is why it is optional. The blind
PLM head scores come from k-fold OOF on the train embeddings
(`train_head_oof.py`, ORIX); the blind model is `../benchmark/neural/cafa6_recon/integrator_pret0.json`.

## Remaining gap to GOAlpha (0.524)

Not models, and not the integrator training population (shown negligible above).
The open lever is a **full-component** blind/pre-t0 integrator: the check above
used the 3 components reproducible for arbitrary train proteins
(`prostt5/esm2_3b/net`); reproducing `diam/foldseek/clean/interpro` for a
validation set (IBEX work) would let the *full* model train on the pre-t0
population and let XGBoost + IA/freq help leak-free. STRING Net-KNN and a
(leakage-clean) literature channel are now done.
