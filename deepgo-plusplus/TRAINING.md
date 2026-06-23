# DeepGO-PlusPlus — training, updating, and data

Everything needed to understand **how the shipped model was trained**, **how to
retrain it at a new UniProt / STRING / CAFA release**, **where every input comes
from and how to (re)download it**, and **how it performs on CAFA6**.

- Quick start / module overview → `README.md`
- Full ablation tables + temporal-integrity deep dive → `RESULTS.md`
- Input release provenance per shipped model → `VERSIONS.md`

---

## 1. How it was trained

DeepGO-PlusPlus is **not a new model architecture**. It is a learned *integrator*
(a frozen per-aspect logistic regression) over the per-component GO scores that
GSPA's existing predictors already produce. The CAFA6 post-mortem found our
original entry failed on **train/test setup and generalisation**, not model
class — and that naive max-merge of components actively *destroyed* the
molecular-function signal (MF f_w 0.464 → 0.153). Replacing max-merge with a
learned stacker is the whole idea.

### 1.0 Faithful CAFA6 reconstruction (the evaluation protocol)

Training and evaluation happen on a faithful offline reconstruction of CAFA6 so
the numbers track the real leaderboard rather than a leaky internal split:

- **t0 = 2026-02-02** (CAFA6 final-submission deadline). Everything used to
  *predict* must predate t0; the test truth is what proteins *gain* after t0.
- `pipeline/build_groundtruth.py` reads a dated GOA filter
  (`testsuperset_exp_annots.tsv`: `accession⇥GO⇥F|P|C⇥evidence⇥YYYYMMDD`) and
  emits the three CAFA6 knowledge-class ground truths:
  - **no-knowledge** — protein had *no* experimental annotation in *any* aspect
    on/before t0 (the novel-protein case; the real public-LB proxy).
  - **limited** — experimental annotation in 1–2 aspects pre-t0.
  - **partial** — all three aspects annotated pre-t0.
  A (protein, aspect) is a test target iff it accumulated an experimental
  annotation in that aspect *after* t0; the ground-truth term set is the
  protein's *complete* experimental set in that aspect (full-evaluation mode;
  `cafaeval` propagates ancestors).
- Scoring is the **official `cafaeval`** (BioComputingUP) with the provided
  `IA.tsv`: IA-weighted max-F (`f_w`), `-norm cafa -prop max`. **The CAFA6
  leaderboard headline is the `no-knowledge` class f_w** (verified: our
  submission's no-knowledge recon 0.359 matches our real LB 0.37749; the 3-class
  mean 0.575 does not). Report and compare on **no-knowledge**, not the 3-class
  mean.

### 1.0a ⚠️ Known bug — GAF last-modified dates contaminate the no-knowledge GT

**Symptom (found 2026-06-23).** 91% of the reconstruction's "no-knowledge"
proteins, and **95% of its "no-knowledge MF" targets (2658 → 143)**, already have
that protein's MF term in the pre-t0 `train_terms.tsv`. So the "novel" targets are
largely *not novel* — their function was already in the training labels. IA-weighted,
**61% of the no-knowledge MF truth is pre-known** (BP 8%, CC 5% — MF-specific).

**Root cause.** `build_groundtruth.py` decides "is this a novel post-t0 target?"
from the GOA `date` column (col 14). That column is the annotation's
**last-modified** date, **not** its creation date. An annotation that existed
before t0 but was touched after t0 (re-curated, re-validated, format-migrated) gets
a post-t0 date and is wrongly counted as a brand-new discovery. Example: `GO:0005515`
on `P76092`/`Q5SQS8`, evidence IPI, date `20260613` (post-t0) — yet the same term
is in the pre-t0 `train_terms`.

**Effect.** Homology/PPI methods that vote pre-t0 labels *retrieve the pre-known
answer*. The `net` (STRING Net-KNN) component is hit hardest: no-knowledge MF
**0.803 → 0.441**, mean **0.475 → 0.347**, once the contamination is removed.
(Consistency check: GOAlpha report their Net-KNN *drops* on the genuinely-novel
private set, 0.29→0.26 — the opposite of our contaminated rise.) The *integrator*
is largely robust (0.532 → 0.521) because it blends `net` with non-leaking signal,
but every standalone no-knowledge-MF and `net` number on the dirty GT is inflated.

**Fix.** Don't trust GAF dates; use the canonical pre-t0 label set instead. A
(protein, aspect) is genuinely no-knowledge iff the protein has **no `train_terms`
entry in that aspect**. `pipeline/build_clean_gt.py` builds the leak-free GT:

```bash
python pipeline/build_clean_gt.py \
  --gt-no gt/gt_no.tsv --train-terms <train_terms.tsv> --obo <go.obo> \
  --out gt/gt_no_clean.tsv --novel-proteins-out gt/gt_no_novelprot.tsv
```
(Aspect comes from the OBO `namespace`, the authoritative source cafaeval uses —
**not** from a go-dag closure, which is order-dependent for terms with cross-aspect
`part_of` ancestors.)

Corrected, leak-free numbers are in [`RESULTS.md` → Corrected results](RESULTS.md#corrected-results-2026-06-23).
**A proper re-freeze should retrain the integrator on the clean GT** (current clean
numbers score a model still *trained* on dirty labels — robust, but not the final
word). A deeper fix would re-derive `testsuperset_exp_annots.tsv` from a *dated GOA
diff* (annotations present in the post-t0 snapshot but absent from a ≤t0 snapshot)
rather than the per-annotation last-modified date.

### 1.1 Components (the evidence, produced upstream)

Each component is a TSV `protein⇥term⇥score`. The integrator first **propagates
every component to GO ancestors by max** so the columns are comparable, then
forms, per aspect, the union of predicted terms as (protein, term) candidates.

| component | signal | where it comes from |
|---|---|---|
| `diam`     | BLAST-KNN homology      | `DiamondPredictor` (DIAMOND vs pre-t0 SwissProt) |
| `foldseek` | structure-KNN          | `FoldSeekPredictor` |
| `clean`    | EC → GO                | `CleanPredictor` |
| `interpro` | domain → GO            | `InterProScanPredictor` |
| `mlp`      | sequence-PLM           | ESM2-650M MLP head (10× ensemble; ORIX/IBEX) |
| `prostt5`  | structure-aware PLM    | ProstT5 MLP head (ORIX) — the complementary lever |
| `net`      | PPI guilt-by-assoc.    | **`make net`** — Net-KNN over STRING v12 |
| `lit`      | literature text-kNN    | **`make lit`** — BM25 over SwissProt names |

`diam/foldseek/clean/interpro/mlp/prostt5` are produced by GSPA's predictors /
GPU PLM heads (§4); `net`/`lit` are CPU-built by this pipeline.

### 1.2 The integrator (per-aspect logistic regression)

For each aspect (MF/BP/CC) and each (protein, term) candidate:

- **features = the component scores only** (one column per component, propagated
  value or 0 if absent). *Nothing else* — no term IA, no term frequency.
- **label = 1** if the term is in the protein's propagated experimental GT for
  that aspect, else 0.
- Fit `StandardScaler` → `LogisticRegression(class_weight='balanced')`.
- The frozen artifact stores, per aspect, the scaler `mean`/`scale` and the
  logistic `coef`/`intercept`. Apply = `sigmoid(Σ coefᵢ·(xᵢ−meanᵢ)/scaleᵢ + b)`
  with the numerically-stable, branch-by-sign sigmoid.

**Why scores-only / logreg (the leak-free choice).** A richer XGBoost variant
with IA + log-frequency features leaked badly — it memorised term base-rates and
hit MF f_w 0.963 on in-distribution CV, which does **not** generalise to novel
proteins. A linear model over component scores *cannot* memorise term identity,
so its cross-validated score is honest. That variant was discarded; `--model
logreg --features scores` is the only configuration `--save-model` accepts (the
guard is enforced and tested).

### 1.3 Honest evaluation (out-of-fold)

The reported generalisation number is **GroupKFold-by-protein out-of-fold**: each
protein's predictions come from a model trained on other proteins, so no protein
informs its own score. `train_integrator.py` writes the OOF predictions to
`<out>/preds/ltr.tsv`; `cafaeval` scores them. The *shipped* JSON is the
full-data fit (all proteins) — applying it in-sample scores 0.489, consistent
with the 0.483 OOF, confirming the deployable artifact is faithful.

### 1.4 Blind, fully CAFA-faithful variant (temporal integrity)

The default integrator's ~8 weights are tuned (GroupKFold) on no-knowledge
*test* labels — which you would not have at submission time. To quantify that
optimism, `pipeline/train_head_oof.py` (GPU, ORIX) produces k-fold OOF PLM-head
scores for the *pre-t0 train* population, and an integrator trained **blind** on
a 25k-protein pre-t0 validation set is compared to the test-trained one:

| integrator trained on … | official 3-class f_w |
|---|---|
| blind pre-t0 validation (never sees test) | **0.561** |
| no-knowledge test labels | 0.566 |

**Test-label tuning is worth only +0.005**, and `net`'s weight *increases* under
blind training — so the result is essentially CAFA-submission-faithful. The only
non-faithful element is `lit` (query-name leakage risk), which is why it ships
optional. Full detail in `RESULTS.md` §"Temporal integrity".

---

## 2. Performance on CAFA6

Faithful reconstruction, official `cafaeval` + `IA.tsv`, t0 = 2026-02-02.

> **⚠️ Corrected 2026-06-23.** Two earlier reporting errors are fixed here; see
> the [Correction banner in `RESULTS.md`](RESULTS.md#-correction-2026-06-23--read-this-before-trusting-any-number-below)
> for the full audit. In short: (1) the CAFA6 headline is the **no-knowledge**
> f_w, *not* the 3-class mean — comparing our 3-class mean (0.647) to GOAlpha's
> no-knowledge 0.524 was apples-to-oranges; (2) the no-knowledge MF GT was
> contaminated by the **GAF last-modified-date bug** (see §1.0a), which inflated
> the `net` component (MF 0.803 → 0.441 on a leak-free GT). The integrator is
> robust to it (0.532 → 0.521 clean), so the honest headline is **level with
> GOAlpha, not ahead**.

**Headline (corrected).** Our real CAFA6 entry (team *HoehndorfLab*) scored
**0.37749, rank 263/2177**; the winner (GOAlpha) **0.524**. With *zero* new
models, learned integration of the same components reaches **no-knowledge
IA-weighted f_w ≈ 0.52** (0.532 on the contaminated GT, **0.521 on the leak-free
GT**) — **level with GOAlpha, not ahead.** (Do not quote the 3-class mean,
~0.63–0.65, against GOAlpha: it is inflated by easy partial-knowledge proteins.)

**Per-component, no-knowledge f_w** (corrected; `clean` = leak-free per-aspect GT):

| component | dirty f_w | **clean f_w** | note |
|---|---|---|---|
| `net`      | 0.475 | **0.347** | MF **0.803 → 0.441**: ~half was retrieving pre-known labels (the GAF-date leak) |
| `mlp`      | 0.449 | **0.472** | strongest PLM head; *rises* on clean (doesn't leak) |
| `prostt5`  | 0.435 | — | structure-aware PLM, complementary to `mlp` |
| `diam`     | 0.389 | — | BLAST-KNN |
| `foldseek` | 0.364 | — | structure-KNN |
| `interpro` | 0.165 | — | weak; redundant with PLM/homology |
| `clean`    | ~0    | — | EC→GO, niche |
| our 2025 submission (max-merge) | 0.359 | 0.454 | max-merge destroyed MF; *rises* on clean |
| **integrator (6+net)** | **0.532** | **0.521** | robust to the leak — blends `net` with non-leaking signal |

**Ablation by knowledge class** (the 6-comp model =
`diam,foldseek,clean,interpro,mlp,prostt5`). **These per-class numbers are on the
contaminated GT** (esp. no-knowledge MF) and the **3-class mean must not be
compared to GOAlpha** — kept for the relative component story only:

| model | no-knowledge | limited | partial | 3-class mean *(not LB-comparable)* |
|---|---|---|---|---|
| 6-comp (baseline)         | 0.489 | 0.630 | 0.768 | 0.629 |
| **6-comp + net** *(ship)* | 0.538 | 0.654 | 0.748 | 0.647 (+0.018) |
| 6-comp + lit + net        | 0.553 | 0.646 | 0.715 | 0.638 |
| 6-comp + lit              | 0.507 | 0.631 | 0.720 | 0.620 |

`net`'s apparent no-knowledge win is **substantially the GAF-date leak** (see the
corrected per-component table above and §1.0a); on a leak-free GT its standalone
mean is 0.345, and its value in the *integrator* is much smaller than the dirty
+0.018/+0.049 figures suggest. `lit` additionally carries a name-leak caveat. See
`RESULTS.md` for the corrected results and the full audit.

**Full leave-one-out / cumulative / aggregator ablation** (reproduce with
`pipeline/ablation.py`; raw numbers in `ablation_no_results.tsv`, deep dive in
`RESULTS.md`). The headlines on no-knowledge (novel proteins):

* **`net` is the only large marginal win** (+0.049 leave-one-out). Every other
  component is within ±0.013; `clean`, `foldseek`, `mlp` sit at or below zero once
  `net` is present — six PLM/homology/domain channels are largely redundant for
  novel proteins. **`clean` actively *decreases* the cumulative score** (−0.001)
  and is kept only because it is near-free and marginally helps partial-knowledge.
* **Gradient-boosted trees do *not* beat the linear stacker.** Leak-free XGBoost
  over the same 6 component scores scores **0.478 vs logreg 0.483**; its extra
  capacity overfits the tiny feature set. Giving XGBoost term-identity features
  (IA + log-freq) **leaks** — MF jumps to a non-credible **0.963** (term base-rate
  memorisation; GroupKFold splits by protein but terms are shared across folds).
  This is why DG++ ships the **linear, scores-only** integrator. The bottleneck is
  signal diversity + generalisation, not aggregator sophistication.

### DeepGO-PlusPlus-Light (no GPU) — corrected 2026-06-23

DG++-Light drops the GPU PLM heads (`mlp`, `prostt5`, `esm2_3b`) and uses CPU
components (DIAMOND, FoldSeek, InterProScan/-LR, STRING Net-KNN, optional BM25
`lit`, and a CPU 1D-CNN `cnn` —`pipeline/build_cnn_component.py`).

> **⚠️ The old DG++-Light claims were leak-driven and are RETRACTED.** They were
> computed on the GAF-date-contaminated no-knowledge GT (§1.0a), and the Light
> panels are exactly the `net`/`net_union`-centric models the leak inflates most.
> On the leak-free GT (full table in
> [`RESULTS.md`](RESULTS.md#corrected-deepgo-plusplus-light-no-gpu--the-beats-gpu-claim-reverses)):

**Panels** (no-knowledge f_w):

| model | dirty mean | **clean mean** | note |
|---|---|---|---|
| full GPU integrator (6+net) | 0.532 | **0.521** | robust to the leak — **wins on novel proteins** |
| `light.json` (diam,foldseek,interpro,net) | 0.550 | **0.488** | best Light panel on clean |
| `light_cnn.json` (cnn,diam,foldseek,interpro,net) | 0.516 | **0.470** | — |
| `light_cpu.json` (diam,interpro,net_union) *(was "recommended")* | 0.564 | **0.464** | most leak-inflated |
| `light_fast.json` (diam,net_union) *(webservice default)* | 0.562 | **0.451** | weakest; net_union-only |

**Full standalone-component ablation** (every component incl. the CNN, scored
directly; dirty → clean; full table + GPU refs in [`RESULTS.md`](RESULTS.md#corrected-deepgo-plusplus-light-no-gpu--the-beats-gpu-claim-reverses)):

| component | dirty | **clean** | | component | dirty | **clean** |
|---|---|---|---|---|---|---|
| `lit` (CPU)        | 0.445 | **0.459** | | `cnn` (CPU)        | 0.209 | 0.268 |
| `diam` (CPU)       | 0.388 | **0.451** | | `interpro` (CPU)   | 0.164 | 0.215 |
| `foldseek` (GPU-gated†)| 0.363 | **0.428** | | `net_bridge` (CPU bridge) | 0.504 | **0.358** ↓ |
| `interpro_lr` (CPU)| 0.233 | 0.325 | | `net_union` (CPU bridge) | 0.506 | **0.353** ↓ |
| `clean` (CPU)      | 0.007 | 0.116 | | `net` (CPU, STRING‡)| 0.475 | **0.347** ↓ |

**Two query-input gates that fail exactly on novel proteins** (so a component being
"CPU" on a benchmark ≠ usable on a fresh sequence):
- **† `foldseek` is not unconditionally CPU.** FoldSeek search is CPU but needs a
  **structure or ProstT5-derived 3Di** per query — a CPU lookup only for AFDB-covered
  proteins; folding/encoding a novel sequence (ESMFold/ProstT5 3Di) needs a **GPU**.
- **‡ plain `net` needs the query to be a STRING node.** Novel proteins are *not* in
  STRING, so plain `net` = 0 for them; the deployable net signal is the **DIAMOND
  homology bridge** `net_union`/`net_bridge` (query → STRING-member homolog → its
  neighbours' labels — CPU). The CAFA6 no-knowledge set is ~87 % STRING members, which
  is why plain `net` still scores 0.347 here.

* **The `net` family is the only group that *drops* on clean** (the leak); every
  other component *rises*.
* **On genuinely-novel proteins the strongest signals are `lit` (0.459), `diam`
  (0.451), `foldseek` (0.428†) — NOT `net` (0.347).**
* **The full GPU model BEATS DG++-Light** (0.521 vs 0.45–0.49); the PLM heads are
  not redundant. DG++-Light trades ~0.03–0.07 f_w for no-GPU deployment.
* The CPU **1D-CNN** is the weakest sequence signal (clean 0.268) but doesn't leak;
  its value stays *coverage* of true orphans, not accuracy.

**Leak-free net-bridge accuracy** (clean-B = proteins ∉ `train_terms` → bridge cannot
self-match; 0/251 in train_terms verified). Split by STRING membership:

| component | all clean-B (251) | STRING-member (214) | non-STRING (37) |
|---|---|---|---|
| `net` (direct) | 0.339 | 0.380 | **0.000** |
| `net_union` | 0.351 | 0.379 | **0.252** |
| `net_bridge` | 0.353 | 0.377 | **0.252** |
| `diam` (ref) | 0.442 | 0.414 | **0.598** |

The bridge's *actual* leak-free accuracy is ~0.35 (0.252 on non-STRING, where direct
`net`=0 — so it genuinely fires). **But `diam` dominates it everywhere, including on
non-STRING proteins (0.598 vs 0.252)** — the one case the bridge exists for. The bridge
just transfers a DIAMOND homolog's *neighbours'* labels; transferring the homolog's
*own* labels (`diam`) is far better. So **the net/PPI family is redundant with `diam`
and a candidate to drop entirely** (n=37 non-STRING small but the gap is large).

**At the next clean re-freeze, redesign DG++-Light around `diam` + `lit`** (strictly-CPU,
any novel sequence) + `cnn` for orphan coverage, `foldseek` only in a
*no-GPU-given-structures/3Di* tier, and **drop the net/PPI family** (dominated by
`diam`). Re-pick the default (`light.json` currently leads, not `light_cpu`/`light_fast`).
Re-pick/re-freeze with `pipeline/build_clean_gt.py` + `pipeline/eval_light.py` +
`train_integrator.py --save-model` (clean set small, n≈282, ±0.02). No predictor code
changes — the component list is in the JSON.

**DONE (2026-06-23) — net-free, leak-free retrain.** Re-trained the Light integrator
without the net family on the clean GT: best net-free panel **0.508**
(`deepgo_plusplus_light_clean.json` = diam,foldseek,interpro,lit), strictly-CPU
**0.500** (`deepgo_plusplus_light_cpu_clean.json` = diam,interpro,lit,cnn) — both
**beat the old net-based Light models (0.464–0.488)** and nearly match the full GPU
model (0.521). Frozen on n≈282 (interim; retrain on a large pre-t0 set for production).
For the **roadmap to a genuinely better CPU predictor** (embedding-kNN third channel,
small-PLM/light-attention head, ProteInfer dilated CNN, PLM distillation, loss
upgrades — with citations), see `RESULTS.md` → "Better CPU predictor — literature-driven
roadmap". GSPA already ships `Esm2CentroidPredictor` (embedding-kNN) and
`ProteInferPredictor` to build the top two on.

---

## 3. How to update / retrain at a new release

Retraining is reproducible and re-runnable for any UniProt / STRING / CAFA
release. The CPU path (index → ground truth → net/lit → train → eval) is driven
by the `Makefile`; the GPU components (§4) are produced upstream.

```bash
cd deepgo-plusplus
cp config.mk.example config.mk          # point at the new release's inputs
$EDITOR VERSIONS.md                      # record the exact release versions
uv sync                                  # install the pinned env (commit uv.lock)

uv run make check-inputs                 # every declared input path resolves?
uv run make index                        # UniProt pass -> text + STRING + taxon index
uv run make gt                           # dated GOA -> gt_{no,limited,partial,all}
uv run make net                          # Net-KNN component (needs STRING_DIR; §4)
# uv run make lit                        # optional literature component
# place the GPU/homology components under build/components/ (§4)
uv run make integrator                   # train + freeze -> models/deepgo_plusplus_integrator.json
uv run make eval                         # official cafaeval on the OOF predictions
uv run make test                         # regression suite stays green
```

To include `net` in the frozen model, set in `config.mk`:
`COMPONENT_LIST = diam,foldseek,clean,interpro,mlp,prostt5,net` (and append
`,lit` with `make lit` for the literature variant).

### Per-release checklist

1. **Decide t0.** For a *CAFA-faithful* re-score keep `T0 = 20260202` and a GO
   ontology **≤ t0**. For a *production* model on current data, set t0 to the new
   evaluation cutoff and use the latest GO release.
2. **Update inputs** in `config.mk` (paths) and **record releases** in
   `VERSIONS.md` (SwissProt release, STRING version, GO release, CAFA artifacts).
3. **Rebuild ground truth** if GOA changed (`make gt`).
4. **Rebuild components**: `make net`/`make lit` for the CPU ones; rerun the
   upstream GPU/homology components (§4) if UniProt changed (DIAMOND/FoldSeek
   reference DBs and PLM embeddings are all release-sensitive).
5. **Retrain + freeze** (`make integrator`) and **evaluate** (`make eval`).
6. **Promote** the new JSON: it already lands in `models/`; point GSPA at it via
   `--deepgo-plusplus-integrator models/<new>.json` or the YAML config.
7. **Re-run tests** (`make test`) and, if any dependency changed, `uv lock` and
   commit the updated `uv.lock`.

### What is and isn't reproducible here

- **Fully reproducible (CPU, this Makefile):** the text/STRING/taxon index, the
  ground truth, the `net` and `lit` components, integrator training, evaluation.
- **Upstream (heavy / GPU, documented in §4):** the `diam/foldseek/clean/interpro`
  homology+domain components and the `mlp/prostt5` PLM heads. Drop their score
  TSVs into `build/components/` before `make integrator`.

---

## 4. Data sources & (re)download

Declare these in `config.mk`. Record the exact release each retrain uses in
`VERSIONS.md`. **Compute placement:** the multi-GB GOA filter belongs on **IBEX**;
GPU embeddings on **ORIX** (pixi); STRING + the CPU pipeline run locally or on
`ws` (the 56-core office box). Avoid the metered laptop link for the big pulls.

### 4.1 UniProt SwissProt flat file → `net`/`lit` text index

Per-accession identification text, `DR STRING` xref, and `OX` taxon.

```bash
# ~250 MB+ ; verify you want current_release vs a pinned release directory
curl -L -o uniprot_sprot.dat.gz \
  https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
# set UNIPROT_SPROT_DAT in config.mk ; then: uv run make index
```

### 4.2 UniProt-GOA → CAFA6 ground truth (`testsuperset_exp_annots.tsv`)

The full GAF is ~11 GB — **filter on IBEX**, don't pull it to a workstation. Keep
test-superset proteins with **experimental** evidence and project to the columns
`build_groundtruth.py` expects (`acc⇥GO⇥F|P|C⇥evidence⇥date`):

```bash
# on IBEX (dm/login): goa_uniprot_all.gaf.gz lives under the EBI GO mirror
#   https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz
# GAF cols: 2=accession 5=GO 7=evidence 9=aspect(F/P/C) 14=date(YYYYMMDD)
EXP='EXP|IDA|IPI|IMP|IGI|IEP|TAS|IC|HTP|HDA|HMP|HGI|HEP'
zcat goa_uniprot_all.gaf.gz \
 | awk -F'\t' -v OFS='\t' -v ids=testsuperset_ids.txt -v e="$EXP" '
     BEGIN{while((getline l<ids)>0) keep[l]=1}
     !/^!/ && ($2 in keep) && ($7 ~ ("^("e")$")) {print $2,$5,$9,$7,$14}' \
 > testsuperset_exp_annots.tsv
# set TESTSUPERSET_ANNOTS in config.mk ; then: uv run make gt
```

`testsuperset_ids.txt` = the CAFA6 test-superset accessions (§4.4). The current
snapshot is staged at `../benchmark/neural/cafa6_recon/testsuperset_exp_annots.tsv`
and on IBEX `/ibex/scratch/projects/c2014/cafa6/recon/`.

### 4.3 STRING per-species links → `net`

The test set spans only ~68 species (~97 % of proteins have a STRING id), so
download **per species**, not the full dump:

```
https://stringdb-downloads.org/download/protein.links.v12.0/{taxid}.protein.links.v12.0.txt.gz
```

`../benchmark/neural/run_net_ws.sh` downloads all needed species with gzip
integrity-check + retry (drops corrupt files) and builds `net` in one step; point
`STRING_DIR` at the result. **Taxid gotcha:** STRING uses *canonical* taxids, not
the testsuperset NCBI taxids — e.g. yeast `559292 → 4932`, E.coli K-12
`83333 → 511145` (others 404). SwissProt `DR STRING` ids already carry the
`4932.`/`511145.` prefix, so save the correct-taxid file under the OX-taxon
filename. A complete, validated mirror (68/68) is at
`../benchmark/neural/cafa6_recon/string/`.

### 4.4 CAFA6 official artifacts (Kaggle)

`IA.tsv`, `train_terms.tsv` (the pre-t0 training labels — *the data we actually
had at submission*), `testsuperset-taxon-list.tsv`, the test-superset FASTA, and
`sample_submission`. Download from the **CAFA6 Kaggle competition** with the
**leechuck** account — only that `kaggle.json` is a CAFA6 participant
(`roberthoehndorf` → HTTP 403):

```bash
kaggle competitions download -c <cafa6-competition-slug>   # leechuck kaggle.json
# unzip into ~/Public/software/cafa6/kaggle-official/
```

Already mirrored at `~/Public/software/cafa6/kaggle-official/` and on IBEX
`/ibex/scratch/projects/c2014/cafa6/data/`.

### 4.5 GO ontology + DAG

- **`go-basic.obo` / `go.obo`** — for `cafaeval`. Pin a release **≤ t0** for a
  faithful re-score: `http://purl.obolibrary.org/obo/go/go-basic.obo`
  (release archives under `.../go/releases/`).
- **`go-dag.tsv`** (`child⇥ancestor`) — the integrator/apply propagation table.
  It is the **transitive closure** of the ontology over `is_a` + `part_of`
  (each term → *all* its ancestors; the pipeline does not close it at runtime).
  Use `~/Public/software/cafa6/go-dag.tsv`, or regenerate it from the matching
  `go-basic.obo` release. Set `OBO` and `DAG` in `config.mk`.

### 4.6 Upstream component scores (GPU / homology)

The six non-net/lit components are produced outside this Makefile and dropped
into `build/components/<component>.tsv[.gz]`:

- **`diam`, `foldseek`, `clean`, `interpro`** — run GSPA's predictors
  (`DiamondPredictor`, `FoldSeekPredictor`, `CleanPredictor`,
  `InterProScanPredictor`) on the testsuperset FASTA with **pre-t0 reference
  DBs**, emitting `protein⇥term⇥score`.
- **`mlp`, `prostt5`** — PLM heads on ORIX (pixi env `gspa-glm`):
  `benchmark/neural/extract_embeddings.py` /
  `extract_{esm2_3b,prostt5}.sbatch` → `train_head.py` (per-aspect MLP,
  frequency-weighted BCE), then predict the testsuperset. Working dir
  `/mnt/data/u/hohndor/cafa6-plm/`. The `mlp` component is the original CAFA6
  10× ESM2-650M ensemble (IBEX `…/cafa6/data/mlp{0..9}_*.th`).

A staged snapshot of all components is at
`../benchmark/neural/cafa6_recon/components/`.
