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
  `IA.tsv`: IA-weighted max-F (`f_w`), `-norm cafa -prop max`. The headline
  metric is the mean over MF/BP/CC, then over the three knowledge classes.

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

**Headline.** Our real CAFA6 entry (team *HoehndorfLab*) scored **0.37749, rank
263/2177**; the winner (GOAlpha) **0.524**. With *zero* new models, learned
integration of the same components recovers **novel-protein (no-knowledge)
IA-weighted f_w 0.359 → 0.483** — a would-be top-10 result — and adding the
network signal lifts the official 3-class mean to **0.647**.

**Per-component, no-knowledge f_w** (why integration, not a single model, wins):

| component | f_w | note |
|---|---|---|
| `net`      | 0.475 | **best single-component MF 0.803** (guilt-by-association) |
| `mlp`      | 0.447 | strongest PLM head |
| `lit`      | 0.445 | MF 0.539 (name-leak caveat) |
| `prostt5`  | 0.435 | **MF 0.466** — best PLM MF, complementary to `mlp` |
| `diam`     | 0.389 | BLAST-KNN |
| `foldseek` | 0.364 | structure-KNN |
| `interpro` | 0.165 | weak; redundant with PLM/homology |
| `clean`    | ~0    | EC→GO, niche |
| our 2025 submission (max-merge) | 0.359 | **below `mlp` alone** — max-merge destroyed MF |

**Ablation** (IA-weighted f_w by knowledge class; the 6-comp model =
`diam,foldseek,clean,interpro,mlp,prostt5`):

| model | no-knowledge | limited | partial | **official 3-class** |
|---|---|---|---|---|
| 6-comp (baseline)         | 0.489 | 0.630 | 0.768 | 0.629 |
| **6-comp + net** *(ship)* | 0.538 | 0.654 | 0.748 | **0.647** (+0.018) |
| 6-comp + lit + net        | 0.553 | 0.646 | 0.715 | 0.638 |
| 6-comp + lit              | 0.507 | 0.631 | 0.720 | 0.620 |

`net` is the clean, leak-free win (ship `deepgo_plusplus_integrator_net.json`).
`lit` raises no-knowledge but drags the 3-class mean (and carries a name-leak
caveat), so it ships as an optional no-knowledge booster
(`deepgo_plusplus_integrator_lit_net.json`). See `RESULTS.md` for the full story,
per-aspect MF/BP/CC splits, and the remaining gap to GOAlpha.

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

### DeepGO-PlusPlus-Light (no GPU)

DG++-Light drops the GPU PLM heads (`mlp`, `prostt5`, `esm2_3b`) and uses CPU
components (DIAMOND, FoldSeek, InterProScan/-LR, STRING Net-KNN, optional BM25
`lit`, and a CPU 1D-CNN `cnn` that replaces the PLM heads —
`pipeline/build_cnn_component.py`). Two CAFA6 findings (no-knowledge f_w; tables
in `RESULTS.md`):

* **No-GPU panels beat the full GPU model (0.532) on novel proteins** — the PLM
  heads are redundant with `net`. But `foldseek` needs a query *structure* (CPU
  lookup if in AlphaFold DB, else GPU folding), and `net` only fires for STRING
  members. The **strictly no-GPU, any-sequence** model is
  `models/deepgo_plusplus_light_cpu.json` = `diam,interpro,net_union` (**0.564**,
  **recommended**); `models/deepgo_plusplus_light.json` = `diam,foldseek,interpro,net`
  (0.550) is for when AFDB structures are available.
* **`net_union` extends `net` to proteins not in STRING** via a DIAMOND homology
  bridge (`build_net_bridge.py`): query → pre-t0 STRING-member homolog → its
  neighbours' pre-t0 labels. Leak-safe (pre-t0 homolog DB; novel queries have no
  pre-t0 labels so can't self-match). Takes the 356 no-STRING no-knowledge proteins
  from f_w 0 → 0.42 and lifts the structure-free panel 0.544 → 0.564.
* **The 1D-CNN does not improve novel-protein f_w** (standalone 0.206; −0.012 in
  panel) — a sequence model trained on pre-t0 data hits the same generalisation
  wall. Its value is *coverage* (it predicts for every protein, incl. orphans), so
  it ships as the separate `models/deepgo_plusplus_light_cnn.json` (0.516).

Re-pick / re-freeze at each release with `pipeline/eval_light.py` +
`train_integrator.py --save-model`; the `cnn` component is rebuilt with
`build_cnn_component.py` (`--save-model`/`--load-model` to apply without
retraining). No predictor code changes — the component list is read from the JSON.

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
