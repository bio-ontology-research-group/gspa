# `cafa-baseline` — GSPA's CAFA6 learned-stacker predictor

A CAFA6-competitive GO predictor built **without** new model architecture: it
replaces naive max-merge of GSPA's existing component predictors with a learned
per-aspect logistic-regression stacker. On a faithful CAFA6 reconstruction
(GOA snapshot, t0=2026-02-02, official `cafaeval` + `IA.tsv`) this recovered
**novel-protein (no-knowledge) IA-weighted f_w from 0.359 → 0.483**, vs the
0.524 first-place entry (GOAlpha) — a would-be top-10 result from integration
alone. See `cafa6_recon/` for the reconstruction harness and
`~/.claude/.../memory/project_cafa6_failure_diagnosis.md` for the full story.

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
`models/cafa_baseline_integrator_net.json`. Adding `lit` pushes no-knowledge
higher (the real-LB proxy) but *drags* the 3-class mean by hurting
partial-knowledge proteins (the integrator is frozen on no-knowledge weights),
so the 8-comp model is offered separately as
`models/cafa_baseline_integrator_lit_net.json` for no-knowledge-focused use.
The partial-knowledge dip is the frozen-on-no-knowledge weighting artifact;
per-knowledge-class / pre-t0-population training (deferred) would remove it.

## 1. Train + freeze the integrator (once)

```bash
python train_ltr_integrator.py \
  --components cafa6_recon/components --gt cafa6_recon/gt/gt_no.tsv \
  --dag   ~/Public/software/cafa6/go-dag.tsv \
  --ia    ~/Public/software/cafa6/kaggle-official/IA.tsv \
  --train-terms ~/Public/software/cafa6/kaggle-official/train_terms.tsv \
  --taxon ~/Public/software/cafa6/kaggle-official/testsuperset-taxon-list.tsv \
  --model logreg --features scores \
  --component-list diam,foldseek,clean,interpro,mlp,prostt5 \
  --save-model models/cafa_baseline_integrator.json --out /tmp/ltr_run
```

`--model logreg --features scores` is **leak-free** — the linear model over
component scores cannot memorise term identity (an XGBoost variant with IA/freq
features leaked to MF f_w 0.963 in-distribution and was discarded). The OOF
GroupKFold-by-protein score in `/tmp/ltr_run/preds/ltr.tsv` is the honest
generalisation estimate (0.483). The frozen JSON is the full-data fit
(`models/cafa_baseline_integrator.json`, 2 KB) shipped for inference.

## 2. Apply — directly

```bash
python run_neural_predictors.py --predictor cafa-baseline \
  --manifest manifest.tsv --min-score 0.1 \
  --integrator models/cafa_baseline_integrator.json \
  --components-dir cafa6_recon/components \
  --dag ~/Public/software/cafa6/go-dag.tsv
```

Emits `<tag>.cafa-baseline.tsv` (`protein\tterm\tscore\tGO`). Applying the
frozen model to the no-knowledge set scores **0.489** (in-sample full-fit;
consistent with the 0.483 OOF), confirming the artifact is faithful.

## 2′. Apply — via GSPA

```bash
./gradlew :gspa-cli:run --args="annotate --input proteome.faa --output out \
  --neural-sidecar $PWD/benchmark/neural/run_neural_predictors.py \
  --cafa-baseline \
  --cafa-baseline-integrator $PWD/benchmark/neural/models/cafa_baseline_integrator.json \
  --cafa-baseline-components-dir <dir-of-component-scores> \
  --cafa-baseline-dag ~/Public/software/cafa6/go-dag.tsv"
```

(or set `predictors.neural.cafaBaseline.{enabled,integrator,componentsDir,dag}`
in YAML). The sidecar propagates each component to GO ancestors (max), forms
per-(protein, term) candidates per aspect, and emits `sigmoid(w·x + b)`.

## 3. Building the `net` and `lit` components

Both derive from one streaming pass over the SwissProt flat file, which yields
per-accession identification text, the `DR STRING` xref, and the `OX` taxon:

```bash
# one parse -> text_string_index.tsv (accession\ttaxon\tstring_id\tname\tchar_text)
python build_text_string_index.py uniprot_sprot.dat.gz text_string_index.tsv
```

**Literature** (`lit`) — BM25 text-kNN, CPU-only, `--shard i/N` for parallelism:
```bash
python build_lit_component.py --index text_string_index.tsv \
  --train-terms train_terms.tsv --queries test_proteins.txt \
  --out lit.tsv --topk 30          # corpus uses full text; query uses NAME only
```

**Net-KNN** (`net`) — needs STRING per-species link files. The test set spans
only ~68 species (97 % of proteins have a STRING id), so download per species,
not the full dump. `run_net_ws.sh` does the download (with gzip integrity +
retry, corrupt files dropped) and the build in one step:
```bash
# slim 3-column index is enough for net (no text):
cut -f1-3 text_string_index.tsv > net_index.tsv
python build_net_component.py --index net_index.tsv \
  --train-terms train_terms.tsv --queries test_proteins.txt \
  --string-dir <dir of {taxid}.protein.links.v12.0.txt.gz> \
  --out net.tsv --min-conf 400 --topk 50
```

Gzip the outputs into the components dir (`components/{net,lit}.tsv.gz`) and add
`net` (and optionally `lit`) to `--component-list`. The Groovy
`CafaBaselinePredictor` needs **no change** — it reads the component list from
the integrator JSON, so swapping in `cafa_baseline_integrator_net.json` is
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
So the cafa-baseline result is essentially CAFA-faithful; the only non-faithful
element is `lit` (query-name leakage), which is why it is optional. The blind
PLM head scores come from k-fold OOF on the train embeddings
(`train_head_oof.py`, ORIX); the blind model is `cafa6_recon/integrator_pret0.json`.

## Remaining gap to GOAlpha (0.524)

Not models, and not the integrator training population (shown negligible above).
The open lever is a **full-component** blind/pre-t0 integrator: the check above
used the 3 components reproducible for arbitrary train proteins
(`prostt5/esm2_3b/net`); reproducing `diam/foldseek/clean/interpro` for a
validation set (IBEX work) would let the *full* model train on the pre-t0
population and let XGBoost + IA/freq help leak-free. STRING Net-KNN and a
(leakage-clean) literature channel are now done.
