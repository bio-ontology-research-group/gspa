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

Per-component no-knowledge f_w: mlp 0.447, prostt5 0.435 (best single-component
**MF** 0.466), diam 0.389, foldseek 0.364, interpro 0.165, clean ~0. ProstT5 is
the only PLM that *complements* the MLP rather than cannibalising it; adding
ESM2-3B on top regresses (correlated-PLM collinearity).

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

## Remaining gap to GOAlpha (0.524)

Not models — component-hunting has hit diminishing returns. The two open levers:
(1) train the integrator on the **pre-t0 population** (re-run components on the
pre-t0 train FASTA) so XGBoost + IA/freq can help leak-free; (2) add the
missing **Net-KNN (STRING PPI)** and a down-weighted literature channel.
