# DeepGO-PlusPlus

GSPA's CAFA6-competitive GO predictor — a **learned per-aspect stacker** over
the heterogeneous evidence GSPA already produces, plus a network (STRING PPI)
and literature channel. It is the productionised result of the CAFA6 post-mortem
(see `RESULTS.md`): replacing naive max-merge with a frozen logistic integrator
recovered novel-protein IA-weighted f_w from **0.359 → 0.483** (vs the 0.524
first-place entry), and adding Net-KNN lifts the official 3-class mean to
**0.647**.

This folder is **self-contained and re-runnable at every UniProt / STRING
release**. The frozen model under `models/` is an artifact; the recipe that
produces it lives here so it can always be rebuilt and audited.

**Full documentation:** [`TRAINING.md`](TRAINING.md) — how it was trained, how to
update/retrain, where every input comes from and how to (re)download it, and the
CAFA6 performance numbers. [`RESULTS.md`](RESULTS.md) — full ablation tables +
temporal-integrity audit. [`VERSIONS.md`](VERSIONS.md) — input-release provenance.

```
deepgo-plusplus/
├── README.md            ← you are here (overview + GSPA wiring)
├── TRAINING.md          ← how it was trained · update/retrain · data (re)download · CAFA6 perf
├── RESULTS.md           ← full ablation tables + temporal-integrity audit
├── VERSIONS.md          ← which external releases each frozen model used
├── Makefile             ← the reproducible retrain DAG  (`make help`)
├── pyproject.toml       ← pinned deps for `uv`  (commit uv.lock for exact repro)
├── config.mk.example    ← input paths — copy to config.mk
├── models/              ← frozen integrators (shipped, retrainable)
│   ├── deepgo_plusplus_integrator.json          (6-comp, default)
│   ├── deepgo_plusplus_integrator_net.json      (+net, recommended)
│   ├── deepgo_plusplus_integrator_lit_net.json  (+lit+net)
│   ├── deepgo_plusplus_light_fast.json          (diam+net_union — webservice default)
│   ├── deepgo_plusplus_light_fast_cnn.json      (+cnn — orphan coverage)
│   ├── deepgo_plusplus_light_cpu.json           (+interpro — strictly no-GPU best)
│   ├── deepgo_plusplus_light_full.json          (+interpro+cnn)
│   ├── deepgo_plusplus_light.json               (no-GPU given AFDB structures)
│   └── deepgo_plusplus_light_cnn.json           (offline 6-comp + cnn, coverage)
├── ablation_no_results.tsv  ← committed ablation numbers (pipeline/ablation.py)
├── pipeline/            ← the retraining + apply scripts
│   ├── build_text_string_index.py   one UniProt pass → text + STRING + taxon
│   ├── build_groundtruth.py         dated GOA → CAFA6 knowledge-class GT
│   ├── build_net_component.py       Net-KNN over STRING PPI
│   ├── build_lit_component.py       BM25 literature text-kNN
│   ├── build_cnn_component.py       CPU 1D-CNN over sequence (DG++-Light)
│   ├── build_net_bridge.py          DIAMOND-bridged Net-KNN for non-STRING proteins
│   ├── extract_sprot_fasta.py       SwissProt .dat → FASTA for given accessions
│   ├── train_integrator.py          train + freeze the per-aspect logreg
│   ├── train_head_oof.py            (GPU) k-fold OOF PLM heads for blind model
│   ├── apply_integrator.py          apply a frozen model to a FASTA
│   ├── apply_net_bridge.py          fast bridge at inference (index lookup, no STRING scan)
│   ├── ablation.py                  LOO / cumulative / GBT-vs-LR ablation
│   └── eval_light.py                pick the best no-GPU panel (DG++-Light)
├── service/             ← CPU webservice (FastAPI + Docker; POST /predict)
└── tests/               ← pytest regression suite (no GPU, no network)
```

## The GSPA module

The predictor is wired into GSPA as `deepgo-plusplus` (legacy alias
`cafa-baseline` still works everywhere):

- **`DeepGoPlusPlusPredictor`** (`gspa-predictors/.../neural/`) — extends
  `AbstractNeuralSidecarPredictor`; delegates to the shared sidecar
  `benchmark/neural/run_neural_predictors.py --predictor deepgo-plusplus`.
- **`GspaConfig.NeuralConfig.deepGoPlusPlus`** — `{enabled, integrator,
  componentsDir, dag, batchSize, minScore}`.
- **CLI** (`annotate`): `--deepgo-plusplus`, `--deepgo-plusplus-integrator`,
  `--deepgo-plusplus-components-dir`, `--deepgo-plusplus-dag` (and the
  `--cafa-baseline*` aliases).

```bash
./gradlew :gspa-cli:run --args="annotate --input proteome.faa --output out \
  --neural-sidecar $PWD/benchmark/neural/run_neural_predictors.py \
  --deepgo-plusplus \
  --deepgo-plusplus-integrator $PWD/deepgo-plusplus/models/deepgo_plusplus_integrator_net.json \
  --deepgo-plusplus-components-dir <dir-of-component-scores> \
  --deepgo-plusplus-dag ~/Public/software/cafa6/go-dag.tsv"
```

The Groovy predictor reads the component list from the integrator JSON, so
swapping a different frozen model (e.g. the `_net` variant) needs no code change.

## Components

Each is a TSV `protein⇥term⇥score` named `<component>.tsv[.gz]` in the
components dir. Six come from GSPA's existing predictors / PLM heads (GPU/IBEX —
**produced upstream, not by this Makefile**); two (net, lit) are CPU-built here.

| component | source | built by |
|---|---|---|
| `diam`     | `DiamondPredictor` (BLAST-KNN homology) | upstream (GSPA) |
| `foldseek` | `FoldSeekPredictor` (structure-KNN)     | upstream (GSPA) |
| `clean`    | `CleanPredictor` (EC→GO)                | upstream (GSPA) |
| `interpro` | `InterProScanPredictor` (domain→GO)     | upstream (GSPA) |
| `mlp`      | ESM2-650M MLP head                      | upstream (ORIX) |
| `prostt5`  | ProstT5 structure-aware head            | upstream (ORIX) |
| `net`      | Net-KNN over STRING v12 PPI             | **`make net`** |
| `lit`      | BM25 literature text-kNN                | **`make lit`** |

Per-component no-knowledge f_w and the full ablation are in
[`TRAINING.md` §2](TRAINING.md) (performance) and `RESULTS.md` (deep dive). How
the integrator is trained and why it is leak-free: [`TRAINING.md` §1](TRAINING.md).

## No-GPU variant — DeepGO-PlusPlus-Light

DG++-Light drops the three GPU PLM heads (`mlp`, `prostt5`, `esm2_3b`) and uses
only **CPU** components — DIAMOND BLAST-KNN, FoldSeek structure-KNN (over
precomputed AlphaFold-DB structures, no GPU folding), InterProScan (raw + `_lr`),
STRING Net-KNN, optional BM25 `lit`, and a CPU **1D-CNN over sequence** (`cnn`,
DeepGOCNN-style, `pipeline/build_cnn_component.py`) intended to *replace* the PLM
heads. These are the alignment/search/lookup signals the CAFA6 winner found
*generalise best on the private test set*. No code change is needed — the
predictor reads its component list from the integrator JSON:

```bash
./gradlew :gspa-cli:run --args="annotate --input proteome.faa --output out \
  --neural-sidecar $PWD/benchmark/neural/run_neural_predictors.py \
  --deepgo-plusplus \
  --deepgo-plusplus-integrator $PWD/deepgo-plusplus/models/deepgo_plusplus_light.json \
  --deepgo-plusplus-components-dir <dir-of-CPU-component-scores> \
  --deepgo-plusplus-dag ~/Public/software/cafa6/go-dag.tsv"
```

Shipped Light models (CAFA6 no-knowledge IA-weighted f_w; all beat the full GPU
model, 0.532, on novel proteins):

| model | components | f_w | use |
|---|---|---|---|
| `deepgo_plusplus_light_cpu.json` | `diam,interpro,net_union` | **0.564** | **recommended** — strictly no-GPU (no structures) *and* works for proteins not in STRING |
| `deepgo_plusplus_light.json` | `diam,foldseek,interpro,net` | 0.550 | no-GPU *given AFDB structures* for foldseek |
| `deepgo_plusplus_light_cnn.json` | `cnn,diam,foldseek,interpro,net` | 0.516 | coverage-first: CPU 1D-CNN covers orphan proteins |

**Findings worth knowing** (full tables in `RESULTS.md`):
1. On *novel* proteins the no-GPU panels **beat** the full GPU model — the PLM
   heads are redundant with `net`.
2. **`foldseek` isn't unconditionally GPU-free**: FoldSeek search is CPU but needs
   a query *structure*; that's a CPU lookup for AFDB-covered proteins, but folding
   a brand-new sequence (ESMFold/ProstT5) needs a GPU. The strictly-no-GPU panel
   drops it (`diam+interpro+net`, 0.544 — costs only 0.006).
3. **`net` only fires for STRING members.** `net_union` (`build_net_bridge.py`)
   bridges via DIAMOND — query → pre-t0 STRING-member homolog → its neighbours'
   labels — taking the 356 no-STRING proteins from f_w 0 → **0.42** and lifting the
   structure-free panel to **0.564** (the recommended `_cpu` model). Leak-safe: the
   homolog DB is pre-t0 and novel queries can't self-match.
4. The **1D-CNN does not improve novel-protein f_w** (standalone 0.206) — same
   generalisation wall; it ships only as the coverage-first `_cnn` variant.

Re-pick / rebuild at a new release with `pipeline/eval_light.py`.

### Webservice (`service/`)

DG++-Light runs as a CPU REST API (`service/`, FastAPI + Docker): `POST /predict`
with a FASTA → JSON GO predictions, ~5 s/protein, no GPU. One DIAMOND search powers
`diam` + the bridged `net`; the 6.1 GB STRING scan is **precomputed once** into
`train_net_index.tsv` so each request is a DIAMOND search + index lookup
(`pipeline/apply_net_bridge.py` — 8 s vs ~13 min, recovering ~99 % of the slow
bridge's f_w). Default model `diam+net_union`; **`?interpro=true`** and **`?cnn=true`**
each opt in an extra component (4 frozen models cover the combinations) — `cnn`
gives a signal to orphan proteins with no homolog. See `service/README.md`.

## Reproducible retrain

> Per-release checklist and the full data-source / (re)download guide:
> [`TRAINING.md` §3–§4](TRAINING.md).

```bash
cd deepgo-plusplus
cp config.mk.example config.mk        # edit input paths + record in VERSIONS.md
uv sync                               # install the pinned environment (uv.lock)

uv run make check-inputs              # confirm every input path resolves
uv run make index                     # UniProt pass → text_string_index + net_index
uv run make gt                        # dated GOA → gt_{no,limited,partial,all}.tsv
uv run make net                       # Net-KNN component  (needs STRING_DIR)
# (place upstream diam/foldseek/clean/interpro/mlp/prostt5 under build/components/)
uv run make integrator                # train + freeze → models/deepgo_plusplus_integrator.json
uv run make eval                      # official cafaeval on the OOF predictions
```

`make all` chains `index → gt → net → integrator`. To include net in the model:
set `COMPONENT_LIST = diam,foldseek,clean,interpro,mlp,prostt5,net` in
`config.mk` (and `make lit` + append `,lit` for the literature variant).

**Why it is leak-free / CAFA-faithful.** `--model logreg --features scores` is a
linear model over component scores only — it cannot memorise term base-rates (an
XGBoost+IA/freq variant leaked to MF f_w 0.963 and was discarded). The OOF score
(`build/ltr_run/preds/ltr.tsv`, GroupKFold by protein) is the honest
generalisation estimate. All inputs are pre-t0; a blind pre-t0-validation
integrator costs only 0.005 vs the test-tuned one (`RESULTS.md`
§"Temporal integrity"). The only non-faithful element is `lit` (query-name
leakage risk), which is why it is optional.

### Ground truth

`build_groundtruth.py` expects `testsuperset_exp_annots.tsv` — the current
`goa_uniprot_all.gaf.gz` filtered to test-superset proteins + experimental
evidence codes (`accession⇥GO⇥F|P|C⇥evidence⇥YYYYMMDD`). Produce it on IBEX
(the full GAF is ~11 GB); see `RESULTS.md`. It applies t0 = 2026-02-02 and emits
the no/limited/partial knowledge-class ground truths.

### STRING download

`benchmark/neural/run_net_ws.sh` downloads the ~68 per-species STRING files with
gzip-integrity + retry (drops corrupt files). Point `STRING_DIR` at the result.

## Apply a frozen model directly (outside GSPA)

```bash
uv run python pipeline/apply_integrator.py \
  --integrator models/deepgo_plusplus_integrator_net.json \
  --components-dir build/components --dag ~/Public/software/cafa6/go-dag.tsv \
  --fasta queries.faa --out preds.tsv
```

This calls the *same* `run_deepgo_plusplus` runner GSPA uses (single source of
truth for the integration math). Add `--three-col` for cafaeval-ready output.

## Tests

```bash
uv run make test          # or: uv run pytest -q tests
```

`tests/` is a pure-CPU, no-network regression suite over tiny synthetic
fixtures: ground-truth knowledge-class logic and the t0 boundary, the net/lit
builders (incl. the literature *name-only* leak guard and corrupt-STRING skip),
the UniProt index parser, integrator training (schema + determinism), the frozen
apply (DAG propagation, sigmoid, the numerically-stable extreme-z regression),
an end-to-end chain, and a schema check over the shipped `models/*.json`. The
Groovy predictor contract is covered by `NeuralSidecarPredictorsSpec` in
`gspa-predictors`.
