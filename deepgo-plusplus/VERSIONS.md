# DeepGO-PlusPlus — input release provenance

Record the exact external release used for each frozen model so a retrain is
reproducible and auditable. Append a new row whenever you re-freeze.

## Shipped models (`models/`)

> **⚠️ f_w corrected 2026-06-23.** The original `OOF no-knowledge f_w` column was on
> the GAF-date-**contaminated** GT (see [`TRAINING.md` §1.0a](TRAINING.md)). Both the
> `dirty` and `clean` (leak-free, `build_clean_gt.py`) values are shown; **compare on
> `clean`**. The 3-class mean is **not** CAFA-LB-comparable and was dropped. The
> Light models' "beats GPU" was leak — see [`RESULTS.md`](RESULTS.md).

| model JSON | components | STRING | GO ontology | t0 | no-knowledge f_w (dirty → **clean**) |
|---|---|---|---|---|---|
| `deepgo_plusplus_integrator.json`          | diam,foldseek,clean,interpro,mlp,prostt5 | — | go-basic 2025-10-10 | 2026-02-02 | 0.483 → ~0.49 |
| `deepgo_plusplus_integrator_net.json`      | + net (STRING v12.0, 2023) | v12.0 (2023) | go-basic 2025-10-10 | 2026-02-02 | 0.532 → **0.521** |
| `deepgo_plusplus_integrator_lit_net.json`  | + lit + net | v12.0 (2023) | go-basic 2025-10-10 | 2026-02-02 | 0.553 → (re-eval) |
| `deepgo_plusplus_light.json`               | diam,foldseek,interpro,net | v12.0 (2023) | go-basic 2025-10-10 | 2026-02-02 | 0.550 → **0.488** (best Light on clean) |
| `deepgo_plusplus_light_cnn.json`           | + cnn (CPU 1D-CNN) | v12.0 (2023) | go-basic 2025-10-10 | 2026-02-02 | 0.516 → **0.470** |
| `deepgo_plusplus_light_cpu.json`           | diam,interpro,net_union | v12.0 (2023) | go-basic 2025-10-10 | 2026-02-02 | 0.564 → **0.464** (most leak-inflated) |
| **`deepgo_plusplus_light_clean.json`** *(net-free, retrained leak-free)* | diam,foldseek,interpro,lit | — | go-basic 2025-10-10 | 2026-02-02 | **0.508** (clean OOF; needs structures for foldseek) |
| **`deepgo_plusplus_light_cpu_clean.json`** *(net-free, strictly-CPU)* | diam,interpro,lit,cnn | — | go-basic 2025-10-10 | 2026-02-02 | **0.500** (clean OOF; any sequence) |
| **`deepgo_plusplus_integrator_tierA.json`** *(genome cascade, homology tier)* | diam,net_union,interpro | v12.0 (2023) | go-basic 2025-10-10 | 2026-02-02 | **0.509** (clean, Tier-A homology proteins) |
| **`deepgo_plusplus_integrator_tierB.json`** *(genome cascade, orphan tier)* | esm2_knn (ESM2-35M kNN) | — | go-basic 2025-10-10 | 2026-02-02 | **0.508** (clean, orphan/no-homolog proteins) |

The genome-cascade tier models (see [`CASCADE.md`](CASCADE.md)) are applied per-protein
by `service/predict.py::cascade()`: one DIAMOND search triages the proteome — homolog →
Integrator-A, orphan → Integrator-B. **Tier-B requires the ESM2-35M train embedding store**
(`emb/train_esm2_35m.npz`, 82,404 × 480 fp16 = the CPU kNN reference; gitignored scratch,
ship as a deployment asset). Rebuild on ORIX:
`extract_embeddings.py --model esm2_35m --fasta train_cascade.fasta` (~5 min, 1×H100), then
`pipeline/build_esm2_knn.py` for the LOO train feature / test transfer.

Reference: full GPU integrator (6+net) = 0.532 → **0.521** clean; it **beats** every
Light panel on the leak-free GT. GOAlpha (CAFA6 1st) no-knowledge ≈ 0.524.

The `cnn` component for the Light-CNN model is a CPU 1D-CNN trained on the pre-t0
SwissProt sequences (`pipeline/build_cnn_component.py`; train FASTA extracted from
`uniprot_sprot.dat.gz` via `extract_sprot_fasta.py`). It is not shipped as weights;
rebuild it (or save weights with `--save-model`) at each release.

The `net_union` component (Light-CPU model) = plain `net` for STRING members +
`build_net_bridge.py` for the rest (DIAMOND vs the pre-t0 train DB → STRING-member
homolog → neighbour labels). Both halves are pre-t0; novel queries can't self-match
the train DB. Rebuild at each release alongside `net`.

> Fill the `2025_xx` SwissProt release placeholders with the precise release the
> next retrain uses (e.g. `2025_04`). All ontology/CAFA inputs are pre-t0, so the
> models are CAFA-submission-faithful (see README §"Temporal integrity").

## Temporal-integrity invariants (must hold at every retrain)

- GO ontology release **on or before t0** (2026-02-02).
- `train_terms.tsv` / `IA.tsv` are the CAFA6 official pre-t0 artifacts.
- STRING release predates t0 (v12.0 is 2023).
- The literature component's *query* text is name-only (identification fields),
  never post-t0 `CC FUNCTION` — see `pipeline/build_lit_component.py`.
