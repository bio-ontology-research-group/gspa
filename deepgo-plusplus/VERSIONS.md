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
| **`deepgo_plusplus_integrator_full_aux.json`** *(full + aux components)* | diam,foldseek,clean,interpro,mlp,prostt5,net,lit,**proteinfer,eggnog,deepfri** | v12.0 (2023) | go-basic 2025-10-10 | 2026-02-02 | **0.541** (clean-A; aux add only +0.001 over the lit-full base 0.540) |
| **`deepgo_plusplus_integrator_cpu_aux.json`** *(CPU-only + aux, novel-computable)* | diam,interpro,cnn,net_union,esm2_knn,**proteinfer,eggnog,deepfri** | v12.0 (2023) | go-basic 2025-10-10 | 2026-02-02 | **0.530** (clean-A; CPU base 0.517 → +0.013, ≈ all from proteinfer; beats old full GPU 0.521) |
| **`deepgo_plusplus_integrator_cpu_lean_mcm.json`** *(hierarchy-aware CPU)* | diam,interpro,**cnn(MCM)**,net_union,esm2_knn,proteinfer | v12.0 (2023) | go-basic 2025-10-10 | 2026-02-02 | **0.524** (clean-A; cnn head retrained with C-HMCNN over is_a∪part_of; bce baseline 0.519) |
| **`deepgo_plusplus_integrator_full_aux_mcm.json`** *(hierarchy-aware FULL)* | diam,foldseek,clean,interpro,**mlp(650M,MCM)**,**prostt5(MCM)**,net,lit,proteinfer,eggnog,deepfri | v12.0 (2023) | go-basic 2025-10-10 | 2026-02-02 | **0.545** (clean-A; all trainable PLM heads retrained with C-HMCNN; matched bce baseline 0.539) |

### Hierarchy-aware heads (C-HMCNN, is_a ∪ part_of) — 2026-06-24

The `_mcm` integrators use components whose **trainable heads were retrained with a
hierarchy-aware loss**: the C-HMCNN Max-Constraint Module (Giunchiglia & Lukasiewicz,
NeurIPS 2020) over the GO **is_a ∪ part_of** DAG, replacing flat BCE. Implementation:
`pipeline/train_head_hmcnn.py` (`--loss mcm`, scalable `scatter_reduce(amax)` MCM — no
(B,n,n) blow-up) for the PLM heads; `pipeline/build_cnn_component.py --loss mcm` for the
CPU 1D-CNN. **Standalone f_w gain (clean-A / clean-B):** prostt5 +0.006/+0.008, esm2_3b
+0.008/+0.008, ESM2-650M +0.008/+0.015, **cnn +0.037/+0.041** (the weakest head gains most
— it leans hardest on the DAG prior). **Integrated:** FULL +0.005 (prostt5 swap) to +0.009
(all heads); CPU +0.005 (cnn swap). The hard max-constraint beat a soft true-path penalty
(`--loss softreg`). Gains concentrate in BP/CC (deep DAG); MF ≈ flat. **Deployment weights**
(gitignored, `models/weights/`, rebuild at release): `cnn_mcm.pt` (DG++-Light cnn component,
loaded by `service/predict.py`), `head_prostt5_mcm.pt` + `head_650m_mcm.pt` (FULL PLM heads,
applied via the `dgpp-head` sidecar runner → `extract_embeddings.py` then
`train_head_hmcnn.py --load-model`). **Pending the existing IBEX dependency:** re-freezing
the genome-cascade *tier* integrators (tierA/tierB) on the pre-t0 orphan population with the
MCM components (same blocker as the production tier re-freeze).

The `_aux` integrators fold in the benchmarked auxiliary components (see
[`CASCADE.md`](CASCADE.md) "Auxiliary components"). **Honest attribution (integrator-
inclusion ablation):** `proteinfer` is the only meaningful lift (+0.013 in the CPU model);
`eggnog`+`deepfri` add +0.002 together (redundant with diam/esm2_knn and cnn); in the FULL
model the aux add only +0.001 (PLM heads already capture it). `deepfri` here is **seq-mode
(CPU, no structure)** — the structural GCN (GPU) is untested. `psortb` excluded
(bacterial-only → kingdom-gate). The CPU-only model still exceeds the previous full GPU
integrator — driven by `esm2_knn` + `proteinfer`, not eggnog/deepfri.

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
