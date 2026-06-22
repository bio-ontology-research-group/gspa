# DeepGO-PlusPlus — input release provenance

Record the exact external release used for each frozen model so a retrain is
reproducible and auditable. Append a new row whenever you re-freeze.

## Shipped models (`models/`)

| model JSON | components | UniProt (SwissProt) | STRING | GO ontology | CAFA train_terms / IA | t0 | OOF no-knowledge f_w |
|---|---|---|---|---|---|---|---|
| `deepgo_plusplus_integrator.json`          | diam,foldseek,clean,interpro,mlp,prostt5 | 2025_xx (pre-t0) | — | go-basic 2025-10-10 | CAFA6 official (pre-t0) | 2026-02-02 | 0.483 |
| `deepgo_plusplus_integrator_net.json`      | + net (STRING v12.0, 2023) | 2025_xx | v12.0 (2023) | go-basic 2025-10-10 | CAFA6 official | 2026-02-02 | 0.538 (3-class 0.647) |
| `deepgo_plusplus_integrator_lit_net.json`  | + lit + net | 2025_xx | v12.0 (2023) | go-basic 2025-10-10 | CAFA6 official | 2026-02-02 | 0.553 |
| `deepgo_plusplus_light.json`               | diam,foldseek,interpro,net (no GPU) | 2025_xx | v12.0 (2023) | go-basic 2025-10-10 | CAFA6 official | 2026-02-02 | **0.550** |
| `deepgo_plusplus_light_cnn.json`           | + cnn (CPU 1D-CNN) | 2025_xx | v12.0 (2023) | go-basic 2025-10-10 | CAFA6 official | 2026-02-02 | 0.516 |

The `cnn` component for the Light-CNN model is a CPU 1D-CNN trained on the pre-t0
SwissProt sequences (`pipeline/build_cnn_component.py`; train FASTA extracted from
`uniprot_sprot.dat.gz` via `extract_sprot_fasta.py`). It is not shipped as weights;
rebuild it (or save weights with `--save-model`) at each release.

> Fill the `2025_xx` SwissProt release placeholders with the precise release the
> next retrain uses (e.g. `2025_04`). All ontology/CAFA inputs are pre-t0, so the
> models are CAFA-submission-faithful (see README §"Temporal integrity").

## Temporal-integrity invariants (must hold at every retrain)

- GO ontology release **on or before t0** (2026-02-02).
- `train_terms.tsv` / `IA.tsv` are the CAFA6 official pre-t0 artifacts.
- STRING release predates t0 (v12.0 is 2023).
- The literature component's *query* text is name-only (identification fields),
  never post-t0 `CC FUNCTION` — see `pipeline/build_lit_component.py`.
