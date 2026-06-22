# DeepGO-PlusPlus-Light webservice

A strictly-CPU GO-prediction REST API. One DIAMOND search of the query against the
pre-t0 train proteins powers both components of the no-GPU model:

- **`diam`** — BLAST-KNN: bit-score-weighted vote of each homolog's pre-t0 GO labels.
- **`net_union`** — homology-bridged STRING Net-KNN: vote each homolog's *precomputed*
  STRING-neighbour label vector (`train_net_index.tsv`). **No STRING files are read
  per request** (the 6.1 GB scan is precomputed once → `pipeline/apply_net_bridge.py`),
  so it serves proteins that are not in STRING at all — the realistic novel-protein
  case. Hold-out benchmark: the bridge recovers ~100 % of direct-STRING f_w
  (`../RESULTS.md`).

No GPU, no structures, nothing external needed for the default (fast) path. Two
components are **optionally** composable on top (each picks the matching frozen
model — all four combinations ship):

| `interpro` | `cnn` | model | adds |
|---|---|---|---|
| false | false | `deepgo_plusplus_light_fast.json` (default) | — |
| false | true  | `deepgo_plusplus_light_fast_cnn.json` | CPU 1D-CNN (orphan coverage) |
| true  | false | `deepgo_plusplus_light_cpu.json` | InterProScan domains |
| true  | true  | `deepgo_plusplus_light_full.json` | both |

`?cnn=true` adds the CPU 1D-CNN over sequence — it gives a signal to **orphan
proteins with no homolog** (where `diam`/`net` are blind), at the cost of bundling
PyTorch (CPU) and a `cnn_model.pt`. `?interpro=true` adds InterProScan domains
(heavy install, minutes/protein; not bundled). Either is rejected with 400 if its
asset isn't configured.

## API

```
GET  /health                      -> {status, models, interpro_available, cnn_available, n_homolog_nodes}
POST /predict  (body: FASTA text) -> {model, n_proteins, predictions: {prot: [{term,name,aspect,score}]}}
     ?interpro=true   add InterProScan domains   (needs DGPP_INTERPROSCAN)
     ?cnn=true        add the CPU 1D-CNN          (needs cnn_model.pt / DGPP_CNN_MODEL)
     ?min_score=0.1   score cutoff
     ?topk=5          homologs voted per query
```

```bash
curl -X POST 'http://localhost:8000/predict?min_score=0.4' \
     -H 'Content-Type: text/plain' --data-binary @query.faa
```

## Run with Docker

The image bundles the code, the frozen model JSONs and the DIAMOND binary. The
large data assets are mounted at `/assets` (built once, reused across releases).

```bash
# 1. build the asset bundle (DIAMOND DB + the precomputed bridge index + dag/obo/labels)
#    train_net_index.tsv is the one-time ~45 min precompute (see make_assets.sh).
service/make_assets.sh /data/dgpp-assets \
    train.fasta train_net_index.tsv train_terms.tsv go-dag.tsv go.obo

# 2. build the image (context = the deepgo-plusplus/ dir)
docker build -f service/Dockerfile -t dgpp-light .

# 3. run (mount the assets read-only)
docker run -d -p 8000:8000 -v /data/dgpp-assets:/assets:ro dgpp-light
```

### Enabling the optional components

- **`?cnn=true`** — train the CNN once with weights saved, drop the `.pt` into the
  asset bundle (`make_assets.sh ... go.obo cnn_model.pt`), and it is auto-detected at
  `/assets/cnn_model.pt` (or set `DGPP_CNN_MODEL`). PyTorch (CPU) is already in the
  image; the model is lazy-loaded only on a `cnn=true` request.
  ```bash
  build_cnn_component.py --train-fasta train.fasta --train-terms train_terms.tsv \
      --dag go-dag.tsv --test-fasta tiny.faa --out /dev/null --save-model cnn_model.pt
  ```
- **`?interpro=true`** — InterProScan is large (~tens of GB) and slow, so it is
  **not** bundled. Make `interproscan.sh` reachable in the container and set
  `-e DGPP_INTERPROSCAN=/opt/interproscan/interproscan.sh`.

## Environment

| var | default | meaning |
|---|---|---|
| `DGPP_ASSETS` | `/assets` | dir with `train_db.dmnd`, `train_net_index.tsv`, `train_terms.tsv`, `go-dag.tsv`, `go.obo` |
| `DGPP_MODELS` | `/app/models` | dir with the frozen integrator JSONs |
| `DGPP_THREADS` | `8` | DIAMOND / InterProScan threads |
| `DGPP_INTERPROSCAN` | — | path to `interproscan.sh` to enable `?interpro=true` |
| `DGPP_CNN_MODEL` | `/assets/cnn_model.pt` | saved CNN weights to enable `?cnn=true` |
| `DGPP_DIAMOND` | `diamond` | DIAMOND binary |

## Notes

- Latency is DIAMOND-bound (~5 s/protein at `--very-sensitive`); batch many
  sequences per request to amortise. The index loads once at startup (~a few s).
- The default model is novel-protein-first (CAFA6 no-knowledge) and needs *some*
  homolog. For true no-homolog **orphans**, use `?cnn=true` — the CPU 1D-CNN gives a
  sequence-based signal where `diam`/`net` are blind (it lowers mean f_w on the
  homolog-rich benchmark but is the only signal for orphans).
- Rebuild assets when the UniProt/STRING/CAFA release changes (`../TRAINING.md`).
