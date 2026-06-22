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

No GPU, no structures, no InterProScan needed for the default (fast) path.
`?interpro=true` switches to the 3-component model (`diam+interpro+net_union`,
f_w 0.564) but requires an InterProScan install (heavy; not bundled).

## API

```
GET  /health                      -> {status, fast_model, full_model, interpro_available, n_homolog_nodes}
POST /predict  (body: FASTA text) -> {model, n_proteins, predictions: {prot: [{term,name,aspect,score}]}}
     ?interpro=true   use the +InterProScan model (needs DGPP_INTERPROSCAN)
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

### Enabling the `interpro=true` path

InterProScan is large (~tens of GB) and slow (minutes/protein), so it is **not**
bundled. To enable it, make `interproscan.sh` reachable in the container and point
`DGPP_INTERPROSCAN` at it (e.g. mount it and `-e DGPP_INTERPROSCAN=/opt/interproscan/interproscan.sh`).

## Environment

| var | default | meaning |
|---|---|---|
| `DGPP_ASSETS` | `/assets` | dir with `train_db.dmnd`, `train_net_index.tsv`, `train_terms.tsv`, `go-dag.tsv`, `go.obo` |
| `DGPP_MODELS` | `/app/models` | dir with the frozen integrator JSONs |
| `DGPP_THREADS` | `8` | DIAMOND / InterProScan threads |
| `DGPP_INTERPROSCAN` | — | path to `interproscan.sh` to enable `?interpro=true` |
| `DGPP_DIAMOND` | `diamond` | DIAMOND binary |

## Notes

- Latency is DIAMOND-bound (~5 s/protein at `--very-sensitive`); batch many
  sequences per request to amortise. The index loads once at startup (~a few s).
- The model is novel-protein-first (CAFA6 no-knowledge). It needs *some* homolog;
  true no-homolog orphans get little signal (the offline `cnn` component is the
  coverage fallback — not wired into the service).
- Rebuild assets when the UniProt/STRING/CAFA release changes (`../TRAINING.md`).
