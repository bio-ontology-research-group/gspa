# GSPA genome-scale annotation service

A small, stateless REST wrapper around `gspa-cli annotate`. It exposes the JVM
genome-scale pipeline — **GFF3 CDS translation → DeepGO-PlusPlus(-Light)
prediction → per-contig genome-scale metrics → optional SAT taxon-consistency /
completeness / coherence enforcement → provenance** — so a web frontend (e.g.
DeepGOWeb's "Genome" tab) can drive it.

This is the genome-scale counterpart to `deepgo-plusplus/service` (which serves
per-protein DG++Light only): here the unit of work is a **genome / metagenome
FASTA + GFF3**, and the response carries the **per-contig quality table** and the
**enforcement / provenance** trail, not just a flat GO list.

## API

```
GET  /health    -> readiness (gspa binary, assets, ontology, sidecar present?)
POST /annotate  (multipart/form-data) -> one genome-scale annotation as JSON
```

`POST /annotate` fields:

| field | type | default | meaning |
|---|---|---|---|
| `genome` | file | — | genome / metagenome nucleotide FASTA (one contig per sequence) |
| `gff3` | file | — | GFF3 paired with the genome; CDS are translated, not re-called |
| `proteins` | file | — | pre-called protein FASTA (alternative to `genome`) |
| `predictor` | str | `light` | `light` (CPU, self-contained) \| `full` \| `none` |
| `metrics_scope` | str | `contig` | `contig` (per contig) \| `genome` (pooled) \| `both` |
| `kingdom` | str | auto | `bacteria` \| `archaea` \| `eukaryote` \| `virus` |
| `mag` | bool | false | input is a MAG (adjusts quality thresholds) |
| `enforce_consistency` | bool | false | run the SAT taxon-constraint enforcement pass |
| `consistency_mode` | str | `remove` | `remove` \| `downrank` \| `flag` \| `minimal-flip` |
| `taxon` | str | — | assert the organism taxon (`bacteria`/`NCBITaxon_2`/...) for precise per-term removal |
| `enforce_completeness` | bool | false | promote missing essential functions onto the best-evidenced protein |
| `enforce_coherence` | bool | false | fix complex singletons + missing has_part partners (needs ELK) |
| `provenance` | bool | true | record the provenance trail + enforcement-actions log |

Response JSON: `{ ok, predictor, metrics_scope, annotations[], per_contig_metrics[],
enforcement_actions[], quality_json_files[], log }`, where `annotations` mirrors
`<genome>_annotations.tsv` (incl. the `provenance` column) and
`per_contig_metrics` mirrors `<genome>_quality_per_contig.tsv`.

```bash
curl -X POST http://localhost:8000/annotate \
     -F genome=@contigs.fna -F gff3=@annotation.gff3 \
     -F predictor=light -F metrics_scope=contig \
     -F enforce_consistency=true -F taxon=bacteria
```

## Run with Docker

```bash
# context = the gspa repo root (needs gspa-*/, benchmark/, deepgo-plusplus/)
docker build -f service/Dockerfile -t gspa-service .

# mount the DG++Light asset bundle (same one DeepGOWeb downloads)
docker run -d -p 8000:8000 -v /data/dgpp-assets:/opt/dgpp_assets:ro gspa-service
```

## Environment

| var | default | meaning |
|---|---|---|
| `GSPA_BIN` | `gspa-cli` | the gspa launcher (on PATH in the image) |
| `DGPP_ASSETS` | `/opt/dgpp_assets` | DG++Light bundle (`train_db.dmnd`, `train_net_index.tsv`, `train_terms.tsv`, `go-dag.tsv`, `go.obo`, ...) |
| `NEURAL_SIDECAR` | `/app/benchmark/neural/run_neural_predictors.py` | sidecar gspa-cli shells out to for DG++Light |
| `GSPA_GO_OWL` | `$DGPP_ASSETS/go.obo` | ontology enabling genome-scale metrics + enforcement |
| `GSPA_THREADS` | `0` | DIAMOND / pipeline threads (0 = gspa picks) |
| `GSPA_TIMEOUT` | `3600` | per-request annotate ceiling (seconds) |
| `GSPA_MAX_UPLOAD_BYTES` | `209715200` | upload size cap (200 MB) |

## Notes

- **CPU-only.** No GPU needed; the same profile DeepGOWeb runs. Latency is
  DIAMOND-bound — batch a whole genome per request, not protein-by-protein.
- **Stateless.** Each request runs in a throwaway temp dir; queueing /
  persistence is the caller's job (DeepGOWeb uses Celery).
- `go.obo` shipped in the asset bundle carries no taxon-constraint axioms, so
  `enforce_consistency` without `taxon` is a safe no-op on it; the bundled NCBI
  constraints + `--taxon bacteria` are what make removal act (see
  `../GENOME_GFF3_ANNOTATION.md`).
