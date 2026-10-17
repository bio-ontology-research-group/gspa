# Web service

Use this workflow to deploy and query the GSPA FastAPI web service (`service/`).
The web service provides a stateless REST interface wrapping `gspa annotate`,
enabling web frontends (such as DeepGOWeb's "Genome" tab) or automated workflows
to submit genome assemblies and receive calibrated functional annotations,
per-contig quality scores, and symbolic enforcement audit logs.

The service is deliberately stateless: each request executes in an isolated
temporary directory, parses the emitted TSV and JSON files into a structured
response, and discards temporary files. Persistence, scheduling, and job queues
can be handled externally (for example using Celery or Redis).

## What you need

| Requirement | Description |
|---|---|
| Container runtime | Docker, Podman, or Apptainer |
| DeepGO-PlusPlus-Light assets | Shared asset directory containing `train_db.dmnd`, `train_net_index.tsv`, `train_terms.tsv`, `go-dag.tsv`, `go.obo` |
| System resources | 4+ CPU cores and 8+ GB RAM (CPU-only, no GPU required) |

## REST API specification

### 1. Health and readiness check (`GET /health`)

Verifies that the GSPA binary is available, the neural sidecar is present, and
the reference asset bundle is mounted.

```bash
curl http://localhost:8000/health
```

Response JSON:

```json
{
  "status": "ok",
  "gspa_bin": "/usr/local/bin/gspa-cli",
  "assets_dir": "/opt/dgpp_assets",
  "assets_present": true,
  "go_owl": "/opt/dgpp_assets/go.obo",
  "go_owl_present": true,
  "neural_sidecar": "/app/benchmark/neural/run_neural_predictors.py",
  "neural_sidecar_present": true
}
```

### 2. Genome annotation (`POST /annotate`)

Accepts multipart form-data containing the genome assembly, paired GFF3, and
configuration parameters.

#### Form fields

| Field | Type | Default | Description |
|---|---|---|---|
| `genome` | File | none | Genome or metagenome nucleotide FASTA |
| `gff3` | File | none | Paired structural GFF3; existing CDS features are translated |
| `proteins` | File | none | Pre-called protein FASTA (alternative to `genome`) |
| `predictor` | String | `light` | Base predictor: `light` (CPU), `full`, or `none` |
| `metrics_scope` | String | `contig` | Evaluation scope: `contig`, `genome`, or `both` |
| `kingdom` | String | auto | Organism kingdom: `bacteria`, `archaea`, `eukaryote`, `virus` |
| `mag` | Boolean | `false` | Identifies input as a MAG, adjusting completeness thresholds |
| `enforce_consistency` | Boolean | `false` | Runs Stage 2A SAT taxon-consistency enforcement |
| `consistency_mode` | String | `remove` | Action for inconsistent calls: `remove`, `downrank`, `flag`, `minimal-flip` |
| `taxon` | String | none | Asserted organism taxon (e.g. `bacteria`, `NCBITaxon_562`) |
| `enforce_completeness`| Boolean | `false` | Promotes missing essential functions onto best-evidenced genes |
| `enforce_coherence` | Boolean | `false` | Prunes orphan subunits of obligate heteromers |
| `provenance` | Boolean | `true` | Records per-call predictor evidence in the output |
| `threshold` | Float | `0.50` | Posterior probability floor for emitted annotations |

#### Response schema

The service returns a consolidated JSON object:

- `annotations`: array of functional calls matching `<genome>_annotations.tsv`.
- `per_contig_metrics`: table of completeness and coherence metrics per contig.
- `enforcement_actions`: audit log of terms relaxed or removed by Stage 2A.
- `inferred_taxon`: lineage inferred from sequence evidence if `--taxon` was omitted.
- `timing`: wall-clock seconds spent in prediction, ontology loading, and metrics.

## Deploying with Docker

### Step 1: build the image

Build the container image from the GSPA repository root:

```bash
docker build -f service/Dockerfile -t gspa-service .
```

### Step 2: run the container

Mount your DeepGO-PlusPlus-Light asset directory into `/opt/dgpp_assets`:

```bash
docker run -d \
  --name gspa-api \
  -p 8000:8000 \
  -v /data/reference/dgpp_assets:/opt/dgpp_assets:ro \
  -e GSPA_THREADS=8 \
  -e GSPA_TIMEOUT=1800 \
  gspa-service
```

### Step 3: deployment with Docker Compose

For automated service orchestration, create `docker-compose.yml`:

```yaml
version: '3.8'

services:
  gspa-service:
    build:
      context: .
      dockerfile: service/Dockerfile
    image: gspa-service:latest
    container_name: gspa-api
    ports:
      - "8000:8000"
    volumes:
      - /data/reference/dgpp_assets:/opt/dgpp_assets:ro
    environment:
      - GSPA_THREADS=16
      - GSPA_TIMEOUT=3600
      - GSPA_MIN_SCORE=0.50
      - GSPA_MAX_UPLOAD_BYTES=209715200
    restart: unless-stopped
```

Start the service:

```bash
docker compose up -d
```

## Querying the service

### Using cURL

```bash
curl -X POST http://localhost:8000/annotate \
  -F genome=@genome.fna \
  -F gff3=@annotation.gff3 \
  -F predictor=light \
  -F metrics_scope=contig \
  -F taxon=bacteria \
  -F enforce_consistency=true \
  -F consistency_mode=remove \
  -o response.json
```

### Using Python

```python
import requests

url = "http://localhost:8000/annotate"

files = {
    "genome": open("contigs.fna", "rb"),
    "gff3": open("annotation.gff3", "rb"),
}

data = {
    "predictor": "light",
    "metrics_scope": "contig",
    "taxon": "NCBITaxon_562",
    "enforce_consistency": "true",
    "consistency_mode": "remove",
    "threshold": "0.50",
}

response = requests.post(url, files=files, data=data)
result = response.json()

print(f"Annotated {len(result['annotations'])} protein-function pairs.")
if result.get("enforcement_actions"):
    print(f"Removed {len(result['enforcement_actions'])} taxon-inconsistent terms.")
```

## Environment configuration

Configure service behavior using environment variables:

| Variable | Default | Purpose |
|---|---|---|
| `GSPA_BIN` | `gspa-cli` | Path to the GSPA launcher binary |
| `DGPP_ASSETS` | `/opt/dgpp_assets` | Directory containing DeepGO-PlusPlus-Light assets |
| `NEURAL_SIDECAR` | `/app/benchmark/neural/run_neural_predictors.py` | Python sidecar for neural base predictors |
| `GSPA_GO_OWL` | `$DGPP_ASSETS/go.obo` | Ontology file used for metrics and enforcement |
| `GSPA_THREADS` | `0` | Execution threads for DIAMOND and pipeline (0 = auto-detect) |
| `GSPA_TIMEOUT` | `3600` | Maximum wall-clock seconds allowed per request |
| `GSPA_MIN_SCORE` | `0.50` | Default score threshold for emitted annotations |
| `GSPA_MAX_UPLOAD_BYTES` | `209715200` | Maximum upload size in bytes (200 MB) |

## Warm inference and performance

The service operates entirely on CPU. Processing time is dominated by DIAMOND
alignment against the reference database:

- An average bacterial isolate (4,000 proteins) completes in 20 to 45 seconds
  on an 8-core CPU.
- METEOR and DeepGO-PlusPlus-Light models remain warm on disk through Linux
  page caching, eliminating model reloading latency across requests.
- For high-volume production deployments, run multiple container replicas behind
  a reverse proxy (such as Nginx, Traefik, or an AWS Application Load Balancer).

## Pitfalls

- **Payload upload limits**: attempting to upload multi-gigabyte metagenomes
  exceeding `GSPA_MAX_UPLOAD_BYTES` returns an HTTP 413 error. Increase the limit
  via the environment variable if processing large metagenomes.
- **Unmounted asset directory**: if `train_db.dmnd` is missing from `DGPP_ASSETS`,
  the `/health` endpoint reports `assets_present: false` and annotations fail.
  Verify volume mount permissions.
- **Client timeout errors**: because genome-scale annotation takes 30 to 60
  seconds, ensure HTTP reverse proxies and client libraries configure read
  timeouts of at least 300 seconds.
