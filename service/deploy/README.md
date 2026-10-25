# GSPA public prediction service

`https://gspa.bio2vec.net` runs on phenotype-reactor (`10.254.146.211:18080`).
Public nginx on borg-server2 terminates TLS and proxies to
`http://phenomebrowser.net/gspa/`. Both frontend load balancers route `/gspa/`
to the VM. They do not run prediction code.

The deployed application is `service.portal:app`, selected with `GSPA_APP`.
The existing synchronous `service.app:app` remains available for internal API
integrations, but is not exposed by this deployment.

## Prediction and results

Every new submission selects the full profile by default. Individual analyses
can be disabled in the upload form. The API uses the same defaults; `disabled`
is an optional comma-separated list of analysis keys.

GFF3 is the recommended input: embed the genome after `##FASTA`, or upload a
matching companion FASTA. Existing CDS models, identifiers, genes and products
are retained. Unannotated prokaryotic genomes use Prodigal. Protein-only inputs
receive function, homology, AMR and annotation-quality analyses; analyses that
require genomic coordinates are explicitly marked inapplicable.

The workflow runs DeepGO++ Light (GO, score floor 0.5), DIAMOND against reviewed
Swiss-Prot (GO/EC and optional product names), candidate operons, KEGG map/module
coverage, AMRFinderPlus, antiSMASH, PSORTb, GSPA annotation quality and bacterial
metabolic reconstruction/FBA. All requested tools must complete; a tool failure
fails the job rather than being presented as zero hits. Supported bacterial
biomass templates are Gram negative/positive. Other domains explicitly report
that metabolic reconstruction needs an appropriate organism model. PSORTb
requires bacteria or archaea. Its terse formatter retains tied best locations
instead of choosing one according to Perl hash iteration order. BGC detection requires prokaryotic coordinates.

Homology product assignment requires ≥40% identity and ≥80% coverage on both
query and subject, and fills missing/generic names. Input products take priority.
Candidate operons use consecutive, same-strand CDS separated by ≤300 bp. Their
labels summarize functional evidence; they are not curated operon names.
KEGG coverage counts distinct EC-mapped GO functions and is not a Boolean
assessment of alternative module branches. AMRFinderPlus runs protein searches;
organism-specific resistance mutations are not tested.

Quality measures annotation completeness, ontology/taxon consistency and
functional coherence. CheckM/BUSCO assembly-quality metrics are not inferred.
FBA uses the exact ModelSEED stoichiometry, Gram-specific biomass, and a named
medium supplemented by the reference's basal uptake set. Download SBML and the
reaction table for exchange bounds, candidate GPRs, fluxes and repairs. Feasible
biomass flux is a conditional model result, not an experimental growth rate.

Results include the standalone HTML browser, original inputs, annotations,
integrated scores, provenance/actions, raw tool outputs, GAF quality input,
SBML/fluxes, configuration, run log, and version/hash manifest. The two reference
mounts contain DG++ Light assets and the full analysis reference collection.

## Operations

- Build the tools image with `docker build -f service/Dockerfile.tools -t gspa-analysis-tools:2026-10-25 .`, then the service with `docker build -f service/Dockerfile -t gspa-service:full-2026-10-25 .`.
- Build PSORTb with `docker build -f service/Dockerfile.psortb -t gspa-psortb:2026-10-25 .`. This private HTTP sidecar has no host port and no Docker socket.
- Configuration: `compose.yaml`, installed as `/opt/gspa/compose.yaml`.
- Persistent references: `/opt/gspa/assets` and `/opt/gspa/resources` (read-only in the containers).
- Persistent jobs: `/opt/gspa/jobs`, owned by UID/GID 10001.
- Limits: 4 CPUs / 12 GiB for GSPA plus 4 CPUs / 3 GiB for PSORTb (four independent batches); one running job, ten active jobs including the queue.
- Uploads: 20 MB per file (41 MB total proxy request limit for paired GFF3/FASTA), at most 20,000 sequences; individual proteins at most 20,000 residues.
- Refuse new uploads below 8 GiB free space. Completed/failed jobs expire after seven days.
- Run exactly one uvicorn process: the worker is in-process, and this is not a distributed queue.
- Queued jobs survive restarts. Interrupted jobs are marked failed and can be resubmitted.
- Each stage has a one-hour timeout. The warm prediction server is disabled so
  timeout/shutdown can terminate the complete job process group.
- Non-root container, read-only root filesystem, no capabilities, 2 GiB temporary filesystem.
- Docker logs rotate at 3 × 10 MB. Container restart policy is `unless-stopped`.
- Health: `GET /health` checks the CLI, DG++ assets and the full analysis reference paths.
- Inspect: `docker logs --tail 100 gspa`, `docker stats --no-stream gspa`.
- This VM has Compose v1: `docker-compose -f /opt/gspa/compose.yaml up -d`.

The public service has no accounts. A random 128-bit job URL grants access to
its files. Input data is not listed publicly, but links can be shared; use
public/non-sensitive sequences. No raw upload filenames become filesystem paths.

## API and verification

- `POST /api/jobs`: multipart `file`, `kind` (`gff3` default / `genome` / `proteins`), optional companion `genome`, `preset`, `gram` (`neg` default / `pos`), optional `disabled`.
- `GET /api/jobs/{id}`: status, stage, counts, timestamps.
- `GET /api/jobs/{id}/browser`, `/download`, `/log`, `/model`, `/fluxes`, `/analyses`: allowlisted artifacts.
- Presets: bacteria, archaea, prokaryote, eukaryote, mag, metagenome, virus.

Validation tests:
`uv run --with-requirements service/requirements.txt --with httpx python -m unittest service.test_portal`.
The ClaimExtractor regression tests cover DeepGO++ source recognition.
End-to-end deployment validation must upload actual sequences, require nonzero
predictions for a known positive sample, and open both the browser and ZIP.
HTTP health alone is not proof of working inference.

## Rollback

Stop only GSPA with `docker stop gspa`; its input/reference volumes remain.
Restore the previous GSPA image and recreate its container when rolling back a
release. Existing phenotype-reactor containers are not part of this Compose
project. Proxy backups on both frontends:
`/root/phenomebrowser.net.before-gspa-2026-10-25`. Prefer removing only the GSPA
location blocks if unrelated routing has changed since that backup. Run
`nginx -t` before every reload. The public vhost and automatically renewed
certificate are named `gspa.bio2vec.net`.

## Reference preparation and pinned examples

The full reference tree contains `diamond/swissprot.dmnd`, `go/go.owl`,
`ec2go/ec2go`, `kegg_pathways.tsv`, `kegg_modules.tsv`, the ModelSEED `fba_data/`
files, and versioned AMRFinderPlus/antiSMASH data under `tool-data/`. Prepare the
Swiss-Prot database with `scripts/gspa_resources.py` so its headers retain
products, gene names, GO and EC evidence. Download/index antiSMASH data before
starting the read-only service. All reference paths must be readable by UID 10001.

Precompute the full GO reasoning cache once using `gspa-cli evaluate` with the
installed `--go-owl` and `--reasoner-cache /opt/gspa-resources/reasoner-cache`,
plus a representative GFF3 and GAF. The expensive `has_part` extraction covers
the entire ontology and is independent of that sample. The cache includes the
ontology SHA-256; a changed ontology invalidates it. Regenerate it when updating
GO, then mount the reference tree read-only. Readiness requires the cache file.
Do not copy a cache from a different GO release. Per-job quality assessment still
runs on that job's integrated annotations.

FBA runs in `/opt/fba` with NumPy 2.2.6, SciPy 1.15.3 and libSBML 5.21.1,
separately from antiSMASH's environment. Each bacterial model is tested with
22 medium conditions, including individual amino-acid uptake closures and a
joint closure of the 20 listed amino acids. These tests use the exported network
without further gap filling. The archive includes `metabolism/conditions.json`.

Curated public examples use `/opt/gspa/jobs/examples.json`, mapping a name such
as `ecoli` to a completed job ID. Mark that job's status with `"example": true`
to exempt it from routine expiration. `/examples/ecoli` serves its browser;
`/examples/ecoli/annotation`, `/sequence`, `/model`, `/fluxes` and `/download`
serve the corresponding original inputs and results. These aliases only expose
the same allowlisted artifacts as ordinary job routes.


## Full-profile validation, 2026-10-25

The public E. coli example is pinned at `/examples/ecoli` (job
`c18870c073ad4e8589b92cade5ca45dc`). Its annotated RefSeq chromosome yielded
4,300 proteins, 2,818 integrated EC annotations, four AMR-classified hits,
six stress hits, one virulence hit, three BGC regions and 3,300 localization
assignments. Six generic product names were replaced by qualifying reviewed
homolog names. There are 902 candidate operons, 577 functionally labelled and
89 with a dominant pathway. Quality assessment found 31/32 essential functions
and no remaining consistency violations at the assessment threshold.

The exported model has 3,366 reactions and 2,991 metabolites. Independent SBML
parsing, flux-bound and S*v checks passed (maximum residual <1e-9). Its 22 nutrient
tests retain biomass feasibility without threonine uptake, but not when all 20
listed amino-acid uptakes are closed together. This run took about 19 minutes.
The standalone report, genome viewer, model download and mobile layout were tested.
A separate three-protein reconstruction exercised six repair reactions and
verified that small nonzero fluxes survive export without breaking mass balance.

Full-profile release image: `sha256:6afe114d3c091e9cb39990cb2123e1707d0fc2f3f5ef58a7c7323ff24834b494`.
The full-profile run used `sha256:176dd900735bf7e21998afc718c83756c3fc2836ec91d0feba724cbf99e9de6a`;
the subsequent change only affects explicit homology/FBA opt-outs, tested through
the public API. PSORTb image: `sha256:f39ebfecc77940355da37c6e17e9f80b125b482eec5f5c2ff159b27095a201d1`.

Compose v1 can fail to recreate an existing container with `KeyError: ContainerConfig`.
After backing up the configuration and confirming no active jobs, remove only the
stopped GSPA application container and rerun `up -d`. Keep the persistent bind
mounts and the previous image for rollback.


## EC report navigation update

Current runtime image: `sha256:214534a9c46e0a378d88dc1203e6edaa6fea4019d5e97d037f3aee76130afbd7`.
The Enzymes (EC) tab and header shortcut expose the protein–EC predictions,
with search, score filtering, evidence sources and a filtered TSV download.
EC annotations appear first in protein details. Direct example link:
`https://gspa.bio2vec.net/examples/ecoli#enzymes`.
The existing full-profile HTML reports and ZIPs were refreshed without changing
predictions. The E. coli view and export were verified with 2,818 EC predictions
across 1,812 proteins; the browser regression suite passed 14 tests.
