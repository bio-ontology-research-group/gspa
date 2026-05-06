# `gspa-nf/` — Nextflow wrapper around the annotation pipeline

A container-based alternative to the JVM CLI. Runs the same external tools
(DIAMOND, HMMER, InterProScan, FoldSeek, eggNOG-mapper, AMRFinderPlus,
antiSMASH, BARRNAP, MINCED, CheckM2, GTDB-Tk, …) as Nextflow processes,
each in its own Singularity / Docker image. Useful when you want a
resumable, HPC-native deployment without installing tools on the worker
nodes.

**Not part of the Gradle build** — it's a sibling pipeline, not a
submodule. `./gradlew build` does not touch this directory.

## When to use Nextflow vs. the JVM CLI

| Scenario | Use |
|---|---|
| Single genome, tools already installed locally | `./gradlew :gspa-cli:run --args='annotate ...'` |
| HPC with Singularity, many genomes, checkpoint/resume | `nextflow run gspa-nf/main.nf` |
| Integration testing the GSPA Groovy code | JVM CLI |
| First-time pipeline on a fresh cluster | Nextflow (containers handle deps) |

The Nextflow path does **not** run the GSPA integration / quality
machinery — it produces the raw per-tool outputs plus a merged TSV. Feed
that into `gspa-cli integrate` or `gspa-cli evaluate` if you want GSPA's
combiner + quality scores.

## Quick start

```bash
# Minimal: DIAMOND + Pfam on a single FASTA using Docker
nextflow run gspa-nf/main.nf \
    --input genome.fna \
    --diamond_db /path/to/uniprot_sprot.dmnd \
    --pfam_db /path/to/Pfam-A.hmm \
    -profile docker
```

```bash
# Batch over a samplesheet
cat > sheet.csv <<EOF
sample_id,fasta
mg1655,/data/mg1655.fna
bsubtilis,/data/bsubtilis.fna
EOF
nextflow run gspa-nf/main.nf --input sheet.csv -c databases.config -profile docker
```

## On unimatrix01 (SLURM + Singularity)

See `UNIMATRIX01.md` for the fully worked example. Short version:

```bash
# One-time: pull containers into the shared Singularity cache
bash gspa-nf/stage_images.sh

# Per-run
nextflow run gspa-nf/main.nf \
    -c gspa-nf/slurm_singularity.config \
    -c gspa-nf/databases.config \
    --input samplesheet.csv \
    --outdir /data/<user>/gspa-nf-out \
    -profile singularity_slurm
```

The SLURM config targets the `debug` partition with GlusterFS-aware
`beforeScript` handling (see `UNIMATRIX01.md` for the GlusterFS quirks
that motivated it).

## Pipeline structure

`main.nf` orchestrates 22 processes across 8 module files:

| Module | Processes | When they run |
|---|---|---|
| `modules/gene_calling.nf` | `PYRODIGAL` | Always — produces proteins + GFF. |
| `modules/similarity.nf` | `DIAMOND_BLASTP`, `MMSEQS2_SEARCH`, `FOLDSEEK` | `--run_diamond` (default on), `--run_mmseqs2`, `--run_foldseek`. |
| `modules/domains.nf` | `HMMSEARCH`, `INTERPROSCAN`, `EGGNOG_MAPPER`, `DBCAN` | `--run_hmmer` (default), `--run_interproscan`, `--run_eggnog`, `--run_dbcan`. |
| `modules/specialized.nf` | `BARRNAP`, `MINCED`, `AMRFINDER`, `ANTISMASH`, `SIGNALP`, `CHECKM2`, `GTDBTK` | Flags and DB paths per tool. |
| `modules/quality.nf` | `MERGE_ANNOTATIONS` | Always — emits the merged per-sample TSV. |
| `modules/integrate.nf` | `BUILD_CLAIMS`, `INTEGRATE` | `--run_integrate` — runs Phase 7 evidence integrator end-to-end (claims.jsonl + integrated posteriors with priors). |
| `modules/neural.nf` | `ESM2_DEEPGOPLUS`, `ESM2_CENTROID`, `CLEAN`, `PROTEINFER` | `--run_esm2_deepgoplus`, `--run_esm2_centroid`, `--run_clean`, `--run_proteinfer`. |
| `modules/ensemble.nf` | `ENSEMBLE_PREDS` | `--run_ensemble` — fuses all enabled neural outputs (`max\|mean\|rank`). |
| `modules/eval.nf` | `EVAL_PGAP` | `--run_eval` + `--truth_dir` — F-max, Smin, EC F-max per predictor. |
| `modules/report.nf` | `MAKE_REPORT` | `--run_report` — emits `<sample>.{html,ttl,jsonld}`. SIO-based RDF vocabulary; configurable per-predictor inputs. |

Every process is toggled by a `--run_<tool>` flag and a DB path; omit
either and the step is skipped (empty-file stub is joined into the merge
so the channel topology stays stable).

## Configuration

- `nextflow.config` — default parameters, container images, profile definitions.
- `databases.config` — DB paths; edit locally or pass `-c your_dbs.config`.
- `slurm_singularity.config` — unimatrix01 SLURM profile with GlusterFS
  workarounds.

`-profile` options: `docker`, `singularity`, `singularity_slurm`, `conda`
(declare these in `nextflow.config` — check there for the current list).

## Helper scripts

- `stage_images.sh` — pre-pulls every referenced Singularity image into
  the shared cache. Handles the GlusterFS "0-byte linkto" case documented
  in `UNIMATRIX01.md` by pulling into `/tmp` first and `cp`-ing across.
- `prepull_singularity.sh` — simpler variant for environments without the
  GlusterFS quirks.
- `run_unimatrix01.sh` — convenience wrapper that stitches the above
  together for a single-command unimatrix run.

## Neural function predictors

Four neural predictors plus an ensemble + eval pass — validated on a
21-genome bacterial panel (`benchmark/RESULTS.md`). All opt-in; nothing
runs by default.

### Predictors

| Predictor | Annotation | Backbone | Container | GPU? | F-max (CAFA-prop, n=17) |
|---|---|---|---|---|---|
| ESM2-DeepGOPlus | GO | frozen ESM2-t33 + FC | `gspa-esm2-stack` | Yes | 0.39 |
| ESM2-centroid | GO + EC | NPZ centroids over SwissProt | `gspa-esm2-stack` | CPU OK | 0.04 |
| CLEAN | EC | ESM2 + contrastive head | `gspa-esm2-stack` | Yes | (EC: 0.85) |
| ProteInfer | GO + EC | shallow CNN over sequence | `gspa-proteinfer-stack` | CPU OK | 0.52 |
| **ensemble-mean** | GO + EC | per-(protein,term) mean over enabled | `gspa-eval-stack` | – | **0.55 / EC 0.88** |

See `benchmark/RESULTS.md` for full per-truth-source numbers including
Smin and the propagated-vs-unpropagated comparison that demonstrates
the CAFA-style propagation effect.

### Database fetch

A versioned manifest (`gspa-nf/database_manifest.tsv`) lists every
artefact, version, URL, sha256, and size. Most artefacts are mirrored on
`https://gspa.bio2vec.net/db/`; **CLEAN is fetched directly from upstream
(`tttianhao/CLEAN`)** because its license does not allow re-hosting. The
`fetch_databases.sh` helper downloads only what you ask for:

```bash
# Pull just ProteInfer + the smoke truth fixture, no GPU stack
bash gspa-nf/fetch_databases.sh \
    --predictor proteinfer --predictor truth \
    --dest ~/gspa-db
```

The script prints the `databases.config` lines to paste in. To use a
different host (e.g. local mirror) override the URL prefix:

```bash
GSPA_DB_BASE=https://my-mirror.example.org/gspa-db \
  bash gspa-nf/fetch_databases.sh --predictor esm2-centroid --dest ~/gspa-db
```

### Run

```bash
# Single neural predictor, CPU
nextflow run gspa-nf/main.nf -profile docker -c databases.config \
    --input genome.fna \
    --run_proteinfer --proteinfer_model_dir ~/gspa-db/proteinfer/v1.0/model

# All four + ensemble + eval + report (HTML/TTL/JSON-LD), GPU
nextflow run gspa-nf/main.nf -profile docker,gpu -c databases.config \
    --input samplesheet.csv \
    --run_esm2_deepgoplus --run_esm2_centroid --run_clean --run_proteinfer \
    --run_ensemble --ensemble_mode mean \
    --run_eval --truth_dir ~/gspa-db/truth --go_obo ~/gspa-db/go_aspect_map/v2026.04/go.obo \
    --run_report
```

### Output formats

Per sample, under `${outdir}/${sample_id}/`:

| Path | Format | Source |
|---|---|---|
| `<predictor>/<sample>.<predictor>.tsv` | TSV (4-col) | each enabled neural process |
| `ensemble/<sample>.ensemble.tsv` | TSV (4-col) | `ENSEMBLE_PREDS` |
| `eval/<sample>.<predictor>.eval.json` | JSON | `EVAL_PGAP` (one per predictor) |
| `report/<sample>.html` | HTML | `MAKE_REPORT` (`--run_report`) |
| `report/<sample>.ttl` | RDF / Turtle | `MAKE_REPORT` |
| `report/<sample>.jsonld` | JSON-LD | `MAKE_REPORT` |
| `<sample>_annotations.tsv` | TSV (6-col, classical) | `MERGE_ANNOTATIONS` |
| `<sample>_annotated.gff3` | GFF3 (classical only) | `MERGE_ANNOTATIONS` |

The RDF model uses a small custom vocabulary at
`https://gspa.bio2vec.net/ns/` layered on
[SIO](http://semanticscience.org/) — every prediction is a
`gspa:FunctionPrediction` (subClassOf `sio:000663` data item) with
`gspa:hasTarget` (UniProt or per-sample protein IRI), `gspa:hasFunction`
(OBO GO IRI or UniProt-enzyme IRI for EC), `gspa:hasScore` (`xsd:float`),
`gspa:hasPredictor`, and `gspa:annotationType`. Every predictor (existing
or future) appears as one `gspa:Predictor` instance. Eval metrics enter as
`gspa:EvalRecord` instances. The TTL and JSON-LD agree triple-for-triple
(verified with `rdflib`).

### Smoke test (no GPU, no big DBs)

```bash
bash gspa-nf/fetch_databases.sh --predictor proteinfer --predictor truth --dest ~/gspa-db
GSPA_DB=~/gspa-db nextflow run gspa-nf/main.nf -profile docker,smoke_proteinfer
ls test/nf-output/smoke/proteinfer/smoke.proteinfer.tsv
ls test/nf-output/smoke/eval/smoke.proteinfer.eval.json
```

### Container images

Three Dockerfiles under `gspa-nf/docker/`:

- `Dockerfile.esm2-stack` (PyTorch + fair-esm) — covers DGP, centroid, CLEAN
- `Dockerfile.proteinfer-stack` (TF 2.15) — covers ProteInfer
- `Dockerfile.eval-stack` (python:3.12-slim) — covers ensemble + eval

Build + push with the recipes in `gspa-nf/docker/README.md`. Images are
published at `docker.io/leechuck/gspa-{esm2,proteinfer,eval}-stack:0.1`
and wired in `nextflow.config` (process `withName:` directives).

### Performance budget (per genome, ~4k proteins)

| Predictor | RTX-4090 | CPU (16 cores) |
|---|---|---|
| ProteInfer | n/a | ~5 min |
| ESM2-centroid | ~3 min | ~10 min |
| ESM2-DeepGOPlus (t33) | ~5 min | ~40 min |
| CLEAN | ~10 min | impractical |
| FoldSeek-homology + ProstT5 | ~3 hr | ~12 hr |
| Ensemble fusion | ~10 min | ~10 min |

## Region predictors (FOSS)

Per-residue predictors that emit `(protein_id, region_start, region_end,
region_type, score)` 5-column TSVs. Each region appears in the
`region/` HTML section and as a `gspa:Region` instance in the RDF/JSON-LD.

| Predictor | Capability | License | Container |
|---|---|---|---|
| `--run_metapredict` | Disorder regions (Metapredict v2) | MIT | `gspa-region-stack` |
| `--run_deepsig` | Sec/Tat signal peptides | GPL-3.0 | `gspa-region-stack` |
| `--run_tmbed` | TM helices via ProtT5 | Apache-2.0 | `gspa-region-stack` |
| `--run_tppred3` | N-terminal targeting peptides | GPL-3.0 | `gspa-region-stack` |

## Localization & secretion

| Predictor | Capability | License | Container |
|---|---|---|---|
| `--run_psortb` | Bacterial subcellular localization (5–8 classes) | GPL-3.0 | upstream `brinkmanlab/psortb_commandline:1.0.4` |

PSORTb output is whole-protein (4-col TSV), so it auto-joins
`ENSEMBLE_PREDS` when `--run_ensemble` is on.

## Term-extra predictors (auto-join the ensemble)

| Predictor | Capability | License | Container |
|---|---|---|---|
| `--run_deepfri` | Sequence-only GO (MF/BP/CC) | BSD-3-Clause | `gspa-tf-stack` |
| `--run_deepec` | EC numbers | **AGPL-3.0** ⚠ | `gspa-tf-stack` |
| `--run_deeparg` | Antimicrobial-resistance gene calls | MIT | `gspa-tf-stack` |

⚠ **DeepEC AGPL clause:** fine to bundle and run locally / on HPC. If
you ever expose GSPA as a hosted web service, the AGPL network-clause
requires publishing source to users of that service. If that's a
problem, omit `--run_deepec` — CLEAN + DIAMOND already cover EC well.

## Site-level predictors

| Predictor | Capability | License | Container |
|---|---|---|---|
| `--run_musitedeep` | PTM sites (phospho-S/T/Y by default) | MIT | `gspa-tf-stack` |
| `--run_scannet` | PPI interface residues (needs structures) | Apache-2.0 | `gspa-struct-stack` |

ScanNet requires per-protein PDB/CIF files under
`${outdir}/${sample_id}/structures/`. Set `--structures_from
{esmfold,afdb}` to provision them, or pre-stage your own.

## License-restricted tools (replaced)

The following tools were previously considered but are NOT included
because their licenses are not OSI-approved Free Software. The chosen
FOSS replacement is in the right column.

| Excluded tool | Restriction | FOSS replacement |
|---|---|---|
| DeepLoc 2 | DTU academic EULA | PSORTb 3.0 |
| MULocDeep | "academic users only" | (no FOSS DL alt; PSORTb covers bacterial) |
| DeepTMHMM | BioLib academic ToS | TMbed |
| IUPred3 | enquiry-walled | Metapredict |
| TargetP 2 | CC BY-NC-SA | TPpred 3 |
| SignalP 6 | DTU academic EULA | DeepSig |
| NetPhos / NetPhosBac | DTU EULA | MusiteDeep |

Existing JVM wrappers for SignalP 6 and DeepTMHMM remain in the
codebase from v1.1.0 for users who have those licenses, but are NOT
wired into Nextflow. Use the FOSS replacements above by default.

## End-to-end with Phase 7 integration (`--run_integrate`)

`MERGE_ANNOTATIONS` writes a merged TSV that is **not** the same as the
GSPA integrated TSV — it's a thin per-tool union without the Bayesian
combiner. To run the full Phase 7 integrator inside Nextflow, set
`--run_integrate` and supply the reference data:

```bash
./gradlew :gspa-cli:shadowJar    # one-time, produces gspa-1.5.0.jar

nextflow run gspa-nf/main.nf \
    -profile docker \
    --input genome.fna \
    --diamond_db /refs/uniprot_sprot.dmnd \
    --pfam_db    /refs/Pfam-A.hmm \
    --run_integrate \
    --gspa_jar   $PWD/gspa-cli/build/libs/gspa-1.5.0.jar \
    --goa        /refs/goa_uniprot_all.gaf.gz \
    --go_owl     /refs/go.owl \
    --ec2go      /refs/ec2go \
    --pathways   /refs/kegg_pathways.tsv
```

This adds two processes after `MERGE_ANNOTATIONS`:
- `BUILD_CLAIMS` — wraps `benchmark/02b_parse_predictors_to_claims.py`,
  emits `${sample_id}_claims.jsonl`.
- `INTEGRATE` — invokes `gspa-cli integrate` on that claims.jsonl with
  the full prior stack, emits `${sample_id}_integrated.tsv`
  (per-(protein, function) posterior probabilities with provenance).

Optional flags: `--pfam2go`, `--theta_file`, `--essential_profile`
(default `bacteria`), `--enable_priors` (comma-separated; default
`essentiality,coherence,gap_filling,genomic_context`).

`--run_integrate` is opt-in. Without it, the pipeline behaves exactly
as before — produces raw per-tool outputs + the merged TSV, leaving
integration to a manual `gspa-cli integrate` invocation.

## Relationship to the JVM side

If you change the JVM predictor wrappers in a way that changes expected
output columns, update the matching Nextflow module under
`modules/` — the two currently duplicate output-format knowledge.
