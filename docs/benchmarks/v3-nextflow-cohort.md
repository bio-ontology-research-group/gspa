# Joint cohort inference through Nextflow

The Nextflow route now collects all genomes before invoking the same two-pass
batch integration used by `gspa integrate --batch-manifest`. Pass one constructs
the joint raw-evidence snapshot and held-out cohort context; pass two scores each
genome with that context. Nextflow implements no parallel scoring algorithm.
This advances recommendation 2. It does not certify all prediction adapters,
the full upstream annotation pipeline or biological accuracy.

## Running existing claims

Use the official Nextflow launcher with Java 21. Do not invoke its framework JAR
directly: the launcher supplies JVM module-access options needed by its cache.
The existing repository-wide Nextflow config uses Groovy syntax; with Nextflow 26
set `NXF_SYNTAX_PARSER=v1`. The new integration entry point and module are tested
with both parsers when using the isolated test config.

```bash
export TMPDIR="$HOME/tmp"
export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$TMPDIR"
export NXF_SYNTAX_PARSER=v1
nextflow run gspa-nf/cohort.nf \
  --cohort_manifest batch.tsv \
  --cohort_resources resources.json \
  --gspa_jar gspa-cli/build/libs/gspa-3.0.0.jar \
  --outdir results -with-trace trace.tsv
```

The TSV accepts the production CLI's batch columns. It requires `genome_id`,
`claims_path` and `genome_layout_path`; optional columns include taxon, essential
profile, operons, gaps, per-genome orthogroups, cluster consensus and FBA medium.
Input paths are resolved relative to this file. Nextflow owns `out_path` and
ignores any supplied output destination. Duplicate/unsafe genome IDs, unknown
columns, missing required inputs and control characters are rejected.

Example `resources.json` (paths resolve relative to this JSON):

```json
{
  "reactionGraphFile": "reference/reactions.tsv",
  "reactionEcAliasesFile": "reference/ec-aliases.tsv",
  "cohortMetadataFile": "context/metadata.json",
  "cohortOrthogroupsFile": "context/orthogroups.tsv",
  "goOwlFile": "reference/go.owl",
  "ec2goFile": "reference/ec2go.txt",
  "pathwaysFile": "reference/pathways.tsv",
  "moduleFiles": ["reference/modules.tsv"],
  "learnedModelFile": "models/model.json",
  "overlayFile": "behavior.json",
  "reasonerCacheDir": "reference/reasoner-cache"
}
```

The first four resources are required. Optional typed resources include diffusion
metabolites, role evidence, interaction edges, reaction directions, theta,
explicit external context, taxonomy, taxon constraints, complex terms and
essential functions. The production validators retain their authority over file
schemas, model contracts, evidence identity and allowed behavior overlays. No
model is enabled unless supplied. The base settings match the existing Nextflow
integration route (`lite=true` and its four-prior default); an overlay can change
behavior, including `lite=false` for full ontology inference. For comparison with
another route, match these settings explicitly as well as the input files.

## Main annotation pipeline

`main.nf --run_integrate --cohort` uses the same process after merging predictor
claims and computing operons. It requires `--cohort_layouts`, `--cohort_metadata`,
`--cohort_orthogroups`, `--reaction_graph`, and `--reaction_ec_aliases`. Layouts are
named `<sample_id>.layout.tsv`; their protein identifiers must match this run's
predicted proteins. Orthogroups are supplied evidence, never guessed from names.
Optional `--cohort_gaps` files are named `<sample_id>.gaps.jsonl`. Cohort inference
without `--run_integrate` fails explicitly. The main entry point has been parsed
in preview mode; this is not an end-to-end rerun of all upstream predictors.

## Staging, caching and output

Every sample file and shared resource is a declared `path` input, staged into an
indexed directory. Identical basenames across genomes and resource roles do not
collide. The task also hashes explicit genome-to-file and role-to-file bindings:
Nextflow hashes file collections without order, so a collection alone would
incorrectly reuse a task after swapping two genomes' claim files. A regression
reproduced that failure before the binding fix and passes afterwards. Deep
content hashing also invalidates same-size, same-timestamp edits.

The optional reasoner cache is copied into the task before use; engine writes
cannot alter its source. Scratch files stay under the task's own `.tmp` directory.
Behavior overlays are staged and shell-quoted, including spaces and apostrophes
in filenames. The staging rules follow the official
[Nextflow process input contract](https://nextflow.io/docs/latest/process.html#multiple-input-files).

One successful task publishes `<outdir>/cohort/`, containing per-genome scores,
provenance, suggestions, cohort context, learned reports and stage ledgers. Feature
exports, when requested, are in `cohort/features/<genome_id>/`. The actual staged
workflow JSON and batch TSV are also preserved. On failure Nextflow fails the
task and retains its work directory; any per-genome failed stage ledgers remain
there. Failed output is not published as a successful result.

## Verification

`benchmark/v3/tests/test_nextflow_cohort.py` executes real Nextflow and the
production JAR on four synthetic genomes with colliding protein IDs and file
basenames. It compares complete prediction TSVs and cohort reports with a direct
CLI batch, checks stage-state parity and feature/provenance exports, and verifies
that the two same-clade genomes are held out together. Other tests exercise cache
reuse, content and binding invalidation, duplicate/missing inputs and failed
cohort ledgers for every genome.

```bash
export GSPA_NEXTFLOW_BIN=/path/to/the/official/nextflow
export GSPA_JAR="$PWD/gspa-cli/build/libs/gspa-3.0.0.jar"
"$HOME/tmp/gspa-v3-testenv/bin/python" -m unittest discover \
  -s benchmark/v3/tests -p test_nextflow_cohort.py -v
```

## Real cohort result (2026-10-26)

The final Nextflow module ran all 18 acquired development genomes in one task,
using the same 484,202 claims, layouts, cohort metadata, reaction references,
ontology and learned logistic model as the previously verified intact batch CLI
run. The full ontology setting and six-prior configuration were matched explicitly.
The staged reasoner cache loaded 19,577 source-verified dependency pairs; its
source remained unchanged. Nextflow reports 3 minutes 16 seconds task duration
and 10.2 GB peak RSS for this cached-reference execution.

`check_nextflow_cohort.py` verified the engine and input hashes, each genome's
staged file bindings, successful execution, and all per-genome outputs:

- All 484,202 prediction rows are byte-identical across the two routes.
- All 18 cohort reports agree in full, including held-out genomes and peers.
- All 484,202 learned candidate records and the learned-model hash agree.
- Stage names/states agree, with cohort inference executed for every genome.

Eight live Nextflow regressions pass. The full benchmark suite ran 482 tests,
with three optional skips and no failures; 15 service configuration tests also
pass. The final main annotation entry point parses successfully in preview mode.
The first real execution rebuilt ontology dependencies; the final run used the
explicit staged cache and final input-binding/empty-input fixes. Only the final
run is certified by the archived parity record.

Evidence is archived in
[`benchmark/v3/results/nextflow-cohort-v1/`](../../benchmark/v3/results/nextflow-cohort-v1/),
including source/input hashes, the direct batch receipt, Nextflow trace, complete
comparison record, test logs and the reproduced pre-fix failures. Bulky production
files remain at the cache paths named in those records. Recheck with:

```bash
"$HOME/tmp/gspa-v3-testenv/bin/python" benchmark/v3/check_nextflow_cohort.py \
  --run "$HOME/.cache/gspa-v3-data/nextflow-cohort-v1" \
  --direct "$HOME/.cache/gspa-v3-data/sequence-loss-v1/runs/intact/out" \
  --out "$HOME/tmp/nextflow-cohort-recheck.json"
```

These are execution and parity checks, not experimentally validated annotations.
Full upstream predictor execution, all remaining route combinations and browser
reconciliation still require their own evidence. The broader twelve-item
acceptance record remains authoritative.
