# gspa annotate

`gspa annotate` runs the per-protein predictors on a genome, writes the
annotations, detects operons and, when a GO OWL file is given, computes
genome-scale quality metrics and optionally enforces consistency, coherence
and completeness on the predictions.

```text
gspa annotate (-i FILE | --genome FILE | --gff3 FILE | --proteins FILE) [options]
```

Settings are merged in this order, later sources overriding earlier ones:
built-in defaults, the kingdom preset selected with `-k`, the workflow preset
selected with `--preset`, the YAML file given
with `-c`, and finally the command-line options on this page. See
[Configuration files](configuration.md).

## Workflow defaults

Use `--preset prokaryote|bacteria|archaea|eukaryote|mag|metagenome`.
`metagenomics` is an alias for `metagenome`. Run `gspa presets` for guidance or
add `--print-config` to print resolved YAML and exit without running tools.
See [workflow presets](../scenarios/presets.md) for dependencies and assumptions.
All existing options remain available. Enforcement flags also accept `--no-`
forms, such as `--no-enforce-consistency`, to disable preset defaults.

## Examples

Bacterial genome assembly, genes called by Pyrodigal or Prodigal:

```bash
gspa annotate -i genome.fna -k bacteria -o out/
```

Reuse the CDS of an existing annotation and compute quality metrics with
taxon-consistency enforcement:

```bash
gspa annotate --genome genome.fna --gff3 genome.gff3 -k bacteria \
  --go-owl go.owl --enforce-consistency --taxon bacteria \
  --metrics-scope both -o out/
```

Metagenome-assembled genome with DeepGO-PlusPlus-Light as the base predictor
and an inferred taxon:

```bash
gspa annotate -i bin_017.fna --mag -k auto \
  --base-predictor light \
  --neural-sidecar "$PWD/benchmark/neural/run_neural_predictors.py" \
  --deepgo-plusplus-light-assets /data/dgpp_light_assets \
  --min-score 0.1 -t 16 \
  --go-owl go.owl --infer-taxon -o out_bin017/
```

Protein FASTA only (single synthetic contig, no operons):

```bash
gspa annotate --proteins proteome.faa -k eukaryote -o out/
```

## Inputs

| Option | Type | Default | Description |
|---|---|---|---|
| `-i`, `--input` | file | none | Primary input: a nucleotide FASTA (genes are called), or a GFF3 file (`.gff`, `.gff3`, optionally gzipped) whose CDS are translated; the GFF3 may carry sequences in an embedded `##FASTA` block. For FASTA input it is equivalent to `--genome`. Takes precedence over `--genome`. |
| `--genome` | file | none | Genome or metagenome nucleotide FASTA, one contig per record. With `--gff3`, the GFF3 CDS are translated against these sequences instead of calling genes. |
| `--gff3`, `--gff` | file | none | GFF3 annotation. With `--genome`, CDS coordinates are reused and translated; with `--proteins`, supplied proteins are joined to contigs through the CDS attributes `protein_id`, `ID` or `locus_tag`; alone, it is the primary input. |
| `--proteins` | file | none | Pre-called protein FASTA. Alone, proteins are annotated on one synthetic contig named `proteins`; operons and per-contig metrics are then not meaningful. With `--gff3`, proteins without a matching CDS are placed on a contig named `unplaced`. |
| `--mag` | flag | off | Treat the input as a MAG: gene calling runs in metagenomic mode and completeness is evaluated with MAG adjustments. Sets `input.type` to `mag`. Does not load the `mag.yaml` preset; pass `-c` for that. |
| `-k`, `--kingdom` | string | none (equivalent to `auto`) | Organism kingdom: `bacteria`, `archaea`, `eukaryote`, `virus` or `auto`. For `bacteria`, `archaea` and `eukaryote` the matching bundled preset is loaded. `virus` and `auto` have no legacy preset (built-in defaults are used). Without `--preset`, the kingdom also selects the essential-function profile for completeness (`auto` uses the bacterial profile) and the SignalP organism group. There is no automatic kingdom detection in this release. |

At least one of `-i`, `--genome`, `--gff3` or `--proteins` is required. The
genome identifier used in output file names is the input file name without
its FASTA or GFF3 extension.

## Configuration and output

| Option | Type | Default | Description |
|---|---|---|---|
| `-c`, `--config` | file | none | YAML configuration file ([schema](configuration.md)). A missing file is an error. |
| `-o`, `--output` | directory | `gspa_output` | Output directory; created if needed. Contents are listed under [annotate output directory](file-formats.md#annotate-output-directory). |
| `--no-provenance` | flag | off | Do not record enforcement provenance: omits the `provenance` column of the annotations TSV and `<genome>_enforcement_actions.tsv`. |
| `--min-score` | double in [0, 1] | none (keep everything) | Drop GO annotations scoring below this value right after prediction, before metrics and output. Also raises the emission threshold of the DeepGO-PlusPlus-Light sidecar to this value so the low-score tail is never computed. `0` keeps everything. |

## Predictors

| Option | Type | Default | Description |
|---|---|---|---|
| `--db`, `--database` | path | none | Sequence-similarity database (DIAMOND or MMseqs2, depending on `predictors.similarity.tool`). |
| `-t`, `--threads` | int | 0 | Thread budget. A value above 0 sets the DIAMOND thread count of DeepGO-PlusPlus-Light (`predictors.neural.deepGoPlusPlusLight.threads`, default 8). It does not affect other predictors. |
| `--disorder` | flag | off | Enable the Metapredict intrinsic-disorder predictor. |
| `--foldseek-centroid-db` | path | none | FoldSeek function-centroid database built with `benchmark/neural/build_foldseek_centroids.py`. Enables the structure predictor. |
| `--foldseek-centroid-mode` | string | `none` (from configuration) | FoldSeek centroid mode: `none`, `go`, `ec` or `both`. Any value other than `none` makes FoldSeek search the centroid database instead of doing homology transfer against `predictors.structure.database`. Giving this option enables the structure predictor. |

### Neural predictors

All neural predictors run through one Python sidecar,
`benchmark/neural/run_neural_predictors.py`, executed with
`predictors.neural.pythonExecutable` (default `python3`).

| Option | Type | Default | Description |
|---|---|---|---|
| `--neural-sidecar` | path | none | Absolute path to `run_neural_predictors.py`. Required by every neural predictor below. |
| `--esm2-deepgoplus` | flag | off | Enable the ESM2-based DeepGO-Plus GO predictor (base ESM2 head only). |
| `--esm2-deepgoplus-checkpoint` | path | none | Trained ESM2 head checkpoint (`.pt`). Required with `--esm2-deepgoplus`. |
| `--esm2-deepgoplus-terms` | path | none | GO vocabulary, one term per line, in the order of the checkpoint's output units. |
| `--esm2-deepgoplus-model` | string | `esm2_t33_650M_UR50D` | ESM2 variant that matches the checkpoint (for example `esm2_t12_35M_UR50D`). |
| `--proteinfer` | flag | off | Enable the ProteInfer GO predictor. |
| `--proteinfer-model-dir` | path | none | ProteInfer model directory. |
| `--clean` | flag | off | Enable the CLEAN EC-number predictor. |
| `--clean-model-dir` | path | none | CLEAN checkpoint directory (contains `CLEAN_infer_fasta.py`). |

### DeepGO-PlusPlus base predictor

| Option | Type | Default | Description |
|---|---|---|---|
| `--base-predictor` | string | `auto` (from configuration) | Select DeepGO-PlusPlus as the base GO predictor. `full` (aliases `deepgo-plusplus`, `deepgoplusplus`, `gpu`): learned stacker over precomputed component scores. `light` (aliases `deepgo-plusplus-light`, `deepgopluspluslight`, `cpu`): self-contained CPU predictor that runs DIAMOND and STRING network k-NN from the FASTA. `none`: use the individual enable flags. `auto`: use Light when its assets are configured, otherwise full when its integrator and components are configured, otherwise neither. `full` and `light` are mutually exclusive: the chosen one is enabled and the other disabled. Other values are rejected. |
| `--deepgo-plusplus`, `--cafa-baseline` | flag | off | Enable DeepGO-PlusPlus (full). `--cafa-baseline` is a legacy alias. |
| `--deepgo-plusplus-integrator`, `--cafa-baseline-integrator` | path | none | Frozen integrator JSON written by `deepgo-plusplus/pipeline/train_integrator.py --save-model`. |
| `--deepgo-plusplus-components-dir`, `--cafa-baseline-components-dir` | path | none | Directory of per-component score files (`<component>.tsv` or `<component>.tsv.gz`). |
| `--deepgo-plusplus-dag`, `--cafa-baseline-dag` | path | none | GO DAG file (`child<TAB>ancestor`) used for true-path propagation. |
| `--deepgo-plusplus-light` | flag | off | Enable DeepGO-PlusPlus-Light. |
| `--deepgo-plusplus-light-assets` | path | none | Asset bundle from `make_assets.sh` (`train_db.dmnd`, `train_net_index.tsv`, `train_terms.tsv`, `go-dag.tsv`, optionally `go.obo` and `cnn_model.pt`). Required for Light; with `--base-predictor auto` it is enough to set this. |
| `--deepgo-plusplus-light-models` | path | `deepgo-plusplus/models` (resolved by the sidecar) | Directory of frozen integrator JSONs. |
| `--deepgo-plusplus-light-interpro` | flag | off | Add the InterProScan domain component (needs `--deepgo-plusplus-light-interproscan`). |
| `--deepgo-plusplus-light-cnn` | flag | off | Add the CPU 1D-CNN component (improves coverage of proteins without homologs; needs PyTorch and `cnn_model.pt`). |
| `--deepgo-plusplus-light-interproscan` | path | none | InterProScan launcher; enables the InterPro component. |

See [DeepGO-PlusPlus base predictor](../scenarios/deepgo-plusplus.md).

## Genome-scale metrics

Metrics and enforcement run only when `--go-owl` is given and at least one GO
annotation was predicted. Concepts are described in
[Genome-scale quality metrics](../concepts/quality-metrics.md).

| Option | Type | Default | Description |
|---|---|---|---|
| `--go-owl` | path | none | GO OWL file. Enables completeness, coherence, consistency and information-content metrics and is required by all enforcement options. Also used to add GO labels and aspects to the annotations TSV. |
| `--no-metrics` | flag | off | Skip metrics even when `--go-owl` is given. Enforcement still runs if requested. |
| `--metrics-scope` | string | `contig` | `contig`: one report per contig (each contig evaluated on its own proteins, appropriate for metagenomes and multi-replicon assemblies). `genome`: one pooled report. `both`: both. |
| `--complex-terms` | path | bundled `complex_terms.tsv` | Complex GO term classification (`term<TAB>h|n|a`) for complex coherence and its enforcement ([format](file-formats.md#complex-terms-tsv)). |

## Taxon consistency

| Option | Type | Default | Description |
|---|---|---|---|
| `--taxon` | string | none | Assert the organism taxon (`NCBITaxon_2`, `NCBITaxon:2`, `2`, or `bacteria`, `archaea`, `eukaryote`, `virus` and their variants). Terms that cannot occur in the lineage are then flagged individually. Unrecognised strings are ignored with a warning. For single-organism inputs. |
| `--infer-taxon` | flag | off | Infer the organism taxon from the predictions through the GO taxon constraints and use it for consistency enforcement. Writes `<genome>_taxon_inference.tsv`. Ignored if `--taxon` is given. |
| `--infer-taxon-min-score` | double in [0, 1] | 0.9 (from configuration) | Only predictions scoring at least this value inform taxon inference. The help text prints `Default null`; the effective default is `quality.consistency.inferTaxonMinScore` = 0.9. |
| `--taxon-constraints` | path | bundled `go-taxon-constraints.tsv` | Taxon constraints as OBO (`.obo`, e.g. `go-computed-taxon-constraints.obo`) or TSV ([format](file-formats.md#taxon-constraint-tsv-and-obo)). A missing file falls back to the bundled constraints with a warning. |
| `--taxonomy` | path | bundled `ncbi-taxon-hierarchy.tsv` | Taxon hierarchy and disjointness file (`Term<TAB>is_a|disjoint_from|union_of<TAB>term`) ([format](file-formats.md#taxonomy-file)). |

## Enforcement

Enforcement passes run in the order consistency, coherence, completeness. All
promotions are checked against the taxon constraints. Every action is recorded
in `<genome>_enforcement_actions.tsv` unless `--no-provenance` is set, and the
annotation outputs are rewritten afterwards. See
[Stage 2A enforcement](../concepts/stage2a-enforcement.md).

| Option | Type | Default | Description |
|---|---|---|---|
| `--enforce-consistency` | flag | off | Act on GO annotations that violate taxon or anatomical constraints. |
| `--enforce-consistency-mode` | string | `remove` (from configuration) | `remove`: delete the annotation. `downrank`: multiply its score by `quality.consistency.downrankFactor` (0.5). `flag`: keep it and record the violation. `minimal-flip`: joint minimum-cost demotion with weighted MaxSAT, resolving co-annotation conflicts. `relax` is also accepted: replace the term with its nearest constraint-satisfying non-root ancestor. Other values behave like `remove`. |
| `--enforce-completeness` | flag | off | For each missing essential function, promote it on the protein with the strongest sub-threshold evidence (minimum score 0.05), recorded with evidence code `ISC`. |
| `--enforce-coherence` | flag | off | Fix singleton subunits of obligate heteromeric complexes (demote the lone protein, or promote a plausible partner) and promote missing `has_part` partner terms. Initialises the ELK reasoner, which is slower and needs more memory. |

## Standard options

| Option | Description |
|---|---|
| `-h`, `--help` | Show help and exit. |
| `-V`, `--version` | Print the version and exit. |

## Notes

- Predictors whose external tool or database is not available are skipped; the
  log lists how many predictors were registered and which ran.
- Each run writes `<genome>_timing.tsv` with wall-clock time per phase
  (`setup`, `prediction`, `filter`, `output`, `ontology_load`,
  `taxon_inference`, `enforcement`, `metrics`).
- Operon files are always written when proteins have coordinates, so that
  `gspa integrate --operons out/operons_for_integrate.tsv` and
  `gspa visualize` can use them directly.
