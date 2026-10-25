# Configuration files

GSPA supports declarative YAML configuration files to configure gene calling,
per-protein functional predictors, genome-scale quality metrics, and symbolic
constraint enforcement. The configuration file is passed to `gspa annotate`
using the `-c` or `--config` option.

## Configuration resolution hierarchy

When `gspa annotate` runs, settings are merged hierarchically in the following
precedence order (later sources override earlier sources):

1. **Built-in class defaults**: Hardcoded default values in the `GspaConfig` model classes.
2. **Kingdom preset**: A legacy profile (`bacteria`, `archaea`, `eukaryote`) loaded via `-k` or `--kingdom`.
3. **Workflow preset**: Optional `--preset` defaults; see [workflow presets](../scenarios/presets.md).
4. **User configuration file**: A YAML file provided via `-c` or `--config`.
5. **Command-line overrides**: Explicit CLI arguments passed to `gspa annotate`.

Only keys actually present in the YAML override inherited settings; a missing config file is an error. `--print-config` prints the merged YAML without starting a run.

Nested maps are merged recursively. Scalar values, lists, and sets from higher-precedence sources replace lower-precedence values.

## Configuration keys and sections

### Input section (`input`)

Configures genome sequence format, taxonomic kingdom, and coordinates.

| Key | Type | Default | Description |
|---|---|---|---|
| `format` | string | `fasta` | Primary input format: `fasta`, `genbank`, or `gff`. |
| `type` | string | `genome` | Input biological type: `genome`, `mag`, or `metagenome`. |
| `kingdom` | string | `auto` | Taxonomic kingdom: `bacteria`, `archaea`, `eukaryote`, `virus`, or `auto`. |
| `proteinFasta` | string | `null` | Path to a pre-called protein FASTA file. Joined to contigs via `gff3` when provided. |
| `gff3` | string | `null` | Path to an existing GFF3 annotation paired with the genome FASTA. CDS features are translated against the genome rather than re-called. |
| `translateAltStart` | boolean | `true` | Emit Methionine for alternative initiator codons (GTG, TTG) when translating CDS features. |
| `requireExistingGenes` | boolean | `false` | Refuse fallback gene calling when no proteins/CDS were loaded; enabled by the eukaryote workflow preset. |
| `proteinsOnly` | boolean | `false` | When true, annotate a bare protein FASTA on a single synthetic contig without chromosomal operon context. |

### Predictor section (`predictors`)

Configures sequence similarity, domain annotation, structural comparison, operon calling, metabolic pathways, and neural models.

#### Similarity (`predictors.similarity`)

| Key | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Enable sequence similarity search. |
| `tool` | string | `diamond` | Alignment engine: `diamond` or `mmseqs2`. |
| `database` | string | `null` | Path to reference alignment database (e.g. UniProt Swiss-Prot `.dmnd`). |
| `evalue` | float | `1e-5` | Maximum alignment expectation value (E-value) threshold. |
| `maxTargetSeqs` | integer | `10` | Maximum number of target hits retained per query protein. |

#### Structure (`predictors.structure`)

| Key | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `false` | Enable structural similarity search. |
| `tool` | string | `foldseek` | Structural search tool (typically `foldseek`). |
| `database` | string | `null` | Path to structure database (e.g. AlphaFold DB / Swiss-Prot). |
| `predictStructures` | boolean | `false` | Predict de novo tertiary structures via ESMFold prior to search. |
| `prostt5Model` | string | `null` | Path to ProstT5 model directory for sequence-to-structure translation. |
| `centroidDb` | string | `null` | Path to FoldSeek function-centroid database. |
| `centroidMode` | string | `none` | FoldSeek centroid mode: `none`, `go`, `ec`, or `both`. |

#### Domains (`predictors.domains`)

| Key | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Enable profile HMM domain annotation. |
| `tool` | string | `interproscan` | Execution tool: `interproscan` or direct `hmmer`. |
| `applications` | list of string | `['Pfam', 'TIGRFAM', 'CDD', 'SUPERFAMILY']` | InterProScan member database analyses to execute. |

#### AntiFam pseudogene filter (`predictors.antifam`)

| Key | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `false` | Enable AntiFam HMM scanning for spurious or pseudogenic ORFs. |
| `database` | string | `null` | Path to concatenated `AntiFam.hmm` database. |
| `useCutGa` | boolean | `true` | Use model-specific gathering thresholds (`--cut_ga`). |
| `evalue` | float | `1e-5` | Fallback E-value threshold when `useCutGa` is false. |

#### Operons (`predictors.operons`)

| Key | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Detect operon structures based on intergenic distance and strand orientation. |
| `maxIntergenicDistance` | integer | `300` | Maximum intergenic distance in base pairs between adjacent co-operonic CDSs. |
| `requireSameStrand` | boolean | `true` | Require adjacent genes to reside on the same strand. |

#### Pathways (`predictors.pathway`)

| Key | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Enable metabolic pathway reconstruction. |
| `tool` | string | `gapseq` | Pathway reconstruction engine (typically `gapseq`). |
| `crossfeeding` | boolean | `false` | Enable crossfeeding community analysis for MAGs and metagenomes. |
| `exportSbml` | boolean | `false` | Export reconstructed networks as SBML models. |
| `runFba` | boolean | `false` | Execute Flux Balance Analysis via COBRApy. |

#### Subcellular localization (`predictors.localization`)

| Key | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `false` | Master switch for subcellular localization predictors. |
| `signalP` | boolean | `true` | Run SignalP signal peptide prediction (requires license). |
| `tmhmm` | boolean | `true` | Run TMHMM transmembrane helix prediction (requires license). |
| `psort` | boolean | `false` | Run legacy PSORT predictor. |
| `deepSig` | boolean | `false` | Run open-source DeepSig signal peptide predictor. |
| `tmbed` | boolean | `false` | Run open-source TMbed ProtT5 transmembrane predictor. |
| `tppred3` | boolean | `false` | Run open-source TPpred3 targeting peptide predictor. |
| `psortb` | boolean | `false` | Run open-source PSORTb bacterial localization predictor. |
| `deepSigKingdom` | string | `gramn` | DeepSig model kingdom: `gramp`, `gramn`, or `euk`. |
| `tppred3Kingdom` | string | `nonplant` | TPpred3 model kingdom: `plant` or `nonplant`. |
| `psortbGram` | string | `negative` | PSORTb target organism type: `positive`, `negative`, or `archaea`. |

#### Intrinsic disorder (`predictors.disorder`)

| Key | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `false` | Enable Metapredict intrinsic disorder prediction. |
| `minRegionLen` | integer | `10` | Minimum amino acid length of an intrinsically disordered region to report. |
| `minScore` | float | `0.5` | Minimum mean disorder probability threshold. |

#### Neural predictors (`predictors.neural`)

Configures deep learning models mediated by `benchmark/neural/run_neural_predictors.py`.

| Key | Type | Default | Description |
|---|---|---|---|
| `sidecarScript` | string | `null` | Path to `run_neural_predictors.py` sidecar script. |
| `pythonExecutable` | string | `python3` | Python binary used to run the neural sidecar. |
| `basePredictor` | string | `auto` | Function predictor backbone: `auto`, `none`, `full`, or `light`. `auto` selects `light` if assets are configured, else `full` if components are configured. |

Sub-configs within `predictors.neural`:

- `deepGoPlusPlusLight`:
  - `enabled` (boolean, default `false`): Run self-contained CPU-based DeepGO-PlusPlus-Light.
  - `assets` (string): Directory with `train_db.dmnd`, `train_net_index.tsv`, `train_terms.tsv`, and `go-dag.tsv`.
  - `modelsDir` (string): Directory containing frozen logistic regression integrator JSON files.
  - `diamond` (string, default `diamond`): Path to DIAMOND executable.
  - `interpro` (boolean, default `false`): Include InterProScan domain component.
  - `cnn` (boolean, default `false`): Include 1D-CNN orphan sequence component.
  - `interproscan` (string): Path to InterProScan launcher.
  - `cnnModel` (string): Path to CNN weights file (`cnn_model.pt`).
  - `threads` (integer, default `8`): DIAMOND worker threads.
  - `topK` (integer, default `5`): Homolog k-nearest neighbors voted per query.
  - `batchSize` (integer, default `16`): Batch size.
  - `minScore` (float, default `0.1`): Minimum score threshold.
- `deepGoPlusPlus`:
  - `enabled` (boolean, default `false`): Run full GPU DeepGO-PlusPlus integrator over precomputed scores.
  - `integrator` (string): Path to frozen integrator model JSON.
  - `componentsDir` (string): Directory of per-component score TSV files.
  - `dag` (string): GO DAG edge hierarchy file (`child\tancestor`).
  - `batchSize` (integer, default `16`): Batch size.
  - `minScore` (float, default `0.1`): Minimum score threshold.
- `esm2DeepGoPlus`:
  - `enabled` (boolean, default `false`): Run ESM-2 representation with DeepGO-Plus classification head.
  - `checkpoint` (string): Path to PyTorch model checkpoint (`.pt`).
  - `terms` (string): GO term vocabulary file (one term per line).
  - `model` (string, default `esm2_t33_650M_UR50D`): Pretrained ESM-2 backbone architecture.
  - `batchSize` (integer, default `16`): Batch size.
  - `minScore` (float, default `0.1`): Minimum score threshold.
- `proteinfer`:
  - `enabled` (boolean, default `false`): Run ProteInfer CNN model.
  - `modelDir` (string): Directory containing ProteInfer SavedModel assets.
  - `batchSize` (integer, default `16`): Batch size.
  - `minScore` (float, default `0.1`): Minimum score threshold.
- `clean`:
  - `enabled` (boolean, default `false`): Run CLEAN contrastive EC number predictor.
  - `modelDir` (string): Directory containing CLEAN checkpoints and scripts.
  - `batchSize` (integer, default `16`): Batch size.
  - `minScore` (float, default `0.1`): Minimum score threshold.
- `esm2Centroid`:
  - `enabled` (boolean, default `false`): Run ESM-2 embedding distance to precomputed functional centroids.
  - `db` (string): NPZ database containing centroid vectors and GO terms.
  - `model` (string, default `esm2_t12_35M_UR50D`): ESM-2 embedding model.
  - `topK` (integer, default `5`): Centroid nearest neighbors to keep.
  - `batchSize` (integer, default `16`): Batch size.
  - `minScore` (float, default `0.2`): Minimum score threshold.

#### Open-source auxiliary predictors (`predictors.foss`)

Configures auxiliary region, term, and site predictors mediated by sidecar scripts.

- `regionSidecar` (string): Path to `run_region_predictors.py`.
- `termSidecar` (string): Path to `run_term_predictors.py`.
- `siteSidecar` (string): Path to `run_site_predictors.py`.
- `deepFri`: `enabled` (boolean, default `false`), `modelDir` (string), `mode` (string, default `seq`), `minScore` (float, default `0.5`).
- `mdf` (metagenomic-deepFRI): `enabled` (boolean, default `false`), `weightsDir` (string), `mdfExecutable` (string), `skipPdb` (boolean, default `true`), `foldcompDb` (string), `threads` (integer, default `4`), `minScore` (float, default `0.1`).
- `deepEc`: `enabled` (boolean, default `false`), `modelDir` (string), `minScore` (float, default `0.5`), `acknowledgeAgpl` (boolean, default `false`).
- `deepArg`: `enabled` (boolean, default `false`), `modelDir` (string), `type` (string, default `prot`), `minScore` (float, default `0.5`).
- `musiteDeep`: `enabled` (boolean, default `false`), `modelDir` (string), `residueTypes` (string, default `Phosphoserine_Phosphothreonine`), `minScore` (float, default `0.5`).
- `scanNet`: `enabled` (boolean, default `false`), `modelDir` (string), `structureDir` (string), `minScore` (float, default `0.5`).

#### Viral predictors (`predictors.viral`)

Configures genomic phage and provirus predictors consuming contig nucleotide FASTA:

- `genomicSidecar` (string): Path to `run_genomic_predictors.py`.
- `genomad`: `enabled` (boolean, default `false`), `dbPath` (string), `sif` (string), `minScore` (float, default `0.7`).
- `checkv`: `enabled` (boolean, default `false`), `dbPath` (string), `sif` (string), `threads` (integer, default `4`), `minScore` (float, default `0.5`).
- `phispy`: `enabled` (boolean, default `false`), `sif` (string), `trainset` (string), `minScore` (float, default `0.5`).

### Quality assessment and constraint enforcement (`quality`)

Controls completeness, coherence, consistency metrics, and deterministic post-prediction repair.

| Key | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Compute genome-scale quality metrics when `goOwlFile` is present. |
| `goOwlFile` | string | `null` | Path to Gene Ontology OWL or OBO file. |
| `scope` | string | `contig` | Evaluation granularity: `contig` (per contig), `genome` (pooled across contigs), or `both`. |
| `provenance` | boolean | `true` | Record provenance trail for every annotation and write an enforcement action log. |

#### Completeness (`quality.completeness`)

| Key | Type | Default | Description |
|---|---|---|---|
| `profile` | string | `auto` | Essential function profile name: `auto`, `bacteria`, `archaea`, `eukaryote`, `virus`, or `skip` (omit completeness scoring). |
| `addTerms` | list of string | `[]` | Additional essential GO terms added to the profile for this run. |
| `removeTerms` | list of string | `[]` | GO terms excluded from the completeness evaluation profile. |
| `customFile` | string | `null` | Path to a custom essential functions TSV file (overrides `profile`). |
| `enforce` | boolean | `false` | Impute missing essential functions onto the strongest candidate protein. |
| `promoteEvidence` | string | `ISC` | Evidence code stamped on imputed essential annotations. |
| `promoteMinScore` | float | `0.05` | Minimum sub-threshold predictor score required for essential function promotion. |

#### Coherence (`quality.coherence`)

| Key | Type | Default | Description |
|---|---|---|---|
| `process` | boolean | `true` | Evaluate process coherence across GO `has_part` relations. |
| `pathway` | boolean | `true` | Evaluate pathway completeness across KEGG and MetaCyc modules. |
| `complex` | boolean | `true` | Evaluate multi-subunit coherence of obligate heteromeric protein complexes. |
| `enforce` | boolean | `false` | Deterministically resolve complex singletons and missing `has_part` partners. |
| `promotePartner` | boolean | `true` | Permit promoting sub-threshold candidate partners for complex singletons. |
| `enforceProcess` | boolean | `true` | Promote missing `has_part` partner terms during coherence enforcement. |
| `promoteEvidence` | string | `ISC` | Evidence code stamped on promoted coherence annotations. |
| `promoteMinScore` | float | `0.05` | Minimum sub-threshold predictor score required to promote a partner subunit. |

#### Consistency (`quality.consistency`)

| Key | Type | Default | Description |
|---|---|---|---|
| `taxonConstraints` | boolean | `true` | Evaluate taxonomic constraint consistency using SAT reasoning. |
| `taxonomyFile` | string | `null` | Path to taxonomy hierarchy and disjointness TSV file (defaults to bundled backbone). |
| `constraintsFile` | string | `null` | Explicit taxon constraints source (`.obo` or TSV) when using `go-basic.obo`. |
| `inferTaxonMinScore` | float | `0.9` | Minimum prediction score required for a GO term to inform organism taxon inference. |
| `enforce` | boolean | `false` | Execute SAT taxon constraint enforcement post-pass. |
| `enforceMode` | string | `remove` | Action for violating annotations: `remove`, `downrank`, `flag`, or `minimal-flip`. |
| `downrankFactor` | float | `0.5` | Multiplicative penalty applied to scores when `enforceMode` is `downrank`. |
| `organismTaxon` | string | `null` | Assert organism NCBI taxon (e.g. `NCBITaxon_2`, `2`, or kingdom name). |

### Output section (`output`)

| Key | Type | Default | Description |
|---|---|---|---|
| `formats` | list of string | `['gff3', 'tsv']` | Annotation file formats to generate. |
| `htmlReport` | boolean | `true` | Generate an interactive HTML quality report. |
| `outputDir` | string | `gspa_output` | Directory where output files are written. |

### Integration section (`integration`)

Controls the Bayesian evidence integration and iterative coordinate ascent engine.

| Key | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `false` | Master switch for multi-evidence integration. |
| `maxIter` | integer | `6` | Maximum coordinate ascent refinement iterations. |
| `epsilon` | float | `0.005` | Convergence tolerance on mean absolute posterior delta across function keys. |
| `damping` | float | `0.5` | Jacobi under-relaxation factor between iterations. |
| `thetaFile` | string | `null` | Path to a learned hyperparameters JSON file (`theta.json`). |
| `defaultsFile` | string | `null` | Path to a YAML defaults override file. |

Sub-configs within `integration`:

- `darkMatter`:
  - `enabled` (boolean, default `false`): Enable dark-matter operon gap suggester.
  - `mode` (string, default `separate`): Output mode: `separate`, `top-q`, or `distributed`.
  - `bfMin` (float, default `10.0`): Minimum Bayes factor required to nominate an operon.
  - `gammaInP` (float, default `50.0`): Likelihood elevation multiplier for in-pathway candidate functions.
  - `coverageThreshold` (float, default `0.9`): Credible set coverage threshold for disjunctive nominations.
  - `satDisambiguation` (boolean, default `false`): Enable SAT-based constraint disambiguation.
  - `useLmContext` (boolean, default `false`): Incorporate genomic language model context.
- `intragenomeCluster`:
  - `enabled` (boolean, default `false`): Cluster proteome and propagate predictions from centroids.
  - `identity` (float, default `0.9`): Sequence identity threshold for clustering (90%).
  - `coverage` (float, default `0.8`): Alignment coverage threshold (80%).
- `outerLoop`:
  - `enabled` (boolean, default `false`): Run outer fixed-point loop between gapseq and dark-matter suggester.
  - `maxIter` (integer, default `5`): Maximum outer loop iterations.
  - `qBase` (float, default `0.5`): Base posterior cutoff for reaction pinning.
  - `qStep` (float, default `0.05`): Posterior step increment per outer cycle.
  - `qCap` (float, default `0.75`): Maximum posterior ceiling for pinning.
  - `tauCover` (float, default `0.5`): Minimum pathway coverage threshold.
  - `pinPromotions` (boolean, default `true`): Prevent demotion of pinned annotations across outer cycles.

## Bundled profiles comparison

These are legacy kingdom/configuration profiles; the new `--preset` workflows are described in [workflow presets](../scenarios/presets.md).

GSPA ships with four legacy preset profiles in `gspa-core/src/main/resources/config/`. Selecting a profile with `-k` or `--kingdom` preconfigures biological assumptions:

| Parameter | `bacteria.yaml` | `archaea.yaml` | `eukaryote.yaml` | `mag.yaml` |
|---|---|---|---|---|
| `input.kingdom` | `bacteria` | `archaea` | `eukaryote` | unset (`auto`) |
| `input.type` | `genome` | `genome` | `genome` | `mag` |
| `predictors.similarity.enabled` | `true` (`diamond`) | `true` (`diamond`) | `true` (`diamond`) | `true` (`diamond`) |
| `predictors.domains.applications` | Pfam, TIGRFAM, CDD, SUPERFAMILY | Pfam, TIGRFAM, CDD, SUPERFAMILY | Pfam, TIGRFAM, CDD, SUPERFAMILY | Pfam, TIGRFAM, CDD, SUPERFAMILY |
| `predictors.operons.enabled` | `true` (300 bp) | `true` (300 bp) | **`false`** | `true` (300 bp) |
| `predictors.pathway.tool` | `gapseq` | `gapseq` | `gapseq` | `gapseq` |
| `predictors.pathway.crossfeeding` | `false` | `false` | `false` | **`true`** |
| `predictors.localization.enabled` | `true` | `true` | `true` | `false` |
| `quality.completeness.profile` | `bacteria` | `archaea` | `eukaryote` | `auto` |
| `quality.coherence` | process, pathway, complex | process, pathway, complex | process, pathway, complex | process, pathway, complex |
| `quality.consistency.taxonConstraints` | `true` | `true` | `true` | `true` |

Key biological distinctions:
- **Eukaryote profile (`eukaryote.yaml`)**: Operon detection is disabled because polycistronic transcription units are absent in canonical eukaryotic nuclear genomes. The completeness profile uses the eukaryotic essential set.
- **MAG profile (`mag.yaml`)**: Metabolic crossfeeding analysis is enabled to model incomplete metabolic networks across metagenomic communities. Completeness defaults to `auto`, which falls back to the bacterial essential set; this does not dynamically infer the domain. Use `--preset mag` for conservative defaults without an assumed essential set.

## Command-line overrides

Arguments passed on the command line override values defined in YAML files or presets:

| CLI Option | Overridden YAML Key |
|---|---|
| `-k`, `--kingdom <name>` | `input.kingdom` |
| `--mag` | `input.type = 'mag'` |
| `-o`, `--output <dir>` | `output.outputDir` |
| `--go-owl <file>` | `quality.goOwlFile` |
| `--no-metrics` | `quality.enabled = false` |
| `--metrics-scope <scope>` | `quality.scope` |
| `--db`, `--database <file>` | `predictors.similarity.database` |
| `--enforce-consistency` | `quality.consistency.enforce = true` |
| `--enforce-consistency-mode <mode>` | `quality.consistency.enforceMode` |
| `--taxon-constraints <file>` | `quality.consistency.constraintsFile` |
| `--taxonomy <file>` | `quality.consistency.taxonomyFile` |
| `--taxon <id>` | `quality.consistency.organismTaxon` |
| `--infer-taxon` | `quality.consistency.taxonConstraints = true` |
| `--infer-taxon-min-score <score>` | `quality.consistency.inferTaxonMinScore` |
| `--enforce-completeness` | `quality.completeness.enforce = true` |
| `--enforce-coherence` | `quality.coherence.enforce = true` |
| `--no-provenance` | `quality.provenance = false` |
| `--disorder` | `predictors.disorder.enabled = true` |
| `--neural-sidecar <path>` | `predictors.neural.sidecarScript` |
| `--base-predictor <mode>` | `predictors.neural.basePredictor` |
| `--esm2-deepgoplus` | `predictors.neural.esm2DeepGoPlus.enabled = true` |
| `--proteinfer` | `predictors.neural.proteinfer.enabled = true` |
| `--clean` | `predictors.neural.clean.enabled = true` |
| `--deepgo-plusplus-light` | `predictors.neural.deepGoPlusPlusLight.enabled = true` |
| `--deepgo-plusplus-light-assets <dir>` | `predictors.neural.deepGoPlusPlusLight.assets` |
| `--min-score <score>` | `predictors.neural.deepGoPlusPlusLight.minScore` (also filters output) |
| `-t`, `--threads <count>` | `predictors.neural.deepGoPlusPlusLight.threads` |

## Complete annotated example

Below is a complete standalone configuration file illustrating all common options:

```yaml
# GSPA comprehensive pipeline configuration example

input:
  format: "fasta"
  type: "genome"
  kingdom: "bacteria"
  translateAltStart: true

predictors:
  similarity:
    enabled: true
    tool: "diamond"
    database: "/data/databases/uniprot_sprot.dmnd"
    evalue: 1.0e-5
    maxTargetSeqs: 10

  structure:
    enabled: false
    tool: "foldseek"
    database: "/data/databases/afdb_swissprot"
    centroidMode: "none"

  domains:
    enabled: true
    tool: "interproscan"
    applications:
      - "Pfam"
      - "TIGRFAM"
      - "CDD"
      - "SUPERFAMILY"

  antifam:
    enabled: true
    database: "/data/databases/AntiFam.hmm"
    useCutGa: true

  operons:
    enabled: true
    maxIntergenicDistance: 300
    requireSameStrand: true

  pathway:
    enabled: true
    tool: "gapseq"
    crossfeeding: false
    exportSbml: false
    runFba: false

  localization:
    enabled: true
    signalP: false
    tmhmm: false
    deepSig: true
    tmbed: true
    deepSigKingdom: "gramn"

  disorder:
    enabled: true
    minRegionLen: 10
    minScore: 0.5

  neural:
    pythonExecutable: "python3"
    sidecarScript: "/opt/gspa/benchmark/neural/run_neural_predictors.py"
    basePredictor: "light"
    deepGoPlusPlusLight:
      enabled: true
      assets: "/data/dgpp_light_assets"
      diamond: "diamond"
      threads: 16
      minScore: 0.1

quality:
  enabled: true
  goOwlFile: "/data/ontologies/go.owl"
  scope: "both"
  provenance: true

  completeness:
    profile: "bacteria"
    enforce: true
    promoteEvidence: "ISC"
    promoteMinScore: 0.05

  coherence:
    process: true
    pathway: true
    complex: true
    enforce: true
    promotePartner: true
    enforceProcess: true

  consistency:
    taxonConstraints: true
    enforce: true
    enforceMode: "remove"
    organismTaxon: "NCBITaxon_562"

output:
  formats:
    - "gff3"
    - "tsv"
  htmlReport: true
  outputDir: "annotated_genome"
```
