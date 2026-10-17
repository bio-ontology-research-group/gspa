# gspa evaluate

`gspa evaluate` measures the genome-scale quality of an existing annotation:
completeness of essential functions, process, pathway and complex coherence,
taxon consistency, and mean information content. It does not run any
predictor. The metrics are defined in
[Genome-scale quality metrics](../concepts/quality-metrics.md).

```text
gspa evaluate -i GENOME.gff3 --go-owl go.owl [-a annotations.gaf] [options]
```

## Example

```bash
gspa evaluate \
  -i genome.gff3 -a predictions.gaf \
  --go-owl go.owl --reasoner-cache cache/ \
  -k bacteria --taxon NCBITaxon_562 \
  --ec2go ec2go.txt --pathways kegg_pathways.tsv --modules kegg_modules.tsv \
  -o reports/genome_quality.json
```

For a quick run without the ELK reasoner (no process coherence):

```bash
gspa evaluate -i genome.gff3 -a predictions.gaf --go-owl go.owl --lite \
  -o reports/genome_quality.tsv
```

## Inputs

| Option | Type | Default | Description |
|---|---|---|---|
| `-i`, `--input` | file | required | GFF3 file with the genome structure. Every CDS feature becomes a protein, identified by its GFF3 `ID`. |
| `-a`, `--annotations` | file | none | GO annotations in GAF 2.x format (plain or `.gz`). Column 2 (`DB_Object_ID`) must match the CDS `ID` in the GFF3. Lines with a `NOT` qualifier and lines with fewer than 15 columns are skipped; the annotation score is derived from the evidence code. Without this option no proteins are created and the metrics are computed on an empty annotation. |
| `--go-owl` | path | required | GO in OWL format. |
| `-k`, `--kingdom` | string | `bacteria` | Selects the essential-function profile (`bacteria`, `archaea`, `eukaryote`/`eukarya`/`eukaryota`, `virus`/`viruses`; anything else uses the bacterial profile). When `--taxon` is not given, the kingdom is also asserted as the organism taxon for the consistency check (`bacteria`, `archaea`, `eukaryote` and `eukaryota` are mapped to NCBI taxa; other strings are passed through unchanged). |
| `--taxon` | string | none | Organism taxon asserted during the consistency check: `NCBITaxon_562`, `NCBITaxon:562`, `562`, or one of `bacteria`, `archaea`, `eukaryote`, `eukaryota`, `fungi`, `yeast`, `metazoa`, `viridiplantae`, `plant`. Overrides the taxon implied by `-k`. |
| `--mag` | flag | off | Evaluate completeness with the MAG adjustment (expected incompleteness of metagenome-assembled genomes). |

## Reference data

| Option | Type | Default | Description |
|---|---|---|---|
| `--complex-terms` | file | bundled `complex_terms.tsv` | Complex GO term classification ([format](file-formats.md#complex-terms-tsv)). Used for complex coherence. A path that does not exist is ignored. |
| `--ec2go` | file | none | EC-to-GO mapping ([format](file-formats.md#ec2go)). With `--pathways` it is used to resolve pathway steps; alone it builds a pathway database from the EC-to-GO map. |
| `--pathways` | file | none | Pathway definitions TSV ([format](file-formats.md#pathways-and-modules-tsv)). Used only together with `--ec2go`. |
| `--modules` | file list (comma-separated, repeatable) | empty | Additional pathway TSVs stacked on `--pathways` (KEGG modules, MetaCyc, BioCyc). Used only when both `--ec2go` and `--pathways` exist. |
| `--taxonomy` | file | bundled hierarchy | Two-column taxonomy hierarchy (`child<TAB>parent`). See the warning below. |

!!! warning "`--taxonomy` in evaluate"
    In this release, passing `--taxonomy` to `evaluate` creates the taxon
    checker before the bundled GO taxon constraints are loaded, and those
    constraints are then not loaded at all. Consistency is effectively not
    checked. Omit `--taxonomy` unless you also verify the violation count.

## Performance

| Option | Type | Default | Description |
|---|---|---|---|
| `--lite` | flag | off | Skip the ELK reasoner. Process coherence is not computed and is reported as `N/A` (-1 in JSON). |
| `--reasoner-cache` | directory | none | Cache for reasoner results (`has_part` pairs); speeds up repeated full runs. |

## Output

| Option | Type | Default | Description |
|---|---|---|---|
| `-o`, `--output` | file | `quality_report.json` | Report path. A name ending in `.json` produces the JSON report ([evaluate JSON](file-formats.md#evaluate-json-report)); any other name produces a one-line summary TSV ([format](file-formats.md#quality-summary-tsv)). |

If the consistency check finds violations, a violations TSV is written next to
the report, named by replacing the report's extension with `_violations.tsv`
(for example `genome_quality_violations.tsv`; see
[format](file-formats.md#violations-tsv)). A summary is also printed to the
console.

## Standard options

| Option | Description |
|---|---|
| `-h`, `--help` | Show help and exit. |
| `-V`, `--version` | Print the version and exit. |
