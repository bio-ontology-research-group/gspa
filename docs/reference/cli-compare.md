# gspa compare

`gspa compare` evaluates several annotation sets for the same genome with the
same quality pipeline and writes them into one report, so that methods can be
compared on completeness, coherence, consistency and information content.

```text
gspa compare -i GENOME.gff3 --go-owl go.owl [options] METHOD1.gaf [METHOD2.gaf ...]
```

## Example

```bash
gspa compare -i genome.gff3 --go-owl go.owl -k bacteria \
  -o comparison.html \
  gspa.gaf eggnog.gaf interproscan.gaf
```

Each positional GAF file is one method. The method name shown in the report is
the file name without `.gaf`, `.tsv`, `.gaf.gz` or `.tsv.gz`.

## Options

| Option | Type | Default | Description |
|---|---|---|---|
| `-i`, `--input` | file | required | GFF3 with the genome structure. CDS feature `ID`s must match column 2 of the GAF files. The GFF3 is re-read for every method. |
| `--go-owl` | path | required | GO in OWL format. |
| `-k`, `--kingdom` | string | `bacteria` | Essential-function profile used for completeness (`bacteria`, `archaea`, `eukaryote`/`eukarya`/`eukaryota`, `virus`/`viruses`; anything else uses the bacterial profile). Unlike `evaluate`, the kingdom is not asserted as the organism taxon. |
| `--lite` | flag | off | Skip the ELK reasoner (no process coherence). |
| `-o`, `--output` | file | `comparison.html` | Output. A name ending in `.html` produces a multi-genome HTML report. A name ending in `.json` writes one `<method>_quality.json` per method into the directory of the given path (the given file name itself is not written). Any other name produces a summary TSV with one row per method ([format](file-formats.md#quality-summary-tsv)). |
| `<gafFiles>...` | files (positional) | none | GAF files to compare, one per method. If none are given, an error message is printed and the command exits with status 0. |
| `-h`, `--help` | | | Show help and exit. |
| `-V`, `--version` | | | Print the version and exit. |

GAF parsing follows [`gspa evaluate`](cli-evaluate.md#inputs): `NOT` lines are
skipped and scores are derived from evidence codes.

See also [Evaluate, compare, report](../scenarios/evaluate-compare-report.md).
