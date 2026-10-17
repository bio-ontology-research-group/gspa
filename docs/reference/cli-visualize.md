# gspa visualize

`gspa visualize` writes a single self-contained HTML file for browsing one
GSPA workspace: integrated annotations with their supporting predictors,
operons, quality metrics, optional sidecar predictor results and an igv.js
genome browser. The command extracts a Python script bundled in the jar
(`visualize/make_viz.py`) and runs it; Python 3.8 or later is required and only
the standard library is used.

```text
gspa visualize -w WORKDIR [options]
```

## Example

```bash
gspa visualize -w workspaces/MR59-6 \
  --go-obo ref/go.obo --ec2go ref/ec2go.txt \
  --fasta workspaces/MR59-6/input/MR59-6_assembly.fa
# writes workspaces/MR59-6/gspa_out/MR59-6_browser.html
```

## Options

| Option | Type | Default | Description |
|---|---|---|---|
| `-w`, `--workdir` | directory | required | Workspace directory (layout below). Must exist. |
| `-r`, `--run-dir` | string | `gspa_out` | Subdirectory of the workspace that holds the integrator outputs. |
| `-g`, `--genome-id` | string | workspace directory name | Genome identifier used to build sidecar file names (`<id>.mdf.tsv`, `<id>.proteinfer.tsv`, ...). |
| `-o`, `--out` | file | `<workdir>/<run-dir>/<genome-id>_browser.html` | HTML output path. (The help text says `<workdir>/<genome>_browser.html`; the code places the file inside the run directory.) |
| `--go-obo` | file | see note | `go.obo` used for GO term names and aspects. |
| `--ec2go` | file | see note | EC-to-GO file used for EC names. |
| `--fasta` | file | `<workdir>/input/<genome-id>_assembly.fa` | Assembly FASTA embedded (base64) for the igv.js genome browser tab. |
| `--python` | string | `python3` | Python interpreter used to run the script. |
| `-h`, `--help` | | | Show help and exit. |
| `-V`, `--version` | | | Print the version and exit. |

!!! note "Defaults for `--go-obo` and `--ec2go`"
    The help text mentions `$GSPA_REF/go.obo` and `$GSPA_REF/ec2go.txt`, but
    the bundled script does not read `GSPA_REF`. When the options are omitted,
    the script falls back to site-specific absolute paths compiled into it.
    Always pass both options explicitly.

The options are passed to the script as environment variables
(`GSPA_WORKDIR`, `GSPA_RUN_DIR`, `GSPA_GENOME_ID`, `GSPA_OUT`, `GSPA_GO_OBO`,
`GSPA_EC2GO`, `GSPA_FASTA`; see [Environment variables](environment.md#launcher-environment-gspa)).
If the script exits with a non-zero status, the command fails.

## Workspace layout

Paths are relative to `--workdir`; `<id>` is the genome identifier and `<run>`
the run directory.

| Path | Required | Content |
|---|---|---|
| `prokka_out/<id>.tsv`, `prokka_out/<id>.gff` | yes | Per-CDS metadata and coordinates (Prokka/Bakta-style TSV and GFF). |
| `<run>/integrated.tsv` | yes | Output of `gspa integrate --out` ([format](file-formats.md#integrated-tsv)). |
| `<run>/provenance.json` | yes | Output of `gspa integrate --provenance` ([format](file-formats.md#provenance-json)). |
| `<run>/operons.tsv` | optional | Operons written by `gspa annotate` ([format](file-formats.md#annotate-output-directory)). |
| `<run>/quality_gspa.json` | optional | JSON quality report from `gspa evaluate` ([format](file-formats.md#evaluate-json-report)). |
| `mdf_out/<id>.mdf.tsv`, `proteinfer_out/<id>.proteinfer.tsv`, `clean_out/<id>.clean.tsv` | optional | Sidecar predictor outputs (`protein_id`, `term`, `score`, `annotation_type`). |
| `ipr_out/<id>.faa.tsv` | optional | InterProScan TSV. |
| `amrfinder_out/<id>.amr.tsv` | optional | AMRFinderPlus results. |
| `antismash_out/regions.js` | optional | antiSMASH biosynthetic gene clusters. |
| `input/checkm.tsv`, `input/busco.tsv`, `input/gtdbtk.tsv`, `input/culture_conditions.txt` | optional | Genome metadata shown in the header. |
| `input/<id>_assembly.fa` | optional | Default location of `--fasta`. |

See [Interactive browser](../scenarios/visualize.md) for a walkthrough.
