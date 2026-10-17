# gspa report

`gspa report` renders HTML from JSON quality reports produced by
`gspa evaluate`, `gspa annotate` (`*_quality.json`) or `gspa compare`
(`.json` output).

```text
gspa report [-o OUTPUT.html] REPORT1.json [REPORT2.json ...]
```

## Example

```bash
gspa report -o quality.html out/*_quality.json
```

## Options

| Option | Type | Default | Description |
|---|---|---|---|
| `-o`, `--output` | file | `quality_report.html` | HTML output. With one input file a single-genome report is written; with several, a multi-genome comparison. |
| `<jsonFiles>...` | files (positional) | none | JSON quality reports ([format](file-formats.md#evaluate-json-report)). If none are given, an error message is printed and the command exits with status 0. |
| `-h`, `--help` | | | Show help and exit. |
| `-V`, `--version` | | | Print the version and exit. |

## Fields read from each JSON file

| JSON field | Used as |
|---|---|
| `genome_id` | Genome label (file name if absent). |
| `assessment_date` | Assessment date. |
| `completeness.score` | Completeness. |
| `completeness.present_functions`, `completeness.missing_functions` | Present and missing essential functions. |
| `coherence.process_coherence`, `pathway_coherence`, `complex_coherence` | Coherence scores. |
| `consistency.consistent` | Consistency flag. |
| `information_content.mean_ic` | Mean information content. |
| `summary.total_proteins`, `summary.annotated_proteins` | Coverage. |

Detail sections of the JSON (named terms, pathway detail, violations) are not
used by this command.

!!! note "Known limitations in 2.1.0"
    Values that are false or zero are treated as missing when the JSON is
    read: `consistent: false` is rendered as consistent, and a coherence score
    of exactly 0 is rendered as not computed. Check the JSON or the
    `evaluate` console summary for these cases.
