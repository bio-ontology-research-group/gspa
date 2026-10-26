# Reading the prediction ledger

The results browser (`gspa visualize`, and the web portal's "Explore predictions") has a
**Prediction ledger** tab. It shows what each integration stage did to the predictions, using the
stage ledger the shared workflow writes next to the integrated table:

```
gspa_out/integrated.tsv                 final (protein, function) posteriors
gspa_out/integrated.tsv.stages.json     per-stage status and probability changes
gspa_out/*suggestions*.tsv              optional: candidate alternatives and not-assigned probability
```

The web portal offers the ledger as **Prediction ledger (JSON)** and in the job status. The browser
also downloads the filtered rows as TSV and the complete ledger as JSON.

## What the ledger does and does not say

- A change is a change to a *prediction*. The ledger never says a changed prediction is correct or
  more accurate. Accuracy needs an independent benchmark with truth labels; nothing here supplies one.
- **Baseline.** Today the baseline is the output of evidence integration, which already includes
  priors and context. It is **not** a sequence-only baseline, and the browser never labels it as one.
  A future ledger can set `baseline_kind: "sequence_only"`; only then does the browser use that label.
- **Missing ledger.** A legacy run without `integrated.tsv.stages.json` shows *unavailable*. It is
  not shown as zero changes or as success. An unreadable ledger is treated the same way.

## Stage status

Each stage is `executed`, `skipped` (with a reason such as `disabled`), `failed` (with the error),
or `not_run`. Only executed stages report change counts. A skipped, failed or not-run stage shows its
reason instead of zeros. `evidence_integration` creates the baseline; its candidates are reported as
"created", not as changes.

## Counting rules

Every count in the tab is computed from the rows listed in the table.

- Counts are **unique (protein, function) pairs**, with unique **proteins** in parentheses. A protein with
  three changed GO terms counts as three functions and one protein.
- Changes are `added`, `removed`, `increased` or `decreased`. Removed and downweighted predictions are
  listed alongside additions.
- The view "Evidence-integrated baseline → final" is the net change. Per-stage views list what that
  stage did; a function changed by two stages appears in both views and once in the net view.
- The status line states whether the ledger's final candidate count equals the exported
  `integrated.tsv` row count and whether every net final probability matches the export. A mismatch is
  reported, not hidden.

## Enzyme classification

Enzyme-related changes use the function itself, not its name:

| Class | Rule |
|---|---|
| Complete EC number | an EC number with four levels (`EC:2.7.1.1`; preliminary `n` numbers included) |
| Partial EC number | an EC number with a missing level (`EC:1.1.-.-`); shown separately from complete EC |
| GO catalytic activity | the GO term is `GO:0003824` or a descendant of it through `is_a` in `go.obo` |
| GO complex | a descendant of `GO:0032991` (protein-containing complex); **not** an enzyme, even when its name says "catalytic" |
| Other GO | the hierarchy is known and the term is neither catalytic nor a complex |
| Unknown | the GO hierarchy is missing, or the term is absent from it |

Without `go.obo` hierarchy, GO terms are **unknown**; the browser does not fall back to name matching.
EC numbers still classify. A protein's enzyme flag in the Proteins tab follows the same rules, using
integrated scores of 0.5 or higher.

## Evidence, alternatives and the not-assigned outcome

- **Per protein.** The detail panel adds a *Prediction trace*: the stage changes for that protein,
  the sources that supported each function, the priors that fired, and other functions of the same
  protein with their scores.
- **Sources and priors** are recorded separately. A prior is context, not an independent
  observation of the protein's function.
- **Alternatives.** Suggestion tables list candidate proteins with probabilities and a separate
  **not-assigned** probability. Candidates join to proteins by the FASTA sequence ID. IDs that are not
  genome proteins are listed as not joined, not matched by approximate name. A suggestion that did not
  report a not-assigned probability shows "not reported", not zero.
- **Role scores** are shown when the ledger or suggestion table carries role metadata. The browser does
  not compute them.
- **Placeholder product names** such as "hypothetical protein" are not functional evidence and are
  noted as such in the trace.
- **Metabolic assumptions.** When the FBA stage ran, its medium and biomass template appear beside it.
  Feasibility is conditional on those assumptions; it is not an observation.

## Tests

`scripts/tests/test_visualize_ledger.py` and `service/test_ledger_summary.py` use small genomes with
interpretable functions: unchanged enzymes beside non-enzyme changes, a catalytic GO term with no
enzyme-like name, a "catalytic complex" that is not an enzyme, removed and downweighted enzymes,
failed and skipped stages, and table counts that must equal the TSV and ledger. Browser counting and
filter code is run under Node against the exported rows.

```bash
python3 -m unittest discover -s scripts/tests -p 'test_visualize_ledger.py'
python3 -m unittest service.test_ledger_summary
```
