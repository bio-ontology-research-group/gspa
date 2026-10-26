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

- **Baseline, schema 2.** A schema-2 ledger names its baseline in `changes_measured_from` (today
  `legacy_prior_baseline`: the posterior after the heuristic evidence combiner and priors) and explains it in
  `baseline_definitions`. The browser shows that definition. For a schema-1 ledger it falls back to
  `baseline_kind`, then to evidence integration. Stage 0 (raw maximum score) is a different baseline and is
  never mixed with it.

## Scoring mode and the learned stage

`scoring_mode` in the ledger says which scoring produced the scores: **Learned** (a frozen model replaced the
heuristic posterior), **Legacy heuristic** (no trained model), or **Unknown** (an older ledger that does not
say). The tab and the portal job status show the mode, and for learned runs the model id and hash. A
calibrated model probability is not an accuracy statement; the tab says so next to the mode.

The learned stage also writes `integrated.tsv.learned.json`. The tab reads it as follows.

- **Marginals, locus and null stay separate.** Each (protein, function, role) candidate has its own
  calibrated marginal. The *locus* probability is conditional within a function/role group, and
  *no candidate correct* is a separate outcome. The browser never sums marginals or derives an assignment
  from them. A `null` locus or no-candidate value is shown as unavailable with the report's reason
  (for example the locus layer was disabled), not as zero.
- **Roles.** Accessory, transporter, regulator, complex-subunit and unspecified roles are shown as roles,
  marked "not an enzyme assignment". Rejected or below-threshold roles are listed as *scored only*; this is
  not an improvement claim.
- **Join.** Candidates join to proteins and to the final table by the canonical function key
  `protein|TYPE|function` (the protein id may itself contain `|`). Proteins that are not in this genome are
  listed as *not joined* and cannot be opened.
- **Trace.** The protein detail panel lists the learned candidates per function with their outcome
  (`rescored`, `added`, `retracted`, `removed_no_projection`, `scored_only`), every supporting source, and
  whether the role was projected into the final role-free annotation. A learned candidate that never reached
  the final EC table is still shown.
- **Trust.** Scores are shown only when the ledger says `learned`, its `learned_integration` stage is
  `executed`, and the report's model hash equals the ledger's. A hash mismatch, a report with no hash while
  the ledger records one, a failed or missing learned stage, or a report beside a ledger that does not say the
  learned stage ran is **not trusted** and the scores are withheld. If the ledger records no hash the report
  is shown labelled *model hash unverified* (never "verified"), and the portal's download link says so. A
  missing report, or one that is not the writer's `gspa.learned-stage-report/1` object (wrong schema,
  non-object `model`, non-list `candidates`/`outcomes`), is **unavailable, not zero**. Non-numeric scalar
  fields inside an otherwise valid report are shown as unavailable rather than crashing the page.
- **Writer-to-browser test.** `LearnedArtifactsWriterSpec` runs the shared workflow on the `ContextFixtures`
  genome with a frozen, hand-set, TEST ONLY model (`intended_use = mechanism_test`) and writes the real
  artifacts; `scripts/tests/test_visualize_learned_e2e.py` reads them through the browser builder and the
  portal summary. It tests that writers and readers agree on the files. It is not a validation of a trained
  model and says nothing about accuracy. The production loader refuses such a model, so it does not go
  through `gspa integrate`.
- The candidate table is filterable by role, projection and text, and paginated (100 per page).

## Stage status

Stages appear in the order the ledger lists them (the run order), with readable labels such as *Cohort
context (two-pass)*, *Learned integration (frozen model)*, *Feature export* and *Metabolic experiments*.
Stages that prepare inputs or write side files (reference preparation, cohort context, metabolic
experiments, feature export) are marked "does not edit predictions" rather than shown with zero counts.

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
- **Role scores** are shown when the ledger, a suggestion table or the learned-stage report carries role
  metadata. The browser does not compute them.
- **Placeholder product names** such as "hypothetical protein" are not functional evidence and are
  noted as such in the trace.
- **Metabolic assumptions.** When the FBA stage ran, its medium and biomass template appear beside it.
  Feasibility is conditional on those assumptions; it is not an observation.

## Tests

`scripts/tests/test_visualize_learned.py` builds schema-2 ledger and learned-report fixtures in the shape of
the JVM writers (a contract test reads the Groovy sources so the field names cannot drift), and
`scripts/tests/test_visualize_browser.py` drives them in real Chromium (`GSPA_BROWSER_TESTS=1`, Playwright).

`scripts/tests/test_visualize_ledger.py` and `service/test_ledger_summary.py` use small genomes with
interpretable functions: unchanged enzymes beside non-enzyme changes, a catalytic GO term with no
enzyme-like name, a "catalytic complex" that is not an enzyme, removed and downweighted enzymes,
failed and skipped stages, and table counts that must equal the TSV and ledger. Browser counting and
filter code is run under Node against the exported rows.

```bash
python3 -m unittest discover -s scripts/tests -p 'test_visualize_ledger.py'
python3 -m unittest service.test_ledger_summary
```
