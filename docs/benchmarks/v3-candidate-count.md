# Candidate-count discontinuity repair

The legacy `IterativeRefiner` previously applied its weak-evidence context ceiling
only when the run contained at least 100 distinct protein/function candidates.
An unrelated 100th candidate could therefore lower an existing candidate's score.
Small mechanism fixtures exercised different scoring rules from whole-genome runs.

The ceiling now applies whenever the candidate has negative likelihood log-odds
and receives a positive total prior contribution, regardless of the number of
other candidates. Strong-evidence candidates and negative prior contributions
retain their previous behavior. The existing rule for runs of 100 or more
candidates is unchanged. The ceiling remains a heuristic; this repair does not
train weights, calibrate confidence, or supply the missing null/role labels.

## Evidence

A new core regression first failed on the old implementation: with identical weak
evidence and a fixed positive prior, the target scored 0.9276544687893495 for 1/99
candidates and 0.5240973543916296 for 100/101/1,000. The repaired implementation
returns the latter score at every size and preserves the convergence iteration.
Three additional cases cover strong positive and negative context contributions.
The failing test output is retained in the evidence archive.

`benchmark/v3/check_candidate_count.py` separately exercises the actual CLI and
shared workflow. It creates synthetic EC claims, a declared metabolic gap and
`alpha_gap=4`, then runs the old and new JARs with only the gap-filling prior.
The extra candidates have a different function and no shared protein identity.
Both batch invocations have content-verified inputs, outputs and stage ledgers.

| Total candidates | Old CLI target score | Repaired CLI target score |
|---|---:|---:|
| 1 | 0.9609545058791815 | 0.5243365144833768 |
| 99 | 0.9609545058791815 | 0.5243365144833768 |
| 100 | 0.5243365144833768 | 0.5243365144833768 |
| 101 | 0.5243365144833768 | 0.5243365144833768 |

The complete 100- and 101-candidate TSVs are byte-identical between builds. No
biological truth labels are used: this is a scoring-invariance regression, not a
measurement of annotation accuracy.

The full JVM suite passes **972 tests**, with two skips and no failures/errors;
the production JAR builds successfully. The benchmark suite passes **523 tests**,
with three optional skips, including all eight live Nextflow checks against the
rebuilt JAR. The learned-context overview now reflects actual development-model
training, the version-2 feature contract, unsupported-input behavior and the
untrained null layer.

## Reproduction and identities

```bash
export TMPDIR="$HOME/tmp"
export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$TMPDIR"
export GSPA_TEST_PYTHON="$HOME/tmp/gspa-v3-testenv/bin/python"
./gradlew test :gspa-cli:shadowJar --no-daemon
"$GSPA_TEST_PYTHON" benchmark/v3/check_candidate_count.py \
  --before-jar "$HOME/.cache/gspa-v3-data/engines/c469273224723bff080efd9440d1c1ae07fc0c943c2cacffa04a98003481f801.jar" \
  --after-jar gspa-cli/build/libs/gspa-3.0.0.jar \
  --out "$TMPDIR/candidate-count-replay"
```

Use a new output directory; existing results are never overwritten. The replay
requires an old build that exhibits the defect and refuses to certify one that
does not reproduce it. Its old/new JAR hashes are recorded in `verification.json`.
The verified repaired JAR hash is
`58dbebf1cf6ad97cd216d676b0cf1221600a97c7350f6f45abf27b54c9669654`.
Both JARs are retained in the local content-addressed engine cache.

Source snapshots, before/after predictions, receipts and test logs are under
`benchmark/v3/results/candidate-count-v1/`. Earlier capacity jobs continue with
their immutable original JAR; they were not switched mid-run. Item 4 remains open
for its training, calibration, negative/null and independent-validation requirements
in `tasks/v3-acceptance.md`.
