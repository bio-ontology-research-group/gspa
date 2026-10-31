# Complete TopEC, MAPred and EnzBERT reference panels

The complete 22-genome TopEC, MAPred and EnzBERT panels show **no EC top-1 gain
from default priors**. Job 10306 scores the 44 TopEC/MAPred cells and job 10313
the 22 EnzBERT cells from the frozen original `price-production-v1` collection.
Together they cover 13,682,870 acquired whole-proteome input claims.
The evaluation compares the original exported scores, production integration
without priors, and production integration with default priors. This evaluates
the original JAR, not completion of the newer compact-audit collection.

| Acquired predictor input | Raw top-1 reference matches | No-prior integration | Default priors | Reference functions present among candidates |
| --- | ---: | ---: | ---: | ---: |
| MAPred | 50/148 (33.8%) | 50/148 | 50/148 | 114/151 |
| TopEC | 5/148 (3.4%) | 5/148 | 5/148 | 5/151 |
| EnzBERT | 53/148 (35.8%) | 53/148 | 53/148 | 96/151 |

Both integration comparisons have zero top-1 gains, losses or changed top-1 IDs.
There are no top-score ties in these three panels. MAPred matches a reference
within its top three on 95/148 proteins; TopEC remains at 5/148 and EnzBERT reaches
76/148. The 2,000-resample genome bootstrap gives descriptive top-1 intervals of
26.5–41.0% for MAPred, 0.65–7.0% for TopEC and 28.0–42.9% for EnzBERT. Paired
top-1 differences are identically zero. These intervals describe
this reused dataset, not independent validation or guaranteed probability coverage.

The result supports two specific limitations. Existing default priors do not
change the top predictions on these direct EC inputs. In addition, improving the
ranking of the existing TopEC candidates alone cannot recover the other 146
reference functions: they were never proposed in these acquired outputs. Neither
observation is a conclusion about all GSPA modules or the predictors' published
benchmarks. Additional candidate sources, joint context, trained confidence and
the remaining predictor panels still require their own matched comparisons.

## Population and identifier policy

All 148 eligible proteins and 151 original reference functions are retained,
including missing candidates. The original 149th entry, `NP_384885`, remains an
explicit exclusion because its mapped identifier `NC_003047.1` is a replicon.
Every eligible protein ID is checked against the pinned genome sequence inventory,
independent of whether any predictor emitted a candidate. These are known enzyme
targets; other proteins lack accuracy labels, so this is not enzyme detection.

Exact EC identifiers are the primary analysis. The sensitivity analysis uses the
[pinned ENZYME audit](v3-ec-identity.md), merging only single-destination transfers
that reach an active entry. Ambiguous, deleted and unknown IDs remain under their
original identities. Projected reference duplicates are merged and candidate
scores use the maximum, never a sum. Original and projected denominators are
reported separately. These three panels contain only resolvable target/reference
ECs. TopEC/MAPred have no merged candidate identities; EnzBERT has three merged
candidate groups across three reference proteins. Max-merging those aliases
changes none of the reported summary metrics. The sensitivity results therefore
match the primary results. Other predictor panels may behave differently.

The evaluator retains preliminary four-level ECs as candidates. Its primary
top-1 metric gives fractional credit across tied highest scores, with worst/best
tie bounds. Lexical top-k, reciprocal rank and lexical gain/loss counts are
supplementary; reciprocal-rank bounds are also included. Raw exported scores have
no imposed probability cutoff and are not native tool-selected calls. At the
fixed diagnostic integrated cutoff of 0.5, MAPred matches 61 reference functions,
leaves 90 unmatched and makes 136 unconfirmed calls; TopEC matches two, leaves
149 unmatched and makes 83 unconfirmed calls. EnzBERT matches 41 functions, misses
110 and makes 34 unconfirmed calls. Both prior settings agree for each predictor.
These nonmatches are not experimentally proven false functions, and 0.5 is not a
validated confidence threshold.

## Execution and verification

`evaluate_price_reference.py` requires every genome for every selected predictor.
An incomplete panel, contradictory status, changed receipt, different arm inputs,
unexpected command difference or changed full artifact fails evaluation. Commands
may differ only in the declared prior setting, arm output directory and private
Java scratch path. The normalized batch manifests must also agree. The reader
hashes and validates every raw/prediction row while retaining only reference-protein
scores in memory. The TopEC/MAPred run rehashes all 792 bound inputs and artifacts;
the EnzBERT run rehashes 506, including the full ledgers in each run, before
publication.

All 15 focused tests pass locally and on node004. They cover missing candidates,
transfer splits, deleted/unknown/preliminary ECs, score ties, max-merge semantics,
the full 22-genome fixture, missing panels, altered ledgers, wrong protein IDs and
different arm inputs. A read-only Sonnet 5.5 review and its dispositions are
retained. Remote evaluation completes with exit 0 in 108.93 seconds and 35,840 KiB
peak RSS; no inference is rerun or changed.

The [portable evidence archive](../../benchmark/v3/results/price-reference-v1/)
contains sources, review, remote tests, registration, full reference score extracts
and reports. Its independent arithmetic checker verifies all 1,776
protein/arm/identity-policy records and all 44 panels without importing the
evaluator. It checks reduced-table arithmetic; re-reading the original large
remote files remains the recorded evaluation's responsibility. Run:

```bash
python3 benchmark/v3/results/price-reference-v1/verify_archive.py
```

Historical Price outcomes were already examined before this policy was written.
No threshold or model was fitted here, but this is development reuse, not a blind
test or a pre-registered accuracy result. CLEAN, other remaining predictors,
ensembles, learned/context arms and independent biological validation remain
unfinished. No whole acceptance item is signed off and no release is made.

Follow-up job 10313 completes the entire EnzBERT panel after original production
array 10229 succeeds. It uses the same frozen evaluator, retains all 148 eligible
proteins and 151 reference assertions, and rehashes 506 bound files covering
11,553,500 acquired claims. All 15 tests pass; execution finishes with exit 0 in
8:45.86 at 34,400 KiB peak RSS. The [separate 41-file EnzBERT archive](../../benchmark/v3/results/price-reference-enzbert-v1/README.md)
independently checks 888 protein/arm/identity-policy records and confirms zero
top-1 identity changes directly from retained scores. Its portable checker has
the same reduced-table scope as the first archive. Run:

```bash
python3 benchmark/v3/results/price-reference-enzbert-v1/verify_archive.py
```

EnzBERT has 55 reference functions absent from its acquired candidate set; these
cannot be recovered by reordering that set alone. Its 53/148 top-1 result does
not establish superiority to another predictor: this report has not tested a
paired between-predictor hypothesis and does not reproduce native tool protocols.
The three completed panels establish a negative result for these default-prior
comparisons. Dense remaining panels, ensembles and learned/context comparisons
are still required.

## Full six-predictor comparison queued

Job **10432** now waits for full collection finalizer **10410** before evaluating
the newer `price-production-v3` collection. It uses the exact same five pinned
evaluator, metric, identifier, test and reference files as the completed panels.
The scope is all **132 cells**, six predictors × 22 genomes, and **1,526,919,081**
acquired input claims. All **148 eligible proteins and 151 reference functions**
remain in each predictor's denominator, including missing predictions. This
compares raw exported rankings, no-prior integration and default-prior integration;
it does not run native predictors or learned/context inference.

The launch requires the original registration and plan hashes, successful
finalizer exit, matching finalizer/collection receipts, and the exact 132-cell
verification inventory. The evaluator subsequently reads and hashes every score
row and checks every bound file. A separate existing arithmetic checker will
verify the reduced tables after scoring. Source and package hashes are checked
before execution. The job has no automatic requeue and is cancelled if its
dependency fails. It does not alter running inference jobs or their frozen code.

All 15 focused evaluator tests pass locally. Sonnet 5.5 found no blocking issue
in the two launch files. Deployment checks confirm the package hashes, log
directory and Python 3.12.3. A read-only metadata preflight checks all **73**
completed cells present at that checkpoint, including matched arm commands,
receipt identities and file availability. That preflight neither rehashes the
large score files nor scores an incomplete panel. Job 10432 remains pending;
**there is no new accuracy result or acceptance sign-off**.

The launch package, review, tests and scheduler receipt are retained in
[`price-reference-full-launch-v1`](../../benchmark/v3/results/price-reference-full-launch-v1/README.md).
