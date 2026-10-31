# Duplicate-strain permutation repair

The [verified control replay](v3-cohort-controls.md) found that exact strain
copies changed learned confidence. A label-free synthetic replay isolates four
affected feature columns: raw peer count, background-window count, permutation
p-value and permutation excess. This repair addresses the randomization mechanism;
raw-count effects and full biological calibration remain open.

## Failure and correction

The old permutation loop drew each genome's matched background independently,
then averaged those draws within the strain cluster. Adding identical copies
therefore reduced null variance even though the independent biological evidence
had not increased. It also consumed extra random draws that changed other units'
sampling. With five copies, the synthetic p-value fell from 0.0995025 to 0.00497512.

The loop now draws one seed per strain cluster and permutation. Each genome in
that cluster replays that stream for its corresponding matched member pools.
Exact copies therefore share their randomization, and adding them cannot consume
other clusters' draws. This couples dependent observations within a strain cluster;
it does not establish a calibrated phylogenetic null for arbitrary heterogeneous
lineages or validate the supplied clade/strain assignments.

The existing duplication test previously excluded every permutation feature from
its invariant checks. It now includes them. Four new cases cover one, three and
five copies under two seeds. All five checks failed against the original engine
and pass after repair, alongside the cohort inference and workflow tests.

## Verification

The rebuilt engine is
`a792849b501e7ad7a9989488ea3b2ce2bd85242de44b0f279b213272bd5a8547`.
All 1,050 JVM tests pass with two optional skips after enabling the four live
HMMER/AntiFam cases. The original full gate had six skips; the separate predictor
rerun closes those four live-tool checks without changing the engine.

The same standalone synthetic replay now gives permutation p-value 0.1044776
for zero, one, three and five copies; permutation excess is also invariant.
Only the two raw-count columns change. Their combined contribution to the fitted
model margin still increases by 0.0153528 for five copies. These are feature-delta
attributions, not calibrated predictions for synthetic proteins.

The benchmark suite passes 538 tests with three optional skips, including eight
live Nextflow checks. The real 59-genome/arm replay completed using the same
prepared inputs and fitted model. Independent hash and stage checks pass for
all 59 outputs. `check_cohort_permutation_replay.py` verifies both sets of receipts
and compares candidate coverage and scores without reading function labels or
fitting a new decision policy.

All 484,202 original candidate identities survive in each arm. Changing the
engine changes 5,124 intact, 4,898 shuffled and 5,183 duplicate-arm prediction
scores (absolute tolerance 1e-12); the maximum change is 0.002215. Both increases
and decreases occur. These are algorithm effects, not measured improvements.

Across all 18 original genomes, duplication's maximum learned-probability change
falls from 0.00196057 to 0.00164723. Both builds compare the same 346,269 numeric
probabilities and 137,933 unsupported pairs, with no missing candidates or support
transitions. These denominators cover all 18 genomes; the earlier control report's
0.00173679 maximum covered only four test genomes. The remaining nonzero change
means full score invariance still fails.

The unchanged model's calibration is not certified for the corrected feature
distribution. Removing count dependence and refitting under a versioned feature
contract remain necessary before broader learned-inference acceptance.

Evidence: `benchmark/v3/results/cohort-permutation-v1/`. The earlier results and
JAR remain unchanged in `benchmark/v3/results/cohort-controls-v2/`.

## Mechanism replay

```bash
export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$HOME/tmp"
java -Xmx2g -cp "$HOME/.cache/gspa-v3-data/engines/a792849b501e7ad7a9989488ea3b2ce2bd85242de44b0f279b213272bd5a8547.jar" \
  groovy.ui.GroovyMain benchmark/v3/check_cohort_duplicate_features.groovy \
  gspa-core/src/test/groovy/gspa/integration/cohort/CohortFixtures.groovy \
  benchmark/v3/results/model-input-support-v1/models/model_logistic.json \
  "$HOME/tmp/cohort-permutation-replay.json"
```

## Holdout-independent permutation seed (`content:v2`)

### Defect

The base seed of each (family, reaction) permutation null was derived from
the configured seed, the family, the reaction and the holdout key of the focal
genome (`seedFor(f, r, hk)`, with `hk` = `c:<clade>` under `CLADE` and
`s:<cluster>` under `STRAIN_CLUSTER`). The committed v3.0.0 source (`498560f`)
shows this form. Two holdouts that select exactly the same peers therefore drew
different Monte Carlo nulls. As a result, permutation p-values and excesses
depended on the holdout label even when the evidence was identical. The
per-cluster replay described above does not remove this dependence, because it
starts from the same seed.

The [WP4 cohort-recurrence benchmark](v3-wp4-cohort-recurrence.md) exposed it.
In R0 no two cohort genomes share a genus, so `STRAIN_CLUSTER` and `CLADE` select
the same peers, and the R0-SC control requires identical features. The first run
failed that control. Sealed numeric columns were identical, but candidate
feature records differed in `cohort.w.perm_excess`. Because WP4 compares a
`STRAIN_CLUSTER` regime with a `CLADE` regime, this noise would have entered its
main contrast through the learned permutation features.

### Correction

`CohortInference` now records `permutationSeed: content:v2` in run and record
provenance. The seed is the first eight bytes of a SHA-256 over the scheme tag,
the configured seed, the family, the reaction and a canonical signature of
exactly what the permutation consumes. Per evidence unit (after the
evidence-identity collapse, sorted by clade and then cluster) the signature holds
its clade id, cluster id, exact weight and the sorted background-pool
descriptors of its members. Members are described by their pools, not by genome
id. The holdout label, the peer order and the choice of representative for a
collapsed exact copy therefore cannot change the seed. Whole-genome identity
fingerprints are deliberately not used, because they also change with unrelated
claims elsewhere in a genome. All fields are length-prefixed.

`gspa-core/src/test/groovy/gspa/integration/cohort/CohortPermutationSeedSpec.groovy`
adds seven cases. They check that:

- the holdout label does not change the null when the peer sets are identical;
- a different peer set gives a different permutation stream;
- genome and claim input order change no output;
- an exact duplicate collapsed by evidence identity keeps the null (three
  cases: the duplicate sorting before or after the original, under `CLADE` and
  `STRAIN_CLUSTER`);
- the configured seed still changes the null.

The cases fail on the old code. The JAR bundle `ec-prefix-identity-v4` was built
from the 10402 snapshot plus this patch (unimatrix01 job 10699). The focused
spec passes 7/7, and the full suite passes 1,480 tests with 0 failures and 6
skips. All WP4 runs were repeated on this bundle; their R0-SC control now passes
exactly.

### Consequence for earlier results

Every engine before `ec-prefix-identity-v4` used the hold-out-keyed seed. That
includes the engine of the repair above and JAR 10402. Results produced with
these engines include WP1 and the earlier cohort control, permutation, context
and replay benchmarks. They have not been rerun with the content-keyed seed.
Whether any of their reported numbers change is **unverified**. The new seed
changes the permutation draws, and with them the permutation features, of every
candidate that has a permutation null, in every holdout regime. Whether that
moves a reported metric is unknown until each result is rerun. In WP4, besides
the R0-SC control, one unconfirmed-call count changed between the two runs. That
observation does not transfer to other benchmarks. Results within one benchmark
remain internally comparable when all their arms ran on the same engine.
