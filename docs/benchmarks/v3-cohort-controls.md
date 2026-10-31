# Verified neighborhood and duplicate-strain controls

The repaired controls run production inference on the original 18-genome cohort,
then freeze thresholds on the three calibration genomes before evaluating the four
test genomes. They replace the historical perturbation runner, which reused an
unverified baseline and selected Fmax thresholds on test labels. This is still a
previously examined development cohort, not a new blind biological test.

## Intervention and execution

`prepare_cohort_controls.py` preserves all 49,463 original proteins and claims.
It permutes protein IDs within each original contig while preserving its coordinate
and strand slots. All 49,447 placed proteins retain their contig membership;
49,418 change position. The 16 unplaced proteins remain unplaced. Preparation
does not read function labels. An independent checker verifies these properties,
FASTA coverage, and all 61 input and 20 output hashes.

The duplicate arm adds five exact copies of training genome GCA_000007005.1. Each
uses the original claims and layout, copies all 2,994 family memberships, and keeps
the same clade (`ncbi_genus_2100760`) and strain cluster (`ncbi_species_2287`).
Original metadata and memberships remain unchanged. This tests exact copies,
not near-identical strains with different assemblies or annotations.

All three arms use the same coverage-aware logistic model, references, clade
holdout and engine. The intact and shuffle arms each produce 18 genomes; the
duplicate arm produces 23. All 59 genome/arm outputs have verified execution
receipts and executed cohort and learned stages. Only the original genomes enter
evaluation: duplicate copies never add truth observations. No metabolic gaps are
supplied, so reaction-context suggestion is disabled and no effect is claimed.

Engine SHA-256: `c0def1656fc17170869b1d8f25b80e1bf433ecbc583794f63c27aff23b1bc821`.
The model is archived at `benchmark/v3/results/model-input-support-v1/models/model_logistic.json`.

## Frozen evaluation and results

The threshold rule requires at least ten calls and 90% support in experimental
reference annotations on intact calibration genomes. It maximizes supported calls,
then support fraction, then threshold. EC selects 0.44 (25/27 supported); GO has
no qualifying threshold and abstains. The policy binds inputs, source and all three
production receipts. Genome and protein identities remain distinct throughout.

| Test measure | Intact | Shuffled neighborhoods | Five strain copies |
|---|---:|---:|---:|
| EC top-1 reference matches | 14/37 | 14/37 | 14/37 |
| Supported EC calls at 0.44 | 6 | 7 | 6 |
| Unverified EC calls at 0.44 | 6 | 5 | 6 |
| New / lost supported EC calls versus intact | 0 / 0 | 2 / 1 | 0 / 0 |
| GO calls under the frozen policy | 0 | 0 | 0 |

The 90% calibration support does not transfer to test genomes. Unmatched calls
are unverified, not experimentally disproved. The shuffle's net supported EC
change is +1, with paired genome-bootstrap 95% interval −3 to +3 (2,000 draws).
Only three test genomes have EC reference labels. These results do not establish
an improvement from intact neighborhoods. GO abstention does not establish GO
prediction invariance; ancestor annotations are also correlated observations.

The duplicate comparison preserves all 156,101 test candidate identities:
114,697 have numeric learned probabilities in both arms and 41,404 are unsupported
in both. There are no support transitions. The maximum probability difference is
**0.0017367852562789055**, so exact score invariance fails even though the reported
EC decisions remain unchanged. The current model uses nonzero coefficients for
raw peer counts, background-window counts and permutation statistics; the engine's
existing duplication unit test explicitly excludes these quantities. A label-free synthetic replay through the production engine and feature bridge
confirms that these four columns change: peer count, background-window count,
permutation excess and permutation p-value. Five identical copies change the
p-value from 0.0995025 to 0.00497512, and their combined contribution to the
fitted model margin increases by 0.018582. This is a mechanism finding, not a
calibrated synthetic prediction. These channels require repair before claiming
confidence is invariant to duplicated evidence. No threshold was retuned
after observing this result.

The subsequent [permutation repair](v3-cohort-permutation.md) preserves this
original result and corrects the duplicate-sensitive randomization. Raw-count
dependence in the fitted model remains unresolved.

## Tests and evidence

Three pre-fix evaluator regressions reproduced null-probability crashes, omitted
candidate coverage, and acceptance of duplicate/malformed candidate records.
The repaired comparator separately reports numeric comparisons, unsupported pairs,
support transitions and missing identities; it does not represent an entirely
unsupported comparison as zero probability change.

Ten focused control tests and eleven sequence-loss regressions pass. The complete
benchmark suite passes 538 tests with three skips, including eight live Nextflow
checks on the frozen engine. No JVM implementation was changed in this control
replay. The last JVM gate remains 1,046 tests, two skips, zero failures/errors.

`benchmark/v3/results/cohort-controls-v2/` contains immutable preparation, policy,
results, receipts, stage coverage, independent checkers, source snapshots, failing
before-fix tests and passing gates. Bulk inputs and outputs remain under
`~/.cache/gspa-v3-data/cohort-controls-v2/`. Historical coordinate-loss results
remain distinct from the [genuine nucleotide-loss control](v3-sequence-loss.md).
Recommendations 1, 6 and 8 remain open.

## Reproduction

Use fresh output directories. Preparation and policy refuse overwrite; changed
sources or inputs invalidate receipt reuse. Keep the recorded engine and model
fixed when reproducing this particular result.

```bash
export TMPDIR="$HOME/tmp"
export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$TMPDIR"
export GSPA_JAR="$HOME/.cache/gspa-v3-data/engines/c0def1656fc17170869b1d8f25b80e1bf433ecbc583794f63c27aff23b1bc821.jar"
PY="$HOME/tmp/gspa-v3-testenv/bin/python"
CONTROL="$HOME/.cache/gspa-v3-data/cohort-controls-replay"
$PY benchmark/v3/run_cohort_controls.py prepare --inputs "$CONTROL/inputs" --source-genome GCA_000007005.1
$PY benchmark/v3/run_cohort_controls.py run --inputs "$CONTROL/inputs" --work "$CONTROL/runs"
$PY benchmark/v3/run_cohort_controls.py calibrate --inputs "$CONTROL/inputs" --work "$CONTROL/runs" --policy "$CONTROL/policy.json"
$PY benchmark/v3/run_cohort_controls.py evaluate --inputs "$CONTROL/inputs" --work "$CONTROL/runs" --policy "$CONTROL/policy.json" --out "$CONTROL/report.json"
$PY benchmark/v3/results/cohort-controls-v2/check_inputs.py "$CONTROL/inputs"
$PY benchmark/v3/results/cohort-controls-v2/check_outputs.py "$CONTROL/runs"
```
