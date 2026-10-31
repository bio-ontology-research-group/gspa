# Controlled nucleotide loss in cohort inference

This addresses the missing-sequence control in recommendations 6 and 8. It does
not complete either recommendation or establish calibrated absence attribution.
The earlier `fragmented_assembly_50pct` result removed genomic coordinates while
keeping all protein claims and orthogroup memberships. It is a coordinate-loss
control, not a sequence-loss experiment.

## Inputs and intervention

`benchmark/v3/prepare_sequence_loss.py` reads the acquired assemblies, protein
FASTAs, coordinate layouts, predictor claims, and cohort metadata/memberships.
It does not read labels or the labelled-protein selection table. Each original
contig is divided into 10,000-base blocks. SHA-256 of seed 42, genome, contig and
block start gives a reproducible order; the lowest half of block hashes in each
genome are removed. No fragment spans two original contigs. Short terminal blocks
remain eligible, so half the blocks need not mean exactly half the nucleotides.

Only CDSs wholly inside a retained block are kept. A CDS crossing any cut is
excluded even if both adjacent fragments survive: these controls reuse complete
reference translations rather than predicting partial proteins. The script writes
the retained nucleotide fragments, protein FASTAs and remapped coordinates, filters
the claims and orthogroups to the surviving proteins, and sets completeness to the
actual retained nucleotide fraction. Family identities remain fixed to isolate
observation loss; this is not de novo family clustering or gene calling.

The 16 originally unplaced proteins remain available, with their original claims,
and are explicitly outside the spatial intervention. Unknown position is not
converted into sequence absence. The preparation checks that each original
protein is either placed once or explicitly unplaced and rejects undeclared
claims or orthogroup members. Outputs are immutable; input/output hashes and every
removed protein's reason are recorded in `preparation.json`.

For the 18-genome cohort, preparation retains 22,457 of 49,463 proteins and 214,706
of 484,202 original claims. The difference between nucleotide coverage and protein
coverage includes CDSs intersecting fragment boundaries.

## Matched production comparison

`benchmark/v3/run_sequence_loss_benchmark.py run` invokes the actual batch CLI
with the same model, references, clade holdout and cohort settings for three arms:

| Arm | Coordinates | Protein claims and family membership | Completeness |
|---|---|---|---|
| `intact` | Original | Original | Original |
| `coordinate_loss` | Surviving fragment coordinates; others explicitly unplaced | Original | 1.0: protein evidence remains complete |
| `sequence_loss` | Same surviving fragment coordinates | Lost proteins removed | Measured retained nucleotide fraction |

This is a cohort comparison. No metabolic gaps are supplied, so reaction-context
inference is not requested and no reaction-context ablation is claimed. Every
batch requires content-verified production outputs and executed cohort/learned
stage ledgers. Previous outputs cannot certify changed inputs or source code.

## Evaluation boundary

The evaluator keeps genome and protein IDs throughout. It does not pool identical
protein sequences: a prediction for a surviving paralogue or another genome must
not cover a removed locus. It rejects sequence-loss outputs that reintroduce a
removed protein.

The separate `calibrate` command chooses a threshold for each annotation type
using only intact-arm calibration genomes. The declared rule requires at least
ten calls and 90% support in the experimental reference annotations, maximizes
supported calls, then supported fraction, then the threshold. Unsupported
thresholds abstain. The policy is immutable and bound to inputs and run receipts
before `evaluate` scores test genomes. GO annotations and predictions propagate
through the frozen ontology.

Test results separate all original labelled proteins, retained labelled proteins,
and removed labelled proteins. The original denominator is preserved. Reports
include supported and unverified calls, new/lost supported calls, reference recall,
and genome-bootstrap intervals. EC top-1 reference agreement and paired changes
are also reported without a threshold; those are descriptive ranking metrics,
separate from the frozen call policy. Unmatched annotations are **unverified**,
not experimentally disproved functions. No false-assignment reduction is inferred
from incomplete positive-only labels. GO ancestor counts are dependent annotations,
not independent biological observations.

This is the previously examined development cohort. Fixed seeds and separated
calibration prevent additional outcome-dependent tuning but do not turn it into
a new blind test. Missing sequence is a known intervention here; recovering that
explanation by the inference model remains a separate acceptance requirement.

## Executed results (2026-10-26)

All 54 genome/arm runs completed with verified execution receipts and executed
cohort and learned stages. The independent input check confirmed 2,606 nucleotide
fragments against exact source slices, 22,457 unchanged protein translations,
214,706 byte-identical retained claims, and 22,457 retained family memberships.
Preparation removed 22,565 proteins in deleted blocks and 4,441 CDSs intersecting
block boundaries. Retained DNA fractions range from 49.87% to 50.57% by genome.

The frozen EC threshold is 0.44: calibration has 25 supported calls out of 27
(92.6%). **That support fraction does not transfer to test genomes.** On the intact
test input, only 6 of 12 calls are reference-supported; the other six remain
unverified. No GO threshold meets the declared calibration requirement, so the
GO call policy abstains for every arm. These results do not establish reliable
90%-precision annotation and were not used to retune the policy.

| Test measure | Intact | Coordinate loss | Sequence loss |
|---|---:|---:|---:|
| EC top-1 reference matches, all 37 labelled proteins | 14/37 | 12/37 | 6/37 |
| EC top-1 reference matches, 13 retained labelled proteins | 7/13 | 6/13 | 6/13 |
| EC top-1 reference matches, 24 removed labelled proteins | 7/24 | 6/24 | 0/24 |
| Supported EC calls at the frozen threshold | 6 | 3 | 1 |
| Unverified EC calls at the frozen threshold | 6 | 1 | 0 |

No deleted protein reappears in the sequence-loss production output. Its decrease
in all-original top-1 reference agreement is 0.2162 (paired genome-bootstrap 95%
interval for the change: −0.4167 to −0.0500). On retained labelled proteins the
change is −0.0769 (−0.2000 to 0.0000). Only three of the four test genomes have
experimental EC labels; these intervals have limited independent-genome support.
No improvement in gene/function recovery or experimentally established reduction
in false assignment is demonstrated.

Eleven targeted preparation/evaluation tests and the full benchmark Python suite
pass: 474 tests, three optional skips. The tests include identical-sequence locus
separation, lost-protein denominators, calibration/test separation, malformed
input rejection, immutable preparation and the actual route's input bindings.
Source/output hashes, preparation, input verification, all production receipts,
stage coverage, frozen policy and evaluation are archived under
[`benchmark/v3/results/sequence-loss-v1/`](../../benchmark/v3/results/sequence-loss-v1/).
Bulky assemblies and production outputs remain in the cache paths named by the
receipts. Recommendations 6 and 8 remain open for the other declared controls,
independent biological validation and inferred absence attribution.

## Reproduction

Use a new input directory for a new preparation. `run` reuses results only after
checking complete execution receipts; `calibrate` and `evaluate` require new
policy/report paths.

```bash
export TMPDIR="$HOME/tmp"
export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$TMPDIR"
PY="$HOME/tmp/gspa-v3-testenv/bin/python"
CACHE="$HOME/.cache/gspa-v3-data/training-cohort"
WORK="$HOME/.cache/gspa-v3-data/sequence-loss-v1"
$PY benchmark/v3/prepare_sequence_loss.py --cache "$CACHE" --out "$WORK/inputs"
$PY benchmark/v3/run_sequence_loss_benchmark.py run --inputs "$WORK/inputs" --work "$WORK/runs"
$PY benchmark/v3/run_sequence_loss_benchmark.py calibrate --inputs "$WORK/inputs" --work "$WORK/runs" --policy "$WORK/policy.json"
$PY benchmark/v3/run_sequence_loss_benchmark.py evaluate --inputs "$WORK/inputs" --work "$WORK/runs" --policy "$WORK/policy.json" --out "$WORK/report.json"
$PY -m unittest discover -s benchmark/v3/tests -p test_sequence_loss.py -v
```
