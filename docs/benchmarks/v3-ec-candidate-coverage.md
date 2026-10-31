# EC candidate coverage and routine full-reference reporting

The missing EC functions in the [context comparison](v3-context-comparison-results.md)
were absent before integration. This controlled comparison uses one DIAMOND
homology-transfer source, not the complete GSPA predictor ensemble. Across all
335 experimental-reference EC assertions, raw claim identities and scores
reconstruct exactly from the fenced reference annotations and retained search
hits. No missing reference candidate is explained by an integration cutoff.

## Complete trace

The diagnostic covers every positive `ECO:0000269` EC assertion in the registered
18-genome cohort: 232 train, 64 calibration and 39 test assertions. Each assertion
is traced through the sequence-identity map, reference annotation inventory,
retained search hits, family filter and raw claim. Both modes are retained, giving
670 trace rows. All 756,568 claims across 36 genome/mode files are scanned and
checked against the previously bound hashes; all consumed inputs are rechecked.
For the full test panel, raw candidate membership also matches every context and
scorer's integrated output in the independently audited production run.

| Split | Reference assertions | Absent from reference | No matching retained hit | Homology emitted | Removed by family filter | Family-excluded emitted |
|---|---:|---:|---:|---:|---:|---:|
| train | 232 | 59 | 43 | 130 | 77 | 53 |
| calibration | 64 | 30 | 2 | 32 | 21 | 11 |
| test | 39 | 19 | 5 | 15 | 10 | 5 |

“Absent from reference” means that the exact EC has no qualifying annotation in
the fenced reference library. It does not mean the activity is absent from the
organism, from all Swiss-Prot records, or from other predictors. The reference
excludes the cohort taxa and cohort/Price genera, and transfers only experimental
annotations. Either experimental GO or EC annotations can qualify an entry for
inclusion; an entry's non-experimental EC annotations are not thereby transferred.
The trace does not determine how much each exclusion rule contributes.

“No matching retained hit” means none of the saved hits carries the reference EC.
It does not distinguish insufficient sensitivity, the e-value threshold, the
retrieval cap or lack of useful homology. DIAMOND retrieves at most ten targets
before the family filter removes hits with at least 30% identity and 80% query
coverage. Of 26,605 queries with hits, 10,910 reach the cap; 3,193 of those contain
only family-level hits and 7,717 retain at least one remote hit. The cap therefore
limits this stress-test search space. These counts do not prove that further
remote hits exist or that a wider search would recover another reference EC.
The original comparison and thresholds remain unchanged.

The biological limitations of the EC references also remain: complex annotations
need component-role interpretation, several assertions share the same complex,
and unconfirmed predictions are not proven errors. These are candidate-coverage
and reference-support results, not isolated-enzyme discovery accuracy.

## Reporting change

`analyze_observed_gap_context_models.py` now uses comparison schema 2 and includes
`full_ec_reference` under each threshold policy. The population comes directly
from the registered labels and genome splits, independently of the emitted
supervision rows. Every reference protein remains, including proteins with no
candidates. Every reference EC remains, including functions never proposed.
Test genomes with zero EC reference coverage are reported explicitly as zero
reference assertions, not as known nonenzymatic genomes.

The new `ec_reference_panel.py` helper reports candidate coverage, reference
matches, unconfirmed calls, missed assertions, reference precision/recall/F1,
top-1/3/5, equal-score ties and all per-protein records. The evaluator uses exact
EC identifiers and the existing mixed GO/EC calibration cutoffs or 0.5. It does
not generate candidates, refit scores or select a better cutoff on test data.
It assumes that reference protein identifiers were validated upstream; requiring
a prediction for each protein would incorrectly remove unannotated proteins.
This diagnostic independently checks every reference identity against the pinned
sequence-to-genome map.

All supervision dependencies, including labels and splits, are now directly bound
in the analysis registration against their original fit registration. Old schema-1
analyses remain preserved and require their original source to replay. New code
refuses an old protocol rather than silently changing its population. A failed
partial report is retained; a retry uses a fresh output directory.

## Validation and remaining work

The old report fails the new full-reference assertion because it omits that output.
Nine focused tests now pass, covering missing proteins/functions, multifunction
references, unconfirmed calls, empty reference genomes, contradictory evidence,
invalid probabilities, genome populations, ties and threshold preservation.
The final threshold regression uses a 0.6 test score with a calibration cutoff of
0.9 and fixed cutoff of 0.5, so the two policies must produce different calls.

Unimatrix job 10224 replays the complete report from the original production
predictions. It verifies and stages all 328 consumed prediction/table/reference
files against the previous bindings, then runs the changed evaluator. All 18
cutoffs, 36 original metric summaries and 48 paired statistics are unchanged.
All 36 expanded EC panels exactly match the independent audit, including every
per-protein record. No prediction model is run or fitted. The production code is
unchanged after that replay; only the cutoff regression test was strengthened.

The [retained archive](../../benchmark/v3/results/ec-candidate-coverage-v1/README.md)
includes the complete trace, source inventories and hits, reference-locus raw
claims, executed source package, full replay outputs, tests, and two Sonnet
reviews with adjudication. Its independent portable verifier recomputes all trace
categories, checks source identities and confirms report parity without importing
the production evaluator. Full raw claim and integrated-TSV files remain at their
bound original locations; extracted claims alone cannot prove extraction
completeness without the original source files.

The next candidate-generation comparison must keep the genus/taxon exclusions,
freeze search breadth before evaluation, and report all reference assertions and
all additional unconfirmed calls. It should compare wider retrieval and genuinely
independent annotation channels separately, with new predictor contracts and
appropriate calibration. A cap increase alone is not a biological improvement.
The present changes repair reporting and diagnose a concrete bottleneck; they do
not establish a context benefit, superiority over CLEAN or whole-item acceptance.
All twelve requirements remain open.
