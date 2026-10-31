# Reaction context modes and complete search acquisition

Status: implementation, mechanism validation, matched real-cohort exports and
model comparison pass their execution checks. Biological benefit is inconclusive. No whole audit requirement is
accepted. [Acceptance record](../../tasks/v3-acceptance.md).

The production engine now exposes universal-network, named-pathway and
organism-restricted context through the shared workflow, CLI and Nextflow.
[Input formats and semantics](../concepts/cohort-context.md) describe their
differences. Candidate generation remains common to all modes.

Pathway alternatives form one evidence group. Overlapping groups are merged,
including indirect alternatives of the focal reaction, which cannot anchor it.
Tests cover local observations, peer states, joint coverage and shuffled
background windows. Organism traversal uses each genome's own set, prevents
paths through excluded intermediates, and keeps target signature availability
separate from the peer union. Known empty inputs remain distinct from missing
inputs. Neither mode supplies a flux-feasibility claim.

`cohort_context_mode` accompanies exported vectors without changing numerical
contract v3. Legacy input-support schema 1 covers universal context only.
Schema 2 binds class coverage to mode as well as source, function, role and
availability. The JVM and Python reject unseen modes; new-mode data cannot
borrow another mode's training or calibration class coverage. Mode survives
table preparation, serialization, training identity and evaluation. Historical
universal-only training identities and support manifests remain unchanged.

## Completed validation

- 1,245 JVM tests pass, with two unrelated optional skips. All four live
  HMMER/AntiFam checks execute in this gate.
- 610 benchmark tests pass, with three optional skips and all ten live
  Nextflow tests executed. The new test checks direct/Nextflow pathway parity
  and invalidates a cached task after pathway mappings change.
- 41 trainer tests and 13 input-support tests pass. The input-support tests
  include mode-conditioned coverage, serialization and model/locus abstention.
- The default real-cohort control preserves all 484,202 candidate identities,
  source identities, numerical features and functional prediction TSV bytes
  relative to the preceding signature-availability engine.

Development failures remain in the evidence: a transitive-alternative leak,
omitted mode in training identity, missing report field, and updated export
fixtures. Test-fixture errors are distinguished from production defects. A
broader gate exposed openpyxl overwriting a supposedly fixed workbook modified
timestamp; the fixture now normalizes the serialized timestamp. Metabolic
production code was unchanged.

## Complete upstream search acquisition

Reference acquisition job 10120 completed with exit 0 after 1:25:23 and verified
46,477 files. Retry array 10131 completed the ten archaeal searches. The combined
audit now verifies all 18 frozen genomes, 49,463 proteins, 175,272 reaction rows,
19,424 pathway rows and 2,247 positive pathway predictions. The original ten
archaeal preflight failures remain in the audit history.

The context-input preparer reads no function labels. It preserves 9,636 named
reaction groups and audits 145 IDs absent from the production graph and 4,732
groups with no usable graph mapping. Organism sets contain mapped strong
non-complex hits and spontaneous reactions. Complex activities require
independent subunit inference and are omitted from those sets. These are
observational hypotheses, not proof of enzyme function or feasible flux.

The first prepared input included a reaction-table ID with empty stoichiometry;
the production loader excludes that row. The batch rejected the inconsistency.
The corrected preparer now matches the production loader's complete set of
34,703 reaction IDs. Both input versions and the failed batch are retained.

## Matched comparison protocol

`prepare_context_modes.py` binds search observations, mappings and metadata to
their hashes. `run_context_modes.py export` runs a default compatibility control
and three matched modes, with the same homology candidates, layouts, family
assignments, holdouts and prepared metadata. Only the context relation changes
between the three matched arms. Every arm must retain candidate/source identities,
non-cohort numerical values and functional TSV bytes.

The separate training stage preserves the preceding 4,262 supervision records,
splits and fitting settings, excluding the mode metadata and changed features
from the equality comparison. Logistic and tree models are fitted separately
per arm; thresholds use calibration only. `check_context_mode_models.py` compares
actual JVM model scores and abstentions with Python on the same exported rows.
Results must be reported as reference-support comparisons on an already examined
cohort, not independent biological validation.

All 72 genome/arm outputs verify: the default compatibility control and the
three matched modes each retain 484,202 candidates. Adding organism metadata
to the universal arm changes three organism-availability/fraction columns.
Pathway and organism modes each change 36 cohort columns relative to that
matched universal arm. Representative reaction selection can also change the
reported specificity of the strongest supporting EC call; this affects
`cohort.w.ec_reaction_count` for 166 pathway-arm and 171 organism-arm candidates and is retained in
the column-level audit. The candidate universe remains unchanged.

Both models are refitted separately for all three modes, with all 4,262
supervision records and splits preserved. Supported test coverage is 462/608
for universal, 463/608 for pathway and 461/608 for organism context. Every null
layer remains untrained: the training data contain zero no-correct-candidate
groups and one single-positive group.

The common support intersection contains 460 rows from four test genomes:

| Context | Logistic AUROC | Logistic Brier | Tree AUROC | Tree Brier |
|---|---:|---:|---:|---:|
| Universal | 0.754933 | 0.123525 | 0.716409 | 0.127080 |
| Pathway | 0.752073 | 0.124997 | 0.704895 | 0.127837 |
| Organism | 0.761921 | 0.122568 | 0.714182 | 0.128049 |

These pooled rows weight genomes unequally. The mean per-genome Brier difference
relative to universal is +0.000670 for logistic pathway (paired 95% bootstrap
interval −0.001199 to +0.002860) and +0.000541 for logistic organism
(−0.003227 to +0.005769). Tree differences are +0.000687
(−0.000750 to +0.001476) and +0.001387 (−0.001151 to +0.004241).
The bootstrap resamples genomes in pairs with 10,000 replicates and seed
20261027; all intervals include zero. This descriptive comparison establishes
no reliable context benefit. No mode is promoted to a new default on these results.

The actual JVM model API and Python agree on support/abstention and calibrated
probabilities for all 25,572 row-model combinations (4,262 rows × six models).
Maximum absolute probability difference is 6.11e−16. This verifies model execution,
not calibration transfer or biological correctness.

The previous reference-overlap audit confirmed 395 identical target proteins,
including 93 test proteins, and was a lower-bound audit of bacterial searches.
The archaeal acquisition does not resolve that limitation. Full reference/family
exclusion, independent enzyme/role/null truth, biological context ablations,
organism feasibility and the remaining original audit requirements stay open.

## Evidence

Archive: `benchmark/v3/results/context-modes-v1/`, with 1,281 payloads
independently decompressed and rehashed (16,270,693,645 original bytes;
959,252,822 stored payload bytes).
The acquisition stage contains the complete retry observations, execution and
scheduler receipts, reference hashes, and the combined audit. Large payloads
use deterministic gzip and are decompressed and rehashed against their original
SHA-256 values. Later stages retain exports, failed attempts, frozen sources,
test logs, model artifacts and comparisons.

Production JAR SHA-256:
`7f9d1d477098087469f90564d79a923ee28e85ba0a010c33840f2db44131671b`.
Preceding baseline JAR:
`6757088967769a6d9204f50b34a7e96ebf42a177e4ae50588db27a562cdc3dd7`.
