# Native CLEAN evidence for the training cohort

Status: acquisition and collection audit are in progress. This is not a completed
accuracy comparison, a fitted model, or acceptance of any whole v3 requirement.

The [completed context comparison](v3-context-comparison-results.md) uses
DIAMOND-derived evidence. Its learned input-support checks cannot establish
support for native CLEAN scores. The next comparison therefore needs actual
CLEAN evidence from the same complete cohort, followed by new production feature
exports and training on the existing training/calibration partitions. No new
candidate can be treated as calibrated merely because its native score lies
between zero and one.

## Registered population and execution

The acquisition keeps all 18 resolved genomes and all 49,463 protein sequences
from the existing pathway-search cohort. Its source manifest SHA-256 is
`b9491693f7b34b29abe99c7f077c44f2b70b3003db258816d487d3badf2e8fa9`.
It reads no evaluation labels or split assignments and retains every protein in
the benchmark denominator. One unsupported input receives an explicit abstention.

[The acquisition runner](../../benchmark/v3/observe_clean_cohort.py) pins the
FASTA files, native source, adapter, model weights, `split100.csv` reference,
ESM weights, and eight package versions. It uses native max-separation selection
with GMM confidence, seed zero and no confidence filtering. These are selected
EC calls, not the earlier dense Price matrices. Zero-confidence and preliminary
EC calls remain present. Native four-decimal score text is retained; no
`calibrated_prob` is invented.

Each GPU job uses its own application directory and opaque protein identifiers.
The audit joins every native CSV token to its typed TSV row, verifies the mapping
against FASTA order, then exports genome-scoped claims. Claims carry the generic
predictor contract required by learned integration and preserve the original
adapter contract separately. The original source and staged files are rehashed
before a completion receipt is written.

Remote snapshots are separate and immutable:

- `clean-cohort-v1`: job 10324 stopped before inference because its test bundle
  omitted two dependency files. The failure evidence is retained.
- `clean-cohort-v2`: the corrected seven-file source bundle. Setup job 10330
  passes eight acquisition and fourteen existing adapter tests. Package
  installation succeeds, but preflight then fails on one unsupported residue
  before inference. Original logs and source remain retained.
- The original dependent jobs 10333, 10334 and 10346 were cancelled after setup
  10330 was confirmed failed. A fresh schema-2 snapshot, `clean-cohort-v3`, reuses
  its successfully installed environment and preserves all original inputs.
  Native jobs still request one GPU, four CPUs and 40 GiB on node004.

The remote prefix is `/data/hohndor/gspa-v3/`. Job submission does not establish
successful inference. Completed native predictions are not yet reported here. Final preparation
identity and reference-overlap counts are recorded below. The initial 31-file checkpoint
records the earlier running state; the terminal outcome above supersedes it.

## Complete collection verification

[The separate collection auditor](../../benchmark/v3/audit_clean_cohort.py)
imports neither the acquisition runner nor its adapter. It requires all 18
genomes and 49,463 proteins, exactly one completed attempt per genome, no
unfinished attempts, unchanged registered inputs and exact artifact inventories.
It independently reads every FASTA identifier, native CSV call, typed TSV row
and exported claim. Missing calls, altered scores, incorrect driver identity,
added probability fields and swapped protein identities cannot be accepted as
a successful collection. Failed attempts are retained in the report.

The initial snapshot passed seventeen auditor and eight acquisition tests locally. Tests include
changed native-to-claim values despite resealed hashes, invented calibration,
missing genomes, multiple successes, unfinished attempts and duplicate/nonfinite
JSON. These are integrity checks, not biological tests. A separate Sonnet 5.5
read-only review led to exact metadata checks, independent reconstruction of
provenance identifiers, parent-population validation and command/manifest checks.
The initial nine-file audit snapshot remains frozen at `clean-cohort-audit-v1`;
its job 10346 was cancelled without running after the upstream preflight failure. A killed attempt without a terminal marker remains an error
until explicit terminal evidence supports its classification; it is not silently
ignored.

The [31-file checkpoint archive](../../benchmark/v3/results/clean-cohort-checkpoint-v1/README.md)
retains source, local and remote unit-test logs, review dispositions, scheduling,
the original setup failure and a timestamped I/O observation. Its complete
checksum inventory and both source manifests pass. It contains no completed
native cohort results; this distinction is explicit in its machine-readable
checkpoint and README.

## Interpretation and next comparison

Preparation reports exact full-sequence overlap with CLEAN's supplied reference
and equality after the native 1,022-residue truncation, including counts per
genome. It does not use reference EC labels to choose cases or exclude matching
proteins. This report is not a complete homolog or pretraining-leakage analysis.
In particular, the existing family exclusion applies to the DIAMOND transfer
reference; adding pretrained CLEAN cannot inherit a claim of family independence.

After the full acquisition passes, combine its unchanged claims with the two
existing DIAMOND transfer regimes. Preserve genome identity, supervision policy,
splits and context ablations; fit supported models on training/calibration only.
Report complete experimental-reference recall, correction and harm counts,
calibration and paired genome uncertainty, together with reference overlap and
the reused-cohort limitation. Unsupported roles and null outcomes must still
abstain. A negative result is a valid result; neither score movement nor extra
candidates alone demonstrates better enzyme-function prediction.

## Multi-source export preparation

The existing cohort checker assumes every candidate has the same single
predictor identity. [Candidate-specific sequence checks](../../benchmark/v3/cohort_sequence_inputs.py)
now support the CLEAN/DIAMOND pair through the training export checker. They
preserve complete input bytes when joining claims, retain zero-score observations,
and verify each candidate's actual source set, reference identities, raw maxima,
source/group availability and supplied role alternatives. They reject a combined
probability substituted for a raw score. An empty successful CLEAN observation
is invalid; a DIAMOND search with no hits remains valid. Genome identity and
mixed bare/UniProt ID rewrites are checked before publication.

Eighteen new tests and six existing training/dispatch tests pass locally after
a Sonnet 5.5 review. Generic pipe-containing IDs remain supported, and multiple
versions of a single source keep the old uniform-identity check. The helper does
not reconstruct GO complex classification from the ontology; it checks default
role consistency and leaves complete feature-schema validation to the existing
table builder. No joint CLEAN/DIAMOND feature export, refit or accuracy comparison
has run yet. Full native acquisition and its independent audit remain prerequisites.

## Unsupported residue and explicit abstention

The complete input scan finds exactly one unsupported residue: `J` at zero-based
position 163 of the 182-residue `AAC23084.2` protein in `GCA_000027305.1`. No other
protein among the 49,463 inputs contains a character outside the installed ESM
alphabet. The schema-2 plan retains that original sequence and records its hash,
position and `unsupported_input` / `abstain` outcome. CLEAN receives the other
49,462 sequences unchanged. The original denominator remains 49,463; the missing
prediction is not a biological negative and is not a claim that the gene has no
function. No J-to-X substitution is made.

The independent auditor reconstructs the partition from original FASTAs, verifies
valid sequence bytes and order, requires unchanged native records for unaffected
genomes, and checks supported plus unsupported counts against each original
proteome and the full cohort. Its report retains original genome hashes and flat
genome-scoped abstention outcomes. Downstream evaluation must load that original
population rather than infer its denominator from claim rows. Thirty focused
checks pass, including an affected-genome native/typed/claim audit and altered
sequence, missing abstention and wrong-count regressions. A further read-only
Sonnet review found the partition consistent and prompted the explicit denominator,
order and report-provenance checks. No native prediction or fitted model is yet
reported from this revised snapshot.

The final fresh schedule is setup 10365, pilot 10366, remainder array 10367 and
independent audit 10368. Setup has passed all 30 focused and 14 adapter tests;
setup completes with exit 0. The frozen plan SHA-256 is
`2ffa72dca65200eb8d21baeadfb8804da846630385ab87d61b5df37272c309d5`.
Pilot 10366 completes successfully with CUDA 12.4 on a Quadro K2200 and all eight
registered package versions matched. All 2,994 original proteins in
`GCA_000007005.1` produce 3,695 native calls, including 801 zero-confidence calls;
none is filtered. Its claim-file SHA-256 is
`e035c3122ab5ec927685e3d99d9575a1fa4f5f00a1e660f3304fd3f559fa2b0e`.
Array 10367 has started the remaining genomes; audit 10368 is pending.
Dependent jobs require successful predecessors;
no frozen v1/v2 source or output was overwritten. This is a live inference run,
not a completed native collection or accuracy result.

Preparation reports **3,730 exact full-sequence matches** to the supplied CLEAN
`split100.csv` reference, distributed across all 18 genomes. The 1,022-residue
comparison finds the same 3,730 matching proteins and no additional ones. These
counts refer to the full original population, before any native predictions are
scored. This documented reference overlap prevents describing the forthcoming
comparison as an untouched sequence test; it is not a complete ESM pretraining
or homolog-leakage analysis. No proteins are removed on the basis of overlap.
