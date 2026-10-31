# Observed pathway gaps in shared integration

Annotation now passes the current pathway search into the shared integration
workflow. The adapter retains both the reaction table and the separate pathway
prediction table, with their hashes. Before this change, the latter table was
discarded and no automatic gaps reached integration.

## Input and eligibility

The workflow accepts current in-memory observations from annotation, or an
exported observation JSON through `integrate --pathway-observations FILE` /
`pathwayObservationsFile`. The published configuration-key list includes the
new setting. Standalone calls to `AnnotationPipeline.integrate` do not reuse an
earlier annotation run's observations. An observation for another genome fails
preparation with a failed ledger.

Automatic derivation requires `target=proteome`, a `SUCCESS` or `SUCCESS_ZERO`
outcome, and both observed tables. Completion means every supplied protein was
processed; biological assembly completeness is separate. Partial, failed,
unavailable, nucleotide and representative-only searches do not supply gap
evidence. Missing pathway tables remain explicitly unavailable. A malformed
observed table fails the adapter; it cannot certify a gap search.

The builder groups every row by pathway and reaction before interpreting missing
hits. A group is eligible only when its pathway was explicitly predicted present,
all hits are weak or absent, and no row indicates a strong hit, missing reference
sequences, spontaneity or a complex. Complex subunit candidates remain in raw
sequence evidence; missing subunits do not become a whole-reaction EC gap.
Any weak hit marks the eligible group as guessed, independent of row order.

EC alternatives and all supplied EC-to-GO alternatives are retained. Partial and
preliminary EC identifiers use the same parser as sequence candidates; local
reference tokens remain in raw observations and are not asserted as standard
ECs. Missing mappings remain null. `--ec2go` supplies the optional mapping file
and its hash is recorded. Pathway IDs are source-qualified opaque identifiers
such as `gapsmith:CHLOROPHYLL-SYN`; the workflow does not assign every identifier
to MetaCyc, as the historical benchmark converter did.

Explicit `--gaps` takes precedence, including an explicitly empty file. That
choice and the number of superseded searches are reported. Candidate gaps remain
hypotheses, not evidence of biological absence or proof of enzyme function.

## Reporting and reuse

Successful integration writes `<out>.gap-context.json`, including availability,
source hashes, eligibility/exclusion counts and reasons for unavailable searches.
When gap context is available, `<out>.gaps.jsonl` contains the normalized records
for replay through the existing single, batch or Nextflow gap-input path.
An observed zero produces an empty file; unavailable input produces no gap file.
Reference preparation records the gap count and availability in the stage ledger.

`ctx.gap.available` now distinguishes a supplied/observed empty set from missing
input. The existing feature columns and contract fingerprint are unchanged.
Outdated artifacts are removed on repeated runs, while explicitly reusing the
prior gap export with the same output prefix preserves that input. A regression
reproduced its deletion before the cleanup guard was added.

Benchmark receipts now hash `--pathway-observations` in both option spellings.
The new regression failed before this addition: changing observations otherwise
reused stale results. Existing workflow-file dependencies still need the runner's
declared dependency mechanism; this change does not claim general recursive
discovery of every path inside arbitrary JSON.

## Verification

The [evidence archive](../../benchmark/v3/results/pathway-gap-bridge-v1/README.md)
contains before/after failures, complete output, source snapshots and checksums.

- The pipeline regression fails before wiring and passes after it. Thirty-three
  added core cases cover eligibility, alternatives, missingness, file/in-memory
  parity, explicit precedence, stale artifacts and genome binding. Six added
  predictor cases cover automatic handoff and malformed pathway tables; the
  existing repeated-failure case also verifies both tables are cleared.
- The full JVM gate passes **1,211 tests**, zero failures/errors, two optional
  skips; all four live HMMER/AntiFam cases run. The benchmark gate passes **573
  tests**, with twelve optional discovery skips. All nine live Nextflow checks
  pass separately, leaving three optional checks unexecuted.
- The preceding JAR fails the new real-output oracle: `raw pathway table lost`.
  The repaired JAR passes **ten actual CLI runs** replaying the verified gapsmith
  output for the acquired training proteome GCA_000012525.1. This is a replay of
  the earlier external search, not ten new searches.

All **13,147 reaction rows**, **1,923 pathway rows** and **3,469 raw candidate
claims** survive. Automatic integration yields exactly the **50 candidate-gap
records across 46 reactions** independently audited earlier, with the deliberate
pathway namespace correction described above; 48 records have a GO mapping.
The oracle is bound to the same reaction table, pathway table and EC mapping.

Annotation, direct observation-file integration and explicit-gap replay have
identical feature vectors and functional prediction bytes. Relative to the old
engine, only three feature columns change: gap availability for 3,469 candidates,
gap listing for 53, and guessed-gap listing for the same 53. Raw claims and
baseline prediction TSVs remain byte-identical. Both current learned models
abstain on all candidates from this unsupported source. Injected zero, missing
pathway table, malformed reaction table and malformed pathway table runs retain
distinct output/status behavior.

The final engine SHA-256 is
`fb7f750e451799a93446d6931a29817ec29adcd74ed00cf24828238fcdf3ee3b`;
the preceding engine is
`0fde94f581cc5b871939681e1d11163662ef1ae6b9194c3bdba945f802850e68`.

## Remaining acceptance

This completes the observed-gap handoff for the verified proteome-search path.
It does not validate gap hypotheses biologically, calibrate their contribution,
or establish the full context ablations. The legacy RLC suggester also needs a
pathway database supplying its EC-to-GO anchor lookup; merely deriving mapped
gap records does not complete that independent input path. Pathway-only,
universal-network and organism-specific comparisons, calibrated role/locus/null
inference and held-out biological validation remain open. No whole audit item is
signed off.
