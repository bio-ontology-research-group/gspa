# Local reaction context without independent cohort peers

The reaction engine already computed each candidate's own neighborhood, but the
production bridge discarded those features whenever the target had no independent
peer genomes. The repair preserves local observations without turning missing
recurrence into negative evidence. It advances item 5; full context and biological
acceptance remain open.

## Behavior and contract

The existing `ExternalContext.Cohort` payload now has an optional `localOnly`
boolean, default false for older inputs. In this mode, the bridge retains local
anchors, matched local background, direct-EC mapping alternatives, family
membership and supplied intrinsic role evidence. All peer statistics remain null;
validation rejects peer measurements, even measured zero, in a local-only payload.
Missing coordinates still withhold all local-window features.

The extractor uses the existing contract-v3 columns: `cohort.w.available` remains
zero, while `cohort.w.local_available` independently records an observed local
window. Peer availability flags and values remain zero. No column was added or
redefined, and the contract fingerprint remains unchanged. Local/intrinsic
features retain their existing definitions. An unsupported availability profile
must not be interpreted as calibrated evidence by a model.

Without peers, an ambiguous EC mapping uses the reaction with the largest local
anchor weight, then reaction ID, as its representative. Other mappings remain
reported alternatives, with null peer support and effective N; evidence is not
summed across mappings. Runs with peers preserve the previous representative rule
and feature values. The report keeps `no_peers=true` and separately counts
`keys_with_local_features`, `keys_with_cohort_features` and `keys_with_context`.
The stage ledger includes the local-feature count.

Feature exports add the derived `available_contexts.local_reaction` boolean.
The Python training-data reader accepts older exports without it, validates the
new boolean against the feature vector, and still rejects unrelated unknown
fields. Existing export fixtures were updated only after verifying every other
field, feature and mechanism-model score was unchanged.

## Verification

[Archived evidence](../../benchmark/v3/results/local-reaction-context-v1/README.md)
includes all test logs, source snapshots, CLI outputs and hashes.

- The initial four focused regressions fail before the repair. The final eleven
  pass: singleton anchors, missing coordinates, curated interaction/role inputs,
  direction, JSON/merge preservation, and six forbidden peer-statistic cases.
  Existing unresolved-metadata and single/reference/batch parity tests now check
  preserved local evidence and absent recurrence.
- The complete JVM gate passes **1,172 tests**, zero failures/errors, two optional
  skips. All four live HMMER/AntiFam cases run. The Python benchmark gate passes
  **572 tests**, with twelve optional discovery skips; the nine live Nextflow
  cases also pass separately, leaving three optional checks unexecuted.
- Three added Python cases verify backward-compatible export reading, rejection
  of a mismatched local availability boolean, and continued rejection of unknown
  fields. Initial export-reader and fixture failures are retained in the archive.

`benchmark/v3/check_local_reaction_context.py` runs eight real production CLI
invocations against frozen inputs. Target genome GCA_000012525.1 contributes
**29,079 identical homology claims** in every arm, along with the same **50
candidate-gap records** from the verified gapsmith observation. The comparison
uses complete existing reaction/EC references, actual layouts and family
assignments, and subsets the existing metadata to the genomes in each arm.
No truth labels enter inference. Two acquired reference genomes exercise peer
compatibility; their inclusion is a software check, not a biologically optimized
cohort design.

| Arm | Candidates | Local window available | Nonzero local anchors | Old/new features |
|---|---:|---:|---:|---|
| Single genome | 29,079 | 1,182 | 20 | Local/intrinsic columns restored |
| Single genome, coordinates withheld | 29,079 | 0 | 0 | Intrinsic EC/family mapping columns restored |
| Same target with two independent peers | 29,079 | 1,182 | 15 | Every column identical |

Candidate identity, all sequence features, and all features outside the cohort
block are identical in every old/new comparison. Baseline functional prediction
TSVs are byte-identical in all three arms. Both current learned models abstain
on all 1,182 local-context candidates in the singleton replay, with null learned
probabilities. These runs also supply genuine gap indicators, so this abstention
does not isolate which availability flag is outside training coverage.

The final engine SHA-256 is
`0fde94f581cc5b871939681e1d11163662ef1ae6b9194c3bdba945f802850e68`;
the preceding engine is
`9c471e8ba6652d519e12bb99ad2b5b186e8e4f8c7cf75a62fd177e3a921370b3`.
An initial replay failed because the check supplied metadata for genomes outside
its selected cohort. The corrected check subsets metadata without changing any
retained values; that failed invocation is preserved.

## Remaining scope

This is evidence preservation through an explicitly requested cohort/context
stage. Automatic pathway-observation-to-gap wiring is still missing. Universal,
pathway-only and organism-specific network comparisons, calibrated scoring of
these local-only profiles, complete role/locus/null supervision, and independent
held-out biological validation remain required. The 50 gaps are hypotheses, not
truth labels. No whole requirement is accepted by this check.
