# Community result ledger and summary validation

Status: repaired and verified through the production workflow. This advances
items 2 and 9; it does not validate biological cross-feeding or connect rescue to
a gene assignment.

The prior `CommunityStage.verify` accepted incomplete results when their model
hashes matched. Three direct reproductions through the immutable production JAR
accepted a missing solver ledger, 41 used calls against a budget of 40, and a
missing organism summary. Twenty-three new regression cases failed against the
old verifier. Four malformed real-CLI sidecar outputs were also reported as
executed by the old engine.

## Repair

The stage now requires the exact model set, unique organism identifiers and
complete organism/conditional-evidence coverage. Solver budget and usage must be
integers, usage must stay within budget, and every indexed call needs a valid
status and finite optimal objective. Experiment and FVA references must account
for every solver call exactly once. Growth records must agree with their solver
objectives, threshold, status and matched monoculture/community identities.
Summaries are checked against those experiments, and conditional records must
retain the same organism classification, model-conditional status, explicit null
probability and no-enzyme-assignment disposition.

Missing data cannot silently become an empty successful report. Budget exhaustion
remains a valid explicit unknown; it is not turned into a biological negative.
This patch changes the Java/Groovy stage validator only. The Python solver,
community model and frozen remote numerical experiment remain unchanged.

## Evidence

The old engine is
`74d9e8fc8c170adcb8d222cc0fb2ba82199b877e63a51efffca40f8dac512405`;
the repaired engine is
`b23a7f6d94eb06c375e63ad344cfbf4c15d7644e33c1848766a36c4287cfd090`.

- **76 focused tests pass**: 29 contract cases and 47 real-sidecar workflow cases.
  They retain valid complete and budget-limited results while rejecting missing,
  malformed, duplicated, contradictory and over-budget records.
- **1,356 JVM tests pass**, with zero failures/errors and six optional skips, after
  forced recompilation and execution of all test tasks.
- **14 production CLI runs** compare both engines on a baseline, a valid community,
  a one-call budget, and four malformed reports. The old engine accepts every
  malformed report. The new engine rejects all four, records a failed stage, and
  publishes no community result.
- Valid functional TSVs, complete candidate-feature files and raw community
  results are byte-identical between engines. Valid community and budget-limited
  runs preserve baseline functional TSVs; no candidate feature is invented.

The first CLI checker incorrectly expected the disabled stage status to be
`not_requested`; the actual contract is `skipped` with reason `disabled`. That
failure and its completed baseline run are retained. The corrected checker reran
all fourteen cases in a fresh directory. The retry regenerated one unused
first-attempt sidecar fixture before archival; its original bytes are unavailable,
and only its expected hash remains. The first attempt had stopped at the baseline
before using that fixture. All final-run inputs and outputs are retained and
verified; the historical omission is recorded explicitly.

Working evidence is in `~/tmp/gspa-community-contract-v1/`, including old/new
source, reproductions, failure/pass logs, complete test XML, registered commands,
input/output hashes and `cli-verification-v2.json`.

## Limits

These are synthetic SBML/LP and result-contract checks. They do not establish
exchange or transporter truth, biological accuracy, or reduced false enzyme
attribution. The community stage still exports organism-specific conditional
explanations with no candidate feature, probability or gene promotion. Complete
request/assumption provenance and independent biological interpretation remain
separate acceptance requirements; these checks are not an independent solver.
The six-shard empirical community study and its frozen analysis remain separate.

The [immutable repair archive](../../benchmark/v3/results/community-result-contract-v1/README.md)
preserves the complete final proof and explicitly records the unused historical
fixture omission. All 731 archive checksums pass.
