# Predictor role handoff into learned candidates

The pathway adapter preserved `candidate_role` in raw claims, but the shared
learned extractor ignored it. An actual gapsmith search on training genome
GCA_000012525.1 supplies 829 complex-subunit candidates; the preceding production
JAR reported all 829 as `UNSPECIFIED`. This is a reproducible item-7 defect,
separate from whether any proposed role is biologically correct.

## Repair

`ContextFeatureExtractor` now uses adapter-supplied roles before its conservative
function defaults. It validates the optional enum string, retains distinct roles
as separate candidates in deterministic enum order, and does not filter weak
scores. Malformed non-null values fail explicitly. Structured
`ExternalContext.roles` remains authoritative; candidate notes retain the
original supplied role and source, including when structured context overrides it.
Missing roles preserve the previous EC/GO/default behavior. Names or function
identifier prefixes are not interpreted as role evidence.

This uses the existing role features and candidate identity contract. It does
not add a probability, change sequence scores, or train a role classifier. The
observed adapter is the only production producer of this metadata field at the
time of this change.

## Verification

The [evidence archive](../../benchmark/v3/results/claim-role-bridge-v1/README.md)
retains before/after tests, source snapshots, complete CLI outputs and checksums.

- Twelve of the initial thirteen mechanism regressions fail before the repair.
  The final fourteen cases pass, including all six roles, weak-score retention,
  competing/duplicate roles, structured override, malformed values, conservative
  defaults, and raw JSONL round-trip feature/provenance parity.
- The old frozen JAR also fails the enhanced CLI oracle on the full observed
  reaction table: `observed role lost during learned candidate construction`.
- The repaired JAR passes seven production CLI runs, using the frozen genuine
  external output. The oracle independently derives roles from reaction-table
  complex status and EC availability. Both logistic and tree reports retain
  2,207 catalytic, 829 complex-subunit and 433 unspecified candidates.
- Against the preceding successful replay, precisely 829 role identities and
  their feature digests change. The other 2,640 feature digests are identical for
  each model. All raw claims and functional prediction TSVs are byte-identical.
  Both models abstain on all 3,469 unsupported-source candidates with null
  learned margins and probabilities.
- The complete JVM suite passes: **1,161 tests**, zero failures/errors, two
  optional skips. All four live HMMER/AntiFam cases run. Nine live Nextflow
  cases also pass using the established 26.04.6 launcher. An initial invocation
  accidentally used the older PATH launcher, which rejected the installed Java
  version; its failure log is preserved.

The tested engine SHA-256 is
`9c471e8ba6652d519e12bb99ad2b5b186e8e4f8c7cf75a62fd177e3a921370b3`.
The prior engine is
`bc9487caf1d2a4b1b787eff2eb40fbad70350f3dda9545aa87c719c99ef2b6f9`.
The actual external search and its provenance are documented in
[the gapsmith observation](v3-gapsmith-observed.md); these CLI runs replay its
outputs and are not new searches.

## Remaining acceptance

This closes the raw-predictor-to-candidate role handoff, not item 7. Independent
role labels, calibrated role/locus/null models, difficult same-operon alternatives
and held-out biological validation remain necessary. Automatic pathway-gap
wiring, matched reaction-context ablations and preservation of local context
when no cohort peers exist also remain open. No whole audit item is signed off.
