# Gene-deletion identity across metabolic workflows

Explicit gene deletions could silently do nothing when the supplied identifier
was an SBML gene label rather than its gene-product ID. Unknown genes were also
treated as harmless deletions. The community route separately accepted unknown
genes and overwrote earlier deletion records for the same organism. These faults
can make a required activity appear dispensable or leave an alternative isozyme
active during a supposed combined knockout.

The shared `resolve_knockouts` function now validates gene/reaction identifiers,
resolves declared labels and evaluates all deleted genes together through the
GPR. Both single-organism and community experiments use it. The community route
validates requested deletions before baseline solves and rejects an unknown or
excluded organism. The SBML loader rejects ambiguous labels, including a label
that collides with another gene's ID. This deliberately rejects ambiguous models
instead of choosing an arbitrary gene. It does not guess identifiers by removing
prefixes; a locus tag absent from the declared ID/label map remains unknown.

All explicit knockout identities are validated before bounds are changed.
This is not a general transaction guarantee for every other malformed medium
modification. An unknown candidate gene retained alongside an explicit reaction
hypothesis remains provenance, not an executed gene deletion.

## Validation

Final Slurm job **10218**, node002, passes **144** experiment/community tests and
**118** affected benchmark tests: **262 tests, no failures or skips**. The
14-method focused suite produces 15 failure events and five errors against the
original files, including subtests, and passes after repair. It covers label/ID
equivalence, unknown identifiers, collisions, mixed direct-reaction/gene
deletions, AND subunits, OR isozymes, the masked-assignment loop and community
experiments with repeated organism records. Unknown community deletions fail
before a solver call.

The [Sco-GEM diagnostic](v3-germicidin-diagnostic.md) adds a real network test.
For all 18 native-ID scenarios, before/after predictions and deleted reactions
are exactly equal, and independent COBRApy/GLPK agrees with the production
HiGHS results. A derived fixture changes only `G_SCO7221`'s declared SBML label
to `SCO7221`; stoichiometry, bounds, GPRs, reaction order and metabolite order
remain identical. Before repair, deleting that label leaves all four germicidin
reactions active and predicts positive synthesis capacity. After repair, it
disables the same four reactions as the native ID, and the two product-capacity
checks exactly reproduce the native-ID deletion. The fixture is an identifier
representation test, not an additional biological strain.

Claude Sonnet 5.5 reviewed the identity patch and its callers. Its review exposed
the community overwrite and ambiguous-label faults, which were then repaired.
It also prompted positive biological controls, explicit provenance for curated
experimental metadata and relabelling the repeated native solve as a state
control rather than complementation.

## Retained execution history

- **10215:** initial identity repair passes 103 experiment and 67 benchmark
  tests; the original eight-method regression run fails as expected.
- **10216:** first real-model run passes all 18 independent scenario comparisons
  and 90 LP checks. Native repetition was initially named `gcs_restored`; the
  final protocol correctly calls it `native_repeat`, a state-control replay.
- **10217:** expanded suite passes 143 of 144 tests. The failing assertion
  expected zero common growth after deleting vitamin synthesis and overflow;
  removing both disposal routes actually makes ATP maintenance infeasible.
  The test was corrected to preserve that physiological distinction. No solver
  setting or model bound was changed to obtain a passing result.
- **10218:** complete final tests, native-ID parity, independent solver checks,
  positive/negative controls and two real-model label regressions pass.

Archive: `benchmark/v3/results/knockout-identity-v1/`; its verifier checks input
and source hashes, all archived bytes, test totals, native-ID parity, solver
agreement, label-fixture invariants and the complete primary observation table.
This repair improves execution correctness. It does not establish biological
cross-feeding attribution or function recovery; all twelve whole items remain open.
