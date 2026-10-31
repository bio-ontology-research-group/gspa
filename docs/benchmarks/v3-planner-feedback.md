# Standard metabolic planner: feedback repair

The standard agent loop previously discarded the measurements needed to reason
about secondary metabolism. Equal growth produced identical prompts whether a
product was feasible, impossible, unassessed or inconclusive. FVA and reaction
knockout checks were also discarded, and prior evidence appeared only as IDs.
Four focused tests reproduce these omissions before the repair (two failures,
two missing-interface/field errors) and pass afterward.

`scripts/metabolic_experiments.py` now projects hypothesis/evidence records,
medium and objective assumptions, menu interventions, growth, FVA, reaction
knockouts and product minimum/maximum flux into a versioned prompt. Missing
checks and solver failures retain their distinct statuses. Arbitrary metadata,
configuration secrets and model matrices are not copied. This is a field
projection, not anonymization: identifiers and allowed text may still contain
prior knowledge. Benchmark builders remain responsible for preventing leakage.
Evidence types are not learned or calibrated probabilities.

Each Qwen selection records the prompt schema and content hash. Prompt-building
errors propagate instead of becoming misleading endpoint failures. The existing
Keio benchmark wrapper accepts the additional engine arguments while preserving
its original prompt; its logs say `custom_unversioned`. Historical benchmark
results have not been recomputed or reinterpreted as validation of this prompt.

## Validation

All model and solver execution ran on unimatrix node002 through Slurm `debug`.
Job **10213** passes **95** experiment tests and **67** affected benchmark tests.
These include product/FVA feedback in the real production loop, explicit checks
against swallowed callback assertions, local prompt failures, and the separate
Keio wrapper. This is a Python-only repair; no JVM production code changed.

Job **10214** compares fixed, random and evidence selection before/after at two
seeds: all six complete output pairs are exactly equal. It also runs both agent
versions against the unimatrix Qwen endpoint using a synthetic SBML network
with an optional secondary product and a biomass-coupled coproduct. Both use the
same menu, four-experiment limit and 60-call solver cap. This is a mechanism
diagnostic, not an established biological benchmark.

| Live diagnostic | Before | After |
|---|---:|---:|
| Successful Qwen requests | 4 | 4 |
| Invalid choices / fallback | 0 / no | 0 / no |
| Solver calls | 21 | 21 |
| Prompt tokens, total | 1,345 | 6,299 |
| Completion tokens, total | 242 | 263 |
| Endpoint seconds, total | 4.568 | 6.604 |

Before the repair, Qwen chose both knockouts and then both product checks. After
the repair it chose the product checks first. All final experiment measurements
are identical, and all four final prompt hashes reproduce from the retained
transcripts. This establishes delivery and use of a richer observation surface;
it does **not** establish more correct functions, greater information gain or
better performance at a restricted budget. Prompt cost increased substantially.
The endpoint model alias was `qwen3.8-27b`; the server's weight revision was not
independently pinned, so these transcripts do not establish future repeatability.

The actual rationales also expose unresolved biological reasoning errors. Qwen
describes feasible flux as upregulation, treats optional/dispensable growth flux
as refuting the secondary-product hypothesis, and describes a growth-only
knockout as testing product secretion. The last interpretation exceeds what that
experiment measured. These rationales remain verbatim in the archive. No result
is promoted to a function assignment or calibrated confidence; item 12 remains
open. Correctly transferring the observations does not make the agent's
interpretation reliable.

## Retained failures and review

- **10211:** missing test import path prevented two test modules from loading;
  the 58 loaded checks passed. The launcher now supplies the test module path.
- **10212:** 94 experiment and 41 benchmark tests pass; the paired diagnostic
  then fails because a new-module model raises a different exception class from
  the old-module handler. Each version now loads its own model instance.
- **10213:** all final tests and six non-agent comparisons pass; the live call
  fails twice with `URLError` and the engine falls back explicitly. This failed
  run is retained, not counted as agent success.
- **10214:** the same frozen source succeeds using the verified unimatrix IPv4
  endpoint, avoiding the relocated container's unavailable DNS configuration.

Claude Sonnet 5.5 reviewed the narrow change without running models. Its findings
led to the benchmark-wrapper fix, assertions outside the callback, explicit
request-count/validity checks, and separating local prompt errors from endpoint
errors. Optional three-argument prompt construction remains for diagnostic
callers; production always supplies context and resolved configuration. Prompt
size control and empirical planner benefit remain unresolved.

Archive: `benchmark/v3/results/planner-feedback-v1/`. Run its
`verify_archive.py` to check the complete byte inventory, source package binding,
test counts, six output pairs, live transcripts, prompt hashes and per-experiment
measurement parity. Sources and failed attempts are retained. All twelve whole
acceptance items remain open.
