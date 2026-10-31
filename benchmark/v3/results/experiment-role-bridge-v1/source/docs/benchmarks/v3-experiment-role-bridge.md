# Raw candidate roles in metabolic experiments

The masked-assignment builder previously read candidate roles only from external
context. A raw `candidate_role: COMPLEX_SUBUNIT` observation could therefore
become a target-alone enzyme assignment. Accessory, regulator and transporter
claims had the same problem. Invalid role spellings silently acquired the EC
default. This violated the existing candidate-role contract and the documented
partner requirement.

The repair uses the shared role resolver with only the retained raw evidence
channels that can contribute to experiment construction. Derived learned,
metabolic and inferred claims cannot override those roles. Explicit structured
context retains its established precedence. Absent role metadata retains the
existing function default. Each proposed hypothesis records raw role observations
and their sources, the resolved roles and the source of the resolution in the
private construction report; they do not enter the planner's public problem.

Complex-subunit evidence requires independently supplied partner clauses, even
when another observation also calls the candidate catalytic. The supported
assignment is an AND clause such as `MASKED_TARGET AND gB`, with an explicit
partner-knockout condition. Noncatalytic-only candidates cannot become enzyme
assignments merely because a support record names partners. Invalid raw role
metadata fails validation instead of falling back to a catalytic default.

The same review found that explicit and guided absence comparisons allowed a
catalytic-plus-subunit candidate to bypass their standalone-activity guard.
Those formats have no supported partner clause. Explicit candidate binding now
rejects that combination; guided construction retains the reaction comparison
with candidate assignment unavailable. This does not remove complex hypotheses
from the richer masked-assignment format.

## Verification

Thirteen regression cases were added. Twelve fail before repair; the explicit
context-precedence case already passes and guards the retained behavior. The
failing cases cover raw noncatalytic roles, raw subunits with and without partner
support, mixed catalytic/subunit observations, derived-role interference,
malformed values, and both absence routes. All 105 focused cases pass after
repair, including the existing learned role bridge and real-sidecar workflow
tests.

`benchmark/v3/run_experiment_role_bridge.py` compares frozen old and new JARs
through the actual CLI. It supplies identical synthetic model, raw claims,
conditions and support records to both versions. It checks that unsupported
enzyme assignments disappear, valid partner assignments remain AND clauses,
the new role provenance survives, and ordinary completed runs preserve annotation
bytes. The unchanged default and supported cases must also retain identical
public problems and phenotype cells. Every executed phenotype is replayed with
the independent COBRApy/GLPK checker.

These checks establish the role handoff and conditional simulation behavior.
The claims, partner support and network are synthetic. They do not establish
biological role accuracy, correct partner identity, calibrated assignment
recovery or whole-item acceptance. Independent role-labelled loci and empirical
validation remain required for items 7 and 12.

The complete JVM gate passes 1,278 tests with zero failures and six existing
optional skips. The final CLI replay completes nine matched old/new cases
(18 invocations), with 90 executed phenotype comparisons agreeing with GLPK.
Five masked cases change from execution to no admissible assignment; the mixed
absence case changes from execution to validation failure. All eight pairs that
complete successfully preserve annotation bytes. The three cases that continue
to execute preserve phenotype cells, and both masked cases preserve the public
problem exactly. Failure publishes no experiment result.

The initial CLI harness omitted the required `completeness_fraction: null` from
its unknown-assembly fixture and stopped before the absence cases. That failure
and the initial harness are retained. The final replay supplies the required
field and runs every case from a fresh directory.

Evidence is archived under `benchmark/v3/results/experiment-role-bridge-v1/`.
The old engine is
`2b07dfbdb166a016f244caccc2324f82be211527c1b740b3ebe9ec146e9c0c10`;
the repaired engine is
`256121146db83406f5e02bb20920f793d722da977404f37aaad2f1b932ef5cf9`.
Both are retained in the local content-addressed engine cache. The runner takes
`--old-jar`, `--new-jar` and a new `--out` directory, using its Python environment
for the sidecars and independent checkers. Exact CLI arguments, source snapshots,
test XML, models, claims, assumptions, outputs and checksums accompany the report.
