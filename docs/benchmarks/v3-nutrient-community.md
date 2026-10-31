# Matched nutrient-dependent coculture evaluation

This extends the [audited nutrient dataset](v3-nutrient-environments.md) to all
1,320 two-strain cocultures: 1,200 primary assays and 120 no-added-nitrogen
controls. The six existing medium/oxygen arms yield 7,920 case/arm plans. No
case is selected using OD, coverage, agreement or predicted rescue.

**Current status:** all production, numerical and frozen endpoint analyses have
completed. The [complete results](v3-nutrient-community-results.md) retain one
original solver unknown and mixed/negative biological outcomes. Cross-feeding
attribution and whole-item acceptance remain unverified.

## Model and assay correspondence

The models and canonical gene mappings are the same registered iML1515 and
iYO844 inputs used in the [single-organism evaluation](v3-nutrient-phenotypes.md).
Each pair has two separate stoichiometric blocks and two separately evaluated
deletion genotypes. The pair is never encoded as one double mutant. All affected
GPR reactions and isozyme bypasses are preserved, along with the model-native
ATP maintenance and internal transport bounds.

The printed and additional-inorganic arms each use oxygen ceilings of 0, 2 and
20. Ingredient bindings and fixed uptake assumptions are reused unchanged from
the single-organism protocol, including unrepresented ingredients and the
explicitly unprinted inorganic sensitivity. Every missing exchange is recorded.

Ports map exact extracellular model species to a shared pool. Pairs contain
strains of the same species/model; intracellular species are never merged.
All model exchange ports are marked **predicted**, not experimentally supported
transport. Environment supply is finite and biomass-weighted; removal has an
explicit ceiling of 1,000 mmol/gDW-community/h for each pool species. Opening a
boundary exchange does not add an internal membrane-transport reaction.

The common-growth formulation fixes biomass fractions at 1:1 and requires both
members to have the same growth rate. Equal inoculum OD is only a motivation for
that assumption: the assay did not measure equal biomass or steady balanced
growth. Its endpoint is bulk OD after 96 hours. A positive endpoint does not
establish that both members grew, identify exchanged metabolites or validate
transporters.

## Frozen comparisons and budgets

Each case/arm receives exactly five production LP calls, in this order:

1. The first deletion strain alone.
2. The second deletion strain alone.
3. The balanced community with permissive exchange uptake: every exchange can
   import up to at least 1,000, constrained by pool supply and partner secretion.
4. The community with native uptake bounds retained for compounds with zero
   environmental supply; supplied compounds retain their opened ports.
5. The community with uptake of every unsupplied compound closed. Sharing
   supplied compounds remains possible, so this is not an absolute prohibition
   of all organism interaction.

The total is 39,600 production LPs. All three community policies are reported;
none is selected using endpoint measurements. Failure and unrun outcomes remain
unknown. Maintenance infeasibility retains its status and predicts no feasible
growth under the stated constraints. No annotation or calibrated probability is
produced. Transfer fluxes from a single optimal solution are conditional
candidates; they are not evidence of unique or observed transfer.

## Implementation and independent checks

`run_nutrient_community.py prepare` reads outcome-free assay inputs and model
registration. It freezes common pool/port definitions and all execution plans,
binding inputs, models and implementation by SHA-256. The production runner
calls `community_metabolism.Analysis.growth` and its LP builder with cached
parsed models. It uses the same solver, organism validation, pool constraints,
budget accounting and transfer records as the production community path.

The six shards are deterministic interleaved partitions of the frozen plan.
`merge` verifies receipts, hashes, exact plan identity/order and complete
coverage before producing a combined run receipt. It rejects missing, duplicate,
modified or mismatched shard outputs. Simulation packages contain no observation
file. Only the separate `evaluate` command joins the reconstructed endpoints.

The independent `check_nutrient_community_cobra.py` uses the existing COBRApy
shared-pool formulation, reparses SBML/GPRs independently and solves with GLPK.
It does not use the production matrix, GPR evaluator or solver. It verifies
every member's disabled-reaction set and every one of the five LP outcomes.
The GLPK basis is reset before each solve, with primal feasibility tolerance
1e-9. Objective comparison uses 1e-6 absolute and relative tolerance; finite
primal values, mass balance and reaction bounds are also checked at 1e-6.
Unknown production solves and numerical disagreement remain failed comparisons.

The evaluator reports bulk-growth agreement and a separate paired-rescue proxy:
bulk coculture OD > 0.02 while both matched monoculture endpoints are at or below
0.02. Missing endpoints remain unknown. Always-positive and always-negative
comparators are retained for both endpoints. These endpoint comparisons do not
become gene-function or transporter truth.

## Verified checkpoint

Seven focused tests pass. They exercise reciprocal rescue, both exchange signs,
separate mutant genotypes, pool/uptake controls, exact five-call budgets,
unchanged parsed-model bounds, independent GLPK agreement, unknown outcomes,
outcome-free preparation/execution, shard/serial parity and rejection of
incomplete or tampered data. The complete benchmark gate passes 645 tests with
13 optional skips. No production JVM or Python solver source changed.

A deterministic preliminary replay uses the first input case of each species
under all six arms: 12 real-model cases and 60 LP comparisons. All agree with
the independent checker. This is preliminary numerical evidence, not proof for
the full catalogue.

The final plan is `nutrient-community-v1/plan-v3` in the local data cache; all
7,920 plan rows have SHA-256
`eaa665f93ea9a2c9b4bbd943b28b441711222c8cc77d5c770ebea59c4a438b46`.
Earlier preparation snapshots have identical plan rows but precede the final
sharding/resource-management implementation; their implementation bindings must
not be used for execution. The immutable remote package is
`/data/hohndor/gspa-v3/nutrient-community-v1`. Its files passed checksum
verification before Slurm submission. Array **10150**, tasks 0–5, runs production
and independent replay with two tasks at a time. Initially each has a two-hour limit,
12 GB RAM, local scratch, persisted logs/results and an early termination handler.
The checkpoint observed 243 and 234 written rows in the first two shards. These
are progress counts, not verified completed results. Current production log
messages label the original interleaved plan position as `completed`; use actual
row counts and the final run receipts for per-shard counts. Some shard indices
do not emit those periodic messages. Direct worker SSH was unavailable;
read-only Slurm steps within the existing allocation supplied the progress check.

The checkpoint package, tests, model inputs, frozen plans, source snapshots and
submission/status receipts are retained under
`benchmark/v3/results/nutrient-community-checkpoint-v1/`. It deliberately contains
no claimed full-run evaluation.

The separate [paired analysis specification](v3-nutrient-community-analysis.md)
is now frozen and tested. It retains unknowns and adds environment/strain
resampling and leave-one-strain-out sensitivity without changing this execution
package. Full numerical verification and matched endpoint results remain pending.
Measured exchanged-metabolite and transporter-knockout data,
independent gene attribution and the other original acceptance gates also remain
open. Earlier negative single-organism results are unchanged.

### Runtime allowance update (2026-10-27)

Independent oracle checks ran slower than production. An attempt to extend array
10150 to three hours succeeded for queued tasks 2–5; Slurm denied increasing the
limit of already-running tasks 0–1, which retain two hours. The attempt returned
exit code 1 with a partial update, so it is not recorded as a complete extension.
The frozen package, plans and per-case solver budgets are unchanged. Terminal and
partial receipts must be inspected before any continuation. The separate cache
record `nutrient-community-v1/runtime-budget-extension.json` retains the request,
response, status snapshot and frozen submission-script hash.

The [first terminal-shard audit and numerical diagnosis](v3-nutrient-community-audit.md)
records the completed first shard's solver unknown, the second shard's timed-out
verification and its separate oracle-only replay. It also validates growth
classifications independently of scalar tolerances. Full collection verification
and the frozen endpoint analysis are still pending.
