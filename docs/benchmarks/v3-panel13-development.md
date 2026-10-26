# Production GO benchmark: 13-genome development panel

All 39 arm/genome cells completed through the production CLI built at `3c75a38`
(Slurm job 9739, 21 min 22 s). The [full report](../../benchmark/v3/results/development-3c75a38/panel13-go-report.json)
records the exact commands, input/output hashes, per-genome results, failures and
paired genome-bootstrap intervals. Its binary SHA-256 is
`65fe3068a17e3b1c50a157335097432633dba6dfff7e9b2dd80abc37a06c028f`.
The remote copy lacks `.git`, so the report's source revision is null; the
commit-to-binary association is recorded here and in the launch record, not
inferred from a non-null runner field.

These are previously inspected development genomes, not an untouched test set.
Exact duplicate claims were removed by the registered preparation procedure.
The arms use identical remaining claims: evidence integration with no priors,
the default heuristic priors, and those priors plus symbolic constraint relaxation.
They do not evaluate the new trained context model, cohort recurrence, or adaptive
FBA stages. Historical auxiliary gap inputs include synthetic fallback gaps;
the real GO labels do not make those auxiliary inputs experimentally validated.

Experimental GO labels cover 12 of the 13 genomes: `rprowazekii` executed but has
no evaluable labels. There are 323 labelled proteins and 3,045 propagated true
protein–GO pairs. GO ancestors are propagated in every arm. Unlabelled predictions
on labelled proteins count as false under the declared incomplete-label convention;
this can penalize correct novel annotations. The comparison universe contains
14,225 pairs and is shared by all arms.

| Production arm | Fmax | F1 at 0.5 | ECE | Corrected / harmed pairs at 0.5 |
| --- | ---: | ---: | ---: | ---: |
| No priors | 0.34638 | 0.34157 | 0.36034 | Baseline |
| Default priors | 0.35021 | 0.34678 | 0.34353 | 78 / 22 |
| Priors + constraint relaxation | 0.35058 | 0.34472 | 0.33470 | 106 / 154 |

Fmax selects its threshold on these evaluation labels and is descriptive, not a
deployable threshold selected independently. At the fixed threshold, the 95%
paired genome-bootstrap interval for the F1 change is −0.00008 to +0.01667 for
priors, and −0.00356 to +0.02000 for priors plus constraints (2,000 replicates).
The constraint arm has slightly higher Fmax but more harmed than corrected
binary decisions. The net correction interval is −198 to +111 pairs. The
pair-level McNemar p-values in the report ignore dependencies among GO ancestors
and proteins; genome-level intervals are the relevant uncertainty here.

This run supports a small development-set change, not a large accuracy gain.
It is not an EC prediction comparison with CLEAN. Independent splits, trained
context arms, more complete labels and the remaining ablations are still required.

The earlier job 9732 timed out while concurrent genomes repeated ontology
extraction. Its outputs remain preserved separately. The successful run uses
single-flight exact inferred-relation extraction and a validated persistent cache;
no asserted-only approximation was substituted to obtain these results.
