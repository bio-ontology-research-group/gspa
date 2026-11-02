# Pre-registered plan: scoring the four absence explanations against known truth, and the effect of false enzyme assignments on coculture predictions

Date: 2026-11-01. Author: Fable agent for Robert Hoehndorf. Frozen before any scored check is run. The
SHA-256 of this file (body above the "Deviations" heading) is recorded in
`benchmark/v3/results/absence-crossfeeding-closure-v1/plan/` and in the report
`docs/benchmarks/v3-absence-crossfeeding-closure.md`. Nothing above that heading is changed after the
hash is taken; anything that has to change is appended there, dated.

This plan addresses the MISSING elements of acceptance rows 8 and 9 in `tasks/v3-acceptance.md`, as
itemised by `tasks/v3-acceptance-audit-2026-11-01.md` (rows 8 and 9). It does not change either bar.

Bars, quoted:

* Row 8: "Gene attribution scored separately from growth/KO accuracy on complete genomes, auxotrophs and
  controlled sequence loss."
* Row 9: "Measured mono/coculture, exchanged metabolites and transporter KOs; separate false enzyme
  assignment effects."

Already shown elsewhere and not repeated here: growth/KO accuracy on auxotrophs
(`docs/benchmarks/v3-nutrient-phenotypes.md`), measured mono/coculture agreement
(`docs/benchmarks/v3-nutrient-community-results.md`), gene attribution on complete genomes
(`docs/benchmarks/v3-offgrowth-discovery.md`), masked recovery on auxotrophs
(`docs/benchmarks/v3-wp5-masked-recovery.md`).

## Reconnaissance already done before freezing (environment, not scored results)

* The production absence engine (`scripts/absence_experiments.py`, schema `gspa.absence_experiments/1`)
  runs on the real `iYO844` model with a real frozen nutrient manifest: 4 hypotheses x 1 condition
  execute, and `external_supply` requires the supplement exchange to be declared in
  `baseline.exchanges` and its transport ports to be real model reactions
  (`absence_experiments.transport_reactions`). Declaring an exchange does not supply it:
  `metabolic_experiments.baseline_bounds` closes uptake of every external exchange that is not in
  `baseline.medium`. One coculture plan costs 2.6-5.2 s of one core through the production community
  engine (`scripts/community_metabolism.py` via `benchmark/v3/run_nutrient_community.simulate`).
* All inputs used below were extracted from the committed archives and verified by SHA-256 into
  `/mnt/data2/leechuck/gspa-work/absence-crossfeeding/inputs/`: `iML1515.xml`
  (`9c772d44...`), `iYO844.xml` (`28e4ca32...`), the frozen nutrient plan
  (`plans.jsonl` `4cedbf05...`, `protocol.json` `feef0885...`), the frozen coculture plan
  (`plans.jsonl` `eaa665f9...`, `contexts.json` `2c5fcb13...`, `protocol.json` `cf5bbda3...`), the
  merged baseline coculture predictions (`d4c481eb...` per shard; merged `analysis-output:merged`),
  and the measured observations (`886326e8...`). The `$HOME` cache tree is being migrated by another
  process today and is deliberately not an input.

## Part A (row 8): is the known explanation in the feasible set?

### A0. Scoring definitions, fixed here

One *comparison* is one production `absence_comparison` run: four hypotheses
(`missing_annotation`, `missing_sequence`, `pathway_loss`, `external_supply`) x the perturbations the
engine reserves, over explicit conditions. For a case with measured conditions:

* A hypothesis is **available** when the engine executes it (`status: completed`); `missing_annotation`
  is unavailable without a retained candidate claim and `external_supply` without a declared supply and
  transport, by the engine's own contract. Unavailable is reported, never scored as a negative.
* The **feasible set** F is the set of available hypotheses whose evaluated growth class equals the
  measured growth class in every scored condition. `infeasible_under_constraints` is mapped to
  `no_growth` (the engine's own rule); `unknown` cells make that condition not comparable and are
  reported.
* **Rank**: hypotheses are ordered by the construction's declared queue order, which is evidence
  support and then stable identity. The rank of the true explanation is reported; no new ordering is
  invented for this benchmark.
* Per explanation class t: **recall** = (cases with truth t whose F contains t) / (cases with truth t);
  **precision** = (cases with truth t whose F contains t) / (cases of any truth whose F contains t);
  **unique identification** = F equals exactly {t}. Intervals: Wilson 95% for binomials, plus a
  clustered bootstrap over (species, strain, environment) for the nutrient family, which is not a
  population estimate and is labelled descriptive.
* Growth agreement is **not** part of these numbers. It stays in the existing reports.

### A1. Nutrient auxotroph deletions (truth: pathway loss / external supply) - local, ws

Inputs: the frozen nutrient plan and measured observations above; models `iML1515` (E. coli EcNR1) and
`iYO844` (B. subtilis 168). Case construction reads no measurement.

Per primary assay of a deletion strain, the focal reaction is the set of reactions the deleted
canonical gene blocks in that model, taken from the frozen registration. Scope rules, fixed:

* exactly one blocked reaction: **primary scope**;
* more than one blocked reaction (E. coli `serA`): **sensitivity scope**, one comparison per blocked
  reaction, reported separately and never pooled into the primary numbers;
* no blocked reaction (E. coli `leuB`, isozyme bypass): **not constructible**, counted and reported,
  never dropped silently.

Two conditions per comparison, both measured: `unsupplemented` (the arm's medium without the
supplement) and `supplemented` (the same medium plus the measured supplement at the frozen rate of the
matched assay). `external_supply` declares that supplement exchange, its real transport ports, and the
same rate. `candidate_function_key` is null, because no retained GSPA EC claim is bound to these two
curated models in this benchmark; `missing_annotation` is therefore unavailable here and that is
reported as a structural limit of this family, not as a negative result.

Truth, fixed before any run: for an unsupplemented deletion assay the focal reaction really is absent,
so truth = `pathway_loss`; for the matched supplemented assay the measured growth is really due to the
supplied amino acid, so truth = `external_supply`. Wild-type assays are negative controls: no absence
explanation is true, and any uniquely selected explanation there is a false positive.

All six frozen medium/oxygen arms are run; none is chosen by its result.

Criteria:

* **A1-pass-1**: every constructible primary case executes the available hypotheses with a verified
  result and an independent LP check that agrees (below), with zero unexplained engine failures.
* **A1-pass-2**: recall of `pathway_loss` on unsupplemented deletion assays and of `external_supply`
  on supplemented deletion assays is >= 0.80 in at least one frozen arm, with its Wilson interval.
* **A1-pass-3**: unique-identification rate and precision are reported for all four explanations side
  by side, with the number of unavailable and non-comparable cells. No threshold: the bar states no
  level, so this is descriptive and a low value is reported as such.
* **A1-pass-4** (control): no wild-type assay yields a uniquely selected explanation. A failure here
  is a defect, not a performance result.

### A2. Controlled sequence loss (truth: missing sequence) - Shaheen

Inputs: the frozen sequence-loss cohort (`benchmark/v3/results/sequence-loss-v1/`, preparation
`fb76c32b...`): 24 removed labelled proteins and 13 retained labelled proteins in the four test
genomes, with the measured retained nucleotide fraction per genome; the gapsmith draft and filled
models of those genomes from `/scratch/project/k10018/gspa-v3/benchmark-models-v1`; the EC-to-reaction
map from the pinned gapsmith seed tables.

Per removed labelled protein whose labelled EC maps to a model reaction, one comparison on the
`sequence_loss` arm input and one on the `intact` arm input, with `assembly.completeness_fraction` set
to that genome's measured retained fraction (sequence-loss arm) or 1.0 (intact arm), and
`candidate_function_key` bound to the retained claim where the arm still has one. There is no measured
phenotype for these genomes, so no growth-matched feasible set exists; what is scored is availability
and rank:

Criteria:

* **A2-pass-1**: on the `sequence_loss` arm, `missing_sequence` is available in every case and
  `missing_annotation` is unavailable in every case whose candidate claim was removed with the
  sequence; on the `intact` arm the same cases have `missing_annotation` available. A failure is a
  defect of the construction.
* **A2-pass-2**: the rank of `missing_sequence` and the available-set composition are reported for
  every case, with the count of cases whose labelled EC maps to no model reaction (reported, not
  dropped).
* **A2-pass-3**: it is stated explicitly that without a measured phenotype no explanation is selected
  and no recall against growth is claimed for this family.

### A3. Price-148 hidden annotations (truth: missing annotation) - Shaheen

Inputs: the 148 Price et al. (2018) reannotated enzymes in 22 bacteria, the per-case masking of
`benchmark/v3/offgrowth/mask.py` (masking audit must pass for every case), the gapsmith draft/filled
models of those genomes, and the real predictor candidate table of the release-candidate Price
collection. Truth: the enzyme activity is present in the organism and experimentally assigned, but
absent from the masked annotation, so truth = `missing_annotation`.

Per case, one comparison with the focal reaction taken from the case's mapped seed reactions, the
candidate bound to the top-scoring retained claim for that enzyme, and the default medium of the
model. No measured phenotype exists per case, so the scored quantities are availability, rank, and
pairwise distinguishability of the four explanations under the model's own growth classes.

Criteria:

* **A3-pass-1**: `missing_annotation` is available in every case with a retained candidate, and the
  number of cases without one is reported.
* **A3-pass-2**: for every case the four explanations are reported side by side with their pairwise
  distinguishability; the fraction of cases in which `missing_annotation` is distinguishable from each
  of the other three is reported with a Wilson interval.
* **A3-pass-3**: if the four explanations are indistinguishable under the model's growth classes in
  most cases, that is reported plainly as the result.

### A4. Independent check for Part A

Every executed comparison is replayed by `benchmark/v3/check_absence_cobra.py` (COBRApy + GLPK, which
imports neither the production parser nor the production solver). **A4-pass**: every comparable cell
agrees within the script's tolerance; disagreements are listed individually and not summarised away.

## Part B (row 9): effect of false enzyme assignments on the coculture predictions

### B1. Injection

`inject_false_assignments.py` takes a model and (rate, seed) and rewrites gene-protein-reaction rules:
it draws, without replacement, `round(rate * n)` reactions that have a GPR, and replaces each drawn
reaction's whole GPR with one single gene drawn uniformly from the model's other gene products
(excluding that reaction's own genes). Draws use a SHA-256 stream keyed by (model SHA-256, rate, seed),
so the perturbation is reproducible from the receipt alone. The receipt records every (reaction, old
rule, new gene) triple and the resulting change in each coculture strain's blocked-reaction set.

Rates: 0.00, 0.01, 0.05, 0.10. Seeds: 1, 2, 3 (rate 0.00 uses one seed; it is the control). Both
models are perturbed independently with the same (rate, seed). No rate or seed is chosen after seeing
a result.

### B2. Rerun

All 7,920 frozen coculture plans (1,320 cocultures x six frozen arms) are rerun per (rate, seed)
through the same production community engine, with one change: each member's knockout set is
recomputed from the perturbed GPR, which is the effect under test. Everything else - medium, pool,
ports, ratio, budgets, growth threshold, solve order - is the frozen protocol. Scoring reuses the
frozen endpoint definitions (bulk growth, paired rescue proxy) and the measured observations.

Criteria:

* **B2-pass-1** (control): at rate 0.00 every one of the 7,920 plans reproduces the archived baseline
  prediction exactly. Any difference stops Part B and is reported as a defect.
* **B2-pass-2**: for each rate, the change in bulk and paired-rescue accuracy against the measured
  endpoints, per species and arm, is reported with the number of prediction flips and the spread over
  the three seeds, next to the unchanged constant baselines.
* **B2-pass-3**: the number of perturbed reactions that actually change a strain's blocked set is
  reported per (rate, seed); if most injected errors touch no coculture prediction, that is the
  reported result and the sensitivity is described as a lower bound for this design.
* **B2-pass-4**: a fraction of the perturbed runs (at least 1 in 20 plans, chosen by a content hash,
  not by outcome) is replayed by the independent COBRApy/GLPK community oracle
  (`benchmark/v3/check_nutrient_community_cobra.py`) and must agree.

### B3. Exchanged metabolites and coculture transporter knockouts

A published dataset with measured exchange and transporter knockouts that fits iML1515/iYO844 or the
existing community engine is searched for, with every attempt recorded (source, URL, HTTP outcome,
licence, what the paper actually measures). If one is obtained within this run with verified
provenance, it is built and scored under a design appended here before any label is read. Otherwise
the attempts and the reason for exclusion are reported, and these two elements stay MISSING.
Unconfirmed labels never count as negatives. **B3-pass**: an obtained dataset is scored, or the search
is reported with enough detail to repeat it.

## Rules carried over

* Defaults unchanged; production code paths used for every scored number. Touched tests are run.
* Heavy compute on Shaheen via Slurm (account k10018, `--exclude=nid00022,nid00497,nid01033`, explicit
  `REPO_ROOT`, Singularity options before the image, `timeout` on every ssh). Frozen jobs
  (`eqa_*`, `g_*`, `gc2*`, unimatrix 10406/10297/10298/10410/10432, PAVS) are not touched. No
  login-node compute. Nothing is written to `/tmp`.
* Outputs: `docs/benchmarks/v3-absence-crossfeeding-closure.md`,
  `benchmark/v3/results/absence-crossfeeding-closure-v1/` with `SHA256SUMS`, one progress line in
  `tasks/v3-supervision-plan-2026-10-28.md`. `tasks/v3-acceptance.md` and `README.md` are not edited;
  proposed row text goes in the final reply only.
* A bar element counts as shown only with real evidence. Elements measured and negative are
  distinguished from elements still missing.

## Deviations and incidents (appended after the checks ran; the plan above is unchanged)

All entries below were written when the change was made, before the affected numbers were read.

* **D1 (2026-11-01, before any Part B arm ran): targeted injection added.** A uniform draw over all
  reactions with a GPR almost never touches the six or nine genes the coculture benchmark actually
  deletes (at 10% of iYO844, one of six). Measuring only that would answer a different question, so a
  second, pre-declared injection mode was added: `targeted` draws only from the reactions whose rule
  mentions a deletion gene of the frozen registration, and replacement genes exclude the deletion
  genes so the wrong assignment is not self-cancelling. Arms: uniform rates 0.01/0.05/0.10 and
  targeted rates 0.1667/0.5/1.0, seeds 1-3 each, plus the uniform rate-0.00 control: 19 arms, each the
  full 7,920 plans. The frozen uniform scheme is reported unchanged beside it.
* **D2 (2026-11-01, at case construction): wild-type controls borrow a focal reaction.** A wild-type
  assay has no deleted gene, so the registered construction (focal reaction = the reaction the deleted
  gene blocks) yields nothing. Each wild-type assay is therefore paired with each primary auxotroph
  focal reaction of the same species, in the same environment and arm, where that reaction really is
  present. Truth stays "no absence explanation applies".
* **D3 (2026-11-01, at scoring): one case = one measured assay.** The plan's A1 defines truth per
  assay, so each of the two conditions of a comparison is scored as its own case with exactly one
  scored condition, against that assay's own measured OD600.
* **D4 (2026-11-02, at genome case construction): SBML gene identifiers normalised.** gapsmith writes
  gene products such as `G_AAM33030.1`; a dot is not a legal SBML SId, so libSBML rejects the file
  (error 10310) and the production parser cannot read it. Only identifiers beginning with `G_` are
  rewritten (dot to double underscore) into a copy; stoichiometry, bounds, objectives and reaction
  identifiers are untouched, and both hashes plus the count of rewritten identifiers are in the build
  receipt. Parts A2 and A3 could not run at all otherwise.
* **D5 (2026-11-02, at genome case construction): extracellular compartment declared.** gapsmith names
  the extracellular compartment `e0`, so `baseline.extracellular_compartments` is set to `["e0"]`
  explicitly and checked against the model rather than left to the engine's default of `["e"]`.
* **D6 (2026-11-02, for the independent check): the committed oracle is used, on a sample.** The
  check is `benchmark/v3/check_absence_cobra.py` exactly as committed (it compares status and growth
  value, not only the class). A first attempt wrote a separate oracle into that path; it was reverted
  with `git checkout` and the sampling wrapper now lives in
  `benchmark/v3/absence_crossfeeding/oracle_sample.py`. Samples are every 60th nutrient case, every
  sequence-loss case, every second Price case, chosen by file order and never by outcome.
* **D7 (2026-11-02, venue): split between the workstation and Shaheen.** Part A and the Part B
  rate-0.00 control ran on the 56-core workstation (strict memory overcommit there forced independent
  sharded processes instead of a worker pool); the 18 perturbed Part B arms ran one Shaheen node each
  (160 shards, jobs 15745622-15745639).
* **D8 (2026-11-02, finding about the inputs, not about the method): one archived cell is unknown.**
  For plan `NUTRIENT2026:Bs:AmSu + Gluc:R+C` in arm `printed:oxygen=0` the archived baseline has
  `no_unsupplied_uptake: null`; every Shaheen arm resolves it to `false`, including arms in which no
  knockout set changed, while the workstation control reproduces the archived `null`. It is an
  environment difference on one of 7,920 plans in a policy that no scored endpoint uses, and flip
  counts for that policy are reported with this constant named.
