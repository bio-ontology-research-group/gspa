# Which absence explanation does the evidence support, and what do wrong enzyme assignments do to the coculture predictions

Date: 2026-11-02. Pre-registered plan: `tasks/v3-row8-row9-absence-crossfeeding-2026-11-01.md`
(body SHA-256 `f04f95466fe042f58eff1a42fb1877630203e0d0e77367e32780db0a66f9ac13`, frozen before any scored
check; the deviations appended to it were each written before the affected numbers were read).
Results and receipts: `benchmark/v3/results/absence-crossfeeding-closure-v1/`.

This closes two measurement gaps. The four-scenario absence comparison existed and was audited as a
mechanism, but had never been scored against cases whose true explanation is known. The coculture
benchmark measured growth agreement, but the effect of a wrong enzyme assignment on those predictions
had never been measured. Growth and knockout agreement themselves are not re-reported here; they stay
in `v3-nutrient-phenotypes.md` and `v3-nutrient-community-results.md`.

## What was scored, and how

A comparison is one production `absence_comparison` run (`scripts/absence_experiments.py`, schema
`gspa.absence_experiments/1`): four explanations — a missing annotation, a missing sequence, a genuine
pathway loss, and external supply of the product — evaluated over explicit conditions. An explanation
is *available* when the engine executes it; it is *feasible* when its evaluated growth class equals
the measured class in the scored condition. Feasibility is conditional on the model, the medium and
the bounds: it is evidence about a reconstruction, not biological proof of a cause.

| Family | Cases | True explanation | Measured phenotype |
|---|---|---|---|
| Nutrient auxotroph deletions (`iML1515`, `iYO844`) | 2,640 primary comparisons (5,280 primary assay cells; 9,504 with the sensitivity scope and the wild-type controls) | pathway loss (unsupplemented), external supply (supplemented) | yes, OD600 per assay |
| Controlled sequence loss (3 of 4 test genomes) | 64 comparisons, 22 with a true explanation | missing sequence | no |
| Price et al. reannotations (22 genomes) | 100 of 149 constructible | missing annotation | no |

## Auxotroph deletions: the true explanation is usually feasible, but rarely alone

Recall is the fraction of cases whose true explanation is in the feasible set; precision is how often
a case in which that explanation is feasible really has it as the truth. Intervals are Wilson 95%.

| Explanation | Cases with this truth | Recall | Precision | Uniquely identified |
|---|---|---|---|---|
| pathway loss | 2,640 | 0.861 [0.847, 0.874] | 0.637 | 0.285 |
| external supply | 2,640 | 0.490 [0.471, 0.509] | 0.465 | 0.000 |
| missing sequence | 0 (never the truth here) | — | 0.000 | — |
| missing annotation | 0 | unavailable in all 9,504 cells | — | — |

A clustered bootstrap over the 440 (species, strain, environment) groups gives 0.831-0.889 for
pathway loss and 0.473-0.508 for external supply; it is descriptive, not a population estimate.

Three things in that table matter more than the point values.

External supply is recovered only where the medium makes it visible. Its recall is 0.864 [0.828,
0.893] and 0.889 [0.856, 0.915] in the two aerated `additional_inorganic` arms, and 0.130 [0.101,
0.164] in `printed:oxygen=0`, where the model cannot grow on the arm's medium at all, so the supplied
amino acid changes nothing it can be seen by. The same comparison, on the same strains, answers well
or not at all depending on the condition it is given.

No case ever identified external supply uniquely, and only 28.5% identified pathway loss uniquely:
with a single binary growth endpoint, disabling a reaction and withdrawing its product usually produce
the same prediction. This is a property of the design, not a solver failure, and it is the honest
ceiling for one-endpoint evidence.

A missing annotation cannot be scored on this family at all. The engine makes that explanation
available only with a retained candidate claim, and no GSPA claim is bound to these two curated
models in this benchmark. That is a structural limit of the family, reported as such, not a negative
result.

The wild-type controls, where the focal reaction really is present, show what the method costs in
false positives: of 2,640 control cells, 166 (6.3% [5.4%, 7.3%]) selected exactly one explanation, and
every one of those 166 is a case where the unperturbed model already disagrees with the measurement.
Where the model reproduces the wild type (1,144 cells), no control ever selects an explanation.

## Sequence loss and complete genomes: availability behaves correctly, growth classes do not separate

For a removed labelled protein the engine made the missing-sequence explanation available in 22 of 22
cases, and made the missing-annotation explanation unavailable in exactly those 22 while keeping it
available in all 10 retained-labelled cases of the same arm and in all 32 intact-arm cases. The
availability contrast that distinguishes "the sequence is gone" from "the annotation is gone" is
therefore produced correctly, by construction, from the arm's own evidence. For the Price enzymes the
missing-annotation explanation was available in all 100 constructible cases; 49 of 149 were not
constructible because the enzyme's mapped seed reactions are absent from that genome's model, which is
reported rather than dropped.

Under the models' own growth classes, however, the explanations do not separate: in 99 of 100 Price
cases and 31 of 32 sequence-loss cases every pair carries the same evaluated growth class, and the
external-supply explanation is unavailable throughout (no supply is declared). Exactly one case per
family is distinguishable. Neither family has a measured phenotype, so nothing is selected and no
recall against growth is claimed for them. The useful conclusion is negative and specific: on whole
genomes these four explanations are not separable by a default-medium growth prediction, and
separating them needs a condition in which they disagree, which is what the auxotroph family supplies.

## Independent check

A sample was replayed through the committed independent oracle, `benchmark/v3/check_absence_cobra.py`
(COBRApy 0.32.1 on GLPK, which compares status and growth value, not only the class): every
sequence-loss case (296 cells), every second Price case (50 cases, 300 cells) and every 60th nutrient
case (102 cases, 612 cells). All 216 sampled cases agree on status and growth class. In one nutrient
case the two solvers' objective values differ by 1.3e-6 relative (0.1024002049 against 0.1024000715)
across four cells, which is just outside the oracle's 1e-6 relative tolerance, so it is reported as a
disagreement rather than hidden: the registered criterion asked for agreement within that tolerance,
and this cell is outside it while agreeing on every class.

## Wrong enzyme assignments barely move the measured coculture endpoints

Each arm replaces the gene-protein-reaction rule of a drawn set of reactions with one wrong gene,
keyed reproducibly by (model hash, mode, rate, seed), and reruns all 7,920 frozen coculture plans
through the production community engine; only the members' knockout sets change. The control at rate
0.00 reproduced all 7,920 archived predictions exactly, and the scoring reproduces the archived
baseline endpoint accuracies exactly (bulk growth 0.4875, paired rescue proxy 0.6125 over 7,200
primary cells), so the perturbed numbers are comparable with the published ones.

| Arm | Reactions perturbed (Bs/Ec) | Plans with a changed knockout set | Plans identical to baseline | Bulk growth | Paired rescue |
|---|---|---|---|---|---|
| control, rate 0.00 | 0/0 | 0 | 7,920 | 0.4875 | 0.6125 |
| uniform 0.01, seeds 1-3 | 9/23 | 0-1,320 | 7,314-7,919 | 0.4875 | 0.6125-0.6364 |
| uniform 0.05, seeds 1-3 | 45/113 | 0-3,696 | 7,217-7,919 | 0.4225-0.4875 | 0.5903-0.6364 |
| uniform 0.10, seeds 1-3 | 90/227 | 0-6,072 | 7,529-7,919 | 0.4875 | 0.5947-0.6125 |
| targeted 0.1667, seeds 1-3 | 1/2 | 1,320-2,640 | 7,529 | 0.4875 | 0.5969-0.6003 |
| targeted 0.5, seeds 1-3 | 3/4 | 5,544-6,336 | 6,128-6,612 | 0.4875 | 0.5775-0.6044 |
| targeted 1.0, seeds 1-3 | 6/9 | 7,920 | 5,055 | 0.4875 | 0.6075 |

The targeted arms are the strong test: at rate 1.0 every reaction of every deleted gene is reassigned
to a gene that is not deleted, so in the model no strain is auxotrophic any more. That changes 2,865
of 7,920 plans in at least one policy — 2,708 monoculture predictions for the first member, 1,662 for
the second — and yet bulk growth agreement does not move at all (0.4875) and paired rescue moves from
0.6125 to 0.6075. Only one arm of nineteen changed the community prediction itself (uniform 0.05 seed
2, 605 `permissive` flips, bulk growth 0.4875 → 0.4225).

The reason is visible in the flip counts: the perturbation moves monoculture and native-uptake
predictions readily, but the permissive community endpoint that carries the bulk-growth comparison is
set by the shared pool and the medium rather than by which reactions a strain has lost. So the
measured sensitivity is a lower bound for this design, and the result to carry forward is about the
benchmark rather than about the engine: these coculture endpoints are weak at detecting enzyme
misassignment, and a sharper design would score monoculture-level or exchange-level endpoints where
the perturbation does land. An unchanged prediction is not evidence that a wrong assignment is
biologically harmless.

One of the 7,920 plans (`NUTRIENT2026:Bs:AmSu + Gluc:R+C`, arm `printed:oxygen=0`) carries a null
`no_unsupplied_uptake` in the archive which every Shaheen arm resolves to false, including arms where
no knockout set changed, while the workstation control reproduces the archived null. It is an
environment difference in a policy no scored endpoint uses; the single flip it contributes to that
policy is named in the receipts.

## Measured exchange and coculture transporter knockouts: still missing, with the search recorded

No usable published dataset was obtained, and the two elements stay MISSING rather than being filled
with something weaker. What was tried, and what it showed:

* Mee et al. 2014 (PNAS 111:E2149, PMC4034247) is the closest design — all 91 pairs of 14 single
  amino-acid auxotrophs of one *E. coli* background, 84 h fold growth, with qPCR relative abundances
  in Dataset S2. Its own text, however, reports no transporter knockouts and no measured secreted
  amino-acid concentrations, so it would not supply either missing element even if obtained. The files
  were not obtainable here in any case: `pnas.org` returned HTTP 403 behind Cloudflare, the PMC
  supplementary paths returned a reCAPTCHA page rather than the spreadsheet, and the Europe PMC
  `supplementaryFiles` service returned an empty archive for this PMCID.
* Harcombe et al. 2014 (Cell Reports 7:1104, PMC4097880) measures the *E. coli*/*S. enterica*
  methionine-acetate exchange that the benchmark would need. It is not open access: `cell.com`
  returned HTTP 403, Europe PMC reports the record as not open access, and the one PMC supplementary
  file (`NIHMS585160-supplement-01.pdf`) is served behind the same reCAPTCHA.
* Pande et al. 2014 (ISME J 8:953) and the Kost-lab follow-ups use deletion-based cross-feeding pairs
  but report no exporter knockouts and no deposited exchange measurements; searches of Dryad and
  Zenodo for a coculture transporter-knockout dataset fitting `iML1515`/`iYO844` returned nothing.

An unconfirmed label is not a negative result, so nothing from these sources was scored. Obtaining
the Mee Dataset S2 and the Harcombe supplement needs either institutional access or a request to the
authors; the exchange and transporter-knockout elements of row 9 should stay marked missing until then.

## Scope and limits

* Feasibility and distinguishability are statements about these reconstructions under these media and
  bounds, not about biology. "No absence explanation is feasible" can mean the model is wrong rather
  than that no explanation holds.
* The nutrient family scores three explanations; a missing annotation is unavailable there for want of
  a retained candidate claim. The genome families score availability, rank and distinguishability
  only, never selection, because they have no measured phenotype per case.
* The Price candidate key is the enzyme's own labelled EC, so that family measures whether the
  comparison keeps a missing annotation available and separable, never whether it can find the right
  candidate.
* The wild-type controls borrow a focal reaction from the matched auxotroph of the same species; a
  wild type has no deletion of its own.
* The false-assignment injection is a deliberately crude model of an enzyme assignment given to the
  wrong gene, and the sensitivity it measures is a lower bound for the endpoints this benchmark scores.
* The 49 Price enzymes whose mapped reactions are absent from their genome's model, the 264 E. coli
  `leuB` cases with no blocked reaction (isozyme bypass), the 792 `serA` comparisons kept in a separate
  sensitivity scope, and the fourth sequence-loss test genome with no EC labels are all reported, and
  none of them is silently dropped.

## Reproducing

```bash
# Part A, auxotrophs (workstation; 12 sharded processes, ~25 min)
python benchmark/v3/absence_crossfeeding/absence_truth.py build  --plan-dir <frozen nutrient plan> \
    --ec-model iML1515.xml --bs-model iYO844.xml --out cases
python benchmark/v3/absence_crossfeeding/absence_truth.py run    --cases cases --out shard --shard-index I --shard-count 12 ...
python benchmark/v3/absence_crossfeeding/absence_truth.py merge  --cases cases --shards shard0 ... --out run
python benchmark/v3/absence_crossfeeding/absence_truth.py score  --cases cases --run run --observations observations.jsonl --out score.json

# Part A, genomes
python benchmark/v3/absence_crossfeeding/genome_absence.py build-sequence-loss --preparation ... --labels ... --models-dir ... --out cases --work normalised
python benchmark/v3/absence_crossfeeding/genome_absence.py build-price --manifest benchmark/manifests/suite6_price149_manifest.tsv --models-dir ... --out cases --work normalised
python benchmark/v3/absence_crossfeeding/genome_absence.py run --cases cases --out run
python benchmark/v3/absence_crossfeeding/genome_absence.py score --cases cases --run run --out score.json

# Independent oracle (committed COBRApy/GLPK script, sampled)
python benchmark/v3/absence_crossfeeding/oracle_sample.py --cases cases/cases.jsonl --every 60 --work tmp --out oracle.json

# Part B, one arm per Shaheen node (jobs 15745622-15745639)
MODE=targeted RATE=1.0 SEED=1 sbatch benchmark/v3/shaheen/crossfeeding-false-assignment.sbatch
python benchmark/v3/absence_crossfeeding/false_assignments.py score --plan-dir <frozen coculture plan> \
    --baseline baseline_predictions.jsonl --arms <19 merged arms> --observations observations.jsonl --out score.json
```

Tests: `benchmark/v3/tests/test_absence_crossfeeding.py` (6), with
`benchmark/v3/tests/test_absence_cobra.py`, `scripts/tests/test_absence_experiments.py`,
`test_nutrient_community.py` and `test_nutrient_phenotypes.py` all passing unchanged.
