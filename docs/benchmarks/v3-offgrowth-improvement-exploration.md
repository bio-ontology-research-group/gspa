# Off-growth enzyme discovery: exploratory analysis of the pilot outputs (2026-11-02)

Status: **exploratory**. Everything below was computed on the stored outputs of the off-growth pilot
(`benchmark/v3/results/offgrowth-pilot-v1/report-combined/summary.json.gz`, the per-genome tars under
`/scratch/project/k10018/gspa-v3/out/offgrowth-v1/genomes/` and the candidate tables next to them) and scored
against the **same labels the pilot used** (Price-148 and the 35 Swiss-Prot experimental enzyme numbers). No
number here is a result; the analysis exists to motivate the pre-registered plan in
`tasks/v3-row12-offgrowth-improvement-2026-11-02.md`, which must be tested on labels that were not used here.

Labels seen in this analysis: Price-148 case proteins (148) and the Swiss-Prot experimental labels (35 on 32
proteins). Labels **not** scored here: the Fitness Browser EC reannotations identified in §7 (their file was
written, their outcomes were not looked at), and the Swiss-Prot labels of the recurrence-cohort genomes.

Extraction: Shaheen job 15754030 (one node, 6 minutes), script
`~/tmp/gspa-v3-offgrowth-explore/shaheen/extract_pertest.py`, outputs copied to
`~/tmp/gspa-v3-offgrowth-explore/pertest/` (41 MB); analysis scripts `explore1.py`–`explore4.py` and
`feba_join.py` in `~/tmp/gspa-v3-offgrowth-explore/`.

## 1. Where the recoveries come from

Of the 180 hidden cases (148 Price, 32 Swiss-Prot), the union of the three targeted tests re-added the true
reaction in 54 and named the hidden protein in 23. The three test types contribute differently.

| Test types that recovered the case | Cases (of 23) |
|---|---:|
| Unblock + dead-end (A2 + A3) | 8 |
| Sole source only (A1) | 4 |
| Sole source + unblock (A1 + A2) | 3 |
| Unblock only (A2) | 3 |
| All three | 2 |
| Dead-end only (A3) | 2 |
| Sole source + dead-end (A1 + A3) | 1 |

The recoveries are concentrated in a few enzyme families that repeat across the seven Pseudomonas-group genomes,
so the 23 are fewer independent findings than they look.

| True enzyme | Re-added | Recovered | Note |
|---|---:|---:|---|
| 1.3.8.1 butyryl-CoA dehydrogenase | 8 | 0 | two true paralogs per genome; the call always goes to a third acyl-CoA dehydrogenase at 0.567 |
| 4.2.1.43 2-dehydro-3-deoxy-L-arabinonate dehydratase | 6 | 2 | |
| 5.3.1.7 mannose isomerase | 5 | 4 | |
| 4.3.1.9 glucosaminate ammonia-lyase | 5 | 0 | hidden protein has no candidate for this EC in 4 of 5 genomes |
| 1.1.1.48 galactose 1-dehydrogenase | 4 | 1 | |
| 3.2.2.3 uridine nucleosidase | 4 | 3 | |
| 3.1.1.15, 1.1.1.203, 4.1.2.19 | 2 each | 2, 2, 1 | |
| 13 further enzymes | 1 each | 8 | |

Per genome, recovery ranges from 0 of 13 (*Shewanella oneidensis*) to 6 of 31 (*Sinorhizobium meliloti*); 11 of
the 22 genomes contribute at least one recovery.

## 2. The sole-source tests fill 150–220 reactions each

Every arm gap-fills the **draft** model (the production gap-filled model is used only to define the growth-path
stratum). The gapsmith default medium contains glucose, lactate, all twenty amino acids and two vitamins, so the
draft was never required to synthesise amino acids. A sole-carbon-source test removes glucose, lactate and the
amino acids, and the linear program then rebuilds amino-acid biosynthesis together with the catabolic route for
the test compound. Fill sizes for one *Phaeobacter* case:

| Test type | Tests | Tests with a fill | Fill size: quartiles (25 %, median, 75 %, 90 %) | Largest |
|---|---:|---:|---|---:|
| A1 sole source | 293 | 240 | 150, 168, 193, 204 | 222 |
| A2 unblock | 1,005 | 646 | 0, 1, 4, 12 | 156 |
| A3 dead end | 421 | 216 | 0, 1, 4, 10 | 131 |

A reaction that comes back inside a 170-reaction fill is a weak call; a reaction that is the whole fill of a test
is a specific one. Across all union-arm gene calls, labelled precision (supported ÷ (supported + harmful)) falls
with the smallest fill in which the reaction was added:

| Smallest fill among the tests that added the reaction | Supported | Harmful | Unconfirmed | Labelled precision |
|---|---:|---:|---:|---:|
| 1 reaction | 99 | 150 | 19,244 | 0.40 |
| 2–3 | 28 | 312 | 24,132 | 0.08 |
| 4–10 | 18 | 73 | 26,398 | 0.20 |
| 11–50 | 0 | 22 | 6,125 | 0.00 |
| > 50 (all sole-source fills) | 15 | 335 | 28,588 | 0.04 |
| not attributable (added only in a test type not extracted) | 30 | 121 | 24,676 | 0.20 |

By test type: A2 calls at fill ≤ 3 have precision 0.22 (117 / 528), A3 calls at fill ≤ 3 have 0.28 (88 / 317),
and A1 calls (all at fill > 10) have 0.10 (55 / 529). All ten sole-source recoveries sit inside fills of
77–194 reactions; in three of them the same reaction was also re-added by an unblock or dead-end test with a fill
of 1–4. Wherever the same enzyme was found in two genomes by different test types (mannose isomerase: fill 168 in
*Phaeobacter* by A1, fill 1 in *Paraburkholderia* and *Pseudomonas simiae* by A2/A3), the small fill was the
specific one.

## 3. The gene-assignment ceiling

For each hidden case the unmasked candidate table (best score per protein and enzyme over the six predictors) was
inspected for the true enzyme number.

| Outcome of the union arm | Cases | Hidden protein is rank 1 for its EC | Rank 2–3 | Rank > 3 | No candidate for its EC | Median own score |
|---|---:|---:|---:|---:|---:|---:|
| Recovered | 22 | 21 | 0 | 0 | 1 | 0.566 |
| Reaction re-added, wrong protein named | 30 | 0 | 5 | 15 | 10 | 0.202 |
| Reaction not re-added | 99 | 32 | 27 | 29 | 11 | 0.552 |
| All (per truth EC row; 186 rows for 180 cases) | 186 | 71 | 37 | 47 | 31 | |

Two facts follow. Recovery by "top-scoring candidate" can never exceed the 71 cases in which the hidden protein is
already the best candidate for its enzyme (38 %); 21 of those 71 were recovered, and 32 more sit in cases whose
reaction no test re-added, so the simulation side, not the ranking, is what limits recall among recoverable cases.
Conversely, in all 30 "right reaction, wrong gene" cases the hidden protein ranked below another protein or had
no candidate at all, so these calls were wrong by construction and no re-weighting of the same ranking can fix
them (the competitive-assignment follow-up confirmed this: 25 of 54).

Whether the true enzyme is the protein's own best enzyme ("mutual best") separates less well than hoped: 10 of 22
recovered cases are mutual best, 1 of 30 wrong-gene cases, 17 of 99 not-re-added cases. As a gate on calls, mutual
best raises labelled precision from 0.12 (not-best, 848 calls on labelled proteins) to 0.25 (355 calls) but would
drop half of the recoveries. The median margin of the hidden protein over the best other protein for its EC is
+0.18 in recovered cases and −0.31 in wrong-gene cases, so a margin rule is informative where the hidden protein
is rank 1.

## 4. Anatomy of the harm

The 1,013 harmful calls of the union arm are **62 distinct (genome, protein, enzyme) decisions**, repeated once per
hidden case of the same genome (31 cases in *S. meliloti*); the 190 supported calls are 23 distinct decisions. The
per-call harm rate of 0.84 therefore rests on 85 distinct decisions, with the genomes that have many cases
weighted most. Of the harmful calls, 445 (44 %) name an enzyme in the same sub-subclass as the protein's true
enzyme (first three EC levels equal: acyl-CoA dehydrogenases, sugar dehydrogenases), 324 share one level and 181
none. Roughly half of the harm is a related-substrate mistake that a growth simulation cannot distinguish, the
other half is a dense predictor giving a labelled protein a top score for an unrelated enzyme.

Score helps little: labelled precision by 0.05-wide score bin peaks at 0.37 in the 0.55–0.60 bin (124 supported,
207 harmful) and is 0 above 0.60 (169 harmful, 0 supported; one predictor's ceiling of about 0.66). The maximum
over six differently calibrated predictors is not a usable confidence.

Agreement across the three test types raises precision from 0.08 (one type) to 0.25 (two) and 0.22 (three). How
many different reactions the same protein is called for in a case does not separate right from wrong.

Calls that appear in **every** case of a genome are what an unmasked run would produce (the hidden case's own call
is by definition absent from the other cases). In the 18 genomes with at least three cases there are a median
660 distinct gene calls and 440 reaction-only calls per genome; the labelled decisions among them are 11
supported and 80 harmful. Restricted to reactions whose smallest fill is one reaction: 128 gene calls per genome,
7 supported, 12 harmful (0.37). The labelled denominators of a deployment-like view are tiny, which is the main
reason the pilot's precision could not be measured usefully.

## 5. Why standard gap-filling and the measured phenotypes gave nothing

Only 2 of 180 cases lie on the default growth path, so biomass gap-filling has no reason to add them. The RB-TnSeq
follow-up could not test the gap-filled models because they grow on none of the sole-source wells; the reason is
the same medium effect as in §2: a model filled for growth on glucose plus amino acids has no amino-acid
biosynthesis to use when a single carbon source is offered.

## 6. What the exploration suggests (hypotheses, not findings)

1. Run sole-source tests **differentially**: first fill the draft once for growth on a minimal medium (the baseline
   fill, shared by all compounds), then fill each compound relative to that model. The specific fill for a compound
   should then be the 1–5 reactions of its catabolic route, and sole-source calls would inherit the precision seen
   for small fills instead of the 0.04–0.10 of 150-reaction fills. The same baseline fill lets the model grow on
   measured sole-source wells, which removes the obstacle that stopped the RB-TnSeq comparison.
2. Gate calls on **fill specificity** (reaction added by a test whose fill is at most F reactions) and on
   **fill stability** (present under two solvers or two weightings; the pilot's replay agreed on the reaction set
   in only 47 % of gap-fills).
3. Assign a gene only when the sequence ranking is itself unambiguous (rank 1 for the reaction's enzyme with a
   margin over rank 2, a score floor, and either mutual best or no existing model assignment); otherwise record a
   **reaction-level** call ("the organism needs this reaction; candidates X, Y") and abstain on the gene.
4. Group reactions of one small fill into **pathway calls** and test whether their candidate genes cluster on the
   chromosome; recoveries such as rhamnose isomerase and rhamnulose-phosphate aldolase (adjacent genes in
   *S. meliloti*) and the Fitness Browser's inositol operon in the same genome suggest that cluster calls can carry
   precision that single-gene calls cannot.
5. Measure precision on **distinct decisions of an unmasked run per genome**, scored against a much larger label
   set (§7), and keep masked per-case runs only for recovery. The two numbers answer different questions and
   should not be mixed.

## 7. A larger label set of the same kind exists and has not been used

The Fitness Browser database (`feba.db`, downloaded 2026-10-31 for the RB-TnSeq follow-up) carries the curated
reannotations of the Price laboratory in tables `Reannotation` (1,886 loci, 41 organisms) and `ReannotationEC`
(1,204 enzyme-number rows). Joining them to the 22 Price proteomes by exact protein sequence
(`~/tmp/gspa-v3-offgrowth-explore/fitness_browser_ec_reannotations_price22.tsv`, 763 rows):

| Set | Proteins | With a complete four-level EC |
|---|---:|---:|
| Fitness Browser EC reannotations matched in the 22 genomes | 763 | 747 |
| of which Price-148 case proteins | 132 | 132 |
| of which Swiss-Prot experimental label proteins | 0 | 0 |
| **of which neither (never used as truth or harm label)** | **631** | **615** |
| Fitness Browser EC reannotations in 17 other organisms (not among the 145 benchmark genomes) | 302 | 293 |

Per genome the new labels number 4 (*Azospira*) to 65 (*Pseudomonas simiae*, *Pseudomonas* N2E3); the
Pseudomonas-group genomes hold about 350 of the 615. The outcomes of the pilot's calls on these proteins were
**not** computed. The recurrence-cohort bacteria among the 145 benchmark genomes (eight genomes: *Synechococcus
elongatus* 7942, *Campylobacter jejuni* NCTC 11168, *Legionella pneumophila*, *Clostridium acetobutylicum*,
*Magnetospirillum magneticum*, *Haemophilus influenzae* Rd, *Acinetobacter* ADP1, *Dickeya dadantii* 3937) carry
79 Swiss-Prot experimental enzyme labels in `/project/k10018/gspa-v3/training-cohort/labels/labels.tsv`; their
gapsmith models exist, their six-predictor candidate tables do not.

## Limits of this analysis

* Scored on the labels the pilot was scored on; any rule read off these tables is tuned to them.
* Fill sizes come from the pilot's weighted linear programs (weights 1 / 100, draft model, GLPK); a differential
  design changes the programs, so the small-fill precision above is an indication, not a prediction.
* Distinct-decision counts treat identical calls across cases of one genome as one decision; the pilot's own
  harm figures count them per case.
