# Completed context comparison: independent audit

The computation is consistent, but this comparison does not establish a broad
benefit from pathway or organism-specific context. Learned scoring improves some
reference-based results, especially after family evidence is excluded, while
candidate coverage remains poor. There is no matched CLEAN comparison here.
No model is promoted and all twelve whole audit items remain open.

## Population and what the numbers mean

This controlled experiment uses one DIAMOND homology-transfer source; it does
not compare the complete GSPA predictor ensemble. The [candidate trace and
reporting repair](v3-ec-candidate-coverage.md) explain the missing functions and
make complete-reference EC panels part of future schema-2 reports.

The original run covers all 18 registered genomes: 11 training, three calibration
and four test genomes. The test genomes are `GCA_000007065.1`, `GCA_000008765.1`,
`GCA_000009985.1` and `GCA_000147055.1`. It compares three contexts, each with
rule, logistic and boosted-tree scoring, in two transfer regimes. UNIVERSAL uses
the universal reaction context, PATHWAY uses pathway grouping and ORGANISM uses
organism-specific context. Family exclusion removes the specified homology
transfer evidence; it does not establish absence of overlap in pathway-search
references or predictor pretraining.

The original test tables contain 608 homology rows (75 positive) and 416
family-excluded rows (34 positive). Their EC subsets contain only 52 rows with
15 positives and 22 rows with five positives. The remaining labels are
**unconfirmed reference nonmatches**, not experimentally established negatives.
Every row, including learned-model fallback, remains in the comparison.

The independent audit also evaluates **all 37 experimental-reference proteins and
all 39 positive EC assertions** in the four test genomes. Exact EC identifiers
and evidence code `ECO:0000269` define this panel. Missing proteins and functions
remain in the denominator. This is the complete panel under that evidence rule,
not every enzyme in these genomes. Neither predicted functions on other proteins
nor inferential-only reference assertions are covered by this expanded panel.
The full candidate matrix is used for ranking within these 37 proteins.
The 39 EC assertions occur in three of the four test genomes: 21 in
`GCA_000007065.1`, 13 in `GCA_000008765.1` and five in `GCA_000147055.1`.
`GCA_000009985.1` contributes no EC assertion under this rule; that is missing
reference coverage, not evidence that its proteins lack enzyme functions.

All selected thresholds were fitted on the mixed GO/EC calibration rows, with
higher cutoffs winning equal-F1 ties. The audit preserves these thresholds and
0.5; it never retunes on the EC panel or test rows. The expanded panel is a
post-hoc descriptive analysis of an already examined cohort. Calibration fitting
and threshold selection reuse the calibration split, so calibration F1 is not an
independent estimate either.

## Original GO and EC test rows combined

M = reference matches, U = unconfirmed calls, FN = missed positive reference rows.
Triples below are M/U/FN; all denominators are fixed within a transfer regime.
The full confusion counts, precision, recall, Brier, ECE, 48 paired contrasts and
per-genome counts remain in the [original results](../../benchmark/v3/results/context-comparison-audit-v1/analysis/results.json).
Four-genome bootstrap intervals describe this panel; they are not evidence of
population significance. The two transfer regimes overlap and are not independent
replications.

| Regime | Context/scorer | Selected cutoff | M/U/FN | F1 | At 0.5: M/U/FN | F1 at 0.5 | Brier |
|---|---|---:|---|---:|---|---:|---:|
| homology | UNIVERSAL/rule | 0.3072 | 23/61/52 | 0.2893 | 3/4/72 | 0.0732 | 0.1081 |
| homology | UNIVERSAL/logistic | 0.3126 | 14/29/61 | 0.2373 | 6/7/69 | 0.1364 | 0.1042 |
| homology | UNIVERSAL/gbdt | 0.2029 | 22/63/53 | 0.2750 | 0/0/75 | 0.0000 | 0.1090 |
| homology | PATHWAY/rule | 0.3072 | 23/61/52 | 0.2893 | 3/4/72 | 0.0732 | 0.1081 |
| homology | PATHWAY/logistic | 0.3133 | 14/29/61 | 0.2373 | 5/6/70 | 0.1163 | 0.1055 |
| homology | PATHWAY/gbdt | 0.2674 | 14/36/61 | 0.2240 | 2/5/73 | 0.0488 | 0.1094 |
| homology | ORGANISM/rule | 0.3072 | 23/61/52 | 0.2893 | 3/4/72 | 0.0732 | 0.1081 |
| homology | ORGANISM/logistic | 0.3752 | 8/15/67 | 0.1633 | 4/7/71 | 0.0930 | 0.1040 |
| homology | ORGANISM/gbdt | 0.1921 | 23/69/52 | 0.2754 | 1/4/74 | 0.0250 | 0.1099 |
| family_excluded | UNIVERSAL/rule | 0.1780 | 6/94/28 | 0.0896 | 0/1/34 | 0.0000 | 0.0808 |
| family_excluded | UNIVERSAL/logistic | 0.2058 | 6/33/28 | 0.1644 | 0/2/34 | 0.0000 | 0.0760 |
| family_excluded | UNIVERSAL/gbdt | 0.2133 | 12/37/22 | 0.2892 | 0/1/34 | 0.0000 | 0.0753 |
| family_excluded | PATHWAY/rule | 0.1780 | 6/94/28 | 0.0896 | 0/1/34 | 0.0000 | 0.0808 |
| family_excluded | PATHWAY/logistic | 0.2029 | 5/29/29 | 0.1471 | 0/1/34 | 0.0000 | 0.0765 |
| family_excluded | PATHWAY/gbdt | 0.2272 | 8/24/26 | 0.2424 | 0/1/34 | 0.0000 | 0.0762 |
| family_excluded | ORGANISM/rule | 0.1780 | 6/94/28 | 0.0896 | 0/1/34 | 0.0000 | 0.0808 |
| family_excluded | ORGANISM/logistic | 0.2373 | 2/4/32 | 0.1000 | 1/2/33 | 0.0541 | 0.0751 |
| family_excluded | ORGANISM/gbdt | 0.2102 | 11/35/23 | 0.2750 | 0/1/34 | 0.0000 | 0.0737 |

With homology evidence, no learned arm exceeds the rule arm's combined F1 at the
selected cutoffs. With family evidence excluded, UNIVERSAL trees achieve the
highest combined F1 (0.2892 versus 0.0896 for rules). Changing context does not
consistently improve the same scorer: ORGANISM trees have a lower Brier score
than UNIVERSAL trees in that regime, but lower F1. A single favorable metric
cannot establish a general context benefit.

## All experimental-reference EC proteins

Only **15/39 reference EC assertions are emitted** with homology evidence and
**5/39** after family exclusion, identically across context modes and scorers.
There are EC candidates of any identity for 30/37 and 13/37 proteins respectively.
Thus, reranking the current candidates cannot exceed 38.5% and 12.8% reference
recall on this panel. A missing candidate is different from a candidate rejected
by the learned model's support checks.

Each row below uses the same 39 positive assertions. F1 and precision penalize
reference nonmatches for comparison purposes; those calls remain biologically
unconfirmed. Top-1 is the number of proteins with a reference match ranked first,
out of all 37, using lexical EC order to resolve equal scores. It does not apply
a score threshold. Ties and top-3/top-5 are retained in the archive.

| Regime | Context/scorer | Selected: M/U/FN | Reference F1 | Recall | At 0.5: M/U/FN | Reference F1 at 0.5 | Top-1 /37 |
|---|---|---|---:|---:|---|---:|---:|
| homology | UNIVERSAL/rule | 7/9/32 | 0.2545 | 0.1795 | 0/1/39 | 0.0000 | 14 |
| homology | UNIVERSAL/logistic | 11/8/28 | 0.3793 | 0.2821 | 4/3/35 | 0.1739 | 14 |
| homology | UNIVERSAL/gbdt | 10/9/29 | 0.3448 | 0.2564 | 0/0/39 | 0.0000 | 13 |
| homology | PATHWAY/rule | 7/9/32 | 0.2545 | 0.1795 | 0/1/39 | 0.0000 | 14 |
| homology | PATHWAY/logistic | 11/8/28 | 0.3793 | 0.2821 | 3/2/36 | 0.1364 | 14 |
| homology | PATHWAY/gbdt | 7/6/32 | 0.2692 | 0.1795 | 1/1/38 | 0.0488 | 12 |
| homology | ORGANISM/rule | 7/9/32 | 0.2545 | 0.1795 | 0/1/39 | 0.0000 | 14 |
| homology | ORGANISM/logistic | 6/3/33 | 0.2500 | 0.1538 | 2/3/37 | 0.0909 | 12 |
| homology | ORGANISM/gbdt | 11/9/28 | 0.3729 | 0.2821 | 0/0/39 | 0.0000 | 11 |
| family_excluded | UNIVERSAL/rule | 1/6/38 | 0.0435 | 0.0256 | 0/0/39 | 0.0000 | 4 |
| family_excluded | UNIVERSAL/logistic | 2/6/37 | 0.0851 | 0.0513 | 0/0/39 | 0.0000 | 4 |
| family_excluded | UNIVERSAL/gbdt | 5/5/34 | 0.2041 | 0.1282 | 0/0/39 | 0.0000 | 4 |
| family_excluded | PATHWAY/rule | 1/6/38 | 0.0435 | 0.0256 | 0/0/39 | 0.0000 | 4 |
| family_excluded | PATHWAY/logistic | 1/4/38 | 0.0455 | 0.0256 | 0/0/39 | 0.0000 | 3 |
| family_excluded | PATHWAY/gbdt | 3/5/36 | 0.1277 | 0.0769 | 0/0/39 | 0.0000 | 4 |
| family_excluded | ORGANISM/rule | 1/6/38 | 0.0435 | 0.0256 | 0/0/39 | 0.0000 | 4 |
| family_excluded | ORGANISM/logistic | 2/2/37 | 0.0930 | 0.0513 | 1/1/38 | 0.0488 | 3 |
| family_excluded | ORGANISM/gbdt | 4/3/35 | 0.1739 | 0.1026 | 0/0/39 | 0.0000 | 5 |

The strongest homology EC reference F1 is 0.3793, with 11 matches and eight
unconfirmed calls; rules produce seven matches and nine unconfirmed calls.
With family evidence excluded, UNIVERSAL trees produce five matches and five
unconfirmed calls (F1 0.2041), versus one match and six unconfirmed calls for
rules (F1 0.0435). These are small absolute counts. The attractive selected-row
EC F1 values—0.6471 for homology UNIVERSAL logistic and 0.6667 for family-excluded
UNIVERSAL trees—fall to 0.3793 and 0.2041 when missing reference functions are
included. The expanded denominator changes the interpretation substantially.

At 0.5, no boosted-tree arm produces a reference match except PATHWAY with
homology (one match); homology UNIVERSAL logistic produces four. Most
family-excluded arms make no EC calls at this cutoff; ORGANISM logistic produces
one match and one unconfirmed call. This illustrates why improvements measured
at a lower benchmark cutoff need not appear as more displayed enzyme calls.
It is not a direct audit of a particular web-interface run, and no default
threshold was changed. Lowering a cutoff alone would not solve the missing
candidate problem or establish calibrated probabilities.

Top-1 ranking also provides no general improvement: homology rules and UNIVERSAL
logistic both match 14/37 proteins, while UNIVERSAL trees match 13/37. Across all
homology arms, top-3 and top-5 match 15/37 proteins. In the family-excluded regime,
top-1 ranges from three to five matches; all arms reach five at top-5. The
archived records retain all 37 proteins, including empty predictions and ties.

## Execution, support and independent verification

Original job 10176 finished at 2026-10-27T21:25:23+03:00 with all seven stages
successful. It fitted all 12 models and completed 216 full-production runs,
4,539,408 candidate/model checks, and 40,236 JVM/Python table-reader comparisons.
Maximum table-reader disagreement was 6.11e-16. Earlier exports cover all 144
registered genome/arm combinations; their failed relocation attempts remain
recorded in the [execution history](v3-observed-gap-context-modes.md).

Independent job 10223 verified **6,981 bound files (87,052,664,741 bytes)**,
including all seven supplementary supervision dependencies against their parent
registration. The registry was already bound there; previous prose suggesting
otherwise was incorrect. It reconstructed all 6,706 scored supervision rows
from production TSVs, reading 6,809,112 arm/candidate rows, and verified unchanged
supervision metadata across contexts. It independently reproduced all 18
calibration thresholds, 36 primary summaries and 48 paired statistics, then
rechecked the 337 consumed bound inputs. Five focused standard-library tests
pass. No model was refitted or rerun for this audit.

Supported production candidates below cover **all 18 genomes**, not just labelled
test rows. Both learned scorers in each mode have the same support counts.
Unsupported candidates retain rule fallback and remain in the evaluation.

| Regime | Context | Supported | All production candidates |
|---|---|---:|---:|
| homology | UNIVERSAL | 343944 | 484202 |
| homology | PATHWAY | 344812 | 484202 |
| homology | ORGANISM | 343164 | 484202 |
| family_excluded | UNIVERSAL | 196621 | 272366 |
| family_excluded | PATHWAY | 194223 | 272366 |
| family_excluded | ORGANISM | 196121 | 272366 |

The [retained archive](../../benchmark/v3/results/context-comparison-audit-v1/README.md)
contains original scores, model files, input references, production audit summaries,
the independent file ledger and report, executed Python sources, job records and
Sonnet review with adjudication. Its portable verifier checks retained hashes,
reference membership, thresholds, original metric arithmetic and all 36 expanded
EC panels. It does not rerun remote models, rehash the remote 87 GB, or reconstruct
rankings from raw production TSVs locally. Those TSVs and the immutable execution
package remain at the registered remote locations.

## Remaining experimental requirements

This comparison has neither a CLEAN arm nor a matched no-gap arm. There are only
two gap-linked training rows per regime, both unconfirmed, and no gap-linked test
rows in the family-excluded table. The null layer remains untrained. Consequently
these results cannot establish missing-enzyme recovery, correct null assignments,
biological calibration or superiority over CLEAN.

The next discriminating experiment must evaluate candidate generation separately
from scoring, retaining all experimentally characterized proteins and functions.
Freeze predictor/reference versions and family/genome splits, compare sequence
baselines (including CLEAN) with matched context additions on the same proteins,
and report candidate recall, rank changes, calibrated calls, omissions and newly
introduced errors. Training must provide meaningful positive and negative/null
examples for the context tasks. A learned reranker cannot recover a function
that its candidate generator never proposes.

The panel is also concentrated: 12 assertions concern EC 1.5.98.3 in the same
genome. Eight of those proteins instead have EC 7.1.1.2 among their emitted
candidates. These are distinct activities: the official definitions specify
[F420-dependent methanophenazine reduction](https://iubmb.qmul.ac.uk/enzyme/EC1/5/98/3.html)
and [NADH-dependent ubiquinone reduction](https://iubmb.qmul.ac.uk/enzyme/EC7/1/1/2.html),
respectively. The [curated EC entry](https://enzyme.expasy.org/EC/1.5.98.3)
associates the F420 activity with multiple Fpo subunits. These observations identify
a useful context/subunit diagnostic, not twelve independent enzyme discoveries
or proof that every isolated subunit catalyses the complete reaction. Before
using it as function-recovery truth, verify the primary assays, component roles
and which evidence may enter training or inference. No EC alias remapping or
post-hoc label change has been made in the reported experiment.

Point 12 additionally needs held-out
condition-dependent phenotypes and function truth, beyond known-model feasibility
checks. These remain acceptance requirements, not gains demonstrated here.
