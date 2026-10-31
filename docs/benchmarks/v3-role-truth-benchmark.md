# Within-locus role benchmark (paralog and subunit discrimination)

Status 2026-10-28: registry (11 cases, 8 accepted; 3 negative-only loci), builder and evaluator exist
and are tested; four integration arms and a simple baseline evaluated on the seven scorable cases; the learned-context arm on the three
development-examined genomes has its three SEQONLY arms evaluated (2026-10-29, below) and its three context arms
being rerun (`tasks/v3-supervision-plan-2026-10-28.md`, WP2d).
Nothing here is a release claim. This document is the evidence trail for acceptance items 7
(retention and roles) and 11 (agent versus simple methods).

## Why a separate truth set

Closed-world labels call every unconfirmed candidate negative. For paralog loci that is exactly the
question under test: which of several homologous proteins carries the activity. The registry
therefore holds only loci whose candidate universe is complete for the function and whose inactive
members are negative by published evidence (inactive paralog, regulator, editing subunit,
different substrate), never by absence of annotation. Each case records the primary reference
(PMID/PMC), the reviewer and the evidence quoted from the archived text.

## Registry (`benchmark/v3/role_truth_cases.jsonl`)

| case | genome | function | active | evidence negatives | status |
|---|---|---|---|---|---|
| PMID15240874-threonylation | S. solfataricus P2 (cohort, train) | GO:0004829 | AAK42622.1 ThrRS catalytic | AAK40712.1 freestanding editing subunit (PMID 15240874) | accepted |
| speh-adometdc-sso | S. solfataricus P2 (cohort, train) | EC 4.1.1.50 | AAK40898.1 SpeH | AAK40855.1 pyruvoyl arginine decarboxylase, "ArgDC had no AdoMetDC activity" (PMC2533785) | accepted |
| coaa-coax-pantothenate-kinase-mtb | M. tuberculosis H37Rv | GO:0004594 | P9WPA7 CoaA | P9WPA1 CoaX, "did not detect hydrolysis of ATP by CoaX" (PMC11626117) | accepted |
| guab-imp-dehydrogenase-mtb | M. tuberculosis H37Rv | GO:0003938 | P9WKI7 GuaB2 | P9WKI3 GuaB1 (GMP reductase, PMC9790621), P9WKI5 GuaB3 (inactive, PMID 21081761) | accepted |
| guab-gmp-reductase-mtb | M. tuberculosis H37Rv | GO:0003920 | P9WKI3 GuaB1 | P9WKI7 GuaB2 | pending (no direct GuaB2 GMPR assay located; PMC9790621 states only the GuaB1 positive) |
| folp-dihydropteroate-synthase-mtb | M. tuberculosis H37Rv | GO:0004156 | P9WND1 FolP1 | P9WNC9 FolP2 "lacks enzymatic DHPS activity" (PMID 18680522) | accepted |
| ipms-dre-tim-family-saci | S. acidocaldarius (cohort) | EC 2.3.3.13 | AAY80302.1 LeuA | AAY80642.1 homocitrate synthase; AAY81610.1 CimA unknown (PMID 31330039) | pending |
| ccmm-carbonic-anhydrase-syn7942 | S. elongatus PCC 7942 (cohort) | GO:0004089 | CcaA | CcmM; three members unknown | pending (universe incomplete) |
| pgm-yhfr-gpmi-bsubtilis | B. subtilis 168 | GO:0004619 | P39773 GpmI | O07617 YhfR/PhoE "little if any PGM activity … no detectable PGM activity" (PMID 10869096) | accepted (drafted by `draft_role_cases.py`) |
| transaldolase-talb-fsa-ecoli | E. coli K-12 | GO:0004801 | P0A870 TalB | P32669 FsaA, P78055 FsaB/TalC "fructose-6-phosphate aldolase activity and not a transaldolase" (PMID 11120740); TalA unknown | accepted (drafted) |
| fgd1-fgd2-f420-g6pdh-mtb | M. tuberculosis H37Rv | GO:0052749 | P9WNE1 Fgd1 (EC 1.1.98.2 ECO:0000269 via ec2go) | P96809 Rv0132c "does not have the annotated FGD activity" (PMC3478262); two unknown | accepted (drafted) |

Negative-only loci (`role_truth_negative_only.jsonl`): S. solfataricus GO:0004089 (AAK41461.1 inactive
carbonic anhydrase; second member unknown), M. maripaludis GO:0004332 (CAF30242.1 no FBP-aldolase
activity) and the M. tuberculosis cutinase-like null group GO:0050525 (Culp1–5, all curated negative,
pending until the primary texts are read). They test abstention, not ranking.

Model-organism cases are a *disclosed development-examined* set: E. coli K-12, B. subtilis 168 and
M. tuberculosis H37Rv were excluded from the training cohort because earlier benchmark panels on
them were inspected during development. They were never used to fit any model.

### Candidate queue from Swiss-Prot negative statements

`mine_negative_function_statements.py` finds 1,324 curated negative statements (526 activity
negatives on 492 proteins) and 227 GO `NOT` annotations across all prokaryotic Swiss-Prot entries;
`propose_negated_functions.py` maps them to GO MF / EC identifiers and joins them to the candidate
universes of a claims file; `draft_role_cases.py` turns complete universes with experimental
positives into pending registry cases. First real run on the legacy claims of the three model
genomes: 2,349 coverage rows, 10 usable statements, 5 drafts (B. subtilis phosphoglycerate mutase
and E. coli transaldolase loci accepted after reading the references; M. tuberculosis cutinase-like
null group and Fgd1/Fgd2 pending; one draft rejected because a broad synonym lifted a specific
negation to a parent term — fixed in the matcher). 951 negative statements concern proteins that no
source proposes for the negated function: correct abstentions by construction, not evidence.

## Results on the enlarged registry (seven scorable cases, three genomes)

Archive `benchmark/v3/results/role-truth-v1/evaluations-registry-v2/` (registry sha256 at the time in
each `summary.json`). The three drafted cases are decided by sequence evidence in every arm (margins
0.4–0.7), so they move the averages but not the conclusion:

| arm | scored | strict top-1 | lenient top-1 | MRR | evidence negatives called | negative-only members called |
|---|---|---|---|---|---|---|
| legacy April-2026 integration | 7 of 11 | 5 | 5 | 0.83 | 4 of 13 | 0 of 7 |
| JAR 10526, priors none / + layout / + operons (identical) | 7 of 11 | 4 | 6 | 0.76 | 4 of 13 | 0 of 7 |
| simple method: best single claim score | 7 of 11 | 4 | 7 | 0.76 | 4 of 13 | 0 of 7 |

## Results: four integration arms on the M. tuberculosis cases (first registry)

Evaluator: `evaluate_role_truth.py`, strict top-1 = active member strictly above every evidence
negative; lenient accepts exact ties. Archive: `benchmark/v3/results/role-truth-v1/evaluations/`.

| arm | scored cases | strict top-1 | lenient top-1 | MRR | negatives called at the active's score | negative-only loci called |
|---|---|---|---|---|---|---|
| legacy April-2026 integration (9 sources, additive priors) | 4 | 2 | 2 | 0.71 | 4 of 9 | 0 of 2 |
| JAR 10526, legacy claims, priors none | 4 | 1 | 3 | 0.58 | 4 of 9 | 0 of 2 |
| + genome layout | 4 | 1 | 3 | 0.58 | 4 of 9 | 0 of 2 |
| + legacy operons, default priors | 4 | 1 | 3 | 0.58 | 4 of 9 | 0 of 2 |
| simple method: best single claim score per candidate (legacy claims; DIAMOND identity, Pfam = 1.0), no integration | 4 | 1 | 4 | 0.58 | 4 of 9 | 0 of 2 |

Per case (margin = active score minus best evidence-negative score):

| case | legacy | JAR 10526 arms |
|---|---|---|
| CoaA vs CoaX | −0.026 (CoaX first) | −0.026 (CoaX first) |
| GuaB2 vs GuaB1/GuaB3 | −0.018 (GuaB1 first, genomic-context bonus on the wrong paralog) | 0.000 (three-way tie) |
| FolP1 vs FolP2 | +0.013 | 0.000 (tie) |
| GuaB1 GMPR (pending) | +0.403 | +0.385 |

Reading: sequence evidence alone cannot separate these paralogs; the current integration leaves
them tied, and neither the legacy additive priors nor the current prior set with layout and operons
breaks the ties in the right direction (the legacy genomic-context bonus broke one in the wrong
direction). The margins are two orders of magnitude smaller than the GMPR case, where sequence
evidence is decisive. Both the ThrRS and SpeH cohort cases are complete evidence-labelled groups
in the 10406 exports and will be reported from that job's production replay.

## Learned-context arm, 21-genome joint integration (jobs 10658 and 10666, 2026-10-29)

The 18 cohort genomes and the three model genomes (M. tuberculosis H37Rv, E. coli K-12, B. subtilis 168) were
integrated jointly with JAR 10402 and the `diamond_homology` models of job 10406. The model-genome claims come
from the role-fenced Swiss-Prot DIAMOND reference (cohort fence plus the Bacillus, Escherichia, Mycobacterium and
Shigella exclusions); SEQONLY means no genome context. Evaluator: `evaluate_role_truth.py --format integrated_tsv
--score-field score`, threshold 0.5, registry sha256 fa8e916a…, negative-only sha256 f2225d73…, on the concatenated
tables of the three model genomes. Archive: `benchmark/v3/results/role-truth-v1/evaluations-learned-arm-v1/`
(`SEQONLY-none` = fenced DIAMOND claims only, no learned model; `SEQONLY-logistic` and `SEQONLY-gbdt` = the 10406
SEQONLY models).

| arm | cases scored of 11 | strict top-1 | lenient top-1 | MRR | mean margin | actives ≥ 0.5 | negatives ≥ 0.5 | negative-only loci called |
|---|---|---|---|---|---|---|---|---|
| SEQONLY-none | 4 | 4 | 4 | 1.00 | 0.190 | 1 of 11 | 0 of 13 | 0 of 7 |
| SEQONLY-logistic | 4 | 4 | 4 | 1.00 | 0.190 | 1 of 11 | 0 of 13 | 0 of 7 |
| SEQONLY-gbdt | 4 | 4 | 4 | 1.00 | 0.190 | 1 of 11 | 0 of 13 | 0 of 7 |

Top-1 counts are over the scored cases; the `scores_sha256` of the three evaluations is identical. Per-case scores
(the registry marks GuaB2, FolP1, B. subtilis GpmI and E. coli TalB as the active members; the ThrRS pair comes from
the evaluation that adds the Saccharolobus cohort table GCA_000007005.1):

| case | member | SEQONLY-none | SEQONLY-logistic | SEQONLY-gbdt |
|---|---|---|---|---|
| GuaB (M. tuberculosis) | GuaB2 (active) | 0.2626 | 0.2626 | 0.2626 |
| | GuaB1 | 0.2066 | 0.2066 | 0.2066 |
| | GuaB3 | 0.1316 | 0.1316 | 0.1316 |
| FolP (M. tuberculosis) | FolP1 (active) | 0.2499 | 0.2499 | 0.2499 |
| | FolP2 | 0.1836 | 0.1836 | 0.1836 |
| PGM (B. subtilis) | GpmI (active) | 0.5400 | 0.5400 | 0.5400 |
| | YhfR | 0.1828 | 0.1828 | 0.1828 |
| transaldolase (E. coli) | TalB (active) | 0.4180 | 0.4180 | 0.4180 |
| | TalC/FsaB (P32669) | 0.1372 | 0.1372 | 0.1372 |
| | Fsa (P78055) | not scored | not scored | not scored |
| ThrRS (S. solfataricus; **in-cohort, not held out**) | catalytic subunit (active) | 0.2358 | 0.1536 | 0.1117 |
| | editing subunit | 0.1831 | 0.0966 | 0.1168 |

TalA (P0A867) has unknown truth in the registry and is not a negative; it does not enter any margin. The remaining
registry cases are not scored on the model genomes: the CoaA/CoaX, Fgd1/Fgd2 and GuaB1 GMP-reductase cases because
the fenced Swiss-Prot DIAMOND reference proposes none of their members, and the cohort-genome cases because their
genomes are not among the model-genome tables.

The central finding is an abstention. On every model-genome candidate the learned models return the claim
unchanged and unsupported: the role-fenced reference has a different predictor-contract `reference_sha256`
(bc9802f1…) from the reference the models were trained against (2ef431fe…, the cohort fence), and the
predictor-contract gate rejects an unseen reference hash by design; `check_predictor_support.py` tests exactly this.
The two learned arms therefore reproduce the unscored projection on the model genomes, which is why their scores
and summaries equal `SEQONLY-none` exactly. The gate behaves as specified. The arm as set up cannot show whether the
learned models help on the model genomes, and these rows are no evidence either way. What they do establish is the
baseline a working learned arm has to beat: sequence evidence alone ranks the active member first in every scored
case, mostly with small margins, and only GpmI crosses the call threshold.

The only rescored cases are the in-cohort ThrRS pair, which are not held out from the models and say nothing about
generalisation. The logistic model keeps the order of the unscored projection at lower scores; the gbdt SEQONLY
model inverts it and ranks the editing subunit above the catalytic subunit.

The proposed remedy for the abstention is to retrain the `diamond_homology` arm on cohort claims computed against
the role-fenced reference, so that all 21 genomes share one predictor contract; the decision is pending. The
alternatives were rejected: using the same reference for the model genomes leaks self-hits, and loosening the gate for
subset references changes the JVM and contradicts the tested abstention contract.

### Context arms (job 10666)

The three context arms (UNIVERSAL, PATHWAY, ORGANISM) of job 10658 failed on a layout wiring defect in
`prepare_role_genome_context.py plan` (the model-genome rows pointed at `role_genomes.tsv` instead of their genome
layouts, and the context step left `unplacedProteins` empty for proteome entries without a RefSeq coordinate). After
the fix they were rerun as job 10666 (`arm-v2.sbatch`, plans `plans-v2`, JAR 10402, the same 10406 models); all twelve
arms exited 0. The nine context arms were evaluated with the same evaluator, registry and threshold and are archived
beside the SEQONLY arms. The first table adds the arm summaries (both evaluations: the three model genomes, and with
the in-cohort Saccharolobus table); the second the candidate scores of the second evaluation. Archive README and
`SHA256SUMS` cover all twelve folders.

| folder | learned model | candidates rescored | scored cases, model genomes | scored cases, plus SSO | strict top-1 (plus SSO) | mean margin, model genomes | mean margin, plus SSO | actives >= 0.5 | negatives >= 0.5 |
|---|---|---|---|---|---|---|---|---|---|
| `SEQONLY-none` | no learned model | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.190 | 0.163 | 1 of 11 | 0 of 13 |
| `SEQONLY-logistic` | 10406 logistic | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.190 | 0.163 | 1 of 11 | 0 of 13 |
| `SEQONLY-gbdt` | 10406 gbdt | not determined | 4 of 11 | 5 of 11 | 4 of 5 | 0.190 | 0.151 | 1 of 11 | 0 of 13 |
| `UNIVERSAL-none` | no learned model | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.175 | 0.150 | 2 of 11 | 1 of 13 |
| `UNIVERSAL-logistic` | 10406 logistic | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.175 | 0.153 | 2 of 11 | 1 of 13 |
| `UNIVERSAL-gbdt` | 10406 gbdt | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.175 | 0.151 | 2 of 11 | 1 of 13 |
| `PATHWAY-none` | no learned model | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.175 | 0.150 | 2 of 11 | 1 of 13 |
| `PATHWAY-logistic` | 10406 logistic | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.175 | 0.153 | 2 of 11 | 1 of 13 |
| `PATHWAY-gbdt` | 10406 gbdt | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.175 | 0.152 | 2 of 11 | 1 of 13 |
| `ORGANISM-none` | no learned model | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.175 | 0.150 | 2 of 11 | 1 of 13 |
| `ORGANISM-logistic` | 10406 logistic | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.175 | 0.153 | 2 of 11 | 1 of 13 |
| `ORGANISM-gbdt` | 10406 gbdt | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.175 | 0.149 | 2 of 11 | 1 of 13 |

| case | member | SEQONLY-none | SEQONLY-logistic | SEQONLY-gbdt | UNIVERSAL-none | UNIVERSAL-logistic | UNIVERSAL-gbdt | PATHWAY-none | PATHWAY-logistic | PATHWAY-gbdt | ORGANISM-none | ORGANISM-logistic | ORGANISM-gbdt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GuaB (M. tuberculosis) | GuaB2 (active) | 0.2626 | 0.2626 | 0.2626 | 0.2626 | 0.2626 | 0.2626 | 0.2626 | 0.2626 | 0.2626 | 0.2626 | 0.2626 | 0.2626 |
|  | GuaB1 | 0.2066 | 0.2066 | 0.2066 | 0.2066 | 0.2066 | 0.2066 | 0.2066 | 0.2066 | 0.2066 | 0.2066 | 0.2066 | 0.2066 |
|  | GuaB3 | 0.1316 | 0.1316 | 0.1316 | 0.1316 | 0.1316 | 0.1316 | 0.1316 | 0.1316 | 0.1316 | 0.1316 | 0.1316 | 0.1316 |
| FolP (M. tuberculosis) | FolP1 (active) | 0.2499 | 0.2499 | 0.2499 | 0.5109 | 0.5109 | 0.5109 | 0.5109 | 0.5109 | 0.5109 | 0.5109 | 0.5109 | 0.5109 |
|  | FolP2 | 0.1836 | 0.1836 | 0.1836 | 0.1836 | 0.1836 | 0.1836 | 0.1836 | 0.1836 | 0.1836 | 0.1836 | 0.1836 | 0.1836 |
| PGM (B. subtilis) | GpmI (active) | 0.5400 | 0.5400 | 0.5400 | 0.5400 | 0.5400 | 0.5400 | 0.5400 | 0.5400 | 0.5400 | 0.5400 | 0.5400 | 0.5400 |
|  | YhfR (negative) | 0.1828 | 0.1828 | 0.1828 | 0.5047 | 0.5047 | 0.5047 | 0.5047 | 0.5047 | 0.5047 | 0.5047 | 0.5047 | 0.5047 |
| transaldolase (E. coli) | TalB (active) | 0.4180 | 0.4180 | 0.4180 | 0.4180 | 0.4180 | 0.4180 | 0.4180 | 0.4180 | 0.4180 | 0.4180 | 0.4180 | 0.4180 |
|  | TalC/FsaB | 0.1372 | 0.1372 | 0.1372 | 0.1372 | 0.1372 | 0.1372 | 0.1372 | 0.1372 | 0.1372 | 0.1372 | 0.1372 | 0.1372 |
| ThrRS (S. solfataricus; in-cohort) | catalytic subunit (active) | 0.2358 | 0.1536 | 0.1117 | 0.2358 | 0.1624 | 0.1295 | 0.2358 | 0.1616 | 0.1285 | 0.2358 | 0.1566 | 0.1259 |
|  | editing subunit | 0.1831 | 0.0966 | 0.1168 | 0.1831 | 0.0967 | 0.0732 | 0.1831 | 0.0952 | 0.0693 | 0.1831 | 0.0934 | 0.0784 |

The evaluation folders do not record how many model-genome candidates each learned model rescored, so that column is
not determined from the archive; the identical scores below are what the abstention predicts.

What the arms show, and what they do not:

1. On the three model genomes the candidate scores are identical across none, logistic and gbdt within every context
   (the `scores_sha256` of the model-genome evaluations is the same for all nine context arms). The learned models
   abstain on every model-genome candidate, the predictor-contract fence already seen for SEQONLY. The context arms
   therefore carry no evidence of a context benefit on held-out role truth: no learned model scored a held-out locus,
   and nothing in this job separates the context kinds UNIVERSAL, PATHWAY and ORGANISM, whose model-genome results
   are identical to one another.
2. The only differences from SEQONLY come from the unscored prior and rescue path, which acts independently of the
   learned model, and they are the same in UNIVERSAL, PATHWAY and ORGANISM: FolP1 (P9WND1) rises above the call
   threshold and is called, and the YhfR negative (O07617) also rises above it and is wrongly called. Every other case
   is unchanged. The first is a correct call of an active member; the second is a false positive that narrows the
   GpmI margin without changing the strict top-1. The net effect on the called counts is one more active and one more
   negative; this is a property of the prior and rescue path, not of the learned models, and the registry is too small
   to say which direction generalises.
3. The in-cohort ThrRS pair is the only locus the learned models rescore. Its two scores differ by model kind (and
   slightly between contexts), and the catalytic subunit stays ranked above the editing subunit in every context arm,
   including the gbdt arms. The SEQONLY gbdt inversion is therefore not reproduced with context. This pair is not
   held out from the models, so the order says nothing about generalisation.

Not reported: any learned-model score on a held-out model-genome candidate, and any claim that a context improves
role-level discrimination.

## Role-fenced retraining (WP2e) and the 21-genome arms under its models (jobs 15690304 and 15690691, 2026-10-29)

The abstention of the 10406 models was removed by retraining the `diamond_homology` arm on cohort claims computed
against the role-fenced Swiss-Prot reference (arm `diamond_homology_rolefence`; Shaheen III job 15690304, JAR 10402,
four contexts times logistic and gbdt), so that the cohort and model genomes share one predictor contract. The twelve
21-genome arms were then rerun with these models (job 15690691, plans `context/plans-rolefence`) and evaluated as
before (registry sha256 fa8e916a…, negative-only sha256 f2225d73…, threshold 0.5). Archive:
`benchmark/v3/results/role-truth-v1/evaluations-rolefence-arm-v1/` (arm folders, receipts, training receipts,
`SHA256SUMS`).

The supervision universe is not that of job 10406. The disclosed comparison is identical in all four contexts:

| context | shared examples | original_only (10406 universe) | new_only (role-fenced) | representation changes |
|---|---|---|---|---|
| SEQONLY | 3,424 | 838 | 366 | 423 |
| UNIVERSAL | 3,424 | 838 | 366 | 423 |
| PATHWAY | 3,424 | 838 | 366 | 423 |
| ORGANISM | 3,424 | 838 | 366 | 423 |

The models' own held-out evaluation on the cohort test split (supported rows only; the row set differs slightly
between contexts, so the rows are not paired):

| context | kind | AUROC | log-loss | supported rows | total rows | test positive rate | n |
|---|---|---|---|---|---|---|---|
| SEQONLY | logistic | 0.740 | 0.419 | 310 | 429 | 0.145 | 310 |
| SEQONLY | gbdt | 0.738 | 0.433 | 310 | 429 | 0.145 | 310 |
| UNIVERSAL | logistic | 0.776 | 0.404 | 309 | 429 | 0.146 | 309 |
| UNIVERSAL | gbdt | 0.669 | 0.403 | 309 | 429 | 0.146 | 309 |
| PATHWAY | logistic | 0.772 | 0.428 | 308 | 429 | 0.146 | 308 |
| PATHWAY | gbdt | 0.708 | 0.395 | 308 | 429 | 0.146 | 308 |
| ORGANISM | logistic | 0.785 | 0.430 | 305 | 429 | 0.144 | 305 |
| ORGANISM | gbdt | 0.663 | 0.407 | 305 | 429 | 0.144 | 305 |

The role-truth summary per arm. Strict top-1 is over the scored cases; the first evaluation is the three model
genomes, the second adds the in-cohort Saccharolobus table (ThrRS pair, not held out):

| arm | scored, model genomes | strict top-1, model genomes | mean margin, model genomes | scored, plus SSO | strict top-1, plus SSO | mean margin, plus SSO | actives >= 0.5 (plus SSO) | case negatives >= 0.5 (plus SSO) |
|---|---|---|---|---|---|---|---|---|
| `SEQONLY-none` | 4 of 11 | 4 of 4 | 0.190 | 5 of 11 | 5 of 5 | 0.162 | 1 of 11 | 0 |
| `SEQONLY-logistic` | 4 of 11 | 4 of 4 | 0.356 | 5 of 11 | 5 of 5 | 0.298 | 2 of 11 | 0 |
| `SEQONLY-gbdt` | 4 of 11 | 4 of 4 | 0.255 | 5 of 11 | 5 of 5 | 0.211 | 2 of 11 | 0 |
| `UNIVERSAL-none` | 4 of 11 | 4 of 4 | 0.175 | 5 of 11 | 5 of 5 | 0.150 | 2 of 11 | 1 |
| `UNIVERSAL-logistic` | 4 of 11 | 4 of 4 | 0.304 | 5 of 11 | 5 of 5 | 0.259 | 1 of 11 | 0 |
| `UNIVERSAL-gbdt` | 4 of 11 | 2 of 4 | 0.107 | 5 of 11 | 3 of 5 | 0.094 | 0 of 11 | 0 |
| `PATHWAY-none` | 4 of 11 | 4 of 4 | 0.175 | 5 of 11 | 5 of 5 | 0.150 | 2 of 11 | 1 |
| `PATHWAY-logistic` | 4 of 11 | 4 of 4 | 0.292 | 5 of 11 | 5 of 5 | 0.248 | 1 of 11 | 0 |
| `PATHWAY-gbdt` | 4 of 11 | 2 of 4 | 0.134 | 5 of 11 | 3 of 5 | 0.118 | 0 of 11 | 0 |
| `ORGANISM-none` | 4 of 11 | 4 of 4 | 0.175 | 5 of 11 | 5 of 5 | 0.150 | 2 of 11 | 1 |
| `ORGANISM-logistic` | 4 of 11 | 4 of 4 | 0.294 | 5 of 11 | 5 of 5 | 0.250 | 1 of 11 | 0 |
| `ORGANISM-gbdt` | 4 of 11 | 2 of 4 | 0.097 | 5 of 11 | 3 of 5 | 0.088 | 0 of 11 | 0 |

Learned-prior coverage of the model-genome integrated tables (counted from the integrated tables, which are not
archived; a row counts when `priors_fired` contains a `learned.` entry):

| context | rows in the three model-genome tables | rows carrying a learned prior, logistic | rows carrying a learned prior, gbdt |
|---|---|---|---|
| SEQONLY | 101,705 | 67,291 | 67,291 |
| UNIVERSAL | 101,705 | 62,444 | 62,444 |
| PATHWAY | 101,705 | 62,088 | 62,088 |
| ORGANISM | 101,705 | 61,649 | 61,649 |

Per-case candidate scores of the second evaluation (GuaB2, FolP1, GpmI, TalB and the ThrRS catalytic subunit are the
registry's active members; YhfR is a registry negative; Fsa is not scored in any arm; `n/s` = not scored):

| case | member | SEQONLY-none | SEQONLY-logistic | SEQONLY-gbdt | UNIVERSAL-none | UNIVERSAL-logistic | UNIVERSAL-gbdt | PATHWAY-none | PATHWAY-logistic | PATHWAY-gbdt | ORGANISM-none | ORGANISM-logistic | ORGANISM-gbdt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PMID15240874-threonylation | AAK40712.1 | 0.1831 | 0.1126 | 0.1211 | 0.1831 | 0.0823 | 0.0777 | 0.1831 | 0.0814 | 0.0742 | 0.1831 | 0.0808 | 0.0712 |
|  | AAK42622.1 (active) | 0.2320 | 0.1809 | 0.1539 | 0.2320 | 0.1615 | 0.1223 | 0.2320 | 0.1567 | 0.1298 | 0.2320 | 0.1564 | 0.1227 |
| guab-imp-dehydrogenase-mtb | P9WKI3 | 0.2066 | 0.1433 | 0.1071 | 0.2066 | 0.1302 | 0.0994 | 0.2066 | 0.1271 | 0.0962 | 0.2066 | 0.1266 | 0.0992 |
|  | P9WKI5 | 0.1316 | 0.0592 | 0.0901 | 0.1316 | 0.0446 | 0.0773 | 0.1316 | 0.0452 | 0.0678 | 0.1316 | 0.0445 | 0.0737 |
|  | P9WKI7 (active) | 0.2626 | 0.2325 | 0.1365 | 0.2626 | 0.1877 | 0.0790 | 0.2626 | 0.1817 | 0.0911 | 0.2626 | 0.1812 | 0.0886 |
| folp-dihydropteroate-synthase-mtb | P9WNC9 | 0.1836 | 0.1132 | 0.1211 | 0.1836 | 0.1000 | 0.0994 | 0.1836 | 0.0984 | 0.0962 | 0.1836 | 0.0978 | 0.0992 |
|  | P9WND1 (active) | 0.2499 | 0.2102 | 0.1539 | 0.5109 | 0.1598 | 0.0706 | 0.5109 | 0.1554 | 0.0847 | 0.5109 | 0.1546 | 0.0773 |
| pgm-yhfr-gpmi-bsubtilis | O07617 | 0.1828 | 0.1122 | 0.1211 | 0.5047 | 0.1090 | 0.0916 | 0.5047 | 0.1070 | 0.0917 | 0.5047 | 0.1064 | 0.0919 |
|  | P39773 (active) | 0.5400 | 0.8392 | 0.5747 | 0.5400 | 0.7728 | 0.3487 | 0.5400 | 0.7487 | 0.3963 | 0.5400 | 0.7538 | 0.3170 |
| transaldolase-talb-fsa-ecoli | P0A870 (active) | 0.4180 | 0.5735 | 0.6199 | 0.4180 | 0.4850 | 0.2830 | 0.4180 | 0.4625 | 0.3174 | 0.4180 | 0.4664 | 0.2601 |
|  | P32669 | 0.1372 | 0.0641 | 0.1159 | 0.1372 | 0.0487 | 0.0639 | 0.1372 | 0.0489 | 0.0712 | 0.1372 | 0.0484 | 0.0653 |
|  | P78055 | n/s | n/s | n/s | n/s | n/s | n/s | n/s | n/s | n/s | n/s | n/s | n/s |

What the arms show, and what they do not:

1. The role-fenced models act on the model genomes. The scores of every scored candidate differ from the `none` arm
   of the same context, and the learned prior is present on most, not all, model-genome rows. The earlier archive
   could only show that the models abstained; this is the first arm in which learned integration acts on held-out role
   truth. The three model genomes are held out from the models, the ThrRS pair is not.
2. SEQONLY: both models keep the active member first in every scored case and widen the margins of PGM and
   transaldolase, while calling no negative. In the logistic arm the other three cases keep small margins; the gbdt
   arm narrows the GuaB and FolP margins relative to logistic. The logistic and gbdt models also push TalB over the
   threshold.
3. Context arms, logistic: the unscored prior path of the `none` arms called the YhfR negative and rescued FolP1;
   the logistic model removes the YhfR call, and it also pushes the correct FolP1 rescue back under the threshold.
   Strict top-1 is kept, and the margins of PGM and transaldolase are wide, but below the SEQONLY logistic margins in
   all three contexts. On held-out role truth the context does not help the logistic model; the registry cannot
   separate UNIVERSAL, PATHWAY and ORGANISM from one another.
4. Context arms, gbdt: under every learned context the gbdt model ranks GuaB1 above
   GuaB2 and FolP2 above FolP1, so strict top-1 drops on the held-out cases, and the mean margin falls below that of
   the arms without a learned model. Under context, gbdt is worse than no model on this registry; SEQONLY gbdt is not.
5. The cohort-test evaluation points the same way on discrimination: with context the logistic AUROC is higher than
   SEQONLY and the gbdt AUROC lower. Log-loss does not follow this pattern (gbdt log-loss is lower under context, and
   logistic log-loss is mixed), the row sets differ, and the test split is small, so the agreement is suggestive of
   (3) and (4) and nothing more.

Limits. The registry has 11 cases, of which five are scorable on these genomes; only one of the five is in-cohort
(ThrRS), and the held-out cases come from three genomes. Every difference above rests on a handful of loci, and no
interval or significance statement is possible. The supervision universe differs from that of job 10406, so these
models are not the 10406 models retrained on the same data. No claim of improved role-level discrimination, or of a
context benefit, follows from these arms.

## Open arms

1. Learned-context arm: done for the role-fenced models (WP2e, above). The twelve arms score the model genomes, and
   the registry shows no benefit of learned context on held-out role truth. What remains is the GuaB1 GMP-reductase
   case (pending in the registry), the cutinase cases, and a larger truth set: with five scorable cases no
   difference between learned arms can be established.
2. Registry growth from the drafter's second run (ec2go-backed positives, specificity guard).
3. Item 11 comparison against the learned-context arm once it exists; the simple-method row above is the first half of that comparison (all four loci tie at the Pfam ceiling of 1.0, so the simple method never ranks a paralog first or second: strict 1/4, lenient 4/4).
