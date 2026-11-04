# Frozen plan: positive-unlabelled training of the learned stage and of the no-candidate outcome (row 4), 2026-11-02

Status: frozen 2026-11-02 before any fit or score was read under it. The text is copied to
`benchmark/v3/results/pu-learning-v1/plan_frozen.md` and its sha256 recorded in
`benchmark/v3/results/pu-learning-v1/PLAN_SHA256`; later changes go under "Deviations", dated, append-only;
the body above that heading is not edited.

Direction (Robert Hoehndorf, 2026-11-02): "check on PU-GO and positive unlabelled learning; also, there are some
confirmed NOT annotations in GOA, so these can potentially be used."

Approval (Robert Hoehndorf, 2026-11-02, relayed by the supervising session): the plan is approved in full,
unconfirmed rows treated as unlabelled, GOA experimental NOT rows and accepted Swiss-Prot negations used for
evaluation only; "Lock in the plans, run implementation with agy, you monitor, control, evaluate." Execution
rules attached to the approval: new files and folders only (the existing trainer `scripts/train_context_model.py`
is imported, not edited, because the complete-locus runner of the other session imports it; no JVM change in
the shared tree); cohort fits on ws; production replays and the Price part on Shaheen (k10018, nodes
nid00022, nid00497, nid01033 excluded); the Price part waits for the other session's frozen contract claims;
other sessions' folders are read only; no commits; no /tmp.

Scope: acceptance row 4 ("Learned inference", Not met). The audit of 2026-11-01 names one missing element:
enough honest labels to train and test the no-candidate (null) outcome. The complete-locus plan of this
morning (`tasks/v3-row4-complete-locus-labels-2026-11-02.md`, frozen, registry frozen 13:36) measured how far
closed-world completeness can go: 941 complete competitions, but 0 all-negative ones, 6 contested usable ones,
0 of them in calibration or test; criterion C7 (null layer trained) failed. The blocker table of that registry
shows why: 84,843 cohort competitions, 231 with a positive, and nearly all of the rest blocked by one or more
unlabelled members. Completeness will not arrive by curation. This plan changes the learning problem instead:
unconfirmed candidates are treated as what they are, unlabelled, and the model is trained with
positive-unlabelled (PU) learning. Experimental negatives (GOA NOT, reviewed Swiss-Prot negations) are kept
for evaluation.

## 1. What PU-GO does (verified against the paper and the code)

Zhapa-Camacho F, Tang Z, Kulmanov M, Hoehndorf R. Predicting protein functions using positive-unlabeled
ranking with ontology-based priors. Bioinformatics 2024;40(Suppl 1):i401–i409, doi:10.1093/bioinformatics/btae237,
PMID 38940168, PMC11211813. Code: github.com/bio-ontology-research-group/PU-GO (`nn.py`, class `PUModel`);
data: Zenodo 10.5281/zenodo.11079885.

- Problem framing: every (protein, GO term) pair without an experimental annotation is unlabelled, not negative.
- Risk estimator: the non-negative PU risk of Kiryo et al. (2017), per GO term i, with a logistic loss
  l(z) = −ln σ(z): R_i = π_i R_i^P+ + max{0, R_i^U− − π_i R_i^P− + β}. The clamp is a plain `relu` in the code
  (no gradient-ascent step). β = γ π with a tuned margin factor γ.
- Ranking form (the one the paper and the README select, `pu_ranking_multi`): the unlabelled term is replaced by
  a pairwise ranking loss, R_i^U− = −1/(|P_i||U_i|) Σ_{x∈P_i} Σ_{y∈U_i} ln σ(g(x) − g(y)); only positives are
  required to outrank unlabelled samples, instead of unlabelled scores being pushed to zero. (In the released
  code the "ranking" term is implemented element-wise rather than pairwise; the paper states the pairwise form.)
- Class prior: not learned. π_i = π_0 · S_i / S_max, where S_i is the fraction of training proteins annotated
  with term i after true-path propagation and π_0 a tuned hyperparameter searched in [1e-4, 1e-3] (selected
  MF 1.43e-4, CC 1.11e-4, BP 7.85e-4). This is the "ontology-based prior": a frequency under propagation,
  scaled; terms near the root get larger priors.
- Model and inputs: ESM2-15B mean embeddings (5,120-d), a two-block MLP with one sigmoid output per term,
  descendants' maxima propagated to ancestors at prediction time; no calibration step; optional mixing with
  DIAMOND scores.
- Reported (Swiss-Prot 2023_03, DIAMOND similarity split, Table 2): Fmax MF 0.532 / BP 0.437 / CC 0.721 for
  PU-GO against 0.508 / 0.428 / 0.709 for the same MLP trained with unlabelled pairs as negatives (MLP(ESM2));
  Smin lower in all three; AUPR higher in MF and CC, lower in BP. The paper's own statement of the benefit:
  "the advantage of PU learning is directly demonstrated when comparing PU-GO to MLP(ESM2) which uses the same
  classifier function as PU-GO but considers unlabeled samples as negatives." Time-split results against
  TALE, SPROF-GO, NetGO3 and ATGO+ are in its Table 3 (not re-read here).
- Not in PU-GO: GO NOT annotations (the code has no qualifier handling); evaluation still scores unannotated
  pairs as negatives (Fmax, AUPR), as every CAFA-style evaluation does.

What transfers to GSPA: the risk estimator (nnPU with logistic loss), the ranking variant, the idea of a prior
tied to annotation frequency, and the paired comparison "same model, PN loss vs PU loss". What does not: the
per-term prior formula (GSPA scores candidate rows across thousands of functions with one model, so the prior
is per candidate row, not per term), the ESM input (GSPA's model reads the context-feature contract, not
sequence embeddings), and the absence of calibration (GSPA must emit probabilities that survive a reliability
check and a harm count).

## 2. The learning problem in GSPA, restated as PU

A row of the training table is one candidate (genome, protein, function, role) with 160 contract-3 features.
Today's label: 1 if Swiss-Prot has an experimental annotation of that function on that protein
(`positive_truth`), 0 if the protein carries experimental annotations but not this one
(`closed_world_unconfirmed`). The current logistic/tree models are fitted with those zeros as negatives
("reference support" target). Candidates on proteins without any experimental annotation (`no_truth_not_negative`,
472,019 rows in the homology export) are left out.

Feasibility counts measured on 2026-11-02 (ws, light Python; scripts in `~/tmp/gspa-v3-pu-design/`):

| quantity | homology table | family-excluded table | full family-excluded export |
|---|---:|---:|---:|
| candidate rows | 4,262 | 2,444 | 272,366 |
| positives (train / calibration / test) | 394 / 72 / 75 | 183 / 29 / 34 | same proteins |
| unlabelled rows on labelled proteins (train / cal / test) | 2,644 / 544 / 533 | 1,499 / 316 / 382 | 4,302 rows on labelled proteins in all |
| labelled proteins | 364 | 235 | 19,295 proteins in the export |
| competitions (genome, type, function, role) | 2,945 | 1,849 | 86,987 |
| competitions with ≥ 2 members | 684 | 396 | 38,415 |
| competitions with ≥ 2 members and a labelled positive (train / cal / test) | 85 / 9 / 17 | 44 / 4 / 7 | — |
| competitions with ≥ 2 members and no labelled positive (train / cal / test) | 398 / 75 / 100 | 240 / 34 / 67 | — |
| competitions with ≥ 1 labelled protein and ≥ 1 unlabelled member | — | — | 2,234 |

Class prior. The labelled rate among candidate rows on labelled proteins is 0.117–0.130 (homology) and
0.082–0.109 (family-excluded). An Elkan–Noto estimate (a group-aware out-of-fold logistic classifier of
labelled vs rest; c = mean score on labelled rows; π̂ = labelled rate / c) gives π̂ = 0.44 (train OOF), 0.28
(calibration), 0.47 (test) on the homology table and 0.51 / 0.36 / 0.48 on the family-excluded table, i.e. an
estimated 18–45 % of the "unconfirmed" rows would be true positives. The estimator is biased upward when the
classifier is weak (c is only 0.17–0.43 here) and when the "selected completely at random" assumption fails,
so these are upper-range values. The plan treats π as an uncertain quantity with a pre-registered sweep
(section 5), not as a fitted constant.

Experimental negatives available (GOA 2026-06-17 inventory `gaf-negative-audit-v1`, experimental codes only,
proteins only; Swiss-Prot negation statements from the complete-locus scan of 356,637 prokaryotic entries):

| stratum | GOA experimental NOT rows (proteins) | Swiss-Prot negation statements (all classes) | negative candidates joined to competitions (complete-locus registry) | competitions with a negative candidate |
|---|---:|---:|---:|---:|
| cohort, 18 genomes (19 taxids) | 12 (7) | 74 | 40 pairs | 8 (all in train genomes) |
| Price, 22 genomes (19 taxids) | 1 (1) | 12 | 1 pair | 0 |
| model genomes E. coli K-12, B. subtilis 168, M. tuberculosis H37Rv | 85 (73) | 464 | 137 pairs | 27 |
| all GOA, experimental, proteins | 4,095 rows in 171 taxa; 2,967 of them human, mouse, Arabidopsis, fly, zebrafish | 1,551 prokaryotic statements | — | — |

Of the 47 negative candidates in those 35 competitions, 10 were reviewed in the complete-locus plan and 7
accepted (`registry/negative_reviews.tsv`); the review queue holds 22 more. The three model genomes hold
three quarters of the usable negatives, and all 27 of their negative-bearing competitions.

Conclusion from the counts: training needs no negatives, and the negatives are too few to train on anyway;
they are enough for a harm count in evaluation if the model-genome stratum is scored (section 4.3).

## 3. Method

### 3.1 Candidate layer (marginal correctness): PU logistic regression (unbiased risk, nnPU clamp as a check), same artifact format

- Features, contract, normalisation, group weighting, five grouped folds, seed 42, family masking, splits:
  unchanged from the family-excluded recipe (`run_family_excluded_training.py`, trainer
  `scripts/train_context_model.py`). The table rows are identical to the current table, so every comparison
  with the current models is paired row for row.
- Positives P: rows with `label_basis = positive_truth`. Unlabelled U1: rows with
  `label_basis = closed_world_unconfirmed`. The single `explicit_negative_truth` row stays out of the fit
  (evaluation only, section 3.3).
- Objective (new `--objective pu` in the trainer; default stays `pn`): the unbiased PU risk of du Plessis et
  al. (2015) with the logistic loss and L2, R = π R_P^+ + R_U^− − π R_P^−, where R_P^± are the weighted mean
  losses of positives labelled 1 / labelled 0 and R_U^− the weighted mean loss of unlabelled rows labelled 0.
  With the logistic loss, l(z,1) − l(z,0) = −z, so R = −π · mean_P(z) + mean_U ln(1 + e^z) + L2: convex, and
  L-BFGS fits it exactly as the PN loss is fitted today. The non-negative correction of Kiryo et al. (2017)
  (the form PU-GO uses) exists to stop flexible models driving the estimated negative-class risk
  R_U^− − π R_P^− below zero; for a 160-weight linear model it is a check, not the fit: the trainer reports
  that quantity on train OOF and on calibration, and if it is negative the nnPU clamp (gradient-ascent step,
  β = 0) is applied and the artifact says so. The second arm is the PU-GO ranking form (`--objective pu_rank`):
  the unlabelled term is the pairwise logistic ranking loss between positive and unlabelled rows **of the same
  genome and function type** (so that a GO positive is never asked to outrank an EC unlabelled row), with
  π R_P^+ kept.
- Prior π: one scalar per fit, from the train split only, by the Elkan–Noto estimator on the OOF labelled-vs-rest
  margins of the PN model; its value and both components (labelled rate, c) are written into the artifact.
  Pre-registered sweep around it (section 5). No per-term prior: with 3,103 functions and 541 positives no
  per-function frequency is estimable, and the PU-GO formula would reduce to the propagation frequency of
  each GO term, which is already a feature (`fn.go_aspect_*`, depth is not a feature and is not added here).
- Output probability: the PU risk minimiser with the right π estimates P(y = 1 | x) directly. Platt scaling on
  the calibration split against PN labels would re-impose the closed-world target, so it is not refitted.
  Instead, on the calibration split, the Elkan–Noto identity P(labelled | x) = c · P(y = 1 | x) is checked in
  ten bins of p̂ (observed labelled rate against c · p̂); the only adjustment allowed is a single temperature
  on the margin, fitted on the calibration split to minimise that PU-adjusted calibration error, with the
  intercept fixed. The artifact keeps the Platt slot (`a` = temperature or 1, `b` = 0) so the JVM reader is
  unchanged, and records `probability_target = experimental_support_under_positive_unlabelled_policy`, π, c
  and the labelled rate.
- GBDT: the same PU risk fits the Newton-boosted trees. Per row: a positive contributes gradient −π/|P| (its
  loss is linear in z) and Hessian 0; an unlabelled row contributes gradient σ(z)/|U| and Hessian
  σ(z)(1 − σ(z))/|U|. All Hessians are non-negative, so the Newton step is as stable as today; the leaf
  regulariser λ = 1 keeps leaves with only positives finite. The negative-class risk check and the clamp apply
  as for the logistic arm.
- Rule arm: unchanged (`cohort_rlc` and the Stage 0/1 rules).

### 3.2 No-candidate (null) layer: the same softmax, trained by PU marginal likelihood

The production null layer (`locus_model`, unchanged JVM form) is a softmax over [null, member 1, …, member n]
with logits [b0 + w·g, s·m_1, …, s·m_n], where m_j are the candidate-layer OOF margins and g the six group
features. Today it is fitted only on competitions whose every member has a label, hence never.

Under PU every competition is usable. A competition is the set of table rows sharing (genome, function type,
function id, role), the same key the production stage uses to form a locus group; no completeness certificate
is needed because unlabelled members are allowed, and `locus_group_id` is not used. Let c be the label
frequency (probability that a true positive member carries an experimental label), estimated in 3.1 from the
train split. For a competition:

- if one member is labelled positive, the target is that member (as now);
- if two or more members are labelled positive, the competition is excluded (isoenzymes; as now) and counted;
- if no member is labelled, the observation is "no label", whose probability under the model is
  P(null) + (1 − c) Σ_j P(j). The trainer maximises this marginal likelihood (the same `locus_nll` with the
  no-label case added) over all competitions with ≥ 2 members of the train split; L2 chosen on the calibration
  split's marginal likelihood as now.

This identifies the null bias and weights because c < 1 is pinned from the candidate layer and not refitted
here (if c were free, P(null) and 1 − c would trade off). Counts available: train 398 no-label and 85
single-positive competitions with ≥ 2 members (homology), 240 and 44 (family-excluded), against the trainer's
default gates of ≥ 20 and ≥ 20; calibration 75/9 and 34/4 against gates ≥ 5 and ≥ 5 (the family-excluded
calibration split has 4 single-positive contested groups, one short of the gate; the gate is not changed,
and if it blocks, the family-excluded null layer is reported as blocked at the gate and the homology-regime
null layer carries the result). The PU-trained null layer is a model of "none of the proposed candidates is
experimentally supported, beyond what the label frequency explains"; its probability target is written into
the artifact in those words.

### 3.3 How each label kind is used, and why there is no circularity

| label kind | training (candidate layer) | training (null layer) | prior / c estimation | evaluation |
|---|---|---|---|---|
| experimental positive (Swiss-Prot GO experimental codes, EC with ECO:0000269; Price-148 on Price) | positive | positive member | yes (train split only) | positive |
| unconfirmed candidate on a labelled protein | unlabelled | unlabelled member | yes (train split only) | unlabelled: never a negative, never a false positive |
| candidate on an unlabelled protein | not used in v1 (optional U2 arm, section 5) | not used in v1 | no | not scored as truth |
| GOA experimental NOT (exact term, no propagation), reviewed and accepted Swiss-Prot activity negations, accepted role-truth negatives | **not used** | **not used** | **not used** | negative: harm count, all-negative competitions, specificity |

Reasons. (a) The negatives are 47 candidate rows in all strata, 7 accepted so far; a model trained on them
could not also be checked by them. (b) PU training is complete without them, so keeping them out gives an
independent check of the PU assumption: if the "selected completely at random" assumption and the prior are
badly wrong, the model will put high probability on experimentally denied functions and the harm count will
show it. (c) The complete-locus registry already defines the negatives with provenance, review and
dispositions; this plan reads that registry and adds nothing to it. A later plan may try positive-negative-
unlabelled training once the review queue is worked through; that is out of scope here.

## 4. Data, strata and splits

4.1 Cohort (primary). The 18 label-rich genomes, the frozen 11/3/4 split disjoint on genome, group, genus and
family, family masking at 30 %/80 %, both claim regimes (homology, family-excluded), the production feature
exports already archived (`family-excluded-training-v1`, `training-cohort`). Nothing is re-mined. The only new
rows are none: the tables are the existing ones with the loss changed.

4.2 Price (22 genomes). The contract-bearing DIAMOND claims being produced by the Price plan
(`tasks/v3-price-contract-learned-2026-11-02.md`, Shaheen job family gspa-price-contract) are reused. Once its
Step C claims are frozen, the PU artifacts are run in its `learned` arm form on the same cells (a second
artifact path, no code change), and its evaluator functions (`row1_calibration.py`) score them with the frozen
row-1 genus split (7 calibration / 15 test). Price-148 labels are positives; everything else unlabelled; the
one Price GOA NOT row is reported but cannot carry a metric.

4.3 Model-genome stratum (evaluation only, never trained on). E. coli K-12 MG1655 (GCF_000005845.2),
B. subtilis 168 (GCF_000009045.1), M. tuberculosis H37Rv (GCF_000195955.2) hold 85 GOA experimental NOT rows
and 137 joined negative pairs in 27 competitions. Their existing claims (`rolefence_diamond`) carry a
predictor contract the models do not accept. Under the family-excluded regime, however, the training
reference needs no re-fencing: every hit at ≥ 30 % identity and ≥ 80 % coverage is discarded, which removes
the proteins' own Swiss-Prot entries and close paralogues exactly as for the cohort genomes. So the stratum is
produced with the unchanged family-excluded adapter and reference (same predictor identity as the
family-excluded models), a feature export for three genomes, and the complete-locus registry's model-set
labels and accepted negatives. It is reported as a held-out, development-examined stratum (the registry's
own words), separate from cohort test. If the predictor-identity gate or the input-support gates (unseen
availability profile, out-of-range score) reject the claims or most candidates, the stratum is reported as
not scorable with the reasons and counts, and the harm count falls back to the 8 cohort competitions (train
genomes, so descriptive only) and the family-excluded regime's single explicit negative.

4.4 Family and genus holdouts. Unchanged frozen split. Added sensitivity, pre-registered: leave-one-genus-out
over the 18 cohort genomes for the PN and PU logistic models (18 refits each, seconds per fit), reporting the
spread of the PU–PN difference; no model is selected from it.

## 5. Arms and comparisons

| arm | candidate layer | null layer | regime |
|---|---|---|---|
| A0 rules | Stage 0, Stage 1, cohort_rlc | none | both |
| A1 PN logistic (current default) | existing frozen artifact | untrained | both |
| A2 PN tree (current) | existing frozen artifact | untrained | both |
| A3 PU logistic | new, π from train OOF | PU marginal likelihood | both |
| A4 PU-rank logistic | new, PU-GO ranking form | PU marginal likelihood | both |
| A5 PU tree | new | PU marginal likelihood | both (or "PN only" if 3.1 applies) |
| S1 prior sweep | A3 at π ∈ {labelled rate (= PN limit), 0.20, π̂_EN, min(2π̂_EN, 0.8)} | — | both |
| S2 unlabelled pool | A3 with U2 = a stratified sample of 10 unlabelled candidates per labelled-protein competition from the full export (`no_truth_not_negative` rows), π re-estimated | — | family-excluded |
| S3 leave-one-genus-out | A1 and A3 | — | family-excluded |

Everything in A3–A5 and S1–S3 is fitted before any test or Price score is read; artifacts are hashed and listed
in `SHA256SUMS` first.

## 6. Metrics

All on supported rows, with the same support gates as now; abstentions counted and reported.

Candidate layer, per split (calibration, test, Price test, model stratum):
1. Supported precision as now (labelled positives among calls at the frozen threshold): a lower bound under PU,
   reported with the PU-adjusted expected precision (mean p̂ among calls) and the labelled-fraction identity
   (observed labelled rate among calls ÷ c).
2. Recall over all positives of the split; positive-vs-unlabelled AUROC (biased for every arm in the same
   way; comparable across arms); per-positive rank among the unlabelled rows of its own competition (the
   PU-GO ranking view): fraction of positives ranked first, mean reciprocal rank.
3. Calibration under PU: reliability table of observed labelled rate vs c · p̂ in ten bins, its expected
   calibration error (PU-ECE), and the PN reliability for reference; Brier against labels after c scaling.
4. Harm and corrections against the PN logistic and against the rules, pairwise on the same rows:
   corrections = positives called only by the arm; harms = experimental negatives called only by the arm
   (model stratum, cohort negatives), and, as the closed-world figure the audit already reports, unconfirmed
   rows called only by the arm, labelled as such and never summed with the negatives.
5. Specificity on experimental negatives: fraction of the 47 (or however many accepted at freeze time) negative
   candidates with p̂ below the threshold, per arm.

Null layer:
6. Trained or not, with the gate counts; on calibration and test competitions with one labelled positive
   (homology: 9 and 17 contested; plus the 812 single-positive groups of the complete-locus registry, mostly
   singletons): P(null) should be low; NLL against the uniform baseline log(n + 1).
7. Deletion test (new, pre-registered): for every test competition with exactly one labelled positive and ≥ 2
   members, remove the positive member and rescore; the honest answer is "no known candidate is correct"
   (a hidden positive may remain, so the figure is a lower bound). Report P(null) before and after, the AUROC
   of P(null) between intact and deleted groups, and the same for the PN models with a naive null
   (1 − max p̂) so the comparison is like for like.
8. On the model stratum's 27 negative-bearing competitions and the registry's one all-negative singleton:
   P(null) and the rank of the negative members.
9. Harm of the null decision in production terms: for each genome run, how many Stage 1 calls the arm removes
   because null wins, split by labelled positive (harm), experimental negative (correction), unconfirmed.

Candidate-count artifact check (bar text): P(null) against log1p(n) on competitions matched by max margin;
the slope must not be the dominant term (reported, not gated).

## 7. Pass criteria (fixed at freeze)

| id | criterion |
|---|---|
| P0 | Plan sha256 recorded before the first fit; all artifacts hashed before any test, Price or model-stratum score |
| P1 | PU logistic loads in the unchanged JVM reader; Python/JVM probabilities agree to 1e-9 on every table row and on one full production genome run (as the current parity checks do) |
| P2 | Null layer trained at the trainer's default gates in at least one regime, with calibration-split marginal NLL below the uniform baseline |
| P3 | Deletion test: AUROC of P(null) between deleted and intact single-positive test competitions ≥ 0.70 for A3, and higher than the PN naive null |
| P4 | Supported test precision of A3 ≥ A1 (currently 0.25) at matched or higher recall, in the family-excluded regime; the PU-adjusted expected precision reported beside it |
| P5 | Harm: on experimental negatives (model stratum plus cohort), A3 calls no more negatives than A1 at the frozen thresholds; absolute specificity reported |
| P6 | PU-ECE of A3 ≤ PU-ECE of A1 on test |
| P7 | Robustness: P3–P5 hold at the central π̂ and at the lower sweep value 0.20; the direction of each PU–PN difference is reported for every π of S1 and every fold of S3 |
| P8 | Price: A3 scores Price-148 EC candidates in the Price plan's cells; supported precision, calibration and harm reported under its frozen split; no criterion is attached to Price values in this plan |

A failed criterion is reported as failed with the measured numbers; nothing is re-selected after a test read.

## 8. Compute and wall-clock

- Trainer change (PU and ranking objectives, PU marginal likelihood for the null layer, PU reliability,
  deletion test, leave-one-genus-out): Python only, in a new module `benchmark/v3/pu_learning/pu_trainer.py`
  that imports the unchanged `scripts/train_context_model.py` (table reader, model classes, folds, artifact
  writer, reference predictor); written by agy-auto from a spec, reviewed here; existing trainer tests must
  keep passing and new tests cover: PU risk equals PN risk at π = labelled rate with c = 1, the clamp, the
  marginal-likelihood no-label case, and the artifact's unchanged JVM-facing fields.
- Fits: each logistic fit is seconds; all arms, sweeps and 18 leave-one-genus-out refits under one hour on ws.
- Production parity and production replays on the 18 cohort genomes (2 artifacts × 18 genomes, as the
  family-excluded verification did): ws, a few hours, or one Shaheen node for one hour.
- Model stratum: DIAMOND family-excluded claims for three proteomes (minutes on ws), three feature exports
  (ws, under an hour), scoring (seconds).
- Price: 22 genomes × 2 extra learned arms on the Price plan's cells, Shaheen k10018, one node-hour class,
  excluding nodes nid00022, nid00497, nid01033; only after the Price plan's claims are frozen.
- Total well under 200 core-hours; nothing here needs a GPU.

## 9. Dependencies on the other session's outputs (read, never edited)

| needed | from | state |
|---|---|---|
| negative label registry with provenance and review dispositions (`registry/negative_reviews.tsv`, `audit/negative_groups.jsonl`, `registry/labels.tsv`), model-set labels | `benchmark/v3/results/complete-locus-labels-v1/` | registry frozen 13:36, evaluation done (C7 failed); usable now |
| rebuilt family-excluded table and restored feature export (272,366 candidates) | `/mnt/data2/leechuck/gspa-v3-data/complete-locus-labels-v1/{fit/table,features}` | present; this plan copies what it uses into its own folder and records hashes |
| contract-bearing DIAMOND claims for the 22 Price genomes, cells and evaluator | `benchmark/v3/results/price-contract-learned-v1/`, `benchmark/v3/price_contract_learned/` | Step C and the pilot cell done; full cells running on Shaheen; Step E (scores) pending. This plan's Price part waits for its `PLAN.sha256`-bound claims and uses its frozen split |
| the 22 pending negative reviews | complete-locus review queue | optional; whatever is accepted at freeze time enters evaluation, the count is recorded |

## 10. Outputs

`docs/benchmarks/v3-row4-pu-learning.md`; `benchmark/v3/results/pu-learning-v1/` (plan copy, PLAN_SHA256,
artifacts, evaluation, SHA256SUMS, LARGE_FILES.tsv for files over 5 MB); code in `benchmark/v3/pu_learning/`
(trainer module with `--objective`, runner, deletion test, evaluation), tests in `benchmark/v3/tests/`
(`test_pu_learning*.py`); the existing trainer and the JVM tree are not edited; raw work in
`/mnt/data2/leechuck/gspa-v3-data/pu-learning-v1/`; one progress line in
`tasks/v3-supervision-plan-2026-10-28.md`. `tasks/v3-acceptance.md` is updated only after the report, by the
supervising session.

## 11. Known risks, stated before any result

- The SCAR assumption (a true positive is labelled with a probability that does not depend on its features)
  is false in a known direction: well-studied enzymes and central metabolism are labelled more often. The
  prior and c are then overestimated for peripheral functions. The sweep (S1), the negatives-only harm count
  and the deletion test are the three independent checks; if they disagree, the plan says so.
- The Elkan–Noto prior is upward-biased when the classifier is weak; the current AUROC is 0.72–0.76. The lower
  sweep value is there for that reason.
- Experimental negatives are few and concentrated in three model genomes; a harm count of this size can show
  a large effect, not a small one. The power statement is written with the realised count at freeze time.
- The null layer learns "beyond what the label frequency explains"; if c is misestimated the null rate is
  misestimated by the same factor. The deletion test does not depend on c and is the primary null criterion.
- The GBDT PU variant may not fit stably; the logistic arm is primary.
- Hidden multi-positive competitions (unlabelled isoenzymes) are misspecified by one-of-n-or-null; reported
  through the excluded-multi-positive count in labelled groups as a proxy.

## Deviations (append only, dated; none before freezing)

- (2026-11-02 15:50, before any fit or score) Implementation delivered by agy-auto (Gemini) from
  `~/tmp/gspa-v3-agy/row4-pu/spec.md` as new files only (`benchmark/v3/pu_learning/{pu_trainer,pu_evaluate,run_pu_learning}.py`,
  three test files); reviewed line by line by the supervising session; the agent's own report misdescribes some
  numerical details (temperature grid, L2 candidates, clipping) and the code, which follows the spec, is authoritative.
- (2026-11-02 15:50, before any fit) Correction of the counts in §3.2: the figures 85/9/17 (homology) and 44/4/7
  (family-excluded) are competitions with ≥ 2 members and **at least one** labelled positive. After excluding
  competitions with two or more labelled positives, as §3.2 and the trainer require, the single-positive counts are
  train/calibration/test 52/3/7 (homology) and 30/1/4 (family-excluded). The calibration gate of ≥ 5 single-positive
  groups is therefore unmet in both regimes, not only in the family-excluded one. P2 is judged exactly as written
  (expected to fail). Added, fixed now: a secondary null-layer fit with the calibration single-positive gate lowered
  to the available count (3 for homology; 1 for family-excluded, reported as uninformative), labelled "gate deviation";
  its deletion-test and NLL figures are reported beside the gated result and do not count towards P2.
- (2026-11-02 15:50) Test T1 reads "PU risk equals PN risk at π = labelled rate with c = 1" as agreement in
  expectation (unbiasedness) at n = 20,000 within 0.02, since the identity is not exact on a finite sample.
- (2026-11-02 15:50) The non-negative correction is implemented as a smooth barrier (softplus with k = 100 plus a
  quadratic penalty on the violation) followed by Kiryo-style ascent steps on the negative-class risk and, if still
  negative, a bisection towards θ = 0; the tree variant scales the ensemble towards the prior log-odds. It applies
  only when the unbiased fit's negative-class risk is negative; the artifact records whether it fired.
- (2026-11-02 15:50) Model-genome stratum (§4.3): family-excluded contract claims for the three genomes were made
  on ws (contracts match the family-excluded model); the feature export ran on Shaheen (job 15759029) without cohort
  inputs because the genomes are not cohort members. Pre-check: 69,789 candidates, 47,781 (68 %) supported by the
  frozen family-excluded logistic model, 450 rows join model-set positives, 5 of the 6 accepted model negatives
  join a candidate (P9WKI3 / GO:0003938 has no claim). The stratum is scored offline by every arm.
- (2026-11-02 15:50) S3 (leave-one-genus-out) refits use the frozen table rows without recomputing family masking
  across the new folds; it is a sensitivity analysis only, as §4.4 states, and the report says so.
- (2026-11-02 17:25) Production replays (P1, metric 9) run on Shaheen with the release-candidate-3 engine
  (`gspa-3.0.0.jar`, sha256 `8391d73f…`), not the archived family-excluded verification engine `74d9e8fc…`, which
  is not on Shaheen; the JVM reader is the same class. A first replay used the cohort bundle's claims, which carry
  no predictor contract, so every candidate abstained (`candidate_predictor_contract_unavailable`); it is discarded
  and the replay was rerun on the contract-bearing family-excluded claims of `predictor-contract-v1` (uploaded with
  hashes). A one-genome replay is impossible because the cohort pass requires all 18 genomes of the cohort
  metadata; the replay therefore covers all 18 genomes, which is what metric 9 needs anyway.
- (2026-11-02 17:25) Price cells: the first submission ran only cell 0 (Slurm splits `--export` values on commas;
  the cell list now uses `+`) and failed because the other session's `parts` folder was not bound inside the
  container; resubmitted with that folder bound read-only. No score was read from the failed run.
