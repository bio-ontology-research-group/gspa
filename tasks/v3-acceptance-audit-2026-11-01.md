# Acceptance-bar audit of the ten "Not met" rows (2026-11-01)

Read-only audit of `tasks/v3-acceptance.md` (bar table, lines 349–362) against the linked reports in
`docs/benchmarks/v3-*.md` and the archives in `benchmark/v3/results/`. Rows 2 and 10 are already "Met" and are
not examined. For every other row the bar text is quoted exactly, split into its elements, and each element is
marked SHOWN, SHOWN-NEGATIVE (measured fairly, result bad) or MISSING. Unconfirmed labels are never treated as
negatives. "Class" follows the brief: (a) every element shown and the "Not met" rests on a performance
expectation the bar does not state; (b) some element missing but obtainable with concrete new work; (c) the bar
states a performance level that was measured fairly and not reached.

## Summary table

| Row | Class | Missing elements (relative to the bar text) | Proposed work |
|---|---|---|---|
| 1. Actual application | b | Calibration-only thresholds on the Price benchmark (0.5 is a declared fixed diagnostic); calibration metrics (Brier/ECE) for the release-candidate run; under the strict reading of "full-system", a learned/cohort/reaction-context arm on the six Price enzyme predictors (the learned stage abstains on them for lack of a predictor contract) | (i) Pre-registered genus-blocked split of the 22 Price genomes into calibration/test, apply the existing "≥10 calls and ≥90% support" rule to the retained scores, report Brier/ECE from the retained score files: about one day, no new inference. (ii) Native CLEAN with the pinned cohort adapter on the 22 Price proteomes (GPU, ~22 jobs) so that the 10406 `clean`/`clean_diamond_homology` models can score; run cohort + reaction-context + learned stages through the batch route and score with the unchanged Price evaluator: 3–5 days |
| 3. Predictors | b | Versioned genuine outputs exist for 9 of 22 adapters (HMMER, InterProScan, AMRFinder, eggNOG v2/v3, gapsmith, antiSMASH [captured, not fresh], DIAMOND, Pyrodigal/Prodigal, AntiFam [synthetic model only]). Missing: FoldSeek, MMseqs2, SignalP, DeepTMHMM, Metapredict, CRISPR, VFDB, gapseq nucleotide route, geNomad/CheckV/PhiSpy sidecars (controlled fixtures only), dbCAN real search (formatter fixture only); nucleotide attribution; Nextflow `main.nf` predictor containers never run | Per adapter, repeat the AMRFinder recipe (`benchmark/v3/check_amrfinder_route.py`): pinned install, run on the two public proteins already used (P18843, P0ADB7) plus the synthetic control, hash everything, before/after route check with annotation, annotation+integration, direct integration, missing executable, malformed output, header-only output and one dropped-protein coverage-loss case. 0.5–1 day each, ~10 adapters, CPU only (FoldSeek/DeepTMHMM optional GPU). No benchmark dataset needed; the bar is about outputs and paths |
| 4. Learned inference | b | "Honest complete locus labels" exist only as one reviewed training competition (ThrRS pair) and 8 accepted role-truth cases on three held-out model genomes; zero complete calibration/test groups in the supervision table; the no-candidate (null) outcome is therefore untrained. All other elements shown | Grow the complete-locus registry from the Swiss-Prot negative-statement miner (526 activity negatives on 492 proteins, 227 GO NOT rows) joined to the DIAMOND candidate universes of the 18 cohort genomes and the 22 Price genomes; target ≥20 complete competitions with ≥5 in calibration/test, each with a quoted primary reference; then fit and evaluate the null layer. 1–2 weeks curation, one day fitting. Data: Swiss-Prot and GOA (established) |
| 5. Reaction/interaction context | b (a under a loose reading) | The three-mode ablation ran on held-out loci with 521 observed gaps, but the held-out labelled loci are not gap-linked (0 gap-linked test positives in either regime), so the bar's conjunction "held-out enzyme loci with actual gaps" has no power. Organism-specific context was tested on 148 gap-linked held-out enzymes (Price) but not pathway-only vs universal | Run the reaction-context stage in pathway-only, universal and organism-specific modes on the 22 Price genomes with the gapsmith observed gaps already built for the flux-certified stage; priors act without the learned model (the gap-filling prior changed 202 scores in the cohort); score with the unchanged Price evaluator. 2–3 days, CPU |
| 6. Joint recurrence | a | None. Clade holdouts, duplicate strains, shuffled neighborhoods, genuine fragmentation and sequence removal, candidate- and genome-level effects are all measured, pre-registered, with passing invariance controls. Result: no gain; the narrowing test is underpowered (not a negative) | None required to mark Met-with-negative-result. Optional: a removal arm without context so that harm versus sequence-only can be evaluated (the report names this gap) |
| 7. Retention and roles | a (b under a strict reading of "same-operon") | Experimentally supported role truth (8 accepted cases, 3 negative-only loci), real predictor scores, recall and locus ranking, no-candidate abstention are all measured. Whether any registry case is a same-operon pair is not stated in the report | Tag each registry case with operon membership from the RefSeq layouts already built for the role genomes; add ≥3 same-operon cases (non-catalytic subunits) from the negative-statement queue. 2–3 days curation. Registry is small (5 scorable held-out cases); state this |
| 8. Absence/environment | b | Gene attribution is scored separately from phenotype agreement on auxotrophs (masked recovery benchmark) and on complete genomes (off-growth pilot), but the four competing explanations (missing annotation / missing sequence / pathway loss / external supply) were never scored against known truth on real genomes, and no attribution was scored under controlled sequence loss | Run the evidence-guided absence construction (exists: `--experiment-mode evidence_guided`) on three existing inputs with known explanations: the sequence-loss cohort (24 removed labelled enzymes plus 39 labelled removals from the close-relative benchmark: truth = missing sequence), the 480 nutrient auxotroph deletions (truth = pathway loss; supplemented = external supply) and Price-148 (truth = missing annotation); score whether the known explanation is in the feasible-pattern set, separately from growth accuracy. About one week, CPU only. Data: training cohort, nutrient dataset, Price |
| 9. Cross-feeding | b | Measured mono/coculture: shown (negative). Exchanged metabolites: missing. Transporter knockouts in coculture: missing (sulfonate transporter knockouts exist only in monoculture). Separate false-enzyme-assignment effect: missing | (i) False-assignment sensitivity: inject controlled wrong reactions into iML1515/iYO844 GPRs (e.g. assign the deleted auxotroph gene's reaction to another gene, or add a Price-style wrong EC) and re-run the 7,920 coculture plans; report change in bulk/rescue predictions. 1–2 days, ~40,000 LPs. (ii) Exchanged metabolites and transporter knockouts: acquire the E. coli / Salmonella enterica methionine–acetate consortium (Harcombe et al. 2014, measured metabolite and knockout data) or the E. coli auxotroph-pair exchange panel (Mee et al. 2014); 1–2 weeks |
| 11. Agent comparison | b | Correctness/harm: not identifiable, because none of 196 real unresolved cases intersects experimental truth (361 assertions inspected, none supports a selected candidate–function pair). Unseen genomes, real cases, frozen policies, abstention and cost are shown | Run the production anomaly extractor (`run_agent_evidence.py`) on the 22 Price release-candidate outputs (reannotated enzymes are by selection likely gaps) and on the three dense-label model genomes used in the role benchmark (disclosed as development-examined); split by genus into calibration/test; rerun the six simple methods and the live agent. 3–5 days; several hundred agent requests. Data: Price 2018, Swiss-Prot |
| 12. Experiments | b | Measured secondary-product labels: missing (Sco-GEM intake accepted zero cases; fitted-flux circularity and strain/name conflicts). All other elements shown; fixed-precision recovery is zero for every experiment arm because none obtained a calibration threshold, with 13 test cases | Secondary products: freeze the Sco-GEM base model (not the flux-fitted one), use published knockout-abolished production of actinorhodin, undecylprodigiosin and CDA (MIBiG knock-out evidence) under the published media, withhold product observations. 2–3 weeks curation. Also enlarge the masked test split (Price RB-TnSeq cases gated by the design conditions) |

## Per-row detail

### Row 1. Actual application

Bar (quoted): "Whole-proteome claims; matched stage/full-system comparisons; frozen family/genome splits;
calibration-only thresholds; calibration, correction/harm counts and paired genome uncertainty."

| Element | Status | Evidence |
|---|---|---|
| Whole-proteome claims | SHOWN | `v3-price-rc-evaluation.md`: 132 cells, 1,526,919,081 claims, every cell's publication verified (job 15708110). `v3-price148.md`: full candidate matrices (`full_candidate_coverage: true`) for all six predictors |
| Matched stage comparisons | SHOWN-NEGATIVE | Integration without priors vs default priors: byte-identical in all 132 cells (priors act on GO only). Flux-certified refinement (`v3-price-stage2b.md`): CLEAN top-1 60→59, top-3 87→80, 10 top calls changed, 0 corrected, 1 harmed; five other predictors unchanged; CAFA F-max change 0.000 [0.000, 0.000] for all six. Taxon fix (`v3-price-taxon-fix.md`): 465→0 forbidden calls, 0 top calls changed. Six pre-registered follow-ups (`v3-stage2b-followups.md`): no Price gain in any arm. The bar asks for the comparison, not for a gain |
| Full-system comparison | ambiguous; MISSING under the strict reading | If "full system" means all default priors plus the metabolic refiner, it is shown (above). If it means the learned, cohort and reaction-context stages as well, it is missing on Price: those stages require a `predictor_contract` and the six enzyme predictors carry none, so the learned stage abstains (`v3-routes-learned-cohort.md`: 146,160 of 146,160 abstain where the contract is absent). The three Price genomes run with learned+cohort used DIAMOND claims and were not scored against the Price reference |
| Frozen family/genome splits | SHOWN with a caveat | For the learned models: 11/3/4 genome split, disjoint on genome, group, genus and family, `family_disjoint_verified: true`; Price audit: 0 genus conflicts, 0 family overlaps between the cohort and the 22 Price organisms (`v3-training-cohort.md` §2, §4). Caveat stated in `v3-price148.md`: overlap between the six external predictors' pretraining data and the Price proteins is unassessed and cannot be fixed here |
| Calibration-only thresholds | MISSING on Price | `v3-price-rc-evaluation.md`: "Calls are integrated scores of at least 0.5, a fixed diagnostic threshold, not a calibrated operating point." Calibration-only threshold selection exists elsewhere (`v3-sequence-loss.md`: EC threshold 0.44 from calibration genomes; observed-gap analysis) but was not applied to Price |
| Calibration metrics | MISSING for the release-candidate run | `price-reference-rc1/*/report.json` contains no Brier/ECE fields. `v3-price-calibration-audit.md` covers the historical runs only (zero-score baseline has lower Brier than 13 of 14 learned combinations; CLEAN top-1 falls 60→2/0 under the historical learned arms) |
| Correction/harm counts | SHOWN | Stage 2B table (0/1 CLEAN, 0/0 others); follow-ups (harm among gated updates 0.97–0.999; competitive assignment harm 0.851); taxon fix (0/0) |
| Paired genome uncertainty | SHOWN | Genome-bootstrap top-1 intervals (e.g. CLEAN 0.405 [0.330, 0.485]); CAFA F-max delta 95% intervals; `summary_exact.tsv` |

Class (b). The current "Not met" reason ("the full system shows no gain") is a performance expectation the bar
does not state; but two bar elements (calibration-only thresholds, calibration metrics) are genuinely missing on
Price and are cheap to obtain, and the strict "full-system" reading needs the contract-bearing CLEAN run. Note:
if "full-system" is read loosely and the two cheap elements are added, the row becomes (a) with negative results.

### Row 3. Predictors

Bar (quoted): "Versioned genuine tool outputs and complete parser/dispatch paths, including injected failures and
coverage loss."

| Element | Status | Evidence |
|---|---|---|
| Versioned genuine tool outputs | PARTIAL | Genuine, hashed, versioned: HMMER 3.4 (`v3-hmmer-observed.md`, synthetic model `GSPA_TEST_001.1`, real hmmbuild/hmmsearch), InterProScan 5.72-103.0 (job 10095, 2 public proteins + control, 4 matches, 3,951 files hashed), AMRFinder 4.2.7 (17-protein test FASTA, 14 calls, DB 2026-08-07.1), eggNOG-mapper 2.1.13 and 3.0.0-beta6 (67 and 26 claims), gapsmith 0.1.0/DIAMOND 2.0.15 (job 10111, 13,147 rows → 3,469 candidates), antiSMASH 8.0.4 (captured output, 3 regions, 70 memberships; fresh execution open), DIAMOND 2.1.11 (`v3-predictor-contract-replay.md`), Pyrodigal 3.7.1 and Prodigal 2.6.3 (`v3-gene-caller-repair.md`, job 10515). Not genuine or absent: dbCAN (version-4 formatter run on synthetic hits; no real search), geNomad/CheckV/PhiSpy ("do not establish real geNomad execution"), AntiFam (synthetic HMM, genuine AntiFam models not run), FoldSeek, MMseqs2, SignalP, DeepTMHMM, Metapredict, CRISPR, VFDB, gapseq. The old `mock-output/` fixtures are explicitly re-labelled "format fixtures" |
| Complete parser/dispatch paths | PARTIAL | Annotation vs direct-integration parity verified for every repaired adapter; gene-caller dispatch and empty-input dispatch repaired. Nucleotide dispatch: `v3-gapsmith-observed.md` ("Nucleotide gapseq hits cannot become protein claims without coordinate attribution; this missing capability is reported as failure"); `v3-eggnog-observed.md` lists nucleotide dispatch as open. Nextflow `main.nf` predictor containers not run (`v3-routes-learned-cohort.md`) |
| Injected failures | SHOWN for repaired adapters | Typed FAILED / UNAVAILABLE / SUCCESS_ZERO / PARTIAL on malformed, missing-executable, header-only and corrupt input (AMRFinder 7 CLI runs; eggNOG; HMMER; geNomad fixture table with 8 outcomes; gene callers 23 invocations) |
| Coverage loss | PARTIAL | gapsmith partial coverage on missing sequences; antiSMASH partial records; AntiFam scan-coverage states; InterProScan "full upstream coverage loss ... unverified"; `v3-amrfinder-observed.md` "full dispatch/coverage-loss paths remain to be certified" |

Class (b). The word "complete" in the bar is a completeness requirement, not a performance level, and it is not met:
about 9 of 22 adapter classes have genuine versioned outputs. "Real-family accuracy" in the status notes is not in
the bar and should not block the row. The work is mechanical and the recipe exists.

### Row 4. Learned inference

Bar (quoted): "Honest complete locus labels; OOF training; family/genus holdouts; logistic/tree/rule comparisons;
calibration and harm; supported predictor distributions."

| Element | Status | Evidence |
|---|---|---|
| Honest complete locus labels | PARTIAL (MISSING for calibration/test) | `v3-locus-truth.md`: polarity bug fixed; 4,262 rows, 541 positive, 3,721 closed-world unconfirmed, 0 explicit negatives, 0 complete locus groups; `locus_model.trained=false`. `v3-gaf-negative-audit.md`: one reviewed complete two-member training competition (ThrRS catalytic vs editing, PMID 15240874). `v3-role-truth-benchmark.md`: 8 accepted role cases with quoted primary evidence on held-out genomes, 3 negative-only loci. No complete calibration or test group; the null outcome cannot be trained or evaluated |
| OOF training | SHOWN | Group-aware five-fold CV, seed 42, calibration-split temperature scaling; thresholds frozen before test (`v3-training-cohort.md` §5; `v3-family-excluded-training.md`) |
| Family/genus holdouts | SHOWN | 11/3/4 split disjoint on genome, group, genus, family; MMseqs2 30%/80% family masking (76 masked training rows); family-excluded regime (hits ≥30%/80% discarded) |
| Logistic/tree/rule comparisons | SHOWN | `v3-context-comparison-results.md`: 18 arms (3 contexts × rule/logistic/gbdt × 2 regimes); `v3-observed-gap-cohort.md` final-score table; `v3-source-context-comparison.md`: 61 arms including native CLEAN |
| Calibration and harm | SHOWN (negative) | Supported-test Brier 0.088–0.123, AUROC 0.72–0.76, precision 0.25–0.50, recall 0.03–0.35; corrections/harms vs no-gap rules 56/24 (logistic), 38/33 (tree), family-excluded 92/5, 91/23 — against closed-world unconfirmed labels, stated as such. Zero-score baseline audit |
| Supported predictor distributions | SHOWN | `v3-predictor-support.md`: predictor-contract and observed-range gate; 40 CLI calls, 1,921 claims, 34,096 JVM/Python comparisons; 544,732 production comparisons with abstention reasons |

Class (b). "Supported test precision 25–29 %" is not a level the bar states. The only missing element is enough
complete locus labels to train and test the null layer. Work is curation, with the mining tools already written.

### Row 5. Reaction/interaction context

Bar (quoted): "Pathway-only versus universal-network versus organism-specific ablations on held-out enzyme loci
with actual gaps/context."

| Element | Status | Evidence |
|---|---|---|
| Three-mode ablation | SHOWN (negative) | `v3-context-comparison-results.md`: UNIVERSAL / PATHWAY / ORGANISM × rule/logistic/gbdt × homology/family-excluded; 12 fits, 216 production runs, independently audited (job 10223, 6,981 files). No consistent context advantage. Close-relative benchmark: no relation advantage. Role benchmark: contexts identical on held-out role loci |
| Held-out enzyme loci | SHOWN | 4 test genomes; expanded panel of all 37 experimental-reference proteins / 39 EC assertions (15/39 emitted with homology, 5/39 family-excluded; best selected arm 11/39 matched) |
| Actual gaps/context | SHOWN, but not on the labelled loci | `v3-observed-gap-cohort.md`: 521 observed gap records from complete proteome searches, 72 genome/arm executions, 202 score changes. Supervision coverage: gap-linked test positives 0 (homology) and 0 (family-excluded); "the family-excluded table-test metrics contain no gap-linked examples". Organism-specific context on gap-linked held-out enzymes: Price flux-certified refinement (0 corrections, 1 harm; 2 of 148 on the growth path) and off-growth pilot (21 of 148 recovered by targeted tests vs 0) |

Reading A (loose): each noun in the bar is present somewhere → (a), negative result. Reading B (strict): the
ablation must be on held-out loci that have actual gaps; the labelled loci have none, so the measurement is
unpowered for exactly the case the item is about → (b). I recommend reading B. The Price set is the natural
gap-linked held-out set and gapsmith models for all 22 genomes already exist (`benchmark/v3/stage2b/`).

### Row 6. Joint recurrence

Bar (quoted): "Clade holdouts; duplicate strains; shuffled neighborhoods; genuinely fragmented/missing sequence;
candidate-level and genome-level effects."

| Element | Status | Evidence (`v3-wp4-cohort-recurrence.md`, `v3-sequence-loss.md`, `v3-cohort-controls.md`) |
|---|---|---|
| Clade holdouts | SHOWN | Genus-level (`CLADE`) vs species-level (`STRAIN_CLUSTER`) holdouts; invariance control passes (54 rows, 0 changed, after the permutation-seed fix); congener peers (6 per family) added for the first time |
| Duplicate strains | SHOWN | Exact copy in the same cluster: 0 change; copy relabelled as new species: inflation up to 0.21–0.67 (reported as the expected inflation); earlier 5× duplication: probabilities exactly invariant for 346,269 candidates after refit |
| Shuffled neighborhoods | SHOWN | Orthogroup shuffle (5 seeds) and pathway shuffle (3 seeds); label-free narrowing not above shuffles for pathway/organism relations; earlier shuffle: net supported-call change +1, interval −3 to +3 |
| Genuinely fragmented/missing sequence | SHOWN | Nucleotide-block removal: 22,457 of 49,463 proteins retained; EC top-1 14/37 → 6/37 (change −0.216, paired interval −0.417 to −0.050); fragmentation tiers 100/25/5 kb: recovery 9→5 (pathway), 10→5 (universal) |
| Candidate-level and genome-level effects | SHOWN | Per-assertion toward/away counts, per-protein entropy changes, per-target recovery, paired McNemar (R2 vs R0: 0 gains, 0–4 losses, p ≥ 0.125), family-block bootstrap |

Class (a). Every element is measured, pre-registered, with passing invariance controls (the first run was
invalidated and repeated). The narrowing test is "not testable" under the power gate (9/8/5 assertions vs 20
required): that is a no-result, not a failure. The bar names no performance level; the row can be marked Met with
the negative result stated. One internal gap: harm versus sequence-only after true removal was not evaluable
(no removal arm without context).

### Row 7. Retention and roles

Bar (quoted): "Real predictor scores and difficult same-operon alternatives; experimentally supported role truth;
recall and locus ranking, including no-candidate cases."

| Element | Status | Evidence (`v3-role-truth-benchmark.md`, `benchmark/v3/results/role-truth-v1/`) |
|---|---|---|
| Real predictor scores | SHOWN | Legacy 9-source claims (DIAMOND, Pfam), role-fenced DIAMOND claims, learned priors on 61,649–67,291 of 101,705 model-genome rows |
| Difficult same-operon alternatives | AMBIGUOUS | Registry cases are paralog pairs (CoaA/CoaX, GuaB1/2/3, FolP1/2, GpmI/YhfR, TalB/FsaA/FsaB, Fgd1/Rv0132c) and one subunit pair (ThrRS catalytic/editing). Margins are small (0.000 to −0.026 with legacy claims). Whether any pair is co-transcribed is not recorded; strict reading → MISSING, loose reading ("difficult within-genome alternatives") → SHOWN |
| Experimentally supported role truth | SHOWN | 8 accepted cases with PMID/PMC-quoted inactive paralogs, 3 pending, 3 negative-only loci; Swiss-Prot negative-statement miner (526 activity negatives) feeds the queue |
| Recall and locus ranking | SHOWN | Strict/lenient top-1, MRR, mean margin per arm; scored 4–7 of 11 cases; candidate-absent cases counted. Role-fenced models on held-out genomes: sequence alone 4/4 strict top-1; logistic keeps 4/4 with wider margins; gbdt under context 2/4 (inverts GuaB and FolP) |
| No-candidate cases | SHOWN | Negative-only loci: 0 of 7 called in every arm (abstention correct) |

Class (a) with the stated ambiguity. The "Not met" reasons ("hard paralog loci stay tied or inverted; no labelled
role training set") are, respectively, a performance expectation and an implementation-column item; neither is in
the bar. Sample size (5 scorable held-out cases) must be stated prominently.

### Row 8. Absence/environment

Bar (quoted): "Gene attribution scored separately from growth/KO accuracy on complete genomes, auxotrophs and
controlled sequence loss."

| Element | Status | Evidence |
|---|---|---|
| Growth/KO accuracy on auxotrophs | SHOWN-NEGATIVE | `v3-nutrient-phenotypes.md`: 1,144 assays × 6 arms (6,864 conditions); best assay agreement 84.8% (B. subtilis, added inorganic, O2 20); paired rescue 20–69% vs constant-rescue 80%; E. coli predicts no growth in all printed-medium primary assays. Bar does not state a level |
| Gene attribution scored separately, auxotrophs | SHOWN | Masked recovery benchmark (`v3-wp5-masked-recovery.md`): function recovery (sequence-only 3/13 = 0.231 [0.082–0.503], harm 0; all experiment arms 0 by frozen policy) reported separately from called-pair phenotype agreement (0.979 for sequence-only, sensitivity scope) |
| Gene attribution, complete genomes | SHOWN-NEGATIVE | Off-growth pilot: reaction re-added (51) vs hidden protein recovered (21/148) vs harm (0.85) reported separately; RB-TnSeq knockout readout AUC 0.51 |
| Gene attribution under controlled sequence loss | MISSING | `v3-sequence-loss.md` scores EC top-1 on retained/removed loci (6/37 after loss; 0/24 removed) but states "recovering that explanation [missing sequence] by the inference model remains a separate acceptance requirement" |
| Competing explanations scored against truth | MISSING | `v3-absence-comparison.md` and `v3-guided-absence.md`: four-scenario comparison runs only on synthetic fixtures (147 and 100 GLPK comparisons); `gene_attribution: unverified`, `selected_explanation: null` by design; "Real cohort comparisons against independent measured absence explanations ... remain open" |

Class (b). The mechanism and three suitable real inputs with known explanations exist; what is missing is running
the comparison on them and scoring the explanation separately from growth.

### Row 9. Cross-feeding

Bar (quoted): "Measured mono/coculture, exchanged metabolites and transporter KOs; separate false enzyme assignment
effects."

| Element | Status | Evidence |
|---|---|---|
| Measured mono/coculture | SHOWN-NEGATIVE | `v3-nutrient-community-results.md`: 1,320 cocultures with matched monocultures, 7,920 plans, 39,599/39,600 LP agreements. Printed medium: B. subtilis bulk 30.67–32.50% vs always-growth 69.33%; E. coli predicts no growth (= never-growth baseline). Added-inorganic B. subtilis bulk 74.33% vs 69.33% (30 corrections, 0 new errors); rescue +12 points but dependence-aware range [−11.7, 45.8] |
| Exchanged metabolites | MISSING | "No measured exchanged metabolite ... is verified here"; transfer fluxes are "conditional candidates" |
| Transporter KOs | PARTIAL/MISSING | `v3-sulfonate-diagnostic.md`: E. coli tauABC/ssuABC deletion panel, 260 cases, 414 LPs, 6 butanesulfonate discrepancies — monoculture only, model already contains the transporter GPRs. No coculture transporter knockout |
| Separate false enzyme assignment effects | MISSING | "No ... false enzyme-assignment effect is verified here" |

Class (b). The first element is measured fairly with a negative result; the other three are absent.

### Row 11. Agent comparison

Bar (quoted): "Unseen genomes/families; real unresolved cases; calibration-frozen decisions; correctness/harm/
abstention and complete inference/tool cost."

| Element | Status | Evidence (`v3-stage3-production-comparison.md`, `stage3-agent-evaluation-v2/validation.json`) |
|---|---|---|
| Unseen genomes/families | SHOWN with caveat | 3 calibration + 4 test genomes held out from model training by genome/genus/family; the report calls it "the previously examined development cohort, not a new blind test" |
| Real unresolved cases | SHOWN | 196 cases from the production anomaly extractor: 182 pathway gaps, 10 borderline conflicts, 4 borderline essential functions; one genome at the 150-case cap |
| Calibration-frozen decisions | SHOWN | Rule fixed before test (≥10 calibration decisions, ≥90% precision); all seven policies immutable; all abstain by rule |
| Correctness/harm | MISSING (not identifiable) | 361 experimental assertions inspected; none supports a selected candidate–function assignment; one selected protein overlaps the labelled set |
| Abstention | SHOWN | All 158 test cases abstained by every method; agent's uncalibrated actions 193 abstain / 3 demote |
| Complete inference/tool cost | SHOWN | Shared production 359.18 s; agent 392 requests, 589.93 s, 730,626 prompt and 71,169 completion tokens; billing rate not inferred |

Supplementary: the masked recovery benchmark does provide a calibration-frozen agent comparison with correctness,
harm (0), abstention and LP/request cost on held-out organisms, but on masked known functions, not real unresolved
cases; 32 of 82 agent runs failed by repeating a visited condition. The historical 40-case comparison is invalid
(answer channels in public IDs). Class (b): the missing element is truth that intersects real unresolved cases.

### Row 12. Experiments

Bar (quoted): "Mask target assignments from ALL evidence/model inputs; held-out organisms and measured media/KO/
secondary-product labels; single-medium, fixed, random, simple evidence and agent under equal solver budgets;
fixed-precision recovery/harm/abstention/cost."

| Element | Status | Evidence (`v3-wp5-masked-recovery.md`, `v3-offgrowth-discovery.md`, `v3-stage2b-followups.md`) |
|---|---|---|
| Mask from ALL inputs | SHOWN | Masking audit: 0 gating alias hits, 162 adjudicated shared identifiers, 70 problems with 0 invariance violations, 82 transcripts with 0 hits; off-growth masking audit passes all 148 cases; shuffled-gene control recovers 0 |
| Held-out organisms | SHOWN | E. coli EcNR1/MC4100 and B. subtilis 168 test; Keio calibration (different organism coverage, flagged); 22 Price bacteria |
| Measured media labels | SHOWN | 44 environments × supplements; Biolog PM1–PM4 (293 wells usable); RB-TnSeq 1,400 condition experiments |
| Measured KO labels | SHOWN | Keio auxotrophs; 480 nutrient deletions; RB-TnSeq fitness (AUC 0.51, chance) |
| Measured secondary-product labels | MISSING | `v3-secondary-source-intake.md`: "No biological cases or function labels have been accepted"; `v3-germicidin-identity.md`: name/formula conflicts; fitted-flux circularity |
| Five arms under equal budgets | SHOWN | minimal medium, fixed panel, random (10 seeds), evidence rule, agent; k ∈ {3, 6} with LP caps; exhaustive and sequence-only references |
| Fixed-precision recovery/harm/abstention/cost | SHOWN-NEGATIVE / underpowered | No experiment arm reached 0.90 precision on ≥10 calibration decisions, so all abstain (recovery 0); sequence-only 0.231; harm 0 everywhere; 13 test cases, no paired difference detectable (p = 0.25); LP and request costs tabulated. Off-growth: 21/148 recovered vs 0, H1/H2 hold, H3 (≥80% kept at 25% of tests) not met: 11 of 21; harm 0.85. Follow-ups: 25/54 vs bar 35; panel 40% vs bar 25% |

The pre-registered criteria (H3, ≥35 of 54, ≤25% of tests) are performance bars of the individual designs, not
of the acceptance bar text, which asks that these quantities be measured at fixed precision. Class (b) because
secondary-product labels are missing; with that element added and the small test split stated, the row would be
(a) with negative results.

## Three most valuable closable items

1. Row 1, cheap half: calibration-only threshold on a pre-registered genus-blocked split of the 22 Price genomes,
   plus Brier/ECE for the release-candidate scores. One day, no new inference; converts two MISSING elements to
   SHOWN. The expensive half (native CLEAN with the pinned cohort adapter on the 22 Price proteomes so that the
   existing `clean`-contract models can score) also closes the strict "full-system" reading and feeds rows 4 and 5.
2. Row 11: anomaly extraction on the 22 Price release-candidate outputs and the three dense-label model genomes,
   then rerun the seven methods. This is the only way to make correctness/harm identifiable on real unresolved
   cases; 3–5 days and a few hundred agent requests.
3. Row 8: run the existing evidence-guided absence comparison on the sequence-loss cohort (truth: missing
   sequence), the nutrient auxotrophs (truth: pathway loss / external supply) and Price-148 (truth: missing
   annotation), scoring the explanation separately from growth. One week, CPU only, no new data. The same
   model-perturbation machinery gives row 9's false-assignment sensitivity in another one to two days.

Rows 6 and 7 need no new work to be marked Met with negative results stated (row 7 with the "same-operon" caveat).
