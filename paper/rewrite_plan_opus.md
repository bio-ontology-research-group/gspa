# Rewrite plan for `paper/main.tex` (and supporting `supplement.tex` changes)

**Reading note.** I loaded `SKILL.md` through the Claude copy of the skill. I could not read `robert-scientific-style.md` or `revision-evidence.md`, because both paths were permission-blocked (the `.gemini` and `.claude` copies). This plan follows SKILL.md and your direction above. If you grant read access to those two files, I'll check the plan against them before drafting.

---

## 0. Things you need to know before the rewrite

1. **Four citation keys are missing from `references.bib`:** `you2018golabeler`, `amarasinghe2020longread`, `logsdon2020longread` and `kolmogorov2020metaflye`. `you2019netgo` and `wang2023netgo3` are present. The `niu2026meteor` entry still says *bioRxiv*; to cite METEOR as PSB 2027, the entry needs the proceedings details. I won't invent bibliographic records. I can draft them for a `grounded-research` check, or you can supply them.
2. **Driver 1 is responsible for most of the multi-predictor Price-149 gain.** The current abstract credits "+0.0528 to +0.1987 F_max" to "multi-predictor refinement". The table shows a different split:
   - GraphEC rises from 18 to 48 at the Stage 0 Noisy-OR step and stays at 48 at Iter 3.
   - EnzBERT peaks at Noisy-OR (62) and falls to 54 at Iter 3.

   The narrative should credit pooling (Driver 1) for this gain. Stages 1/2B contribute the extra gains for MAPred (58 → 73), DeepProZyme (54 → 59) and CLEAN (63 → 66). This makes the three-driver story more convincing, not less.
3. **The zero-flux contrast is not specific to GSPA.** "0.00 zero-flux additions" holds by construction. The threshold positive-cost LP also has zero zero-flux additions (81.31 per genome, all flux-carrying). Against MILP, the scientific contrast is *positive vs. negative cost*. GSPA's own contribution is parsimony (5.7 vs. 81.3 additions per genome) plus locus attribution. The Results should say exactly that.
4. **The Stage 3 comparison favors supervised classifiers.** 8-feature supervised classifiers with the guardrail reach F₁ 0.92–0.97 against the agent's 0.91. On LRO-50, 3-scalar classifiers reach 90–94% against 82.7%. The defensible claim is: zero-shot, no labels, and the only configuration that proposes ancestor relaxations at scale (13.7 vs. ≤2). The abstract should not present 34 → 83% as a standalone headline.
5. **Price-149 uses a reduced, script-level form of Stages 1 and 2B** (neighbor window, fixed coefficients, outer damping schedule), not the CLI pipeline. Methods must state this as the protocol used. The full formulas move to the Supplement.
6. **Table count.** The paper has eight tables. The suite overview (`tab:suites`) is the eighth, and its last column consists entirely of file paths. I recommend moving it to the Supplement as the provenance table (S-Table "Benchmark-to-file map"), which leaves Tables 1–7 in the main text, as your direction assumes.
7. **The figures need regenerating.** Figure panels show `A1…A9`, `E1…E5`, `M1…M4`, "Config 1–4" and "pillar" labels. Updating the captions alone leaves the codes visible, so `generate_figures.py` and `generate_figure6.py` need relabeling. Also, Fig. 3b and Fig. 4c (current file `fig6`) both show role-stratified LRO; I suggest dropping the panel from Fig. 4 and using that slot for the 1,000-genome structural audit or the BGC result.
8. **Suite numbers ("Suite 1–8") read as internal identifiers.** I recommend removing them from main-text headings, prose and captions; the Supplement keeps the suite ↔ question mapping.

---

## 1. Title

**Proposed:** *Combining per-protein predictors with whole-genome context and metabolic reasoning improves protein function annotation*

**Minimal-change alternative:** *Genome-scale context, stoichiometric refinement, and neuro-symbolic reasoning for proteome function annotation*

I recommend the first: it names the three drivers and states the finding. Update `pdftitle` to match.

---

## 2. Abstract (five sentences, about 180 words, four headline numbers)

1. **Problem and Driver 1 premise.** Per-protein predictors combine sequence homology, domain profiles, protein language models and structural similarity into accurate calibrated scores. They score each protein independently of the genome that encodes it. As a result, divergent pathway enzymes stay unannotated, subunits of heteromeric complexes are called in isolation, and homology transfer assigns developmental terms to clades that lack the anatomy.
2. **Timeliness, GSPA and Driver 2.** Long-read sequencing and modern assemblers now recover complete genomes routinely. GSPA therefore starts from a calibrated ensemble of per-protein predictors and refines it using the whole genome: operons, core essential functions, complex membership, and GO taxon and Uberon/PO/CL anatomical constraints.
3. **Driver 3.** GSPA goes beyond homology in three steps:
   - it couples reaction-local genomic context with parsimonious flux balance analysis to attribute gap-filling reactions to specific loci, competing against an unassigned-gene hypothesis;
   - it resolves the fewer than 2.5% of loci that stay ambiguous with a tool-using language model;
   - it accepts a model proposal only after verification against the ontologies.
4. **Findings (four numbers).**
   - Across bacteria, archaea, viruses and eukaryotes, GSPA removes all taxonomic and anatomical violations, including **1,009** in six eukaryotic proteomes, while keeping protein-centric accuracy.
   - Genomic context raises recovery of held-out complex subunits from **56% to 94%**.
   - Negative-cost gap-filling leaves **84%** of added reactions without flux. GSPA adds only flux-carrying reactions, attributed to loci, and improves enzyme prediction on fitness-validated enzymes and BacDive phenotypes.
5. **Implication and availability.** The genome, not the isolated sequence, should be the unit of functional annotation; code, constraint tables and benchmarks are available at the repository.

**Cut from the abstract:** all CIs, the Price-149 per-predictor Top-1 lists, the FVA/knockout certificate counts, the BacDive F₁/AUC values, the 3-seed means, "Stage 0/1/2A/2B" labels (keep only the concepts), the model name, and "GSPA v2.1.0" (the version belongs in Methods and availability).

---

## 3. Introduction (six paragraphs)

**P1 – Problem definition.**
- **Role:** define genome-scale annotation as a computational problem.
- **Inputs:** the proteome, gene coordinates/strand, organism lineage, and reference ontologies and biochemistry.
- **Output:** a joint GO/EC assignment.
- **Constraints:** the five invariants: essential processes present; pathways stoichiometrically connected; obligate complexes complete; taxon constraints satisfied; anatomy entailed by developmental terms present in the clade.
- **Edits:** keep the formal notation but drop the four GO IDs; one example term (*heart development* in a nematode) makes the point.
- **Ends on:** a per-protein score can be correct protein by protein yet violate every one of these invariants.

**P2 – Concept 1: per-protein inference and ensembles (sets up Driver 1).**
- Describe four representations: alignment to annotated homologues; profile HMM domains; learned sequence and PLM embeddings; structural alphabets/fold similarity. Cite each as a compact citation group, not a list of tool names.
- Ensembles: GOLabeler and NetGO/NetGO 3 showed that learning-to-rank across complementary channels beats any single channel, because the channels err on different proteins. This is the base GSPA builds on.
- Failure: all of these methods treat the protein as conditionally independent of the rest of the proteome. Below about 30% identity, pathway enzymes become holes, domain hits give singleton subunits, and distant transfers bring clade-inappropriate terms.
- Remove the 13-tool CAFA sentence. CAFA appears only as the evaluation setting, with its citations.

**P3 – Concept 2: genomic context, and why the shift is timely (sets up Driver 2).**
- Classical context signals: conserved neighborhoods, co-directional operons, phylogenetic profiles, STRING-type association networks. Genomic language models learn the same co-occurrence implicitly.
- What these signals capture: local co-occurrence on a contig.
- Gap: local context does not impose cellular invariants — essential-function coverage, complex stoichiometry, or formal taxon/anatomy axioms.
- Timeliness: context-based inference was limited by fragmented short-read contigs; long-read sequencing (amarasinghe2020longread, logsdon2020longread) and assemblers such as metaFlye (kolmogorov2020metaflye) now yield complete or near-complete genomes and MAGs, so whole-genome context is available for most new genomes. This sentence motivates the shift of the unit of inference.

**P4 – Concept 3: systems-level and ontological constraints (sets up Driver 3 and the Stage 2A gap).**
- Metabolic reasoning, in order:
  - pathway hole fillers (Green & Karp) used pathway membership to rank candidates;
  - constraint-based reconstruction (ModelSEED, CarveMe, gapseq, GapMind) uses stoichiometry to add reactions that restore biomass flux;
  - METEOR coupled gap-filling with enzyme prediction. Its negative-cost reaction-indicator MILP rewards adding high-scoring reactions, which brings in reactions that carry no flux, and it boosts every candidate gene of a reaction uniformly, so paralogs are inflated.
- Ontological constraints: GO taxon constraints (Deegan) and the GO–Uberon/PO/CL logical definitions (Mungall) formalize which functions a clade can have. Existing post-filters delete a violating term together with its valid ancestors, losing correct general annotations.
- Final sentence: tool-using language agents can integrate heterogeneous evidence on hard cases, and formal ontology verification is needed to keep their proposals consistent.

**P5 – GSPA mapped to the gaps.**
- Introduce GSPA directly, then one sentence per driver, each naming the gap it closes:
  - **Driver 1, Stage 0:** calibrated pooling of correlation-grouped channels, or ingestion of any single predictor.
  - **Driver 2, Stage 1 priors and Stage 2A:** essentiality, complex/pathway coherence, operon coupling, and taxon/anatomy priors; deterministic enforcement that relaxes a forbidden term to its nearest valid ancestor instead of deleting it.
  - **Driver 3, Stage 1 RLGC, Stage 2B and Stage 3:**
    - RLGC and operon coupling;
    - role-aware DNF complex/isozyme aggregation;
    - positive-cost lexicographic pFBA, which by construction adds only flux-carrying reactions;
    - competitive locus responsibility against an explicit unassigned-gene hypothesis;
    - guarded LM reasoning on the fewer than 2.5% of loci that remain ambiguous.
- State the advances over METEOR explicitly: positive cost, locus-level attribution, genomic context.
- Model name moves to Methods.

**P6 – Evaluation, impact, availability.**
- One sentence per evidence type:
  - temporal and prospective CAFA-style holdouts, and reference-genome panels across four domains (accuracy and invariants);
  - leave-reaction-out and leave-subunit-out recovery (masking experiments);
  - curated GEMs, RB-TnSeq fitness-validated enzymes and BacDive phenotypes (experimental anchors).
- One sentence of supported impact: annotation of complete genomes and MAGs that is consistent and auditable.
- Availability: URL and release.
- No numbers in this paragraph; the current P4 number dump is deleted.

---

## 4. Results

**Proposed reorganization.** The Results follow the drivers in order. Dataset definitions move to Methods, so each Results paragraph opens on the question.

- 2.1 Combining complementary predictors establishes the per-protein baseline (Driver 1)
- 2.2 Whole-genome priors restore structural invariants without loss of accuracy (Driver 2)
- 2.3 Operon context completes heteromeric complexes (Driver 2)
- 2.4 Anatomical constraints and ancestor relaxation in multicellular eukaryotes (Driver 2)
- 2.5 Genomic context and stoichiometry localize missing enzymes (Driver 3)
- 2.6 Flux-parsimonious repair and locus attribution vs. negative-cost gap-filling (Driver 3)
- 2.7 Fitness-validated enzymes and experimental phenotypes (Drivers 1 + 3)
- 2.8 Guarded language-model reasoning on residual ambiguity (Driver 3)

Each paragraph below gives (a) the question, (b) how it was tested, (c) the numbers to keep and what to cut, (d) what it means for the drivers, and (e) the bridge to the next paragraph.

### 2.1 Driver 1: per-protein baseline (Fig. 2d; Supp.)

**¶1 – Reference panels**
- (a) Does pooling complementary channels beat each channel and production annotators under the same scoring?
- (b) Ancestor-propagated micro and CAFA F_max on the 13- and 21-genome panels against PGAP, InterProScan, metagenomic-deepFRI, gLM/gLM2 and each Stage 0 channel.
- (c) Keep:
  - 3-channel GSPA vs. PGAP on 11 prokaryotes: micro F_max 0.84 vs. 0.78;
  - "exceeds every constituent channel on the 21-genome panel" (best single channel ESM-2 kNN, 0.858 vs. 0.869).

  Cut: the per-channel strict-holdout values, the 7- vs. 3-channel splits and S_min. They go to Fig. 2d and a Supplementary table, which should also list the strict leaf scores.
- (d) Channels are complementary, and the ensemble is a strong starting distribution; this is Driver 1.
- (e) Do these gains hold on proteins without close homologues?

**¶2 – Temporal holdout and the twilight zone**
- (a) Does the baseline, and the genome-level refinement, hold up on proteins annotated after the training cutoff, including those below 30% identity?
- (b) CAFA temporal holdout, overall and in the below-30%-identity subset.
- (c) Keep:
  - Stage 0 macro F_max 0.693 vs. the best published comparator (ProteInfer 0.614);
  - Stage 1+2 lifts the below-30% subset from 0.672 to 0.678.

  Cut: per-aspect values for each iteration, and the other comparators (to the table). The BP-peaks-at-Iter-1 observation goes to Methods as the stopping-rule justification.
- (d) The ensemble carries bulk accuracy, and genome context helps where homology is weakest. This is a bridge between Drivers 1 and 2.
- (e) Do the same patterns hold for enzymes and for single predictors?

**¶3 – Enzyme holdout and multi-predictor pooling**
- (a) Is pooling also the main source of gain for EC prediction?
- (b) Holdout-125, and the Noisy-OR column of Price-149.
- (c) Keep:
  - Holdout-125 F_max: Stage 0 pooling 0.438 vs. CLEAN 0.436 alone, with Top-5 rising further after refinement (54 → 55);
  - GraphEC Top-1 18 → 48 after pooling (Table 3).

  Cut: the Top-3/Top-5/MRR for every method, and the METEOR v8 holdout values.
- (d) Pooling complementary EC predictors gives the largest single step; Driver 1 is necessary.
- (e) Per-protein pooling leaves genome-level invariants unchecked, which is the subject of 2.2.

**Prospective CAFA6 check:** keep it here as ¶4 or move it to 2.8 (my recommendation is 2.8).
- (a) Does the pipeline generalize prospectively when only isolated sequences are available?
- (c) Keep: macro F_max 0.543 vs. 0.496 (PLM+3Di MLP) and vs. 0.513 (unguarded agent).
- Place it in 2.8, because the comparison isolates the guardrail and the agent.

### 2.2 Driver 2: whole-genome priors (Table 1, Fig. 2a–c)

**¶1 – Ablation of the priors**
- (a) Which genome-scale priors repair which invariant, and at what cost to accuracy?
- (b) Cumulative ablation on nine reference genomes, with hyperparameters frozen from *E. coli* K-12; report the eight held-out genomes.
- (c) Keep:
  - full Stage 1 raises pathway coherence 0.50 → 0.51 and complex coherence 0.86 → 0.89 (Wilcoxon p = 0.008 for pathway);
  - protein-centric CAFA F_max stays 0.804 at baseline and after refinement.

  Cut: the s.d.s, every single-prior intermediate value, per-genome increments, and the per-genome violation counts (to Table 1 and the Supplement).
- (d) The priors move structural invariants that per-protein scoring cannot see, while accuracy is unchanged; this is Driver 2.
- (e) The taxon/anatomy priors were the one block that also moved accuracy. Is that because they suppress transferred errors?

**¶2 – Perturbation test**
- (a) Do the constraints suppress cross-clade transfers specifically?
- (b) 900 injected eukaryotic and clade-restricted claims in nine prokaryotes.
- (c) Keep: 122 of the 125 claims that survive calibration are suppressed (97.6%). Cut: the sub-splits (55/55, 67/70) and the calibration-only figure.
- (d) Ontology constraints catch the errors that calibration lets through.
- (e) Does this scale from nine genomes to a thousand, and across domains?

**¶3 – Scale and domains**
- (a) Do these effects hold at RefSeq scale and outside bacteria?
- (b) 1,000-genome structural audit and the 88-proteome, four-domain panel.
- (c) Keep:
  - taxon-violating predictions per genome fall from 23,384 (Stage 0) to 0 (DeepFRI: about 107,000);
  - complex coherence rises from 66% to 100%;
  - one sentence that the bootstrap CIs exclude zero for Bacteria and Viruses.

  Cut: all per-domain BP/CC/iteration values and the runtime ranges (runtime moves to a single Methods/Supplement sentence).
- (d) Organism-wide consistency is enforceable at genome scale and is cheap (seconds per genome).
- (e) Coherence counts whether complexes are complete, but does genomic context identify the *right* subunit? That question leads to 2.3.

### 2.3 Driver 2: complexes (Table 5, Fig. 5a,b)

**¶1 – Leave-subunit-out recovery**
- (a) Can genomic context recover a weakly supported subunit among matched decoys?
- (b) Leave-subunit-out recovery: 175 trials in six prokaryotes; the true subunit's evidence is attenuated to 0.22–0.38, against five matched decoys (definition moves to Methods).
- (c) Keep:
  - Hit@1 0.56 (baseline and coherence prior) → 0.94 with operon context;
  - "the coherence prior alone raises recall but not ranking".

  Cut: the precision/recall/F₁@0.5 cascade and the 145-trial repeats (to Table 5).
- (d) Knowing that a complex is present is not enough to identify its subunits; local architecture identifies which protein fills the role. This is the clearest evidence for Driver 2.
- (e) Where operons do not co-locate subunits, context helps less.

**¶2 – *Synechocystis***
- (a) What happens when subunits are dispersed across the chromosome?
- (b) The *Synechocystis* row of Table 5.
- (c) Keep: Hit@1 0.64 with operon context, rising to 0.68 with Stages 2–3; that is the full set.
- (d) This identifies the regime for Driver 3 (non-local reasoning), which 2.8 takes up.
- (e) Complexes are one organism-wide constraint. Anatomy in multicellular eukaryotes is the other.

### 2.4 Driver 2: eukaryotic anatomy (Table 6, Fig. 5c,d)

**¶1 – Endogenous violations**
- (a) Do real eukaryotic predictions, without any injection, violate anatomical constraints?
- (b) Six eukaryotes, unperturbed predictions.
- (c) Keep:
  - 1,009 endogenous violations, concentrated in *C. elegans* and *D. melanogaster*;
  - all removed, with mean CAFA F_max 0.710 → 0.712.

  Cut: per-species F_max deltas.
- (d) Homology transfer from mammals and plants produces these violations routinely; Driver 2 corrects them.
- (e) Which constraint layer is needed?

**¶2 – Which layer is needed**
- (a) Are GO taxon constraints sufficient?
- (b) Stress test that compares the constraint layers.
- (c) Keep:
  - direct GO taxon constraints leave 44% of violations unresolved;
  - the soft prior leaves high-confidence violations in invertebrates (559 in *Drosophila*).

  Cut: per-species counts and the two different "mean elimination" definitions.
- (d) Anatomical axioms (Uberon/PO/CL) are required, and deterministic enforcement must follow soft priors.
- (e) Enforcement can delete or relax. Which is better?

**¶3 – Relaxation vs. deletion**
- (a) Is ancestor relaxation better than deletion?
- (b) Transfer-evaluated F_max and recall against unmodified Swiss-Prot.
- (c) Keep:
  - relaxation beats deletion in all four multicellular animals and plants (largest gain: *Drosophila*, +0.016 F_max);
  - relaxation retains 95–99% of annotation information content.

  Cut: per-species IC values and recall deltas. The yeast/moss precision trade-off is reported as one positive sentence ("in unicellular yeast and in moss, deletion gives higher precision because the injected terms fall on housekeeping proteins"). The developmental-boost 502/618 result goes to the Supplement or a single clause.
- (d) Relaxation keeps the valid general annotation; this answers the naive-filter gap set up in the Introduction.
- (e) Driver 2 enforces consistency. Recovering a missing function needs Driver 3.

### 2.5 Driver 3: localizing missing enzymes (Table 2, Fig. 3)

**¶1 – Leave-reaction-out with embeddings masked**
- (a) Can genomic and pathway context recover a masked enzyme when homology and embedding donors are removed?
- (b) Leave-reaction-out in *E. coli*. Track A masks homologues and embedding donors with cosine similarity above 0.5 (definitions move to Methods; the "All"/"Pool" distinction is stated once).
- (c) Keep: Hit@1 2% for sequence ranking and for METEOR's uniform boost, 23% with RLGC + coherence. Cut: GBDT, cosine-masked, Pool and MRR values (to Table 2).
- (d) Reaction-local context localizes enzymes that sequence similarity misses; METEOR's uniform boosting gives no ranking signal (Driver 3 vs. METEOR).
- (e) How does context compare with protein language models?

**¶2 – Embeddings unmasked**
- (a) What does context add once embeddings are allowed?
- (b) Track B (homology masked, embedding centroids available) and *B. subtilis* out-of-sample.
- (c) Keep:
  - PLM centroids reach Hit@1 0.62; the full model reaches 0.64;
  - *B. subtilis*, frozen parameters: 0.43, compared with 0.11 under Track A.

  Cut: the CIs, the 143-reaction intersection, and the stratum breakdown (Strata A–C go to Table 2 with a single sentence: "gains hold across core, catabolic and secondary strata").
- (d) Stated positively: in a well-annotated proteome, embeddings carry most of the signal, and context resolves ties between paralogs. Context adds most when embedding neighbors are absent, which is the case for MAGs.
- (e) Ranking a gene for a known gap is one task. Deciding which reactions to add is another.

**¶3 – Biosynthetic gene clusters and MAGs**
- (a) Does context-aware refinement avoid over-promotion outside primary metabolism, and what does it propose for uncultured genomes?
- (b) The 106-BGC audit, then 97 MAGs (an illustration, labeled as such).
- (c) Keep:
  - BGC F₁ 0.74 → 0.90 at 99.6% non-BGC specificity;
  - in the MAGs, 568 high-confidence candidate fills (about 6 per MAG).

  Cut: the core/tailoring split and the completeness percentages.
- (d) The approach scales to uncultivated lineages; these fills are hypotheses for experimental testing.
- (e) Which reactions should be added to a draft network, and which genes are responsible?

### 2.6 Driver 3: repair vs. METEOR (Fig. 4a; Supp. Table S5)

**¶1 – Flux of added reactions**
- (a) Do negative-cost formulations add reactions that carry no flux?
- (b) Audit of 108 bacterial reconstructions under three formulations.
- (c) Keep:
  - MILP adds 313 reactions per genome, of which 84% carry zero flux;
  - positive-cost formulations add only flux-carrying reactions;
  - GSPA's lexicographic pFBA adds 5.7 per genome, against 81 for threshold LP.

  Cut: the active/imbalanced/dead-end counts, the identity check, runtimes, and the FVA/knockout certificate counts (to the Supplement, with one clause in Methods).
- (d) The sign of the cost term causes the zero-flux additions (METEOR's limitation). Parsimony plus attribution is GSPA's contribution.
- (e) Is a parsimonious network still accurate against curated models?

**¶2 – Curated BiGG models (short)**
- (b) Six curated BiGG GEMs.
- (c) Keep one comparison: within the ModelSEED namespace, EC F₁ 0.585 for METEOR vs. 0.584 for Reconstructor. State that CarveMe's higher score reflects that it shares the BiGG namespace. If GSPA's own F₁ on these GEMs is not reported, say so to me rather than in the text. **Decision for you:** keep this paragraph only if it carries a GSPA result; otherwise move it to the Supplement.

**¶3 – METEOR flip filters on 95 genomes (moved from 2.1)**
- (a) Does stoichiometric acceptance turn METEOR's top-1 EC changes into net gains?
- (b) The single-predictor stoichiometric filter and the three-channel consensus filter.
- (c) Keep:
  - unfiltered METEOR changes lose EnzBERT Top-1 at every tier;
  - the consensus filter yields net gains in all 12 predictor × tier combinations.

  Cut: all per-tier counts, p-values and F_max values. Describe the filters in words; their definitions go to Methods and the Supplement. Remove `gspa_single` and `blind_*_accept` entirely.
- (e) Leads to experimental validation in 2.7.

### 2.7 Experimental anchors (Tables 3 and 4, Fig. 4b,d)

**¶1 – Price-149 benchmark**
- (a) Do refinements improve predictions on fitness-validated enzymes?
- (b) Single-predictor vs. multi-predictor refinement.
- (c) Keep:
  - MAPred Top-1 50 → 56 single-predictor and 73 multi-predictor;
  - explicit credit split: pooling gives most of the multi-predictor gain (GraphEC 18 → 48), and Stages 1/2B add more for MAPred, DeepProZyme and CLEAN.

  Cut: all dataset macro/micro F_max values, iteration lists, McNemar p-values and stratum numbers (to Table 3 and the Supplement).
- (d) Drivers 1 and 3 compound.
- (e) Do these gains also appear for whole-organism phenotypes?

**¶2 – Gap-filling objective on BacDive diagnostic enzymes**
- (a) Does the gap-filling objective affect phenotype accuracy?
- (b) MILP on BacDive diagnostic enzymes and urease.
- (c) Keep: MILP lowers precision on diagnostic enzymes (0.644 → 0.619); complex-aware calling on urease gives F₁ 0.76. Cut: the v6/a025/v7 parameterizations (to the Supplement, and rename them there).
- (e) Leads to the blind and supervised tracks.

**¶3 – Blind and supervised tracks**
- (a) Does GSPA match or exceed reconstruction tools without labels, and does it add to supervised classifiers?
- (b) The three BacDive tracks.
- (c) Keep:
  - blind: same binary calls as gapseq, ROC-AUC 0.81 → 0.86 from calibrated scores;
  - genus-blocked hybrid: F₁ 0.857 → 0.871 over the best classifier (Traitar);
  - phenotype-anchored (30% of assays revealed): held-out carbon F₁ 0.79 → 0.89.

  Cut: everything else (to Table 4).
- (d) Calibrated genome-level probabilities add information that binary reconstruction calls lack.
- (e) Some ambiguity remains after deterministic inference, which leads to 2.8.

### 2.8 Driver 3: guarded LM reasoning (Table 7, Fig. 6)

**¶1 – Benchmark of residual cases**
- (a) How many loci stay ambiguous after deterministic inference, and can they be resolved without labels?
- (b) Trigger rate on 88 proteomes; the 120-anomaly benchmark (its construction, masking and decoys move to Methods).
- (c) Keep: 1.8% of proteins are triggered. Cut: the per-genome anomaly counts.
- (e) Leads to how well the agent resolves them.

**¶2 – Agent vs. baselines**
- (b) The agent compared with a threshold rule, supervised classifiers and the guardrail ablation.
- (c) Keep:
  - agent F₁ 0.91 with 0 accepted violations; removing the guardrail admits violations and drops F₁ to 0.85;
  - a supervised classifier given the agent's post-tool features reaches 0.97.

  Cut: the action counts per seed, latency details (to the table) and the per-category breakdown.
- (d) Stated positively: supervised models suffice for binary promotion when labels exist. The guarded agent works zero-shot and is the only configuration that proposes ancestor relaxations (about 14 vs. ≤2).
- (e) Leads to ambiguous metabolic gaps.

**¶3 – Ambiguous gaps and prospective test**
- (b) LRO-50, then CAFA6 (prospective).
- (c) Keep:
  - Hit@1 34% → 83% on LRO-50;
  - CAFA6 macro F_max 0.543 vs. 0.513 for an unguarded agent.

  Cut: per-seed values and the MRR.
- (d) The ontology guardrail makes LM reasoning both safe and useful.

**¶4 – Desert isolates (illustration, two sentences)**
- (c) Keep: pathway coherence 0.49 → 0.62, with candidate osmoprotection and oxidative-stress repairs. Cut: the CDS counts and essential-completeness values.

---

## 5. Discussion (five paragraphs) and Limitations

**D1 – Main claim, framed by the three drivers.** The genome is the right unit of inference. The Results attribute each improvement to a driver: ensemble accuracy (Driver 1), structural invariants (Driver 2), and localization of missing functions (Driver 3). No numbers in this paragraph.

**D2 – Driver 1 relative to GOLabeler/NetGO.**
- GSPA's Stage 0 is in the same multi-source line of work; it accounts for most protein-centric accuracy and the largest single Price-149 step.
- What GSPA adds is a calibrated *probability distribution*, not a ranking, so genome-level priors can update it.
- Rephrase "homology methods are not obsolete" positively: "homology and embedding similarity supply the prior that whole-genome inference updates."

**D3 – Driver 2 relative to genomic-context methods and ontology filters.**
- Local context (operons, gLMs) gives a ranking signal (the subunit result).
- Formal axioms give global consistency (the eukaryotic result).
- Relaxation retains the annotations that naive filters delete.
- Timeliness: with complete assemblies now routine, these constraints apply to most new genomes.

**D4 – Driver 3 relative to METEOR and CarveMe/gapseq.**
- State plainly that the sign of the gap-filling cost determines whether zero-flux reactions are added.
- Parsimony plus competitive responsibility turns a network repair into gene-level evidence.
- Name the advances over METEOR: positive cost, attribution against an unassigned-gene hypothesis, context.

**D5 – When an LM agent is warranted.** Labeled anomalies favor supervised classifiers. New clades and ancestor-relaxation decisions favor the guarded agent. The guardrail is what makes either acceptable.

**Limitations (one ordinary paragraph or two).** Each limitation is phrased as what the design establishes, followed by a named next test. No forbidden negation forms.
1. On well-annotated proteomes, Stage 0 determines most accuracy, and Stages 1–3 act on about 15% of loci. The next test is a benchmark weighted toward MAGs and novel phyla.
2. CAFA6 and Holdout-125 supply isolated sequences, so Driver 2 cannot operate on them. Linking targets to their assemblies would allow the full pipeline to be tested prospectively.
3. Leave-reaction-out recovery is bounded by Stage 0 candidate recall. Embedding backbones saw UniRef50 during pretraining; time-split PLMs would remove that exposure.
4. Leave-subunit-out decoys come from unrelated loci, which favors operon context by design. Same-operon paralog decoys would give a stricter test. *(This limitation is not in the current draft; I recommend adding it.)*
5. Steady-state FBA leaves out regulation, transporters and cofactors. Expression priors and transporter models are the extension.
6. The anatomical bridge covers 2,715 GO terms with logical definitions.
7. Hyperparameters were set on *E. coli* K-12 and frozen; the held-out results support transfer.
8. MAG and isolate predictions are hypotheses for expression or fitness assays.

---

## 6. Methods (restructured: science in the main text, code identifiers in the Supplement)

**Rule for this section:** every formal object keeps its equation and parameter values. Groovy/Python class, file and function names, line numbers, enum labels, CLI flags and TSV paths move to a new Supplement section, "Implementation and provenance map", plus `results/PROVENANCE_MAP.md`.

| # | Subsection | Contents in main text | Moves to Supplement |
|---|---|---|---|
| 4.1 | Problem formulation and notation | Notation for 𝒫, 𝒦, τ, 𝒯, 𝓜, 𝒪 and **Y**; four-stage overview | — |
| 4.2 | Benchmarks and datasets *(new; collects the definitions now scattered in Results)* | One paragraph per data family: (i) CAFA temporal, CAFA6 post-cutoff, Holdout-125, 95-genome EC; (ii) the 9/13/21/88/1,000-genome and 6-eukaryote panels, with selection criteria; (iii) leave-reaction-out (strict/full, strata, Track A/B, All vs. Pool), leave-subunit-out (attenuation range, decoy matching), 108 reconstructions, 6 BiGG GEMs, 106 BGCs, 97 MAGs; (iv) Price-149, BacDive panels, 120-anomaly and LRO-50 construction (symbol masking and hash ordering described in words) | Suite ↔ question ↔ file table (the former `tab:suites`) |
| 4.3 | Stage 0 | Platt calibration equation; reliability scaling; correlation-group maximum; Noisy-OR (Eq. 1); log-odds clipping; ancestor propagation | Per-tool (a_e, b_e) table; reliability-weight table; correlation-group membership table |
| 4.4 | Stage 1 | Damped update (Eq. 2); convergence and rollback rule; six priors, each defined mathematically with α values; RLGC (Eq. 3) | Class names; `theta.json` key list (Table S1 already covers it) |
| 4.5 | Constraint resources | Counts of rules, axioms, classes and properties; extraction sources and releases | Resource file names and the extraction script |
| 4.6 | Stage 2A | Forbidden set = taxon constraints ∪ SAT unsatisfiable core; relaxation (Eq. 4); singleton attenuation | CLI flags, mode strings |
| 4.7 | Stage 2B | Universal network construction; DNF aggregation (Eq. 5); positive-cost lexicographic pFBA (both steps, tolerances); certificate definitions in a single sentence; responsibility (Eq. 6); recoverability indicators | Price-149 reduced protocol (neighbor window, coefficient set, outer damping Eq. 7, multi-predictor blend), named "Price-149 refinement protocol" with one main-text sentence pointing to it |
| 4.8 | Stage 3 | Anomaly categories and triggers; tool functions described in plain words; guardrail resources; model checkpoint, serving, context window, temperature, seeds | Tool schemas, function names |
| 4.9 | BacDive protocols | Track 1 pooling (Eq. 8); genus-blocked CV in words ("deterministic hash of genus name"); hybrid and Track 2 equations | Classifier hyperparameters; hash expressions; feature list |
| 4.10 | Metrics and statistics | F_max (micro and CAFA), WF_max, S_min; homology-masking protocol; CIs and tests; GAEF definition (Eq. 9), presented as an internal diagnostic | — |
| 4.11 | Calibration and compute | All hyperparameters fit on *E. coli* K-12 and frozen; iteration stopping justified by the temporal holdout; hardware | — |

**Data and code availability:** repository, release tag and archive DOI (if one is minted). Point to the Supplement for the file inventory instead of listing scripts in the main text.

**Supplement additions:**
- an identifier → label map (A1–A9 → "Stage 0 baseline … + Stage 2A relaxation"; E1–E5; M1–M4; Config 1–4; METEOR v6/a025/v7/v8; the two acceptance filters);
- the moved per-species, per-iteration and per-tier numbers;
- strict leaf vs. ancestor-propagated scores;
- the Price-149 protocol;
- the provenance table.

---

If you approve the structure and settle the three decisions above (CAFA6 placement, the BiGG paragraph, the suite table going to the Supplement), I'll draft `main.tex` in that order. I'd start with Results and Methods, then write the Introduction and Abstract against the finished Results.
