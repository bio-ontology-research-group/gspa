# Structural and narrative rewrite plan

No files have been edited. This plan covers [main.tex](/home/leechuck/Public/software/gspa/paper/main.tex) and the supporting reorganization of [supplement.tex](/home/leechuck/Public/software/gspa/paper/supplement.tex).

## Global editorial architecture

The rewritten paper should make one argument:

> GSPA improves functional annotation through three cumulative drivers: complementary per-protein evidence establishes a strong baseline; whole-genome architecture adds information unavailable from isolated proteins; and systems-level plus neuro-symbolic reasoning resolves the smaller set of cases that remain ambiguous.

The manuscript should apply these rules throughout:

- Treat homology, domains, protein language models, and structural similarity as the essential starting distribution, not as obsolete competitors.
- Organize prior work around representations and operations, not lists of tools.
- Give every Results paragraph a question, test, finding, interpretation, and bridge.
- Retain at most one principal numerical contrast, or two tightly linked values, per empirical finding in prose.
- Put confidence intervals, per-species values, secondary metrics, iteration traces, and complete comparator results in Tables 1–7, Figures 1–6, or supplementary tables.
- Replace internal configuration labels with scientific names. For example:
  - `A1_noisy_or` → “Stage 0 baseline”
  - `A7_stage1_all` → “all Stage 1 priors”
  - `A9_stage2_rel` → “Stage 2A semantic relaxation”
  - `M3` → “coherence and operon context”
  - `E5` → “semantic ancestor relaxation”
  - `Config 3` → “guarded multi-turn agent”
- Remove filenames, paths, line numbers, function/class names, configuration keys, serialized action labels, and code-specific field names from all main-text prose and captions.
- Move the current benchmark inventory, now Table 8, to the supplement. The main manuscript will then contain the seven scientific result tables referenced in the editorial direction.
- Rewrite all captions as short explanations of the scientific comparison. Dataset provenance belongs in the supplement and `results/PROVENANCE_MAP.md`.
- Keep exact equations and essential parameter values in Methods, but describe them as properties of scientific objects rather than implementation objects.
- Move the Empty Quarter extremophile analysis to a supplementary case study. It is illustrative and lacks experimental validation, so it should not carry the main causal argument.
- Use “Supplementary” only with a specific destination, such as “Supplementary Table S4” or “Supplementary Methods S2.”

---

# Title

Exact proposed title:

> **Genome-scale context, stoichiometric flux modelling, and neuro-symbolic reasoning improve proteome function annotation**

Also update the PDF metadata and supplementary title to this exact sentence-case form.

---

# Abstract

Use exactly five sentences.

1. **Problem and unit of inference.**  
   Explain that sequence, domain, language-model, and structural predictors provide powerful protein-level evidence, but isolated-protein scoring cannot use the pathways, complexes, chromosomal organization, physiology, and taxon-specific anatomy encoded by a complete genome.

2. **Method and three drivers.**  
   Introduce GSPA as a staged framework combining:  
   (i) calibrated complementary per-protein predictors;  
   (ii) whole-genome pathway, complex, operon, taxonomic, and anatomical information; and  
   (iii) stoichiometric flux modelling followed by guarded language-model reasoning for the residual ambiguous loci.

3. **Driver 1 and Driver 2 headline result.**  
   State that the ensemble provides the strong baseline and that genome architecture supplies information beyond that baseline, using only the leave-subunit-out contrast: complex-subunit Hit@1 rises from **56.0% to 93.7%** with coherence and operon context. Point implicitly to the main results rather than enumerating species or secondary metrics.

4. **Driver 3 headline results.**  
   State that positive-cost pFBA reduces the **83.72% zero-flux padding** produced by negative-cost gap filling to **zero**, while guarded reasoning, applied to **fewer than 2.5% of loci**, raises ambiguous reaction-level Hit@1 from **34.0% to 82.7%**. Do not add FVA, knockout, latency, MRR, or seed-specific values here.

5. **Conclusion and scope.**  
   Conclude that accurate protein-level evidence, complete genomic architecture, and selectively applied systems reasoning are complementary levels of functional inference, and state that the software, constraint resources, and reproducible benchmarks are openly available.

Numbers to remove from the Abstract include all per-species eukaryotic changes, Price-149 predictor-by-predictor counts, BacDive track-specific values, FVA/knockout certificate counts, MAG counts, confidence intervals, and iteration results.

---

# 1. Introduction

## Paragraph 1: Functional annotation as a joint genome problem

**Role:** Define the scientific problem and its outputs.

- Begin with protein function annotation as the assignment of GO terms and EC functions to the protein-coding genes of an organism.
- Explain that a complete annotation must satisfy both local evidence and organism-level relationships.
- Introduce pathways, obligate complexes, taxonomic constraints, and anatomy as examples of relationships among annotations.
- Avoid the current five-item invariant catalogue. Use two or three examples in connected prose.
- End with the limitation that conventional evaluation and prediction normally treat proteins as independent instances.

## Paragraph 2: Why whole-genome inference is now timely

**Role:** Establish the technological premise for Driver 2.

- Explain that individual genes could be predicted from fragmented short-read assemblies even when chromosomal architecture was unavailable.
- Contrast that historical setting with long-read sequencing and modern genome/metagenome assemblers, which increasingly recover complete or near-complete chromosomes.
- Cite `amarasinghe2020longread`, `logsdon2020longread`, and `kolmogorov2020metaflye`.
- Explain the consequence: gene order, strand, operon structure, neighborhood conservation, and genome-wide pathway composition are now routinely available inputs.
- Bridge: this expanded input changes functional annotation from isolated classification to structured inference.

## Paragraph 3: Driver 1, complementary per-protein predictors

**Role:** Explain the first conceptual pillar of the state of the art.

- Describe the representations supplied by sequence homology, domain profiles, protein language models, and structural similarity.
- Explain how multi-source ensembles such as GOLabeler and NetGO combine evidence channels and improve over any single representation.
- Cite `you2018golabeler`, `you2019netgo`, and `wang2023netgo3`.
- State positively that homology and similarity provide the strongest initial probability distribution for most well-characterized proteins.
- Identify the remaining gap: all these channels generally score each protein independently of the rest of its proteome.
- Do not list current GO or EC tools individually. Use at most one example per representation class.

## Paragraph 4: Driver 2, genomic context and chromosomal architecture

**Role:** Explain what local genome context adds and where it stops.

- Describe operon prediction, conserved gene neighborhoods, phylogenetic profiles, and genomic language models as representations of local co-occurrence and chromosomal organization.
- Explain why these signals distinguish paralogs and support weakly predicted members of pathways and heteromeric complexes.
- State the boundary positively: these approaches model local or comparative context, while global cellular physiology, complex stoichiometry, and formal taxonomic/anatomical consistency require additional structures.
- Bridge from local genome architecture to whole-cell constraints.

## Paragraph 5: Driver 3, stoichiometric systems biology and METEOR

**Role:** Explain metabolic reconstruction and the specific advance over the authors’ prior work.

- Introduce constraint-based metabolic reconstruction as a representation of biochemical reactions and steady-state mass balance.
- Describe pathway hole filling and reconstruction systems conceptually, mentioning CarveMe and gapseq only as examples.
- Present METEOR as the authors’ initial effort to move enzyme inference beyond homology, accepted for the PSB 2027 “Biological Molecular Function: Beyond Homology” session.
- Explain the METEOR operation precisely: reaction-level indicator optimization with negative addition costs, followed by uniform transfer of reaction support to candidate genes.
- State the two observed consequences motivating GSPA: zero-flux reaction padding and indistinguishable boosting of paralogous candidate loci.
- Cite `niu2026meteor`; update its bibliographic status when final proceedings metadata are available.

## Paragraph 6: Formal ontology constraints and selective language-model reasoning

**Role:** Complete Driver 3 beyond metabolism.

- Explain that GO taxon constraints and GO cross-products with Uberon, PO, and CL encode organismal restrictions that similarity scores cannot guarantee.
- Explain why deletion-based post-filters lose valid ancestral biological-process information.
- Introduce tool-using language models as a way to synthesize heterogeneous evidence at genuinely ambiguous loci.
- State that formal symbolic verification supplies the acceptance criterion for those proposals.
- Avoid a catalogue of agent systems; cite one or two representative examples.

## Paragraph 7: GSPA mapped directly to the three gaps

**Role:** State the contribution.

- Stage 0 implements Driver 1 by calibrating and pooling complementary per-protein evidence.
- Stages 1 and 2A implement Driver 2 by integrating essentiality, pathways, complexes, reaction-local genome context, operon coupling, and organism-wide ontology constraints.
- Stages 1 and 2B implement Driver 3 by combining genomic context with role-aware reaction mechanisms, positive-cost lexicographic pFBA, and competitive locus responsibility against an unassigned-gene hypothesis.
- Stage 3 completes Driver 3 by invoking guarded tool-using reasoning only for residual ambiguities.
- Refer to Figure 1.
- Keep model checkpoint, parameter, and constraint-table counts out of this paragraph.

## Paragraph 8: Evaluation and availability

**Role:** Define what the paper tests.

- Summarize the evaluation questions, not the entire result set:
  1. Does evidence integration improve over its constituent channels?
  2. Does genome context improve pathway and complex annotation while preserving protein-level accuracy?
  3. Does positive-cost flux modelling avoid reaction padding and assign responsibility to loci?
  4. Do ontology constraints preserve valid ancestral annotations?
  5. Does guarded reasoning help on the residual ambiguous cases?
- Mention the breadth of the benchmarks across prokaryotes, viruses, and multicellular eukaryotes without listing all cohorts.
- End with the versioned public repository and release.
- Do not repeat numerical results from the Abstract.

---

# 2. Results

## 2.1 Complementary predictors establish the protein-level baseline

### Paragraph 1

- **(a) Question/Hypothesis:** Does combining complementary per-protein representations improve over any constituent channel?
- **(b) How tested:** Compare the calibrated multi-channel ensemble with sequence, domain, PLM, neural, and structural channels under the same ancestor-propagated protocol on the 21-genome panel.
- **(c) Load-bearing numbers:** Report **micro \(F_{\max}=0.8705\)** for the seven-channel ensemble versus **0.8580** for the strongest individual channel; point to Figure 2d. Cut all other channel-by-channel micro, CAFA, EC, and strict-holdout values.
- **(d) Meaning:** Supports Driver 1: complementary predictors establish a stronger initial distribution than any single representation.
- **(e) Bridge:** Ask whether this advantage persists under temporal and homology-restricted evaluation.

### Paragraph 2

- **(a) Question/Hypothesis:** Does the ensemble remain useful when close sequence transfer is restricted?
- **(b) How tested:** Use the CAFA temporal holdout and its below-30%-identity subset, with donor masking applied to the searchable reference resources.
- **(c) Load-bearing numbers:** Report the twilight-zone change from **0.6723 to 0.6775** after genome-scale refinement; point to Figure 2d and the appropriate supplementary table. The Stage 0 overall score may be stated without further decimals if needed. Cut aspect-by-aspect values and all three iteration traces.
- **(d) Meaning:** Similarity remains the foundation, while contextual inference produces its clearest gains where direct similarity is weakest.
- **(e) Bridge:** Move from isolated targets to complete proteomes, where additional genome-scale signals become available.

### Paragraph 3

- **(a) Question/Hypothesis:** Is the ensemble contribution reproducible across enzyme predictors and masking regimes?
- **(b) How tested:** Summarize the post-cutoff enzyme holdout and 95-genome homology-masked evaluation.
- **(c) Load-bearing numbers:** State only that multi-predictor acceptance produced positive Top-1 gains in **all 12 predictor-by-masking comparisons**; point to Supplementary Tables S3–S4. Cut every per-predictor gain, discordant count, \(p\)-value, Top-3/Top-5 value, and masking-tier sequence.
- **(d) Meaning:** Driver 1 is a general ensemble effect, not a result tied to one protein model.
- **(e) Bridge:** Having established the baseline, test what whole-genome reasoning adds.

## 2.2 Whole-genome priors improve structural coherence while preserving accuracy

### Paragraph 1

- **(a) Question/Hypothesis:** Can genome-scale priors improve proteome structure without sacrificing the strong Stage 0 baseline?
- **(b) How tested:** Apply the full ablation sequence to nine reference genomes, with one development genome and eight out-of-sample genomes.
- **(c) Load-bearing numbers:** Report CAFA \(F_{\max}\) as **0.8041 at Stage 0 and 0.8037 after full refinement**; point to Table 1 and Figure 2a–b. Use this pair to establish preservation. Cut standard deviations, experimental-only variants, Wilcoxon values, and per-genome changes from prose.
- **(d) Meaning:** Drivers 2 and 3 refine structural consistency around a baseline established by Driver 1 instead of replacing it.
- **(e) Bridge:** Examine which biological structures improve.

### Paragraph 2

- **(a) Question/Hypothesis:** Which genome-level invariants respond to whole-proteome inference?
- **(b) How tested:** Ablate essentiality, pathway/complex coherence, gap filling, operon context, and ontology constraints.
- **(c) Load-bearing numbers:** Report the composite structural score increase from **0.8459 to 0.8643** and point to Table 1 and Figure 2a. Leave pathway, complex, and essentiality components in the table.
- **(d) Meaning:** The combined priors add complementary structural information, supporting Driver 2.
- **(e) Bridge:** Structural improvement must also withstand deliberately introduced cross-clade errors.

### Paragraph 3

- **(a) Question/Hypothesis:** Do formal genome-wide constraints suppress high-scoring cross-clade transfers?
- **(b) How tested:** Inject cross-clade developmental and taxon-restricted claims, then compare calibrated Stage 0 with the constrained stages.
- **(c) Load-bearing numbers:** Report suppression of **122 of 125** active post-calibration violations; point to Table 1 and Figure 2c. Cut raw-score ranges, term-by-term examples, the 897/900 alternative denominator, and direct-versus-cross-product subdivisions.
- **(d) Meaning:** Driver 2 adds organism-level information that no protein-local score can guarantee.
- **(e) Bridge:** Test the same principle at population scale.

### Paragraph 4

- **(a) Question/Hypothesis:** Do the structural gains scale beyond the small reference panel?
- **(b) How tested:** Audit 1,000 bacterial genomes and the four-kingdom cohort using the same frozen settings.
- **(c) Load-bearing numbers:** Report MetaCyc coherence increasing from **36.71% to 66.75%**, with **zero surviving taxonomic violations per genome**; point to Figure 2d and a supplementary benchmark table. Cut the six-comparator catalogue, domain-by-domain confidence intervals, runtime ranges, memory, and iteration-specific values.
- **(d) Meaning:** Whole-genome inference scales and improves the properties it was designed to enforce.
- **(e) Bridge:** The next analysis isolates the mechanistic contribution of local genomic architecture.

## 2.3 Genomic architecture resolves metabolic gaps

### Paragraph 1

- **(a) Question/Hypothesis:** Can reaction-local genome context recover a masked enzyme when sequence ranking is insufficient?
- **(b) How tested:** Mask known gene–reaction links in the strict leave-reaction-out benchmark and rank candidate loci using sequence-only evidence, genomic context, and combined models.
- **(c) Load-bearing numbers:** Report Hit@1 increasing from **1.95% for sequence-only ranking to 23.05% with reaction-local genome context**; point to Table 2 and Figure 3a. Cut confidence intervals, pool-conditioned metrics, MRR, Hit@3, and intermediate rerankers.
- **(d) Meaning:** Direct evidence for Driver 2: chromosomal neighborhood identifies pathway loci that isolated sequence scores cannot resolve.
- **(e) Bridge:** Determine how genomic context interacts with a strong PLM representation.

### Paragraph 2

- **(a) Question/Hypothesis:** Are protein-level similarity and genome context substitutes or complements?
- **(b) How tested:** Compare PLM reaction-centroid similarity alone with the combined PLM, pathway-coherence, and genomic-context model.
- **(c) Load-bearing numbers:** Report Hit@1 of **62.11% for PLM similarity and 63.67% for the combined model**; point to Table 2 and Figure 3a. Cut the full candidate-pool intersection analysis and all secondary ranking metrics.
- **(d) Meaning:** Driver 1 carries most cases in a well-annotated organism, while Driver 2 resolves ties and paralog ambiguity.
- **(e) Bridge:** Test whether the combined inference transfers across species and biological roles.

### Paragraph 3

- **(a) Question/Hypothesis:** Does the combined model generalize beyond the development organism?
- **(b) How tested:** Apply frozen settings to the *B. subtilis* leave-reaction-out panel and evaluate across core, catabolic, and peripheral reaction strata.
- **(c) Load-bearing numbers:** Report **43.06% unconditional Hit@1** in the out-of-sample species; point to Table 2 and Figure 3b–c. Leave stratum-specific values, pool-conditioned performance, and all confidence intervals in the displays.
- **(d) Meaning:** Genome context and protein representations remain complementary across species.
- **(e) Bridge:** Apply the same inference to incompletely characterized genomes.

### Paragraph 4

- **(a) Question/Hypothesis:** Can the method generate experimentally testable annotations in uncultured genomes?
- **(b) How tested:** Apply the fixed model to 97 high-quality metagenome-assembled genomes.
- **(c) Load-bearing numbers:** Report **568 high-confidence candidate assignments** and pathway completeness increasing from **0.511 to 0.656**; point to Figure 3d. Cut the lower-confidence count, per-MAG means, phylum lists, and threshold catalogue.
- **(d) Meaning:** Drivers 1 and 2 produce ranked hypotheses where curated annotations are sparse.
- **(e) Bridge:** Candidate pathway completion still needs a physiological criterion that distinguishes useful reactions from network padding.

The antiSMASH specificity result should move to a short supplementary Results paragraph unless it is needed by the target venue. Its full three-metric recital currently interrupts the metabolic argument.

## 2.4 Positive-cost flux modelling assigns metabolic responsibility to loci

### Paragraph 1

- **(a) Question/Hypothesis:** Does negative-cost reaction selection add reactions that participate in the inferred flux state?
- **(b) How tested:** Audit reaction additions across 108 reconstructions under negative-cost MILP, positive-cost repair, and lexicographic pFBA.
- **(c) Load-bearing numbers:** Report that **83.72%** of negative-cost additions carry zero flux, compared with **0%** under GSPA’s positive-cost pFBA retention criterion; point to Figure 6a and Supplementary Table S5. Keep counts of added reactions, FVA certificates, knockout certificates, dead ends, and runtime out of prose.
- **(d) Meaning:** Driver 3 converts gap filling from reward-driven reaction accumulation into physiologically active repair.
- **(e) Bridge:** Active reactions must then be assigned to the most plausible gene locus.

### Paragraph 2

- **(a) Question/Hypothesis:** Does competitive locus responsibility improve enzyme annotation?
- **(b) How tested:** Evaluate single-predictor and multi-predictor refinement on the Price-149 fitness-supported enzyme benchmark.
- **(c) Load-bearing numbers:** Use one representative single-predictor result, **MAPred Top-1 from 50 to 56**, and the largest ensemble \(F_{\max}\) gain, **+0.1987**; point to Table 3 and Figure 6b. Cut all per-iteration values, macro/micro/per-protein variants, predictor-by-predictor counts, and \(p\)-values.
- **(d) Meaning:** Role-aware reaction aggregation and competitive gene responsibility improve locus attribution, supporting Driver 3.
- **(e) Bridge:** Test whether those annotations predict observed whole-organism phenotypes.

The six curated BiGG-model comparison should remain in Supplementary Table S5 and receive at most one sentence here, because database namespace compatibility complicates its interpretation.

## 2.5 Systems-level annotations reproduce experimental phenotypes

### Paragraph 1

- **(a) Question/Hypothesis:** Do the inferred pathway probabilities correspond to experimentally measured growth?
- **(b) How tested:** Compare blind carbon-source predictions with reconstruction pipelines on BacDive assays.
- **(c) Load-bearing numbers:** State that the continuous GSPA score preserves the leading binary \(F_1\) of **0.7952** while reaching **ROC-AUC 0.8567**; point to Table 4 and Figure 6d. Cut precision, recall, balanced accuracy, runtimes, and all comparator rows from prose.
- **(d) Meaning:** Driver 3 yields calibrated physiological confidence while retaining the established binary decisions.
- **(e) Bridge:** Ask whether systems probabilities complement supervised phenotype models across unseen genera.

### Paragraph 2

- **(a) Question/Hypothesis:** Does combining mechanistic reconstruction with supervised genomic phenotype prediction improve out-of-genus transfer?
- **(b) How tested:** Use genus-blocked cross-validation for carbon and diagnostic-enzyme assays.
- **(c) Load-bearing numbers:** Report hybrid \(F_1\) of **0.8710 for carbon growth** and **0.8713 for enzyme activity**; point to Table 4. Cut the three baseline \(F_1\) values, per-fold means, standard deviations, and other metrics.
- **(d) Meaning:** Mechanistic and statistical representations provide complementary information.
- **(e) Bridge:** Evaluate whether limited phenotype observations can personalize the organism-level model.

### Paragraph 3

- **(a) Question/Hypothesis:** Can observed activities in one part of a strain’s phenotype profile improve predictions for its unobserved activities?
- **(b) How tested:** Reveal 30% of assays and evaluate only the held-out 70%.
- **(c) Load-bearing numbers:** Report \(F_1\) of **0.8871 for held-out carbon sources** and **0.8545 for held-out enzyme activities**; point to Table 4. Keep the hybrid results and all secondary metrics in the table.
- **(d) Meaning:** Whole-organism observations constrain the remaining pathway hypotheses, extending Driver 3 from genotype to physiology.
- **(e) Bridge:** Metabolic pathways are one form of coordinated function; obligate protein complexes provide an independent test of genomic architecture.

The current paragraph on negative-cost MILP across ten diagnostic ECs should be compressed to one setup sentence or moved to the supplement. Its precision/recall sequences duplicate the stoichiometric audit.

## 2.6 Genome context completes obligate heteromeric complexes

### Paragraph 1

- **(a) Question/Hypothesis:** Can complex membership and operon structure recover weakly supported subunits while rejecting matched decoys?
- **(b) How tested:** Attenuate one subunit in each leave-subunit-out trial and rank it against five matched non-complex proteins.
- **(c) Load-bearing numbers:** Report Hit@1 increasing from **0.560 to 0.937** after adding complex coherence and operon context; point to Table 5 and Figure 4a. Cut recall, precision, \(F_1\), confidence intervals, and per-species results.
- **(d) Meaning:** Driver 2 resolves paralogs by locating protein evidence within a chromosomal assembly.
- **(e) Bridge:** Some complexes are dispersed and therefore remain unresolved by operon proximity.

### Paragraph 2

- **(a) Question/Hypothesis:** Can selective reasoning help when valid complex partners occupy distant loci?
- **(b) How tested:** Analyze the dispersed photosynthetic and respiratory complexes of *Synechocystis* after deterministic genome-context inference.
- **(c) Load-bearing numbers:** Report Hit@1 increasing from **0.640 to 0.680** after guarded Stage 3 reasoning; point to Table 5 and Figure 4b. Cut all other per-species values.
- **(d) Meaning:** Driver 3 complements local architecture in the specific cases where physical proximity is uninformative.
- **(e) Bridge:** Whole-organism context also constrains which biological processes are possible in a clade.

## 2.7 Cross-ontology reasoning preserves valid ancestral biology

### Paragraph 1

- **(a) Question/Hypothesis:** Do unconstrained similarity-based predictions create taxonomic or anatomical contradictions in unperturbed eukaryotic proteomes?
- **(b) How tested:** Evaluate six eukaryotes against direct GO constraints and GO cross-products with Uberon, PO, and CL.
- **(c) Load-bearing numbers:** Report **1,009 endogenous violations reduced to zero**, alongside the mean CAFA \(F_{\max}\) gain of **+0.0021**; point to Table 6 and Figure 4c. Cut every species-specific violation count and \(F_{\max}\) change.
- **(d) Meaning:** Driver 2 improves formal organismal consistency while preserving predictive accuracy.
- **(e) Bridge:** Determine why direct GO constraints and soft penalties are insufficient.

### Paragraph 2

- **(a) Question/Hypothesis:** Which representation is needed to enforce anatomy-linked restrictions?
- **(b) How tested:** Compare direct GO taxon constraints, soft cross-ontology penalties, and deterministic cross-ontology enforcement under cross-clade perturbation.
- **(c) Load-bearing numbers:** Report that direct GO constraints leave **642 violations**, while deterministic GO–anatomy enforcement leaves **zero**; point to Table 6 and Figure 4c. Cut species-specific counts and the weighted/unweighted suppression alternatives.
- **(d) Meaning:** Anatomical restrictions become available only after connecting GO processes to the structures they entail.
- **(e) Bridge:** Once an invalid specific term is found, the remaining question is how much valid meaning can be retained.

### Paragraph 3

- **(a) Question/Hypothesis:** Does semantic ancestor relaxation preserve useful information better than deleting an invalid term?
- **(b) How tested:** Compare removal and nearest-valid-ancestor relaxation against unchanged Swiss-Prot annotations.
- **(c) Load-bearing numbers:** Report an \(F_{\max}\) advantage of **+0.0006 to +0.0157** across multicellular animals and vascular plants; point to Table 6 and Figure 4d. Move all information-content values, recall differences, yeast/moss exceptions, and species-by-species scores into the table.
- **(d) Meaning:** Driver 3 applies ontology semantics constructively by preserving the valid ancestral component of a transferred annotation.
- **(e) Bridge:** Deterministic inference resolves most loci; the final analysis concerns the small residual set.

## 2.8 Guarded language-model reasoning resolves residual ambiguity

### Paragraph 1

- **(a) Question/Hypothesis:** Is language-model reasoning sufficiently selective to be a final-stage method?
- **(b) How tested:** Apply deterministic anomaly triggers across the four-kingdom cohort.
- **(c) Load-bearing numbers:** Report a mean trigger rate of **1.84%**, remaining **below 2.5% in every domain**; point to Figure 5d. Cut the interquartile range and per-domain values.
- **(d) Meaning:** Stage 3 is a residual decision layer, not the source of bulk annotation accuracy.
- **(e) Bridge:** Evaluate whether tool use and the symbolic guardrail improve those selected cases.

### Paragraph 2

- **(a) Question/Hypothesis:** Does guarded multi-turn reasoning resolve anomalies while maintaining formal validity?
- **(b) How tested:** Use the de-identified 120-anomaly benchmark with matched decoys and cross-clade conflicts across three seeds.
- **(c) Load-bearing numbers:** Report mean \(F_1=\)**0.9115** with **zero accepted ontology violations**; point to Table 7 and Figure 5a–b. Cut seed-specific values, action counts, tool-call counts, and latency.
- **(d) Meaning:** Driver 3 combines flexible evidence synthesis with deterministic acceptance semantics.
- **(e) Bridge:** Test candidate ranking on an independent ambiguity task.

### Paragraph 3

- **(a) Question/Hypothesis:** Can the agent distinguish candidate loci that remain tied after deterministic reaction-local scoring?
- **(b) How tested:** Re-rank five candidates in 50 de-identified ambiguous reaction gaps.
- **(c) Load-bearing numbers:** Report Hit@1 increasing from **34.0% to 82.7%**; point to Table 7 and Figure 5c. Cut Hit@3, MRR, seed-specific results, and McNemar details.
- **(d) Meaning:** Tool-mediated reasoning adds value specifically after Drivers 1 and 2 have narrowed the candidate set.
- **(e) Bridge:** Compare this flexible zero-shot method with a cheaper supervised alternative.

### Paragraph 4

- **(a) Question/Hypothesis:** When is an agent preferable to a tabular classifier?
- **(b) How tested:** Compare pre-tool classifiers, post-tool supervised classifiers, fast one-step reasoning, and guarded multi-turn reasoning.
- **(c) Load-bearing numbers:** Report the best supervised upper-bound \(F_1\) of **0.9720** versus **0.9115** for zero-shot guarded reasoning; point to Table 7. Leave all classifier-specific values and runtimes in the table.
- **(d) Meaning:** Labeled in-distribution cases support a cheaper classifier; zero-shot guarded reasoning is useful when labels are unavailable or the output requires structured ancestor relaxation.
- **(e) Bridge:** This tradeoff motivates the limitations and deployment guidance in the Discussion.

---

# 3. Discussion

## Paragraph 1: Synthesis around the three drivers

- Restate the principal conclusion without repeating a metric list.
- Explain the cumulative dependency:
  1. per-protein evidence supplies candidate functions;
  2. genome architecture relates those candidates across loci;
  3. systems and symbolic models select assignments compatible with physiology and organismal biology.
- Emphasize that no stage independently accounts for all gains.

## Paragraph 2: Driver 1 in relation to GOLabeler and NetGO

- Position Stage 0 within the multi-source protein-ensemble lineage established by GOLabeler and NetGO.
- State that the present results reinforce the continuing value of homology, domains, PLMs, and structure.
- Explain that GSPA changes the subsequent unit of inference, from one protein to an organism-wide annotation state.
- Avoid claims that GSPA supersedes protein predictors.

## Paragraph 3: Driver 2 and the transition to complete genomes

- Relate operon coupling, reaction-local context, essentiality, and complex coherence to comparative-genomics and genomic-language-model approaches.
- Explain that long-read assemblies make these signals widely available.
- Distinguish local genomic organization from whole-genome constraints: local neighborhoods identify likely loci, while the proteome-wide model tests pathway, complex, taxonomic, and anatomical compatibility.
- Interpret the LRO and LSO results as evidence for this shift in inference scale.

## Paragraph 4: Driver 3 and the advance beyond METEOR

- Present METEOR as the direct conceptual predecessor.
- Explain the three advances:
  1. reaction-local genomic context and co-directional operon coupling;
  2. role-aware DNF complex/isozyme aggregation, positive-cost lexicographic pFBA, and competitive locus responsibility;
  3. guarded reasoning on residual cases.
- Use the zero-flux-padding result as the single headline comparison.
- Explain that explicit unassigned-gene responsibility prevents physiological necessity from being converted automatically into support for every paralog.

## Paragraph 5: Ontology semantics and guarded reasoning

- Explain that ontology constraints act as generative correction rules through semantic relaxation, not only as rejection filters.
- Relate this to the preservation of ancestral developmental functions.
- Explain why the symbolic guardrail and language model have different roles: the model proposes a synthesis; the ontology determines admissibility.
- State that the small trigger rate keeps this expensive reasoning focused.

## Paragraph 6: What the evaluation establishes

- State positively that the study supports improved annotation coherence, enzyme ranking, complex recovery, formal consistency, and agreement with bacterial phenotypes across the tested benchmarks.
- Distinguish experimentally anchored results from computational hypotheses in MAGs and extremophile isolates.
- End with the practical implication: complete genomes allow annotation systems to reason over cellular organization, not only protein resemblance.

---

# 3.1 Limitations

Frame each limitation as an established boundary followed by a concrete next test.

## Paragraph 1: Dependence on Stage 0 candidate recall

- State that genome-scale inference can rerank or validate functions represented in its candidate set.
- Note that uncovered true enzymes bound unconditional LRO recovery.
- Propose remote structural search, improved fold-level candidate generation, and time-split PLM evaluation as the next tests.

## Paragraph 2: Availability and quality of genomic architecture

- State that operon and neighborhood evidence is strongest in complete or near-complete assemblies.
- Explain that fragmented MAGs supply only partial local context.
- Propose stratifying future evaluations by assembly continuity and explicitly propagating uncertainty from contig boundaries.

## Paragraph 3: Limits of stoichiometric physiology

- State that steady-state flux captures mass-balanced reaction feasibility under specified media.
- Name transcriptional regulation, transporter specificity, enzyme kinetics, cofactors, and environmental variation as additional determinants.
- Propose transcriptomic priors, transporter models, and multi-condition flux experiments.

## Paragraph 4: Coverage of cross-ontology axioms

- State that anatomical enforcement is limited to GO terms with available logical links to Uberon, PO, or CL and to the taxonomic constraints represented in those resources.
- Propose expanding logical definitions and evaluating additional clades.
- Note that improvements were strongest where ground truth retains tissue-level ancestral processes.

## Paragraph 5: Stage 3 evaluation and deployment

- State that the agent was evaluated on selected, de-identified ambiguity sets and a small number of random seeds.
- Explain that post-tool classifiers are preferable for stable, labeled, in-distribution binary decisions.
- Propose prospective expert adjudication, clade-held-out evaluation, calibration of abstention, and measurement of compute cost.
- Keep the comparison positive; do not describe the agent as universally superior.

## Paragraph 6: Experimental status of novel predictions

- State that MAG and extremophile outputs are ranked computational hypotheses.
- Name heterologous expression, targeted knockout, complementation, and substrate-specific growth experiments as validation routes.
- This paragraph replaces disclaimer-style language with the next experiment required.

---

# 4. Methods

## 4.1 Problem definition and staged inference

### Paragraph 1

- Define proteins, genomic coordinates, the target taxon, GO terms, metabolic reactions, and cross-ontology axioms.
- Define the output as a proteome-wide posterior over protein–function assignments.
- Keep one compact notation paragraph and Figure 1.
- Move the current eight-suite inventory table to Supplementary Table S1.

### Paragraph 2

- Describe the dependency among stages:
  - Stage 0 produces calibrated protein-level probabilities.
  - Stage 1 updates them from genome-wide evidence.
  - Stage 2A enforces symbolic consistency.
  - Stage 2B tests metabolic assignments through flux and locus responsibility.
  - Stage 3 examines residual anomalies.
- State which stages are optional and which data they require.

## 4.2 Stage 0: calibrated complementary protein evidence

### Paragraph 1

- Define the evidence channels by representation: sequence similarity, domain profiles, protein-language-model embeddings, learned sequence predictors, and structural similarity.
- State that GSPA may accept either one predictor or several predictors.
- Cite public tools only as examples.

### Paragraph 2

- Define channel-specific calibration and reliability weighting.
- Give the calibration equation and clipping interval.
- Keep exact fitted/default values in Supplementary Methods S2 unless they determine a reported result.

### Paragraph 3

- Explain correlation groups as scientific dependence classes.
- Define within-group maximum aggregation and across-group Noisy-OR.
- Retain the Noisy-OR equation in the main Methods.
- Explain that the result is the anchored Stage 0 log-odds matrix.

**Move to supplement:** exact evidence-type enum mapping, calibration lookup names, every default coefficient, class/function names, configuration files, and source locations.

## 4.3 Stage 1: whole-proteome log-odds refinement

### Paragraph 1

- Define the iterative update, its anchored Stage 0 term, the six prior contributions, damping, convergence, and rollback.
- Keep the principal update equation.
- Explain why anchoring prevents recursive context from replacing the protein evidence.

### Paragraph 2

- Define essentiality and pathway/complex coherence as organism-wide coverage objects.
- Explain how missing essential categories, incomplete pathways, and incomplete obligate complexes generate candidate-specific log-odds changes.
- Describe singleton attenuation scientifically, without implementation labels.

### Paragraph 3

- Define taxonomic and anatomical priors.
- Explain the distinction between a soft log-odds penalty in Stage 1 and deterministic enforcement in Stage 2A.

### Paragraph 4

- Define the metabolic-gap prior and reaction-local genomic context.
- Introduce the reaction neighborhood, anchor genes, genomic distance kernel, strand agreement, and paralog commitment.
- Retain the RLGC equation if space allows; move its full coefficient inventory to Supplementary Methods S3.

### Paragraph 5

- Define co-directional operon coupling for pathway steps and complex subunits.
- Explain how the operon prior differs from the broader RLGC density field.

**Move to supplement:** all internal prior module names, filenames, line references, configuration keys, fixed thresholds, and exact implementation crosswalk.

## 4.4 Cross-ontology taxonomic and anatomical bridge

### Paragraph 1

- Describe the source ontology releases and versions.
- Define the logical chain from GO process to anatomical class to permitted taxon.
- Report the sizes of the resulting constraint resources only once, preferably in Supplementary Table S3.

### Paragraph 2

- Describe extraction, closure, obsolete-term handling, and quality-control checks.
- Give release dates, versions, and access dates required for reproduction.

**Move to supplement:** extraction script, generated resource filenames, repository paths, row-level schema, and checksums or manifest references.

## 4.5 Stage 2A: deterministic consistency and semantic relaxation

### Paragraph 1

- Define the forbidden-term set produced by direct taxon constraints and entailed anatomical restrictions.
- Explain how the system evaluates both a term and its relevant ancestors.

### Paragraph 2

- Define the nearest valid non-root ancestor set and the relaxation update.
- Keep the semantic-relaxation equation.
- Explain that removal is used when no valid non-root ancestor exists.

### Paragraph 3

- Define deterministic pruning of an isolated obligate-complex singleton.
- State the scientific rule and threshold, with the full parameter source in the supplement.

**Move to supplement:** internal modes, provenance-state fields, report labels, command-line switches, and implementation classes.

## 4.6 Stage 2B: role-aware reaction support

### Paragraph 1

- Define the sanitized universal network, balanced reactions, media, and reaction mechanisms.
- Distinguish isozyme disjunctions from multi-subunit conjunctions.

### Paragraph 2

- Define top-candidate isozyme aggregation, the bottleneck gate for required roles, and maximum aggregation across alternative mechanisms.
- Retain the DNF reaction-support equation.

**Move to supplement:** database import paths, internal data structures, exact sanitation workflow, and complete universal-network statistics.

## 4.7 Stage 2B: positive-cost lexicographic pFBA

### Paragraph 1

- Define core and candidate reactions and strictly positive addition costs.
- Explain why positive costs select only reactions needed for the required biomass flux.

### Paragraph 2

- Define the two lexicographic steps:
  1. minimize repair cost while meeting biomass;
  2. minimize total flux within the optimal repair-cost tolerance.
- Define the retained repair set as non-core reactions carrying pFBA flux.
- Keep tolerances needed to reproduce the reported zero-padding claim.

### Paragraph 3

- Define the FVA, knockout, and perturbation-stability certificates.
- Put their exact thresholds in the main Methods only if they are reported in a main table; otherwise place them in Supplementary Methods S4.

## 4.8 Stage 2B: competitive gene-locus responsibility

### Paragraph 1

- Define candidate loci for each required reaction role.
- Introduce the unassigned-gene hypothesis as a competing explanation.
- Define responsibility from sequence odds, genomic context, synteny/coherence, and PLM similarity.
- Retain the responsibility equation.

### Paragraph 2

- Explain how responsibility and flux certificates become a damped log-odds update.
- Describe condition-specific substrate probing for catabolic reactions.
- Keep the general scientific formulation in the main text.

**Move to supplement:** Price-149-specific formulas, neighbor-window sizes, evaluation-only coefficients, internal “single” and “multi” acceptance rules, and exact implementation paths.

## 4.9 Stage 3: selective guarded reasoning

### Paragraph 1

- Define the four anomaly classes in scientific language: pathway gaps, incomplete complexes, missing essential categories, and borderline ontology conflicts.
- Define the trigger intervals and selection cap.

### Paragraph 2

- Identify the public model checkpoint, decoding settings, maximum turns, and seeds needed for reproduction.
- Describe the five tool capabilities by what evidence they return, not by function names.

### Paragraph 3

- Define the permitted resolution outcomes in prose: promote, relax to a valid ancestor, demote, or abstain.
- Explain that every proposed result is checked deterministically against the ontology and taxon constraints before acceptance.

**Move to supplement:** tool schema names, JSON structures, action enum strings, serving alias, host name, implementation script, prompt text, and output schema.

## 4.10 Benchmark design and leakage controls

### Paragraph 1

- Give the development/evaluation split: *E. coli* K-12 for parameter calibration and frozen settings elsewhere.
- Summarize the eight suites through their scientific questions.
- Point to Supplementary Table S1 for the full cohort and stage matrix.

### Paragraph 2

- Define donor-database homology masking across sequence, embedding, and structural resources.
- Distinguish the conservative cosine-masked LRO track from the track retaining pretrained PLM reference centroids.
- State pretrained-model exposure as a design property of the latter track.

### Paragraph 3

- Define the CAFA temporal, prospective/post-cutoff, enzyme holdout, and multi-genome evaluations.
- Put complete cohort compositions and data-file mappings in the supplement.

## 4.11 LRO, LSO, Price-149, and discovery cohorts

### Paragraph 1

- Define strict and role-stratified leave-reaction-out evaluation, candidate-pool coverage, and unconditional versus coverage-conditioned metrics.

### Paragraph 2

- Define leave-subunit-out construction, attenuation, matched decoys, and zero-leakage formatting.

### Paragraph 3

- Define Price-149 host mapping, single-predictor versus multi-predictor refinement, and the experimental fitness-supported endpoint.

### Paragraph 4

- Describe the antiSMASH and MAG cohorts as specificity and hypothesis-generation applications.
- Move full selection criteria and per-cohort statistics to the supplement.

## 4.12 BacDive phenotype protocols

### Paragraph 1

- Define the carbon-source and diagnostic-enzyme panels and the three reconstruction inputs.
- Describe blind probability pooling and cross-assay coupling.

### Paragraph 2

- Define genus-blocked cross-validation and ensure that calibration, threshold selection, and classifier fitting occur within training folds.
- Give the hybrid operation conceptually.

### Paragraph 3

- Define the 30%-revealed/70%-held-out design and strain-specific Bayesian anchoring.
- Keep complete feature vectors and estimator hyperparameters in Supplementary Methods S6.

## 4.13 Evaluation metrics and statistical analysis

### Paragraph 1

- Define \(F_{\max}\), weighted \(F_{\max}\), micro \(F_{\max}\), \(S_{\min}\), Hit@\(k\), MRR, and phenotype metrics by the tasks where they apply.
- Avoid mixing protein-centric, pair-level, and per-protein quantities under one unqualified “\(F_{\max}\).”

### Paragraph 2

- State confidence-interval procedures, paired tests, bootstrap units, multiple-run summaries, and stopping rules.

### Paragraph 3

- Define the genome-scale structural audit and its components.
- State explicitly and positively that it measures targeted annotation invariants and is reported alongside held-out experimental evaluations.
- Retain its equation; move implementation details to the supplement.

---

# Supporting changes to `paper/supplement.tex`

## Proposed supplementary structure

1. **S1. Benchmark suites and manuscript provenance**
   - Move the current main-text Table 8 here as Supplementary Table S1.
   - Include suite, cohort, input channels, active stages, comparators, canonical result files, and generator workflow.
   - Add a main-label-to-internal-code crosswalk for every renamed configuration.

2. **S2. Stage 0 calibration and hyperparameters**
   - Retain all calibration coefficients, correlation-group membership, reliability weights, clipping values, and configuration keys.
   - Include the exact code files and public commands needed to reproduce Stage 0.

3. **S3. Stage 1 and ontology-resource implementation**
   - Retain the 18-parameter table.
   - Add exact prior-module mappings, resource filenames, ontology versions, extraction scripts, and validation checks.

4. **S4. Stoichiometric refinement implementation**
   - Retain full reaction sanitation statistics, pFBA parameters, certificate thresholds, responsibility coefficients, and Price-149-specific evaluation formulas.
   - Preserve exact Python/Groovy source mapping here.

5. **S5. Extended protein and genome benchmark results**
   - Keep Holdout-125, 95-genome masking, CAFA aspect-specific results, four-kingdom breakdowns, confidence intervals, and significance tests.

6. **S6. Extended metabolic and BacDive results**
   - Keep the curated BiGG comparison, complete stoichiometric audit, role-stratified Price-149 results, diagnostic-EC analysis, and per-fold BacDive metrics.

7. **S7. Stage 3 prompts, tools, guardrail, and multi-seed results**
   - Preserve exact tool names, schemas, prompt template, action codes, seeds, model-serving details, classifier features, latency, and per-seed results.

8. **S8. MAG and Empty Quarter case studies**
   - Move the extremophile narrative from the main Results.
   - Clearly label MAG and isolate assignments as experimentally testable predictions.

## Provenance requirements

Update `results/PROVENANCE_MAP.md` so that it maps:

- each human-readable main-text method label to its exact internal configuration;
- every Table 1–7 row and Figure 1–6 panel to the canonical TSV columns;
- every supplementary table to its code, configuration, and data files;
- all main-text headline numbers to exact source rows and columns.

The main manuscript’s Data and Code Availability sections should name the versioned release, repository, archive/DOI, and provenance manifest without listing internal paths. Exact paths remain in Supplementary Methods S1 and `PROVENANCE_MAP.md`.

---

# Bibliography work required during implementation

The current bibliography contains `you2019netgo`, `wang2023netgo3`, and `niu2026meteor`, but the requested keys for GOLabeler and long-read/genome-assembly context are absent. The rewrite should add and verify:

- `you2018golabeler`
- `amarasinghe2020longread`
- `logsdon2020longread`
- `kolmogorov2020metaflye`

The METEOR entry should be updated from its present preprint-only description to the accepted PSB 2027 record once authoritative proceedings metadata are available.