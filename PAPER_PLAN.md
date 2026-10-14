# GSPA Publication & Research Plan
**Target Venues:** *Nature Methods* / *Nucleic Acids Research* (Methods) / *Bioinformatics*  
**Working Title:** *Genome-scale probabilistic inference, ontological constraint optimization, and neuro-symbolic agentic completion improve protein function prediction across the tree of life*

---

## 1. Central Research Question & Thesis

### Research Question
> **Can explicit genome-scale biological context and formal ontological constraints—spanning metabolic networks, multi-subunit protein complexes, operon/chromosomal organization, essential cellular repertoires, and taxonomic as well as anatomical/developmental structure boundaries—measurably improve both protein-level function prediction accuracy (CAFA metrics) and system-level biological plausibility (GAEF metrics) over protein-in-isolation predictors, and how much additional recovery can an optional, neuro-symbolically guarded agentic layer achieve on unresolved genome-scale anomalies?**

### Central Thesis
All standard protein function predictors (sequence homology, domain scanners, protein language models, structure GNNs) and even standard genome annotation pipelines (NCBI PGAP, Bakta, Prokka, eggNOG-mapper, InterProScan) operate under a **reductionist independence assumption**:
$$P(\mathbf{Y}_{1:N} \mid \mathbf{X}_{1:N}) = \prod_{i=1}^N P(Y_i \mid X_i)$$
where each protein $i \in \{1,\dots,N\}$ in a genome is annotated in isolation. While Molecular Function (MF) often depends on a single sequence/structure fold, **Biological Process (BP)** and **Cellular Component (CC)** are emergent properties of interacting protein collectives within a specific organismal body plan and lineage.

**GSPA solves this problem** through a **three-stage genome-scale inference architecture** (with Stage 3 as an optional, separately evaluated layer):
1. **Stage 1 — Soft Bayesian Context Integration (`IterativeRefiner` + `ReactionLocalContextSuggester`):**
   $$P(\mathbf{Y}_{1:N} \mid \mathbf{X}_{1:N}, \mathcal{G}, \mathcal{O}) \propto \prod_{i=1}^N P_{\text{base}}(Y_i \mid X_i) \cdot \Psi_{\text{taxon}}(\mathbf{Y}, T_{\mathcal{G}}) \cdot \Psi_{\text{anat}}(\mathbf{Y}, \mathcal{A}_{\mathcal{G}}) \cdot \Psi_{\text{pwy}}(\mathbf{Y}, \mathcal{M}_{\mathcal{G}}) \cdot \Psi_{\text{cplx}}(\mathbf{Y}, \mathcal{C}) \cdot \Psi_{\text{ctx}}(\mathbf{Y}, \mathcal{L}_{\mathcal{G}}) \cdot \Psi_{\text{ess}}(\mathbf{Y}, \mathcal{R}_{T_{\mathcal{G}}})$$
   Combines calibrated multi-modal protein evidence via correlation-group Noisy-OR and iterative Jacobi-damped Bayesian priors, plus reaction-local genomic context (RLGC) for metabolic dark matter.
2. **Stage 2 — Hard Logical & Anatomical Constraint Optimization (`SatConsistencyChecker` + `TaxonConsistencyMaxSat` + `AnatomicalConstraints`):**
   Infers or validates the organism's clade $T_{\mathcal{G}}$, enforces 100% taxonomic and anatomical/developmental consistency via propositional SAT and weighted-MaxSAT minimum-cost repair, and enforces complex/process coherence.
3. **Stage 3 (Optional, Separately Evaluated) — Lightweight Neuro-Symbolic Agentic Anomaly Completion (`gspa-agent`):**
   Acts selectively on the **~20–60 unresolved genome-scale anomalies** surfaced by Stages 1–2 (unfilled essential functions, 1–2 step pathway/module holes with top-$k$ RLGC candidates, incomplete heteromeric complexes, and high-score taxon/anatomical conflicts). Driven by local open-weight Qwen models (**Small Qwen 7B/8B** and **`Qwen 3.8 27B` on `unimatrix01` via `pi-unimatrix`**), the agent inspects sub-threshold domain architectures, localization, syntenic neighborhoods, and local BM25 literature, while every proposed edit is strictly verified by the ELK/SAT4J/Uberon guardrail.

---

## 2. Codebase Analysis: What Exists vs. What Needs to Be Built/Fixed

### 2.1 Strengths Already Implemented in the Repository
1. **Multi-modal Predictor Ecosystem ([`gspa-predictors`](gspa-predictors/)) & Learned Stacker ([`deepgo-plusplus`](deepgo-plusplus/)):**
   - Wrappers for ~16 tools (DIAMOND, MMseqs2, HMMER/Pfam, InterProScan, FoldSeek, mDeepFRI, DeepFRI, ProteInfer, CLEAN, DeepEC, ESM2-DeepGOPlus, gapseq, OperonPredictor, etc.).
   - `DeepGO-PlusPlus` & `DeepGO-PlusPlus-Light`: CAFA6-competitive per-aspect logistic stacker combining homology, InterPro domains, homology-bridged STRING PPI (`net_union`), and 1D-CNN (`f_w ≈ 0.52–0.56` on CAFA6 no-knowledge proteins).
2. **Phase 7 Bayesian Integration Engine ([`gspa-core/src/main/groovy/gspa/integration`](gspa-core/src/main/groovy/gspa/integration)):**
   - Correlation-group collapse + log-odds Noisy-OR ([`EvidenceCombiner.groovy`](gspa-core/src/main/groovy/gspa/integration/EvidenceCombiner.groovy)).
   - Jacobi-damped fixed-point refinement ([`IterativeRefiner.groovy`](gspa-core/src/main/groovy/gspa/integration/IterativeRefiner.groovy)) with full per-annotation provenance ([`ClaimProvenance.groovy`](gspa-core/src/main/groovy/gspa/integration/ClaimProvenance.groovy)).
   - Six modular priors: `EssentialityPrior`, `CoherencePrior`, `ConsistencyPrior`, `GapFillingPrior`, `GenomicContextPrior`, `HomologyTransferPrior`.
3. **GAEF Quality Metrics & Post-hoc Constraint Enforcers ([`gspa-core/src/main/groovy/gspa/metrics`](gspa-core/src/main/groovy/gspa/metrics)):**
   - Full implementation of Tawfiq et al.'s Completeness, Coherence (process `has_part` via ELK, pathway, heteromeric complex), and Consistency (SAT4J with UNSAT core explanation in [`SatConsistencyChecker.groovy`](gspa-core/src/main/groovy/gspa/ontology/SatConsistencyChecker.groovy)).
   - Genome-wide taxonomic inference from predictions ([`TaxonInference.groovy`](gspa-core/src/main/groovy/gspa/metrics/TaxonInference.groovy)) across 106 NCBI clades.
   - Weighted-MaxSAT minimum-cost consistency repair ([`TaxonConsistencyMaxSat.groovy`](gspa-core/src/main/groovy/gspa/ontology/TaxonConsistencyMaxSat.groovy)), `CoherenceEnforcer`, and `CompletenessEnforcer`.
4. **Metabolic Dark-Matter Suggesters ([`gspa-core/src/main/groovy/gspa/integration/suggester`](gspa-core/src/main/groovy/gspa/integration/suggester)):**
   - Operon Bayes-factor suggester ([`DarkMatterSuggester.groovy`](gspa-core/src/main/groovy/gspa/integration/suggester/DarkMatterSuggester.groovy)).
   - Phase 12 Reaction-Local Genomic Context suggester ([`ReactionLocalContextSuggester.groovy`](gspa-core/src/main/groovy/gspa/integration/suggester/ReactionLocalContextSuggester.groovy)) + PLM catalyst centroid ranker (`benchmark/leave_reaction_out/`), which achieves **63.7% hit@1** (MRR 0.660) on *E. coli* Leave-Reaction-Out and **43.1% hit@1** on held-out *B. subtilis* (up from 1.8% with operon-only baseline).

### 2.2 Critical Gaps & Additions to Build for the Paper
1. **Unify `gspa annotate` and `gspa integrate` into the Cohesive 2+1 Stage Pipeline:**
   - Currently, `gspa annotate` runs predictors (e.g. `DeepGO-PlusPlus-Light`) and applies post-hoc heuristic enforcers (`ConsistencyEnforcer`, `CoherenceEnforcer`, `CompletenessEnforcer`), whereas `gspa integrate` runs the Bayesian `IterativeRefiner` on `claims.jsonl` without `TaxonInference`, without `TaxonConsistencyMaxSat`, and without complex coherence in `CoherencePrior`.
   - **Fix:** Unify them so that both entrypoints execute Stage 1 (Bayesian Context Integration) $\to$ Stage 2 (SAT/MaxSAT Taxon & Anatomical Constraint Optimization), with `--enable-agent` invoking Stage 3 (Optional Neuro-Symbolic Agentic Anomaly Completion).
2. **Add Complex Subunit Coherence to `CoherencePrior`:**
   - [`Coherence.groovy`](gspa-core/src/main/groovy/gspa/metrics/Coherence.groovy) evaluates 3 sub-metrics (process `has_part`, pathway, and heteromeric complexes), and [`CoherenceEnforcer.groovy`](gspa-core/src/main/groovy/gspa/metrics/CoherenceEnforcer.groovy) enforces heteromeric complexes, but [`CoherencePrior.groovy`](gspa-core/src/main/groovy/gspa/integration/prior/CoherencePrior.groovy) only includes process and pathway terms—it ignores protein complexes (`complex_terms.tsv` / Complex Portal subunit stoichiometries).
   - **Fix:** Extend `CoherencePrior` to boost missing subunits when companion subunits of an obligate heteromeric complex are present in the genome, and penalize isolated single-subunit claims when no other subunit has support in the genome.
3. **Add Anatomical & Developmental Structure Constraints (Multicellular Eukaryotes):**
   - Implement `AnatomicalConstraints.groovy` and `AnatomicalPrior.groovy` / `AnatomicalEnforcer.groovy` in `gspa-core`:
     - Programmatically extract GO $\leftrightarrow$ Uberon / Plant Ontology (PO) / Fungal Anatomy cross-product axioms from `go-plus.owl` (`results_in_development_of`, `results_in_morphogenesis_of`, `results_in_formation_of`, `occurs_in`, `part_of`) and combine with Uberon/PO taxon constraints (`only_in_taxon`, `never_in_taxon`, `present_in_taxon`) and anatomical hierarchy (`part_of`, `develops_from`).
     - Enforce two rules on eukaryotic genomes (*C. elegans*, *D. melanogaster*, *D. rerio*, *M. musculus*, *A. thaliana*, *P. patens*, *S. cerevisiae*):
       1. **Body-Plan Consistency:** Suppress (or relax to nearest valid ancestor) any GO biological process $f$ whose target anatomical structure $a$ (`f SubClassOf R some a`) does not exist in the organism's clade $T_{\text{org}}$ (e.g., *heart development*, *lung morphogenesis*, or *neural crest migration* in *C. elegans* or *D. melanogaster*; *flower development* in *P. patens* or *S. cerevisiae*).
       2. **Developmental / Structural Coherence:** If an organism possesses structure $a$, specific morphogenetic steps of sub-structures $a' \sqsubseteq_{\text{part\_of}} a$ require core developmental/structural annotations of $a$ to be present in the genome, boosting weak homologs of canonical developmental regulators and demoting isolated false-positive organ calls.
4. **Implement the Optional Stage 3 Neuro-Symbolic Agentic Anomaly Completion Layer (`gspa-agent`):**
   - **Design Principle:** Never run an LLM on the whole proteome. Stages 1–2 resolve 97–99% of proteins deterministically and emit an `anomalies.jsonl` file containing **~20–60 structured genome-scale anomalies per genome**, each paired with $k \le 5$ pre-filtered candidate proteins:
     - *Type 1:* Unfilled essential functions ($0\text{--}5$/genome) + top-$k$ near-ancestor candidates.
     - *Type 2:* Nearly-complete pathways / KEGG modules ($10\text{--}25$/genome, $\ge 75\%$ complete, 1–2 missing steps) + top-$k$ RLGC candidate ORFs.
     - *Type 3:* Incomplete obligate heteromeric complexes ($5\text{--}15$/genome) + candidate partner ORFs or singleton demotion candidates.
     - *Type 4:* High-confidence taxon or anatomical constraint conflicts ($5\text{--}15$/genome) where the agent must choose between *relaxing to a valid conserved ancestor GO term*, *retaining as a bona fide host-interaction effector*, or *removing as a spurious homology transfer*.
   - **Models Evaluated:**
     - **Small Qwen** (e.g., `Qwen2.5-7B-Instruct` / `Qwen3-8B`) for lightweight single-node execution.
     - **`Qwen 3.8 27B` on `unimatrix01`** (invoked via `pi-unimatrix` / local OpenAI-compatible server on `unimatrix01`).
   - **5 Fast Local Inspection & Verification Tools:**
     1. `inspect_candidate_architecture(protein_id)`: Sub-threshold Pfam/InterPro/SUPERFAMILY/CDD domains, active sites, TMbed helices, DeepSig signal peptides, PSORTb localization, CLEAN/DeepEC calls, length/pI.
     2. `inspect_genomic_neighborhood(protein_id, window_kb=10)`: Contig coordinates, strand orientation, intergenic distances (bp), and flanking gene annotations.
     3. `search_local_uniprot_bm25(query, top_k=5)`: Queries the offline SwissProt/PubMed BM25 index (`deepgo-plusplus/pipeline/build_lit_component.py`) for twilight-zone homologs and known **non-orthologous gene displacements (NODs)**.
     4. `inspect_ontology_neighborhood(go_term, organism_taxon)`: Returns GO definition, `has_part` / Uberon anatomical dependencies, parent/child terms, and which ancestor terms are valid in `organism_taxon`.
     5. `commit_verified_action(action_type, protein_id, go_term, rationale)`: Structured output (`PROMOTE_GAP_CANDIDATE`, `RELAX_TO_VALID_ANCESTOR`, `KEEP_AS_HOST_INTERACTION`, `DEMOTE_ORPHAN_SUBUNIT`, `ABSTAIN`).
   - **Hard Neuro-Symbolic Guardrail:** Every proposed action is immediately validated by `SatConsistencyChecker`, `AnatomicalConstraints`, and `GoOntology`. Any invalid GO ID or taxon/anatomical violation is rejected (with 1-step UNSAT core feedback for self-correction), guaranteeing **0% consistency violations** and complete auditability in `<genome>_enforcement_actions.tsv`.
5. **Fix Confounding & IEA Circularity in Existing Benchmarks (`benchmark/RESULTS.md`):**
   - Evaluate with true-path GO ancestor propagation, strict experimental truth (`EXP, IDA, IPI, IMP, IGI, IEP` + CAFA temporal holdout), and unconfounded step-by-step ablations.

---

## 3. Experimental Design & Benchmark Plan

To satisfy both the **CAFA community** (*Nature Methods* / *NAR* reviewers) and **systems biology / genome annotation** users, our evaluation has **five pillars** (with Pillar 5 dedicated specifically to evaluating the optional Stage 3 Agentic Layer):

### Pillar 1: Strict CAFA-Style Protein-Centric & Term-Centric Benchmark (`EXP` Truth + Temporal Holdout)
Prove that **adding genome-scale priors and constraints (Stages 1–2) and optional agentic completion (Stage 3) to SOTA protein predictors improves standard CAFA metrics ($F_{\max}$, IA-weighted $F_{\max}^w$, and $S_{\min}$)** across all three GO aspects (MF, BP, CC), with the largest gains in **BP** and **CC**.

* **Dataset A (21 Well-Annotated Reference Proteomes with High Experimental GOA Coverage):**
  * **Prokaryotes & Archaea (15 genomes):** *E. coli* K-12 MG1655, *B. subtilis* 168, *P. aeruginosa* PAO1, *M. tuberculosis* H37Rv, *S. enterica* LT2, *S. aureus* N315, *V. cholerae* N16961, *C. crescentus* CB15, *S. coelicolor* A3(2), *Synechocystis* PCC 6803, *H. pylori* 26695, *M. genitalium* G37, *M. jannaschii* DSM 2661, *P. furiosus* DSM 3638, *H. volcanii* DS2.
  * **Eukaryotes (6 model genomes spanning Fungi, Plants, Invertebrates, Vertebrates):** *S. cerevisiae*, *S. pombe*, *C. elegans*, *D. melanogaster*, *A. thaliana*, *M. musculus* (or *D. rerio*).
* **Dataset B (CAFA Temporal Holdout):**
  * Proteins from target genomes that gained new experimental GO annotations (`EXP, IDA, IPI, IMP, IGI, IEP`) after cutoff $t_0$ (using the leak-free `deepgo-plusplus` pre-$t_0$ training snapshot).
* **Strict Leave-Clade-Out / Leave-Proteome-Out Controls:**
  * Exclude target species and stratify at $<90\%$, $<50\%$, and $<30\%$ sequence identity.
* **Base Predictors Tested ($\text{Base}$ vs. $+\text{Stage 1}$ vs. $+\text{Stage 2}$ vs. $+\text{Stage 3 Agent}$):**
  1. `DIAMOND + Pfam + InterProScan` (Classical homology/domain stack)
  2. `DeepGO-PlusPlus-Light` (Fast CPU stacker)
  3. `DeepGO-PlusPlus` (Full PLM + Structure + Network stacker)
  4. `ProteInfer`
  5. `DeepFRI / mDeepFRI`
  6. External baseline comparators: **NCBI PGAP**, **Bakta**, **eggNOG-mapper v2**.
* **Metrics Reported (200-iteration bootstrap 95% CIs):**
  * CAFA protein-centric $F_{\max}$ (true-path ancestor propagated)
  * **Information-Accretion weighted $F_{\max}^w$** (official `cafaeval` + `IA.tsv`)
  * **Semantic Distance $S_{\min}$** (Remaining Uncertainty $ru$ vs. Misinformation $mi$)
  * Term-centric AUC-PR

### Pillar 2: GAEF Genome-Scale Plausibility Benchmark (Tawfiq et al. Framework)
Across the 21-genome prokaryote/eukaryote panel AND a **1,000-genome / MAG diversity panel**:
* Measure **Base Predictor vs. +Stage 1 vs. +Stage 2 vs. +Stage 3 (Agent)**:
  1. **Completeness:** % of kingdom-essential GO terms covered.
  2. **Process Coherence:** % of triggered GO `has_part` dependencies satisfied.
  3. **Pathway Coherence:** Mean coverage of triggered KEGG Pathways, KEGG Modules, and MetaCyc pathways.
  4. **Complex Coherence:** % of obligate heteromeric complexes (Complex Portal / GO CC) with $\ge 2$ distinct subunits annotated.
  5. **Taxonomic & Anatomical Consistency:** % of proteins and genomes with zero SAT4J taxon-constraint or Uberon/PO anatomical-constraint violations.

### Pillar 3: Five Concrete Biological Case Studies ("Where Genome Knowledge Solves What Sequence Alone Cannot")
1. **Vignette 1 — Metabolic Model Gap-Filling & "Dark Matter" Enzyme Recovery (`gapseq` + RLGC + PLM + Optional Agent):**
   - **Experiment:** Leave-Reaction-Out (LRO) benchmark across 10 diverse prokaryotic genomes (expanding `benchmark/leave_reaction_out/`). Mask the true enzyme's sequence homology claims, create a metabolic gap, and evaluate gene recovery (`hit@1`, `hit@3`, `hit@5`, MRR) across:
     - Operon/pathway heuristic baseline (1.8% `hit@1`)
     - Stage 1 RLGC context-only (22.9% `hit@1` on *E. coli*)
     - Stage 1 RLGC + PLM catalyst centroid (63.7% `hit@1` on *E. coli*, 43.1% on held-out *B. subtilis*)
     - Stage 3 Agentic re-ranking of top-5 RLGC candidates using sub-threshold domain folds, localization, and BM25 literature.
   - Also evaluate downstream impact on **FBA metabolic model viability / carbon source phenotype prediction**.
2. **Vignette 2 — Pathway Completion via Operon & Module Coherence:**
   - **Experiment:** Leave-One-Enzyme-Out on curated KEGG Modules / biosynthetic operons (`trp`, `his`, cobalamin, peptidoglycan) where one enzyme's sequence score is attenuated to borderline ($0.20 \le p \le 0.45$).
3. **Vignette 3 — Obligate Heteromeric Protein Complex Completion & Orphan Subunit Pruning:**
   - **Experiment:** Evaluate on Complex Portal gold-standard complexes across *E. coli*, *S. cerevisiae*, and *M. musculus* (ATP synthase, Complex I, 20S proteasome, SecYEG, RNA polymerase), measuring both **orphan subunit false-positive removal** and **Leave-Subunit-Out recovery**.
4. **Vignette 4 — Species/Clade-Based Taxon Inference & Decontamination (Single Genomes & Metagenomic Contigs):**
   - **Experiment:** Run `TaxonInference` + `ConsistencyPrior` / `TaxonConsistencyMaxSat` across 21 genomes and multi-contig metagenomes (`mock_community`), quantifying how removing the ~15–22% taxonomically impossible calls reduces CAFA Misinformation ($mi$) and boosts $F_{\max}^w$.
5. **Vignette 5 — Anatomical & Developmental Structure Constraints in Multicellular Eukaryotes:**
   - **Experiment:** Evaluate GO Biological Process predictions across 6 eukaryotes with contrasting body plans (*S. cerevisiae*, *C. elegans*, *D. melanogaster*, *A. thaliana*, *P. patens* / *C. reinhardtii*, *M. musculus* / *D. rerio*).
   - Quantify how GO–Uberon–PO logical constraints eliminate impossible organ/tissue development calls on conserved signaling proteins (e.g., *heart development*, *lung morphogenesis*, or *neural tube closure* in *C. elegans* or *Drosophila*; *flower development* in moss/yeast) while preserving/relaxing to the true conserved developmental/signaling ancestor.

### Pillar 4: Dedicated Evaluation of the Optional Stage 3 Agentic Layer (Small Qwen vs. Qwen 3.8 27B on `unimatrix01`)
Because Stage 3 is an **optional layer**, we evaluate it in a standalone section and dedicated figure/table:
1. **Accuracy & Plausibility Lift ($\Delta$ Stage 3 over Stages 1+2):**
   - Measure marginal improvement in:
     - **LRO Metabolic Gap-Filling (`hit@1`, `hit@3`, MRR):** Re-ranking top-5 RLGC candidates when PLM/homology is ambiguous or when the gap is filled by a **non-orthologous gene displacement (NOD)**.
     - **Leave-Subunit-Out Complex Recovery:** Identifying missing complex subunits via sub-threshold domain/localization/synteny inspection.
     - **Semantic Ancestor Relaxation vs. Hard Deletion on Taxon/Anatomical Conflicts:** Measure CAFA $F_{\max}^w$ and $S_{\min}$ when the agent relaxes forbidden organ/taxon terms (e.g. `heart development` $\to$ `regulation of cell fate specification` in *C. elegans*) vs. Stage 2's hard deletion.
2. **Model Scaling Comparison (`Small Qwen 7B/8B` vs. `Qwen 3.8 27B` via `pi-unimatrix` on `unimatrix01`):**
   - Compare decision accuracy, tool-use precision, abstention calibration, and wall-clock latency/throughput (seconds per genome, anomalies resolved per minute) between the small Qwen model and `Qwen 3.8 27B`.
3. **Ablation of the Neuro-Symbolic Guardrail (`Agent without SAT/ELK Guardrail` vs. `Agent with SAT/ELK Guardrail`):**
   - Quantify how often an unguarded LLM agent proposes invalid GO IDs or taxon/anatomically inconsistent terms, proving that the **neuro-symbolic SAT/OWL guardrail is essential** to guarantee 0% constraint violations.

---

## 4. Concrete Implementation Roadmap

### Phase 1: Core Unification & New Priors/Constraints (`gspa-core` & `gspa-cli`)
1. **Fix & Verify Gradle Test Suite:**
   - Ensure `./gradlew clean test` passes 100% cleanly across `gspa-core`, `gspa-predictors`, and `gspa-cli`.
2. **Add Anatomical & Developmental Constraint Engine (`gspa-core`):**
   - Programmatically extract `go-anatomical-constraints.tsv`, `uberon-po-taxon-constraints.tsv`, and `uberon-po-hierarchy.tsv` from official `go-plus.owl`, `uberon.owl`, and `po.owl` (strictly validating all ontology IDs programmatically).
   - Implement `AnatomicalConstraints.groovy`, `AnatomicalPrior.groovy`, and ancestor-relaxation support in `ConsistencyEnforcer.groovy`.
3. **Add Complex Subunit Coherence to `CoherencePrior`:**
   - Extend [`CoherencePrior.groovy`](gspa-core/src/main/groovy/gspa/integration/prior/CoherencePrior.groovy) to incorporate Complex Portal / heteromeric complex subunit co-occurrence boosts and singleton penalties.
4. **Unify `gspa integrate` and `gspa annotate`:**
   - Wire `TaxonInference` and `ConsistencyPrior` + `TaxonConsistencyMaxSat` directly into the unified pipeline so Stages 1 and 2 run seamlessly.

### Phase 2: Implement Optional Stage 3 Agentic Layer (`gspa-agent`)
1. **Anomaly Exporter (`AnomalyExtractor.groovy`):**
   - After Stages 1–2, emit `<genome>_anomalies.jsonl` containing the ~20–60 unresolved items (missing essentials, 1–2 hole pathways with top-5 RLGC candidates, incomplete complexes, and taxon/anatomical conflict terms) with their pre-extracted candidate context.
2. **Agent Sidecar (`benchmark/agent/run_gspa_agent.py` + `AgentCompletionRunner.groovy`):**
   - Implement the 5 deterministic local tools (`inspect_candidate_architecture`, `inspect_genomic_neighborhood`, `search_local_uniprot_bm25`, `inspect_ontology_neighborhood`, `commit_verified_action`).
   - Support both **`pi-unimatrix` / `Qwen 3.8 27B` on `unimatrix01`** and **local Small Qwen (`Qwen2.5-7B` / `Qwen3-8B`)** via OpenAI-compatible structured tool calling.
   - Wire the **SAT4J + ELK + Anatomical Guardrail** loop so every committed action is verified before updating `<genome>_annotations.tsv` and `<genome>_enforcement_actions.tsv`.
   - Expose via `--enable-agent` (`--agent-model qwen-small | qwen3.8-27b-unimatrix`) in `gspa annotate` and `gspa integrate` (default: off).

### Phase 3: Benchmark Execution
1. **Run Pillar 1 & Pillar 2 (CAFA + GAEF Ablations across 15 Prokaryotes + 6 Eukaryotes):**
   - Compare `Base` $\to$ `+ Stage 1 (Bayesian Priors)` $\to$ `+ Stage 2 (SAT/MaxSAT + Anatomy)` $\to$ `+ Stage 3 (Optional Agent)` across 5 base predictors and against PGAP / Bakta / eggNOG-mapper v2.
2. **Run Pillar 3 (5 Biological Vignettes):**
   - 10-genome Leave-Reaction-Out (LRO) metabolic gap-filling + FBA phenotype validation.
   - Leave-One-Enzyme-Out pathway recovery & Leave-Subunit-Out Complex Portal recovery.
   - 6-eukaryote Anatomical & Developmental constraint benchmark.
3. **Run Pillar 4 (Stage 3 Agentic Layer Evaluation):**
   - Compare `Small Qwen` vs. `Qwen 3.8 27B (unimatrix01)` on LRO re-ranking, complex completion, and ancestor relaxation vs. hard deletion, plus guarded vs. unguarded agent ablation and runtime profiling.

### Phase 4: Manuscript Figures & Writing
1. **Figure 1:** Three-Stage GSPA Architecture (Stage 1 Bayesian Context Integration $\to$ Stage 2 SAT/MaxSAT Taxon & Anatomical Optimization $\to$ Optional Stage 3 Neuro-Symbolic Agentic Anomaly Resolver).
2. **Figure 2:** Main CAFA ($F_{\max}^w$, $S_{\min}$) & GAEF (Completeness, Coherence, Consistency) benchmark across 15 prokaryotes and 6 eukaryotes + stepwise ablation waterfall.
3. **Figure 3:** Metabolic Dark-Matter Recovery (10-genome Leave-Reaction-Out benchmark + FBA phenotype accuracy + operon locus vignette).
4. **Figure 4:** Pathway & Obligate Heteromeric Complex Completion + Cross-Clade Taxon & Eukaryotic Anatomical/Developmental Decontamination Heatmap (*S. cerevisiae*, *C. elegans*, *Drosophila*, *Arabidopsis*, moss, mouse).
5. **Figure 5:** Dedicated Evaluation of the Optional Stage 3 Agentic Layer (`Small Qwen` vs. `Qwen 3.8 27B` on `unimatrix01`: LRO/complex lift, ancestor relaxation vs. deletion, guarded vs. unguarded hallucination rate, and runtime/cost scaling).
