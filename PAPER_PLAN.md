# GSPA Publication & Research Plan (Upgraded for *Nature Methods* / *Nature Machine Intelligence*)
**Target Venues:** *Nature Methods* / *Nature Machine Intelligence*  
**Working Title:** *Genome-scale probabilistic inference, stoichiometric flux-guided locus attribution, and neuro-symbolic constraint optimization solve protein function annotation and microbial phenotype prediction*

---

## 1. Central Research Question & Thesis

### Research Question
> **Can explicit genome-scale biological context, stoichiometric flux-balance constraints with locus-specific credit assignment, and formal ontological reasoning—spanning metabolic networks, multi-subunit protein complexes, operon/chromosomal organization, essential cellular repertoires, and taxonomic as well as anatomical/developmental boundaries—measurably improve both protein-level function prediction accuracy (CAFA metrics) and organism-level phenotype prediction (reproducing experimentally measured microbial phenotypes in BacDive) over protein-in-isolation predictors and existing genome-scale reconstruction pipelines?**

### Central Thesis
All standard protein function predictors (sequence homology, domain scanners, protein language models, structure GNNs) and standard genome annotation pipelines (NCBI PGAP, Bakta, Prokka, eggNOG-mapper, InterProScan) operate under a **reductionist independence assumption**:
$$P(\mathbf{Y}_{1:N} \mid \mathbf{X}_{1:N}) = \prod_{i=1}^N P(Y_i \mid X_i)$$
where each protein $i \in \{1,\dots,N\}$ in a genome is annotated in isolation. Conversely, existing metabolic network post-processors (such as METEOR) and naive gap-fillers operate at the **reaction or EC level without protein-level identifiability**, suffering from proteome-size Noisy-OR inflation, zero-flux free-rider selection, inability to distinguish isozyme `OR` from complex subunit `AND`, and paralog-blind feedback that shifts all proteins uniformly without attributing a required metabolic reaction to a specific genomic locus.

**GSPA solves both halves of this problem** through a **three-stage genome-scale inference architecture** coupled to a dedicated **GPR-aware, flux-certified FBA gap-filling and locus-attribution method (`StoichiometricLocusRefiner`)**:
1. **Stage 1 — Soft Bayesian Context Integration & Sparse GPR Hypothesis Generation (`IterativeRefiner` + `ReactionLocalContextSuggester` + `PLMCatalystCentroidRanker`):**
   $$P(\mathbf{Y}_{1:N} \mid \mathbf{X}_{1:N}, \mathcal{G}, \mathcal{O}) \propto \prod_{i=1}^N P_{\text{base}}(Y_i \mid X_i) \cdot \Psi_{\text{taxon}}(\mathbf{Y}, T_{\mathcal{G}}) \cdot \Psi_{\text{anat}}(\mathbf{Y}, \mathcal{A}_{\mathcal{G}}) \cdot \Psi_{\text{pwy}}(\mathbf{Y}, \mathcal{M}_{\mathcal{G}}) \cdot \Psi_{\text{cplx}}(\mathbf{Y}, \mathcal{C}) \cdot \Psi_{\text{ctx}}(\mathbf{Y}, \mathcal{L}_{\mathcal{G}}) \cdot \Psi_{\text{ess}}(\mathbf{Y}, \mathcal{R}_{T_{\mathcal{G}}})$$
   Calibrates multi-modal protein evidence, prevents proteome-size saturation via sparse top-$k$ / Poisson-truncated role priors, models complete Gene–Protein–Reaction (GPR) mechanisms (`OR` across isozymes, `AND` across obligate complex subunits), and generates candidate protein–role hypotheses $h \in \mathcal{H}_j$ (including an explicit unresolved/null hypothesis $h_\varnothing$) using genomic neighborhood density (RLGC), PLM catalyst centroids, and domain/localization compatibility.
2. **Stage 2 — Joint Logical, Anatomical & Stoichiometric Constraint Optimization (`SatConsistencyChecker` + `TaxonConsistencyMaxSat` + `AnatomicalConstraints` + `StoichiometricLocusRefiner`):**
   - **Logical & Anatomical Consistency:** Enforces 100% taxonomic and eukaryotic body-plan/developmental consistency (GO–Uberon–PO cross-products) via propositional SAT and weighted-MaxSAT minimum-cost repair / ancestor relaxation.
   - **Stoichiometric Flux-Certified Gap-Filling & Locus Attribution:** Solves a positive-cost, flux-certified repair MILP followed by lexicographic parsimonious FBA (pFBA) and Flux Variability / Knockout Necessity certificates, then maps flux-required reaction repairs back to specific genomic loci via Bayesian posterior responsibility $r_j(h)$.
3. **Stage 3 (Optional, Separately Evaluated) — Lightweight Neuro-Symbolic Agentic Anomaly Completion (`gspa-agent`):**
   Acts selectively on the **~20–60 high-uncertainty genome-scale anomalies** surfaced by Stages 1–2 (unresolved FBA repair candidates where $r_j(h_\varnothing)$ or paralog entropy is high, 1–2 step pathway holes, incomplete obligate complexes, and borderline taxon/anatomical conflicts). Driven by local open-weight Qwen models (`Qwen 3.8 27B` on `unimatrix01` / small Qwen), every agent proposal is verified by BOTH the SAT/OWL ontology guardrail and Stage 2 stoichiometric feasibility before acceptance.

---

## 2. Diagnostic of METEOR & Existing Gap-Fillers vs. GSPA's Clean FBA Integration

### 2.1 Five Mathematical & Biological Flaws in Naive FBA/MILP Refinement (METEOR Audit)
Our audit of `/home/leechuck/Documents/papers/meteor/` (`draft.tex`, `supplementary.tex`, `MESSAGE_TO_KEXIN.md`, `REVIEW_METEOR.txt`) and GSPA's existing `GapFillingPrior.groovy`, synthesized with independent analysis from `codex-auto` (`gpt-6-astra`), identified five fundamental flaws in naive MILP enzyme refinement:

1. **The Zero-Flux Free-Rider Problem (`K07`):**
   METEOR minimizes $\sum_j c_j y_j + \lambda \sum_j |v_j|$ with $c_j = -\log\frac{w_j+\epsilon}{1-w_j+\epsilon}$ and $L_j y_j \le v_j \le U_j y_j$. Because $y_j = 1$ permits $v_j = 0$, every reaction with $w_j > 0.5$ has $c_j < 0$ and is trivially set to $y_j^\star = 1$ with $v_j^\star = 0$, even if it is completely disconnected from the metabolic network, stoichiometrically dead-ended, or mass-imbalanced. "Selected" ($y_j=1$) is conflated with "flux-carrying" ($|v_j| > 0$) and "biomass-essential."
2. **Proteome-Size Noisy-OR Saturation (`K09`):**
   Aggregating $w_j = 1 - \prod_{p \in \mathcal{P}, e \in M(j)} (1 - E_{p,e})$ across all $N \approx 4{,}000$ proteins in a bacterial proteome means a tiny background score $E_{p,e} = 0.0002$ across 4,000 proteins yields $w_j = 1 - (0.9998)^{4000} = 0.551 > 0.5 \implies c_j < 0$. This causes massive false-positive reaction activation (`10,003` reactions per genome in DeepProZyme, of which `4,086` are mass-imbalanced).
3. **Paralog Blindness & Uniform Feedback (`K08`, `K15`):**
   METEOR's feedback update $E'_{p,e} = E_{p,e} \pm \alpha E_{p,e}(1 - E_{p,e})$ and GSPA's existing `GapFillingPrior` (`logOddsBoost` independent of `proteinId`) apply the exact same directional shift to *every* protein in the proteome for EC $e$:
   - It preserves within-EC protein ranking, so it can **never disambiguate competing paralogs**.
   - Because $E_{p,e}(1 - E_{p,e}) = 0$ at $E_{p,e} = 0$, it **cannot lift zero scores or recover "dark matter" genes** missing from the baseline predictor's output vocabulary.
   - If a reaction maps to multiple alternative ECs $M(j) = \{e_1, e_2\}$, activating $y_j = 1$ blindly boosts all ECs in $M(j)$ across all proteins, even when only one catalytic mechanism is used.
4. **Conflation of Isozyme `OR` with Complex Subunit `AND`:**
   Mapping reactions to flat EC sets ignores Gene–Protein–Reaction (GPR) logic: alternative isozymes act disjunctively (`OR`), whereas obligate heteromeric enzyme complexes (e.g., ATP synthase, NADH dehydrogenase, pyruvate dehydrogenase, urease $\alpha\beta\gamma$) require a conjunction (`AND`) of distinct catalytic and structural subunits.
5. **Conflation of Reconstruction Medium with Diagnostic Phenotype Probing (`K15`, `K21`):**
   Enforcing growth on a single static medium misses catabolic, respiratory, and secondary pathways uncoupled from that medium, while opening all exchanges post hoc (as in METEOR's downstream builder) invalidates standalone FBA testing. Reconstructing a draft model under permitted viability assumptions must be strictly separated from **diagnostic flux probing** across external test conditions (such as BacDive carbon, nitrogen, energy, or fermentation assays).

---

### 2.2 Clean Mathematical Formulation: `StoichiometricLocusRefiner` (Synthesizing Our Model + Astra's Design)

We replace protein-blind reaction boosting with a four-part probabilistic & stoichiometric framework that separates **(i) reaction presence**, **(ii) condition-specific flux support**, **(iii) task necessity**, and **(iv) locus-specific protein responsibility**:

#### Part 1: Sparse, GPR-Aware Mechanism & Reaction Evidence (Solving Flaws 2 & 4)
Let $p \in \mathcal{P}$ index proteins in genome $\mathcal{G}$, $r \in \mathcal{R}$ index specific catalytic or structural subunit roles, $j \in \mathcal{J}$ index mass- and charge-balanced biochemical reactions (preserving substrate, cofactor, compartment, and transport specificity rather than collapsing distinct reactions by EC), and $k \in \mathcal{K}_j$ index alternative enzyme mechanisms catalyzing reaction $j$.

1. **Bounded Candidate Retrieval + Sparse Role Prior:**
   For each role $r$, rather than taking an all-proteome Noisy-OR over $N \approx 4{,}000$ proteins, we construct a bounded candidate set $C_r \subseteq \mathcal{P}$ ($|C_r| \le K_{\max} = 15$) from the union of:
   - Sequence/PLM/domain hits above a calibrated noise floor, and
   - Top-ranked dark-matter candidates from `ReactionLocalContextSuggester` (operon/chromosomal density + reaction-graph neighbor anchors) and `PLMCatalystCentroidRanker` (ESM2 catalytic centroid similarity + InterPro/Pfam domain compatibility).
   We place a truncated Poisson sparsity prior $K_r \sim \text{TruncPoisson}(\lambda_r, K_{\max})$ on the number of true loci encoding role $r$ in a genome (learned from reference genomes, typically $\lambda_r \in [1.0, 1.5]$), eliminating proteome-size score inflation by construction.
2. **Explicit Boolean GPR Mechanism Scoring:**
   Let $x_{p,r} \in \{0,1\}$ indicate that protein $p$ performs role $r$. A catalytic mechanism $k \in \mathcal{K}_j$ requires all of its obligate subunit roles $R_{j,k}$ (`AND`), each of which can be supplied by any isozyme candidate $p \in C_r$ (`OR`), while reaction $j$ is available if any complete mechanism $k \in \mathcal{K}_j$ is present (`OR`):
   $$G_{j,k}(x) = \bigwedge_{r \in R_{j,k}} \left( \bigvee_{p \in C_r} x_{p,r} \right), \qquad y_j(x) = \bigvee_{k \in \mathcal{K}_j} G_{j,k}(x)$$
   In probability space, with calibrated protein–role posteriors $p_{p,r} = P(x_{p,r}=1 \mid X)$ from GSPA Stage 1:
   $$P(r \text{ supplied} \mid X) = 1 - \prod_{p \in C_r} (1 - p_{p,r}), \qquad P(G_{j,k}=1 \mid X) = \prod_{r \in R_{j,k}} P(r \text{ supplied} \mid X)$$
   $$w_j = P(y_j = 1 \mid X) = 1 - \prod_{k \in \mathcal{K}_j} \big(1 - P(G_{j,k}=1 \mid X)\big)$$
   *(accounting for shared subunits across mechanisms via exact small-factor marginalization when $|\bigcup_k R_{j,k}| \le 10$).*

#### Part 2: Positive-Cost, Flux-Certified Repair MILP + Lexicographic pFBA (Solving Flaws 1 & 5)
We Partition reactions into the high-confidence draft core $B = \{j \in \mathcal{J} : w_j \ge \tau_{\text{core}}\}$ and the repair candidate pool $\mathcal{J}_{\text{cand}} = \mathcal{J} \setminus B$.
1. **Strictly Positive Repair Costs (No Zero-Flux Free Riders):**
   Only additions $a_j \in \{0,1\}$ for $j \in \mathcal{J}_{\text{cand}}$ incur a strictly positive repair cost:
   $$c_j = c_{\min} - \log(w_j + \epsilon) + d_j > 0$$
   where $c_{\min} > 0$ is a positive floor and $d_j \ge 0$ penalizes reactions lacking a high-confidence genomic candidate ($h_\varnothing$, database-only orphan repairs) or reactions with thermodynamic/directionality uncertainty (excluding mass-imbalanced reactions and energy-generating cycles from $\mathcal{J}_{\text{cand}}$).
2. **Explicit Flux-Support Certificates ($b^+_{jc}, b^-_{jc}$):**
   For each permitted reconstruction condition $c \in \mathcal{C}_{\text{rec}}$ (e.g., minimal/defined viable growth medium for prototrophs or clade-specific core Auxotrophy-rescued medium, **never** leaking held-out test phenotypes), we decompose flux $v_{j,c} = v^+_{j,c} - v^-_{j,c}$ with direction binary indicators $b^+_{j,c}, b^-_{j,c} \in \{0,1\}$:
   $$\mathbf{S}\mathbf{v}^{(c)} = \mathbf{0}, \qquad \mathbf{t}_c^\top \mathbf{v}^{(c)} \ge \tau_c$$
   $$\epsilon_{\text{flux}} b^+_{j,c} \le v^+_{j,c} \le U^+_{j,c} b^+_{j,c}, \qquad \epsilon_{\text{flux}} b^-_{j,c} \le v^-_{j,c} \le U^-_{j,c} b^-_{j,c}, \qquad b^+_{j,c} + b^-_{j,c} \le y_j$$
   where $y_j = 1$ for $j \in B$ and $y_j = a_j$ for $j \in \mathcal{J}_{\text{cand}}$, and $\epsilon_{\text{flux}} = 10^{-5} \gg \text{solver tolerance}$. We enforce:
   $$a_j \le \sum_{c \in \mathcal{C}_{\text{rec}}} (b^+_{j,c} + b^-_{j,c}) \qquad \forall j \in \mathcal{J}_{\text{cand}}$$
   **Mathematical Guarantee:** Every repaired reaction $a_j^\star = 1$ is mathematically guaranteed to carry non-zero flux $|v_{j,c}^\star| \ge \epsilon_{\text{flux}}$ in at least one witness condition $c$ and cannot be a disconnected free rider.
3. **Two-Stage Lexicographic Optimization + Alternative Optima Enumeration:**
   - **Step 2a (Minimum-Cost Mechanism Repair):** Solve $\min_{\mathbf{a}, \mathbf{b}, \mathbf{v}} \sum_{j \in \mathcal{J}_{\text{cand}}} c_j a_j$ to optimality (recording solver status, MIP gap, and runtime).
   - **Step 2b (Lexicographic pFBA Witness):** Fix $\sum_{j \in \mathcal{J}_{\text{cand}}} c_j a_j \le (1 + \delta_{\text{opt}}) Z^\star$ and minimize $\sum_{c, j} |v_{j,c}|$ to obtain a parsimonious flux witness.
   - **Step 2c (Necessity & Ensemble Stability Certificates):** Run Flux Variability Analysis (FVA), single-reaction knockout ($\max \mathbf{t}_c^\top \mathbf{v}$ with $y_j = 0$), and $M \le 5$ integer-cut alternative repair solutions to classify every added reaction into four distinct, auditable tiers:
     1. `WITNESS_FLUX_ACTIVE` (carries flux $\ge \epsilon_{\text{flux}}$ in optimal pFBA witness),
     2. `FVA_COUPLED` (non-zero flux range under FVA),
     3. `MODEL_ESSENTIAL` (knockout in reconstructed model drops task below $\tau_c$),
     4. `ENSEMBLE_ESSENTIAL` (present across all near-optimal repair solutions within $\delta_{\text{opt}}$).

#### Part 3: Backward Locus Credit Assignment via Bayesian Responsibility (Solving Flaw 3)
When reaction $j \in \mathcal{J}_{\text{cand}}$ is repaired ($a_j^\star = 1$), how do we update protein-level annotations without METEOR's uniform proteome shift?
For each repaired reaction $j$, we define the hypothesis space $\mathcal{H}_j = \{h_1, \dots, h_M, h_\varnothing\}$, where each $h_m = (k, \{p_r\}_{r \in R_{j,k}})$ is a specific mechanism $k$ assigned to specific candidate protein loci $p_r \in C_r$, and $h_\varnothing$ is the **explicit null/unresolved hypothesis** ("reaction $j$ is required stoichiometrically, but the encoding gene is outside $C_r$, non-orthologous, or unsequenced").
Each hypothesis $h_m$ has a local prior log-odds score combining:
$$\ell(h_m \mid X) = \sum_{r \in R_{j,k}} \Big( \underbrace{\beta_{\text{seq}} s_{\text{seq}}(p_r, r)}_{\text{Stage 1 calibrated log-odds}} + \underbrace{\beta_{\text{rlgc}} s_{\text{rlgc}}(p_r, j)}_{\text{operon/synteny + neighbor density}} + \underbrace{\beta_{\text{plm}} s_{\text{plm}}(p_r, r)}_{\text{ESM2 catalyst centroid cosine}} + \underbrace{\beta_{\text{dom}} s_{\text{dom}}(p_r, r)}_{\text{sub-threshold Pfam/active-site}} - \underbrace{\gamma_{\text{commit}} c_{\text{commit}}(p_r)}_{\text{commitment to orthogonal roles}} \Big)$$
while the unresolved hypothesis $h_\varnothing$ has a calibrated threshold score $\ell(h_\varnothing) = \tau_\varnothing$.
The **posterior locus responsibility** of hypothesis $h \in \mathcal{H}_j$ given the stoichiometric repair evidence $D_{\text{flux}}$ is:
$$r_j(h) = \frac{\exp(\ell(h \mid X)) \cdot P(D_{\text{flux}} \mid h)}{\sum_{h' \in \mathcal{H}_j} \exp(\ell(h' \mid X)) \cdot P(D_{\text{flux}} \mid h')}$$
- **Locus Identifiability Property:** Only candidate proteins $p \in C_r$ receive a posterior update proportional to their marginal responsibility $r_{j,p} = \sum_{h \ni p} r_j(h) \cdot \omega_{\text{ens}}(j)$ (weighted by reaction $j$'s ensemble stability $\omega_{\text{ens}}(j) \in [0,1]$).
- **Dark-Matter Creation:** Because $C_r$ includes candidates retrieved by RLGC genomic context and PLM catalyst centroids even when baseline sequence predictors emitted zero ($s_{\text{seq}} = \text{missing}$), GSPA can **promote true dark-matter genes from zero**, whereas if all candidates are weak, $r_j(h_\varnothing)$ dominates and the reaction is recorded as a *model-level stoichiometric gap repair without a false-positive gene call*.
- **No Double-Counting:** Inferred responsibilities remain revisable across outer iterations rather than pinned as hard floors, and Stage 2 SAT/MaxSAT ontological consistency is re-verified jointly so no metabolic gene promotion ever violates GO taxonomic or anatomical constraints.

---

## 3. Rigorous Benchmark Architecture, Metabolic Role Stratification & Multi-Channel Homology Masking

A central methodological contribution of this work (developed jointly with Astra's independent audit in [`astra_benchmark_design_proposal.md`](file:///home/leechuck/.gemini/antigravity-cli/brain/037ff17e-e13e-4016-9d69-584e0c2b3ec4/scratch/astra_benchmark_design_proposal.md)) is **delineating what each benchmark actually measures, what is biologically possible to recover under each constraint regime, and how to prevent indirect homology leakage across model organisms**.

### 3.1 Diagnostic of Overfitting & Structural Limits in Existing Benchmarks (`Price-149`, `CAFA`, and `LRO`)

1. **What `Price-149` Actually Tests (and Why Static Minimal-Biomass FBA Cannot Solve It Alone):**
   - **Biological Composition:** Price et al. (2018, *Nature*) discovered functions for bacterial genes of unknown function using RB-TnSeq mutant fitness profiling across **hundreds of distinct sole carbon, nitrogen, and stress conditions** (reannotating 238 transport proteins and 218 catabolic proteins). The `Price-149` subset (curated by ProteInfer and popularized by CLEAN, spanning 149 proteins and 56 EC classes across 22 bacterial genomes) is enriched for **alternative-substrate catabolism, peripheral assimilation, and transport-coupled pathways**—which are *primary metabolism under a specific nutrient condition*, **not constitutive minimal-glucose biomass synthesis** and not secondary natural-product biosynthesis.
   - **Why Single-Medium Biomass MILP Fails on Catabolic Enzymes:** In a minimal-glucose model, an enzyme initiating degradation of L-fucose, D-galactonate, myo-inositol, or D-glucosamine carries **zero flux ($v_j = 0$)** and has **zero knockout effect on glucose biomass**. A single-condition biomass MILP supplies no affirmative stoichiometric evidence for that enzyme; in METEOR (`price149_6baseline.tsv`), apparent gains came from negative-cost zero-flux free-riders ($w_j > 0.5$), whereas true enzymes with baseline score $< 0.5$ were actively suppressed (causing top-$k$ regressions in EnzBERT and GraphEC).
   - **What *Can* Legitimately Recover `Price-149` Enzymes:**
     1. **Genomic Context & Pathway/Module Coherence (Stage 1 RLGC + CoherencePrior):** In bacteria, inducible catabolic enzymes almost invariably reside in **substrate-specific operons** alongside their cognate ABC/MFS transporters, transcriptional regulators, and adjacent catabolic steps. Even without knowing the growth medium, chromosomal neighborhood density (`RLGC`) + pathway neighbor coherence provides strong locus-specific evidence.
     2. **Condition-Specific FBA Task Probing (Stage 2 `StoichiometricLocusRefiner`):** When the RB-TnSeq growth substrate (the carbon or nitrogen source on which the mutant had a fitness defect) is supplied as the condition-specific FBA task $\Phi_c(y)$, the missing catabolic reaction becomes **flux-essential**, allowing `StoichiometricLocusRefiner` to repair the reaction and attribute it to the operon candidate!
   - **Post-2018 Database Contamination in `Price-149`:** Because Price et al. was published in May 2018, many of these 149 experimental discoveries were subsequently absorbed into UniProt/SwissProt, RefSeq, KEGG, InterPro (`TIGRFAM`/`NCBIfam`/`PANTHER`/`HAMAP`), and post-2018 PLM training sets. Evaluating modern checkpoints on `Price-149` without multi-channel homology masking (`H`) or pre-2018 temporal snapshots (`T`) measures retrospective database memorization rather than discovery.

2. **Why Naive "Leave-Species-Out" on Model Organism Proteomes Leaks Massively:**
   - Dense experimental GO/EC annotations (`EXP, IDA, IPI, IMP, IGI, IEP`) are concentrated in ~15 canonical model organisms (*E. coli* K-12, *B. subtilis* 168, *S. cerevisiae*, *S. pombe*, *C. elegans*, *D. melanogaster*, *A. thaliana*, *M. musculus*, *P. aeruginosa*, *M. tuberculosis*).
   - Over decades of biocuration, experimental annotations discovered in *E. coli* K-12 or *S. cerevisiae* have **propagated outward** into:
     1. **Sister species & genera in SwissProt** (e.g., *Shigella*, *E. coli* O157:H7, *Salmonella enterica* sharing 85–99% identity, where transferred annotations were elevated to reviewed SwissProt entries),
     2. **Domain-to-GO curated mapping rules (`InterPro2GO`, `HAMAP`, `TIGRFAM`, `PANTHER`, `FunFam`)** whose seed alignments and functional rules were defined directly on the model organism's protein,
     3. **STRING PPI networks (`net`, `net_bridge`, `net_union`)** where transferred interolog edges or self-matched bridge homologs return the model organism's own labels,
     4. **SwissProt protein names and descriptions (`lit` BM25)** which explicitly state the curated enzyme name,
     5. **Supervised PLM heads and EC catalyst centroids (`ESM2-MLP`, `ProstT5`, `ProteInfer`, `CLEAN`, `PLMCatalystCentroidRanker`)** trained on full SwissProt/RefSeq panels containing the test organism or its >70% identity orthologs (as demonstrated by our `leave_reaction_out/BASELINE_REPORT.md` audit, where `mg1655` in-panel scored `0.637` hit@1 vs. held-out `bsubtilis` `0.431`).
   - **The GAF Last-Modified Date Bug in CAFA Reconstructions:** As proven in our `deepgo-plusplus/TRAINING.md` audit (§1.0a), filtering GAF files by column 14 (`date`) mistakes re-curated pre-$t_0$ annotations for novel post-$t_0$ discoveries (contaminating 95% of `no-knowledge` MF targets and inflating `STRING net` MF from `0.441` to `0.803`). True temporal evaluation requires **strict archive snapshot diffing**: $\Delta E(p,a) = E_{t_1}(p,a) \setminus E_{t_0}(p,a)$ after resolving accession aliases, isoforms, and ontology term merges.

---

### 3.2 Two Orthogonal Generalization Axes & The Multi-Channel Homology Masking Protocol (`H` & `T`)

We explicitly separate two scientific questions that previous papers conflate:
* **Axis 1 — Temporal Discovery (`T`, Prospective CAFA Protocol):** *Given all biological knowledge publicly available at snapshot $t_0$ (including close homologs in other species), can we predict new experimental annotations $\Delta E(p,a) = E_{t_1}(p,a) \setminus E_{t_0}(p,a)$ acquired between $t_0$ and $t_1$?*
* **Axis 2 — Sequence & Clade Generalization (`H`, Multi-Channel Homology Masking Protocol):** *Given a genome whose close sequence relatives are withheld across EVERY evidence channel simultaneously, can genome-scale context, stoichiometry, and ontological constraints recover functions as sequence similarity drops from $>70\%$ to $50\%$, $30\%$, or clade-holdout?*

#### Exact Multi-Channel Masking Specification (`Protocol H`)
For every outer evaluation fold on a target genome $\mathcal{G}_{\text{test}}$ at sequence-identity ceiling $\theta \in \{30\%, 50\%, 70\%\}$ (measured over aligned residues with $\ge 80\%$ query/target coverage **plus** an $\ge 80\%$ catalytic-domain local alignment screen to catch domain fusions) and taxonomic block $B_{\text{tax}} \in \{\text{species}, \text{genus}, \text{family}\}$:
We define the **Forbidden Donor Set** $\mathcal{F}(\mathcal{G}_{\text{test}}, \theta, B_{\text{tax}})$ containing all proteins in $\mathcal{G}_{\text{test}}$, all proteins in clade $B_{\text{tax}}$, and all reference proteins sharing $\ge \theta$ sequence or catalytic-domain identity with any test target. **Every upstream predictor and prior channel must enforce $\mathcal{F}$ before inference, calibration, or integration:**

| Evidence Channel | Strict Multi-Channel Masking Rule (`Protocol H`) | Permissive Deployment Track (`Protocol D`) |
|---|---|---|
| **1. Sequence Homology (`DIAMOND`, `MMseqs2`)** | Exclude all target records, exact duplicates, accession aliases, and any reference protein in $\mathcal{F}(\mathcal{G}_{\text{test}}, \theta, B_{\text{tax}})$ prior to BLAST/DIAMOND search. | Exclude only exact test proteome accessions (`Leave-Proteome-Out`). |
| **2. Domain & Family Rules (`InterProScan`, `Pfam`, `HAMAP`, `PANTHER`, `TIGRFAM`)** | Separate **structural domain presence** (allowed as an unlabelled architectural feature) from **curated function-transfer rules** (`InterPro2GO`, `HAMAP`/`TIGRFAM` EC rules). In strict `H`, disable direct family-to-GO/EC rule transfer if the family seed/rule provenance includes $\mathcal{F}$, or map domains only through fold-retrained non-leaking calibrators. | Full `InterPro2GO` and member-database functional rules enabled. |
| **3. PPI Networks (`STRING net`, `net_bridge`, `net_union`)** | Remove query species nodes from STRING (simulating a newly sequenced genome); in `net_bridge`, exclude all bridge donors in $\mathcal{F}$ and restrict STRING to experimental/co-expression channels (excluding `textmining` and `database`/KEGG-transfer edges). | Standard `net_union` excluding self-labels. |
| **4. Literature / Text (`BM25 lit`)** | **Disabled** in strict sequence-only `H` (or restricted to sequence-derived hits outside $\mathcal{F}$ with query gene names/descriptions stripped). | Pre-$t_0$ SwissProt text index with test species removed. |
| **5. Supervised PLM Heads & Catalyst Centroids (`ESM2-MLP`, `ProstT5`, `CLEAN`, `PLMCatalystCentroidRanker`)** | Rebuild all EC/GO catalyst centroids (`PLMCatalystCentroidRanker`), kNN indices, and supervised linear/MLP heads ** strictly excluding $\mathcal{F}$ from the training/centroid pool**. (Where external frozen checkpoints like `ProteInfer` or `DeepProZyme` cannot be retrained for every $\theta$, evaluate them on their documented pre-cutoff splits and report their exact training-set sequence-identity overlap per target). | Frozen public release checkpoints; centroids built on `Leave-Genome-Out` panel. |
| **6. Genomic Context & Pathway Priors (`RLGC`, `CoherencePrior`, `StoichiometricLocusRefiner`)** | Test genome coordinates (`.gff`) are allowed. For **Leave-Reaction-Out (LRO)**, strip all functional claims (sequence, domain, PLM, cached posteriors) on the target reaction's true encoding gene(s) (and all organism-specific isozyme realizations for reaction-recovery tasks), while retaining non-target neighbor annotations inferred under `Protocol H`. | Target reaction claims masked; neighbor claims from `Protocol D`. |
| **7. Stage 3 Agent (`gspa-agent`)** | Agent tools may only query the `Protocol H`-masked local evidence tables and ontology constraints; no web search, no external database lookup, and no gene symbol / description strings. Compare directly against a non-agentic GBDT/linear reranker given the exact same masked feature table. | Same local tools over `Protocol D` evidence. |

---

### 3.3 Delineated Benchmark Matrix: What Each Benchmark Tests & What Is Biologically Possible (Oracle Ceilings)

No single benchmark tests every aspect of function and phenotype prediction. We organize the evaluation around **four distinct, non-overlapping scientific claims** across **eight delineated benchmark suites**:

| Benchmark Suite | Ground Truth & Cohort Size | Leakage Control (`T` / `H` / `P`) | Which GSPA Mechanisms Can Legitimately Fire? | Biological & Mathematical Recoverability Ceiling (What Is Possible vs. Impossible) |
|---|---|---|---|---|
| **Suite 1: Prospective CAFA Temporal Holdout (`NK` & `LK`)** | True archive snapshot diff $\Delta E = E_{t_1} \setminus E_{t_0}$ (`EXP, IDA, IPI, IMP, IGI, IEP` only; GAF date-bug fixed). Reported separately for **Protein-NK** ($n \approx 251\text{--}500$), **Aspect-NK**, and **LK completion**. | **Strict `T`**: All databases, STRING, SwissProt text, and integrator weights frozen at $\le t_0$. Out-of-time evaluation on $(t_0, t_1]$. | • Multi-modal stacker (`DeepGO-PlusPlus`)<br>• `CoherencePrior` (process `has_part` + Complex Portal `AND`)<br>• Stage 2 `TaxonConsistencyMaxSat` + `AnatomicalConstraints` (eukaryote body plans)<br>• `RLGC` & `FBA` *only* for targets mapped to complete assemblies. | • **Possible:** Correcting cross-clade/anatomical false positives ($mi \downarrow$), boosting borderline BP/CC terms via complex/pathway partners.<br>• **Impossible:** FBA cannot affect non-metabolic signaling/structural proteins; incomplete experimental GOA means precision is bounded by unassayed true functions. |
| **Suite 2: Sequence-Blocked Reference Proteomes (25 Genomes)** | **19 Prokaryotes/Archaea + 6 Eukaryotes** (*S. cerevisiae*, *S. pombe*, *C. elegans*, *D. melanogaster*, *A. thaliana*, *M. musculus*). Primary: experimental GOA + SwissProt EC. Diagnostic: full GOA. | **Strict `H`** at $\theta \in \{30\%, 50\%, 70\%\}$ + Genus/Family holdout across all channels (`DIAMOND`, domain rules, `STRING`, `PLM` centroids). | • All Stage 1 priors (`Essentiality`, `Coherence`, `Complex AND`, `GenomicContext`, `Anatomical`)<br>• Stage 2 SAT/MaxSAT + `StoichiometricLocusRefiner`<br>• Stage 3 Agent on borderline anomalies. | • **At $>70\%$ ID:** Homology already saturates MF; GSPA mainly removes taxon/anatomical violations and improves BP/CC coherence.<br>• **At $<30\text{--}50\%$ ID:** Sequence homology degrades sharply; operon synteny, complex subunit co-occurrence, and FBA repair provide the primary recovery lift. |
| **Suite 3: Multi-Genome LRO — Stratum A (`Core Biomass Essential`)** | **15 diverse prokaryotic genomes** (~150–250 reactions/genome): reactions in the high-confidence GSMM that carry mandatory flux ($v_j \ge \epsilon_{\text{flux}}$) and whose knockout abolishes minimal-medium biomass ($\Phi_{\text{bio}}(y_{-j}) < \tau_{\text{bio}}$). | **Strict `H`** ($\theta \in \{30\%, 50\%, 70\%\}$) + complete masking of all functional claims on the true encoding gene(s) $p^*$ (and isozymes). Fold-rebuilt PLM centroids. | • **Stage 2 `StoichiometricLocusRefiner` (Minimal-Biomass FBA Repair)** detects the exact missing reaction $j$<br>• **Backward Locus Attribution $r_j(h)$** (`RLGC` operon density + `PLM` centroid + domain compatibility)<br>• Stage 3 Agent disambiguation. | • **Reaction Ceiling:** Bounded by whether reaction $j$ is non-redundant for biomass in the universal network (`Task Informativeness Ceiling`).<br>• **Locus Ceiling:** FBA identifies *that* reaction $j$ is needed ($a_j^\star=1$); `RLGC` + `PLM` + domains identify *which gene* $p^*$ encodes it (`Candidate Recall` $\times$ `Conditional Ranking`). |
| **Suite 4: Multi-Genome LRO — Stratum B (`Condition-Specific Catabolic & Transport`)** | **Same 15 genomes** (~100–200 reactions/genome): alternative carbon/nitrogen/energy assimilation and electron-acceptor reactions that carry **zero flux on minimal glucose** ($v_j = 0$) but are essential on a documented substrate $c$. | **Strict `H`** + target gene masking. Evaluated under two sub-settings:<br>• **B1:** Minimal glucose medium only<br>• **B2:** Cognate substrate condition $c$ supplied as FBA task. | • **Under B1 (Glucose FBA):** FBA is provably blind ($v_j=0$); ONLY `RLGC` operon neighborhood (co-clustered transporters/regulators) + `CoherencePrior` + `PLM` can recover $p^*$.<br>• **Under B2 (Substrate-$c$ FBA):** `StoichiometricLocusRefiner` activates $a_j^\star=1$ and boosts $r_j(h)$. | **Crucial Methodological Proof:** Demonstrates empirically that minimal-glucose FBA (`B1`) acts as a rigorous negative control for catabolic enzymes, whereas supplying the cognate growth substrate (`B2`) unlocks stoichiometric repair on top of operon context! |
| **Suite 5: Multi-Genome LRO — Stratum C (`Secondary & Peripheral Metabolism`)** | **Same 15 genomes + MIBiG/antiSMASH characterized clusters**: secondary metabolite (PKS/NRPS/terpene/siderophore) and uncoupled peripheral enzymes. | **Strict `H`** + exclusion of gene-cluster labels that directly state the EC. | • `RLGC` biosynthetic gene cluster (BGC) co-localization<br>• `CoherencePrior` (module completeness)<br>• Domain architecture (`inspect_candidate_architecture` in Stage 3). | **Impossible by Biomass FBA:** Secondary metabolites are not required for biomass; claiming biomass-MILP recovery here (as METEOR did with 0.68–1.00 FPR) is an artifact. Recovery is strictly bounded by **BGC synteny, module coherence, and domain/PLM specificity**. |
| **Suite 6: `Price-149` Re-Audited Enzyme Benchmark** | **All 149 RB-TnSeq proteins (56 ECs)** mapped to their exact source genomes, locus tags, and **Price et al. 2018 mutant-fitness growth substrates**, stratified into Core vs. Condition-Catabolic vs. Peripheral. | Crossed **`T` (pre-2018 snapshot)** and **`H` ($\theta \in \{30\%, 50\%, 70\%\}$)** vs. contemporary baselines (`CLEAN`, `DeepProZyme`, `EnzBERT`, `GraphEC`, `DeepGO-PlusPlus`). | • `RLGC` operon/transporter synteny + `CoherencePrior`<br>• **Condition-specific `StoichiometricLocusRefiner`** using the documented RB-TnSeq carbon/nitrogen source as task $\Phi_c$. | Separates (1) baseline PLM/sequence discrimination, (2) operon/module context lift, and (3) substrate-specific FBA repair lift, while reporting all 149 proteins with zero silent sample drops. |
| **Suite 7: `BacDive` Track 1 (Blind Genome-to-Phenotype Prediction — Flagship Impact)** | **1,000-strain primary RefSeq panel** (and full ~5,000-strain QC panel from `11,374` genome-linked strains in `bacdive_dump.json`): experimentally measured `+` and `-` carbon, nitrogen, energy, assimilation, fermentation, acid production, and 101 EC assays. | **Strict `P` (100% of test-strain phenotypes hidden)** + **Genus-Blocked & Family-Blocked splits** + **`Protocol H` / `D`**. Reconstruction is frozen *before* any test assay is probed. | • Stage 1 + Stage 2 (`StoichiometricLocusRefiner` under clade-viable minimal prior)<br>• **Readout A:** Diagnostic FBA flux probe $\Phi_c(y)$ per substrate<br>• **Readout B:** GPR catalytic + transporter completeness. | • **Strict Prediction–Intervention Boundary:** Never gap-fill toward the test assay.<br>• **Ceiling:** Bounded by ChEBI-to-exchange reaction coverage, transporter annotation, and transcriptional regulation (when a gene is present in the genome but unexpressed in the assay). |
| **Suite 8: `BacDive` Track 2 (Partial-Phenotype Assisted Reconstruction & Gene Attribution)** | **Same held-out BacDive strains**: reveal 30% of assay blocks (grouped by experimental panel/substrate family to prevent correlated leakage), predict the hidden 70% of assays. | **Strict `P` on hidden 70% blocks** + panel-level block splitting. | • `StoichiometricLocusRefiner` uses the revealed 30% positive growth/utilization tasks $\mathcal{C}_{\text{rec}}$ to repair missing transport/catabolic steps and assign gene responsibilities $r_j(h)$, then probes the hidden 70%. | Tests whether **phenotype-conditioned FBA gap-filling + locus attribution** transfers to unrevealed substrates/enzymes beyond Track 1, stratified by metabolic network distance to revealed substrates. |

---

### 3.4 Three Distinct Problem Formulations for Leave-Reaction-Out (LRO) & Unbiased Strata
Following our joint review with Astra, we eliminate two subtle biases in the Phase 12 LRO prototype (`benchmark/leave_reaction_out/BASELINE_REPORT.md`):
1. **Remove Favorable Pre-Filtering from the Main LRO Cohort:**
   - Earlier LRO runs filtered for $\ge 2$ local catalyzed neighbors within 20 kb and $\le 3$ paralogs (`BASELINE_REPORT.md` line 6), which artificially favors local synteny methods and excludes dispersed pathways and large paralog families.
   - **Fix:** The main 15-genome LRO benchmark includes **all verified single-gene and multi-subunit reactions** (dispersed pathways, 0–1 local neighbors, $\ge 4$ paralogs, and heteromeric complexes), and reports neighbor count ($0\text{--}1$, $2$, $3\text{--}4$, $\ge 5$), paralog count ($1$, $2\text{--}3$, $\ge 4$), and sequence identity ($<30\%, 30\text{--}50\%, 50\text{--}70\%, >70\%$) as **post-hoc analysis strata** rather than inclusion filters!
2. **Distinguish Three LRO Task Modes Explicitly:**
   - **Mode LRO-1 (Known-Reaction Locus Assignment):** The missing reaction identity $j$ is given (e.g., from a pathway hole); all functional annotations on its true encoding protein(s) $p^*$ are masked. Measures pure **candidate retrieval (`Recall@K`) and locus ranking (`Hit@1`, `Hit@5`, `MRR`)**.
   - **Mode LRO-2 (Blind Stoichiometric Reaction & Locus Recovery):** Both the reaction $j$ and its encoding protein(s) $p^*$ are removed from the draft model. `StoichiometricLocusRefiner` must first **discover that reaction $j$ is missing** via positive-cost flux-certified FBA repair ($a_j^\star = 1$) under minimal biomass (Stratum A) or cognate substrate probes (Stratum B), and then **attribute $j$ to the correct locus $p^*$** via $r_j(h)$ (or abstain via $h_\varnothing$ when $p^*$ is not in $C_r$).
   - **Mode LRO-3 (Heteromeric Complex Subunit Completion):** Mask one or more obligatory subunits $r \in R_{j,k}$ of a multi-subunit enzyme complex (`AND` GPR rule) and measure whether GPR-aware `CoherencePrior` + `StoichiometricLocusRefiner` recovers the exact partner locus.
3. **Five-Step Unconfounded Ablation Ladder on Identical Folds:**
   Every LRO and `Price-149` experiment compares five retrained/calibrated models on the exact same genomes, candidate sets, and homology masks:
   1. `Sequence / Domain / PLM Only` (no genomic context, no FBA),
   2. `Context Only (RLGC + Operon + Pathway/Complex Coherence)` (trained without PLM/sequence features—never by zeroing features in a PLM-heavy tree),
   3. `Combined Stage 1 (Sequence/PLM + Context + Coherence)`,
   4. `Combined Stage 1 + Stage 2 StoichiometricLocusRefiner (Flux-Certified FBA + Responsibility r_j(h))`,
   5. `Combined Stage 1 + Stage 2 + Stage 3 Neuro-Symbolic Agent` (compared directly against a non-agentic reranker with identical evidence access).

---

### 3.5 Five Explicit Recoverability Ceilings ("What Is Possible" vs. Methodological Limits)
In both the Methods and Results, we decompose every benchmark gap against five explicit, measurable ceilings rather than claiming an ill-defined "topological oracle":
1. **Vocabulary Ceiling ($\mathcal{C}_{\text{vocab}}$):** Fraction of ground-truth GO/EC annotations present in the predictor/ontology output vocabulary.
2. **Stoichiometric Representation Ceiling ($\mathcal{C}_{\text{rxn}}$):** Fraction of assayed enzymes or BacDive substrates mapped to a mass-balanced reaction (and ChEBI exchange/transporter) in the universal reaction universe.
3. **Task Informativeness Ceiling ($\mathcal{C}_{\text{task}}$):** Fraction of target reactions for which the declared FBA condition $c$ (minimal biomass in Stratum A, or substrate probe $c$ in Stratum B / BacDive) creates a non-zero flux differential ($\Phi_c(y_{+j}) > \Phi_c(y_{-j})$). By mathematical necessity, **FBA can only provide affirmative evidence for reactions inside $\mathcal{C}_{\text{task}}$**; outside $\mathcal{C}_{\text{task}}$ (e.g., Stratum C secondary metabolism or unprobed catabolism), recovery relies 100% on sequence, domain, PLM, and operon/pathway context.
4. **Candidate Retrieval Ceiling ($\mathcal{C}_{\text{cand}}$):** Fraction of true loci $p^*$ included in the bounded candidate pool $C_r$ (`Recall@K`).
5. **Locus Identifiability Ceiling ($\mathcal{C}_{\text{id}}$):** Because two candidates $p_1, p_2 \in C_r$ that enable the exact same reaction $j$ have identical stoichiometric likelihood $P(D_{\text{flux}} \mid h_1) = P(D_{\text{flux}} \mid h_2)$, stoichiometry raises the probability that *some* candidate in $C_r$ catalyzes $j$ (lowering $r_j(h_\varnothing)$), while **ranking between $p_1$ and $p_2$ is governed strictly by the prior odds ratio $\frac{q_j(h_1 \mid X)}{q_j(h_2 \mid X)}$ from local genomic context (`RLGC`), PLM centroids, and domain architecture**.

---

## 4. Updated Six-Pillar Experimental Architecture for the Manuscript

1. **Pillar 1 — Dual-Axis Protein Function Benchmark: Prospective CAFA Temporal Holdout (`Suite 1`) & Multi-Channel Sequence-Blocked Proteomes (`Suite 2`):**
   - **Suite 1 (`Protocol T`):** Leak-free archive-diff CAFA6 temporal evaluation (`Protein-NK`, `Aspect-NK`, `LK`) with the GAF last-modified date bug fixed.
   - **Suite 2 (`Protocol H` at 30%, 50%, 70% ID + Genus/Family holdout, and `Protocol D`):** 25 reference proteomes (19 Bacteria/Archaea + 6 Eukaryotes: *S. cerevisiae*, *S. pombe*, *C. elegans*, *D. melanogaster*, *A. thaliana*, *M. musculus*).
   - Compare 8 base predictors (`DIAMOND`, `InterProScan`, `eggNOG-mapper v2`, `Bakta`, `PGAP`, `CLEAN`, `DeepProZyme`, `DeepGO-PlusPlus`) before and after GSPA Stages 1, 2, and 3 using `fmax_micro`, IA-weighted $F_{\max}^w$, $S_{\min}$, validation-fixed operating point $F_1$/MCC, and genome/clade-clustered bootstrap 95% CIs.
2. **Pillar 2 — Mathematical Invariants & Role-Stratified Mechanism Recovery (`Suites 3–6`: `Core Biomass` vs. `Catabolic` vs. `Secondary` LRO + `Price-149`):**
   - **Mathematical Invariant Audit (109-genome panel):** Prove `0.0%` zero-flux free riders (`StoichiometricLocusRefiner` vs. `METEOR`), proteome-size decoy invariance, and identical-medium FBA witness preservation.
   - **Role-Stratified 15-Genome LRO (`Suites 3, 4, 5`) & `Price-149` (`Suite 6`):** Evaluate `Mode LRO-1` (known reaction), `Mode LRO-2` (blind reaction + locus recovery), and `Mode LRO-3` (complex subunit `AND` recovery) across Stratum A (`Core Biomass`), Stratum B (`Condition-Specific Catabolic`, contrasting Glucose-only `B1` vs. Cognate-Substrate `B2`), and Stratum C (`Secondary/Peripheral`), plus all 149 `Price-149` proteins under `Protocol T x H`.
3. **Pillar 3 — Flagship External Impact: Reproducing Experimental Microbial Phenotypes in BacDive (`Suites 7 & 8`, `/home/leechuck/Public/software/bacdive/bacdive_dump.json`):**
   - Evaluate **Track 1 (Blind Genome-to-Phenotype Prediction)** and **Track 2 (30%-Revealed Phenotype-Assisted Reconstruction & Gene Attribution)** across the curated **1,000-strain RefSeq panel** (and ~5,000-strain extended panel) spanning carbon, nitrogen, energy, assimilation, fermentation, acidification, and 101 EC enzyme assays under Genus-Blocked and Family-Blocked splits.
   - Report both **all eligible assays** and the **mechanistically representable subset ($\mathcal{C}_{\text{rxn}}$)** using Dual Readouts (Readout A: condition-specific FBA flux probe; Readout B: GPR catalytic + transporter completeness), comparing against `CLEAN`, `DeepProZyme`, `EnzBERT`, `DeepGO-PlusPlus`, `PGAP`, `Bakta`, `eggNOG-mapper v2`, `CarveMe`, `gapseq`, `METEOR`, and a `Taxonomic Prevalence Baseline` (with clade-clustered bootstrap CIs).
4. **Pillar 4 — GAEF Genome-Scale Plausibility & Eukaryotic Body-Plan Decontamination:**
   - Evaluate Completeness, Process/Pathway/Complex Coherence, and Taxonomic + Anatomical/Developmental Consistency across the 25-genome prokaryote/eukaryote panel, demonstrating 100% elimination of cross-clade and impossible organ/tissue development calls (*heart development* in yeast/worm, *flower development* in bacteria/fungi) via MaxSAT + ancestor relaxation.
5. **Pillar 5 — Dedicated Evaluation of the Optional Stage 3 Neuro-Symbolic Agentic Layer (`gspa-agent` with `Qwen 3.8 27B` on `unimatrix01`):**
   - Evaluate Stage 3 under the exact same `Protocol H` / `Protocol P` evidence restrictions on the 120-anomaly benchmark, high-entropy LRO gaps, and ambiguous BacDive Track 2 repairs—comparing `Qwen 3.8 27B` directly against a non-agentic GBDT/linear reranker with identical feature access and ablating the Symbolic + Stoichiometric Guardrail (`0.0%` violation rate with guardrail vs. `22.5%` without).
6. **Pillar 6 — Concise Secondary Discovery Vignette: Prospective Application to Novel Extremophile Isolates (Empty Quarter):**
   - Present a compact 1-paragraph / 1-panel prospective case study applying the frozen, BacDive-validated GSPA pipeline to novel desert extremophile isolates from the Empty Quarter pangenome catalog, framing predicted stress-response, carbon-scavenging, and biosynthetic capabilities strictly as testable biological hypotheses while BacDive carries the statistical generalization claim.

---

## 5. Concrete Implementation Roadmap

### Phase 1: Build Immutable Benchmark Manifests & Multi-Channel Homology Masking Engine (`benchmark/manifests/`)
1. **Immutable Manifest & Provenance Builder:**
   - Generate versioned, checksummed manifests for `Suite 1` (snapshot-diff CAFA temporal `NK`/`LK`), `Suite 2` (25-proteome panel with $\theta \in \{30\%, 50\%, 70\%\}$ and genus/family Forbidden Donor Sets $\mathcal{F}$), `Suites 3–5` (15-genome LRO without neighbor/paralog pre-filtering, annotated with `Stratum A/B/C` and cognate substrates), `Suite 6` (`Price-149` mapped to assemblies, locus tags, and RB-TnSeq substrates), and `Suites 7–8` (`BacDive` deduplicated strain–assembly–assay–condition tuples with genus/family blocks).
2. **Multi-Channel Masking Pipeline (`apply_protocol_H.py`):**
   - Implement automated per-fold filtering of $\mathcal{F}(\mathcal{G}_{\text{test}}, \theta, B_{\text{tax}})$ across DIAMOND, InterPro/domain rules, STRING bridge donors, and fold-specific `PLMCatalystCentroidRanker` / supervised heads so `Protocol H` is 100% leak-free and reproducible.

### Phase 2: Implement `StoichiometricLocusRefiner` (Clean FBA-Based Gap-Filling & Locus Attribution)
1. **Sparse GPR & Reaction-Mechanism Compiler (`gspa-core/src/main/groovy/gspa/integration/fba/`):**
   - Map ModelSEED / BiGG / MetaCyc reactions with verified ChEBI stoichiometry, mass/charge balance filters, compartmentation, and explicit `OR`-of-`AND` GPR roles (`GprMechanismCompiler.groovy`).
   - Implement truncated-Poisson / top-$k$ candidate pooling (`SparseRoleAggregator.groovy`) combining Stage 1 calibrated posteriors, `ReactionLocalContextSuggester` (RLGC), and fold-masked `PLMCatalystCentroidRanker`, plus the explicit unresolved hypothesis $h_\varnothing$.
2. **Flux-Certified Positive-Cost Repair Solver & Lexicographic pFBA (`StoichiometricLocusRefiner.groovy` + Python HiGHS worker):**
   - Formulate the positive-cost repair MILP with flux-support binaries $b^+_{jc}, b^-_{jc}$ ($\epsilon_{\text{flux}} = 10^{-5}$) using open-source **HiGHS** (deterministic seed, sparse matrix construction $<500\text{ MB}$ RSS), followed by lexicographic pFBA, FVA, single-reaction knockout necessity checks, and top-$M$ integer-cut alternative repair enumeration.
   - Implement backward Bayesian responsibility attribution $r_j(h)$ (`LocusResponsibilityUpdater.groovy`) to update specific candidate proteins $p \in C_r$ and couple cleanly with Stage 2 `TaxonConsistencyMaxSat` / `AnatomicalConstraints`.

### Phase 3: Execute the 8 Benchmark Suites (`Suites 1–8`) & Update Manuscript (`paper/main.tex`)
1. **Run `Suites 3–6` (Role-Stratified 15-Genome LRO, Mathematical Invariants vs. METEOR, and `Price-149`) on `unimatrix01` / `ibex`.**
2. **Run `Suites 7–8` (BacDive Track 1 Blind & Track 2 Assisted Phenotype Reproduction across all baselines) on `unimatrix01` / `ibex`.**
3. **Run `Suites 1–2` (Snapshot-Diff Temporal CAFA & Multi-Channel Sequence-Blocked 25 Proteomes) and `Pillar 5` (Stage 3 Agent vs. Non-Agentic Reranker).**
4. **Update Manuscript, Figures & Standalone Release (`v2.1.0`):**
   - Integrate the new mathematical formulation, the 8-suite results with explicit recoverability ceilings, and the flagship BacDive phenotype reproduction figures into `paper/main.tex` and tag `v2.1.0`.

