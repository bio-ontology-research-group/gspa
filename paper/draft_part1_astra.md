## Plan confirmation

I agree with the consolidated narrative, subsection order, and paragraph-level progression. The three drivers form a coherent argument: complementary protein-level evidence establishes the posterior, genome-scale architecture and ontological constraints refine it, and stoichiometric plus guarded neuro-symbolic reasoning resolves the remaining structured ambiguities.

I would make two minor refinements:

1. The Results allocation contains 22 paragraphs, although the plan states 21. The allocation is \(3+4+2+3+3+2+2+3=22\). I retain all 22 because each addresses a distinct planned analysis.
2. The planned abstract contains more than four headline values. The draft below retains four numerical values: the complex-recovery comparison, the zero-flux fraction, and the final ambiguous-gap recovery.

## Title

```latex
\title{\textbf{\sffamily Genome-scale context, stoichiometric refinement, and neuro-symbolic reasoning for proteome function annotation}}
```

## Abstract

```latex
\begin{abstract}
\noindent
Current protein function predictors combine sequence homology, domain profiles, protein language models, and structural similarity, while complete genome assemblies also provide operon organization, complex composition, and metabolic connectivity that can inform each prediction. We present Genome-Scale Protein Annotation (GSPA), a staged framework that calibrates and pools complementary protein-level predictors, propagates whole-genome constraints across pathways and obligate complexes, enforces taxonomic and anatomical consistency, and couples reaction-local genomic context with positive-cost parsimonious flux balance analysis, competitive locus attribution, and ontology-guarded language-model reasoning. Across bacterial, archaeal, viral, and eukaryotic proteomes, whole-genome inference preserves the protein-level baseline, resolves taxonomic and anatomical contradictions, and raises held-out heteromeric subunit recovery from $56.0\%$ to $93.7\%$. Positive-cost flux refinement removes the $83.7\%$ zero-flux reaction padding produced by negative-cost gap-filling, improves fitness-validated enzyme identification and experimental phenotype prediction, and raises ambiguous metabolic-gap Hit@1 to $82.7\%$. These results show that protein similarity, genome architecture, and systems-level constraints provide complementary layers of functional inference; software, constraint resources, and benchmarks are available at \url{https://github.com/bio-ontology-research-group/gspa}.
\end{abstract}
```

## Section 1: Introduction

### Intro P1

```latex
Automated protein function annotation assigns biological functions to the protein-coding genes of an organism and supports comparative genomics, metabolic reconstruction, and systems biology \citep{radivojac2013large,zhou2019cafa}. Given proteins $\mathcal{P}=\{p_1,\ldots,p_N\}$, their chromosomal coordinates and strand organization $\mathcal{K}$, the organism's taxonomic lineage $\tau_{\text{tax}}$, and reference biochemical and ontological knowledge, genome-scale annotation infers a joint assignment of Gene Ontology (GO) terms and Enzyme Commission (EC) numbers across $\mathcal{P}$ \citep{thomas2019gene}. A coherent joint assignment represents core cellular functions, connects metabolic reactions into routes that support cellular physiology, completes obligate heteromeric complexes, and assigns terms consistent with the organism's taxonomic lineage and anatomical body plan. These properties relate the function of each protein to the organization of the complete proteome.
```

### Intro P2

```latex
The increasing availability of complete genomes makes this joint formulation broadly applicable. Short-read assemblies often divided genomes into contigs that retained individual open reading frames while fragmenting operons, repeated regions, and chromosome-scale reaction inventories. Long-read sequencing now spans many of these regions and supports complete or near-complete genome assembly \citep{amarasinghe2020longread,logsdon2020longread}. Repeat-aware genome and metagenome assemblers extend this capability to cultivated isolates and metagenome-assembled genomes \citep{kolmogorov2020metaflye}. Gene order, co-directional operons, pathway complements, and complex composition can therefore serve as routine inputs to protein function annotation.
```

### Intro P3

```latex
Protein-level prediction methods derive complementary representations of a sequence. Alignment methods transfer functions from homologues \citep{buchfink2021sensitive,li2021refseq,schwengers2021bakta}; profile hidden Markov models identify conserved domains \citep{eddy2011accelerated,jones2014interproscan}; protein language models encode sequence context and evolutionary regularities \citep{lin2023evolutionary,elnaggar2021prottrans,kulmanov2024deepgose,yu2023enzyme}; and structure comparison identifies related folds \citep{gligorijevic2021structure,vankempen2024foldseek}. Multi-source systems such as GOLabeler and NetGO combine these representations and established their value in CAFA evaluations \citep{you2018golabeler,you2019netgo,wang2023netgo3,radivojac2013large,jiang2016expanded,zhou2019cafa}. Their complementary errors yield a more informative initial ranking than an individual representation. This ranking provides the foundation for conserved proteins, while independent per-protein scoring leaves divergent pathway enzymes, isolated complex subunits, and distant cross-clade transfers for genome-scale refinement.
```

### Intro P4

```latex
Comparative genomics adds information from the chromosomal context of a protein. Conserved neighborhoods, co-directional operons, and phylogenetic profiles associate genes that participate in the same pathway or molecular complex \citep{dandekar1998conservation,overbeek1999use,pellegrini1999assigning,price2005novel,szklarczyk2023string}. Genomic language models and neighborhood-aware predictors learn related representations from contiguous gene windows \citep{hwang2024genomic,cornman2024omg,toibazar2024context,wiatrak2025contextualised}. These representations can distinguish paralogues and identify co-regulated partners because the candidate locus is interpreted within its local gene neighborhood. A local window covers one part of genome organization; core essentiality, dispersed complex stoichiometry, and formal taxonomic or anatomical consistency require constraints over the complete proteome.
```

### Intro P5

```latex
Constraint-based metabolic reconstruction supplies a second source of genome-scale evidence. Pathway Tools, ModelSEED, CarveMe, gapseq, and GapMind place predicted enzymes in stoichiometric reaction networks and use flux balance analysis to identify reactions needed for steady-state growth \citep{green2004bayesian,seaver2021modelseed,machado2018fast,zimmermann2021gapseq,price2020gapmind}. METEOR was our initial approach for coupling metabolic gap-filling to protein-level enzyme prediction through reaction-indicator mixed-integer linear programming \citep{niu2026meteor}. Its reaction-level objective assigns negative addition costs to candidate reactions above a score threshold, which rewards the inclusion of reactions irrespective of their steady-state flux. Uniform transfer of reaction support to candidate proteins also distributes evidence across paralogues. These properties motivate flux-certified reaction selection and competitive attribution of each selected reaction to a genomic locus.
```

### Intro P6

```latex
Ontological constraints define which functional assignments are admissible for an organism. GO taxon constraints encode clade restrictions directly \citep{deegan2010formalization}, while logical definitions connect GO biological processes to anatomical classes in Uberon, the Plant Ontology, and the Cell Ontology \citep{mungall2012uberon,cooper2013plant,diehl2016cell}. These cross-ontology relations allow an anatomical restriction to constrain a process term whose GO record carries no direct taxon rule. Deletion-based post-filters remove a forbidden term together with biological information represented by its valid ancestors, whereas semantic relaxation can retain the nearest admissible process. After these deterministic operations, tool-using language agents can synthesize domain, neighborhood, pathway, and ontology evidence for residual ambiguous loci \citep{wang2025_geneagent,caufield2024_spires}. Deterministic symbolic verification can then test every proposed assignment against the formal constraints.
```

### Intro P7

```latex
We developed Genome-Scale Protein Annotation (GSPA) to combine these sources of evidence in their dependency order (Figure~\ref{fig:architecture}). Stage 0 calibrates and pools complementary per-protein predictors. Stage 1 updates their posteriors using whole-proteome priors for essential functions, pathways, complexes, taxonomic consistency, anatomy, metabolic gaps, and operon context, and Stage 2A applies deterministic constraint enforcement with semantic ancestor relaxation. Stage 2B aggregates isozyme and subunit mechanisms, selects positive-cost repairs by lexicographic parsimonious flux balance analysis \citep{lewis2010omic}, and assigns reaction responsibility competitively among candidate loci and an unassigned-gene hypothesis. Stage 3 applies selective tool-mediated language-model reasoning under a symbolic ontology guardrail. We evaluate whether these stages improve protein-level prediction, structural proteome coherence, complex and pathway completion, experimentally supported enzyme and phenotype prediction, and resolution of residual ambiguities across bacteria, archaea, viruses, and eukaryotes. Software, constraint resources, and benchmark manifests are available at \url{https://github.com/bio-ontology-research-group/gspa}.
```

## Section 2: Results

### 2.1 Complementary predictors establish the protein-level baseline

#### Results 2.1 P1

```latex
We first asked whether calibrated integration of complementary protein-level predictors improves annotation accuracy under a common ancestor-propagated evaluation. We compared GSPA with its individual sequence, domain, structure, and protein-language-model channels and with production annotation pipelines across the reference panels. On the prokaryotic panel, GSPA reached micro $F_{\max}=0.8421$, compared with $0.7820$ for NCBI PGAP (Figure~\ref{fig:pillar12}d). The expanded neural panel showed the same ordering across the additional channels and genomic language models (Supplementary Table~S2). Complementary representations therefore establish a more informative protein-level posterior for the genome-scale stages.
```

#### Results 2.1 P2

```latex
We next tested this posterior on proteins annotated after the training cutoff, with particular attention to the sequence-identity twilight zone. Stage 0 provided the leading overall result on the CAFA temporal holdout, while the first whole-genome refinement iteration raised macro $F_{\max}$ in the below-$30\%$ identity subset from $0.6723$ to $0.6775$ (Supplementary Table~S3). Protein-level pooling therefore accounts for most performance on conserved sequences, and genome context contributes where direct similarity provides less decisive evidence.
```

#### Results 2.1 P3

```latex
We then evaluated whether the same integration stabilizes EC prediction under temporal and donor-homology masking. On the post-cutoff enzyme benchmark, Stage 0 pooling reached $F_{\max}=0.4381$, and genome-scale refinement moved additional correct enzymes into the top-ranked candidate set (Supplementary Table~S4). Across the bacterial homology-masking panel, multi-channel consensus combined with stoichiometric filtering produced positive net Top-1 gains in every predictor-by-masking comparison, reaching 155 additional correct CLEAN assignments in the largest gain (Supplementary Table~S5). Multi-predictor agreement therefore provides a reliable enzyme posterior for subsequent structural refinement.
```

### 2.2 Whole-genome priors improve structural coherence while preserving accuracy

#### Results 2.2 P1

```latex
Having established the protein-level posterior, we tested whether whole-genome priors improve structural coherence while preserving protein-centric accuracy. We performed a systematic ablation across bacterial and archaeal reference genomes, using one genome for hyperparameter calibration and the remaining genomes for frozen evaluation. Full refinement retained mean CAFA $F_{\max}$ at approximately $0.80$ and increased the composite structural audit score to $0.8643$ (Figure~\ref{fig:pillar12}a,b and Table~\ref{tab:ablation}). Whole-genome priors therefore concentrate their changes on structurally constrained and lower-confidence loci around the conserved protein-level backbone.
```

#### Results 2.2 P2

```latex
We separated the priors to determine which proteome property each one changes. Essentiality increased the representation of core functions, pathway and gap priors connected metabolic routes, complex coherence completed multi-subunit assemblies, and operon context linked these changes to candidate loci. Together, these priors raised mean pathway coherence from $0.4975$ to $0.5146$ (Table~\ref{tab:ablation}). The targeted responses support the interpretation that each prior updates the genome property it represents.
```

#### Results 2.2 P3

```latex
We next tested the taxonomic and anatomical priors under controlled cross-clade perturbation. High-scoring eukaryotic developmental and clade-restricted GO claims were added to the prokaryotic proteomes before calibration and genome-scale enforcement. The ontological stages suppressed $97.6\%$ of the injected claims that remained above the Stage 0 decision threshold (Figure~\ref{fig:pillar12}c and Table~\ref{tab:ablation}). Formal organism-level constraints therefore identify cross-clade transfers that retain high protein-level similarity scores.
```

#### Results 2.2 P4

```latex
We then examined whether these structural effects persist at population scale and across domains of life. Frozen GSPA parameters were applied to a large RefSeq bacterial cohort and to curated bacterial, archaeal, eukaryotic, and viral proteomes. In the bacterial cohort, whole-genome refinement raised mean MetaCyc pathway coherence from $36.71\%$ to $66.75\%$ while eliminating the detected taxonomic violations (Supplementary Table~S6). Accuracy gains also remained positive for the bacterial and viral cohorts in the multi-domain analysis. Genome-scale refinement therefore applies consistently across diverse proteome sizes and taxonomic groups.
```

### 2.3 Operon architecture completes obligate heteromeric complexes

#### Results 2.3 P1

```latex
We asked whether complex stoichiometry and operon organization can identify a weakly supported subunit among domain-matched proteins. In leave-subunit-out trials, the true subunit and matched decoys received attenuated protein-level evidence, while the remaining partners defined the incomplete complex. Complex coherence identified the missing role, and operon context localized that role to the co-directional candidate locus, raising Hit@1 from $0.5600$ to $0.9371$ (Figure~\ref{fig:complex_anatomy}a and Table~\ref{tab:complex_lso}). Genome-wide complex composition and local chromosome architecture therefore provide distinct and complementary information.
```

#### Results 2.3 P2

```latex
The remaining errors were concentrated in complexes whose subunits occupy dispersed loci. We examined this setting in \textit{Synechocystis} sp.\ PCC 6803, where photosynthetic and respiratory complex genes frequently occur in separate chromosomal regions. Operon-aware inference reached Hit@1 of $0.6400$, and selective domain and pathway reasoning raised it to $0.6800$ (Figure~\ref{fig:complex_anatomy}b and Table~\ref{tab:complex_lso}). Local synteny resolves contiguous assemblies, while selective evidence synthesis contributes to complexes distributed across the chromosome.
```

### 2.4 Cross-ontology anatomical reasoning preserves valid ancestral biology in eukaryotes

#### Results 2.4 P1

```latex
We next tested whether ordinary multi-channel predictions contain anatomical contradictions in unperturbed eukaryotic proteomes. We evaluated six species against logical constraints connecting GO processes to Uberon, Plant Ontology, and Cell Ontology anatomy. Stage 0 contained $1{,}009$ endogenous taxonomic or anatomical violations, and deterministic cross-ontology enforcement reduced this count to zero while preserving mean protein-centric accuracy (Figure~\ref{fig:complex_anatomy}c and Table~\ref{tab:anatomy}). Anatomical reasoning therefore addresses errors produced by high-scoring cross-clade homology transfer.
```

#### Results 2.4 P2

```latex
We compared direct GO taxon rules, soft anatomical priors, and deterministic cross-ontology enforcement to identify the representation required for this correction. Direct GO constraints left $44.2\%$ of stress-test violations unresolved because many restrictions occur on anatomical classes linked through logical definitions. Deterministic enforcement across GO, Uberon, Plant Ontology, and Cell Ontology eliminated all detected violations (Figure~\ref{fig:complex_anatomy}c and Table~\ref{tab:anatomy}). Anatomical validity therefore depends on following the cross-ontology definitions that connect biological processes to clade-restricted structures.
```

#### Results 2.4 P3

```latex
Once a term is identified as anatomically inadmissible, its valid ancestors may still represent conserved biology. We compared deletion with relaxation to the nearest valid non-root GO ancestor against the original Swiss-Prot annotations. Semantic relaxation retained $95.1\%$--$98.7\%$ of the unconstrained annotation information content and improved micro $F_{\max}$ over deletion in all evaluated multicellular animals and vascular plants (Figure~\ref{fig:complex_anatomy}d and Table~\ref{tab:anatomy}). Cross-ontology constraints can therefore refine an assignment constructively by preserving its admissible biological-process ancestry.
```

### 2.5 Genomic architecture and language models localize missing pathway enzymes

#### Results 2.5 P1

```latex
We then asked whether reaction-local genomic context can identify a missing enzyme after sequence and embedding donors for its reaction have been masked. In strict leave-reaction-out trials, a uniform reaction-level boost preserved the sequence-based ordering of candidate paralogues. Adding operon kernel density and pathway coherence raised unconditional Hit@1 from $1.95\%$ to $23.44\%$ (Figure~\ref{fig:metabolism}a and Table~\ref{tab:metabolism}). Genomic proximity to pathway anchor genes therefore supplies locus-specific evidence that a reaction-level score cannot provide.
```

#### Results 2.5 P2

```latex
We next measured how genomic context complements protein-language-model reaction centroids. On \textit{E.~coli}, centroid similarity alone recovered $62.11\%$ of held-out enzymes at Rank 1, and adding reaction-local genomic context plus pathway coherence raised recovery to $63.67\%$ by resolving ties among paralogues (Figure~\ref{fig:metabolism}a,b and Table~\ref{tab:metabolism}). The combined model also transferred with frozen parameters to \textit{B.~subtilis}. Protein language models therefore provide the main candidate representation when reaction centroids are available, while genome architecture selects among candidates with similar representations.
```

#### Results 2.5 P3

```latex
Core flux on a specified medium covers only part of microbial metabolism, so we also evaluated local context in biosynthetic gene clusters and uncultured genomes. Operon synteny raised overall biosynthetic-cluster $F_1$ to $0.8959$ while retaining high specificity outside annotated clusters. Applied to metagenome-assembled genomes, the same context-aware framework identified 568 high-confidence candidates for missing pathway steps (Figure~\ref{fig:metabolism}d and Supplementary Table~S7). Local chromosome organization therefore extends locus attribution to secondary pathways and uncultivated lineages.
```

### 2.6 Positive-cost flux modeling eliminates reaction padding and attributes responsibility to gene loci

#### Results 2.6 P1

```latex
Reaction-local context identifies candidate loci, but a metabolic repair must also carry flux in the reconstructed network. We audited negative-cost mixed-integer gap-filling, positive-cost linear repair, and lexicographic parsimonious flux balance analysis across bacterial reconstructions. Under negative-cost gap-filling, $83.72\%$ of added reactions carried zero steady-state flux; positive-cost lexicographic refinement selected a mean of 5.70 flux-carrying, mass-balanced additions per genome (Figure~\ref{fig:fba_bacdive}a and Supplementary Table~S8). A positive addition cost therefore selects a parsimonious repair set supported by the flux solution.
```

#### Results 2.6 P2

```latex
We then tested whether assigning these flux-supported reactions competitively to gene loci improves experimentally supported enzyme identification. On the fitness-validated bacterial enzyme benchmark, single-predictor genomic and stoichiometric refinement increased MAPred Top-1 recovery from 50 to 56 proteins. Multi-predictor pooling followed by locus-specific refinement increased the corresponding recovery from 58 to 73 proteins (Figure~\ref{fig:fba_bacdive}b and Table~\ref{tab:price149}). Complementary protein predictors and competitive stoichiometric attribution therefore contribute additive information to enzyme discovery.
```

### 2.7 Systems-level annotations reproduce experimental bacterial phenotypes at scale

#### Results 2.7 P1

```latex
We next asked whether pathway-coupled and complex-aware posteriors reproduce experimental bacterial phenotypes under blind evaluation. GSPA pooled metabolic reconstructions with cross-assay pathway coupling and was compared with individual reconstruction systems on carbon-source growth and diagnostic enzyme assays. For carbon growth, the continuous GSPA posterior raised ROC-AUC from the leading binary reconstruction value of $0.8082$ to $0.8567$ while preserving its decision-level $F_1$ (Figure~\ref{fig:fba_bacdive}d and Table~\ref{tab:bacdive}). Systems-level coupling therefore converts pathway calls into calibrated physiological probabilities.
```

#### Results 2.7 P2

```latex
We then tested whether these probabilities complement supervised phenotype models across unseen genera and whether partial observations from a strain inform its remaining assays. Under genus-blocked cross-validation, the hybrid of GSPA posteriors and DiASPora reached carbon-growth $F_1=0.8710$. When a subset of each strain's phenotypes was supplied as evidence, the unsupervised pathway model reached held-out carbon-growth $F_1=0.8871$ (Table~\ref{tab:bacdive} and Supplementary Table~S9). Mechanistic pathway posteriors therefore complement statistical phenotype features and propagate observed traits to related unassayed pathways.
```

### 2.8 Selective guarded language-model reasoning resolves residual ambiguity

#### Results 2.8 P1

```latex
The deterministic stages leave a small set of pathway, complex, essentiality, and ontology anomalies that require evidence synthesis. Across the multi-domain proteome cohort, deterministic triggers selected less than $2.5\%$ of proteins in every domain. On the de-identified anomaly benchmark, the guarded multi-turn agent reached mean $F_1=0.9115$, with every accepted assignment passing the symbolic ontology checks (Figure~\ref{fig:agent}a,b,d and Table~\ref{tab:agent}). Selective invocation therefore confines language-model reasoning to residual cases while formal verification controls admissibility.
```

#### Results 2.8 P2

```latex
We tested whether tool-mediated reasoning can distinguish enzymes that remain tied after reaction-local genomic scoring. Across ambiguous leave-reaction-out gaps, the agent inspected domain architecture, pathway state, operon neighbors, and ontology relations and raised Hit@1 from $34.0\%$ to a three-seed mean of $82.7\%$ (Figure~\ref{fig:agent}c and Table~\ref{tab:agent}). Applied to prospective CAFA6 targets, guarded inference also achieved the leading macro $F_{\max}$ among the evaluated sequence and agent baselines (Supplementary Table~S9). Detailed evidence synthesis therefore contributes after deterministic priors have reduced the candidate set.
```

#### Results 2.8 P3

```latex
Finally, we compared the zero-shot agent with supervised tabular classifiers to identify their respective operating regimes. The guarded agent achieved $F_1=0.9115$, compared with $0.8283$ for the leading classifier trained on pre-tool summary features (Table~\ref{tab:agent}). Classifiers trained on labeled post-tool features provided the leading in-distribution binary promotion results and the lowest latency, while the zero-shot agent synthesized unlabelled evidence and proposed structured ancestral GO relaxations. Supervised classifiers therefore suit repeated binary decisions with representative labels, and guarded tool-mediated reasoning suits new clades and ontology operations that require traversal and term selection.
```

## Proposed figure captions

### Figure 1

```latex
\caption{GSPA integrates calibrated protein-level evidence with whole-genome and systems-level constraints. Stage 0 pools sequence, domain, structure, and protein-language-model predictions. Stage 1 updates proteome-wide posteriors using essentiality, pathway, complex, taxonomic, anatomical, metabolic, and operon evidence. Stage 2 applies semantic ancestor relaxation and flux-certified competitive locus attribution. Stage 3 uses selective tool-mediated reasoning, followed by deterministic ontology verification.}
```

### Figure 2

```latex
\caption{Protein-level integration and whole-genome prior ablation. (\textbf{a}) Structural coherence across successive evidence layers. (\textbf{b}) Per-genome change in the structural audit score. (\textbf{c}) Suppression of injected cross-clade claims. (\textbf{d}) Uniformly evaluated protein-function accuracy for GSPA, individual evidence channels, genomic language models, and production annotation systems.}
```

### Figure 3

```latex
\caption{Completion of obligate heteromeric complexes and enforcement of eukaryotic anatomical constraints. (\textbf{a}) Leave-subunit-out recovery after adding complex and operon evidence. (\textbf{b}) Recovery of dispersed complex subunits in \textit{Synechocystis}. (\textbf{c}) Violations remaining after direct, soft, and deterministic constraint enforcement. (\textbf{d}) Functional accuracy after deletion or semantic ancestor relaxation.}
```

### Figure 4

```latex
\caption{Localization of missing metabolic enzymes. (\textbf{a}) Leave-reaction-out recovery using sequence evidence, protein-language-model reaction centroids, pathway coherence, and reaction-local genomic context. (\textbf{b}) Recovery across metabolic roles and homology strata. (\textbf{c}) Transfer to \textit{Bacillus subtilis} and guarded re-ranking of ambiguous candidates. (\textbf{d}) Candidate pathway repairs in uncultured metagenome-assembled genomes.}
```

### Figure 5

```latex
\caption{Flux-certified locus refinement and experimental validation. (\textbf{a}) Reaction additions and flux activity under negative-cost gap-filling, positive-cost repair, and lexicographic parsimonious flux balance analysis. (\textbf{b}) Enzyme identification on the fitness-validated bacterial benchmark. (\textbf{c}) Locus recovery across metabolic roles and homology strata. (\textbf{d}) Prediction of experimental carbon-growth and diagnostic-enzyme phenotypes.}
```

### Figure 6

```latex
\caption{Evaluation of selective guarded language-model reasoning. (\textbf{a}) Resolution of de-identified pathway, complex, essentiality, and ontology anomalies. (\textbf{b}) Accepted ontology violations with and without symbolic verification. (\textbf{c}) Rank of the correct enzyme before and after agentic re-ranking. (\textbf{d}) Tool use and latency across evaluated species.}
```

## Proposed table captions

### Table 1

```latex
\caption{Systematic ablation of protein-level pooling, whole-genome priors, and deterministic symbolic enforcement across nine bacterial and archaeal reference genomes. Values summarize protein-centric accuracy, essential-function completeness, pathway and complex coherence, taxonomic consistency, structural audit scores, and suppression of cross-clade perturbations.}
```

### Table 2

```latex
\caption{Leave-subunit-out recovery of obligate heteromeric complex subunits across six prokaryotic genomes. Methods add complex coherence, operon context, singleton pruning, and selective guarded reasoning to the calibrated protein-level baseline.}
```

### Table 3

```latex
\caption{Taxonomic and anatomical constraint enforcement across six eukaryotes. The table compares direct GO taxon constraints, soft anatomical priors, deterministic deletion, and semantic ancestor relaxation using surviving violations, functional accuracy, recall, and annotation information content.}
```

### Table 4

```latex
\caption{Leave-reaction-out metabolic enzyme recovery in \textit{Escherichia coli}, out-of-sample recovery in \textit{Bacillus subtilis}, and re-ranking of ambiguous metabolic gaps. Results distinguish unconditional recovery from recovery among reactions whose true enzyme enters the candidate pool.}
```

### Table 5

```latex
\caption{Single-predictor and multi-predictor refinement on fitness-validated bacterial enzymes. Results compare each baseline with protein-level pooling, stoichiometric probing, competitive locus attribution, and full genome-scale refinement.}
```

### Table 6

```latex
\caption{Prediction of experimental carbon-source growth and diagnostic enzyme phenotypes in BacDive. The table reports blind pathway-based prediction, genus-blocked supervised evaluation, hybrid prediction, and phenotype-assisted inference on held-out assays.}
```

### Table 7

```latex
\caption{Resolution of de-identified genome-scale anomalies and ambiguous metabolic gaps. The table compares deterministic refinement, posterior thresholding, supervised tabular classifiers, fast agentic inference, guarded multi-turn reasoning, and symbolic-verification ablations.}
```