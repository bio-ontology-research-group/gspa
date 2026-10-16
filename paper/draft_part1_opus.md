## 1. Plan confirmation

I agree with the consolidated narrative. The three drivers (pooled per-protein evidence, whole-genome constraints, and systems-level and guarded reasoning) give the paper one argument. The 21-paragraph Results structure follows that argument, and every paragraph opens from a question the previous one left open.

When I checked the plan against `main.tex`, several key-number bullets claimed more than the data support. The draft below follows the data in those places:

1. **2.1 P3, the 95-genome EC filter.** The filter that gives significant gains in all 12 predictor–tier combinations is the three-predictor consensus filter. It accepts a top-1 change when the other two predictors support it. The plan calls it "consensus with stoichiometric filtering", which mixes it up with the single-predictor filter. On Holdout-125, pooling barely beats CLEAN (F_max 0.4381 vs 0.4360), so the draft gives that comparison.
2. **2.5 P1–P2.** Track A (cosine-masked, 11.33%) and the configuration without language-model information (coherence + RLGC, 23.44%) are separate settings; the plan merges them. The data also don't support "RLGC drives recovery in twilight-zone strata": RLGC alone scores 0% Hit@1 in the ≤50% identity stratum. What the data do show is that RLGC alone does best on catabolic reactions.
3. **2.6 P2, Price-149.** The plan calls the gains "additive", but that holds only for MAPred, DeepProZyme and CLEAN. For EnzBERT (62 after pooling vs 54 after full refinement) and GraphEC (48 vs 48), pooling accounts for the gain. The draft says so.
4. **2.7 P1, urease.** The 0.5610 and 0.7626 F_1 values come from different panels (5,026 vs 871 assays), and 0.7626 is gapseq's own value. The draft keeps only the comparison within one panel. It also says that gapseq's 0.8082 is a single-point AUC from binary calls.
5. **2.8 P3, the agent vs supervised classifiers.** On the 50 ambiguous metabolic gaps, cross-validated classifiers reach 90–96% Hit@1, above the agent's 82.7%. The draft states this. The agent's advantage holds on the 120-anomaly F_1 and on ancestral relaxations.
6. **2.4 P1 and 2.2 P4.** For Eukarya and Archaea, the confidence interval on the accuracy change includes zero. The draft reports a mean change there instead of "improves".

There are also some structural changes:

7. **Supplementary table numbering.** With the plan's numbering, the first citations would come out of order (S4 and S5 cited before S3). I propose this order: S1 benchmark suites, S2 hyperparameters, S3 Holdout-125, S4 95-genome EC, S5 cross-ontology constraint resources, S6 108-genome audit, S7 BacDive per-fold, S8 Stage 3 per-seed results, S9 Empty Quarter. Methods 4.1–4.10 would need the same renumbering.
8. **Figure 5 panel c (role-stratified LRO).** It duplicates Figure 4b and would be cited before Figure 5 is introduced. I dropped it, so Figure 5 becomes three panels: audit, Price-149, BacDive.
9. **No "Driver 1/2/3" labels in the prose.** They stay as planning labels; the text names the scientific content.
10. **Missing or venue-specific citations.** antiSMASH has no `.bib` entry. Qwen3.8 has none either (`yang2024qwen25` is Qwen2.5). I left METEOR's PSB 2027 acceptance to the bibliography entry. All other keys used below exist in `references.bib`.

---

## 2. Draft LaTeX (Part 1)

### Title

```latex
\title{Genome-scale context, stoichiometric refinement, and neuro-symbolic reasoning for proteome function annotation}
```

### Abstract (5 sentences, about 190 words)

```latex
\begin{abstract}
Protein function predictors combine sequence homology, domain profiles, protein language models, and structural similarity, but score each protein independently of the complete chromosomes that long-read sequencing and modern assembly now routinely recover.
We present Genome-Scale Protein Annotation (GSPA), a staged framework that pools complementary per-protein predictors, refines their output with whole-genome constraints from operons, essential functions, heteromeric complexes, and taxonomic and anatomical ontology axioms, and couples reaction-local genomic context with positive-cost parsimonious flux balance analysis, gene-locus attribution, and ontology-guarded language-model reasoning.
Across bacterial, archaeal, viral, and eukaryotic proteomes, whole-genome inference removes all detected taxonomic and anatomical violations, outperforms deletion of violating terms in multicellular animals and vascular plants, and raises recovery of held-out complex subunits (Hit@1) from 56.0\% to 93.7\%.
Positive-cost flux repair retains exclusively flux-carrying reaction additions, whereas negative-cost gap-filling adds 83.7\% zero-flux reactions, and improves the identification of fitness-validated enzymes and the prediction of experimental growth and enzyme phenotypes; guarded language-model reasoning on fewer than 2.5\% of proteins raises recovery of ambiguous missing enzymes from 34.0\% to 82.7\%.
Protein similarity, genome architecture, and systems-level constraints therefore provide complementary evidence for function annotation; GSPA is available at \url{https://github.com/bio-ontology-research-group/gspa}.
\end{abstract}
```

### Section 1: Introduction

**Intro P1: annotation as a joint, genome-scale problem**

```latex
\section{Introduction}

Genome projects require the function of every protein encoded in a genome, and these annotations are the basis of comparative genomics, metabolic modeling, and systems biology \citep{radivojac2013large,zhou2019cafa}. We consider the following problem. The input consists of the protein-coding genes $\mathcal{P} = \{p_1, \dots, p_N\}$ of an organism, their chromosomal coordinates and strand organization $\mathcal{K}$, the organism's taxonomic lineage $\tau_{\text{tax}}$, and reference biochemical and ontological knowledge bases. The output is an assignment of Gene Ontology (GO) \citep{thomas2019gene} terms and Enzyme Commission (EC) numbers to the proteins in $\mathcal{P}$. Because the proteins of an organism function together, this assignment is subject to constraints that range over many proteins at once. Core functions required for viability in the organism's domain of life, such as translation and the initiation of DNA replication, must be represented. Metabolic pathways must form stoichiometrically connected routes that support biomass synthesis \citep{osterman2003missing,karp2018how}. Obligate heteromeric complexes, such as ATP synthase or urease, must be annotated with all of their subunits \citep{meldal2022complexportal,giurgiu2019corum}. Finally, assigned terms must agree with the organism's taxon and with the anatomical structures present in its clade \citep{deegan2010formalization,mungall2012uberon}. An annotation that is accurate for each protein individually can therefore still be incoherent at the level of the organism.
```

**Intro P2: why whole-genome inference is timely**

```latex
These organism-level constraints can be evaluated only when the complete gene complement and gene order of a genome are known. Short-read sequencing of isolates and metagenomes produced fragmented assemblies from which open reading frames were called and annotated individually, and operon structure and complete reaction inventories were often unavailable. Long-read sequencing \citep{amarasinghe2020longread,logsdon2020longread} and repeat-graph assemblers for genomes and metagenomes \citep{kolmogorov2020metaflye} now routinely recover complete or near-complete chromosomes, both for cultivated isolates and for metagenome-assembled genomes (MAGs) of uncultured organisms. Gene order, co-directional operons, and the genome-wide complement of pathways and complexes are therefore available as standard inputs to function annotation.
```

**Intro P3: per-protein similarity and multi-source ensembles**

```latex
Most function prediction methods compute a representation of an individual protein and transfer annotations from proteins with similar representations. Sequence alignment identifies homologs of known function \citep{buchfink2021sensitive,li2021refseq,schwengers2021bakta}; profile hidden Markov models detect conserved domains \citep{eddy2011accelerated,jones2014interproscan}; protein language models embed sequences in a space in which functionally related proteins are close \citep{lin2023evolutionary,elnaggar2021prottrans,kulmanov2024deepgose,yu2023enzyme}; and structure-based methods compare folds \citep{gligorijevic2021structure,vankempen2024foldseek}. Because these representations capture different properties of a protein and make different errors, combining them improves prediction. GOLabeler \citep{you2018golabeler} and NetGO \citep{you2019netgo,wang2023netgo3} combine homology, domain, motif, and language-model evidence in a learned ranking, and they performed strongly in the Critical Assessment of Functional Annotation (CAFA) \citep{radivojac2013large,jiang2016expanded,zhou2019cafa}. For conserved proteins, these methods provide reliable predictions, and they remain the starting point of any annotation. However, they score each protein $p_i$ independently of the rest of the proteome $\mathcal{P} \setminus \{p_i\}$. As a consequence, divergent pathway enzymes fall below detection thresholds and leave pathway gaps, a domain match to one subunit of a heteromeric complex produces an isolated subunit call, and distant homology transfers developmental processes to clades that lack the corresponding anatomy.
```

**Intro P4: genomic context and chromosomal architecture**

```latex
Comparative genomics uses information that exists at the level of the genome. Conserved gene neighborhoods, co-directional operons, and phylogenetic co-occurrence profiles link proteins that act in the same pathway or complex \citep{dandekar1998conservation,overbeek1999use,pellegrini1999assigning,price2005novel,szklarczyk2023string}. Genomic language models \citep{hwang2024genomic,cornman2024omg} and neighborhood-aware predictors \citep{toibazar2024context,wiatrak2025contextualised} learn representations of a protein within a window of contiguous genes, which allows them to distinguish paralogs and to identify co-regulated partners. These representations cover a local segment of the chromosome. Constraints that span the whole genome, such as the presence of core essential functions, the stoichiometry of complexes whose subunits are encoded at distant loci, and the taxonomic and anatomical admissibility of GO terms, lie outside this segment.
```

**Intro P5: metabolic reconstruction and the limitations of METEOR**

```latex
Constraint-based metabolic reconstruction applies a different genome-scale constraint, the requirement that the metabolic network supports growth. Pathway Tools \citep{green2004bayesian}, ModelSEED \citep{seaver2021modelseed}, CarveMe \citep{machado2018fast}, gapseq \citep{zimmermann2021gapseq}, and GapMind \citep{price2020gapmind} build a stoichiometric reaction network from predicted enzymes and use flux balance analysis to identify reactions that are missing for steady-state growth. In previous work, we developed METEOR \citep{niu2026meteor}, which couples this gap-filling step with protein-level enzyme prediction: a mixed-integer linear program (MILP) selects reactions through binary indicators, and each selected reaction increases the scores of its candidate genes. This formulation has two limitations. First, METEOR assigns a negative objective cost to every candidate reaction whose score exceeds an addition threshold; the solver therefore improves its objective with each such reaction regardless of flux, and it adds hundreds of reactions per genome that carry zero steady-state flux. Second, METEOR transfers the support of a reaction uniformly to all of its candidate genes, which raises the scores of paralogs together and leaves operon neighborhoods, complex roles, and competition among candidate loci unused.
```

**Intro P6: cross-ontology axioms and guarded agentic reasoning**

```latex
Bio-ontologies formalize which functions are admissible in a clade. GO taxon constraints restrict terms to taxa or exclude them from taxa \citep{deegan2010formalization}, and logical definitions link GO biological processes to anatomical structures in the Uberon anatomy ontology \citep{mungall2011cross,mungall2012uberon}, the Plant Ontology (PO) \citep{cooper2013plant}, and the Cell Ontology (CL) \citep{diehl2016cell}. Through these definitions, a process such as \textit{flower development} inherits the taxon restrictions of the anatomical structure in which it occurs. Existing post-filters delete a violating term, and with it the ancestral process annotations that remain valid for the organism. A remaining class of errors concerns individual loci at which the available evidence is ambiguous. Tool-using language-model agents \citep{wang2025_geneagent,caufield2024_spires} can combine heterogeneous evidence, such as domain architectures, genomic neighbors, and pathway context, to resolve such cases, and a deterministic symbolic verifier can check their proposals against the same ontology axioms.
```

**Intro P7: GSPA, the evaluation questions, and availability**

```latex
Here we present Genome-Scale Protein Annotation (GSPA), a framework that addresses these limitations in four stages (Figure~\ref{fig:architecture}). Stage 0 calibrates the scores of complementary per-protein predictors and pools them into a probability for each protein--term pair, grouping correlated predictors so that they are counted once. Stage 1 updates these probabilities iteratively across the whole proteome with six genome-scale priors: essential functions, pathway and complex coherence, taxonomic consistency, anatomical constraints, metabolic gap-filling with reaction-local genomic context (RLGC), and operon coupling. Stage 2A resolves the remaining taxonomic and anatomical violations by relaxing each forbidden term to its nearest valid GO ancestor. Stage 2B addresses the two limitations of METEOR: it aggregates gene support through the isozyme and subunit roles of each reaction, repairs the network with a positive-cost lexicographic parsimonious flux balance analysis (pFBA) \citep{lewis2010omic}, and attributes each repaired reaction to gene loci that compete with an explicit hypothesis that no annotated gene catalyzes it. Stage 3 applies a tool-using language model to the fewer than 2.5\% of proteins that remain ambiguous after Stage 2, and a deterministic ontology guardrail rejects every proposal that violates the constraints. Stage 0 therefore contributes complementary per-protein evidence, Stages 1 and 2A contribute whole-genome constraints, and Stages 2B and 3 contribute systems-level and guarded reasoning beyond homology. We evaluate GSPA on bacteria, archaea, viruses, and eukaryotes (Supplementary Table~S1) and address five questions: whether pooling improves on individual predictors; whether whole-genome priors improve organism-level coherence while preserving accuracy; whether genome architecture and ontology axioms complete complexes and correct anatomical errors; whether positive-cost flux modeling with locus attribution improves enzyme identification and phenotype prediction; and when guarded language-model reasoning is useful. GSPA v2.1.0, the constraint resources, and the benchmarks are available at \url{https://github.com/bio-ontology-research-group/gspa}.
```

### Figure 1 caption

```latex
\caption{\textbf{Architecture of GSPA.} Stage 0 calibrates and pools per-protein predictions from sequence homology, domain profiles, protein language models, and structural similarity, or accepts the output of a single predictor. Stage 1 updates the probabilities of all protein--term pairs with six genome-scale priors that use genomic coordinates, draft metabolic networks, 26,911 GO taxon constraints, and 2,715 logical definitions linking GO terms to Uberon, PO, and CL classes. Stage 2A relaxes terms that violate taxonomic or anatomical constraints to their nearest valid GO ancestors and removes isolated complex subunits; Stage 2B repairs the metabolic network with positive-cost parsimonious flux balance analysis and attributes repaired reactions to gene loci. Stage 3 applies a tool-using language model with a symbolic ontology guardrail to the fewer than 2.5\% of proteins that remain ambiguous.}
```

---

### Section 2: Results

```latex
\section{Results}
\subsection{Complementary predictors establish the protein-level baseline}
```

**2.1 P1: pooling vs individual predictors and production pipelines**

```latex
We first asked whether pooling complementary per-protein predictors gives better annotations than each predictor alone and than production annotation pipelines. We compared GSPA, which pools the predictors in Stage 0 before refining them, with individual sequence, domain, language-model, and structure predictors, with NCBI PGAP \citep{li2021refseq}, InterProScan \citep{jones2014interproscan}, and metagenomic-deepFRI \citep{bezshapkin2026scaling}, and with the genomic language models gLM and gLM2 \citep{hwang2024genomic,cornman2024omg}. All methods were scored with the same ancestor-propagated protocol on a 13-genome and a 21-genome reference panel (Figure~\ref{fig:pillar12}d; all predictors in Supplementary Section~S3), and all GSPA hyperparameters were calibrated on \textit{Escherichia coli} K-12 and kept fixed for all other genomes (Supplementary Table~S2). On the 11 prokaryotes of the 13-genome panel, GSPA with three sequence and domain predictors reaches a micro $F_{\max}$ of 0.8421, compared with 0.7820 for PGAP. On the 21-genome panel, GSPA with all seven predictors reaches 0.8705, above the strongest individual predictor, ESM-2 nearest-neighbor transfer (0.8580). Combining complementary predictors therefore gives a stronger starting distribution than any of them alone. Because the reference panels consist of well-annotated genomes, we next tested whether this advantage holds for proteins annotated after the training data were fixed.
```

**2.1 P2: temporal holdout and the twilight zone**

```latex
We used a CAFA temporal holdout of 16,125 protein--aspect evaluations on proteins that received experimental annotations after the training cutoff, including 3,440 evaluations of proteins with less than 30\% sequence identity to the training data (Supplementary Section~S3). Stage 0 pooling reaches a macro $F_{\max}$ of 0.6929, compared with 0.6138 for ProteInfer \citep{sanderson2023proteinfer}, the strongest single predictor in this comparison. Stages 1 and 2 keep overall accuracy at the same level and raise $F_{\max}$ on the low-identity subset from 0.6723 to 0.6775. Pooling therefore determines most of the protein-level accuracy, and genome context contributes most where sequence similarity is weak. Enzyme annotation raises a related question: gap-filling methods propose changes to the top-ranked EC number of a protein, and we next asked which of these changes to accept.
```

**2.1 P3: post-cutoff enzymes and homology-masked EC consensus**

```latex
We evaluated EC number prediction in two settings. On 125 Swiss-Prot enzymes annotated after the training cutoff of the predictors (Supplementary Table~S3), pooling CLEAN \citep{yu2023enzyme}, EnzBERT \citep{buton2023enzbert}, and DeepProZyme \citep{kim2023dpz} in Stage 0 reaches an $F_{\max}$ of 0.4381, slightly above the strongest single predictor (CLEAN, 0.4360), and Stages 1 and 2 raise the fraction of enzymes with the correct EC number among the top five candidates to 44.0\%. Across 69,666 proteins from 95 bacterial genomes, we masked donor homologs at four sequence-identity thresholds (100\%, 70\%, 50\%, and 30\%) and examined the changes to top-ranked EC numbers proposed by METEOR (Supplementary Table~S4). Accepted without filtering, these changes reduce top-1 accuracy for DeepProZyme and EnzBERT. When GSPA accepts a change only if the two other predictors support the new EC number more strongly than the original one, the net number of correct top-1 predictions increases significantly in all twelve combinations of predictor and identity threshold, by up to 155 proteins for CLEAN ($p < 10^{-25}$, McNemar test). Agreement among complementary predictors therefore identifies which gap-filling proposals to trust. With this per-protein baseline in place, we next isolated what whole-genome priors add.
```

```latex
\subsection{Whole-genome priors improve organism-level coherence and preserve protein-level accuracy}
```

**2.2 P1: accuracy preservation and structural gain**

```latex
Whole-genome priors change the probability of an annotation according to the annotations of other proteins in the same genome, and such updates could lower the accuracy established in Stage 0. We tested this in an ablation of nine configurations on nine bacterial and archaeal genomes, with \textit{E.~coli} K-12 as calibration genome and the other eight genomes evaluated out of sample (81 runs; Figure~\ref{fig:pillar12}a,b and Table~\ref{tab:ablation}). The full Stage 1 and 2 model keeps protein-centric CAFA $F_{\max}$ at the Stage 0 level (0.8037 compared with 0.8041), slightly increases pair-level micro $F_{\max}$, and raises a composite structural audit score \citep{toonsi2025gaef} from 0.8459 to 0.8643 ($p = 0.0039$, Wilcoxon signed-rank test). About 85\% of annotated proteins in these genomes have high-confidence Swiss-Prot homologs, so Stage 0 determines most protein-level predictions, and the priors act on the remaining divergent and structurally constrained loci. We next asked which organism-level properties each prior changes.
```

**2.2 P2: contribution of individual priors**

```latex
Each Stage 1 prior encodes one organism-level property: the essentiality prior the presence of core functions, the coherence prior the completeness of pathways and complexes, the gap-filling and operon priors the connectivity of pathways, and the consistency and anatomical priors the admissibility of terms in the organism's taxon (Methods). Added individually, the essentiality, coherence, gap-filling, and operon priors change pathway coherence by small amounts; combined in the full Stage 1 model, they raise mean pathway coherence from 0.4975 to 0.5146 ($p = 0.0078$) and complex coherence from 0.8636 to 0.8852 (Table~\ref{tab:ablation}). The consistency and anatomical priors, which read the taxon constraints and cross-ontology definitions in Supplementary Table~S5, remove all 25 taxonomic violations present at baseline. Each prior therefore improves the property it encodes, and the priors together improve pathway connectivity beyond any single prior. These violations were already present in the predictions; we next tested whether the same constraints suppress errors introduced by cross-clade homology transfer.
```

**2.2 P3: suppression of injected cross-clade errors**

```latex
To simulate erroneous homology transfer, we added 900 eukaryotic developmental and clade-restricted GO annotations with raw scores between 0.55 and 0.85 to the nine prokaryotic proteomes (Figure~\ref{fig:pillar12}c and Table~\ref{tab:ablation}). Stage 0 calibration lowers most of these scores below the decision threshold of 0.50, and 125 injected annotations remain above it. The taxonomic and anatomical constraints suppress 122 of these 125 annotations (97.6\%). Organism-level constraints therefore remove cross-clade errors that pass protein-level calibration. Because the nine genomes are well-studied reference organisms, we next tested whether these effects extend to a large and taxonomically diverse set of genomes.
```

**2.2 P4: 1,000 bacterial genomes and 88 proteomes from four domains**

```latex
We applied GSPA with fixed hyperparameters to 1,000 RefSeq bacterial genomes and to 88 proteomes from Bacteria, Archaea, Eukarya, and viruses (Supplementary Section~S3). In the 1,000 bacterial genomes, Stage 0 makes a mean of 23,384 taxon-violating protein--term predictions per genome; Stages 1 and 2 remove all of them and raise MetaCyc pathway coherence from 36.7\% to 66.8\%, at 5.4~s per genome. In the 88 proteomes, the 95\% bootstrap confidence interval of the per-genome accuracy change lies above zero for Bacteria and viruses and includes zero for Archaea and Eukarya. Whole-genome refinement therefore scales to large genome collections, enforces taxonomic consistency in every genome, and improves accuracy in the two largest groups. Complex coherence records whether all partners of an assembly are annotated, and we next asked whether genome architecture identifies the gene that encodes a missing partner.
```

```latex
\subsection{Operon architecture completes obligate heteromeric complexes}
```

**2.3 P1: leave-subunit-out recovery**

```latex
Complex coherence identifies an assembly that lacks a subunit, but a subunit with weak homology competes with other proteins that match the same domains. We tested whether complex stoichiometry and operon context distinguish the true subunit from such competitors in 175 leave-subunit-out trials across six prokaryotic genomes. In each trial, we lowered the evidence for one subunit of a curated complex into the twilight-zone range and added five decoy proteins with comparable domain scores from unrelated loci, so that a random ranking places the true subunit first with probability $1/6$ (Figure~\ref{fig:complex_anatomy}a and Table~\ref{tab:complex_lso}). The complex coherence prior raises the probability of all domain-matching candidates equally and keeps Hit@1 at the baseline value of 0.560. Adding operon context raises Hit@1 to 0.937 (0.931 in the five genomes other than the calibration genome). Complex stoichiometry therefore identifies which assembly lacks a subunit, and operon architecture identifies the locus that encodes it. In some organisms, complex subunits are encoded far apart on the chromosome, and we examined such a case next.
```

**2.3 P2: dispersed complexes in *Synechocystis***

```latex
In the cyanobacterium \textit{Synechocystis} sp.\ PCC 6803, subunits of the photosystems, phycobilisomes, and NADH dehydrogenase are often encoded at separate loci (Figure~\ref{fig:complex_anatomy}b and Table~\ref{tab:complex_lso}). In the 25 trials in this genome, operon context raises Hit@1 from 0.52 to 0.64, a smaller gain than in the other five genomes. Adding Stage 2 and the Stage 3 agent, which inspects the domain architecture and pathway context of each candidate, recovers one further subunit (Hit@1 of 0.68) and raises Hit@1 across all six genomes to 0.943. Operon context therefore resolves complexes encoded in contiguous operons, and Stage 3 adds evidence when subunits are dispersed. In multicellular eukaryotes, the anatomy of the organism provides a further organism-level constraint, which we evaluated next.
```

```latex
\subsection{Cross-ontology anatomical constraints correct eukaryotic annotations and preserve valid ancestral terms}\label{sec:euk_results}
```

**2.4 P1: endogenous violations in eukaryotic proteomes**

```latex
We first asked whether unperturbed per-protein predictions contain anatomical and taxonomic errors in eukaryotes. We applied the seven-predictor Stage 0 model to 6,155 proteins from \textit{Saccharomyces cerevisiae}, \textit{Caenorhabditis elegans}, \textit{Drosophila melanogaster}, \textit{Physcomitrium patens}, \textit{Arabidopsis thaliana}, and \textit{Mus musculus}, and checked the predictions against GO taxon constraints and the logical definitions linking GO to Uberon, PO, and CL (Figure~\ref{fig:complex_anatomy}c and Table~\ref{tab:anatomy}). Stage 0 produces 1,009 violations, almost all in \textit{C.~elegans} and \textit{D.~melanogaster}, where processes are transferred from mammalian and flowering-plant homologs. GSPA removes all 1,009 violations, and the mean CAFA $F_{\max}$ across the six species changes from 0.7096 to 0.7117. Homology transfer therefore assigns organ-level processes of vertebrates and flowering plants to invertebrates, and we next asked which constraints detect these errors.
```

**2.4 P2: why cross-ontology axioms and deterministic enforcement are needed**

```latex
To compare constraint representations, we added 1,546 cross-clade developmental annotations to the six proteomes, which yields 1,452 violations, and removed them in three ways (Figure~\ref{fig:complex_anatomy}c and Table~\ref{tab:anatomy}). Direct GO taxon constraints remove 55.8\% of the violations and leave 642, because the taxon restrictions of these organ-level terms are defined on the linked Uberon, PO, and CL classes. The soft anatomical prior of Stage 1 removes almost all violations in four species and leaves 201 in \textit{C.~elegans} and 559 in \textit{D.~melanogaster}, where the homology scores of the transferred terms exceed 0.95. Deterministic enforcement in Stage 2A, which reads the anatomical restrictions through the cross-ontology definitions, removes all violations. Complete enforcement of anatomical constraints therefore combines the links from GO to anatomy ontologies with a deterministic step after the probabilistic update. Such enforcement can either delete a forbidden term or replace it with a valid ancestor, and we compared these two options next.
```

**2.4 P3: ancestor relaxation vs deletion**

```latex
When a protein of an organism without a heart is annotated with \textit{heart development}, its participation in a more general developmental process may still be correct. Semantic ancestor relaxation replaces the forbidden term with its nearest ancestor that is valid in the organism's taxon, whereas deletion removes the term together with the ancestors it implies. We compared both strategies against unmodified Swiss-Prot annotations (Figure~\ref{fig:complex_anatomy}d and Table~\ref{tab:anatomy}). Relaxation retains 95.1\% to 98.7\% of the mean information content of the unconstrained annotations and achieves a higher micro $F_{\max}$ than deletion in all four multicellular animals and vascular plants, with the largest difference in \textit{D.~melanogaster} ($+0.0157$). In yeast and moss, the injected organ-level terms fall on metabolic and housekeeping proteins; there, deletion achieves the higher $F_{\max}$, and relaxation achieves the higher recall and information content. Relaxation therefore preserves valid ancestral biology in organisms whose proteins participate in tissue and organ development. The constraints considered so far act on terms and complexes; we next turned to metabolism, where a missing enzyme must be located in the genome.
```

```latex
\subsection{Genomic context and protein language models locate missing pathway enzymes}
```

**2.5 P1: leave-reaction-out recovery without language-model information**

```latex
Metabolic gap-filling identifies the reactions a genome needs, but it leaves open which gene catalyzes each reaction. We tested whether reaction-local genomic context (RLGC), which scores a candidate gene by its chromosomal proximity to genes of adjacent pathway reactions, locates the enzyme of a reaction when sequence homologs are masked. In 256 leave-reaction-out trials in \textit{E.~coli} K-12, we removed the gene--reaction association of one reaction at a time, masked sequence homologs of the true enzyme, and ranked all candidate proteins without protein language-model information (Figure~\ref{fig:metabolism}a and Table~\ref{tab:metabolism}). Sequence-only ranking places the true enzyme first for 1.95\% of reactions, and METEOR's reaction-level support leaves this ranking unchanged because it adds the same increment to every candidate of a reaction. Combining RLGC with the Stage 1 pathway coherence prior raises Hit@1 to 23.44\%. The chromosomal neighborhoods of pathway partners therefore distinguish candidate genes to which reaction-level support assigns equal scores. Protein language models provide a second, alignment-free source of evidence, and we next asked how it combines with genomic context.
```

**2.5 P2: language-model centroids and genomic context across roles and species**

```latex
In a second setting, we retained ESM-2 embedding centroids of Swiss-Prot enzymes for each reaction and masked sequence alignments to the held-out reaction (Figure~\ref{fig:metabolism}a--c and Table~\ref{tab:metabolism}). In \textit{E.~coli}, the centroids alone place the true enzyme first for 62.11\% of reactions, and adding RLGC and pathway coherence raises Hit@1 to 63.67\% by separating paralogs with similar embeddings. Across metabolic roles, RLGC alone performs best on alternative catabolic reactions, whose genes are often clustered in operons (Figure~\ref{fig:metabolism}b). With hyperparameters fixed on \textit{E.~coli}, the combined model reaches a Hit@1 of 43.06\% on 144 reactions in \textit{Bacillus subtilis} 168. Language-model embeddings therefore provide most of the recovery when reference enzymes for a reaction exist, and genomic context resolves ties between paralogs. Both analyses concern primary metabolism in reference organisms; we next applied genomic context to secondary metabolism and to uncultured genomes.
```

**2.5 P3: biosynthetic gene clusters and 97 uncultured MAGs**

```latex
Biosynthetic gene clusters produce secondary metabolites and carry zero flux during growth on minimal glucose, so flux-based gap-filling assigns their genes zero support. On 106 biosynthetic gene clusters curated with antiSMASH, GSPA scores these loci by cluster synteny and raises the $F_1$ of cluster gene recovery from 0.7395 to 0.8959, at 99.6\% specificity on genes outside the clusters (Supplementary Section~S6). As an application to genomes without reference annotations, we ran GSPA on 97 MAGs of uncultured organisms from 16 bacterial and archaeal phyla (Figure~\ref{fig:metabolism}d). GSPA proposes 568 candidate genes for missing pathway steps at a posterior probability of at least 0.50, which raises mean draft pathway completeness from 0.511 to 0.656. These candidates are ranked hypotheses for experimental testing. Gap-filling in these analyses adds reactions to complete pathways, and we next asked which reactions such repairs add and whether they carry flux.
```

```latex
\subsection{Positive-cost flux repair adds flux-carrying reactions and attributes them to gene loci}
```

**2.6 P1: 108-genome stoichiometric audit**

```latex
Gap-filling adds reactions to a draft metabolic network until the network can produce biomass. We audited the added reactions in 108 bacterial reconstructions under three formulations: METEOR's negative-cost MILP, a positive-cost linear program, and the two-step lexicographic pFBA of Stage 2B (Figure~\ref{fig:fba_bacdive}a and Supplementary Table~S6). Negative-cost MILP adds a mean of 313 reactions per genome, of which 83.7\% carry zero steady-state flux. With a positive addition cost, every added reaction carries flux, and lexicographic pFBA reduces the number of added reactions from 81.3 under the positive-cost linear program to 5.7 per genome, all of them mass-balanced. A negative cost rewards the solver for each added reaction, so the MILP accumulates reactions regardless of flux; a positive cost selects a small repair in which every reaction contributes to growth. Each added reaction must still be assigned to a gene, and we next tested whether competitive attribution of reactions to gene loci improves enzyme identification.
```

**2.6 P2: Price-149 fitness-validated enzymes**

```latex
We evaluated enzyme identification on 148 proteins from 22 bacterial genomes whose functions were established by mutant fitness assays \citep{price2018mutant}, using six whole-genome enzyme predictors and two predictors applied to target sequences (Figure~\ref{fig:fba_bacdive}b and Table~\ref{tab:price149}). Of these proteins, 55 act in catabolic pathways that are inactive during growth on glucose, and glucose-based MILP gap-filling lowers $F_{\max}$ for three predictors. Applied to the output of a single predictor, GSPA raises the number of correct top-ranked EC numbers for MAPred \citep{rong2025mapred} from 50 to 56. Pooling all predictors in Stage 0 gives the largest gains for weak predictors, for example from 18 to 48 for GraphEC \citep{song2024accurately}, and full multi-predictor refinement raises MAPred from 58 after pooling to 73 (49.3\%). For EnzBERT and GraphEC, pooling accounts for most of the improvement (Table~\ref{tab:price149}). Pooling and locus-level stoichiometric refinement therefore both contribute to the identification of experimentally validated enzymes. These benchmarks evaluate individual proteins; we next tested whether GSPA's annotations predict the growth and enzyme phenotypes of whole organisms.
```

```latex
\subsection{Systems-level annotations predict experimental bacterial phenotypes}
```

**2.7 P1: blind phenotype prediction in BacDive**

```latex
BacDive records experimentally measured carbon-source utilization and diagnostic enzyme activities for thousands of bacterial strains \citep{reimer2022bacdive,schober2025bac}. We compared blind GSPA predictions, made with zero phenotype labels, with CarveMe, ModelSEED, gapseq, and METEOR on 2,332 carbon-growth assays in 529 strains, 17,496 enzyme assays in 3,017 type strains, and 50,075 assays of 10 diagnostic enzymes in 6,313 strains (Figure~\ref{fig:fba_bacdive}c and Table~\ref{tab:bacdive}). METEOR's gap-filling lowers precision on the diagnostic enzymes, including the heteromeric urease, for which it raises the scores of isolated subunits. GSPA pools the calibrated calls of the three reconstruction tools with pathway co-occurrence across assays. It reproduces the binary carbon-growth calls of gapseq, the most accurate reconstruction tool ($F_1 = 0.7952$), and its continuous probabilities reach a ROC-AUC of 0.8567, compared with 0.8082 for the single operating point of gapseq's binary calls. Calibrated pathway probabilities therefore rank strains by their capacity for growth more accurately than binary reconstructions, using zero training labels. We next asked whether these probabilities add information to supervised phenotype classifiers and whether the observed phenotypes of a strain inform its remaining predictions.
```

**2.7 P2: out-of-genus transfer and phenotype anchoring**

```latex
In five-fold cross-validation with each genus confined to a single fold, adding GSPA probabilities to the DiASPora classifier \citep{koblitz2025diaspora} raises out-of-genus $F_1$ to 0.8710 for carbon growth and 0.8713 for diagnostic enzymes, above all metabolic and supervised baselines (Table~\ref{tab:bacdive} and Supplementary Table~S7). In a second protocol, we revealed 30\% of each strain's assays and predicted the remaining 70\%. Conditioning GSPA's pathway probabilities on the revealed assays in an unsupervised Bayesian update raises held-out carbon-growth $F_1$ from 0.7922 to 0.8871, and combining the conditioned probabilities with DiASPora reaches 0.9225. Mechanistic pathway coupling and statistical phenotype models therefore contribute complementary information, and observed phenotypes propagate through the metabolic network to unassayed traits. After all deterministic stages, a small fraction of loci remains ambiguous, and we examined these loci last.
```

```latex
\subsection{Guarded language-model reasoning resolves residual ambiguity}
```

**2.8 P1: selectivity and resolution of 120 anomalies**

```latex
Stage 3 applies a language model to proteins that deterministic triggers flag after Stage 2: unresolved pathway gaps, incomplete complexes, borderline essential functions, and borderline ontology conflicts. In the 88 proteomes, these triggers select a mean of 1.84\% of proteins per proteome, and fewer than 2.5\% in every domain. To evaluate Stage 3 under controlled conditions, we built a benchmark of 120 anomalies from seven genomes, masked all gene symbols, shuffled candidate order, and ran the agent (Qwen3.8-27B with five read-only genomic and ontology tools and the symbolic ontology guardrail) with three random seeds (Figure~\ref{fig:agent}a,b,d, Table~\ref{tab:agent}, and Supplementary Table~S8). The guarded agent reaches a mean $F_1$ of 0.9115 with zero accepted ontology violations; removing the guardrail admits a cross-clade violation and lowers $F_1$ to 0.8515. Stage 3 therefore acts on a small fraction of loci, and the guardrail ensures that its outputs satisfy the ontology constraints. Missing metabolic enzymes form a second class of residual ambiguity, which we evaluated next.
```

**2.8 P2: tied metabolic candidates and prospective CAFA6 targets**

```latex
After RLGC scoring, 50 leave-reaction-out gaps in \textit{E.~coli} K-12 retain several candidates with similar scores. Re-ranking the top five candidates of each gap with Stage 3 raises Hit@1 from 34.0\% to a three-seed mean of 82.7\% (Figure~\ref{fig:agent}c and Table~\ref{tab:agent}). On 769 prospective CAFA6 targets, which are supplied as isolated sequences, GSPA with Stages 0, 2A, and 3 reaches a macro $F_{\max}$ of 0.5433, compared with 0.4957 for a classifier on language-model and structure-token embeddings and 0.5130 for an unguarded agent based on a 14-billion-parameter model. Once deterministic priors have narrowed the candidates, inspection of domain architectures and operon neighbors can therefore resolve ties between paralogs. A supervised classifier can learn from the same features that the tools return, and we compared the agent with such classifiers next.
```

**2.8 P3: when a zero-shot agent is preferable to a supervised classifier**

```latex
We trained logistic regression, random forest, gradient-boosted tree, and multilayer perceptron classifiers by five-fold cross-validation, using either three scores available before tool use or eight numerical features returned by the tools, the latter combined with the guardrail (Table~\ref{tab:agent}). On the 120 anomalies, the zero-shot agent ($F_1 = 0.9115$) outperforms classifiers trained on the three pre-tool scores ($F_1$ from 0.7500 to 0.8283), and classifiers trained on the eight tool features reach up to 0.9720. On the 50 metabolic gaps, cross-validated classifiers reach a Hit@1 between 90.0\% and 96.0\%, above the agent. On conflict cases, the agent proposes a mean of 13.7 valid ancestral GO relaxations and the classifiers propose at most two, because relaxation requires traversal of the ontology. Supervised classifiers are therefore preferable for binary promotion and ranking when labeled anomalies from the same distribution are available, and the zero-shot guarded agent is suited to clades that lack such labels and to proposing ancestral relaxations. An application of GSPA to five extremophile isolates from the Empty Quarter desert is described in Supplementary Section~S8 and Supplementary Table~S9.
```

---

### Figure and table captions

**Figure 2**

```latex
\caption{\textbf{Per-protein pooling and whole-genome priors on reference genomes.} (\textbf{a}) Structural audit score, essential function completeness, complex coherence, and pathway coherence across the nine ablation configurations, averaged over nine bacterial and archaeal genomes. (\textbf{b}) Gain in structural audit score from the Stage 0 baseline to the full Stage 1 and 2 model for each genome; the asterisk marks the calibration genome, \textit{E.~coli} K-12. (\textbf{c}) Injected cross-clade annotations remaining above the decision threshold in each genome under the Stage 0 baseline and with taxonomic and anatomical constraints (100 injected annotations per genome). (\textbf{d}) Ancestor-propagated micro and protein-centric $F_{\max}$ of GSPA, external annotation pipelines, individual Stage 0 predictors, and genomic language models on the 13-genome and 21-genome panels.}
```

**Table 1** (suggested row labels: *Stage 0 baseline; + essentiality prior; + coherence prior; + gap-filling prior; + operon context prior; + consistency and anatomical priors; all Stage 1 priors; Stage 1 + Stage 2A removal; Stage 1 + Stage 2A relaxation*)

```latex
\caption{\textbf{Ablation of Stage 1 priors and Stage 2A enforcement on nine reference genomes.} Rows add priors to the Stage 0 baseline one at a time, combine all six Stage 1 priors, and add Stage 2A enforcement by removal or by relaxation. Values are means over \textit{E.~coli} K-12 (calibration genome), \textit{E.~coli} O157:H7, \textit{B.~subtilis} 168, \textit{P.~aeruginosa} PAO1, \textit{M.~tuberculosis} H37Rv, \textit{H.~pylori} 26695, \textit{M.~genitalium} G37, \textit{Synechocystis} sp.\ PCC 6803, and \textit{M.~jannaschii} DSM 2661; the last row gives the mean over the eight genomes other than the calibration genome. $F_{\max}$ is computed against all GO annotations and, in parentheses, against experimental annotations. Consistency is $\max(0, 1 - 0.01V)$ for $V$ violations, with the total number of violations in parentheses. Perturbation suppression is the fraction of the 125 injected annotations above threshold after calibration that the constraints remove (fraction of all 900 injected annotations in parentheses). Removal and relaxation give identical values because the nearest valid ancestor of each violating prokaryotic term is the ontology root.}
```

**Figure 3**

```latex
\caption{\textbf{Completion of heteromeric complexes and enforcement of anatomical constraints in eukaryotes.} (\textbf{a}) Leave-subunit-out recovery in 175 trials across six prokaryotic genomes, each ranking the true subunit against five decoys (chance Hit@1 $= 0.167$): Hit@1, recall and $F_1$ at threshold 0.5, and genome complex coherence for the Stage 0 baseline, the complex coherence prior, coherence with operon context, and the full model with Stages 2 and 3. (\textbf{b}) Subunit recovery $F_1$ for each genome; the asterisk marks the calibration genome. In \textit{Synechocystis} sp.\ PCC 6803, many complex subunits are encoded at dispersed loci. (\textbf{c}) Violations remaining in six eukaryotic proteomes under the Stage 0 baseline, direct GO taxon constraints, the Stage 1 anatomical prior, and Stage 2A enforcement. (\textbf{d}) Micro $F_{\max}$ against unmodified Swiss-Prot annotations for the Stage 0 baseline, Stage 2A deletion, and Stage 2A ancestor relaxation.}
```

**Table 2**

```latex
\caption{\textbf{Leave-subunit-out recovery of 175 heteromeric complex subunits in six prokaryotic genomes.} Each trial ranks one subunit with attenuated evidence against five decoys with comparable domain scores from unrelated loci (chance Hit@1 $= 1/6$). Precision, recall, and $F_1$ at threshold 0.5 are computed per trial and averaged over trials. Rows report all six genomes (\textit{E.~coli} K-12, \textit{E.~coli} O157:H7, \textit{B.~subtilis} 168, \textit{M.~tuberculosis} H37Rv, \textit{P.~aeruginosa} PAO1, and \textit{Synechocystis} sp.\ PCC 6803), the five genomes other than the calibration genome \textit{E.~coli} K-12, and \textit{Synechocystis} alone. Complex coherence across all six genomes is the unweighted mean over genomes.}
```

**Table 3** (suggested column labels: *Stage 0 baseline; direct taxon constraints; Stage 1 anatomical prior; Stage 2A deletion; Stage 2A relaxation*)

```latex
\caption{\textbf{Anatomical and taxonomic constraint enforcement in six eukaryotes.} Violation columns give the violations remaining under the Stage 0 baseline (stress test, with unperturbed predictions in parentheses), direct GO taxon constraints, the Stage 1 anatomical prior, and Stage 2A enforcement. The stress test adds 1,546 cross-clade developmental annotations to 7,655 proteins; the unperturbed set contains 6,155 proteins. Micro $F_{\max}$, recall, and mean information content (bits) compare the Stage 0 baseline, Stage 2A deletion, and Stage 2A relaxation against unmodified Swiss-Prot annotations. The developmental rescue column gives the fraction of species-appropriate developmental annotations with attenuated scores that the Stage 1 anatomical prior raises above the decision threshold. The last row reports protein-centric CAFA $F_{\max}$ on the unperturbed predictions.}
```

**Figure 4**

```latex
\caption{\textbf{Recovery of masked metabolic enzymes and candidate pathway repairs in uncultured genomes.} (\textbf{a}) Hit@1, Hit@3, and mean reciprocal rank in 256 leave-reaction-out trials in \textit{E.~coli} K-12 (calibration genome) for sequence-only ranking, homology- and embedding-masked GSPA, RLGC alone, pathway coherence with RLGC, ESM-2 reaction centroids alone, and the full model combining RLGC, coherence, and centroids. (\textbf{b}) Top-1 recovery by metabolic role: core biomass ($n = 307$), alternative catabolic ($n = 49$), secondary and peripheral ($n = 58$), and core biomass reactions whose enzyme shares at most 50\% identity with any annotated homolog ($n = 23$). (\textbf{c}) Hit@1 on 144 reactions in \textit{B.~subtilis} 168 with hyperparameters fixed on \textit{E.~coli}, and on 50 ambiguous reactions before and after Stage 3 re-ranking. (\textbf{d}) Mean number of candidate pathway repairs per genome at posterior thresholds 0.30 and 0.50 in 97 MAGs of uncultured organisms from 16 bacterial and archaeal phyla.}
```

**Table 4**

```latex
\caption{\textbf{Leave-reaction-out recovery of metabolic enzymes in \textit{E.~coli} K-12, \textit{B.~subtilis} 168, and 50 ambiguous reactions.} \textit{E.~coli} K-12 served as the calibration genome. Strict reactions are catalyzed by monomeric or homo-oligomeric enzymes with a single primary EC number; the role-stratified set adds multi-EC and peripheral reactions. ``All'' scores every held-out reaction and counts reactions whose true enzyme is absent from the candidate pool as misses; ``Pool'' is restricted to reactions whose true enzyme is in the candidate pool ($n_{\text{cov}}$). In the homology- and embedding-masked setting, sequence homologs and reaction-centroid donors with cosine similarity above 0.50 are masked; in the embedding setting, Swiss-Prot ESM-2 reaction centroids are retained and sequence alignments to the held-out reaction are masked.}
```

**Figure 5** (three panels; the role-stratified panel is merged into Figure 4b)

```latex
\caption{\textbf{Stoichiometric repair, enzyme identification, and phenotype prediction.} (\textbf{a}) Mean numbers of added reactions, flux-carrying added reactions, zero-flux added reactions, and unbalanced active reactions per genome in 108 bacterial reconstructions under negative-cost MILP (METEOR), a positive-cost linear program, and the lexicographic pFBA of Stage 2B (symmetric logarithmic scale). (\textbf{b}) Enzyme $F_{\max}$ on 148 fitness-validated bacterial enzymes for six predictors alone, with single-predictor GSPA, and with multi-predictor GSPA. (\textbf{c}) Prediction of carbon-source growth (2,332 assays in 529 strains) and diagnostic enzyme activity (17,496 assays in 3,017 type strains) in BacDive under blind prediction, genus-blocked cross-validation, and prediction of 70\% held-out assays from 30\% revealed assays.}
```

**Table 5**

```latex
\caption{\textbf{Enzyme identification on 148 fitness-validated bacterial enzymes from 22 genomes.} For each predictor, columns give the number of correct top-ranked EC numbers (per-protein $F_{\max}$ in parentheses) for the predictor alone, after METEOR's MILP gap-filling, after Stage 2B refinement under glucose and with the local-context reaction probe, after single-predictor GSPA (iterations 1 to 3), after Stage 0 pooling of all predictors, and after multi-predictor GSPA (iteration 3). The last column gives the change in correct top-ranked EC numbers (and in $F_{\max}$) for single-predictor and multi-predictor GSPA relative to the predictor alone. DCGM and ECrecer predict from target sequences and were evaluated on those sequences alone.}
```

**Table 6**

```latex
\caption{\textbf{Prediction of carbon-source growth and diagnostic enzyme activity in BacDive.} Blind methods use zero phenotype labels; supervised classifiers and the hybrid model are evaluated by five-fold cross-validation with each genus confined to one fold (pooled out-of-fold metrics). In the phenotype-anchored protocol, 30\% of each strain's assays are revealed and the remaining 70\% are evaluated. CarveMe, ModelSEED, and gapseq produce binary calls, so their ROC-AUC equals balanced accuracy; GSPA and the supervised classifiers produce continuous scores. The middle rows compare DeepProZyme with METEOR gap-filling on 10 diagnostic enzymes, on urease, and on catabolic assays in 108 genomes. Runtime is per genome after base reconstruction (end-to-end runtime in parentheses).}
```

**Figure 6**

```latex
\caption{\textbf{Evaluation of the Stage 3 guarded agent.} (\textbf{a}) Precision, recall, $F_1$, and resolution accuracy on 120 anomalies (30 pathway gaps, 30 incomplete complexes, 30 borderline essential functions, and 30 borderline ontology conflicts) for deterministic Stages 1 and 2, a lowered posterior threshold with the guardrail, supervised classifiers, a single-turn agent, the multi-turn agent without the guardrail, and the multi-turn agent with the guardrail (mean $\pm$ s.d.\ over three seeds). (\textbf{b}) Rate of accepted taxonomic and anatomical violations. (\textbf{c}) Rank of the true enzyme for 50 ambiguous metabolic gaps before and after Stage 3 re-ranking. (\textbf{d}) Wall-clock time per anomaly across seven species with eight concurrent requests.}
```

**Table 7** (suggested row labels: *Deterministic Stages 1–2; Threshold 0.35 + guardrail; Logistic regression / GBDT (3 pre-tool scores); Single-turn agent; Multi-turn agent, guardrail removed; Multi-turn agent + guardrail (seed 42; 3-seed mean); LR / RF / GBDT / MLP (8 tool features + guardrail)*)

```latex
\caption{\textbf{Comparison of the Stage 3 agent with deterministic, threshold, and supervised baselines on 120 anomalies and 50 ambiguous metabolic gaps.} Gene symbols are masked and candidate order is shuffled. The anomalies comprise 54 promotable gaps, 36 decoys with matched scores, and 30 ontology conflicts (15 genuine conflicts and 15 adversarial traps). Promotion and relaxation counts give true-positive promotions and valid ancestral GO relaxations. Supervised classifiers are trained by five-fold cross-validation on three pre-tool scores or on eight tool-derived features combined with the guardrail. Wall-clock time is per anomaly with eight concurrent requests (single-request latency in parentheses).}
```

---

### Author note: what needs your decision

- **Supplementary renumbering (item 7):** if you accept it, the Methods section and `supplement.tex` need the same renumbering.
- **Figure 5:** it needs to be regenerated with three panels, and all figures and tables should replace internal row labels with the suggested names.
- **Missing citations:** add `.bib` entries for antiSMASH and Qwen3.8, or drop those names.
- **GOLabeler/NetGO:** please confirm you are comfortable with "performed strongly in CAFA" (Intro P3).
