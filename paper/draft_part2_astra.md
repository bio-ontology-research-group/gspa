## Section 3: Discussion

```latex
\section{Discussion}

% D1
GSPA formulates proteome annotation as a sequence of inferences that begins with individual proteins and extends to the complete genome. Three drivers account for its performance. Calibrated pooling of sequence homology, domain profiles, protein language models, and structural similarity establishes the protein-level probability distribution. Whole-genome architectural and ontological constraints then relate each candidate annotation to operons, chromosomal neighbourhoods, core functions, heteromeric complexes, taxonomy, and anatomy. Finally, stoichiometric modelling and selective guarded reasoning use pathway requirements and structured genomic evidence to resolve loci that remain ambiguous. The evaluation shows that these layers provide complementary information: protein similarity supplies the initial candidates, genome structure relates those candidates to other genes, and systems-level constraints determine which joint assignments form a coherent proteome annotation.

% D2
Stage 0 extends the multi-source ensemble strategy established by GOLabeler \citep{you2018golabeler} and NetGO \citep{you2019netgo,wang2023netgo3}. These methods showed that sequence, domain, motif, and learned representations capture different aspects of protein function and make partially independent errors. GSPA applies the same premise to sequence alignment, profile HMMs, protein language models, and structural similarity, calibrates their scores, and pools evidence after grouping correlated channels. The resulting probability distribution remains the foundation of the framework because conserved proteins often carry sufficient evidence at the sequence or structure level. Stages 1--3 preserve that foundation while updating proteins for which pathway position, complex membership, taxonomic admissibility, or genomic context provides additional information.

% D3
The value of genome-scale inference increases as complete or near-complete chromosomes become routine products of long-read sequencing and modern assembly \citep{amarasinghe2020longread,logsdon2020longread,kolmogorov2020metaflye}. GSPA uses local chromosomal context to identify co-directional operons and measure the density of genes assigned to neighbouring reactions. This local evidence ranks candidate loci within a genomic region. Genome-wide constraints perform a different operation: they assess coverage of core functions and pathways, require compatible subunits for obligate heteromeric complexes, and test GO annotations against taxonomic and anatomical axioms derived from Uberon, PO, and CL. Semantic ancestor relaxation then maps a clade-inappropriate term to its most specific admissible biological ancestor. Local context and organism-wide constraints therefore connect gene order to the functional and anatomical organization of the complete proteome.

% D4
METEOR was our initial effort to use metabolic gap-filling for protein annotation beyond sequence similarity \citep{niu2026meteor}. The present study identifies two properties of its reaction-indicator MILP formulation that affect gene-level inference: negative addition costs favour reaction accumulation, including reactions with zero steady-state flux, and uniform reaction support distributes evidence across competing paralogs. GSPA connects each metabolic gap to reaction-local genomic context and represents catalytic mechanisms as disjunctions of isozymes and conjunctions of required subunits. Positive-cost lexicographic pFBA selects a parsimonious flux-carrying repair, and competitive locus responsibility distributes its support among candidate genes and an explicit unassigned-gene hypothesis. Cognate-substrate probing extends this inference to catabolic functions outside minimal-glucose growth. Improvements on fitness-validated enzymes and experimental growth and enzyme phenotypes show that these operations connect reaction-level requirements to specific gene loci and organism-level physiology.

% D5
Stage 3 addresses the residual cases for which deterministic stages produce a pathway gap, incomplete complex, uncovered essential function, or borderline ontological conflict. Restricting language-model reasoning to the small selected fraction of proteins makes tool-mediated evidence synthesis practical at proteome scale. The language model compares domain architecture, operon neighbourhoods, pathway or complex state, protein evidence, and candidate ancestral terms, while the symbolic guardrail verifies every proposed annotation against GO, taxonomic, and anatomical constraints. Supervised tabular classifiers provide rapid binary promotion and candidate ranking when labelled anomalies from the target distribution are available. Zero-shot guarded reasoning addresses new clades and structured ancestor selection, where the output includes an ontological resolution in addition to a binary decision. These complementary operating regimes allow labelled data, genomic evidence, and formal semantics to determine the appropriate completion method.
```

## Subsection 3.1: Limitations

```latex
\subsection{Limitations}

% L1
The study design establishes that Stage 0 accounts for most protein-centric accuracy on well-annotated reference proteomes, while genome-scale stages concentrate their updates on twilight-zone proteins, pathway gaps, incomplete complexes, and taxonomic conflicts. In leave-reaction-out evaluation, unconditional recovery also depends on whether Stage 0 places the true enzyme in the candidate pool. Isolated-sequence benchmarks use Stage 0, Stage 2A, and Stage 3 because their inputs lack chromosomal coordinates. A next evaluation can expand candidate generation with remote structural fold search, train or select protein language models under a strict temporal split, stratify leave-subunit-out decoys by same-operon paralogy, and link prospective benchmark targets to their source genome assemblies. This design would measure how additional candidate recall and complete genomic context affect the loci currently governed by protein-level evidence.

% L2
Stage 2B establishes locus refinement from steady-state, mass-balanced flux connectivity over mapped universal reactions under specified media, and Stage 2A establishes semantic enforcement for GO terms with explicit Uberon, PO, or CL logical definitions. These scopes define two model extensions. Transcriptomic expression can provide condition-specific reaction priors, and automated transporter substrate models can connect extracellular compounds to the intracellular network. Additional cross-ontology logical definitions can extend anatomical enforcement to more biological processes and clades. Evaluating these extensions across matched transcriptomic, growth, transport, and anatomical annotation datasets would determine their contribution beyond the current stoichiometric and ontological coverage.

% L3
The Stage 3 evaluation establishes zero-shot disambiguation and ancestral GO relaxation under deterministic verification, while supervised classifiers establish a faster option for binary decisions supported by labelled in-distribution anomalies. A prospective deployment study can route each anomaly according to training-set coverage, decision type, and calibrated uncertainty, and then compare accuracy, latency, and curator effort across the two routes. The pathway repairs proposed for uncultured MAGs and Empty Quarter extremophile isolates provide a separate experimental programme: heterologous expression, targeted metabolomics, and mutant fitness assays can test the predicted enzymes and growth-supporting reactions. These experiments would convert ranked computational hypotheses into direct functional evidence.
```

## Section 4: Methods

```latex
\section{Methods}

\subsection{Problem formulation and four-stage architecture}

% 4.1
Let $\mathcal{P}=\{p_1,\ldots,p_N\}$ denote the protein-coding genes of a target genome with NCBI Taxonomy identifier $\tau_{\text{tax}}$. Let $\mathcal{K}=\{(c_i,s_i,e_i,\zeta_i)\}_{i=1}^{N}$ contain the contig $c_i$, start coordinate $s_i$, end coordinate $e_i$, and strand $\zeta_i\in\{+,-\}$ of each gene. Let $\mathcal{T}$ denote the Gene Ontology directed acyclic graph, with $t\sqsubseteq a$ indicating that term $t$ is a subtype or part of ancestor $a$. The metabolic reaction network is $\mathcal{M}$, and $\mathcal{O}$ contains the taxonomic and cross-ontology axioms used to test annotations. GSPA estimates a proteome-wide annotation matrix $\mathbf{Y}\in\{0,1\}^{N\times|\mathcal{T}|}$ and its posterior probabilities conditional on $\mathcal{K}$, $\tau_{\text{tax}}$, $\mathcal{M}$, and $\mathcal{O}$.

Stage 0 calibrates and combines per-protein evidence. Stage 1 updates these probabilities from genome-wide essentiality, pathway, complex, taxonomic, anatomical, metabolic, and operon evidence. Stage 2A applies deterministic ontological enforcement and semantic ancestor relaxation, while Stage 2B uses stoichiometric flux requirements to update candidate enzyme loci. Stage 3 selects residual anomalies for tool-mediated language-model reasoning and verifies proposed changes with the same symbolic constraints used in Stage 2A (Figure~\ref{fig:architecture}). Supplementary Table~S1 lists the benchmark cohorts, evidence channels, active stages, and evaluation outputs.

\subsection{Stage 0: Multi-channel evidence calibration and correlation-grouped Noisy-OR pooling}

% 4.2
Stage 0 accepts claims from sequence alignment, profile HMM domain annotation, InterPro signatures, ESM-2 and ProtT5 embeddings, sequence-based neural predictors, and FoldSeek or ProstT5 structural similarity \citep{buchfink2021sensitive,eddy2011accelerated,jones2014interproscan,lin2023evolutionary,elnaggar2021prottrans,vankempen2024foldseek,heinzinger2024prostt5}. A claim $c$ assigns function $f$ to protein $p$ with raw score $s(c)\in[0,1]$ and evidence channel $e(c)$. Platt calibration maps the raw score to
$\tilde q(c)=\operatorname{clip}\!\left(\sigma_{\log}(a_{e(c)}s(c)+b_{e(c)}),\varepsilon_{\mathrm{cal}},1-\varepsilon_{\mathrm{cal}}\right)$,
where $\varepsilon_{\mathrm{cal}}=10^{-4}$. Multiplication by channel reliability $r_{e(c)}\in(0,1]$ gives $q(c)=r_{e(c)}\tilde q(c)$. Supplementary Section~S1 gives the channel-specific calibration parameters, reliability weights, and calibration procedure.

Claims for a protein--function pair are assigned to correlation groups $\mathcal{G}_{\mathrm{corr}}(p,f)$ representing homology, protein-sequence models, structure, genomic context, genomic language models, inferred context, motifs, localization, region features, viral features, and domain-specific predictors. Within each group $g$, Stage 0 retains the largest calibrated score and an optional second-opinion contribution from the second-largest score. It then combines correlation groups by Noisy-OR:
\begin{equation}
  q_g(p,f)=\operatorname{clip}\!\left(
  \max_{c\in g}\bigl(r_{e(c)}\tilde q(c)\bigr)
  +\beta_{\mathrm{so}}q_{g,(2)}(p,f),
  10^{-3},0.999\right),
  \qquad
  P_0(p,f)=\operatorname{clip}\!\left(
  1-\prod_{g\in\mathcal{G}_{\mathrm{corr}}(p,f)}\bigl(1-q_g(p,f)\bigr),
  10^{-6},1-10^{-6}\right).
  \label{eq:noisyor}
\end{equation}
The default second-opinion coefficient is $\beta_{\mathrm{so}}=0$. Stage 0 converts $P_0(p,f)$ to the anchored log-odds
$L^{(0)}(p,f)=\operatorname{clip}\!\left(\log[P_0(p,f)/(1-P_0(p,f))],-12,12\right)$
and propagates GO annotations upward through \texttt{is\_a} and \texttt{part\_of} relations.

\subsection{Stage 1: Whole-proteome iterative Bayesian log-odds refinement}

% 4.3, update
Let $\mathbf{P}^{(k-1)}$ denote the annotation probabilities after iteration $k-1$. Stage 1 computes six prior log-odds shifts
$B_m(p,f\mid\mathbf{P}^{(k-1)},\mathcal{K},\tau_{\text{tax}},\mathcal{M},\mathcal{O})$
and adds them to the anchored Stage 0 log-odds. A damped Jacobi update then gives
\begin{equation}
  \begin{aligned}
    \widetilde L^{(k)}(p,f)
      &=\operatorname{clip}\!\left(
        L^{(0)}(p,f)+
        \sum_{m=1}^{6}w_m
        B_m\!\left(p,f\mid\mathbf{P}^{(k-1)},\mathcal{K},
        \tau_{\text{tax}},\mathcal{M},\mathcal{O}\right),
        [-12,12]\right),\\
    L^{(k)}(p,f)
      &=(1-\gamma_{\mathrm{iter}})L^{(k-1)}(p,f)
        +\gamma_{\mathrm{iter}}\widetilde L^{(k)}(p,f),
        \qquad
    P^{(k)}(p,f)=\sigma_{\log}\!\left(L^{(k)}(p,f)\right).
  \end{aligned}
  \label{eq:stage1}
\end{equation}
We use $w_m=1$, $\gamma_{\mathrm{iter}}=0.50$, and at most $K_{\max}=6$ iterations. Iteration stops when the mean absolute probability change over changed protein--function pairs falls below $\epsilon=5\times10^{-3}$. An increase in this mean change over two consecutive iterations restores the last stable iterate. Active pinned floors are applied after damping.

The six priors are defined as follows:
\begin{enumerate}
  \item $B_1$ is the essentiality prior. Let $\mathcal{F}_{\mathrm{ess}}(\tau_{\mathrm{tax}})$ be the domain-specific essential GO profile and let $\mathcal{U}_{\mathrm{ess}}(\tau_{\mathrm{tax}})$ contain categories with posterior at most $0.50$ for every protein. Define
  $\mathcal{D}_{\mathrm{uncov}}(\tau_{\mathrm{tax}})
  =\bigcup_{e\in\mathcal{U}_{\mathrm{ess}}(\tau_{\mathrm{tax}})}
  \{t\in\mathcal{T}:t\sqsubseteq e\}$.
  Then
  $B_1(p,f)=\alpha_{\mathrm{ess}}
  \mathbb{I}[f\in\mathcal{D}_{\mathrm{uncov}}(\tau_{\mathrm{tax}})]$,
  with $\alpha_{\mathrm{ess}}=1.5$ and maximum boost multiplier $2.0$.

  \item $B_2$ is the pathway and heteromeric-complex coherence prior. For a process \texttt{has\_part} relation with satisfied fraction $c_{\mathrm{proc}}$, a missing part receives
  $w_{\mathrm{proc}}(f)=\max(0,1-c_{\mathrm{proc}})$.
  For a triggered pathway $W$ with coverage $c_W\in(0,1)$, a missing required term receives
  $w_{\mathrm{path}}(f)=\max_{W\ni f}\max(0.25,1-c_W)$.
  For an obligate heteromeric complex, let $N_{\mathrm{strong}}(f)$ count proteins with $P^{(0)}(p,f)>0.50$ and let $N_{\mathrm{cand}}(f)$ count proteins with
  $\max(P^{(0)}(p,f),P^{(k-1)}(p,f))\geq0.08$.
  A candidate partner receives $w_{\mathrm{cplx}}=0.85$ when
  $N_{\mathrm{strong}}=1$ and $N_{\mathrm{cand}}\geq2$, and $w_{\mathrm{cplx}}=0.5525$ when $N_{\mathrm{strong}}\geq2$.
  An isolated strong singleton with $N_{\mathrm{cand}}\leq1$ receives $w_{\mathrm{cplx}}=-0.75$.
  Therefore,
  $B_2(p,f)=\alpha_{\mathrm{coh}}
  [\max(w_{\mathrm{proc}}(f),w_{\mathrm{path}}(f))
  +w_{\mathrm{cplx}}(p,f)]$, with $\alpha_{\mathrm{coh}}=1$.

  \item $B_3$ is the taxonomic consistency prior. Let
  $\mathcal{A}^{(k-1)}=\{t\in\mathcal{T}:\max_pP^{(k-1)}(p,t)>0.50\}$.
  Direct GO taxon rules and the unsatisfiable core of the joint constraint model define
  $\mathcal{V}_{\mathrm{conf}}^{(k-1)}
  =\mathcal{T}_{\mathrm{forb}}^{\mathrm{TC}}(\tau_{\mathrm{tax}})
  \cup\operatorname{Core}_{\mathrm{SAT}}(\mathcal{A}^{(k-1)})$.
  A claim whose term or ancestor lies in this set receives
  $B_3(p,f)=-\alpha_{\mathrm{cons}}$, with
  $\alpha_{\mathrm{cons}}=3.0$.

  \item $B_4$ is the anatomical prior. Let
  $\mathcal{T}_{\mathrm{forb}}^{\mathrm{anat}}(\tau_{\mathrm{tax}})$
  contain GO terms whose anatomical definitions conflict with the target taxon.
  Let $\mathcal{U}_{\mathrm{act}}(\tau_{\mathrm{tax}})$ contain valid anatomical structures supported by at least two proteins with
  $P^{(0)}(p,f)\geq0.50$. Then
  \[
  B_4(p,f)=
  \begin{cases}
    -3.5, &
    (\{f\}\cup\operatorname{Anc}(f))
    \cap\mathcal{T}_{\mathrm{forb}}^{\mathrm{anat}}(\tau_{\mathrm{tax}})
    \neq\emptyset,\\
    +0.65, &
    \operatorname{Anat}(f)\cap
    \mathcal{U}_{\mathrm{act}}(\tau_{\mathrm{tax}})
    \neq\emptyset,\\
    0, & \text{otherwise.}
  \end{cases}
  \]

  \item $B_5$ is the metabolic gap and reaction-local genomic-context prior. A direct match between candidate function $f$ and an open metabolic gap receives $B_5(p,f)=1.5$; a low-homology proposed gap receives $B_5(p,f)=1.05$. For a target reaction $R$, reaction-graph neighbours within two hops provide anchor weights
  $w_{\mathrm{anc}}(p',R)=
  \max_g\alpha_{\mathrm{rxn}}^{d(g)}P(p',g)$
  over neighbours with $P(p',g)\geq\tau_{\mathrm{anc}}=0.30$, where
  $\alpha_{\mathrm{rxn}}=0.50$. With Gaussian bandwidth
  $h_{\mathrm{bp}}=5{,}000$ bp, the reaction-local density and candidate score are
  \begin{equation}
    \begin{aligned}
      D_{\mathrm{RLGC}}(p,R)
        &=\sum_{p':|\Delta\mathrm{bp}(p,p')|\leq4h_{\mathrm{bp}}}
          w_{\mathrm{anc}}(p',R)
          \exp\!\left(
          -\frac{\Delta\mathrm{bp}(p,p')^2}{2h_{\mathrm{bp}}^2}\right),\\
      S_{\mathrm{RLGC}}(p,R)
        &=\log D_{\mathrm{RLGC}}(p,R)
          +\beta_{\mathrm{div}}\operatorname{Div}(p)
          -\gamma_{\mathrm{com}}\operatorname{Commit}(p)
          -\delta_{\mathrm{self}}P(p,f_R)
          +\varepsilon_{\mathrm{str}}\operatorname{Strand}(p)
          -\zeta_{\mathrm{gap}}
          \log_{10}\!\left(\max(1,\Delta\mathrm{bp}_{\min}(p))\right).
    \end{aligned}
    \label{eq:rlgc}
  \end{equation}
  We use $\beta_{\mathrm{div}}=\gamma_{\mathrm{com}}=1$,
  $\delta_{\mathrm{self}}=2$, $\varepsilon_{\mathrm{str}}=0.2$,
  and $\zeta_{\mathrm{gap}}=0.1$. The method evaluates at most 50 candidates with
  $D_{\mathrm{RLGC}}\geq0.01$ and normalizes the leading 20 scores by softmax.
  Cross-genome orthogroup co-occurrence, Pfam compatibility, and ESM-2 reaction-centroid similarity provide additional ranking evidence where specified.

  \item $B_6$ is the co-directional operon prior. For an operon $O$ with
  $|O|\geq2$ and pathway $W$, let
  $k_O(W)=|\{p'\in O:\max_{f'\in\operatorname{Req}(W)}
  P^{(k-1)}(p',f')>0.50\}|$ and $c_O(W)=k_O(W)/|O|$.
  For $p\in O$ and $f\in\operatorname{Req}(W)$,
  \[
  B_6(p,f)=
  \max_{W\ni f:k_O(W)\geq1}
  \left(\alpha_{\mathrm{ctx}}c_O(W)
  \alpha_{\mathrm{gap\_ctx}}^{\mathbb{I}[f\in\mathcal{G}_{\mathrm{gap}}]}\right)
  +w_{\mathrm{op\_cplx}}
  \mathbb{I}\!\left[
  f\in\mathcal{C}_{\mathrm{het}}\wedge
  \exists p'\in O\setminus\{p\}:P^{(k-1)}(p',f)>0.50
  \right],
  \]
  with $\alpha_{\mathrm{ctx}}=1$,
  $\alpha_{\mathrm{gap\_ctx}}=1.5$, and
  $w_{\mathrm{op\_cplx}}=0.95$.
\end{enumerate}
All configurable Stage 1 and Stage 2B parameters were calibrated on the \textit{E.~coli} K-12 MG1655 development proteome and frozen for the remaining benchmarks. Supplementary Table~S2 gives the full hyperparameter specification.

\subsection{Cross-ontology taxonomic and anatomical constraint bridge}

% 4.4
We constructed five linked resources from the official releases of GO, Uberon \citep{mungall2012uberon}, Plant Ontology \citep{cooper2013plant}, Cell Ontology \citep{diehl2016cell}, and NCBI Taxonomy. The first contains 26,911 direct GO taxon constraints across 18,159 GO terms, comprising 13,819 \texttt{only\_in\_taxon} and 13,092 \texttt{never\_in\_taxon} rules. The second contains 2,715 equivalence or subclass axioms linking 2,674 GO terms to 1,294 anatomical classes through 30 relations, including development, morphogenesis, formation, maturation, localization, transport, and organization relations. The third contains 8,331 taxon constraints over 1,247 Uberon, PO, and CL classes. The fourth contains 7,722 hierarchy edges among 2,855 anatomical classes, and the fifth contains 1,818 hierarchy edges among 1,305 taxonomic clades. These resources allow a GO term to inherit taxonomic restrictions through both its GO ancestors and the anatomical classes in its logical definition. Supplementary Table~S3 reports resource versions, counts, relations, and construction rules.

\subsection{Stage 2A: Deterministic symbolic enforcement and semantic ancestor relaxation}

% 4.5
Stage 2A identifies forbidden GO terms from direct taxon constraints, taxonomic restrictions entailed through Uberon, PO, or CL definitions, and unsatisfiable cores in the joint constraint model. The operation includes the upward GO closure of every candidate term. For a forbidden term $t$, let
$\mathcal{A}_{\mathrm{valid}}(t,\tau_{\mathrm{tax}})$
be the set of most specific taxonomically and anatomically valid ancestors of $t$, excluding the three GO aspect roots. An ancestor belongs to this set when it is valid in $\tau_{\mathrm{tax}}$ and has no valid descendant on the path from $t$.

Before applying consistency penalties, Stage 2A reconstructs
$P_{\mathrm{unpen}}(p,t)$ from $L^{(0)}(p,t)$ and all positive or unrelated Stage 1 prior contributions. It transfers this probability to each minimal valid ancestor and removes the forbidden term:
\begin{align}
  P_{\mathrm{Stage2A}}(p,a)
    &=\max\!\left(
      P^{(K)}(p,a),
      \gamma_{\mathrm{relax}}P_{\mathrm{unpen}}(p,t)
      \right),
      \qquad
      a\in\mathcal{A}_{\mathrm{valid}}(t,\tau_{\mathrm{tax}}),
      \nonumber\\
  P_{\mathrm{Stage2A}}(p,t)&=0.
  \label{eq:relax}
\end{align}
We use $\gamma_{\mathrm{relax}}=1$. A forbidden term with an empty valid-ancestor set is removed. Stage 2A also examines each obligate heteromeric complex. When exactly one protein reaches $0.50$ and the proteome contains zero partner candidates reaching $0.08$, the singleton posterior is set to $0.45$. Supplementary Section~S2 describes constraint reasoning, ancestor selection, and enforcement reports.

\subsection{Stage 2B: Role-aware reaction support and positive-cost lexicographic pFBA}

% 4.6, aggregation
Stage 2B uses a sanitized universal network
$\mathcal{U}_{\mathrm{clean}}$ derived from ModelSEED. The source network contains 47,880 reactions and 33,978 metabolites; 14,243 reactions have EC mappings and 11,826 enzymatic reactions are mass and charge balanced. Pool and imbalanced reactions receive penalty $w_{\mathrm{imb}}(j)=100$. Basal inorganic and cofactor exchanges define the remaining environmental bounds.

A reaction $j$ can have several alternative catalytic mechanisms
$\mathcal{M}_{\mathrm{dnf}}(j)$. Each mechanism $m$ requires a set of roles or subunits
$\mathcal{R}_{\mathrm{sub}}(m)$, and each role $u$ has candidate locus scores
$s_{p,u}\in(0,1)$. Stage 2B combines the three leading isozyme candidates for a role by Noisy-OR, combines required subunits by a gated geometric mean, and selects the strongest alternative mechanism:
\begin{equation}
  \begin{aligned}
    Q(u)
      &=1-\prod_{p\in\operatorname{Top-}M_{\mathrm{iso}}
        (\mathcal{C}_{\mathrm{cand}}(u))}
        (1-s_{p,u}),\\
    S_{\mathrm{mech}}(m)
      &=
      \begin{cases}
        \left(
        \prod_{u\in\mathcal{R}_{\mathrm{sub}}(m)}Q(u)
        \right)^{1/D_m},
        & \min_u Q(u)\geq\theta_{\mathrm{sub}},\\
        \min_{u\in\mathcal{R}_{\mathrm{sub}}(m)}Q(u),
        & \text{otherwise,}
      \end{cases}\\
    S_{\mathrm{rxn}}(j)
      &=\max_{m\in\mathcal{M}_{\mathrm{dnf}}(j)}
        S_{\mathrm{mech}}(m).
  \end{aligned}
  \label{eq:gpr_agg}
\end{equation}
Here $M_{\mathrm{iso}}=3$, $D_m=|\mathcal{R}_{\mathrm{sub}}(m)|$, and
$\theta_{\mathrm{sub}}=0.25$.

Reactions with $S_{\mathrm{rxn}}(j)\geq\tau_{\mathrm{core}}=0.50$ form the draft core $\mathcal{R}_{\mathrm{core}}$, together with open medium exchanges and spontaneous reactions. Candidate repairs belong to the active solver domain when
$S_{\mathrm{rxn}}(j)\geq10^{-3}$ or the reaction belongs to the universal flux-support skeleton. Their positive cost is
$c_j^{\mathrm{add}}=\lambda_{\mathrm{add}}(1-S_{\mathrm{rxn}}(j))
+w_{\mathrm{imb}}(j)+\epsilon_{\mathrm{par}}$,
where $\lambda_{\mathrm{add}}=1$ and $\epsilon_{\mathrm{par}}=10^{-3}$.

For medium $\mathcal{E}_{\mathrm{med}}$, the biomass target is
$v_{\mathrm{bio}}^{\mathrm{target}}
=\max(\alpha_{\mathrm{bio}}v_{\mathrm{bio}}^{\max},0.1)$
with $\alpha_{\mathrm{bio}}=0.10$. The first linear programme minimizes candidate-reaction cost subject to
$\mathbf{S}\mathbf{v}=\mathbf{0}$ and
$v_{\mathrm{bio}}\geq v_{\mathrm{bio}}^{\mathrm{target}}$.
The second minimizes total non-biomass flux while constraining repair cost to
$(1+\delta_{\mathrm{opt}})C^*+10^{-9}$, where
$\delta_{\mathrm{opt}}=10^{-6}$ and $C^*$ is the optimum from the first programme. We solved both programmes with HiGHS at feasibility tolerance $10^{-9}$. The repair set
$\mathcal{J}_{\mathrm{rep}}(\mathcal{E}_{\mathrm{med}})
=\{j\notin\mathcal{R}_{\mathrm{core}}:
|v_j^{\mathrm{pFBA}}|>10^{-6}\}$
therefore contains flux-carrying reactions.

For each repair, flux-variability analysis, reaction knockout, and five cost perturbations provide
$\pi_j^{\mathrm{FBA}}
=0.5\mathbb{I}[j\in\mathcal{R}_{\mathrm{FVA}}]
+0.3\mathbb{I}[j\in\mathcal{R}_{\mathrm{KO}}]
+0.2\omega_{\mathrm{ens}}(j)$.
The flux-variability certificate uses a cost bound of $1.01C^*+10^{-7}$; the knockout certificate uses $1.10C^*+10^{-6}$. Cost perturbations use independent uniform changes of $\pm5\%$, with
$\omega_{\mathrm{ens}}=1$ for activity in at least $80\%$ of solves,
$0.65$ for lower positive frequency, and $0$ for zero frequency.
Supplementary Section~S4 gives network sanitization, media, solver bounds, and certificate procedures, and Supplementary Table~S6 reports the stoichiometric audit.

\subsection{Stage 2B: Competitive gene-locus responsibility and substrate probing}

% 4.7, responsibility
For each reaction $j$ and required role $u$, Stage 2B distributes reaction support among candidate loci and an unassigned-gene hypothesis $h_\varnothing$ with prior mass $\eta_{\mathrm{null}}=0.05$. Let
$o_{p,u}=s_{p,u}/(1-s_{p,u})$ be the sequence odds,
$\operatorname{BF}_{\mathrm{operon}}(p,u)$ the operon Bayes factor,
$C_{\mathrm{syn}}(p,u)$ the pathway-synteny score, and
$S_{\mathrm{ESM2}}(p,u)$ the protein-language-model similarity. Candidate weights and normalized responsibilities are
\begin{align}
  w_{j,u}(p)
    &=o_{p,u}^{\beta_1}
      \operatorname{BF}_{\mathrm{operon}}(p,u)^{\beta_2}
      \exp\!\left(
      \beta_3 C_{\mathrm{syn}}(p,u)
      +\beta_4 S_{\mathrm{ESM2}}(p,u)\right),
      \nonumber\\
  r_{j,u}(p)
    &=\frac{w_{j,u}(p)}
      {\eta_{\mathrm{null}}+
      \sum_{p'\in\mathcal{C}_{\mathrm{cand}}(u)}
      w_{j,u}(p')}.
  \label{eq:responsibility}
\end{align}
We use $\beta_1=\beta_2=1$ and $\beta_3=\beta_4=1.5$.

The stoichiometric contribution to a locus is
$\Delta L_{\mathrm{FBA}}(p,u)
=\gamma_{\mathrm{FBA}}
\sum_{j\in\mathcal{J}_{\mathrm{rep}}(u)}
\pi_j^{\mathrm{FBA}}r_{j,u}(p)(1-S_{\mathrm{rxn}}(j))$,
with $\gamma_{\mathrm{FBA}}=2.2$. Iterative refinement applies step weight
$\lambda_t=\lambda_0\gamma_{\mathrm{damp}}^{t-1}$,
where $\lambda_0=0.70$ and $\gamma_{\mathrm{damp}}=0.60$.
For catabolic enzymes, the method repeats flux probing under the cognate substrate and combines locus responsibility with one-hop and two-hop non-currency-metabolite neighbours and co-directional operon support. Supplementary Section~S4 specifies the three-pass outer update, local gene window, candidate gates, and `Price-149` parameterization.

For each held-out reaction $j$ and true gene $p^*$, five indicators describe structural recoverability. $\mathcal{C}_{\mathrm{struct}}(j)=1$ when $j$ is present and mass balanced in $\mathcal{U}_{\mathrm{clean}}$. $\mathcal{C}_{\mathrm{task}}(j,\mathcal{E}_{\mathrm{med}})=1$ when blocking $j$ requires positive repair flux under medium $\mathcal{E}_{\mathrm{med}}$. $\mathcal{C}_{\mathrm{unique}}(j)=1$ when an unmasked isozyme cannot carry the complete flux. $\mathcal{C}_{\mathrm{operon}}(p^*,j)=\mathbb{I}[\operatorname{BF}_{\mathrm{operon}}(p^*,j)>1]$. Finally, $\mathcal{C}_{\mathrm{complex}}(p^*,j)=1$ when a multisubunit catalyst has partner subunits present. These indicators define the core-growth, cognate-substrate, and secondary-metabolism evaluation strata.

\subsection{Stage 3: Selective neuro-symbolic annotation completion}

% 4.8
Stage 3 selects four types of post-Stage-2 anomaly. A pathway gap has pathway completeness at least $0.50$ and candidate posteriors in $[0.05,0.50)$; the controlled agent benchmark uses pathways with completeness at least $0.65$ and the five leading candidates in $[0.15,0.49]$. An incomplete heteromeric complex has one active anchor at posterior at least $0.50$ and at least one sub-threshold partner candidate; benchmark complexes have partner coverage of at least $50\%$. A borderline essentiality case contains an uncovered member of $\mathcal{F}_{\mathrm{ess}}(\tau_{\mathrm{tax}})$ with sub-threshold candidates. A borderline conflict contains a posterior in $[0.35,0.62]$ whose taxonomic or anatomical restriction permits refinement to an ancestral term. Genome-wide extraction retains at most eight candidates per anomaly and 150 anomalies per genome.

Candidate descriptions exclude evaluation labels and gene symbols and present every candidate with the same evidence fields and deterministic order. A locally hosted FP8-quantized 27-billion-parameter instruction model from the Qwen open-weight family \citep{yang2024qwen25} receives a 32,768-token context, temperature $0$, at most 450 output tokens per turn, and at most six turns. We evaluated seeds 42, 43, and 44. The model can inspect an operon neighbourhood, inspect a pathway or complex gap, query protein evidence, check symbolic constraints, and propose a resolution. These operations are read-only until the final proposal.

The symbolic guardrail accepts a proposal after checking the complete GO term set, 1,305 NCBI Taxonomy clades, 26,911 direct GO taxon constraints, and 2,715 GO--Uberon/PO/CL anatomical constraints. Promotion requires a valid term for the target taxon; relaxation requires a valid ancestor selected by the Stage 2A criterion. Supplementary Section~S7 gives the prompts, operation schemas, de-identification procedure, and acceptance rules, while Supplementary Table~S8 gives the multi-seed evaluation.

\subsection{Experimental phenotype prediction protocols in BacDive}

% 4.9, cohorts and coupling
We evaluated 2,332 carbon-source growth assays from 529 ProTraits-matched bacterial genomes and 17,496 diagnostic enzyme assays from 3,017 BacDive type-strain assemblies. For strain $s$ and assay $y$, let
$g(s,y)$, $m(s,y)$, and $c(s,y)$ be binary calls from gapseq v1.2 \citep{zimmermann2021gapseq}, ModelSEED v2.6.1 \citep{seaver2021modelseed}, and CarveMe v1.5.2 \citep{machado2018fast}. For carbon assays,
$B(s,y)=0.55g(s,y)+0.30m(s,y)+0.15c(s,y)$;
for enzyme assays,
$B(s,y)=0.55g(s,y)+0.28m(s,y)+0.17c(s,y)$.
After column centring,
$\widetilde B_{s,y}=B(s,y)-\overline B_y$,
we define
$\Psi_0(y,y')=\max(0,\widetilde B_{:,y}^{\top}\widetilde B_{:,y'}/
[\|\widetilde B_{:,y}\|_2\|\widetilde B_{:,y'}\|_2+10^{-9}])$
and set $\Psi_0(y,y)=0$. Row normalization gives
$\Psi(y,y')=\Psi_0(y,y')/
[\sum_{y''}\Psi_0(y,y'')+10^{-9}]$.
The unsupervised cross-assay coupling score is
$\psi(s,y)=\sum_{y'}\Psi(y,y')[B(s,y')-0.5]$.

In blind Track 1, each binary call is mapped to a calibrated tool probability. Carbon assays use
$p_{\mathrm{gap}}\in\{0.12,0.88\}$,
$p_{\mathrm{ms}}\in\{0.22,0.80\}$, and
$p_{\mathrm{cm}}\in\{0.30,0.72\}$.
Enzyme assays use
$p_{\mathrm{gap}}\in\{0.14,0.86\}$,
$p_{\mathrm{ms}}\in\{0.24,0.78\}$, and
$p_{\mathrm{cm}}\in\{0.28,0.74\}$.
The blind posterior is
\begin{equation}
  P_{\mathrm{Track1}}(s,y)
    =\sigma_{\log}\!\left(
      w_{\mathrm{gap}}\operatorname{logit}(p_{\mathrm{gap}}(s,y))
      +w_{\mathrm{ms}}\operatorname{logit}(p_{\mathrm{ms}}(s,y))
      +w_{\mathrm{cm}}\operatorname{logit}(p_{\mathrm{cm}}(s,y))
      +w_\psi\psi(s,y)
      \right).
  \label{eq:bacdive_t1}
\end{equation}
The weights are $(0.56,0.28,0.16,0.65)$ for carbon and
$(0.55,0.27,0.18,0.55)$ for enzymes, with decision threshold $0.50$.

We evaluated Traitar \citep{weimann2016traitar}, a random-forest baseline, DiASPora \citep{koblitz2025diaspora,schober2025bac}, and a GSPA--DiASPora hybrid by five-fold genus-blocked cross-validation. Every genus belongs to one fold, and calibration, prevalence estimates, threshold selection, and classifier fitting use the training folds. Baseline features contain the three metabolic-model calls and an assay indicator. GSPA features additionally contain training-fold calibrated model probabilities, assay prevalence, mean strain-level pathway support, cross-assay coupling, the blind Track 1 posterior, a pooled posterior, the complete strain call vector, and the assay indicator. The hybrid combines the out-of-fold classifier probability with the pooled GSPA posterior using weights $0.88$ and $0.12$.

Track 2 reveals $30\%$ of each strain's assay rows and evaluates the remaining $70\%$. For carbon assays, the anchored posterior combines the pooled score with revealed substrate prevalence and revealed-assay co-occurrence:
$P_{\mathrm{Track2}}(s,y)=
\sigma_{\log}[0.72\operatorname{logit}(p_{\mathrm{pool}}(s,y))
+0.28\operatorname{logit}(s_{b,\mathrm{emp}}(y))
+0.85\psi_{\mathrm{rev}}(s,y)]$.
For enzyme assays, it combines the pooled score with the revealed genus--EC posterior:
$P_{\mathrm{Track2}}(s,y)=
\sigma_{\log}[0.70\operatorname{logit}(p_{\mathrm{pool}}(s,y))
+0.30\operatorname{logit}(g_{e,\mathrm{emp}}(s,y))]$.
Thresholds are selected on the revealed partition. The Track 2 hybrid uses
$0.65P_{\mathrm{ML2}}+0.35P_{\mathrm{Track2}}$.
Supplementary Section~S6 gives phenotype mapping and calibration details, and Supplementary Table~S7 gives per-fold results.

\subsection{Benchmarks, homology masking, statistical testing, and structural audit score}

% 4.10
The benchmark collection comprises temporal GO and EC holdouts, homology-masked bacterial enzyme panels, reference-genome ablations, multicellular anatomical perturbations, leave-reaction-out and leave-subunit-out experiments, fitness-validated bacterial enzymes, BacDive growth and enzyme phenotypes, uncultured MAGs, controlled Stage 3 anomalies, and prospective sequence targets. Hyperparameters were selected on \textit{E.~coli} K-12 MG1655 and frozen for out-of-sample genomes. Supplementary Table~S1 gives cohort sizes, inclusion criteria, active stages, and primary comparisons; Supplementary Sections~S3--S8 give benchmark-specific construction and provenance.

Donor-homology masking excludes reference proteins above sequence-identity cutoffs
$\eta_{\mathrm{id}}\in\{100\%,70\%,50\%,30\%\}$ or within the designated genus or family of a test protein. The same exclusions are applied to alignment databases, protein-language-model neighbour indices, structural templates, and probability calibration sets. Neural enzyme predictors use their published homology-split procedures with frozen released backbone weights. Leave-reaction-out Track A excludes reaction centroids or candidates with cosine similarity above $0.50$, while Track B uses frozen ESM-2 reaction centroids. These two tracks separate recovery from local genome context from recovery supported by pretrained sequence representations.

We evaluated GO predictions with protein-centric $F_{\max}$, information-accretion-weighted $WF_{\max}$, pair-level micro $F_{\max}$, and minimum semantic distance $S_{\min}$ \citep{radivojac2013large,jiang2016expanded,clark2013information}. Thresholds range from $0.01$ to $0.99$ in increments of $0.01$; the anatomical perturbation study uses $0.10$ to $0.90$ in increments of $0.10$. Ranking experiments report Hit@$k$ and mean reciprocal rank, and binary phenotype experiments report precision, recall, $F_1$, balanced accuracy, and ROC-AUC as applicable. We calculated $95\%$ bootstrap or exact Clopper--Pearson confidence intervals and used two-sided paired Wilcoxon or exact McNemar tests for paired comparisons. Supplementary Section~S5 specifies sampling units, confidence intervals, and multiplicity handling.

The Genome-Scale Annotation Evaluation Framework summarizes structural properties at $\theta_{\mathrm{eval}}=0.50$. Essential-function completeness is $C_{\mathrm{ess}}$. Functional coherence $C_{\mathrm{func}}$ is the mean of the available pathway, complex, and process coherence scores. Taxonomic and anatomical consistency is
$C_{\mathrm{cons}}=\max(0,1-0.01V(\mathcal{P}))$, where
$V(\mathcal{P})$ counts formal violations. Functional resolution is
$C_{\mathrm{res}}=\min(1,\operatorname{meanIC}/10)$.
With default weights
$(w_{\mathrm{ess}},w_{\mathrm{func}},w_{\mathrm{cons}},w_{\mathrm{res}})
=(0.30,0.40,0.30,0)$,
the score is
\begin{equation}
  \operatorname{GAEF}(\mathcal{P})
    =\frac{
      \sum_{d\in\{\mathrm{ess},\mathrm{func},\mathrm{cons},\mathrm{res}\}:
      C_d\neq\mathrm{null}}
      w_dC_d
    }{
      \sum_{d:C_d\neq\mathrm{null}}w_d
    }.
  \label{eq:gaef}
\end{equation}
When completeness, coherence, and consistency are available, this expression becomes
$0.30C_{\mathrm{ess}}+0.40C_{\mathrm{func}}+0.30C_{\mathrm{cons}}$.
GAEF measures the structural invariants targeted by Stages 1 and 2 and is reported alongside held-out functional and phenotype evaluations.
```

## Data availability and code availability

```latex
\section*{Data availability}

The ontology releases, cross-ontology constraints, benchmark manifests, frozen evaluation results, and machine-readable provenance are available in the GSPA repository at \url{https://github.com/bio-ontology-research-group/gspa}, release v2.1.0. The repository includes a versioned checksum manifest for the benchmark and derived data artifacts. Source ontology licenses and attribution are provided with each redistributed resource.

\section*{Code availability}

GSPA v2.1.0, the executable release, benchmark workflows, stoichiometric refinement implementation, neuro-symbolic agent configuration, and figure-generation workflow are available under the GPL-3.0-or-later license at \url{https://github.com/bio-ontology-research-group/gspa}. The tagged release contains the environments and instructions required to reproduce the analyses and figures.
```

I applied the `write-edit-scientific-paper` skill to preserve the derivational premises, organize Methods by formal objects and operations, and keep limitations as positive statements paired with specific next experiments.