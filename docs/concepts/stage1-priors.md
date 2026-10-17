# Genome-scale priors (Stage 1)

Stage 1 refines the Stage 0 probabilities using what the rest of the genome implies about each protein. It adds six prior log-odds shifts, each derived from one organism-level invariant, to the fixed Stage 0 anchor and iterates until the proteome-wide posterior stops changing. Because every prior is recomputed from the current posterior at each iteration, a change in one protein's annotation (for example a newly confident subunit of a complex) can propagate to other proteins in the next iteration.

## Inputs

| Input | Option | Used by |
|---|---|---|
| Evidence claims (Stage 0) | `--claims` | all |
| GO ontology (OWL) | `--go-owl` | essentiality, coherence, consistency, anatomical |
| Essential function profile | `--essential-profile` (`bacteria`, `archaea`, `eukaryote`, `virus`) or `--essential-functions FILE` | essentiality |
| Pathway definitions | `--pathways`, `--modules` | coherence, operon context |
| Obligate complex classification | `--complex-terms` (bundled default) | coherence, operon context, Stage 2A |
| Organism taxon | `--taxon`, or `--infer-taxon` | consistency, anatomical |
| GO taxon constraints and NCBI lineage | `--taxon-constraints`, `--taxonomy` (bundled defaults) | consistency |
| Metabolic gaps | `--gaps` | gap filling, operon context |
| Operons | `--operons` | operon context |
| Genome layout and reaction graph | `--genome-layout`, `--reaction-graph` | reaction-local genomic context |

A prior whose inputs are missing returns zero and has no effect, so the same command works on a bare proteome and on a fully contextualised genome. `--lite` skips the ELK reasoner; this disables process coherence and the expansion of essential terms to their descendants.

## Update rule

Let \( L^{(0)}(p, f) \) be the Stage 0 anchor and \( \mathbf{P}^{(k-1)} \) the proteome-wide posterior after iteration \( k-1 \). At iteration \( k \) each prior \( B_m \) is recomputed from \( \mathbf{P}^{(k-1)} \), the genome, the taxon, and the reference knowledge, and the posterior log-odds are updated with Jacobi damping:

\[
\tilde L^{(k)}(p, f) = \operatorname{clip}\Bigl(L^{(0)}(p, f) + \sum_{m} w_m\, B_m\bigl(p, f \mid \mathbf{P}^{(k-1)}\bigr),\ -12,\ 12\Bigr),
\]

\[
L^{(k)}(p, f) = (1 - \gamma)\, L^{(k-1)}(p, f) + \gamma\, \tilde L^{(k)}(p, f), \qquad P^{(k)}(p, f) = \sigma\bigl(L^{(k)}(p, f)\bigr).
\]

The iteration stops when the mean absolute probability change over the entries that changed falls below \( \epsilon \), or after \( K_{\max} \) iterations. If the mean change increases on two consecutive iterations, the refiner rolls back to the last stable iterate. Implementation: `gspa-core/.../integration/IterativeRefiner.groovy`.

Throughout, an annotation is *active* when its probability exceeds the decision threshold 0.50.

## The six priors

### Essentiality prior \( B_1 \)

Each domain of life has a core of functions every organism must encode. Let \( \mathcal{F}_{\text{ess}} \) be the GO terms of the selected essential profile and \( \mathcal{U} \subseteq \mathcal{F}_{\text{ess}} \) those not covered by any Stage 0 annotation (after propagation to ancestors). Every candidate claim whose term is an uncovered essential term or one of its descendants receives

\[
B_1(p, f) = \alpha_{\text{ess}}\ \mathbb{I}\Bigl[f \in \textstyle\bigcup_{e \in \mathcal{U}} \{t : t \sqsubseteq e\}\Bigr].
\]

The prior pushes weak candidates for missing core machinery upward; it never lowers a score.

### Pathway and complex coherence prior \( B_2 \)

This prior combines three coherence signals.

- *Process coherence.* For GO `has_part` dependencies whose container process is annotated but whose required part is missing, the missing part receives \( w_{\text{proc}} = \max(0, 1 - c_{\text{proc}}) \), where \( c_{\text{proc}} \) is the fraction of triggered dependencies that are satisfied genome-wide.
- *Pathway coherence.* For every pathway \( W \) with at least one member annotated, each missing required step receives \( w_{\text{path}} = \max(0.25,\ 1 - c_W) \), where \( c_W \) is the fraction of steps covered at Stage 0. A term in several pathways takes the maximum.
- *Obligate heteromeric complexes.* For each complex term classified as obligate heteromeric, let \( N_{\text{strong}} \) be the number of proteins with \( P^{(0)} > 0.50 \) and \( N_{\text{cand}} \) the number of further candidates with \( P^{(0)} \ge 0.08 \) or \( P^{(k-1)} \ge 0.08 \). If exactly one subunit is strong and candidates exist, each candidate receives \( w_{\text{cplx}} = w_{\text{part}} \); if two or more subunits are strong, candidates receive \( 0.65\, w_{\text{part}} \); if a single strong subunit has no candidate partner anywhere in the genome, that subunit receives \( -w_{\text{sing}} \).

\[
B_2(p, f) = \alpha_{\text{coh}}\Bigl(\max\bigl(w_{\text{proc}}(f),\ w_{\text{path}}(f)\bigr) + w_{\text{cplx}}(p, f)\Bigr).
\]

The complex term tells GSPA *which* assembly is incomplete; on its own it boosts every candidate with a matching domain equally. Distinguishing the true subunit needs local genomic context (the operon prior below).

### Taxonomic consistency prior \( B_3 \)

When an organism taxon is asserted (`--taxon`) or inferred (`--infer-taxon`), GO terms forbidden in its lineage by direct GO taxon constraints (`only_in_taxon`, `never_in_taxon`) are collected once. In addition, at every iteration the set of active terms is checked with a SAT solver (SAT4J) over taxon constraints and the NCBI lineage, and the terms in any unsatisfiable core are added to the conflict set \( \mathcal{V}^{(k-1)} \). Any claim whose term, or one of whose ancestors, lies in the conflict set receives

\[
B_3(p, f) = -\alpha_{\text{cons}}.
\]

The penalty is soft by default: a strong likelihood can survive it, which is the intended behaviour for horizontal transfer or contamination. Setting `consistency_hard_filter` in the theta file replaces the penalty by \(-1000\).

### Anatomical prior \( B_4 \)

Many organ-level processes (for example *heart development*) carry no direct GO taxon constraint; their restriction lives on the anatomical structure they act on. GSPA's cross-ontology bridge (see [Stage 2A](stage2a-enforcement.md#the-cross-ontology-constraint-bridge)) entails a GO term as anatomically forbidden in a taxon when the term, or a GO ancestor, is linked to an Uberon, Plant Ontology, or Cell Ontology class whose anatomical ancestors carry a taxon constraint incompatible with the lineage. Let \( \mathcal{T}_{\text{forb}} \) be the forbidden terms and \( \mathcal{U}_{\text{act}} \) the anatomical structures supported by at least two proteins with Stage 0 probability \( \ge 0.50 \). Then

\[
B_4(p, f) = \begin{cases}
-\alpha_{\text{anat}} & f \in \mathcal{T}_{\text{forb}},\\
+\alpha_{\text{boost}} & f \text{ is admissible and linked to a structure in } \mathcal{U}_{\text{act}} \text{ (or its parent)},\\
0 & \text{otherwise.}
\end{cases}
\]

The positive branch rescues species-appropriate developmental terms that Stage 0 attenuated. `anatomical_hard_filter` in the theta file makes the penalty \(-1000\).

### Metabolic gap-filling prior \( B_5 \)

Given metabolic gaps (`--gaps`, one JSON record per gap with the missing EC number and/or GO term, typically from gapseq), every candidate claim whose EC number or GO term matches an open gap receives

\[
B_5(p, f) = \begin{cases} 0.7\,\alpha_{\text{gap}} & \text{gap flagged as a low-homology gapseq guess},\\ \alpha_{\text{gap}} & \text{otherwise.}\end{cases}
\]

Low-homology guesses are marked with `"gapseq_guessed": true` in the gaps file.

This prior says which *functions* the network is missing. It does not say which *gene* supplies them; that is the job of the operon prior, reaction-local genomic context, and Stage 2B.

### Co-directional operon prior \( B_6 \)

Within an operon \( O \) (at least two co-directional genes, from `--operons`), let \( k_O(W) \) be the number of member proteins with an active term required by pathway \( W \). A candidate claim \( (p, f) \) of an operon member, with \( f \) required by \( W \), receives

\[
B_6(p, f) = \max_{W \ni f,\ k_O(W) \ge 1} \Bigl(\alpha_{\text{ctx}}\, \frac{k_O(W)}{|O|}\, \alpha_{\text{gap\_ctx}}^{\,\mathbb{I}[f \text{ fills a metabolic gap}]}\Bigr)
\;+\; \alpha_{\text{ctx}}\, w_{\text{op\_cplx}}\ \mathbb{I}\bigl[f \text{ obligate complex} \wedge \exists p' \in O \setminus \{p\} : P^{(k-1)}(p', f) > 0.50\bigr].
\]

The second term is what lets GSPA pick the true subunit of a complex: among several domain-matched candidates, the one encoded next to an active partner subunit is boosted. In the leave-subunit-out benchmark, adding this prior to the coherence prior is what raised top-ranked subunit recovery well above the Stage 0 baseline (see [Benchmarks](../benchmarks/index.md)). Combined with \( B_2 \) and \( B_5 \), it also lets genome-wide gap signals land on the specific candidate genes in co-directional clusters.

## Reaction-local genomic context

Reaction-local genomic context (RLGC) ranks candidate genes for a missing reaction \( R \) by their position relative to genes that catalyse reactions near \( R \) in the reaction graph. It does not need operon calls.

1. Collect the reactions within \( k \) hops of \( R \) in the reaction graph (breadth-first, currency metabolites removed by a degree percentile).
2. *Anchor genes* are proteins \( p' \) with posterior \( \ge \tau_{\text{anc}} \) for the GO term \( g \) of a neighbouring reaction at distance \( d(g) \). Each anchor has weight \( w(p', R) = \max_g \alpha^{d(g)}\, P(p', g) \).
3. A Gaussian kernel density over base-pair distance, truncated at \( 4h \):

\[
D(p, R) = \sum_{p' : |\Delta \text{bp}(p, p')| \le 4h} w(p', R)\, \exp\Bigl(-\frac{\Delta \text{bp}(p, p')^2}{2h^2}\Bigr).
\]

4. Candidates near the density peaks are scored

\[
S(p, R) = \log D(p, R) + \beta\, \mathrm{Div}(p) - \gamma\, \mathrm{Commit}(p) - \delta\, P(p, f_R) + \varepsilon\, \mathrm{Strand}(p) - \zeta \log_{10} \max\bigl(1, \Delta \text{bp}_{\min}(p)\bigr),
\]

where \( \mathrm{Div} \) rewards support from several distinct neighbouring reactions, \( \mathrm{Commit} \) penalises genes already confidently committed to other functions, the self term excludes genes that already carry \( f_R \), \( \mathrm{Strand} \) rewards strand agreement with the anchors, and the last term penalises large intergenic distance.

5. The top candidates are converted to probabilities with a softmax. If the top candidate exceeds 0.5 a single-gene suggestion is emitted; otherwise a disjunctive suggestion over the smallest set of candidates covering 90% of the mass.

RLGC is enabled with `--rlc-suggester` and needs `--gaps`, `--genome-layout`, and `--reaction-graph` (plus `--reaction-ec-aliases` and optionally `--diffusion-mets` for currency metabolites). Its ranked suggestions are written with `--suggestions-out`; `--features-out` writes the per-candidate feature vector. In the leave-reaction-out benchmarks the RLGC score was combined with orthogroup co-occurrence, Pfam domain compatibility, and (in one track) ESM-2 reaction-centroid similarity; that benchmark scoring lives in the benchmark scripts, not in `gspa integrate`. Implementation: `gspa-core/.../integration/suggester/ReactionLocalContextSuggester.groovy`. The [metabolic gaps and dark matter scenario](../scenarios/metabolic-gaps-dark-matter.md) shows a complete run.

## Hyperparameters

All defaults below are the published values, calibrated on *E. coli* K-12 MG1655 and frozen for every other organism. Theta keys are set in the JSON file passed to `--theta` ([schema](../reference/theta.md)).

| Symbol | Role | Default | Control |
|---|---|---|---|
| \( \alpha_{\text{ess}} \) | essentiality boost | 1.50 | theta `alpha_ess` |
| \( \alpha_{\text{coh}} \) | coherence scale | 1.00 | theta `alpha_coh` |
| \( w_{\text{part}} \) | complex partner boost | 0.85 | theta `alpha_complex_boost` |
| \( w_{\text{sing}} \) | isolated-subunit penalty | 0.75 | theta `alpha_complex_penalty` |
| | complex candidate gate | 0.08 | fixed |
| \( \alpha_{\text{cons}} \) | taxonomic conflict penalty | 3.00 | theta `alpha_cons` |
| \( \alpha_{\text{anat}} \) | anatomical penalty | 3.50 | theta `alpha_anat_penalty` |
| \( \alpha_{\text{boost}} \) | active-structure boost | 0.65 | theta `alpha_anat_boost` |
| \( \alpha_{\text{gap}} \) | gap-filling boost | 1.50 | theta `alpha_gap` |
| \( \alpha_{\text{ctx}} \) | operon pathway scale | 1.00 | theta `alpha_ctx` |
| \( \alpha_{\text{gap\_ctx}} \) | operon gap multiplier | 1.50 | theta `alpha_gap_ctx` |
| \( w_{\text{op\_cplx}} \) | co-operonic partner bonus | 0.95 | fixed |
| \( w_m \) | per-prior weight | 1.0 | theta `prior_weights.<name>` |
| \( \gamma \) | damping | 0.50 | theta `damping` |
| \( K_{\max} \) | maximum iterations | 6 | theta `max_iter` |
| \( \epsilon \) | convergence tolerance | 0.005 | theta `epsilon` |
| \( k \) | RLGC reaction-graph radius | 2 | `--rlc-radius-k` |
| \( \alpha \) | RLGC hop decay | 0.5 | `--rlc-alpha` |
| \( h \) | RLGC kernel bandwidth (bp) | 5000 | `--rlc-kernel-bandwidth` |
| \( \tau_{\text{anc}} \) | RLGC anchor threshold | 0.3 | `--rlc-anchor-threshold` |
| | RLGC currency-metabolite percentile | 99.0 | `--rlc-currency-pct` |
| \( \beta, \gamma, \delta, \varepsilon, \zeta \) | RLGC score coefficients | 1.0, 1.0, 2.0, 0.2, 0.1 | fixed |

Priors are selected with `--enable-priors` (comma-separated; default `essentiality,coherence,consistency,anatomical,gap_filling,genomic_context`). A seventh prior, `homology_transfer`, is available for multi-genome runs: it pulls a protein's claim toward the consensus posterior of its orthogroup (needs `--orthogroups` and `--cluster-consensus`; see [Cross-genome evidence](../scenarios/cross-genome.md)).

## Outputs and audit trail

The contribution of every prior that fired on a pair \( (p, f) \) is stored under its name (`essentiality`, `coherence`, `consistency`, `anatomical`, `gap_filling`, `genomic_context`) in the provenance of that pair. The output TSV lists them in the `priors_fired` column (for example `coherence:0.850,genomic_context:0.950`) next to `likelihood_logodds`, `final_logodds`, `posterior_prob`, and `convergence_iter`; the `--provenance` JSON contains the same information as a map.

## Limitations

- Priors act on candidates proposed by Stage 0; they re-rank and re-weight but do not invent new protein-function pairs (RLGC suggestions are emitted separately and are not folded into the posterior by `gspa integrate`).
- The priors measure the same invariants that the composite audit score measures, so gains in that score are partly by construction. The published evaluation therefore pairs it with held-out GO, EC, leave-reaction-out, leave-subunit-out, and phenotype benchmarks ([Quality metrics](quality-metrics.md)).
- Operon and RLGC signals are strongest where pathway and complex genes are clustered. In genomes with dispersed loci (for example the photosystem and phycobilisome genes of *Synechocystis* sp. PCC 6803) their contribution is smaller.
- Hyperparameters were calibrated on one bacterium. They transferred to other bacteria, archaea, viruses, and eukaryotes in the published benchmarks, but organisms far from that regime may benefit from recalibration.
