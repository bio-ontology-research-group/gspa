# Genome-scale quality metrics

GSPA reports two kinds of measurements. *Genome-scale quality metrics* describe properties of a whole annotation set that need no ground truth: is the essential core present, are pathways and complexes coherent, and is the set consistent with the organism's taxonomy and anatomy? *Protein-level accuracy metrics* compare predictions with reference annotations and are used in the benchmarks. This page defines both and explains how they relate.

The genome-scale metrics are computed by [`gspa evaluate`](../reference/cli-evaluate.md) on any GFF3/GAF annotation (from GSPA or another pipeline), by `gspa annotate` at the end of a run unless `--no-metrics` is given, and compared across annotation sets by [`gspa compare`](../reference/cli-compare.md). The [evaluate, compare, report scenario](../scenarios/evaluate-compare-report.md) walks through a comparison. Implementation: `gspa-core/.../metrics/` and `gspa-core/.../model/QualityReport.groovy`.

All genome-scale metrics are computed on the annotations active at the decision threshold (\( P \ge 0.50 \)), after propagating each annotation to its GO ancestors along `is_a` and `part_of`.

## Completeness

Let \( \mathcal{F}_{\text{ess}} \) be the essential GO terms of the profile chosen with `--kingdom` (bundled profiles for bacteria, archaea, eukaryotes, and viruses). An essential term is *present* if any protein carries the term or one of its descendants. Completeness is

\[
C_{\text{ess}} = \frac{|\{ e \in \mathcal{F}_{\text{ess}} : e \text{ present} \}|}{|\mathcal{F}_{\text{ess}}|},
\]

reported overall and per functional category, together with the lists of present and missing essential terms.

## Coherence

Coherence measures whether functions that depend on each other are annotated together. It has three components, each in \([0, 1]\):

- **Process coherence** \( C_{\text{proc}} \): of the GO `has_part` dependencies triggered by the annotation set (a process is annotated whose definition requires a part), the fraction whose required part is also annotated. Computing it requires the ELK reasoner; with `--lite` it is not evaluated.
- **Pathway coherence** \( C_{\text{path}} \): the mean completeness of the pathways that are triggered (at least one step annotated), using the pathway definitions given with `--pathways` and `--modules`.
- **Complex coherence** \( C_{\text{cplx}} \): of the complex terms that are annotated, the fraction that are coherent. An obligate heteromeric complex is coherent only if at least two distinct proteins carry it; homomeric and ambiguous complexes are coherent when present. The classification comes from the bundled `complex_terms.tsv` or `--complex-terms`.

A component with nothing triggered scores 1.0.

## Consistency

The annotation set is checked for formal violations: terms forbidden in the organism's lineage by direct GO taxon rules or by the GO-Uberon/PO/CL anatomical bridge (when a taxon is given with `--taxon`), and sets of terms that are jointly unsatisfiable under taxon constraints (checked with a SAT solver). With \( V \) violations,

\[
C_{\text{cons}} = \max(0,\ 1 - 0.01\,V).
\]

Each violation is reported with an explanation, so the consistency score doubles as a list of annotations to inspect.

## Composite audit score

The composite score of the Genome-Scale Annotation Evaluation Framework (GAEF) combines the three dimensions, averaging the coherence components that were evaluated:

\[
C_{\text{func}} = \operatorname{mean}\{ C_{\text{proc}}, C_{\text{path}}, C_{\text{cplx}} \text{ that were evaluated} \},
\qquad
\text{GAEF} = \frac{\sum_{d} w_d\, C_d}{\sum_{d} w_d},
\quad (w_{\text{ess}}, w_{\text{func}}, w_{\text{cons}}) = (0.30, 0.40, 0.30),
\]

where the sums run over the dimensions that were evaluated. When all three are available this is \( 0.30\,C_{\text{ess}} + 0.40\,C_{\text{func}} + 0.30\,C_{\text{cons}} \). A dimension that could not be evaluated (for example coherence in a run without pathway definitions and without the reasoner) neither raises nor lowers the score. The paper's formulation includes a functional-resolution term based on mean information content with weight 0; the code reports mean information content separately and does not include it in the composite.

!!! warning "The composite score measures what Stages 1 and 2 optimise"
    Completeness, coherence, and consistency are exactly the invariants that the Stage 1 priors and Stage 2A enforcement target. A higher GAEF score after GSPA refinement shows that refinement did what it was designed to do; it is not independent evidence that the annotations are more correct. Use it as an internal audit of an annotation set, and to compare annotation sets from different pipelines on the same genome. For accuracy, use held-out reference annotations and the protein-level metrics below. The published evaluation reports GAEF only alongside independent held-out GO, EC, leave-reaction-out, leave-subunit-out, fitness, and phenotype benchmarks.

## Metagenome-assembled genomes

A MAG that is 70% complete will be missing essential genes and pathway steps for reasons that have nothing to do with annotation quality, so raw completeness and coherence understate the annotation quality of incomplete bins. `gspa evaluate --mag` marks the input as a MAG. The library divides completeness by the estimated genome completeness (capped at 1.0) when a CheckM or CheckM2 completeness value is attached to the genome, and it contains a further adjuster (`gspa-core/.../metrics/MagAdjuster.groovy`) that corrects coherence for the expected fraction of missing genes, flags violations that may come from contamination, and assigns quality tiers (high: completeness \( \ge 90\% \) and contamination \( < 5\% \); medium: \( \ge 50\% \) and \( < 10\% \)). In the current command-line interface there is no option to supply CheckM estimates, so these adjustments are available only through the library; interpret raw completeness for MAGs with the estimated genome completeness in mind. See the [metagenome-assembled genomes scenario](../scenarios/mag-metagenome.md).

## Other reported quantities

The quality report also contains the mean information content of the annotations (a proxy for specificity), the number and fraction of proteins with at least one annotation, and annotation counts per source.

## Protein-level accuracy metrics

The benchmarks compare predictions with reference annotations using the standard CAFA metrics, computed on ancestor-propagated annotations over a threshold grid:

- **Protein-centric \( F_{\max} \)**: precision averaged over proteins with at least one prediction above the threshold, recall averaged over all proteins with reference annotations, and the maximum harmonic mean over thresholds.
- **Micro \( F_{\max} \)**: precision and recall pooled over all (protein, term) pairs.
- **\( WF_{\max} \)** and **\( S_{\min} \)**: information-accretion-weighted \( F_{\max} \), and the minimum semantic distance combining remaining uncertainty and misinformation (lower is better).
- **Hit@\(k\)** and **mean reciprocal rank** for candidate-ranking tasks such as finding the gene for a held-out reaction or complex subunit.

Exact definitions, thresholds, and statistical tests are given in [Benchmarks: metrics](../benchmarks/index.md#metrics).

## Limitations

- Completeness depends on the essential profile. The bundled profiles are broad domain-level cores; organisms with reduced genomes (endosymbionts, obligate parasites) can legitimately lack some essential terms.
- Pathway coherence depends on the pathway database and is only as meaningful as the pathway definitions supplied.
- The consistency score saturates at 0 after 100 violations, so it separates clean annotations from problematic ones but does not rank heavily contaminated annotation sets.
