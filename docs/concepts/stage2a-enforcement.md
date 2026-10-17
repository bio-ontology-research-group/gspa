# Symbolic enforcement (Stage 2A)

Stage 1 penalises taxonomically and anatomically inadmissible terms softly, so a very strong likelihood can survive. Stage 2A makes these constraints hard. After Stage 1 has converged it removes every remaining violation of the organism's taxonomic lineage and anatomical body plan, and by default it does so by *relaxing* each forbidden term to its most specific admissible GO ancestors instead of deleting it. It also attenuates isolated subunits of obligate heteromeric complexes.

Stage 2A runs when `gspa integrate` is given `--enforce-stage2`. Taxon-level enforcement needs an organism taxon (`--taxon` or `--infer-taxon`); the bundled constraint tables are then loaded automatically.

## The cross-ontology constraint bridge

Many restrictions on GO terms are not stated on the GO term itself. *Heart development* (GO:0007507) is only meaningful in organisms that have a heart, but the taxon restriction is attached to the Uberon class *heart*, not to the GO process. GSPA therefore ships five constraint tables, extracted from the OBO Foundry releases of GO, Uberon, Plant Ontology (PO), Cell Ontology (CL), and NCBI Taxonomy, and bundled in `gspa-core/src/main/resources/`:

| Table | Rows | Content |
|---|---|---|
| `taxon-constraints/go-taxon-constraints.tsv` | 26,911 | direct GO taxon rules: 13,819 `only_in_taxon`, 13,092 `never_in_taxon`, over 18,159 GO terms |
| `anatomical-constraints/go-anatomical-constraints.tsv` | 2,715 | equivalence and subclass axioms linking 2,674 GO terms to 1,294 Uberon, PO, and CL classes over 30 relations (for example `results_in_development_of`, `results_in_morphogenesis_of`, `occurs_in`) |
| `anatomical-constraints/uberon-po-taxon-constraints.tsv` | 8,331 | taxon rules on 1,247 anatomical classes (4,767 `only_in_taxon`, 3,564 `never_in_taxon`) |
| `anatomical-constraints/uberon-po-hierarchy.tsv` | 7,722 | `is_a` and `part_of` edges over 2,855 anatomical classes |
| `taxon-constraints/ncbi-taxon-hierarchy.tsv` | 1,818 | lineage edges over 1,305 NCBI Taxonomy clades |

A GO term \( t \) is *anatomically forbidden* in taxon \( \tau \) when \( t \), or a GO ancestor of \( t \), is linked to an anatomical class \( u \) such that \( u \) or one of its anatomical ancestors carries an `only_in_taxon` rule that excludes the lineage of \( \tau \), or a `never_in_taxon` rule that includes it. A term is *taxonomically forbidden* when it, or an ancestor, violates a direct GO taxon rule for the lineage. Without an asserted taxon, GSPA falls back to checking each protein's own set of annotations for mutually unsatisfiable terms with a SAT solver.

Alternative tables can be supplied with `--taxon-constraints` and `--taxonomy`; see [Bundled resources](../reference/bundled-resources.md).

## Enforcement

Let \( \mathcal{F}(\tau) \) be the union of taxonomically and anatomically forbidden terms. For every protein \( p \) and every GO annotation \( t \) of \( p \) with \( t \in \mathcal{F}(\tau) \) or an ancestor of \( t \) in \( \mathcal{F}(\tau) \), Stage 2A applies the action selected by `--enforce-consistency-mode`:

| Mode | Action |
|---|---|
| `relax` (default) | semantic ancestor relaxation (below) |
| `remove` | delete the annotation |
| `downrank` | multiply the score by 0.5 and tag it `taxon-inconsistent` |
| `flag` | keep the score, tag it `taxon-inconsistent` |
| `minimal-flip` | resolve co-annotation conflicts jointly by weighted MaxSAT, dropping the lower-weight term of each conflicting pair until the set is consistent |

### Semantic ancestor relaxation

Deleting a forbidden term also discards the valid biology implied by its ancestors. An invertebrate ortholog of a vertebrate organ-development regulator that receives *heart development* most likely still participates in a broader developmental or morphogenetic process. Relaxation keeps that information.

For a forbidden term \( t \), Stage 2A walks \( \mathrm{Anc}(t) \) along `is_a` and `part_of`, excluding the three GO aspect roots, and collects the minimal set \( \mathcal{A}_{\text{valid}}(t, \tau) \) of most specific ancestors that are both taxonomically and anatomically admissible in \( \tau \). It reconstructs the support the term had before the Stage 1 conflict penalties, \( P_{\text{unpen}}(p, t) = \sigma\bigl(L^{(0)}(p, t) + \sum_{m} B_m^{+}(p, t)\bigr) \), where the sum keeps all prior contributions except negative consistency and anatomical penalties, and transfers it:

\[
P_{2A}(p, a) = \max\bigl(P^{(K)}(p, a),\ \gamma_{\text{relax}}\, P_{\text{unpen}}(p, t)\bigr)\quad \forall a \in \mathcal{A}_{\text{valid}}(t, \tau), \qquad P_{2A}(p, t) = 0,
\]

with \( \gamma_{\text{relax}} = 1.0 \). If no admissible non-root ancestor exists, the annotation is removed. An ancestor that the protein already carries keeps the larger of its own score and the transferred score.

### Obligate complex singleton pruning

For each GO term classified as an obligate heteromeric complex in the complex classification (`--complex-terms`, bundled default), if exactly one protein in the genome has \( P \ge 0.50 \) and no other protein reaches the candidate gate \( P \ge 0.08 \), that single subunit cannot be part of a complete assembly in this genome. Stage 2A lowers its probability to 0.45, just below the decision threshold. The Stage 1 coherence prior applies a soft version of the same logic; this step makes it final.

## Hyperparameters

| Parameter | Default | Control |
|---|---|---|
| Enforcement mode | `relax` | `--enforce-consistency-mode` |
| Relaxation factor \( \gamma_{\text{relax}} \) | 1.0 | fixed |
| Downrank factor | 0.5 | fixed |
| Singleton prune: active threshold, candidate gate, attenuated value | 0.50, 0.08, 0.45 | fixed |
| Use the anatomical bridge | on when the `anatomical` prior is enabled | `--enable-priors` |

## Outputs and audit trail

Every enforcement action is recorded with the protein, term, reason (for example `anatomical-inconsistent for NCBITaxon_7227: ...`), the evidence basis, and the scores before and after. `--enforcement-actions-out FILE` writes them as a TSV with the columns `dimension`, `action`, `contig`, `protein_id`, `term`, `term_label`, `reason`, `basis`, `score_before`, `score_after` (see [Input and output formats](../reference/file-formats.md)). In the integrated output, relaxed ancestors carry the provenance key `stage2_relax`, pruned singletons carry `stage2_singleton_prune`, and the supporting claims of the original forbidden term are kept on the relaxed ancestor.

## Limitations

- Anatomical enforcement covers only the 2,715 GO terms that currently have Uberon, PO, or CL logical definitions in the bridge. Processes without such a definition are checked only against direct GO taxon rules.
- The bundled NCBI Taxonomy table covers 1,305 clades. For organisms whose lineage is not covered, supply a fuller hierarchy with `--taxonomy`.
- In unicellular organisms and non-vascular plants, where an injected organ-level term usually lands on a housekeeping protein, deletion gave slightly higher precision than relaxation in the published evaluation, while relaxation gave higher recall and information content. Use `--enforce-consistency-mode remove` if precision matters most for such genomes.
- Enforcement is deterministic and complete with respect to the tables; it is not a check of biological truth. A genuine horizontal transfer of a lineage-restricted function will be relaxed or removed.

!!! note "Version note"
    In the GSPA 2.1.0 release jar, a property-name error in `gspa integrate` could leave the SAT checker unwired, so that Stage 2A enforcement was silently skipped. This page describes the intended behaviour, which is what the current source implements.
