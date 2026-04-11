package gspa.integration

import groovy.transform.Canonical
import gspa.integration.suggester.Suggestion
import gspa.model.Annotation
import gspa.model.AnnotationSet
import gspa.model.AnnotationType

/**
 * Result of evidence integration: a final set of annotations with full
 * provenance for every posterior assignment.
 *
 * The {@link #annotations} field holds the concrete annotations written
 * back to the genome (scores = posterior probability). The
 * {@link #provenance} map records, for each function key, how the posterior
 * was built: which claims contributed, what the likelihood log-odds was,
 * which priors fired, and the final posterior.
 *
 * Phase 8's dark-matter suggestions are attached to
 * {@link #suggestions} as a separate channel — they do NOT enter
 * {@link #annotations} unless the user enables one of the collapse modes.
 */
@Canonical
class IntegratedAnnotationSet {

    /** Final annotations, one per (protein, function) that survived integration. */
    AnnotationSet annotations = new AnnotationSet()

    /** Provenance indexed by function key ("proteinId|type|functionId"). */
    Map<String, ClaimProvenance> provenance = new LinkedHashMap<>()

    /** Phase 8 suggestions; empty until the suggester runs. */
    List<Suggestion> suggestions = []

    void put(ClaimProvenance prov, Annotation ann) {
        provenance[prov.functionKey] = prov
        annotations.add(ann)
    }

    ClaimProvenance provenanceFor(String proteinId, AnnotationType type, String functionId) {
        provenance["${proteinId}|${type}|${functionId}".toString()]
    }

    int size() { annotations.size() }

    @Override
    String toString() {
        "IntegratedAnnotationSet(${annotations.size()} annotations, ${provenance.size()} provenance entries, ${suggestions.size()} suggestions)"
    }
}
