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

    /** Annotations indexed by the same biological identity as provenance. */
    Map<String, Annotation> annotationsByFunction = new LinkedHashMap<>()

    /** Includes proteins whose predictions were subsequently removed. */
    Set<String> processedProteinIds = new LinkedHashSet<>()

    void put(ClaimProvenance prov, Annotation ann) {
        provenance[prov.functionKey] = prov
        processedProteinIds.add(prov.proteinId)
        Annotation previous = annotationsByFunction.put(prov.functionKey, ann)
        if (previous != null) annotations.annotations.removeAll { it.is(previous) }
        ann.score = prov.finalProbability
        ann.metadata = new LinkedHashMap<>(ann.metadata ?: [:])
        ann.metadata.protein_id = prov.proteinId
        ann.metadata.function_key = prov.functionKey
        annotations.add(ann)
    }

    /**
     * Provenance is the authoritative score store. Rebuild dependent views after
     * a stage adds, removes, relaxes or re-scores a prediction. Never associate
     * genes and annotations by their list positions: two isozymes may carry the
     * same EC, while relaxation may remove or add GO terms between stages.
     */
    void synchronize(IntegrationState state = null) {
        Map<String, Annotation> next = new LinkedHashMap<>()
        Map<String, Double> posteriors = new LinkedHashMap<>()
        provenance.each { String key, ClaimProvenance prov ->
            String[] parts = IntegrationState.splitFunctionKey(key)
            if (parts == null) throw new IllegalArgumentException("Invalid function key: ${key}")
            if (!Double.isFinite(prov.finalProbability) || prov.finalProbability < 0d || prov.finalProbability > 1d) {
                throw new IllegalArgumentException("Invalid probability for ${key}: ${prov.finalProbability}")
            }
            prov.functionKey = key
            prov.proteinId = parts[0]
            processedProteinIds.add(parts[0])
            double bounded = Math.max(1e-12d, Math.min(1d - 1e-12d, prov.finalProbability))
            prov.finalLogOdds = Math.log(bounded / (1d - bounded))
            def claim = prov.supportingClaims ? prov.supportingClaims.first() : null
            Annotation ann = annotationsByFunction[key] ?: new Annotation(
                type: AnnotationType.valueOf(parts[1]), value: parts[2],
                source: 'integrated', evidence: 'IEA',
                goAspect: claim?.goAspect, evidenceType: claim?.evidenceType)
            ann.score = prov.finalProbability
            ann.metadata = new LinkedHashMap<>(ann.metadata ?: [:])
            ann.metadata.protein_id = parts[0]
            ann.metadata.function_key = key
            next[key] = ann
            posteriors[key] = prov.finalLogOdds
        }
        annotationsByFunction = next
        annotations = new AnnotationSet(new ArrayList<Annotation>(next.values()))
        if (state != null) state.updatePosteriors(posteriors)
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
