package gspa.metrics

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Genome
import gspa.model.Protein
import gspa.ontology.GoOntology
import gspa.ontology.SatConsistencyChecker
import groovy.util.logging.Slf4j

/**
 * Optional completeness enforcement: for each MISSING essential function,
 * promote it (an imputed annotation) onto the protein with the strongest
 * sub-threshold evidence pointing into that function's region of the GO DAG.
 *
 * A truly-missing essential F has no protein annotated with F or a descendant
 * (else it would be present), so the only available evidence is a near-ancestor
 * of F. The candidate protein is the one whose best annotation among F's
 * ancestors (within {@link #maxHops} hops up) has the highest score, decayed by
 * hop distance; the imputed score is that decayed value. If no candidate clears
 * {@link #promoteMinScore}, the gap is left unfilled (recorded, never
 * fabricated from nothing). Terms the organism cannot carry are never promoted.
 *
 * This is the lightweight, fully-explainable analogue of the integration-layer
 * {@code EssentialityPrior} / {@code DarkMatterSuggester}. Every promotion is
 * recorded with its basis for provenance.
 */
@Slf4j
class CompletenessEnforcer {

    GoOntology goOntology
    /** Optional: skip promoting any essential the asserted organism cannot carry. */
    SatConsistencyChecker satChecker
    EnforcementReport report

    String promoteEvidence = 'ISC'
    double promoteMinScore = 0.05
    int maxHops = 1
    double hopDecay = 0.5

    static class Result {
        int promoted = 0
        int unfillable = 0
        @Override
        String toString() { "CompletenessEnforcer.Result(promoted=${promoted}, unfillable=${unfillable})" }
    }

    Result enforce(Genome genome, Collection<String> missingEssentials) {
        def result = new Result()
        if (!missingEssentials) return result
        missingEssentials.each { String f ->
            if (organismForbids(f)) { result.unfillable++; return }
            def cand = PromotionCandidates.best(genome, f, goOntology, maxHops, hopDecay, promoteMinScore)
            if (cand == null) { result.unfillable++; return }
            promote(cand.protein as Protein, f, cand.score as double, cand.basis as String, result)
        }
        log.info("Completeness enforcement: promoted ${result.promoted} missing essential(s); " +
            "${result.unfillable} left unfilled (no evidence / taxon-forbidden)")
        result
    }

    private void promote(Protein p, String f, double score, String basis, Result result) {
        def ann = new Annotation(type: AnnotationType.GO, value: f, score: score,
            source: 'gspa-enforce:completeness', evidence: promoteEvidence)
        p.annotations.add(ann)
        result.promoted++
        report?.record(new EnforcementAction(
            dimension: 'completeness', action: 'promote',
            contigId: p.contig?.id, proteinId: p.id, term: f,
            reason: 'missing essential function', basis: basis,
            scoreBefore: 0.0d, scoreAfter: score), ann)
    }

    private boolean organismForbids(String f) {
        if (satChecker?.organismTaxon == null) return false
        !satChecker.check([f] as Set).consistent
    }
}
