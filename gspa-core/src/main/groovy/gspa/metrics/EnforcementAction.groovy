package gspa.metrics

import groovy.transform.Canonical

/**
 * A single enforcement action taken on a predicted annotation, recorded for
 * provenance. Promotions and demotions/removals are both captured here so the
 * output makes clear HOW each function ended up assigned (or un-assigned) and
 * on what basis (which predictor / partner / constraint).
 */
@Canonical
class EnforcementAction {

    /** 'consistency' | 'completeness' | 'coherence' */
    String dimension

    /** 'remove' | 'downrank' | 'flag' | 'promote' */
    String action

    String contigId
    String proteinId

    /** The GO term acted upon. */
    String term

    /** Human-readable label of {@link #term} (rdfs:label), for display. */
    String termLabel

    /** Human-readable reason (e.g. 'taxon-inconsistent for NCBITaxon_2', 'missing essential function'). */
    String reason

    /**
     * Evidential basis for the action: the predictor + score the decision rested
     * on, a partner annotation, or the constraint that fired. Examples:
     * {@code 'diamond@0.31'}, {@code 'partner GO:0006415 on protein NC_..._042'},
     * {@code 'organism NCBITaxon_2'}.
     */
    String basis

    double scoreBefore = 0.0
    double scoreAfter = 0.0

    /** One provenance-trail line for the affected annotation. */
    String provenanceLine() {
        String b = basis ? "; basis=${basis}" : ''
        "${dimension}:${action}(${reason}${b})"
    }
}
