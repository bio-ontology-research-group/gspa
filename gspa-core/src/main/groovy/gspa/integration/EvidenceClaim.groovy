package gspa.integration

import groovy.transform.Canonical
import groovy.transform.builder.Builder
import groovy.transform.builder.SimpleStrategy
import gspa.model.Annotation
import gspa.model.AnnotationType

/**
 * A single piece of evidence supporting an annotation claim.
 *
 * Claims are the input to the evidence combination machinery
 * ({@code EvidenceCombiner}, {@code PriorEngine}, {@code IterativeRefiner}).
 * Each claim represents "predictor X says protein P has function F with raw
 * score S"; multiple claims for the same (protein, function) are combined
 * into a posterior probability by the combiner.
 *
 * Claims are produced from {@link gspa.model.Annotation} instances by
 * {@link ClaimExtractor}, which fills in the evidence type and calibrated
 * probability from per-source calibration tables.
 */
@Canonical
@Builder(builderStrategy = SimpleStrategy, prefix = '')
class EvidenceClaim {

    /** The protein this claim is about. */
    String proteinId

    /** Functional annotation type (GO, EC, KEGG, AMR, CAZYME, etc.). */
    AnnotationType functionType

    /** Functional identifier (e.g., GO:0006412, EC:2.7.1.1, blaCTX-M-15). */
    String functionId

    /** GO aspect (MF, BP, CC) if applicable, null otherwise. */
    String goAspect

    /** Evidence type (homology, structure, context, LM, domain-specific, ...). */
    EvidenceType evidenceType

    /** Producing predictor name, for provenance. */
    String source

    /** Raw predictor score in [0, 1] (already normalized by the predictor). */
    double rawScore

    /**
     * Calibrated probability in [0, 1]. Produced by
     * {@link CalibrationTable} from the raw score + source combination.
     * This is what enters the combiner as p = P(F | this evidence).
     */
    double calibratedProb

    /** Free-form metadata: target accession, e-value, database, etc. */
    Map<String, Object> metadata = [:]

    /** Back-reference to the original {@link Annotation} for auditability. */
    Annotation origin

    /**
     * Optional context key (e.g., operon id, pathway id) that the prior
     * engine can use to group claims.
     */
    String contextKey

    /** Key used to group claims referring to the same (protein, function). */
    String functionKey() {
        "${proteinId}|${functionType}|${functionId}".toString()
    }
}
