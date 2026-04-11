package gspa.model

import groovy.transform.Canonical
import groovy.transform.builder.Builder
import groovy.transform.builder.SimpleStrategy
import gspa.integration.EvidenceType

/**
 * Represents a single functional annotation assigned to a protein.
 * Annotations can come from different sources (predictors) and have
 * different types (GO, EC, KEGG, Pfam, COG, etc.).
 */
@Canonical
@Builder(builderStrategy = SimpleStrategy, prefix = '')
class Annotation {

    /** The annotation type */
    AnnotationType type

    /** The annotation value (e.g., GO:0006412, PF00001, EC:2.7.1.1) */
    String value

    /** Confidence score from the predictor (0.0 - 1.0) */
    double score = 0.0

    /** Which predictor produced this annotation */
    String source

    /** Evidence code (e.g., IEA, ISS, ISO for GO annotations) */
    String evidence

    /** Free-form metadata */
    Map<String, Object> metadata = [:]

    /** For GO annotations: the GO aspect (MF, BP, CC) */
    String goAspect

    /**
     * Optional evidence type classification used by the Phase 7 integration
     * layer. Predictors may set this directly; otherwise {@code ClaimExtractor}
     * falls back to a source→type lookup table.
     */
    EvidenceType evidenceType

    boolean isGO() { type == AnnotationType.GO }
    boolean isEC() { type == AnnotationType.EC }
    boolean isPfam() { type == AnnotationType.PFAM }
}
