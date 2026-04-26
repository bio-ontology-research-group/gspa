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

    /**
     * Region-level annotations: 1-based inclusive residue indices within the
     * protein sequence. Null for whole-protein annotations (the default).
     * Populated by predictors that localise a feature to a sub-sequence
     * (e.g. SignalP cleavage site, DeepTMHMM helix, Metapredict disorder
     * stretch). {@link #regionType} is a free-form tag such as
     * {@code "disorder"}, {@code "signal_peptide"}, {@code "tm_helix"}.
     */
    Integer regionStart
    Integer regionEnd
    String regionType

    /**
     * Genomic-region annotations (v1.3): when {@link #contigId} is set, the
     * annotation refers to a region on the chromosome / contig rather than
     * to a single protein. Populated by predictors that detect prophages,
     * plasmids or whole-contig viral classifications (geNomad, CheckV,
     * PhiSpy). {@link #genomicStart} / {@link #genomicEnd} are 1-based
     * inclusive contig coordinates. Existing per-protein consumers ignore
     * these fields when {@code contigId == null}.
     */
    String contigId
    Integer genomicStart
    Integer genomicEnd

    boolean isGO() { type == AnnotationType.GO }
    boolean isEC() { type == AnnotationType.EC }
    boolean isPfam() { type == AnnotationType.PFAM }
    boolean hasRegion() { regionStart != null && regionEnd != null }
    boolean hasGenomicRegion() { contigId != null && genomicStart != null && genomicEnd != null }
}
