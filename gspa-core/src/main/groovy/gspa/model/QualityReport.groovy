package gspa.model

import groovy.transform.builder.Builder
import groovy.transform.builder.SimpleStrategy

/**
 * Aggregated quality assessment report for a genome's annotations.
 * Contains completeness, coherence, and consistency results.
 */
@Builder(builderStrategy = SimpleStrategy, prefix = '')
class QualityReport {

    String genomeId
    String assessmentDate

    /** Completeness: fraction of essential functions present (0.0 - 1.0) */
    double completeness = 0.0

    /** Which essential functions are present */
    Set<String> presentEssentialFunctions = []

    /** Which essential functions are missing */
    Set<String> missingEssentialFunctions = []

    /** Process coherence: fraction of has_part dependencies satisfied (0.0 - 1.0) */
    double processCoherence = 0.0

    /** Pathway coherence: average pathway completeness (0.0 - 1.0) */
    double pathwayCoherence = 0.0

    /** Complex coherence: fraction of coherent complexes (0.0 - 1.0) */
    double complexCoherence = 0.0

    /** Whether the annotation set is taxonomically consistent */
    boolean consistent = true

    /** List of consistency violations with explanations */
    List<ConsistencyViolation> violations = []

    /** Unsatisfied (C, F) has_part pairs from process-coherence — what's
     * required but missing. Each entry: key = required (typically C term
     * implied by the annotation), value = the F term that wasn't found. */
    List<Map.Entry<String, String>> incoherentProcessPairs = []

    /** Per-pathway coherence detail. Each entry: id, name, completeness,
     * required_count, present_terms, missing_terms. Used by the GAEF
     * report writer to name the under-covered pathways. */
    List<Map<String, Object>> incoherentPathways = []

    /** GO ID → human-readable name lookup. Populated by QualityScorer when
     * a GoOntology is available, then read by QualityReportWriter so the
     * JSON output names every essential / pathway / process pair instead
     * of just printing bare IDs. Empty map if no ontology was provided. */
    Map<String, String> goLabels = [:]

    /** Mean information content of annotations */
    double meanIC = 0.0

    /** Number of proteins with at least one annotation */
    int annotatedProteinCount = 0

    /** Total number of proteins */
    int totalProteinCount = 0

    /** Per-source annotation counts */
    Map<String, Integer> annotationCountBySource = [:]

    /** Composite quality score (weighted combination) */
    double getCompositeScore() {
        double coherence = (processCoherence + pathwayCoherence + complexCoherence) / 3.0
        double consistencyScore = consistent ? 1.0 : Math.max(0.0, 1.0 - violations.size() * 0.01)
        0.3 * completeness + 0.4 * coherence + 0.3 * consistencyScore
    }

    double getAnnotationCoverage() {
        totalProteinCount > 0 ? annotatedProteinCount / (double) totalProteinCount : 0.0
    }
}
