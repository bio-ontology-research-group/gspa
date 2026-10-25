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

    /** Pathway-map coverage: mean enzyme coverage of the KEGG maps and modules that have any enzyme (0.0 - 1.0).
     *  A breadth measure, not a completeness test. */
    double pathwayCoherence = 0.0

    /** Module completeness: mean fraction of assessable required steps present over the KEGG modules whose
     *  evidence activates them (0.0 - 1.0); -1.0 when module definitions were not loaded or no module was activated. */
    double moduleCoherence = -1.0d
    int modulesAssessed = 0
    int modulesActivated = 0
    int modulesComplete = 0

    /** Per-module detail (all assessed modules): id, name, class, activated, activation_score, completeness,
     *  step counts and per-step status. */
    List<Map<String, Object>> moduleDetail = []

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

    /**
     * Composite quality score (weighted combination). Coherence components are
     * stored negative when they were not evaluated (e.g. ELK-free "lite" runs,
     * which report them as N/A); folding a negative in as a real 0..1 value would
     * wrongly punish the score. So only dimensions actually computed contribute,
     * and the weights are renormalised over them — an un-evaluated dimension
     * neither helps nor penalises.
     */
    double getCompositeScore() {
        List<Double> vals = []
        List<Double> wts = []

        vals << completeness; wts << 0.3d

        double consistencyScore = consistent ? 1.0 : Math.max(0.0, 1.0 - violations.size() * 0.01)
        vals << consistencyScore; wts << 0.3d

        // Module completeness replaces map coverage in the composite when it was assessed.
        double pathwayPart = moduleCoherence >= 0.0d ? moduleCoherence : pathwayCoherence
        List<Double> coh = [processCoherence, pathwayPart, complexCoherence].findAll { it >= 0.0d }
        if (!coh.isEmpty()) {
            vals << (coh.sum() / coh.size()); wts << 0.4d
        }

        double wsum = wts.sum()
        if (wsum <= 0) return 0.0
        double acc = 0.0
        vals.eachWithIndex { v, i -> acc += wts[i] * v }
        acc / wsum
    }

    double getAnnotationCoverage() {
        totalProteinCount > 0 ? annotatedProteinCount / (double) totalProteinCount : 0.0
    }
}
