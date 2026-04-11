package gspa.integration

import gspa.config.EssentialFunctions
import gspa.model.AnnotationType
import gspa.model.Genome
import gspa.ontology.GoOntology
import gspa.ontology.GoReasoner
import gspa.ontology.PathwayDatabase
import gspa.ontology.SatConsistencyChecker

/**
 * Mutable state threaded through the iterative refinement loop.
 *
 * <p>Phase 7.1/7.2 held only the genome reference and the current posterior
 * log-odds map. Phase 7.3 extends this with the ontology / metric handles
 * that the concrete priors need, plus the pre-computed operon groupings
 * and metabolic gaps that genomic-context and gap-filling priors read.</p>
 *
 * <p>All extension fields are optional — a prior that needs
 * {@code GoReasoner} but finds it {@code null} returns zero boost,
 * effectively no-oping itself. This keeps the integrator usable even
 * before the caller has wired every handle.</p>
 */
class IntegrationState {

    /** Maximum-posterior decision threshold for "what's currently annotated".
     *  This is the only threshold used internally; it reflects the Bayesian
     *  MAP rule (a term is considered "present" when its posterior > 0.5),
     *  not a tunable hyperparameter. Plan §A.1 of Phase 8. */
    static final double MAP_THRESHOLD = 0.5d

    /** The genome being integrated. */
    Genome genome

    /** Current posterior log-odds by function key. */
    Map<String, Double> posteriorLogOdds = new LinkedHashMap<>()

    // --- Ontology / metric handles (optional) ---

    GoOntology goOntology
    GoReasoner goReasoner
    SatConsistencyChecker satConsistencyChecker
    EssentialFunctions essentialFunctions
    PathwayDatabase pathwayDatabase

    // --- Pre-computed inputs from Phase 3 / 5 predictors ---

    /** Operons, each as a list of protein IDs (from OperonPredictor). */
    List<List<String>> operons = []

    /** Metabolic gaps (from gapseq). */
    List<MetabolicGap> metabolicGaps = []

    // --- Taxonomy (for ConsistencyPrior) ---

    /** NCBI taxonomy lineage of the genome, parent → child. Optional. */
    Map<String, String> taxonomyParentMap

    IntegrationState(Genome genome) {
        this.genome = genome
    }

    /** Get current posterior log-odds for a function key; 0 if absent. */
    double get(String functionKey) {
        posteriorLogOdds.getOrDefault(functionKey, 0.0d)
    }

    /** Current posterior probability for a function key via sigmoid. */
    double getProbability(String functionKey) {
        sigmoid(posteriorLogOdds.getOrDefault(functionKey, 0.0d))
    }

    /** Set / overwrite a posterior log-odds value. */
    void set(String functionKey, double logOdds) {
        posteriorLogOdds[functionKey] = logOdds
    }

    /** Replace the whole posterior map in one shot (used by refiner). */
    void updatePosteriors(Map<String, Double> newPosteriors) {
        posteriorLogOdds = newPosteriors
    }

    int size() { posteriorLogOdds.size() }

    /**
     * Parse a function key "proteinId|type|functionId" into its parts.
     * @return [proteinId, type, functionId] or null if the key is malformed.
     */
    static String[] splitFunctionKey(String key) {
        if (key == null) return null
        String[] parts = key.split('\\|', 3)
        return parts.length == 3 ? parts : null
    }

    /**
     * Build a function key.
     */
    static String functionKey(String proteinId, AnnotationType type, String functionId) {
        "${proteinId}|${type}|${functionId}".toString()
    }

    /**
     * Set of GO terms currently annotated anywhere in the genome with
     * posterior > {@link #MAP_THRESHOLD}. Not propagated through the
     * hierarchy — priors that need the propagated set should call
     * {@link #currentlyAnnotatedGoTermsPropagated()}.
     */
    Set<String> currentlyAnnotatedGoTerms(double threshold = MAP_THRESHOLD) {
        Set<String> out = new LinkedHashSet<>()
        for (Map.Entry<String, Double> entry : posteriorLogOdds.entrySet()) {
            if (sigmoid(entry.value) <= threshold) continue
            String[] parts = splitFunctionKey(entry.key)
            if (parts == null) continue
            if (parts[1] == 'GO') out << parts[2]
        }
        out
    }

    /**
     * As {@link #currentlyAnnotatedGoTerms} but with upward closure via
     * the GO hierarchy (is_a + part_of). Returns the raw set if no
     * {@code goOntology} is wired. Used for completeness/coherence.
     */
    Set<String> currentlyAnnotatedGoTermsPropagated(double threshold = MAP_THRESHOLD) {
        Set<String> raw = currentlyAnnotatedGoTerms(threshold)
        if (goOntology == null) return raw
        goOntology.propagateAnnotations(raw)
    }

    /** Set of function keys for a given protein. */
    Set<String> functionKeysForProtein(String proteinId) {
        String prefix = proteinId + '|'
        Set<String> out = new LinkedHashSet<>()
        for (String key : posteriorLogOdds.keySet()) {
            if (key.startsWith(prefix)) out << key
        }
        out
    }

    /** Convenience: check if a protein has any annotation above the MAP threshold. */
    boolean hasAnyStrongAnnotation(String proteinId, double threshold = MAP_THRESHOLD) {
        for (String key : functionKeysForProtein(proteinId)) {
            if (sigmoid(posteriorLogOdds[key]) > threshold) return true
        }
        false
    }

    private static double sigmoid(double x) {
        if (x >= 500.0) return 1.0d
        if (x <= -500.0) return 0.0d
        1.0d / (1.0d + Math.exp(-x))
    }
}
