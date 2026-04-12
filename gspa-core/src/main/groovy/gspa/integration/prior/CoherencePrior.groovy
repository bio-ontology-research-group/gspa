package gspa.integration.prior

import gspa.integration.IntegrationState
import gspa.integration.Prior
import gspa.metrics.Coherence
import gspa.metrics.ProcessCoherenceResult
import gspa.ontology.PathwayGraph
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Boosts candidate claims that would close a triggered-but-unsatisfied
 * process dependency, or complete a partially-covered pathway.
 *
 * <p>Plan §A.4 of Phase 7. Uses
 * {@link Coherence#evaluateProcessCoherence(Set)} to find the unsatisfied
 * {@code (C, F)} has_part pairs, then for every candidate claim
 * {@code (protein, F)} returns {@code +alphaCoh * (1 - fracAnnotated)} —
 * where {@code fracAnnotated} reflects how much of the triggering process
 * is already covered (closer to completion → bigger boost).</p>
 *
 * <p>Pathway coherence: for each triggered pathway with partial coverage,
 * boost candidates whose GO term is in the pathway's required set,
 * scaled by {@code (1 - pathwayCoverage)}.</p>
 *
 * <p>No-ops cleanly when either {@code goReasoner} or
 * {@code pathwayDatabase} is missing from state.</p>
 */
class CoherencePrior implements Prior {

    private static final Logger log = LoggerFactory.getLogger(CoherencePrior)

    /** Per-match boost magnitude (log-odds). Learned by the benchmark. */
    double alphaCoh = 1.0d

    /**
     * Cached per iteration: unsatisfied F terms that would complete a
     * process dependency, each paired with the triggering C term and the
     * fraction already satisfied (so the boost can be weighted by
     * (1 - f)).
     */
    private Map<String, Double> processMissingTerms = new LinkedHashMap<>()

    /** Cached per iteration: pathway GO terms needed (term → (1 - coverage)). */
    private Map<String, Double> pathwayMissingTerms = new LinkedHashMap<>()

    @Override
    String name() { 'coherence' }

    @Override
    Set<String> inputs() { ['goReasoner', 'pathwayDatabase', 'posteriors'] as Set }

    @Override
    void beginIteration(IntegrationState state) {
        processMissingTerms = new LinkedHashMap<>()
        pathwayMissingTerms = new LinkedHashMap<>()

        Set<String> annotated = state.goReasoner != null
            ? state.currentlyAnnotatedGoTermsPropagated()
            : state.currentlyAnnotatedGoTerms()

        // --- Process coherence: needs GoReasoner for has_part pairs ---
        if (state.goReasoner != null) try {
            def coherence = new Coherence(state.goOntology, state.goReasoner)
            ProcessCoherenceResult result = coherence.evaluateProcessCoherence(annotated)
            int triggered = result.triggered
            int satisfied = result.satisfied
            double fracAnnotated = triggered == 0 ? 1.0 : (satisfied / (double) triggered)
            double weight = Math.max(0.0d, 1.0d - fracAnnotated)

            for (Map.Entry<String, String> unsat : result.unsatisfied) {
                // Each entry: C present but F missing. Boost F.
                String fTerm = unsat.value
                if (fTerm != null) {
                    processMissingTerms.merge(fTerm, weight, { a, b -> Math.max(a as double, b as double) })
                }
            }
        } catch (Exception ex) {
            log.warn("CoherencePrior: process-coherence check failed: ${ex.message}")
        }

        // --- Pathway coherence: for each pathway, find missing required terms ---
        if (state.pathwayDatabase != null) {
            try {
                for (PathwayGraph pathway : state.pathwayDatabase.pathways.values()) {
                    Set<String> required = pathway.getRequiredGoTerms()
                    if (required == null || required.isEmpty()) continue
                    // Is this pathway "triggered"? Use any-member-annotated as a soft trigger.
                    boolean triggeredPw = required.any { annotated.contains(it) }
                    if (!triggeredPw) continue
                    int requiredCount = required.size()
                    int annotatedCount = required.count { annotated.contains(it) }
                    double coverage = (double) annotatedCount / requiredCount
                    double weight = Math.max(0.0d, 1.0d - coverage)
                    if (weight <= 0.0d) continue
                    for (String term : required) {
                        if (annotated.contains(term)) continue
                        pathwayMissingTerms.merge(term, weight, { a, b -> Math.max(a as double, b as double) })
                    }
                }
            } catch (Exception ex) {
                log.warn("CoherencePrior: pathway check failed: ${ex.message}")
            }
        }

        log.debug("CoherencePrior: ${processMissingTerms.size()} process-missing, ${pathwayMissingTerms.size()} pathway-missing terms")
    }

    @Override
    double logOddsBoost(String proteinId, String functionKey, IntegrationState state) {
        String[] parts = IntegrationState.splitFunctionKey(functionKey)
        if (parts == null || parts[1] != 'GO') return 0.0d
        String goId = parts[2]

        double processWeight = processMissingTerms.getOrDefault(goId, 0.0d)
        double pathwayWeight = pathwayMissingTerms.getOrDefault(goId, 0.0d)
        // Take the max of the two signals; they are not independent.
        double weight = Math.max(processWeight, pathwayWeight)
        if (weight <= 0.0d) return 0.0d
        alphaCoh * weight
    }
}
