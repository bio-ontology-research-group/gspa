package gspa.integration

import gspa.model.Genome

/**
 * Mutable state threaded through the iterative refinement loop.
 *
 * Phase 7.1 holds only the genome reference and the current posterior
 * log-odds map (keyed by function key "proteinId|type|functionId"). Phase
 * 7.3 extends this with ontology handles ({@code GoOntology},
 * {@code GoReasoner}, {@code SatConsistencyChecker}), the essential
 * function list, pathway DB, operon assignments — everything the priors
 * need to look up.
 */
class IntegrationState {

    /** The genome being integrated. */
    Genome genome

    /** Current posterior log-odds by function key. */
    Map<String, Double> posteriorLogOdds = new LinkedHashMap<>()

    IntegrationState(Genome genome) {
        this.genome = genome
    }

    /** Get current posterior log-odds for a function key; 0 if absent. */
    double get(String functionKey) {
        posteriorLogOdds.getOrDefault(functionKey, 0.0d)
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
}
