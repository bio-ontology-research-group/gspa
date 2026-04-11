package gspa.integration

/**
 * Computes per-claim prior contributions to the posterior log-odds.
 *
 * <p>Phase 7.1 ships an empty stub so the refiner compiles. Phase 7.3
 * implements the concrete priors: {@code EssentialityPrior},
 * {@code CoherencePrior}, {@code ConsistencyPrior}, {@code GapFillingPrior},
 * {@code GenomicContextPrior}.</p>
 */
class PriorEngine {

    List<Prior> priors = []

    /** Per-prior strength hyperparameters, learned by the benchmark. */
    Map<String, Double> lambda = new LinkedHashMap<>()

    /**
     * Sum of per-prior log-odds contributions for a (protein, function).
     * Returns 0 when no priors are registered.
     */
    double totalBoost(String proteinId, String functionKey, IntegrationState state) {
        if (priors.isEmpty()) return 0.0
        double sum = 0.0
        for (Prior p : priors) {
            double l = lambda.getOrDefault(p.name(), 1.0d)
            sum += l * p.logOddsBoost(proteinId, functionKey, state)
        }
        sum
    }

    boolean isEmpty() { priors.isEmpty() }
}
