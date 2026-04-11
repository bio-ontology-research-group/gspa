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

    /**
     * Per-prior breakdown for one (protein, function). Used to populate
     * {@link ClaimProvenance#priorContributions}.
     */
    Map<String, Double> perPriorBoost(String proteinId, String functionKey, IntegrationState state) {
        Map<String, Double> out = new LinkedHashMap<>()
        for (Prior p : priors) {
            double l = lambda.getOrDefault(p.name(), 1.0d)
            double raw = p.logOddsBoost(proteinId, functionKey, state)
            if (raw != 0.0) out[p.name()] = l * raw
        }
        out
    }

    /**
     * Notify every prior that a new refinement iteration is starting.
     * Priors can recompute cached per-iteration state here.
     */
    void beginIteration(IntegrationState state) {
        for (Prior p : priors) {
            p.beginIteration(state)
        }
    }

    void register(Prior prior, double strength = 1.0) {
        priors << prior
        lambda[prior.name()] = strength
    }

    boolean isEmpty() { priors.isEmpty() }
}
