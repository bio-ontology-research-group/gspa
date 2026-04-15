package gspa.integration.promotion

import gspa.integration.GapKey
import gspa.integration.IntegrationState
import gspa.integration.suggester.SingletonSuggestion

/**
 * Decides which {@link SingletonSuggestion}s the outer loop should promote
 * this iteration.
 *
 * <p>Three concrete strategies live in this package:</p>
 * <ul>
 *   <li>{@link AllAboveThresholdStrategy} (default) — batch-commit every
 *       singleton whose {@code q} exceeds the iteration's rising threshold.
 *       Produces multiple promotions per gap (one per candidate operon)
 *       which the current benchmark shows to be noisy.</li>
 *   <li>{@link GreedyStrategy} — rank by per-assignment log-posterior,
 *       commit the best non-conflicting batch (≤ one per gap, ≤ one per
 *       protein). Always picks the globally-best move first, leaving
 *       conflicting lower-scoring moves for subsequent outer iterations.</li>
 *   <li>{@link MaxSatStrategy} — weighted MaxSAT (SAT4J) with hard
 *       "at most one per gap" + "at most one per protein" clauses and
 *       soft clauses weighted by log-posterior. Solves the joint
 *       assignment exactly.</li>
 * </ul>
 */
interface PromotionStrategy {

    /**
     * @param candidates      all SingletonSuggestion at or above the outer
     *                        loop's rising q threshold, with their
     *                        gap not yet in {@code state.closedGaps}
     * @param state           current integration state
     * @param iteration       1-based outer-loop iteration counter (for logging)
     * @return the subset of candidates to promote this iteration
     */
    List<SingletonSuggestion> select(
        List<SingletonSuggestion> candidates,
        IntegrationState state,
        int iteration)
}

/**
 * Static helpers used by strategy implementations. Groovy 4 doesn't
 * allow static methods in interfaces, so they live in this companion
 * class.
 */
class PromotionHelpers {

    /**
     * Per-assignment log-posterior the strategies rank by. Pulls from
     * {@code proteinScores[proteinId].totalLogOdds} which the suggester
     * already computes.
     */
    static double logPosteriorOf(SingletonSuggestion ss) {
        def dec = ss.proteinScores?.get(ss.proteinId)
        dec == null ? Double.NEGATIVE_INFINITY : dec.totalLogOdds
    }

    /** Gap key identifying the (pathwayId, reactionId) a suggestion targets. */
    static GapKey gapKey(SingletonSuggestion ss) {
        new GapKey(pathwayId: ss.pathwayId, reactionId: ss.reactionId, goTerm: ss.functionId)
    }
}
