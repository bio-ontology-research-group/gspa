package gspa.integration.promotion

import gspa.integration.GapKey
import gspa.integration.IntegrationState
import gspa.integration.suggester.SingletonSuggestion
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Greedy best-first promotion strategy (option a from the Phase 10 design
 * discussion).
 *
 * <p>Sorts all candidates by log-posterior descending and greedily commits
 * the highest-scoring non-conflicting subset. Two suggestions conflict if
 * they share:</p>
 * <ul>
 *   <li>a gap (pathwayId, reactionId) — at most one protein per gap, and</li>
 *   <li>a protein — at most one new promotion per protein per iteration,
 *       to prevent one protein being assigned many fresh pathway functions
 *       at once</li>
 * </ul>
 *
 * <p>This is a conflict-free-batch greedy (not strict one-per-iteration
 * best-first). With {@code maxPerIteration = 1} it degrades to strict
 * best-first. The default batches as many non-conflicting commits as
 * possible to keep the outer loop tractable on large genomes while
 * retaining the property that every committed promotion is the best
 * currently-available move for its gap.</p>
 */
class GreedyStrategy implements PromotionStrategy {

    private static final Logger log = LoggerFactory.getLogger(GreedyStrategy)

    /**
     * Maximum number of commits per outer iteration. Setting this to 1
     * yields strict greedy best-first; {@link Integer#MAX_VALUE} is
     * conflict-free batch greedy (default).
     */
    int maxPerIteration = Integer.MAX_VALUE

    @Override
    List<SingletonSuggestion> select(
            List<SingletonSuggestion> candidates,
            IntegrationState state,
            int iteration) {
        if (candidates == null || candidates.isEmpty()) return []

        // Sort by log-posterior descending.
        List<SingletonSuggestion> sorted = new ArrayList<>(candidates)
        sorted.sort { a, b -> Double.compare(
            PromotionHelpers.logPosteriorOf(b) as double,
            PromotionHelpers.logPosteriorOf(a) as double) }

        Set<GapKey> usedGaps = new HashSet<>()
        Set<String> usedProteins = new HashSet<>()
        List<SingletonSuggestion> out = new ArrayList<>()

        for (SingletonSuggestion ss : sorted) {
            if (out.size() >= maxPerIteration) break
            GapKey gk = PromotionHelpers.gapKey(ss)
            if (usedGaps.contains(gk)) continue
            if (usedProteins.contains(ss.proteinId)) continue
            out.add(ss)
            usedGaps.add(gk)
            usedProteins.add(ss.proteinId)
        }

        log.info("GreedyStrategy iter=${iteration}: ${candidates.size()} candidates → ${out.size()} conflict-free commits")
        out
    }
}
