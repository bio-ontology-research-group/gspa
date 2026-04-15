package gspa.integration.promotion

import gspa.integration.IntegrationState
import gspa.integration.suggester.SingletonSuggestion

/**
 * Current Phase 10 default: commit every candidate unchanged.
 *
 * Kept as the backward-compatible baseline so existing benchmark
 * numbers remain reproducible. New runs should prefer
 * {@link GreedyStrategy} or {@link MaxSatStrategy} — both empirically
 * emit one promotion per gap rather than one-per-(gap, operon).
 */
class AllAboveThresholdStrategy implements PromotionStrategy {

    @Override
    List<SingletonSuggestion> select(
            List<SingletonSuggestion> candidates,
            IntegrationState state,
            int iteration) {
        new ArrayList<>(candidates ?: [])
    }
}
