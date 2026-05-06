package gspa.integration.promotion

import gspa.integration.GapKey
import gspa.integration.IntegrationState
import gspa.integration.suggester.SingletonSuggestion
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Beam-search promotion strategy.
 *
 * <p>For each gap, restrict to the top-{@link #candidatesPerGap} candidates
 * by log-posterior. Process gaps in priority order (highest best-candidate
 * log-posterior first). Maintain a beam of up to {@link #beamWidth} partial
 * states; each state tracks committed (gap, protein) pairs and a cumulative
 * score. At each gap, branch each state over the top-k candidates whose
 * protein is not already used in that state, prune to the top-W by score.</p>
 *
 * <p>Unlike {@link GreedyStrategy}, which commits the best candidate per
 * gap irrevocably, this strategy can commit a locally-sub-optimal candidate
 * for gap G1 if doing so frees up a higher-scoring protein for gap G2,
 * yielding a better global score.</p>
 *
 * <p>With {@code beamWidth = 1} the strategy degrades to a deterministic
 * single-path traversal that is approximately (but not exactly) greedy —
 * it differs from {@link GreedyStrategy} only in that beam-1 processes
 * gaps in priority order with at-most-one-per-protein enforcement, which
 * is functionally identical to greedy on most inputs.</p>
 *
 * <p>Score of a state = cumulative log-posterior across its commits.
 * Final output = commits from the highest-scoring state in the final beam.</p>
 */
class BeamSearchStrategy implements PromotionStrategy {

    private static final Logger log = LoggerFactory.getLogger(BeamSearchStrategy)

    /** Maximum states retained in the beam at any depth. */
    int beamWidth = 5

    /** Top-k candidates per gap explored as branching options. */
    int candidatesPerGap = 3

    /** Partial state in the beam. */
    private static class PartialState {
        List<SingletonSuggestion> commits = []
        Set<String> usedProteins = new HashSet<>()
        double score = 0.0d

        PartialState extend(SingletonSuggestion cand, double logPosterior) {
            PartialState next = new PartialState()
            next.commits = new ArrayList<>(this.commits)
            next.commits.add(cand)
            next.usedProteins = new HashSet<>(this.usedProteins)
            next.usedProteins.add(cand.proteinId)
            next.score = this.score + logPosterior
            next
        }
    }

    @Override
    List<SingletonSuggestion> select(
            List<SingletonSuggestion> candidates,
            IntegrationState state,
            int iteration) {
        if (candidates == null || candidates.isEmpty()) return []

        // Group by gap. Keep only the top-k candidates per gap.
        Map<GapKey, List<SingletonSuggestion>> byGap = [:].withDefault { [] }
        for (SingletonSuggestion ss : candidates) byGap[PromotionHelpers.gapKey(ss)] << ss

        Map<GapKey, List<SingletonSuggestion>> topK = [:]
        for (Map.Entry<GapKey, List<SingletonSuggestion>> e : byGap) {
            List<SingletonSuggestion> sorted = new ArrayList<>(e.value)
            sorted.sort { a, b -> Double.compare(
                PromotionHelpers.logPosteriorOf(b) as double,
                PromotionHelpers.logPosteriorOf(a) as double) }
            topK[e.key] = sorted.take(candidatesPerGap)
        }

        // Process gaps in priority order: highest best-candidate log-posterior first.
        List<GapKey> gapOrder = new ArrayList<>(topK.keySet())
        gapOrder.sort { g1, g2 -> Double.compare(
            PromotionHelpers.logPosteriorOf(topK[g2][0]) as double,
            PromotionHelpers.logPosteriorOf(topK[g1][0]) as double) }

        // Seed the beam with one empty state.
        List<PartialState> beam = [new PartialState()]

        for (GapKey gk : gapOrder) {
            List<PartialState> newBeam = new ArrayList<>()
            for (PartialState s : beam) {
                boolean extended = false
                for (SingletonSuggestion cand : topK[gk]) {
                    if (s.usedProteins.contains(cand.proteinId)) continue
                    double lp = PromotionHelpers.logPosteriorOf(cand)
                    newBeam.add(s.extend(cand, lp))
                    extended = true
                }
                // If every candidate's protein is already used, the state
                // can't extend on this gap — carry it forward (implicitly
                // "skip this gap").
                if (!extended) newBeam.add(s)
            }
            // Prune: keep top-beamWidth by score.
            newBeam.sort { a, b -> Double.compare(b.score, a.score) }
            beam = newBeam.size() > beamWidth ? newBeam.subList(0, beamWidth) : newBeam
        }

        PartialState best = beam[0]
        log.info("BeamSearchStrategy iter=${iteration}: ${candidates.size()} candidates → " +
            "${gapOrder.size()} gaps × top-${candidatesPerGap} → beam ${beamWidth} → " +
            "${best.commits.size()} commits, score=${String.format(Locale.ROOT, '%.3f', best.score)}")
        best.commits
    }
}
