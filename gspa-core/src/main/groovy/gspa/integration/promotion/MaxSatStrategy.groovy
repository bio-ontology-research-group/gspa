package gspa.integration.promotion

import gspa.integration.GapKey
import gspa.integration.IntegrationState
import gspa.integration.suggester.SingletonSuggestion
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * MaxSAT-style promotion strategy with coherence-bonus coupling — now
 * implemented via greedy + 2-opt local search instead of SAT4J.
 *
 * <p>The class retains its original name for CLI/config backward
 * compatibility. Internally the algorithm is:</p>
 *
 * <ol>
 *   <li><b>Pre-filter</b>: keep only the top-{@link #candidatesPerGap}
 *       candidates per gap by log-posterior.</li>
 *   <li><b>Greedy initial solution</b>: repeatedly commit the candidate
 *       with highest <em>effective score</em>
 *       {@code logPosterior(c) + coherenceBonusWeight · |committed-in-pathway(c)|}
 *       that doesn't conflict with existing commits (same gap, same
 *       protein). Stop when no non-conflicting candidate has positive
 *       effective score.</li>
 *   <li><b>2-opt local search</b>: for each (committed c₁, uncommitted c₂)
 *       pair where swapping c₁ for c₂ respects the hard constraints,
 *       compute the joint delta (individual log-posteriors + all coherence
 *       pair changes). If strictly positive, apply the swap. Repeat until
 *       a pass finds no improvement.</li>
 * </ol>
 *
 * <p>On the rprowazekii benchmark this replaces ~5 min of SAT4J wall time
 * with sub-second local search and produces the same joint-MAP assignments
 * on tractable inputs (SAT4J was timing out anyway).</p>
 *
 * <p>Hard constraints:</p>
 * <ul>
 *   <li>At most one commit per gap.</li>
 *   <li>At most {@link #maxPerProtein} commits per protein (default 1).</li>
 * </ul>
 */
class MaxSatStrategy implements PromotionStrategy {

    private static final Logger log = LoggerFactory.getLogger(MaxSatStrategy)

    /** Max fresh pathway-function assignments per protein per outer iteration. */
    int maxPerProtein = 1

    /**
     * Coherence bonus: added once to a candidate's effective score for
     * every previously-committed candidate in the same pathway (targeting
     * a different gap). Zero (default) disables coherence coupling and the
     * strategy delegates to {@link GreedyStrategy}.
     */
    double coherenceBonusWeight = 0.0d

    /**
     * Per-gap top-k pre-filter. Limits the problem size before the local
     * search; default 3 keeps the 2-opt loop tractable on real inputs.
     */
    int candidatesPerGap = 3

    /**
     * Maximum 2-opt local-search passes. Each pass iterates all
     * (committed, uncommitted) pairs once. 20 is plenty for convergence
     * on problems up to ~500 candidates.
     */
    int maxLocalSearchPasses = 20

    /** Below this effective score a candidate is not worth committing. */
    double minEffectiveScore = 0.0d

    @Override
    List<SingletonSuggestion> select(
            List<SingletonSuggestion> candidates,
            IntegrationState state,
            int iteration) {
        if (candidates == null || candidates.isEmpty()) return []

        // Fast path: no coherence coupling → per-gap argmax = greedy.
        if (coherenceBonusWeight <= 0.0d) {
            return new GreedyStrategy().select(candidates, state, iteration)
        }

        // --- Pre-filter: top-k candidates per gap by log-posterior. ---
        Map<GapKey, List<SingletonSuggestion>> byGap = [:].withDefault { [] }
        for (SingletonSuggestion ss : candidates) {
            byGap[PromotionHelpers.gapKey(ss)] << ss
        }
        List<SingletonSuggestion> pool = new ArrayList<>()
        for (List<SingletonSuggestion> group : byGap.values()) {
            if (group.size() <= candidatesPerGap) {
                pool.addAll(group)
            } else {
                List<SingletonSuggestion> sorted = new ArrayList<>(group)
                sorted.sort { a, b -> Double.compare(
                    PromotionHelpers.logPosteriorOf(b) as double,
                    PromotionHelpers.logPosteriorOf(a) as double) }
                pool.addAll(sorted.subList(0, candidatesPerGap))
            }
        }

        // --- Working state ---
        Set<SingletonSuggestion> committed = new LinkedHashSet<>()
        Set<GapKey> usedGaps = new HashSet<>()
        Set<String> usedProteins = new HashSet<>()
        // pathway → committed candidates (for coherence-pair counting).
        Map<String, Set<SingletonSuggestion>> pathwayCommits = [:].withDefault { new LinkedHashSet<>() }

        // --- Phase 1: greedy initial solution. ---
        while (true) {
            SingletonSuggestion best = null
            double bestScore = minEffectiveScore
            for (SingletonSuggestion ss : pool) {
                if (committed.contains(ss)) continue
                if (usedGaps.contains(PromotionHelpers.gapKey(ss))) continue
                if (usedProteins.contains(ss.proteinId)) continue
                double eff = effectiveScore(ss, pathwayCommits)
                if (eff > bestScore) { bestScore = eff; best = ss }
            }
            if (best == null) break
            commit(best, committed, usedGaps, usedProteins, pathwayCommits)
        }

        int initialCommits = committed.size()

        // --- Phase 2: 2-opt local search. ---
        int passes = 0
        while (passes < maxLocalSearchPasses) {
            boolean improved = twoOptPass(pool, committed, usedGaps, usedProteins, pathwayCommits)
            passes++
            if (!improved) break
        }

        log.info("MaxSatStrategy iter=${iteration}: ${candidates.size()} candidates → " +
            "${pool.size()} after top-${candidatesPerGap} pre-filter → " +
            "${initialCommits} greedy commits → ${committed.size()} after ${passes} 2-opt passes " +
            "(coherence=${coherenceBonusWeight})")
        new ArrayList<>(committed)
    }

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------

    /** Effective score for committing s given current state. */
    private double effectiveScore(SingletonSuggestion ss,
                                  Map<String, Set<SingletonSuggestion>> pathwayCommits) {
        double base = PromotionHelpers.logPosteriorOf(ss)
        if (ss.pathwayId == null) return base
        // Coherence pairs earned against each already-committed candidate
        // in the same pathway at a different gap.
        Set<SingletonSuggestion> same = pathwayCommits[ss.pathwayId]
        if (same == null || same.isEmpty()) return base
        int pairs = 0
        GapKey myGap = PromotionHelpers.gapKey(ss)
        for (SingletonSuggestion o : same) {
            if (PromotionHelpers.gapKey(o) != myGap) pairs++
        }
        base + coherenceBonusWeight * pairs
    }

    /** Net change in the global objective from swapping out `c1` and in `c2`. */
    private double swapDelta(SingletonSuggestion c1, SingletonSuggestion c2,
                             Map<String, Set<SingletonSuggestion>> pathwayCommits) {
        // Individual score change.
        double delta = PromotionHelpers.logPosteriorOf(c2) - PromotionHelpers.logPosteriorOf(c1)

        if (coherenceBonusWeight > 0.0d) {
            // Coherence pairs lost by removing c1.
            int pairsLost = 0
            if (c1.pathwayId != null) {
                GapKey g1 = PromotionHelpers.gapKey(c1)
                for (SingletonSuggestion o : pathwayCommits.getOrDefault(c1.pathwayId, Collections.emptySet())) {
                    if (o != c1 && PromotionHelpers.gapKey(o) != g1) pairsLost++
                }
            }
            // Coherence pairs gained by adding c2 (against everyone ELSE after we remove c1).
            int pairsGained = 0
            if (c2.pathwayId != null) {
                GapKey g2 = PromotionHelpers.gapKey(c2)
                for (SingletonSuggestion o : pathwayCommits.getOrDefault(c2.pathwayId, Collections.emptySet())) {
                    if (o == c1) continue    // c1 is being removed
                    if (PromotionHelpers.gapKey(o) != g2) pairsGained++
                }
            }
            delta += coherenceBonusWeight * (pairsGained - pairsLost)
        }
        delta
    }

    /** Net change in the global objective from removing `c1` entirely. */
    private double removeDelta(SingletonSuggestion c1,
                               Map<String, Set<SingletonSuggestion>> pathwayCommits) {
        double delta = -PromotionHelpers.logPosteriorOf(c1)
        if (coherenceBonusWeight > 0.0d && c1.pathwayId != null) {
            GapKey g1 = PromotionHelpers.gapKey(c1)
            int pairsLost = 0
            for (SingletonSuggestion o : pathwayCommits.getOrDefault(c1.pathwayId, Collections.emptySet())) {
                if (o != c1 && PromotionHelpers.gapKey(o) != g1) pairsLost++
            }
            delta -= coherenceBonusWeight * pairsLost
        }
        delta
    }

    /** Try each possible swap / remove / add; apply the first improving one. */
    private boolean twoOptPass(List<SingletonSuggestion> pool,
                               Set<SingletonSuggestion> committed,
                               Set<GapKey> usedGaps,
                               Set<String> usedProteins,
                               Map<String, Set<SingletonSuggestion>> pathwayCommits) {
        // 1. Try swap: uncommit c1, commit c2 (where c2 uses c1's protein or gap slot).
        for (SingletonSuggestion c1 : new ArrayList<>(committed)) {
            GapKey g1 = PromotionHelpers.gapKey(c1)
            for (SingletonSuggestion c2 : pool) {
                if (committed.contains(c2)) continue
                if (c1 == c2) continue
                GapKey g2 = PromotionHelpers.gapKey(c2)
                // After the swap, the new state must respect AM1/protein, AM1/gap.
                // c2's gap must not be held by any OTHER commit.
                if (g1 != g2 && usedGaps.contains(g2)) continue
                // c2's protein must not be held by any OTHER commit.
                if (c1.proteinId != c2.proteinId && usedProteins.contains(c2.proteinId)) continue
                double delta = swapDelta(c1, c2, pathwayCommits)
                if (delta > 1e-9) {
                    uncommit(c1, committed, usedGaps, usedProteins, pathwayCommits)
                    commit(c2, committed, usedGaps, usedProteins, pathwayCommits)
                    return true
                }
            }
        }
        // 2. Try pure add (no conflict): commit c if effective score positive.
        for (SingletonSuggestion c : pool) {
            if (committed.contains(c)) continue
            if (usedGaps.contains(PromotionHelpers.gapKey(c))) continue
            if (usedProteins.contains(c.proteinId)) continue
            double eff = effectiveScore(c, pathwayCommits)
            if (eff > minEffectiveScore + 1e-9) {
                commit(c, committed, usedGaps, usedProteins, pathwayCommits)
                return true
            }
        }
        // 3. Try pure removal: if a commit's contribution is net-negative, drop it.
        for (SingletonSuggestion c : new ArrayList<>(committed)) {
            double delta = removeDelta(c, pathwayCommits)
            if (delta > 1e-9) {
                uncommit(c, committed, usedGaps, usedProteins, pathwayCommits)
                return true
            }
        }
        false
    }

    private void commit(SingletonSuggestion ss,
                        Set<SingletonSuggestion> committed,
                        Set<GapKey> usedGaps,
                        Set<String> usedProteins,
                        Map<String, Set<SingletonSuggestion>> pathwayCommits) {
        committed << ss
        usedGaps << PromotionHelpers.gapKey(ss)
        usedProteins << ss.proteinId
        if (ss.pathwayId) pathwayCommits[ss.pathwayId] << ss
    }

    private void uncommit(SingletonSuggestion ss,
                          Set<SingletonSuggestion> committed,
                          Set<GapKey> usedGaps,
                          Set<String> usedProteins,
                          Map<String, Set<SingletonSuggestion>> pathwayCommits) {
        committed.remove(ss)
        usedGaps.remove(PromotionHelpers.gapKey(ss))
        usedProteins.remove(ss.proteinId)
        if (ss.pathwayId) pathwayCommits[ss.pathwayId].remove(ss)
    }
}
