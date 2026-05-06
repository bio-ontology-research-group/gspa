package gspa.integration

import gspa.integration.promotion.AllAboveThresholdStrategy
import gspa.integration.promotion.PromotionStrategy
import gspa.integration.suggester.DarkMatterSuggester
import gspa.integration.suggester.SingletonSuggestion
import gspa.model.AnnotationType
import gspa.ontology.PathwayDatabase
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Phase 10 outer fixed-point loop over {@link IterativeRefiner} +
 * {@link DarkMatterSuggester}.
 *
 * <p>Each outer iteration:</p>
 * <ol>
 *   <li>Runs the suggester against the current state.</li>
 *   <li>Promotes singleton suggestions whose {@code q} beats a rising
 *       threshold {@code q_base + q_step·k}.</li>
 *   <li>Pins each promoted claim's posterior floor in {@code state} so
 *       {@link IterativeRefiner#clip} cannot push it below.</li>
 *   <li>Marks the gap (pathway, reaction) as closed in {@code state}.</li>
 *   <li>Re-runs the inner refiner with the augmented claim list.</li>
 *   <li>Recomputes the gap set via the strategy
 *       ({@code GapseqRescorer} when wired; falls back to
 *       {@link GapRecomputer} coverage-only).</li>
 * </ol>
 *
 * <p>Terminates when no new promotions are emitted AND the gap set is
 * stable across two consecutive iterations. Hard cap {@code maxIter}.</p>
 *
 * <p>Cascade safety: if {@code |new_promotions|} grows in two consecutive
 * iterations, the loop halts and rolls back to the last non-increasing
 * state.</p>
 */
class OuterIterativeRefiner {

    private static final Logger log = LoggerFactory.getLogger(OuterIterativeRefiner)

    /** Strategy for recomputing the gap set after promotions. */
    interface GapSource {
        /**
         * Update the state's gap list in place (or return a new one).
         * Return value is the updated set used for convergence testing.
         *
         * @return a set of {@link GapKey} representing the current open gaps
         */
        Set<GapKey> recompute(IntegrationState state, List<EvidenceClaim> promotionClaims)
    }

    /** Default {@link GapSource} that uses {@link GapRecomputer} coverage-only. */
    static class CoverageGapSource implements GapSource {
        double tauCover = 0.5d
        PathwayDatabase pathwayDb
        @Override Set<GapKey> recompute(IntegrationState state, List<EvidenceClaim> promotionClaims) {
            PathwayDatabase db = pathwayDb ?: state.pathwayDatabase
            GapRecomputer.recompute(state, db, tauCover)
        }
    }

    IterativeRefiner innerRefiner
    DarkMatterSuggester suggester = new DarkMatterSuggester()
    GapSource gapSource = new CoverageGapSource()
    PromotionStrategy promotionStrategy = new AllAboveThresholdStrategy()

    int maxIter = 5
    double qBase = 0.5d
    double qStep = 0.05d
    double qCap = 0.75d

    /** Pin policy: if false, promoted claims don't get floors. */
    boolean pinPromotions = true

    OuterIterativeRefiner(IterativeRefiner innerRefiner) {
        this.innerRefiner = innerRefiner
    }

    /**
     * Result of one {@link #refine} call: the final integrated set plus
     * a trace of each outer iteration's activity.
     */
    static class OuterResult {
        IntegratedAnnotationSet integrated
        int outerIterationsRun
        List<Integer> promotedPerIter = []
        List<Integer> gapsPerIter = []
        boolean fixedPointReached
        boolean cascadeRolledBack
    }

    /**
     * Run the outer loop on the given initial claim set.
     *
     * @param initialClaims the claim set from predictors (unchanged through iterations)
     * @param state         state pre-populated with pathways/operons/gaps/ontology handles
     * @return {@link OuterResult}
     */
    OuterResult refine(List<EvidenceClaim> initialClaims, IntegrationState state) {
        OuterResult result = new OuterResult()

        // Iteration 0: Phase-7 integration with the initial claim set.
        IntegratedAnnotationSet integrated = innerRefiner.refine(initialClaims, state)
        suggester.suggest(state, integrated)

        Set<GapKey> currentGaps = new LinkedHashSet<>(collectGaps(state))
        Set<GapKey> prevGaps = new LinkedHashSet<>(currentGaps)

        List<EvidenceClaim> promotionClaims = []
        int prevPromoted = Integer.MAX_VALUE
        int risingStreak = 0

        // Snapshot for cascade rollback: keep the last non-increasing state.
        Map<String, Double> safePosteriors = new LinkedHashMap<>(state.posteriorLogOdds)
        Map<String, Double> safePinnedFloors = new LinkedHashMap<>(state.pinnedFloors)
        Set<GapKey> safeClosedGaps = new LinkedHashSet<>(state.closedGaps)
        IntegratedAnnotationSet safeIntegrated = integrated

        for (int k = 1; k <= maxIter; k++) {
            double qThreshold = Math.min(qCap, qBase + qStep * k)

            // Collect singleton candidates that beat the rising threshold
            // AND whose gap isn't already closed. The promotion strategy
            // then decides which of these to actually commit this iteration.
            List<SingletonSuggestion> eligible = []
            for (def sug : integrated.suggestions) {
                if (!(sug instanceof SingletonSuggestion)) continue
                SingletonSuggestion ss = sug as SingletonSuggestion
                if (ss.q < qThreshold) continue
                GapKey gk = new GapKey(pathwayId: ss.pathwayId, reactionId: ss.reactionId, goTerm: null)
                if (state.isGapClosed(gk)) continue
                eligible << ss
            }
            List<SingletonSuggestion> newPromotions = promotionStrategy.select(eligible, state, k)

            log.info("Outer iter ${k}: qThreshold=${String.format(Locale.ROOT, '%.3f', qThreshold)}, " +
                "suggestions=${integrated.suggestions?.size() ?: 0}, new promotions=${newPromotions.size()}, open gaps=${currentGaps.size()}")

            result.promotedPerIter << newPromotions.size()
            result.gapsPerIter << currentGaps.size()

            // Fixed point: no new promotions AND the gap set is stable.
            if (newPromotions.isEmpty() && currentGaps == prevGaps) {
                result.fixedPointReached = true
                result.outerIterationsRun = k - 1
                break
            }

            // Cascade guard: |new_promotions| rising two iterations in a row → roll back.
            if (newPromotions.size() > prevPromoted) {
                risingStreak++
                if (risingStreak >= 2) {
                    log.warn("Cascade detected at outer iter ${k} (promoted ${prevPromoted} → ${newPromotions.size()}); rolling back to last safe state")
                    state.updatePosteriors(safePosteriors)
                    state.pinnedFloors = safePinnedFloors
                    state.closedGaps = safeClosedGaps
                    integrated = safeIntegrated
                    result.cascadeRolledBack = true
                    result.outerIterationsRun = k - 1
                    break
                }
            } else {
                risingStreak = 0
                // This iteration was non-increasing — safe to snapshot.
                safePosteriors = new LinkedHashMap<>(state.posteriorLogOdds)
                safePinnedFloors = new LinkedHashMap<>(state.pinnedFloors)
                safeClosedGaps = new LinkedHashSet<>(state.closedGaps)
                safeIntegrated = integrated
            }
            prevPromoted = newPromotions.size()

            // Apply promotions: pin floors, mark gaps closed, add claims.
            for (SingletonSuggestion ss : newPromotions) {
                GapKey gk = new GapKey(pathwayId: ss.pathwayId, reactionId: ss.reactionId, goTerm: ss.functionId)
                state.markGapClosed(gk)

                ClaimKey ck = new ClaimKey(
                    proteinId: ss.proteinId,
                    functionType: ss.functionType ?: AnnotationType.GO,
                    functionId: ss.functionId,
                )
                double qClamped = Math.max(1.0e-6d, Math.min(1.0d - 1.0e-6d, ss.q))
                double floorLogOdds = Math.log(qClamped / (1.0d - qClamped))
                if (pinPromotions) state.setPinnedFloor(ck, floorLogOdds)

                promotionClaims << buildPromotionClaim(ss, qClamped)
            }

            // Re-run Phase 7 with augmented claims (initial + all promotions so far).
            List<EvidenceClaim> augmented = new ArrayList<>(initialClaims.size() + promotionClaims.size())
            augmented.addAll(initialClaims)
            augmented.addAll(promotionClaims)
            integrated = innerRefiner.refine(augmented, state)
            suggester.suggest(state, integrated)

            // Recompute gap set via the strategy (coverage or gapseq-rescoring).
            prevGaps = currentGaps
            currentGaps = new LinkedHashSet<>(gapSource.recompute(state, promotionClaims))

            result.outerIterationsRun = k
        }

        result.integrated = integrated
        result
    }

    /** Gap set as seen by the state at iteration 0 (before any promotion). */
    private static Set<GapKey> collectGaps(IntegrationState state) {
        Set<GapKey> out = new LinkedHashSet<>()
        if (state.metabolicGaps == null) return out
        for (MetabolicGap g : state.metabolicGaps) {
            out << new GapKey(pathwayId: g.pathwayId, reactionId: g.reactionId, goTerm: g.goTerm)
        }
        out
    }

    /** Convert a promoted suggestion into an EvidenceClaim carrying DARK_MATTER type. */
    private static EvidenceClaim buildPromotionClaim(SingletonSuggestion ss, double prob) {
        new EvidenceClaim(
            proteinId: ss.proteinId,
            functionType: ss.functionType ?: AnnotationType.GO,
            functionId: ss.functionId,
            goAspect: null,
            evidenceType: EvidenceType.DARK_MATTER,
            source: 'dark_matter',
            rawScore: prob,
            calibratedProb: prob,
            metadata: [
                pathwayId: ss.pathwayId,
                reactionId: ss.reactionId,
                suggestionScore: ss.suggestionScore,
            ] as Map<String, Object>,
        )
    }
}
