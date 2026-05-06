package gspa.integration

import gspa.model.Annotation
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Runs the evidence-integration fixed-point loop.
 *
 * <p>Phase 7.3 adds the full multi-iteration loop:</p>
 * <ol>
 *   <li>Seed with likelihood-only posteriors.</li>
 *   <li>For each iteration: ask each prior to recompute its cached state,
 *       then recompute every posterior as
 *       {@code L_new = L_likelihood + sum_k lambda_k * prior_k.boost()}.</li>
 *   <li>Apply Jacobi-style under-relaxation with {@code damping} to
 *       prevent oscillation when non-monotone priors (Consistency)
 *       flip claims.</li>
 *   <li>Stop when mean {@code |Δp|} over all function keys falls below
 *       {@code epsilon}, or when {@code maxIter} is reached, or when
 *       divergence is detected (two successive iterations with
 *       increasing {@code Δp}).</li>
 * </ol>
 *
 * <p>When {@link #priorEngine} is empty, the loop collapses to a single
 * likelihood pass — the Phase 7.1/7.2 behaviour — and returns immediately
 * after the first iteration.</p>
 */
class IterativeRefiner {

    private static final Logger log = LoggerFactory.getLogger(IterativeRefiner)

    EvidenceCombiner combiner
    PriorEngine priorEngine = new PriorEngine()

    int maxIter = 6
    double epsilon = 0.005
    double damping = 0.5

    IterativeRefiner(EvidenceCombiner combiner) {
        this.combiner = combiner
    }

    /**
     * Refine a set of claims into posterior annotations.
     */
    IntegratedAnnotationSet refine(List<EvidenceClaim> claims, IntegrationState state) {
        Map<String, List<EvidenceClaim>> byKey = ClaimExtractor.groupByFunctionKey(claims)
        log.info("Refining ${claims.size()} claims across ${byKey.size()} (protein, function) groups (maxIter=${maxIter}, damping=${damping})")

        // Precompute likelihood contribution per group — invariant across iterations.
        Map<String, Double> likelihood = new LinkedHashMap<>()
        for (Map.Entry<String, List<EvidenceClaim>> entry : byKey.entrySet()) {
            likelihood[entry.key] = combiner.combineLikelihood(entry.value)
        }

        // Seed posteriors from likelihood alone.
        Map<String, Double> posteriorLogOdds = new LinkedHashMap<>(likelihood)
        state.updatePosteriors(posteriorLogOdds)

        int iterationsRun = 0
        double prevDelta = Double.POSITIVE_INFINITY
        int consecutiveIncreases = 0
        Map<String, Double> lastStableSnapshot = new LinkedHashMap<>(posteriorLogOdds)

        // Always run at least one iteration with priors, since the seed
        // above is likelihood-only. If there are no priors, the loop
        // exits after the first pass because delta will be 0.
        int maxIterations = priorEngine.isEmpty() ? 1 : maxIter

        for (int iter = 0; iter < maxIterations; iter++) {
            priorEngine.beginIteration(state)

            Map<String, Double> newLogOdds = new LinkedHashMap<>()
            for (Map.Entry<String, List<EvidenceClaim>> entry : byKey.entrySet()) {
                String key = entry.key
                double lLik = likelihood[key]
                double lPri = priorEngine.totalBoost(entry.value.first().proteinId, key, state)
                double lNew = clip(lLik + lPri, key, state)

                double lOld = posteriorLogOdds.getOrDefault(key, lLik)
                // Jacobi-style under-relaxation: new = (1-d) * old + d * computed.
                double lDamped = (1.0d - damping) * lOld + damping * lNew
                // Re-apply the pin floor after damping. Damping can drag a
                // floored lNew back toward an unfloored lOld (seed = raw
                // likelihood on iteration 1), which would silently undo the
                // pin. Clamp to [floor, lMax] if a floor exists.
                Double floor = state.pinnedFloors != null ? state.pinnedFloors[key] : null
                if (floor != null && lDamped < floor) {
                    lDamped = Math.min(combiner.lMax, (double) floor)
                }
                newLogOdds[key] = lDamped
            }

            double delta = meanAbsDelta(newLogOdds, posteriorLogOdds)
            log.debug("iter ${iter}: mean |Δp| = ${String.format(Locale.ROOT, '%.5f', delta)}")

            posteriorLogOdds = newLogOdds
            state.updatePosteriors(posteriorLogOdds)
            iterationsRun = iter + 1

            if (delta < epsilon) {
                log.info("Converged at iteration ${iter} (Δp=${String.format(Locale.ROOT, '%.5f', delta)} < ${epsilon})")
                break
            }

            // Divergence detection: two consecutive iterations with rising delta → roll back.
            if (delta > prevDelta) {
                consecutiveIncreases++
                if (consecutiveIncreases >= 2) {
                    log.warn("Divergence detected at iteration ${iter}; rolling back to last stable state")
                    posteriorLogOdds = lastStableSnapshot
                    state.updatePosteriors(posteriorLogOdds)
                    break
                }
            } else {
                consecutiveIncreases = 0
                lastStableSnapshot = new LinkedHashMap<>(posteriorLogOdds)
            }
            prevDelta = delta
        }

        buildResult(byKey, likelihood, posteriorLogOdds, state, iterationsRun)
    }

    /**
     * Build the final IntegratedAnnotationSet with provenance.
     */
    private IntegratedAnnotationSet buildResult(
            Map<String, List<EvidenceClaim>> byKey,
            Map<String, Double> likelihood,
            Map<String, Double> posteriorLogOdds,
            IntegrationState state,
            int iterations) {
        IntegratedAnnotationSet out = new IntegratedAnnotationSet()
        for (Map.Entry<String, List<EvidenceClaim>> entry : byKey.entrySet()) {
            String key = entry.key
            List<EvidenceClaim> group = entry.value
            double lLik = likelihood[key]
            double lFinal = posteriorLogOdds[key]
            double pFinal = sigmoid(lFinal)

            Map<String, Double> priorContribs =
                priorEngine.isEmpty()
                    ? Collections.<String, Double>emptyMap()
                    : priorEngine.perPriorBoost(group.first().proteinId, key, state)

            ClaimProvenance prov = new ClaimProvenance(
                functionKey: key,
                proteinId: group.first().proteinId,
                supportingClaims: new ArrayList<EvidenceClaim>(group),
                priorContributions: priorContribs,
                likelihoodLogOdds: lLik,
                finalLogOdds: lFinal,
                finalProbability: pFinal,
                convergenceIter: iterations,
            )

            Annotation ann = buildAnnotation(group.first(), pFinal)
            out.put(prov, ann)
        }
        out
    }

    private static Annotation buildAnnotation(EvidenceClaim claim, double posteriorProb) {
        new Annotation(
            type: claim.functionType,
            value: claim.functionId,
            score: posteriorProb,
            source: 'integrated',
            evidence: 'IEA',
            goAspect: claim.goAspect,
            evidenceType: claim.evidenceType,
        )
    }

    /**
     * Clip log-odds to the combiner's [lMin, lMax] range, then raise to
     * the pinned floor if one is set for this function key. The floor
     * preserves the "closed stays closed" invariant across Phase 10 outer
     * iterations: a DarkMatter-promoted claim cannot be driven below its
     * promoted posterior by later prior adjustments.
     */
    private double clip(double logOdds, String functionKey, IntegrationState state) {
        double clamped = Math.min(combiner.lMax, Math.max(combiner.lMin, logOdds))
        if (state != null && state.pinnedFloors != null) {
            Double floor = state.pinnedFloors[functionKey]
            if (floor != null && floor > clamped) return Math.min(combiner.lMax, floor)
        }
        clamped
    }

    private static double meanAbsDelta(Map<String, Double> a, Map<String, Double> b) {
        if (a.isEmpty() && b.isEmpty()) return 0.0d
        Set<String> keys = new LinkedHashSet<>()
        keys.addAll(a.keySet())
        keys.addAll(b.keySet())
        double sum = 0.0d
        for (String k : keys) {
            double pA = sigmoid(a.getOrDefault(k, 0.0d))
            double pB = sigmoid(b.getOrDefault(k, 0.0d))
            sum += Math.abs(pA - pB)
        }
        sum / keys.size()
    }

    private static double sigmoid(double x) {
        if (x >= 500.0) return 1.0d
        if (x <= -500.0) return 0.0d
        1.0d / (1.0d + Math.exp(-x))
    }
}
