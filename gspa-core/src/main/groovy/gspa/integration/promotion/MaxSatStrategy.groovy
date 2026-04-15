package gspa.integration.promotion

import gspa.integration.GapKey
import gspa.integration.IntegrationState
import gspa.integration.suggester.SingletonSuggestion
import org.sat4j.core.Vec
import org.sat4j.core.VecInt
import org.sat4j.maxsat.WeightedMaxSatDecorator
import org.sat4j.pb.SolverFactory
import org.sat4j.specs.ContradictionException
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Weighted-MaxSAT promotion strategy (option c from the Phase 10 design
 * discussion).
 *
 * <p>Variables: one boolean {@code x_i} per {@link SingletonSuggestion} —
 * "commit this suggestion".</p>
 *
 * <p>Hard clauses:</p>
 * <ul>
 *   <li>For each gap G, at most one {@code x_i} where suggestion i targets G
 *       (cardinality ≤ 1).</li>
 *   <li>For each protein p, at most {@link #maxPerProtein} {@code x_i} where
 *       suggestion i assigns to p (default 1: no protein gets more than one
 *       fresh pathway function per outer iteration).</li>
 * </ul>
 *
 * <p>Soft clauses: for each suggestion i, one unit clause {@code x_i}
 * with weight = shifted log-posterior. Weight discretization: SAT4J's
 * {@code WeightedMaxSatDecorator} requires non-negative integer weights,
 * so we shift the log-posteriors by {@code −min(logPosteriors)} and scale
 * by {@link #weightScale} (default 1000 → 3-digit precision on log-odds).</p>
 *
 * <p>Solving this yields the globally-optimal per-iteration assignment
 * that is consistent with the hard constraints. For the typical benchmark
 * problem size (≤ a few thousand candidates) SAT4J converges in
 * milliseconds.</p>
 */
class MaxSatStrategy implements PromotionStrategy {

    private static final Logger log = LoggerFactory.getLogger(MaxSatStrategy)

    /** Max fresh pathway-function assignments per protein per iteration. */
    int maxPerProtein = 1

    /** Precision multiplier when discretizing log-posterior weights for SAT4J. */
    int weightScale = 1000

    /**
     * Solver timeout in seconds. If SAT4J doesn't converge in time we
     * fall back to {@link GreedyStrategy}.
     */
    int timeoutSeconds = 30

    @Override
    List<SingletonSuggestion> select(
            List<SingletonSuggestion> candidates,
            IntegrationState state,
            int iteration) {
        if (candidates == null || candidates.isEmpty()) return []

        int n = candidates.size()

        // ---- Compute weights (non-negative integers). ----
        double[] scores = new double[n]
        double minScore = Double.POSITIVE_INFINITY
        double maxScore = Double.NEGATIVE_INFINITY
        for (int i = 0; i < n; i++) {
            double s = PromotionHelpers.logPosteriorOf(candidates[i])
            if (!Double.isFinite(s)) s = -1e6d
            scores[i] = s
            if (s < minScore) minScore = s
            if (s > maxScore) maxScore = s
        }
        double shift = (minScore < 0) ? -minScore + 0.1d : 0.1d
        long[] weights = new long[n]
        for (int i = 0; i < n; i++) {
            double shifted = scores[i] + shift
            long w = Math.max(1L, Math.round(shifted * weightScale))
            weights[i] = w
        }

        // ---- Group by gap / protein for cardinality constraints. ----
        Map<GapKey, List<Integer>> byGap = [:].withDefault { [] }
        Map<String, List<Integer>> byProtein = [:].withDefault { [] }
        for (int i = 0; i < n; i++) {
            byGap[PromotionHelpers.gapKey(candidates[i])] << i
            byProtein[candidates[i].proteinId] << i
        }

        // ---- Build the weighted MaxSAT problem. ----
        WeightedMaxSatDecorator solver
        try {
            solver = new WeightedMaxSatDecorator(SolverFactory.newDefault())
            solver.setTopWeight(new BigInteger(String.valueOf(weightScale * 1_000_000L)))
            solver.newVar(n)                           // 1-indexed: vars 1..n
            solver.setExpectedNumberOfClauses(byGap.size() + byProtein.size() + n)
            solver.setTimeout(timeoutSeconds)

            // Hard: at most one commit per gap.
            for (List<Integer> group : byGap.values()) {
                if (group.size() < 2) continue
                int[] lits = new int[group.size()]
                for (int k = 0; k < group.size(); k++) lits[k] = group[k] + 1
                solver.addAtMost(new VecInt(lits), 1)
            }

            // Hard: at most maxPerProtein commits per protein.
            for (List<Integer> group : byProtein.values()) {
                if (group.size() <= maxPerProtein) continue
                int[] lits = new int[group.size()]
                for (int k = 0; k < group.size(); k++) lits[k] = group[k] + 1
                solver.addAtMost(new VecInt(lits), maxPerProtein)
            }

            // Soft: unit clause x_i with weight[i].
            for (int i = 0; i < n; i++) {
                VecInt clause = new VecInt([i + 1] as int[])
                solver.addSoftClause(weights[i], clause)
            }
        } catch (ContradictionException ce) {
            log.warn("MaxSatStrategy iter=${iteration}: contradiction building problem; falling back to greedy")
            return new GreedyStrategy().select(candidates, state, iteration)
        }

        // ---- Solve. ----
        try {
            boolean hasSol = false
            while (solver.admitABetterSolution()) {
                hasSol = true
                solver.discardCurrentSolution()
            }
            if (!hasSol) {
                log.warn("MaxSatStrategy iter=${iteration}: no solution; falling back to greedy")
                return new GreedyStrategy().select(candidates, state, iteration)
            }
        } catch (ContradictionException ce) {
            // Expected when the optimum is reached.
        } catch (Exception ex) {
            log.warn("MaxSatStrategy iter=${iteration}: solver exception (${ex.class.simpleName}: ${ex.message}); falling back to greedy")
            return new GreedyStrategy().select(candidates, state, iteration)
        }

        int[] model = solver.model()
        if (model == null) {
            return new GreedyStrategy().select(candidates, state, iteration)
        }
        // SAT4J model: positive var = true, negative = false. Literal indices
        // are 1..n; we map back to candidates[literal-1].
        List<SingletonSuggestion> out = new ArrayList<>()
        for (int lit : model) {
            if (lit > 0 && lit <= n) out.add(candidates[lit - 1])
        }
        log.info("MaxSatStrategy iter=${iteration}: ${n} candidates → ${out.size()} commits (hard: ≤1/gap, ≤${maxPerProtein}/protein)")
        out
    }
}
