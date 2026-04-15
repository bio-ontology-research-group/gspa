package gspa.integration.promotion

import gspa.integration.GapKey
import gspa.integration.IntegrationState
import gspa.integration.suggester.SingletonSuggestion
import org.sat4j.core.Vec
import org.sat4j.core.VecInt
import org.sat4j.maxsat.WeightedMaxSatDecorator
import org.sat4j.pb.PseudoOptDecorator
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
     * Solver timeout in seconds. On timeout MaxSAT falls back to
     * {@link GreedyStrategy}. Empirically SAT4J can't handle the
     * coherence-reified problem at real-benchmark scale (370 candidates
     * × O(pathways²) pair clauses) in any reasonable time; 60s is enough
     * to find a first feasible solution on small problems and fail
     * quickly on infeasible-within-budget ones.
     */
    int timeoutSeconds = 60

    /**
     * Restrict to the top-k candidates per gap before building the MaxSAT
     * problem. Without this, real benchmark inputs (thousands of candidates)
     * exceed SAT4J's capacity. Default 3 — tight enough for sub-second
     * solves on ~300-candidate inputs.
     */
    int candidatesPerGap = 3

    /**
     * Coherence bonus: reward weight for jointly committing two candidates in
     * the same pathway (but different gaps). Encoded as a pairwise aux
     * variable {@code y_ij ↔ (x_i ∧ x_j)} with a soft clause on y_ij.
     *
     * <p>Zero (default) disables coherence coupling and MaxSAT reduces to
     * per-gap argmax = greedy. Positive values let MaxSAT prefer "complete
     * this pathway more" over "best individual score per gap", which can
     * pick locally-sub-optimal candidates that close more of a pathway
     * jointly.</p>
     */
    double coherenceBonusWeight = 0.0d

    /**
     * Cap on the number of pairwise coherence clauses per pathway.
     * Empirically 50 was still too much for SAT4J — cut to 10.
     */
    int coherenceBonusPairCap = 10

    @Override
    List<SingletonSuggestion> select(
            List<SingletonSuggestion> candidates,
            IntegrationState state,
            int iteration) {
        if (candidates == null || candidates.isEmpty()) return []

        // Fast path: with no coherence coupling the MaxSAT problem
        // decomposes to per-gap argmax + per-protein AM1 — which is
        // exactly what GreedyStrategy computes. Skip SAT4J entirely;
        // it was taking 5 min per iteration on 370-candidate inputs
        // even without coherence clauses.
        if (coherenceBonusWeight <= 0.0d) {
            return new GreedyStrategy().select(candidates, state, iteration)
        }

        // Pre-filter: keep only the top-k candidates per gap by log-posterior.
        // Without this, SAT4J MaxSAT routinely times out on ~1500-candidate
        // problems. With candidatesPerGap=5 a typical problem shrinks by 5-10×.
        Map<GapKey, List<SingletonSuggestion>> perGapGroups = [:].withDefault { [] }
        for (SingletonSuggestion ss : candidates) {
            perGapGroups[PromotionHelpers.gapKey(ss)] << ss
        }
        List<SingletonSuggestion> prefiltered = new ArrayList<>()
        for (List<SingletonSuggestion> group : perGapGroups.values()) {
            if (group.size() <= candidatesPerGap) {
                prefiltered.addAll(group)
            } else {
                List<SingletonSuggestion> sorted = new ArrayList<>(group)
                sorted.sort { a, b -> Double.compare(
                    PromotionHelpers.logPosteriorOf(b) as double,
                    PromotionHelpers.logPosteriorOf(a) as double) }
                prefiltered.addAll(sorted.subList(0, candidatesPerGap))
            }
        }
        candidates = prefiltered

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

        // ---- Group by gap / protein / pathway. ----
        Map<GapKey, List<Integer>> byGap = [:].withDefault { [] }
        Map<String, List<Integer>> byProtein = [:].withDefault { [] }
        Map<String, List<Integer>> byPathway = [:].withDefault { [] }
        for (int i = 0; i < n; i++) {
            byGap[PromotionHelpers.gapKey(candidates[i])] << i
            byProtein[candidates[i].proteinId] << i
            if (candidates[i].pathwayId) byPathway[candidates[i].pathwayId] << i
        }

        // Pre-compute which candidate pairs would carry a coherence bonus.
        // A pair (i,j) qualifies if both belong to the same pathway AND
        // target DIFFERENT gaps. Per-pathway we cap the pair count at
        // coherenceBonusPairCap and take the top by sum(scores[i]+scores[j]).
        List<int[]> coherencePairs = []
        if (coherenceBonusWeight > 0.0d) {
            for (List<Integer> pwGroup : byPathway.values()) {
                if (pwGroup.size() < 2) continue
                List<int[]> pwPairs = []
                for (int a = 0; a < pwGroup.size(); a++) {
                    for (int b = a + 1; b < pwGroup.size(); b++) {
                        int ia = pwGroup[a], ib = pwGroup[b]
                        if (PromotionHelpers.gapKey(candidates[ia]) == PromotionHelpers.gapKey(candidates[ib])) continue
                        pwPairs << ([ia, ib] as int[])
                    }
                }
                if (pwPairs.size() > coherenceBonusPairCap) {
                    pwPairs.sort { p1, p2 -> Double.compare(
                        scores[p2[0]] + scores[p2[1]],
                        scores[p1[0]] + scores[p1[1]]) }
                    pwPairs = pwPairs.subList(0, coherenceBonusPairCap)
                }
                coherencePairs.addAll(pwPairs)
            }
        }

        // ---- Build the weighted MaxSAT problem. ----
        int nAux = coherencePairs.size()
        int totalVars = n + nAux
        long bonusClauseWeight = Math.max(1L, Math.round(coherenceBonusWeight * weightScale))

        WeightedMaxSatDecorator solver
        try {
            solver = new WeightedMaxSatDecorator(SolverFactory.newDefault())
            solver.setTopWeight(new BigInteger(String.valueOf(weightScale * 1_000_000L)))
            solver.newVar(totalVars)                   // 1-indexed: vars 1..n primary, n+1..n+nAux aux
            solver.setExpectedNumberOfClauses(byGap.size() + byProtein.size() + n + nAux * 4)
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

            // Coherence bonus (optional): for each (i, j) pair in the same pathway
            // targeting different gaps, add aux var y = x_i ∧ x_j and a soft
            // clause rewarding y. This encodes "committing both together is
            // worth an extra bonusClauseWeight beyond their individual sum".
            if (nAux > 0) {
                for (int p = 0; p < nAux; p++) {
                    int i = coherencePairs[p][0]
                    int j = coherencePairs[p][1]
                    int y = n + 1 + p                   // aux var id
                    // y → x_i
                    solver.addHardClause(new VecInt([-y, i + 1] as int[]))
                    // y → x_j
                    solver.addHardClause(new VecInt([-y, j + 1] as int[]))
                    // (x_i ∧ x_j) → y   ≡   ¬x_i ∨ ¬x_j ∨ y
                    solver.addHardClause(new VecInt([-(i + 1), -(j + 1), y] as int[]))
                    // Soft: y with bonus weight.
                    solver.addSoftClause(bonusClauseWeight, new VecInt([y] as int[]))
                }
            }
        } catch (ContradictionException ce) {
            log.warn("MaxSatStrategy iter=${iteration}: contradiction building problem; falling back to greedy")
            return new GreedyStrategy().select(candidates, state, iteration)
        }

        // ---- Wrap with PseudoOptDecorator to get the IOptimizationProblem
        // interface (admitABetterSolution / discardCurrentSolution / model).
        // WeightedMaxSatDecorator alone doesn't expose these. ----
        PseudoOptDecorator opt = new PseudoOptDecorator(solver)
        opt.setTimeout(timeoutSeconds)

        // Each admitABetterSolution()+discard iteration finds a strictly-better
        // solution; capture model() BEFORE discardCurrentSolution() because
        // the latter invalidates the current assignment.
        int[] bestModel = null
        try {
            while (opt.admitABetterSolution()) {
                bestModel = opt.model()
                opt.discardCurrentSolution()
            }
        } catch (ContradictionException ce) {
            // Expected at optimum.
        } catch (Exception ex) {
            log.warn("MaxSatStrategy iter=${iteration}: solver exception (${ex.class.simpleName}: ${ex.message}); falling back to greedy")
            return new GreedyStrategy().select(candidates, state, iteration)
        }
        if (bestModel == null) {
            log.warn("MaxSatStrategy iter=${iteration}: no solution; falling back to greedy")
            return new GreedyStrategy().select(candidates, state, iteration)
        }
        int[] model = bestModel
        // SAT4J model: positive var = true, negative = false. Literal indices
        // are 1..n; we map back to candidates[literal-1].
        List<SingletonSuggestion> out = new ArrayList<>()
        for (int lit : model) {
            if (lit > 0 && lit <= n) out.add(candidates[lit - 1])
        }
        log.info("MaxSatStrategy iter=${iteration}: ${n} candidates + ${nAux} coherence aux vars " +
            "(bonus weight=${coherenceBonusWeight}) → ${out.size()} commits " +
            "(hard: ≤1/gap, ≤${maxPerProtein}/protein)")
        out
    }
}
