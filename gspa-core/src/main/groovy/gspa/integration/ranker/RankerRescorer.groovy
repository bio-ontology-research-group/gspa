package gspa.integration.ranker

import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.integration.suggester.DisjunctiveSuggestion
import gspa.integration.suggester.SingletonSuggestion
import gspa.integration.suggester.Suggestion
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Phase 12 M3: apply a learned ranker to per-candidate feature vectors
 * and re-score suggestions.
 *
 * <p>The feature vector we pass to the ranker must match the order in
 * which it was trained (from the RLC suggester's <code>--features-out</code>
 * TSV). For the M3 initial model the columns are:</p>
 * <pre>
 * log_density, diversity, commitment, self, strand_consistency,
 * log_intergenic, n_anchors, n_nbr_gos, n_candidates, rlc_rank,
 * rlc_score, rlc_q
 * </pre>
 *
 * <p>The suggester's <code>proteinScores</code> map already holds
 * <code>piR = density</code> and <code>q</code>. For other features we
 * expect the caller to pre-populate them into a parallel
 * <code>candidateFeatures</code> map on the IntegratedAnnotationSet —
 * or (simpler, used here) to have emitted the features to disk already,
 * in which case a separate Python pass does the rescoring.</p>
 *
 * <p>For the in-JVM path we currently only have access to q and density
 * in-memory. A full in-JVM rescorer will hook into the suggester's
 * feature computation directly in a follow-up; for M3 initial eval we
 * can call the ranker externally from Python via the LRO harness.</p>
 */
class RankerRescorer {

    private static final Logger log = LoggerFactory.getLogger(RankerRescorer)

    Ranker ranker

    /** Re-score each suggestion's candidates using the ranker. */
    IntegratedAnnotationSet rescore(IntegrationState state,
                                    IntegratedAnnotationSet integrated,
                                    Map<String, double[]> candFeatures) {
        if (ranker == null || integrated.suggestions == null) return integrated
        for (Suggestion s : integrated.suggestions) {
            List<String> pids = (s instanceof SingletonSuggestion)
                ? [((SingletonSuggestion) s).proteinId]
                : (s instanceof DisjunctiveSuggestion ? ((DisjunctiveSuggestion) s).proteinIds : [])
            if (pids.size() < 2) continue
            Map<String, Double> adj = new LinkedHashMap<>()
            for (String p : pids) {
                double[] fv = candFeatures[p]
                if (fv == null) continue
                adj[p] = ranker.score(fv)
            }
            if (adj.isEmpty()) continue
            // softmax over ranker scores
            double maxS = adj.values().max()
            double z = 0.0d
            Map<String, Double> qs = [:]
            for (Map.Entry<String, Double> e : adj) {
                double v = Math.exp(e.value - maxS)
                qs[e.key] = v
                z += v
            }
            for (Map.Entry<String, Double> e : qs) e.value = e.value / z
            if (s instanceof DisjunctiveSuggestion) {
                def ds = (DisjunctiveSuggestion) s
                List<String> ordered = qs.entrySet().sort { a, b -> Double.compare(b.value, a.value) }*.key
                ds.proteinIds = ordered
                ds.qValues = ordered.collect { qs[it] }
            }
        }
        integrated
    }
}
