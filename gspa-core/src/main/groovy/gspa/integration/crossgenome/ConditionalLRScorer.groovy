package gspa.integration.crossgenome

import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.integration.suggester.DisjunctiveSuggestion
import gspa.integration.suggester.SingletonSuggestion
import gspa.integration.suggester.Suggestion
import groovy.transform.CompileStatic
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Phase 12 M2: re-score suggestions from the Reaction-Local Context
 * suggester using the cross-genome conditional-LR statistic per
 * (orthogroup, reaction).
 *
 * <p>The scorer does NOT introduce new candidates — it only re-weights
 * existing suggestions by multiplying each candidate's score by
 * {@code LR(orthogroup(p), R)^λ}. Candidates whose orthogroup has
 * insufficient cross-genome evidence (small {@code nSigTotal} or CI
 * overlapping LR=1) keep their within-genome score.</p>
 *
 * <p>When the re-weighting changes the rank order, singleton-vs-
 * disjunctive decisions are re-made via the same top-q &gt; 0.5 rule.</p>
 */
class CrossGenomeReScorer {

    private static final Logger log = LoggerFactory.getLogger(CrossGenomeReScorer)

    /** Exponent on LR. λ = 0 disables; λ = 1 is the natural prior. */
    double lambda = 1.0d

    /** Min {@code nSigTotal} to trust the LR. */
    int minSupport = 3

    /** Keep only LRs whose 90% CI excludes 1.0 (log-CI excludes 0). */
    boolean requireCredible = true

    /** Singleton promotion threshold on re-weighted top q. */
    double singletonThreshold = 0.5d

    /** Re-scored suggestions score cap (matches DM convention). */
    double maxSuggestionScore = 0.85d

    /**
     * Re-score suggestions in-place on {@code integrated.suggestions}.
     * Returns the same list with updated scores / membership / kinds.
     */
    IntegratedAnnotationSet rescore(IntegrationState state,
                                    IntegratedAnnotationSet integrated,
                                    ReactionLocusCatalog catalog) {
        if (integrated.suggestions == null || integrated.suggestions.isEmpty()) {
            return integrated
        }
        if (catalog == null || catalog.size() == 0) {
            log.info("CrossGenomeReScorer: empty catalog; pass-through")
            return integrated
        }
        Map<String, String> ortho = state.orthogroupMap
        if (ortho == null) {
            log.info("CrossGenomeReScorer: no orthogroup map; pass-through")
            return integrated
        }

        // Build EC → reactions lookup from the reaction graph (for bridging
        // gap reaction ids to the catalog's SEED reaction ids).
        def rxnGraph = state.reactionGraph
        int rescored = 0
        List<Suggestion> out = []
        for (Suggestion s : integrated.suggestions) {
            String rxn = s.reactionId
            if (!rxn) { out << s; continue }

            // Bridge: if direct catalog lookup fails, map gap-reaction to all
            // SEED-equivalent reactions via EC. We resolve candidate rxn IDs
            // once per suggestion.
            List<String> equivalentRxnIds = [rxn]
            if (rxnGraph != null) {
                // Find ECs for this reaction (direct or via gap metadata
                // — simplest is to iterate graph reactions and find any
                // with matching id; if that fails, also try EC-bridge via
                // metabolic gap table).
                Set<String> ecs = rxnGraph.ecsForReaction(rxn)
                if (ecs.isEmpty()) {
                    // Try matching by MetabolicGap
                    for (def gap : state.metabolicGaps) {
                        if (gap.reactionId == rxn && gap.ecNumber) {
                            ecs = Collections.singleton(gap.ecNumber)
                            break
                        }
                    }
                }
                for (String ec : ecs) {
                    for (String altRxn : rxnGraph.reactionsForEc(ec)) {
                        if (altRxn != rxn) equivalentRxnIds << altRxn
                    }
                }
            }

            // Collect (proteinId, baseScore) list from the suggestion.
            List<Tuple2<String, Double>> baseList = collectProteinScores(s)
            if (baseList.isEmpty()) { out << s; continue }

            // Compute LR-adjusted scores.
            List<Tuple2<String, Double>> adjusted = []
            boolean anyRescored = false
            for (Tuple2<String, Double> tup : baseList) {
                String pid = tup.v1
                double base = tup.v2
                String og = ortho[pid]
                double lr = 1.0d
                if (og != null) {
                    // Max LR across all equivalent reaction IDs.
                    ReactionLocusCatalog.Entry best = null
                    for (String r : equivalentRxnIds) {
                        ReactionLocusCatalog.Entry e = catalog.get(og, r)
                        if (e != null && e.nSigTotal >= minSupport) {
                            if (!requireCredible ||
                                    Math.abs(e.logLR()) > e.logLRCiHalfWidth()) {
                                if (best == null || e.lr() > best.lr()) best = e
                            }
                        }
                    }
                    if (best != null) {
                        lr = best.lr()
                        anyRescored = true
                    }
                }
                double adj = base * Math.pow(Math.max(lr, 1e-9), lambda)
                adjusted << new Tuple2(pid, adj)
            }
            if (!anyRescored) { out << s; continue }
            rescored++

            // Renormalise and pick singleton vs disjunctive.
            adjusted.sort { a, b -> Double.compare(b.v2, a.v2) }
            double sum = 0.0d
            for (Tuple2<String, Double> t : adjusted) sum += t.v2
            if (sum <= 0.0d) { out << s; continue }

            double qTop = adjusted[0].v2 / sum
            String topId = adjusted[0].v1

            if (qTop > singletonThreshold) {
                SingletonSuggestion ss = new SingletonSuggestion(
                    proteinId: topId,
                    q: qTop,
                    functionId: s.functionId,
                    functionType: s.functionType,
                    pathwayId: s.pathwayId,
                    reactionId: s.reactionId,
                    operonId: (s.operonId ?: 'rlgc') + '+cg',
                    bayesFactor: s.bayesFactor,
                    suggestionScore: Math.min(maxSuggestionScore, qTop),
                    proteinScores: s.proteinScores,
                    provenance: "${s.provenance} | CG rescore λ=${String.format(Locale.ROOT, '%.2f', lambda)}, qTop=${String.format(Locale.ROOT, '%.2f', qTop)}".toString(),
                )
                out << ss
            } else {
                // Disjunctive with new q values. Keep proteins until cumulative qsum ≥ 0.9.
                List<String> ids = []
                List<Double> qs = []
                double cum = 0.0d
                for (Tuple2<String, Double> t : adjusted) {
                    double q = t.v2 / sum
                    ids << t.v1
                    qs << q
                    cum += q
                    if (cum >= 0.9d) break
                }
                // Concentration score (same shape as RLGC's disjunctive branch)
                double H = 0.0d
                for (double q : qs) if (q > 0.0d) H -= q * Math.log(q)
                double uniform = Math.log((double) ids.size())
                double concentration = uniform > 0 ? (1.0d - H / uniform) : 0.0d
                DisjunctiveSuggestion ds = new DisjunctiveSuggestion(
                    proteinIds: ids,
                    qValues: qs,
                    functionId: s.functionId,
                    functionType: s.functionType,
                    pathwayId: s.pathwayId,
                    reactionId: s.reactionId,
                    operonId: (s.operonId ?: 'rlgc') + '+cg',
                    bayesFactor: s.bayesFactor,
                    suggestionScore: Math.min(maxSuggestionScore, concentration),
                    proteinScores: s.proteinScores,
                    provenance: "${s.provenance} | CG rescore λ=${String.format(Locale.ROOT, '%.2f', lambda)}, k=${ids.size()}, cov=${String.format(Locale.ROOT, '%.2f', cum)}".toString(),
                )
                out << ds
            }
        }
        integrated.suggestions = out
        log.info("CrossGenomeReScorer: rescored ${rescored} / ${out.size()} suggestions (λ=${lambda})")
        integrated
    }

    private static List<Tuple2<String, Double>> collectProteinScores(Suggestion s) {
        if (s instanceof SingletonSuggestion) {
            SingletonSuggestion ss = (SingletonSuggestion) s
            return [new Tuple2<String, Double>(ss.proteinId, (Double) ss.q)]
        }
        if (s instanceof DisjunctiveSuggestion) {
            DisjunctiveSuggestion ds = (DisjunctiveSuggestion) s
            List<Tuple2<String, Double>> out = []
            for (int i = 0; i < ds.proteinIds.size(); i++) {
                out << new Tuple2<String, Double>(ds.proteinIds[i], (Double) ds.qValues[i])
            }
            return out
        }
        []
    }
}
