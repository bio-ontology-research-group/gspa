package gspa.integration.prior

import gspa.integration.IntegrationState
import gspa.integration.MetabolicGap
import gspa.integration.Prior
import gspa.model.AnnotationType
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Boosts candidate claims whose function matches a gapseq-identified
 * metabolic gap.
 *
 * <p>Plan §A.4 of Phase 7. For each gap {@code (pathway P, reaction R)}
 * with EC number {@code E} (and mapped GO term {@code F}), any candidate
 * claim whose functionId equals {@code E} or {@code F} earns
 * {@code +alphaGap}. Gapseq's own gap-fill guesses (marked
 * {@code gapseqGuessed = true}) get a reduced {@code +alphaGap * 0.7}.</p>
 *
 * <p>The metabolic gap list is a pre-computed input on
 * {@link IntegrationState#metabolicGaps}; the caller (pipeline) populates
 * it after running {@code GapseqPredictor}.</p>
 */
class GapFillingPrior implements Prior {

    private static final Logger log = LoggerFactory.getLogger(GapFillingPrior)

    /** Per-match boost magnitude (log-odds). Learned by the benchmark. */
    double alphaGap = 1.5d

    /** Discount for gapseq's own gap-fill guesses. */
    double gapseqGuessedFactor = 0.7d

    private static final Set<String> GENERIC_CATALYTIC_ROOTS = [
        'GO:0003824', 'GO:0016491', 'GO:0016740', 'GO:0016787',
        'GO:0016829', 'GO:0016853', 'GO:0016874'
    ] as Set

    /** Cached per iteration: function_id → strongest boost magnitude (fallback when claimCounts is empty). */
    private Map<String, Double> boostByFunctionId = new LinkedHashMap<>()

    /** Cached per iteration: functionKey -> targeted log-odds adjustment when claimCounts is populated. */
    private Map<String, Double> targetedKeyBoosts = Collections.emptyMap()

    @Override
    String name() { 'gap_filling' }

    @Override
    Set<String> inputs() { ['metabolicGaps'] as Set }

    @Override
    void beginIteration(IntegrationState state) {
        boostByFunctionId = new LinkedHashMap<>()
        targetedKeyBoosts = Collections.emptyMap()
        if (state.metabolicGaps == null || state.metabolicGaps.isEmpty()) return

        Set<String> allGapGoTerms = new LinkedHashSet<>()
        Map<String, List<Map.Entry<String, Double>>> candsByFunc = null
        if (state.claimCounts != null && !state.claimCounts.isEmpty()) {
            candsByFunc = [:].withDefault { [] }
            for (String fkey : state.posteriorLogOdds.keySet()) {
                String[] parts = IntegrationState.splitFunctionKey(fkey)
                if (parts == null) continue
                if (parts[1] != AnnotationType.EC.name() && parts[1] != AnnotationType.GO.name()) continue
                double likeProb = state.getLikelihoodProbability(fkey)
                if (likeProb >= 0.18d && likeProb <= IntegrationState.MAP_THRESHOLD) {
                    candsByFunc[parts[2]] << new AbstractMap.SimpleImmutableEntry<>(fkey, likeProb)
                }
            }
        }

        Map<String, Double> topBoosts = new LinkedHashMap<>()
        for (MetabolicGap gap : state.metabolicGaps) {
            double mag = gap.gapseqGuessed ? (alphaGap * gapseqGuessedFactor) : alphaGap
            Set<String> gapGoTerms = new LinkedHashSet<>()
            Set<String> gapMatchFuncs = new LinkedHashSet<>()
            if (gap.ecNumber) {
                for (String rawEc : gap.ecNumber.split('[/,;\\s]+')) {
                    String ec = rawEc.trim()
                    if (!ec) continue
                    updateMax(ec, mag)
                    gapMatchFuncs << ec
                    String mappedGo = state.pathwayDatabase?.ec2go?.get(ec)
                    if (mappedGo) gapGoTerms << mappedGo
                }
            }
            if (gap.goTerm) {
                for (String rawGo : gap.goTerm.split('[,;\\s]+')) {
                    String go = rawGo.trim()
                    if (go) gapGoTerms << go
                }
            }
            for (String go : gapGoTerms) {
                updateMax(go, mag)
                gapMatchFuncs << go
                allGapGoTerms << go
                if (state.goOntology != null && !GENERIC_CATALYTIC_ROOTS.contains(go)) {
                    for (String desc : state.goOntology.getDescendants(go)) {
                        updateMax(desc, mag)
                        gapMatchFuncs << desc
                        allGapGoTerms << desc
                    }
                }
            }
            if (candsByFunc != null && !gapMatchFuncs.isEmpty()) {
                List<Map.Entry<String, Double>> cands = []
                for (String fid : gapMatchFuncs) {
                    List<Map.Entry<String, Double>> list = candsByFunc.get(fid)
                    if (list != null) cands.addAll(list)
                }
                if (!cands.isEmpty()) {
                    state.sortCandidatesForRescue(cands)
                    String fk = cands[0].key
                    double pVal = cands[0].value as double
                    if (state.candidateRescueTier(fk, pVal) >= 1) {
                        double bVal = Math.max(1.15d * mag, 0.14d - state.getLikelihoodLogOdds(fk))
                        topBoosts[fk] = Math.max(topBoosts.getOrDefault(fk, 0.0d), bVal)
                    }
                }
            }
        }
        if (state.claimCounts != null && !state.claimCounts.isEmpty()) {
            Map<String, Double> boosts = state.computeProteinSpecificityPrior(null, allGapGoTerms, 0.38d * alphaGap, 0.35d * alphaGap)
            for (Map.Entry<String, Double> te : topBoosts.entrySet()) {
                boosts[te.key] = Math.max(boosts.getOrDefault(te.key, 0.0d), te.value)
            }
            targetedKeyBoosts = boosts
        }
        log.debug("GapFillingPrior: ${boostByFunctionId.size()} functions associated with metabolic gaps")
    }

    @Override
    double logOddsBoost(String proteinId, String functionKey, IntegrationState state) {
        if (state.claimCounts != null && !state.claimCounts.isEmpty()) {
            return targetedKeyBoosts.getOrDefault(functionKey, 0.0d)
        }
        if (boostByFunctionId.isEmpty()) return 0.0d
        String[] parts = IntegrationState.splitFunctionKey(functionKey)
        if (parts == null) return 0.0d
        String type = parts[1]
        String funcId = parts[2]
        if (type != AnnotationType.EC.name() && type != AnnotationType.GO.name()) return 0.0d
        return boostByFunctionId.getOrDefault(funcId, 0.0d)
    }

    private void updateMax(String key, double value) {
        boostByFunctionId.merge(key, value, { a, b -> Math.max(a as double, b as double) })
    }
}
