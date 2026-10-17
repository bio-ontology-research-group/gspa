package gspa.integration.prior

import gspa.integration.IntegrationState
import gspa.integration.Prior
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Boosts candidate claims whose function is a (transitive) descendant of
 * an essential function that is currently uncovered by the genome.
 *
 * <p>Plan §A.4 of Phase 7. The essential function list comes from
 * {@link gspa.config.EssentialFunctions}; uncoveredness is assessed against
 * the genome's current propagated annotation set (MAP posterior &gt; 0.5).
 * For every candidate (protein, f) where {@code f} is a GO subclass of an
 * uncovered essential, the prior returns {@code +alphaEss} log-odds.</p>
 *
 * <p>Capped at {@code 2 * alphaEss} per claim so one claim cannot be
 * boosted arbitrarily by many uncovered essentials.</p>
 *
 * <p>If the state is missing {@code goReasoner} or
 * {@code essentialFunctions}, the prior no-ops (returns 0).</p>
 */
class EssentialityPrior implements Prior {

    private static final Logger log = LoggerFactory.getLogger(EssentialityPrior)

    /** Per-descendant boost magnitude (log-odds). Learned by the benchmark. */
    double alphaEss = 1.5d

    /** Maximum total boost per claim, as a multiple of alphaEss. */
    double maxBoostMultiplier = 2.0d

    /** Cached per-iteration: union of GO descendants of every uncovered essential function (fallback when claimCounts is empty). */
    private Set<String> boostableDescendants = Collections.emptySet()

    /** Cached per-iteration: functionKey -> log-odds adjustment when claimCounts is populated. */
    private Map<String, Double> keyBoosts = Collections.emptyMap()

    /** Cached per-iteration: number of boost sources, for per-claim capping. */
    private int uncoveredCount = 0

    @Override
    String name() { 'essentiality' }

    @Override
    Set<String> inputs() { ['essentialFunctions', 'goOntology', 'goReasoner', 'posteriors'] as Set }

    @Override
    void beginIteration(IntegrationState state) {
        boostableDescendants = Collections.emptySet()
        keyBoosts = Collections.emptyMap()
        uncoveredCount = 0

        if (state.essentialFunctions == null) return

        // Currently covered: use propagated set if goOntology or reasoner is available,
        // otherwise fall back to the raw annotated set.
        Set<String> annotated = (state.goOntology != null || state.goReasoner != null)
            ? state.baselineAnnotatedGoTermsPropagated()
            : state.baselineAnnotatedGoTerms()
        Set<String> essentials = state.essentialFunctions.getGoTerms()
        Set<String> uncovered = new LinkedHashSet<>()
        for (String e : essentials) {
            if (!annotated.contains(e)) uncovered << e
        }
        uncoveredCount = uncovered.size()

        Map<String, Set<String>> descByUncovered = new LinkedHashMap<>()
        Set<String> descendants = new LinkedHashSet<>()
        for (String u : uncovered) {
            Set<String> uSet = new LinkedHashSet<>()
            uSet << u
            if (state.goReasoner != null) {
                try {
                    uSet.addAll(state.goReasoner.getSubClasses(u, false))
                } catch (Exception ex) {
                    log.warn("EssentialityPrior: getSubClasses failed for ${u}: ${ex.message}")
                }
            } else if (state.goOntology != null) {
                try {
                    uSet.addAll(state.goOntology.getDescendants(u))
                } catch (Exception ex) {
                    log.warn("EssentialityPrior: getDescendants failed for ${u}: ${ex.message}")
                }
            }
            descByUncovered[u] = uSet
            descendants.addAll(uSet)
        }
        if (state.goReasoner == null && state.goOntology == null && !uncovered.isEmpty()) {
            log.info("EssentialityPrior: no GoOntology/GoReasoner, boosting ${uncovered.size()} exact essential terms (no descendant expansion)")
        }
        boostableDescendants = descendants

        if (state.claimCounts != null && !state.claimCounts.isEmpty()) {
            Set<String> allEssDomain = new LinkedHashSet<>(essentials)
            if (state.goReasoner != null) {
                for (String e : essentials) {
                    try { allEssDomain.addAll(state.goReasoner.getSubClasses(e, false)) } catch (Exception ignored) {}
                }
            } else if (state.goOntology != null) {
                for (String e : essentials) {
                    try { allEssDomain.addAll(state.goOntology.getDescendants(e)) } catch (Exception ignored) {}
                }
            }
            Map<String, Double> boosts = state.computeProteinSpecificityPrior(null, allEssDomain, 0.30d * alphaEss, 0.30d * alphaEss)
            for (Map.Entry<String, Set<String>> ue : descByUncovered.entrySet()) {
                Set<String> uSet = ue.value
                List<Map.Entry<String, Double>> cands = []
                for (String fkey : state.posteriorLogOdds.keySet()) {
                    String[] parts = IntegrationState.splitFunctionKey(fkey)
                    if (parts == null || parts[1] != 'GO') continue
                    if (!uSet.contains(parts[2])) continue
                    double likeProb = state.getLikelihoodProbability(fkey)
                    if (likeProb >= 0.16d && likeProb <= IntegrationState.MAP_THRESHOLD) {
                        cands << new AbstractMap.SimpleImmutableEntry<>(fkey, likeProb)
                    }
                }
                state.sortCandidatesForRescue(cands)
                int limit = Math.min(1, cands.size())
                for (int i = 0; i < limit; i++) {
                    String fk = cands[i].key
                    double bVal = Math.max(Math.min(maxBoostMultiplier * alphaEss, 1.30d * alphaEss), 0.14d - state.getLikelihoodLogOdds(fk))
                    boosts[fk] = Math.max(boosts.getOrDefault(fk, 0.0d), bVal)
                }
            }
            keyBoosts = boosts
        }
        log.debug("EssentialityPrior: ${uncovered.size()} uncovered essentials, ${descendants.size()} boostable descendants")
    }

    @Override
    double logOddsBoost(String proteinId, String functionKey, IntegrationState state) {
        if (state.claimCounts != null && !state.claimCounts.isEmpty()) {
            return keyBoosts.getOrDefault(functionKey, 0.0d)
        }
        if (boostableDescendants.isEmpty()) return 0.0d
        String[] parts = IntegrationState.splitFunctionKey(functionKey)
        if (parts == null || parts[1] != 'GO') return 0.0d
        if (!boostableDescendants.contains(parts[2])) return 0.0d

        return Math.min(maxBoostMultiplier * alphaEss, alphaEss)
    }
}
