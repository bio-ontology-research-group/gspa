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

    /**
     * Cached per-iteration: union of GO descendants of every uncovered
     * essential function. If a candidate's functionId is in this set, it
     * earns a boost.
     */
    private Set<String> boostableDescendants = Collections.emptySet()

    /** Cached per-iteration: number of boost sources, for per-claim capping. */
    private int uncoveredCount = 0

    @Override
    String name() { 'essentiality' }

    @Override
    Set<String> inputs() { ['essentialFunctions', 'goOntology', 'goReasoner', 'posteriors'] as Set }

    @Override
    void beginIteration(IntegrationState state) {
        boostableDescendants = Collections.emptySet()
        uncoveredCount = 0

        if (state.essentialFunctions == null) return
        if (state.goReasoner == null) {
            log.debug('EssentialityPrior disabled: no GoReasoner in state')
            return
        }

        // Currently covered: propagated up through the GO hierarchy.
        Set<String> annotatedPropagated = state.currentlyAnnotatedGoTermsPropagated()
        Set<String> essentials = state.essentialFunctions.getGoTerms()
        Set<String> uncovered = new LinkedHashSet<>()
        for (String e : essentials) {
            if (!annotatedPropagated.contains(e)) uncovered << e
        }
        uncoveredCount = uncovered.size()

        Set<String> descendants = new LinkedHashSet<>()
        for (String u : uncovered) {
            descendants << u
            try {
                descendants.addAll(state.goReasoner.getSubClasses(u, false))
            } catch (Exception ex) {
                log.warn("EssentialityPrior: getSubClasses failed for ${u}: ${ex.message}")
            }
        }
        boostableDescendants = descendants
        log.debug("EssentialityPrior: ${uncovered.size()} uncovered essentials, ${descendants.size()} boostable descendants")
    }

    @Override
    double logOddsBoost(String proteinId, String functionKey, IntegrationState state) {
        if (boostableDescendants.isEmpty()) return 0.0d
        String[] parts = IntegrationState.splitFunctionKey(functionKey)
        if (parts == null || parts[1] != 'GO') return 0.0d
        if (!boostableDescendants.contains(parts[2])) return 0.0d

        // Cap per-claim boost regardless of how many essentials the
        // descendant could be a subclass of. We always emit +alphaEss here;
        // the cap is enforced at the PriorEngine level indirectly because
        // each claim participates in at most one "descendants" match.
        return Math.min(maxBoostMultiplier * alphaEss, alphaEss)
    }
}
