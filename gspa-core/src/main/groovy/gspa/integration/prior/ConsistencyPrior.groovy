package gspa.integration.prior

import gspa.integration.IntegrationState
import gspa.integration.Prior
import gspa.model.ConsistencyViolation
import gspa.ontology.ConsistencyResult
import gspa.ontology.SatConsistencyChecker
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Downweights candidate claims that participate in a taxon-constraint
 * violation, as identified by
 * {@link SatConsistencyChecker#check(Set)}'s UNSAT core.
 *
 * <p>Soft downweight by default: a claim whose GO term is flagged returns
 * {@code -alphaCons} log-odds, which subtracts from the likelihood but
 * does not remove the claim outright. Strong likelihoods can survive one
 * soft penalty — the right behaviour for HGT and contamination cases.
 * Hard-filter mode is available behind the config flag
 * {@code integration.consistency.hardFilter: true}, in which case
 * {@code alphaCons = 1000} is used (effectively removing the claim).</p>
 *
 * <p>The SAT check runs once per refinement iteration over the currently
 * annotated set (posterior > 0.5, the MAP decision), so the per-claim
 * lookup in {@link #logOddsBoost} is a cheap set query.</p>
 */
class ConsistencyPrior implements Prior {

    private static final Logger log = LoggerFactory.getLogger(ConsistencyPrior)

    /** Soft penalty in log-odds. Learned by the benchmark. */
    double alphaCons = 3.0d

    /** If true, effectively remove the claim entirely. */
    boolean hardFilter = false

    /** Cached per iteration: GO terms that appear in any SAT UNSAT core. */
    private Set<String> conflictingTerms = Collections.emptySet()

    @Override
    String name() { 'consistency' }

    @Override
    Set<String> inputs() { ['satConsistencyChecker', 'posteriors'] as Set }

    @Override
    void beginIteration(IntegrationState state) {
        conflictingTerms = Collections.emptySet()
        SatConsistencyChecker checker = state.satConsistencyChecker
        if (checker == null) {
            log.debug('ConsistencyPrior disabled: no SatConsistencyChecker in state')
            return
        }

        Set<String> annotated = state.currentlyAnnotatedGoTerms()
        if (annotated.isEmpty()) return

        try {
            ConsistencyResult result = checker.check(annotated)
            if (result.consistent) return
            Set<String> terms = new LinkedHashSet<>()
            for (ConsistencyViolation v : result.violations ?: []) {
                if (v.involvedGoTerms != null) terms.addAll(v.involvedGoTerms)
            }
            conflictingTerms = terms
            log.debug("ConsistencyPrior: ${terms.size()} conflicting GO terms from ${result.violations?.size() ?: 0} violations")
        } catch (Exception ex) {
            log.warn("ConsistencyPrior: SAT check failed: ${ex.message}")
        }
    }

    @Override
    double logOddsBoost(String proteinId, String functionKey, IntegrationState state) {
        if (conflictingTerms.isEmpty()) return 0.0d
        String[] parts = IntegrationState.splitFunctionKey(functionKey)
        if (parts == null || parts[1] != 'GO') return 0.0d
        if (!conflictingTerms.contains(parts[2])) return 0.0d
        double penalty = hardFilter ? -1000.0d : -alphaCons
        return penalty
    }
}
