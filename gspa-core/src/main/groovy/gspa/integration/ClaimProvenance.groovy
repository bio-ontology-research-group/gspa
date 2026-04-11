package gspa.integration

import groovy.transform.Canonical
import groovy.transform.builder.Builder
import groovy.transform.builder.SimpleStrategy

/**
 * Per-(protein, function) provenance produced by the integrator.
 *
 * Records the supporting claims, the likelihood contribution from the
 * combiner, any prior contributions (filled in by {@code PriorEngine} in
 * Phase 7.3), and the final posterior log-odds / probability.
 */
@Canonical
@Builder(builderStrategy = SimpleStrategy, prefix = '')
class ClaimProvenance {

    /** Function key, e.g. "proteinId|GO|GO:0006412". */
    String functionKey

    /** Protein ID extracted from the function key, stored for fast lookup. */
    String proteinId

    /** Claims that fed into this posterior. */
    List<EvidenceClaim> supportingClaims = []

    /**
     * Per-prior log-odds contribution. Populated in Phase 7.3; empty in
     * Phase 7.1/7.2 where priors are absent.
     */
    Map<String, Double> priorContributions = [:]

    /** Log-odds of the likelihood combination (Noisy-OR output). */
    double likelihoodLogOdds = 0.0

    /** Final log-odds including prior contributions. */
    double finalLogOdds = 0.0

    /** Final posterior probability in [0, 1]. */
    double finalProbability = 0.0

    /** Iteration at which convergence was reached. */
    int convergenceIter = 0
}
