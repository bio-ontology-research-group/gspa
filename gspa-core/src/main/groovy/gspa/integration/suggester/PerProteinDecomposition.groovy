package gspa.integration.suggester

import groovy.transform.Canonical
import groovy.transform.builder.Builder
import groovy.transform.builder.SimpleStrategy

/**
 * Per-protein decomposition of the Phase 8 suggester's posterior
 * calculation for a single (gap function f_R, candidate protein p) pair.
 * Exposes the individual terms so users can audit how a suggestion was
 * generated.
 *
 * <p>Layer 2 of plan §A.1: {@code L_R(p) = L_lik + L_op + L_lm}, then
 * {@code π_R(p) = sigmoid(L_R)}, then Layer 3 softmax to {@code q}.</p>
 */
@Canonical
@Builder(builderStrategy = SimpleStrategy, prefix = '')
class PerProteinDecomposition {

    /** Protein ID. */
    String proteinId

    /**
     * Existing-predictor log-odds for {@code (p, f_R)}. Taken from the
     * integrator's converged posterior; defaults to a low "absent" value
     * (around -6 log-odds, i.e., probability ~0.0025) when no claim
     * exists, so fully-dark proteins aren't given 50/50 priors.
     */
    double likelihoodLogOdds

    /**
     * Operon-membership contribution: {@code log(BF / (1 + BF))} where
     * BF is the Bayes factor that the operon participates in the gap's
     * pathway. Same value for every operon member.
     */
    double operonLogOdds

    /**
     * Genomic-LM contribution (Phase 9). Set to 0 when no LM is wired.
     */
    double lmLogOdds

    /** Total log-odds = lik + op + lm. */
    double totalLogOdds

    /** sigmoid(totalLogOdds). */
    double piR

    /**
     * Softmax over operon members: {@code q = π_R(p) / Σ π_R(p')}.
     * Interpretable as "probability that p is the one that performs R,
     * given the gap constraint".
     */
    double q
}
