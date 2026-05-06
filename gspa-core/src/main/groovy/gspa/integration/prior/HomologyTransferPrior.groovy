package gspa.integration.prior

import gspa.integration.IntegrationState
import gspa.integration.Prior
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Phase 11 cross-genome homology-transfer prior.
 *
 * <p>For a protein {@code p} in orthogroup {@code O}, boosts the
 * log-odds of claim {@code (p, f)} toward the cluster's consensus
 * posterior probability of {@code f} across all orthologs — capped,
 * and only fires when (a) the consensus is strong
 * (≥ {@link #minConsensus}) and (b) the current claim is materially
 * weaker than the consensus (gap ≥ {@link #minDelta}).</p>
 *
 * <p>Upstream requirements on {@link IntegrationState}:</p>
 * <ul>
 *   <li>{@code orthogroupMap}: protein-id → cluster-id</li>
 *   <li>{@code orthogroupConsensus}: "clusterId|type|functionId" →
 *       P(function present in the cluster), computed from a prior
 *       pass that integrates all genomes separately.</li>
 * </ul>
 * Both must be non-null for the prior to fire; otherwise it returns 0
 * log-odds, acting as a no-op. This lets the same integrator code run
 * in single-genome mode (no orthogroups wired) and cross-genome mode.</p>
 *
 * <p>The prior runs in the main {@link gspa.integration.IterativeRefiner}
 * inner loop (once per refinement iteration). It is deliberately
 * asymmetric: consensus &gt; current → positive boost; consensus &lt;
 * current → zero. We don't DOWN-weight claims against the consensus
 * because a well-supported individual claim may be a true lineage-
 * specific annotation that looks like an outlier in the cluster; we
 * do not want homology transfer to erase valid evidence.</p>
 */
class HomologyTransferPrior implements Prior {

    private static final Logger log = LoggerFactory.getLogger(HomologyTransferPrior)

    /** Max log-odds the prior can add per claim. Learned by benchmark. */
    double alpha = 1.5d

    /** Minimum consensus probability (0..1) to fire. Below this, no boost. */
    double minConsensus = 0.5d

    /** Minimum gap (consensus - current_prob) required before firing. Avoids no-op or tiny pushes. */
    double minDelta = 0.05d

    /** Absolute cap on the log-odds boost per claim. */
    double maxBoost = 3.0d

    @Override String name() { 'homology_transfer' }

    @Override Set<String> inputs() { ['orthogroups', 'orthogroupConsensus', 'posteriors'] as Set }

    @Override
    void beginIteration(IntegrationState state) {
        // Consensus is pre-computed upstream; nothing to recompute.
    }

    @Override
    double logOddsBoost(String proteinId, String functionKey, IntegrationState state) {
        if (state.orthogroupMap == null || state.orthogroupConsensus == null) return 0.0d

        String clusterId = state.orthogroupMap[proteinId]
        if (clusterId == null) return 0.0d

        // functionKey is "proteinId|type|functionId"; we need just "type|functionId" for the consensus key
        String[] parts = IntegrationState.splitFunctionKey(functionKey)
        if (parts == null) return 0.0d
        String consensusKey = "${clusterId}|${parts[1]}|${parts[2]}"
        Double consensus = state.orthogroupConsensus[consensusKey]
        if (consensus == null || consensus < minConsensus) return 0.0d

        double currentProb = state.getProbability(functionKey)
        if (consensus <= currentProb + minDelta) return 0.0d

        // Transfer in log-odds space: push toward the consensus, damped by alpha.
        double consensusLogOdds = logit(consensus)
        double currentLogOdds = state.get(functionKey)
        double diff = consensusLogOdds - currentLogOdds
        if (diff <= 0.0d) return 0.0d

        double boost = alpha * diff
        return Math.max(0.0d, Math.min(maxBoost, boost))
    }

    private static double logit(double p) {
        double clipped = Math.min(1.0d - 1.0e-9d, Math.max(1.0e-9d, p))
        Math.log(clipped / (1.0d - clipped))
    }
}
