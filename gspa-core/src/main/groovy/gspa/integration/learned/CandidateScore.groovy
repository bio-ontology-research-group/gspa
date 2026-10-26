package gspa.integration.learned

/**
 * Learned per-candidate result. {@link #probability} is the calibrated probability that
 * THIS candidate's function is correct for THIS protein in THIS role: an independent
 * binary quantity. It is not normalised across proteins, functions or roles that can
 * biologically coexist, and it is not a locus-ranking probability (see {@link LocusAssignment}).
 *
 * <p>{@link #aboveDisplayThreshold} is advisory only: scoring returns every candidate
 * regardless of it, and it never gates retention.
 */
class CandidateScore {
    final CandidateFeatures features
    /** Pre-calibration logit of the frozen model. */
    final double margin
    /** Calibrated probability of correctness. */
    final double probability
    final double displayThreshold
    final boolean aboveDisplayThreshold
    final String modelId

    CandidateScore(CandidateFeatures features, double margin, double probability,
                   double displayThreshold, String modelId) {
        this.features = features
        this.margin = margin
        this.probability = probability
        this.displayThreshold = displayThreshold
        this.aboveDisplayThreshold = probability >= displayThreshold
        this.modelId = modelId
    }

    String candidateKey() { features.candidateKey() }

    @Override
    String toString() { String.format(Locale.ROOT, 'CandidateScore(%s p=%.4f)', candidateKey(), probability) }
}
