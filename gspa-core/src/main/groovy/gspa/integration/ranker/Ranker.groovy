package gspa.integration.ranker

/**
 * Phase 12 M3+: per-candidate ranking interface.
 *
 * <p>Implementations score a single (protein, reaction) pair given a
 * flat feature vector. The ranker is agnostic to feature construction;
 * {@link RankerFeatures} owns that.</p>
 *
 * <p>M3's {@link GbdtRanker} loads a LightGBM text model and scores
 * deterministically. M5's {@code NeuralRanker} (later) implements the
 * same interface via ONNX runtime.</p>
 */
interface Ranker {
    /**
     * @param featureVector in the same index order emitted by
     *        {@link RankerFeatures#featureNames()}.
     * @return a scalar score (higher = better candidate)
     */
    double score(double[] featureVector)

    /** Number of features expected. */
    int featureDim()

    /** Human-readable feature names, in index order. */
    List<String> featureNames()
}
