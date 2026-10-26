package gspa.integration.learned

/**
 * Conditional locus-ranking distribution for one reaction-locus group: the probability that
 * each supplied candidate is THE responsible locus, or that none of them is (the explicit
 * null outcome). Probabilities sum to one within the group and nowhere else.
 *
 * <p>This is a different quantity from {@link CandidateScore#probability}: that one is the
 * calibrated marginal correctness of a candidate function and may be high for several
 * candidates at once (paralogues, isoenzymes, multi-reaction proteins). A singleton group is
 * not normalised to 1: its null probability stays whatever the trained null layer assigns.
 */
class LocusAssignment {
    /** candidate key -> P(this candidate is the responsible locus | exactly one of these or none) */
    final Map<String, Double> candidateProbabilities
    /** P(no supplied candidate is correct). */
    final double nullProbability
    final String modelId

    LocusAssignment(Map<String, Double> candidateProbabilities, double nullProbability, String modelId) {
        this.candidateProbabilities = Collections.unmodifiableMap(candidateProbabilities)
        this.nullProbability = nullProbability
        this.modelId = modelId
    }

    String toString() { "LocusAssignment(null=${nullProbability}, ${candidateProbabilities})" }
}
