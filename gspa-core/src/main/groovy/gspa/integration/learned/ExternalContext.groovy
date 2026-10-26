package gspa.integration.learned

/**
 * Optional, structured context the integration state itself does not carry. Everything
 * here is supplied by the caller from a real analysis (reaction-graph neighbourhood,
 * cohort recurrence, FBA, predictor run manifest); the extractor never fabricates it.
 * Anything absent is reported to the model as unavailable (indicator 0), not as a
 * measured zero.
 *
 * <p>All maps are keyed by the legacy function key {@code proteinId|TYPE|functionId}
 * (see {@link gspa.integration.IntegrationState#functionKey}); they apply to every role of
 * that function.
 */
class ExternalContext {

    /**
     * Predictor completion outcomes by predictor name (PredictorOutcome names:
     * UNAVAILABLE, FAILED, PARTIAL, SUCCESS_ZERO, SUCCESS). Null = not supplied.
     */
    Map<String, String> predictorStatus = null

    /** Neighbourhood features computed from genome context + reaction graph. */
    Map<String, Neighborhood> neighborhood = [:]

    /** Cross-genome recurrence from a cohort analysis. */
    Map<String, Cohort> cohort = [:]

    /** Conditional FBA results. */
    Map<String, FbaConditional> fba = [:]

    /**
     * Explicit roles from structured evidence (e.g. a transporter annotation). When
     * absent, roles are derived conservatively by {@link ContextFeatureExtractor}.
     */
    Map<String, List<CandidateRole>> roles = [:]

    /** External assembly completeness estimate in percent (e.g. CheckM2), or null. */
    Double assemblyCompleteness = null

    static class Neighborhood {
        double supportMax
        int nSupporting
        double conflictMax
        /** Reaction-graph distance to the nearest supported neighbour reaction; null = unavailable. */
        Double reactionDistance
    }

    static class Cohort {
        int nGenomes
        int nWithOrthogroup
        int nWithFunctionSeq
        int nWithContext
        int nIncomplete
    }

    static class FbaConditional {
        int nConditions
        double requiredFraction
        /** Relative growth change on knockout; null = not computed. */
        Double knockoutGrowthEffect
        /** True when the tested objective itself required the candidate (not independent evidence). */
        boolean objectiveRequiresCandidate
    }
}
