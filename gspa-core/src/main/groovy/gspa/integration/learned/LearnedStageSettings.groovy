package gspa.integration.learned

import gspa.model.AnnotationType

/**
 * Typed settings of {@link LearnedIntegrationStage}. There are no probability floors, boosts or
 * count-dependent caps here: the only numeric knob is the optional retraction threshold, and it
 * is off unless set.
 */
class LearnedStageSettings {

    /** The frozen model. Required: without it the stage fails ({@link ModelUnavailableException}); nothing is defaulted. */
    LearnedContextModel model

    /**
     * Remove a candidate's legacy annotation when its learned probability is below this value
     * (0..1). Null (default) never retracts: scores are replaced but every candidate stays. The
     * removal is always recorded in the report with the score that caused it.
     */
    Double retractBelow = null

    /**
     * Add candidates that exist only in structured context ({@link ExternalContext#candidateKeys}),
     * with no raw claim, scoring them with no sequence evidence. They are weak by construction and
     * are not filtered by score.
     */
    boolean addContextCandidates = true

    /**
     * Add claimed functions that the preceding stages did not retain (for example because they were
     * demoted). Off by default: the stage rescores what the integration retained; turn it on when
     * the learned stage replaces the evidence-combination stage and starts from an empty set.
     */
    boolean addUnretainedClaims = false

    /** Compute the conditional locus ranking with explicit null where the model has a trained locus layer. */
    boolean computeLocusAssignments = true

    /** Keep the full feature vector of every candidate in the report (large); the exact-bits digest is always kept. */
    boolean retainFeatureVectors = false

    /** Restrict rescoring to these function types; null = all types. Other keys are left untouched. */
    Set<AnnotationType> types = null

    void validate() {
        if (model == null) {
            throw new ModelUnavailableException('no trained model supplied; the learned stage cannot run (no defaults exist)')
        }
        if (retractBelow != null && (!Double.isFinite(retractBelow) || retractBelow < 0d || retractBelow > 1d)) {
            throw new IllegalArgumentException("retractBelow must be in [0,1], got ${retractBelow}")
        }
        if (types != null && types.isEmpty()) throw new IllegalArgumentException('types must be null (all) or non-empty')
    }

    boolean inScope(AnnotationType t) { types == null || types.contains(t) }
}
