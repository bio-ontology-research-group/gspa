package gspa.predictor

/**
 * Implemented by predictors that can report a typed completion status for
 * their most recent run, so the pipeline can tell a failed run apart from a
 * run that legitimately found nothing.
 */
interface StatusReportingPredictor {

    /** Status of the last predictBatch/predictGenome call; null before the first run. */
    PredictorRunStatus getLastRunStatus()

    /** Why {@code isAvailable()} is false (null when available). */
    String unavailableReason()
}
