package gspa.predictor

/**
 * Terminal state of one predictor in one pipeline run. The distinction between
 * {@link #SUCCESS_ZERO} and {@link #FAILED}/{@link #UNAVAILABLE}/{@link #PARTIAL}
 * matters biologically: "ran on all input and found nothing" is evidence of
 * absence, whereas a missing, crashed or partially fed tool says nothing about
 * the proteins/contigs it did not cover.
 */
enum PredictorOutcome {
    /** Requested tool was not executed because this input has no eligible records. */
    SKIPPED,
    /** Tool/database not installed or not configured; it never ran. */
    UNAVAILABLE,
    /** Tool started but crashed, timed out, threw, or produced output that violates its schema. */
    FAILED,
    /**
     * Tool completed, but some input was not analysed (e.g. contigs without
     * sequence, CDS with invalid coordinates, rows outside the contig). Results
     * are valid for the covered input only; absence is not evidence.
     */
    PARTIAL,
    /** Tool completed normally on all input and reported no annotations. */
    SUCCESS_ZERO,
    /** Tool completed normally on all input and reported at least one annotation. */
    SUCCESS

    /** True only for runs whose absence of a call is meaningful (full coverage). */
    boolean isSuccessful() { this == SUCCESS || this == SUCCESS_ZERO }
}
