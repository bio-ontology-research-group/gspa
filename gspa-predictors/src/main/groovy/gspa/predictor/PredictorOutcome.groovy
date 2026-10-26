package gspa.predictor

/**
 * Terminal state of one predictor in one pipeline run. The distinction between
 * {@link #SUCCESS_ZERO} and {@link #FAILED}/{@link #UNAVAILABLE} matters
 * biologically: "ran and found nothing" is evidence of absence, whereas a
 * missing or crashed tool says nothing about the proteins.
 */
enum PredictorOutcome {
    /** Tool/database not installed or not configured; it never ran. */
    UNAVAILABLE,
    /** Tool started but crashed, timed out, threw, or produced unusable output. */
    FAILED,
    /** Tool completed normally and reported no annotations. */
    SUCCESS_ZERO,
    /** Tool completed normally and reported at least one annotation. */
    SUCCESS

    boolean isSuccessful() { this == SUCCESS || this == SUCCESS_ZERO }
}
