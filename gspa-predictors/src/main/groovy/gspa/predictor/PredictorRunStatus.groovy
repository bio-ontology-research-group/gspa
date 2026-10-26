package gspa.predictor

import groovy.transform.Canonical

/**
 * Structured completion record for one predictor run. Serialised into the
 * machine-readable run manifest by {@link PredictorRunManifest}.
 */
@Canonical
class PredictorRunStatus {
    String predictor
    PredictorOutcome outcome
    /** Annotations returned by the predictor (before assignment to proteins). */
    int annotationCount = 0
    /** Human-readable reason; set for UNAVAILABLE / FAILED. */
    String detail
    /** Process exit code when an external tool was executed, else null. */
    Integer exitCode
    long elapsedMillis = 0L

    static PredictorRunStatus unavailable(String predictor, String detail) {
        new PredictorRunStatus(predictor: predictor, outcome: PredictorOutcome.UNAVAILABLE, detail: detail)
    }

    static PredictorRunStatus failed(String predictor, String detail, Integer exitCode = null) {
        new PredictorRunStatus(predictor: predictor, outcome: PredictorOutcome.FAILED, detail: detail, exitCode: exitCode)
    }

    /** SUCCESS or SUCCESS_ZERO depending on the number of annotations. */
    static PredictorRunStatus completed(String predictor, int annotationCount, Integer exitCode = null) {
        new PredictorRunStatus(predictor: predictor,
            outcome: annotationCount > 0 ? PredictorOutcome.SUCCESS : PredictorOutcome.SUCCESS_ZERO,
            annotationCount: annotationCount, exitCode: exitCode)
    }

    Map<String, Object> toMap() {
        [predictor: predictor, outcome: outcome.name(), annotation_count: annotationCount,
         detail: detail, exit_code: exitCode, elapsed_ms: elapsedMillis] as Map<String, Object>
    }
}
