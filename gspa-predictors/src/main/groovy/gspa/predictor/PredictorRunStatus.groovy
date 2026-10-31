package gspa.predictor

import groovy.transform.Canonical

/**
 * Structured completion record for one predictor run. Serialised into the
 * machine-readable run manifest by {@link PredictorRunManifest}.
 */
@Canonical
class PredictorRunStatus {

    /** At most this many skipped-input ids are listed in the manifest. */
    static final int MAX_LISTED_SKIPPED = 50

    String predictor
    PredictorOutcome outcome
    /** Annotations returned by the predictor (before assignment to proteins). */
    int annotationCount = 0
    /**
     * Human-readable reason for UNAVAILABLE / FAILED / PARTIAL; for SUCCESS and
     * SUCCESS_ZERO an optional interpretation caveat (e.g. what an empty result
     * does and does not mean for this tool).
     */
    String detail
    /** Process exit code when an external tool was executed, else null. */
    Integer exitCode
    long elapsedMillis = 0L
    /** Ids of inputs (contigs, proteins) the tool did not cover; non-empty implies PARTIAL. */
    List<String> skippedInputs = []

    static PredictorRunStatus unavailable(String predictor, String detail) {
        new PredictorRunStatus(predictor: predictor, outcome: PredictorOutcome.UNAVAILABLE, detail: detail)
    }

    static PredictorRunStatus skipped(String predictor, String detail) {
        new PredictorRunStatus(predictor: predictor, outcome: PredictorOutcome.SKIPPED, detail: detail)
    }

    static PredictorRunStatus failed(String predictor, String detail, Integer exitCode = null) {
        new PredictorRunStatus(predictor: predictor, outcome: PredictorOutcome.FAILED, detail: detail, exitCode: exitCode)
    }

    /** SUCCESS or SUCCESS_ZERO depending on the number of annotations. */
    static PredictorRunStatus completed(String predictor, int annotationCount, Integer exitCode = null,
                                        String note = null) {
        new PredictorRunStatus(predictor: predictor,
            outcome: annotationCount > 0 ? PredictorOutcome.SUCCESS : PredictorOutcome.SUCCESS_ZERO,
            annotationCount: annotationCount, exitCode: exitCode, detail: note)
    }

    /** Completed but not on all input: PARTIAL, whatever the annotation count. */
    static PredictorRunStatus partial(String predictor, int annotationCount, String detail,
                                      List<String> skipped, Integer exitCode = null) {
        new PredictorRunStatus(predictor: predictor, outcome: PredictorOutcome.PARTIAL,
            annotationCount: annotationCount, detail: detail, exitCode: exitCode, skippedInputs: skipped ?: [])
    }

    Map<String, Object> toMap() {
        [predictor: predictor, outcome: outcome.name(), annotation_count: annotationCount,
         detail: detail, exit_code: exitCode, elapsed_ms: elapsedMillis,
         skipped_input_count: skippedInputs.size(),
         skipped_inputs: skippedInputs.take(MAX_LISTED_SKIPPED)] as Map<String, Object>
    }
}
