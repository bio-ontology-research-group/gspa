package gspa.integration.learned

/**
 * The learned stage rejected its input or a computed value. Thrown only before the integration
 * state and annotation set are touched (or after they were restored): a failed stage leaves them
 * exactly as they were.
 */
class LearnedStageException extends IllegalStateException {
    LearnedStageException(String message) { super(message) }
    LearnedStageException(String message, Throwable cause) { super(message, cause) }
}
