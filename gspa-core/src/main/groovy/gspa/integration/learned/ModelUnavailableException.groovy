package gspa.integration.learned

/**
 * No usable trained artifact (absent, wrong schema/contract, not frozen, mechanism-test
 * fixture in a production load, or an optional component such as the null-outcome layer
 * that was never trained). Callers must fail or skip the learned stage; there are no
 * built-in default weights, boosts or null probabilities to fall back on.
 */
class ModelUnavailableException extends IllegalStateException {
    ModelUnavailableException(String message) { super(message) }
}
