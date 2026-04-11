package gspa.integration

/**
 * Per-source calibration of raw predictor scores to calibrated probabilities.
 *
 * Predictors emit scores on their own scales (DIAMOND pident / 100, HMMER
 * bit score / 100, eggNOG −log10(evalue) / 50, etc.) which do not directly
 * correspond to "probability that the annotation is correct". The combiner
 * needs probabilities to feed Noisy-OR; Platt-style logistic calibration
 * maps raw → calibrated:
 *
 * <pre>
 *   calibrated = sigmoid( a * raw + b )
 * </pre>
 *
 * Per-source coefficients {@code (a, b)} can be registered here; the default
 * is identity calibration ({@code a = 4}, {@code b = -2}) which gives
 * {@code calibrated = sigmoid(4 * raw - 2)} — a sigmoid centred at
 * {@code raw = 0.5} with gentle slope, preventing log-odds infinities at
 * the endpoints and preserving the monotone ordering of raw scores.
 *
 * Coefficients are loaded from defaults at construction time. Phase 7.5
 * will overwrite them with values learned by the benchmark.
 */
class CalibrationTable {

    static final double EPSILON = 1.0e-4

    /** Platt coefficients indexed by predictor source name. */
    final Map<String, double[]> coefficients = new HashMap<>()

    /** Fallback coefficients if a source is not registered. */
    double[] defaultCoefficients = [4.0d, -2.0d] as double[]

    CalibrationTable() {
        loadDefaults()
    }

    /** Register a calibration curve for a predictor source. */
    void register(String source, double a, double b) {
        coefficients[source] = [a, b] as double[]
    }

    /**
     * Map a raw score to a calibrated probability. Clamped away from
     * {0, 1} by {@link #EPSILON} so the downstream log-odds conversion
     * never produces infinities.
     */
    double calibrate(String source, double rawScore) {
        double[] ab = coefficients.getOrDefault(source, defaultCoefficients)
        double p = sigmoid(ab[0] * rawScore + ab[1])
        Math.max(EPSILON, Math.min(1.0 - EPSILON, p))
    }

    private static double sigmoid(double x) {
        1.0 / (1.0 + Math.exp(-x))
    }

    /**
     * Default Platt coefficients per known predictor. These are provisional
     * defaults chosen so that the highest plausible raw score for each tool
     * maps to ~0.9 and the lowest plausible raw score maps to ~0.2. Phase 7.5
     * replaces these with learned values.
     */
    private void loadDefaults() {
        // Similarity tools: raw = pident / 100 ∈ [0.3, 1.0] (min_id filters).
        // At raw=0.3 → ~0.2; at raw=0.9 → ~0.90.
        register('diamond',  6.0d, -3.0d)
        register('mmseqs2',  6.0d, -3.0d)

        // Domain/HMM tools: raw = bitscore / 100 is already roughly calibrated.
        register('hmmer',    4.0d, -1.5d)
        register('pfam',     4.0d, -1.5d)
        register('interproscan', 5.0d, -2.0d)

        // Structure similarity: FoldSeek. Raw uses TM-score-like quantity,
        // slightly more conservative than sequence.
        register('foldseek', 5.0d, -2.0d)

        // Orthology: eggNOG-mapper emits rawScore = min(1, -log10(evalue)/50)
        // clamped to >=0.4. Gentle curve — the raw score is already tuned.
        register('eggnog-mapper', 3.5d, -1.3d)
        register('eggnog',        3.5d, -1.3d)

        // Localization: binary-ish with moderate confidence.
        register('signalp',     3.0d, -1.0d)
        register('deeptmhmm',   3.0d, -1.0d)

        // Operon / context: weak signal by itself, conservative calibration.
        register('operon',      3.0d, -1.2d)

        // Pathway / metabolic: gapseq output.
        register('gapseq',      4.0d, -1.8d)

        // Domain-specific tools output high-confidence gene hits.
        register('amrfinder',   8.0d, -4.0d)
        register('dbcan',       5.0d, -2.0d)
        register('antismash',   4.0d, -1.5d)
    }
}
