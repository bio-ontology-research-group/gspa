package gspa.integration.prior

import gspa.integration.IntegrationState
import gspa.integration.MetabolicGap
import gspa.integration.Prior
import gspa.model.AnnotationType
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Boosts candidate claims whose function matches a gapseq-identified
 * metabolic gap.
 *
 * <p>Plan §A.4 of Phase 7. For each gap {@code (pathway P, reaction R)}
 * with EC number {@code E} (and mapped GO term {@code F}), any candidate
 * claim whose functionId equals {@code E} or {@code F} earns
 * {@code +alphaGap}. Gapseq's own gap-fill guesses (marked
 * {@code gapseqGuessed = true}) get a reduced {@code +alphaGap * 0.7}.</p>
 *
 * <p>The metabolic gap list is a pre-computed input on
 * {@link IntegrationState#metabolicGaps}; the caller (pipeline) populates
 * it after running {@code GapseqPredictor}.</p>
 */
class GapFillingPrior implements Prior {

    private static final Logger log = LoggerFactory.getLogger(GapFillingPrior)

    /** Per-match boost magnitude (log-odds). Learned by the benchmark. */
    double alphaGap = 1.5d

    /** Discount for gapseq's own gap-fill guesses. */
    double gapseqGuessedFactor = 0.7d

    /** Cached per iteration: function_id → strongest boost magnitude. */
    private Map<String, Double> boostByFunctionId = new LinkedHashMap<>()

    @Override
    String name() { 'gap_filling' }

    @Override
    Set<String> inputs() { ['metabolicGaps'] as Set }

    @Override
    void beginIteration(IntegrationState state) {
        boostByFunctionId = new LinkedHashMap<>()
        if (state.metabolicGaps == null || state.metabolicGaps.isEmpty()) return

        for (MetabolicGap gap : state.metabolicGaps) {
            double mag = gap.gapseqGuessed ? (alphaGap * gapseqGuessedFactor) : alphaGap
            if (gap.ecNumber) updateMax(gap.ecNumber, mag)
            if (gap.goTerm) updateMax(gap.goTerm, mag)
        }
        log.debug("GapFillingPrior: ${boostByFunctionId.size()} functions associated with metabolic gaps")
    }

    @Override
    double logOddsBoost(String proteinId, String functionKey, IntegrationState state) {
        if (boostByFunctionId.isEmpty()) return 0.0d
        String[] parts = IntegrationState.splitFunctionKey(functionKey)
        if (parts == null) return 0.0d
        String type = parts[1]
        String funcId = parts[2]
        if (type != AnnotationType.EC.name() && type != AnnotationType.GO.name()) return 0.0d
        return boostByFunctionId.getOrDefault(funcId, 0.0d)
    }

    private void updateMax(String key, double value) {
        boostByFunctionId.merge(key, value, { a, b -> Math.max(a as double, b as double) })
    }
}
