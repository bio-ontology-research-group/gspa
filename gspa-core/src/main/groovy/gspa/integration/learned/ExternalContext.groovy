package gspa.integration.learned

import gspa.integration.IntegrationState
import gspa.model.AnnotationType

/**
 * Optional, structured context the integration state itself does not carry. Everything
 * here is supplied by the caller from a real analysis (reaction-graph neighbourhood,
 * cohort recurrence, FBA, predictor run manifest); the extractor never fabricates it.
 * Anything absent is reported to the model as unavailable (indicator 0), not as a
 * measured zero.
 *
 * <p>All maps are keyed by the legacy function key {@code proteinId|TYPE|functionId}
 * (see {@link gspa.integration.IntegrationState#functionKey}); they apply to every role of
 * that function.
 *
 * <p>A key that appears in any of the maps (or in {@link #candidates}) names a candidate even
 * if no raw predictor claim exists for it: a weak, context-only candidate is scored with no
 * sequence evidence rather than ignored (see {@link #candidateKeys}).
 */
class ExternalContext {

    private static final Set<String> STATUS_NAMES = ['UNAVAILABLE', 'FAILED', 'PARTIAL', 'SUCCESS_ZERO', 'SUCCESS'] as Set

    /**
     * Predictor completion outcomes by predictor name (PredictorOutcome names:
     * UNAVAILABLE, FAILED, PARTIAL, SUCCESS_ZERO, SUCCESS). Null or empty = not supplied.
     */
    Map<String, String> predictorStatus = null

    /** Neighbourhood features computed from genome context + reaction graph. */
    Map<String, Neighborhood> neighborhood = [:]

    /** Cross-genome recurrence from a cohort analysis. */
    Map<String, Cohort> cohort = [:]

    /** Conditional FBA results. */
    Map<String, FbaConditional> fba = [:]

    /**
     * Explicit roles from structured evidence (e.g. a transporter annotation). When
     * absent, roles are derived conservatively by {@link ContextFeatureExtractor}.
     */
    Map<String, List<CandidateRole>> roles = [:]

    /** Additional candidate function keys proposed by structured context (no claim required). */
    Set<String> candidates = [] as Set

    /** External assembly completeness estimate in percent (e.g. CheckM2), or null. */
    Double assemblyCompleteness = null

    static class Neighborhood {
        double supportMax
        int nSupporting
        double conflictMax
        /** Reaction-graph distance to the nearest supported neighbour reaction; null = unavailable. */
        Double reactionDistance
    }

    static class Cohort {
        int nGenomes
        int nWithOrthogroup
        int nWithFunctionSeq
        int nWithContext
        int nIncomplete
    }

    static class FbaConditional {
        int nConditions
        double requiredFraction
        /** Relative growth change on knockout; null = not computed. */
        Double knockoutGrowthEffect
        /** True when the tested objective itself required the candidate (not independent evidence). */
        boolean objectiveRequiresCandidate
    }

    /** Every function key named by the structured context, in a deterministic order. */
    Set<String> candidateKeys() {
        Set<String> keys = new TreeSet<>()
        [neighborhood?.keySet(), cohort?.keySet(), fba?.keySet(), roles?.keySet(), candidates].each { Collection<String> c ->
            if (c != null) keys.addAll(c)
        }
        keys
    }

    /**
     * Reject values that would produce a non-finite or impossible feature. Called before any
     * extraction so that a bad value fails the whole run instead of one candidate half-way.
     */
    void validate() {
        candidateKeys().each { String k -> requireFunctionKey(k) }
        if (predictorStatus != null) {
            predictorStatus.each { String name, String st ->
                if (!STATUS_NAMES.contains(st)) throw new IllegalArgumentException("predictorStatus[${name}] = ${st} is not one of ${STATUS_NAMES}")
            }
        }
        neighborhood?.each { String k, Neighborhood n ->
            if (n == null) throw new IllegalArgumentException("neighborhood[${k}] is null")
            if (!Double.isFinite(n.supportMax) || !Double.isFinite(n.conflictMax) || n.supportMax < 0d || n.conflictMax < 0d || n.nSupporting < 0) {
                throw new IllegalArgumentException("neighborhood[${k}] has a non-finite or negative value")
            }
            if (n.reactionDistance != null && (!Double.isFinite(n.reactionDistance) || n.reactionDistance < 0d)) {
                throw new IllegalArgumentException("neighborhood[${k}].reactionDistance must be finite and >= 0 (or null)")
            }
        }
        cohort?.each { String k, Cohort c ->
            if (c == null) throw new IllegalArgumentException("cohort[${k}] is null")
            if (c.nGenomes < 0 || c.nWithOrthogroup < 0 || c.nWithFunctionSeq < 0 || c.nWithContext < 0 || c.nIncomplete < 0 ||
                c.nWithOrthogroup > c.nGenomes || c.nWithFunctionSeq > c.nWithOrthogroup || c.nWithContext > c.nWithOrthogroup ||
                c.nIncomplete > c.nGenomes) {
                throw new IllegalArgumentException("cohort[${k}] counts are negative or inconsistent (a subset larger than its superset)")
            }
        }
        fba?.each { String k, FbaConditional f ->
            if (f == null) throw new IllegalArgumentException("fba[${k}] is null")
            if (f.nConditions < 0 || !Double.isFinite(f.requiredFraction) || f.requiredFraction < 0d || f.requiredFraction > 1d) {
                throw new IllegalArgumentException("fba[${k}] has an invalid condition count or requiredFraction outside [0,1]")
            }
            if (f.knockoutGrowthEffect != null && !Double.isFinite(f.knockoutGrowthEffect)) {
                throw new IllegalArgumentException("fba[${k}].knockoutGrowthEffect is not finite")
            }
        }
        roles?.each { String k, List<CandidateRole> rs ->
            if (rs != null && rs.any { it == null }) throw new IllegalArgumentException("roles[${k}] contains null")
        }
        if (assemblyCompleteness != null && (!Double.isFinite(assemblyCompleteness) || assemblyCompleteness < 0d || assemblyCompleteness > 100d)) {
            throw new IllegalArgumentException("assemblyCompleteness must be a percentage in [0,100], got ${assemblyCompleteness}")
        }
    }

    /** Parse and validate a function key; returns [proteinId, type, functionId]. */
    static String[] requireFunctionKey(String key) {
        String[] parts = IntegrationState.splitFunctionKey(key)
        if (parts == null) throw new IllegalArgumentException("invalid function key '${key}' (expected proteinId|TYPE|functionId)")
        try { AnnotationType.valueOf(parts[1]) } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException("function key '${key}' has unknown type '${parts[1]}'")
        }
        if (!parts[2]) throw new IllegalArgumentException("function key '${key}' has an empty function id")
        parts
    }
}
