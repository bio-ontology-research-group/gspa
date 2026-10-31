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

    private static final Set<String> STATUS_NAMES = ['SKIPPED', 'UNAVAILABLE', 'FAILED', 'PARTIAL', 'SUCCESS_ZERO', 'SUCCESS'] as Set

    /**
     * Predictor completion outcomes by predictor name (PredictorOutcome names:
     * SKIPPED, UNAVAILABLE, FAILED, PARTIAL, SUCCESS_ZERO, SUCCESS). Null or empty = not supplied.
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

    /**
     * Cross-genome recurrence. The integer counts are the legacy, unweighted summary (all zero =
     * not supplied). {@link #weighted} carries the phylogeny-weighted, background-matched
     * features of a joint cohort analysis (see {@code gspa.integration.contextbridge}); a null
     * field there means "unavailable", never zero.
     */
    static class Cohort {
        gspa.integration.cohort.ContextMode contextMode = gspa.integration.cohort.ContextMode.UNIVERSAL
        /** Only within-genome reaction and intrinsic evidence; peer statistics must be absent. */
        boolean localOnly = false
        int nGenomes
        int nWithOrthogroup
        int nWithFunctionSeq
        int nWithContext
        int nIncomplete
        /** Weighted cohort features; null when no joint cohort analysis supplied them. */
        Weighted weighted = null
        /** Reaction the weighted features were computed for (an EC may map to several reactions). */
        String reactionId = null
        /** Other reactions the same EC call maps to: reported, never summed into the evidence. */
        List<ReactionAlternative> reactionAlternatives = []
    }

    /** One further reaction hypothesis of the same candidate; descriptive, not an additional evidence term. */
    static class ReactionAlternative {
        String reactionId
        Double supportFraction
        Double effectiveN
    }

    /**
     * Typed weighted cohort features, one nullable field per contract column family. Fractions are
     * in [0,1]; counts are non-negative. Values copy the cohort engine's output without rounding.
     */
    static class Weighted {
        // peer set and uncertainty
        Double effectiveN                  // summed weight of informative independent units
        boolean insufficientPeers          // peers exist, none informative
        boolean adjustedAvailable          // clade/cluster-adjusted outputs exist (metadata known)
        Double phylogenyKnownFraction
        Double nPeerGenomes
        Double coverage                    // share of peer clades that were informative
        Double supportFraction             // Jeffreys-smoothed, weighted
        Double supportFractionLo90         // null unless adjustedAvailable
        Double supportFractionHi90
        Double supportFractionRaw          // unsmoothed weighted fraction
        Double supportUnknownAsNegative
        Double supportUnknownAsSupport
        // observation states of peers (weighted fractions)
        Double fracSupported
        Double fracNegative
        Double fracUnknownContext
        Double fracTrueAbsent
        Double fracUnobservedAbsent
        Double fracPathwayWithoutFamily
        // matched same-genome background
        Double bgRate
        Double bgWindows
        Double enrichmentLog
        Double supportMinusBg
        // complementary partial contexts across lineages
        Double complementarity
        Double complementaryUnannotated    // unsupported label only; never evidence of loss
        Double jointCoverage
        // within-genome permutation null
        Double permP
        Double permExcess                  // observed minus null mean
        // paralogy / copy number
        Double multicopyFraction
        // organism-specific reaction plausibility of peers
        Double peerOrganismSetKnownFraction
        Double peerReactionInOrganismFraction
        // the candidate's own window (cohort engine's local features)
        boolean layoutAvailable
        Double localAnchorWeight
        Double localAnchorCount
        Double localAnchorCurated
        Double localAnchorUpstream
        Double localAnchorDownstream
        Double localContextPositive
        Double localWindowTruncated
        Double localWindowUnknown
        Double localEnrichmentLog
        // how the candidate arose and how ambiguous it is
        boolean ecDirect
        Double ecReactionCount
        Double nReactionAlternatives
        Double nFamilyAlternatives
        boolean familyKnown
        Double membershipWeight
        // raw protein-intrinsic role evidence (missing is not "non-catalytic")
        boolean roleEvidenceMissing = true
        Map<CandidateRole, Double> roleEvidence = [:]
    }

    static class FbaConditional {
        int nConditions
        /** Fraction of the tested conditions in which the candidate's reaction carried required flux. */
        double requiredFraction
        /** Fraction in which the reaction could carry flux without being required (optional flux); null = not reported. */
        Double optionalFraction
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
        cohort?.each { String k, Cohort c -> if (c.weighted != null) validateWeighted(k, c.weighted) }
        cohort?.each { String k, Cohort c ->
            if (c.contextMode == null) throw new IllegalArgumentException("cohort[${k}] contextMode is required")
            if (!c.localOnly) return
            if (c.weighted == null || [c.nGenomes, c.nWithOrthogroup, c.nWithFunctionSeq, c.nWithContext, c.nIncomplete].any { it != 0 } ||
                c.weighted.insufficientPeers || c.weighted.adjustedAvailable ||
                CohortFeatureColumns.peerProperties().any { c.weighted.getProperty(it) != null } ||
                c.reactionAlternatives.any { it.supportFraction != null || it.effectiveN != null }) {
                throw new IllegalArgumentException("cohort[${k}] localOnly requires weighted local evidence and absent peer statistics")
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
            if (f.optionalFraction != null && (!Double.isFinite(f.optionalFraction) || f.optionalFraction < 0d ||
                f.optionalFraction > 1d || f.requiredFraction + f.optionalFraction > 1d + 1e-9d)) {
                throw new IllegalArgumentException("fba[${k}].optionalFraction is outside [0,1] or exceeds 1 - requiredFraction")
            }
        }
        roles?.each { String k, List<CandidateRole> rs ->
            if (rs != null && rs.any { it == null }) throw new IllegalArgumentException("roles[${k}] contains null")
        }
        if (assemblyCompleteness != null && (!Double.isFinite(assemblyCompleteness) || assemblyCompleteness < 0d || assemblyCompleteness > 100d)) {
            throw new IllegalArgumentException("assemblyCompleteness must be a percentage in [0,100], got ${assemblyCompleteness}")
        }
    }

    private static final Set<String> WEIGHTED_FRACTIONS = ['phylogenyKnownFraction', 'coverage', 'supportFraction',
        'supportFractionLo90', 'supportFractionHi90', 'supportFractionRaw', 'supportUnknownAsNegative',
        'supportUnknownAsSupport', 'fracSupported', 'fracNegative', 'fracUnknownContext', 'fracTrueAbsent',
        'fracUnobservedAbsent', 'fracPathwayWithoutFamily', 'bgRate', 'jointCoverage', 'complementaryUnannotated',
        'permP', 'multicopyFraction', 'peerOrganismSetKnownFraction', 'peerReactionInOrganismFraction',
        'localContextPositive', 'localWindowTruncated', 'localWindowUnknown', 'membershipWeight'] as Set
    private static final Set<String> WEIGHTED_COUNTS = ['effectiveN', 'nPeerGenomes', 'bgWindows', 'localAnchorWeight',
        'localAnchorCount', 'localAnchorCurated', 'localAnchorUpstream', 'localAnchorDownstream', 'ecReactionCount',
        'nReactionAlternatives', 'nFamilyAlternatives'] as Set
    private static final Set<String> WEIGHTED_SIGNED = ['enrichmentLog', 'supportMinusBg', 'complementarity', 'permExcess'] as Set

    private static void validateWeighted(String key, Weighted w) {
        WEIGHTED_FRACTIONS.each { String n ->
            Double v = w.getProperty(n) as Double
            if (v != null && (!Double.isFinite(v) || v < 0d || v > 1d)) {
                throw new IllegalArgumentException("cohort[${key}].weighted.${n} must be finite and in [0,1], got ${v}")
            }
        }
        WEIGHTED_COUNTS.each { String n ->
            Double v = w.getProperty(n) as Double
            if (v != null && (!Double.isFinite(v) || v < 0d)) {
                throw new IllegalArgumentException("cohort[${key}].weighted.${n} must be finite and >= 0, got ${v}")
            }
        }
        WEIGHTED_SIGNED.each { String n ->
            Double v = w.getProperty(n) as Double
            if (v != null && !Double.isFinite(v)) throw new IllegalArgumentException("cohort[${key}].weighted.${n} is not finite")
        }
        if (!w.adjustedAvailable && (w.supportFractionLo90 != null || w.supportFractionHi90 != null || w.permP != null)) {
            throw new IllegalArgumentException("cohort[${key}].weighted carries adjusted intervals/permutation p without adjustedAvailable")
        }
        w.roleEvidence?.each { CandidateRole r, Double v ->
            if (v == null || !Double.isFinite(v) || v < 0d || v > 1d) {
                throw new IllegalArgumentException("cohort[${key}].weighted.roleEvidence[${r}] must be finite and in [0,1]")
            }
        }
        if (w.roleEvidenceMissing && w.roleEvidence) {
            throw new IllegalArgumentException("cohort[${key}].weighted reports role evidence as missing but lists values")
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
