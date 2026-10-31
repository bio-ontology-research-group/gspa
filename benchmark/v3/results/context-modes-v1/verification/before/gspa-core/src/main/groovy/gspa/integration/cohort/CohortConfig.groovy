package gspa.integration.cohort

import groovy.transform.CompileStatic

/**
 * Parameters of cohort inference. None of them is a decision threshold on
 * "how often" a pattern must recur; that is left to calibration downstream.
 */
@CompileStatic
class CohortConfig {
    /** Half-width (bp) of the neighbourhood window around a gene midpoint. */
    int halfWidthBp = 5000
    /** Minimum baseline confidence for a neighbouring gene's claim to act as an anchor. */
    double anchorMinConfidence = 0.5d
    /** Distinct anchor reactions required for a window to count as context-positive. */
    int minAnchors = 1
    /** Network hops used to define the reaction signature (curated edges are always hop 1). */
    int maxHops = 1
    /** Distance decay alpha^d for anchor weights. */
    double hopAlpha = 0.5d
    /** Matched background windows may differ from the focal window by this many genes. */
    int densityTolerance = 1
    /** Membership weights below this are ignored. */
    double minMembershipWeight = 0.05d
    /** Total weight of one independent clade, split over its strain clusters. */
    double cladeCap = 1.0d
    /** Effective weight needed for a context reaction to count as jointly covered. */
    double jointMinEffective = 0.5d
    /** Permutations for the matched-window null (0 disables). */
    int permutations = 200
    long seed = 20261026L
    HoldoutLevel holdout = HoldoutLevel.STRAIN_CLUSTER
    /** Cut height for clade extraction from a supplied tree (max leaf-to-leaf distance within a clade). */
    double treeCutDistance = 0.05d
    /** Whether to compute descriptive all-genome aggregates in CohortResult.aggregates. */
    boolean includeAggregates = true
    /** Maximum per-record reaction alternatives materialized on FeatureRecord.alternatives (n_reaction_alternatives feature stays exact). */
    int maxReactionAlternatives = 32
    /** Maximum EC calls per protein expanded into reaction candidates for cohort recurrence (keeps all ties at cutoff). */
    int maxEcCallsPerProtein = 128

    void validate() {
        if (halfWidthBp <= 0) throw new IllegalArgumentException('halfWidthBp must be positive')
        if (!(anchorMinConfidence >= 0d && anchorMinConfidence <= 1d)) throw new IllegalArgumentException('anchorMinConfidence must be in [0,1]')
        if (minAnchors < 1) throw new IllegalArgumentException('minAnchors must be >= 1')
        if (maxHops < 1) throw new IllegalArgumentException('maxHops must be >= 1')
        if (!(hopAlpha > 0d && hopAlpha <= 1d)) throw new IllegalArgumentException('hopAlpha must be in (0,1]')
        if (densityTolerance < 0) throw new IllegalArgumentException('densityTolerance must be >= 0')
        if (!(minMembershipWeight >= 0d && minMembershipWeight <= 1d)) throw new IllegalArgumentException('minMembershipWeight must be in [0,1]')
        if (!(cladeCap > 0d) || Double.isInfinite(cladeCap)) throw new IllegalArgumentException('cladeCap must be positive and finite')
        if (!(jointMinEffective >= 0d)) throw new IllegalArgumentException('jointMinEffective must be >= 0')
        if (permutations < 0) throw new IllegalArgumentException('permutations must be >= 0')
        if (!(treeCutDistance >= 0d)) throw new IllegalArgumentException('treeCutDistance must be >= 0')
    }

    Map<String, Object> toMap() {
        Map<String, Object> m = new LinkedHashMap<>()
        m.put('halfWidthBp', halfWidthBp)
        m.put('anchorMinConfidence', anchorMinConfidence)
        m.put('minAnchors', minAnchors)
        m.put('maxHops', maxHops)
        m.put('hopAlpha', hopAlpha)
        m.put('densityTolerance', densityTolerance)
        m.put('minMembershipWeight', minMembershipWeight)
        m.put('cladeCap', cladeCap)
        m.put('jointMinEffective', jointMinEffective)
        m.put('permutations', permutations)
        m.put('seed', seed)
        m.put('holdout', holdout.name())
        m.put('treeCutDistance', treeCutDistance)
        m
    }
}
