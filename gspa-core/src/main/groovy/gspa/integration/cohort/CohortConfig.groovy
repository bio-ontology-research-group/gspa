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
    /** Pseudo-mass of the UNKNOWN role when resolving roles from evidence. */
    double unknownRolePseudo = 0.5d
    /** Permutations for the matched-window null (0 disables). */
    int permutations = 200
    long seed = 20261026L
    HoldoutLevel holdout = HoldoutLevel.STRAIN_CLUSTER
    /** Cut height for clade extraction from a supplied tree (max leaf-to-leaf distance within a clade). */
    double treeCutDistance = 0.05d

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
        m.put('unknownRolePseudo', unknownRolePseudo)
        m.put('permutations', permutations)
        m.put('seed', seed)
        m.put('holdout', holdout.name())
        m.put('treeCutDistance', treeCutDistance)
        m
    }
}
