package gspa.integration.cohort

import groovy.transform.CompileStatic

/**
 * Resolves a protein's role distribution from protein-intrinsic evidence
 * only. Context (operon, neighbourhood, cohort recurrence) is intentionally
 * not an argument: it cannot make a protein catalytic.
 *
 * <p>p_r = s_r / (sum_s + u), p_unknown = u / (sum_s + u) with pseudo-mass u,
 * so a protein with no evidence is UNKNOWN with probability 1.</p>
 */
@CompileStatic
class RoleResolver {

    static class Resolution {
        Map<Role, Double> probabilities = new LinkedHashMap<>()
        /** Argmax role if its probability is at least 0.5, otherwise UNKNOWN. */
        Role resolved = Role.UNKNOWN
        /** True when catalytic probability rests on sequence/EC evidence (always true if > 0). */
        boolean catalyticEvidenceIndependentOfContext
    }

    static Resolution resolve(Map<Role, Double> scores, double ecDirectConfidence, double unknownPseudo) {
        Map<Role, Double> s = new EnumMap<Role, Double>(Role)
        for (Role r : Role.values()) {
            if (r != Role.UNKNOWN) s.put(r, 0d)
        }
        if (scores != null) {
            scores.each { Role r, Double v ->
                if (r != Role.UNKNOWN && v != null && v > 0d) s.put(r, Math.min(1d, v))
            }
        }
        if (ecDirectConfidence > s.get(Role.CATALYTIC).doubleValue()) s.put(Role.CATALYTIC, Math.min(1d, ecDirectConfidence))
        double total = 0d
        s.values().each { total += it }
        double denom = total + unknownPseudo
        Resolution out = new Resolution()
        Role best = Role.UNKNOWN
        double bestP = unknownPseudo / denom
        for (Role r : Role.values()) {
            double p = r == Role.UNKNOWN ? unknownPseudo / denom : s.get(r).doubleValue() / denom
            out.probabilities.put(r, p)
            if (r != Role.UNKNOWN && p > bestP) { best = r; bestP = p }
        }
        out.resolved = bestP >= 0.5d ? best : Role.UNKNOWN
        out.catalyticEvidenceIndependentOfContext = s.get(Role.CATALYTIC).doubleValue() > 0d
        out
    }
}
