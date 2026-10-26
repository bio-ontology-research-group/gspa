package gspa.integration.cohort

import groovy.transform.CompileStatic

/**
 * Raw role evidence for one protein, passed through for the learned stage.
 *
 * <p>This class deliberately computes no role probability, no null/unknown
 * distribution and no argmax or confidence threshold: that is the learned
 * stage's job. It lists each piece of protein-intrinsic evidence with its
 * source, states whether evidence is missing, and enumerates the role
 * alternatives that remain possible. A direct EC call is reported as raw
 * support that the catalytic role is possible; its confidence is the
 * baseline-snapshot confidence unchanged and is not a probability of role.
 * Context (operon, neighbourhood, cohort recurrence) is never an input.</p>
 */
@CompileStatic
class RoleEvidenceSummary {

    static class Item {
        Role role
        String source
        /** Raw score as supplied (0..1), unchanged. */
        double value
    }

    final List<Item> items = new ArrayList<>()
    /** True when the protein has no protein-intrinsic role evidence at all (missing, not "non-catalytic"). */
    boolean missing = true
    /** Roles with any evidence, in enum order, plus UNKNOWN which is always an open possibility. */
    final List<Role> alternatives = new ArrayList<>()

    /** Largest raw value per role (absent when the role has no evidence). */
    Map<Role, Double> maxByRole() {
        Map<Role, Double> out = new LinkedHashMap<>()
        for (Item i : items) {
            Double prev = out.get(i.role)
            if (prev == null || i.value > prev.doubleValue()) out.put(i.role, i.value)
        }
        out
    }

    static RoleEvidenceSummary summarize(CohortInput.RoleEvidence evidence, Collection<CohortInput.CandidateClaim> claims) {
        RoleEvidenceSummary s = new RoleEvidenceSummary()
        if (evidence != null) {
            for (Role r : Role.values()) {
                Double v = evidence.scores.get(r)
                if (v != null) s.items.add(new Item(role: r, source: evidence.source ?: 'unspecified', value: v.doubleValue()))
            }
        }
        for (CohortInput.CandidateClaim c : claims) {
            if (c.ecDirect) {
                s.items.add(new Item(role: Role.CATALYTIC, source: 'ec-direct:' + c.ecNumber + '@' + (c.provenance ?: 'unspecified'),
                        value: c.baselineConfidence))
            }
        }
        s.missing = s.items.isEmpty()
        Set<Role> seen = new LinkedHashSet<>()
        for (Item i : s.items) seen.add(i.role)
        for (Role r : Role.values()) if (seen.contains(r) && r != Role.UNKNOWN) s.alternatives.add(r)
        s.alternatives.add(Role.UNKNOWN)
        s
    }
}
