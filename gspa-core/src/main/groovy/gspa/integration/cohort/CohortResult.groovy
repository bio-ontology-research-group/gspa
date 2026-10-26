package gspa.integration.cohort

import groovy.transform.CompileStatic

/**
 * Output of cohort inference: per-(gene, reaction) feature records and
 * per-(family, reaction) aggregates, with run provenance.
 *
 * <p>Feature maps omit unavailable values (missing, not zero); explicit
 * {@code *_available} / {@code *_unknown} flags are present as 0/1.</p>
 */
@CompileStatic
class CohortResult {

    /** Another family or reaction hypothesis for the same gene; never dropped. */
    static class Alternative {
        /** "family" or "reaction". */
        String kind
        String familyKey
        String reactionId
        double weight
        Map<String, Double> features = new LinkedHashMap<>()
    }

    static class FeatureRecord {
        String genomeId
        String proteinId
        String reactionId
        /** Primary (highest-weight) family; null when the protein has no family assignment. */
        String familyKey
        double membershipWeight
        double baselineConfidence
        Map<String, Double> features = new LinkedHashMap<>()
        Map<String, String> states = new LinkedHashMap<>()
        Map<Role, Double> roleProbabilities = new LinkedHashMap<>()
        Role resolvedRole = Role.UNKNOWN
        List<Alternative> alternatives = new ArrayList<>()
        Map<String, Object> provenance = new LinkedHashMap<>()
    }

    static class Aggregate {
        String familyKey
        String reactionId
        /** Always "none": aggregates are descriptive and must not be used to score a member genome. */
        String holdout = 'none'
        int nCandidateGenes
        Map<String, Double> features = new LinkedHashMap<>()
        Map<String, Object> provenance = new LinkedHashMap<>()
    }

    final List<FeatureRecord> records = new ArrayList<>()
    final List<Aggregate> aggregates = new ArrayList<>()
    final Map<String, Object> provenance = new LinkedHashMap<>()

    FeatureRecord record(String genomeId, String proteinId, String reactionId) {
        records.find { it.genomeId == genomeId && it.proteinId == proteinId && it.reactionId == reactionId }
    }

    Aggregate aggregate(String familyKey, String reactionId) {
        aggregates.find { it.familyKey == familyKey && it.reactionId == reactionId }
    }
}
