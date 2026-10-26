package gspa.integration.learned

import gspa.model.AnnotationType

/**
 * Feature vector of one candidate: protein x function x role. The role is part of the
 * identity, so the same protein/function with two roles yields two candidates.
 */
class CandidateFeatures {

    final String proteinId
    final AnnotationType functionType
    final String functionId
    final CandidateRole role
    final FeatureContract contract
    final double[] values
    /** Human-readable notes on what was counted / ignored (audit trail, not a model input). */
    final List<String> provenance

    CandidateFeatures(String proteinId, AnnotationType functionType, String functionId, CandidateRole role,
                      FeatureContract contract, double[] values, List<String> provenance = []) {
        if (values.length != contract.size()) {
            throw new IllegalArgumentException("feature vector has ${values.length} values, contract ${contract.id} has ${contract.size()}")
        }
        this.proteinId = proteinId
        this.functionType = functionType
        this.functionId = functionId
        this.role = role
        this.contract = contract
        this.values = values
        this.provenance = provenance
    }

    /** Legacy 3-part key, role-agnostic. */
    String functionKey() { "${proteinId}|${functionType}|${functionId}".toString() }

    /** Candidate identity including the role. */
    String candidateKey() { "${functionKey()}|${role.name()}".toString() }

    double get(String featureName) {
        int i = contract.indexOf(featureName)
        if (i < 0) throw new IllegalArgumentException("unknown feature ${featureName}")
        values[i]
    }

    Map<String, Double> asMap() {
        Map<String, Double> m = new LinkedHashMap<>()
        contract.names().eachWithIndex { String n, int i -> m[n] = values[i] }
        m
    }

    @Override
    String toString() { "CandidateFeatures(${candidateKey()})" }
}
