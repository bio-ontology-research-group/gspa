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
    /** Exact normalized names of usable baseline sources, before pooling into feature columns. */
    final List<String> sequenceSources
    final List<Map> sequenceContracts
    /** Preserved observations; outside the numeric contract and unsupported by current models. */
    final List<Map> qualityEvidence
    final gspa.integration.cohort.ContextMode cohortContextMode

    CandidateFeatures(String proteinId, AnnotationType functionType, String functionId, CandidateRole role,
                      FeatureContract contract, double[] values, List<String> provenance = [], List<String> sequenceSources = null,
                      List<Map> qualityEvidence = [],
                      gspa.integration.cohort.ContextMode cohortContextMode = gspa.integration.cohort.ContextMode.UNIVERSAL,
                      List<Map> sequenceContracts = null) {
        if (cohortContextMode == null) throw new IllegalArgumentException('cohort context mode is required')
        this.cohortContextMode = cohortContextMode
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
        this.sequenceSources = sequenceSources == null ? null : Collections.unmodifiableList(new ArrayList<>(sequenceSources))
        this.sequenceContracts = sequenceContracts == null ? null : PredictorContract.identities(sequenceContracts)
        this.qualityEvidence = Collections.unmodifiableList(new ArrayList<>(qualityEvidence))
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

    /** SHA-256 over the exact IEEE-754 bits of the vector, in contract order (full-precision audit identity). */
    String digest() {
        def md = java.security.MessageDigest.getInstance('SHA-256')
        def buf = java.nio.ByteBuffer.allocate(8 * values.length)
        for (double v : values) buf.putLong(Double.doubleToLongBits(v))
        md.update(buf.array())
        md.digest().collect { String.format('%02x', it & 0xff) }.join()
    }

    /** Fail unless every value is finite; a non-finite feature can never be scored or exported. */
    void requireFinite() {
        for (int i = 0; i < values.length; i++) {
            if (!Double.isFinite(values[i])) throw new IllegalArgumentException("non-finite value for feature ${contract.names()[i]} of ${candidateKey()}")
        }
    }

    Map<String, Double> asMap() {
        Map<String, Double> m = new LinkedHashMap<>()
        contract.names().eachWithIndex { String n, int i -> m[n] = values[i] }
        m
    }

    @Override
    String toString() { "CandidateFeatures(${candidateKey()})" }
}
