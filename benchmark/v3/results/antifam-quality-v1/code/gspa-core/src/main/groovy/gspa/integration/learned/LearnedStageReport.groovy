package gspa.integration.learned

import gspa.integration.EvidenceClaim
import gspa.model.AnnotationType

/**
 * Audit trail of one learned-stage run: the model identity, one record per scored candidate
 * (role-specific, never merged) and one outcome per function key (what happened to the legacy,
 * role-free annotation).
 */
class LearnedStageReport {

    /** One scored (protein, function, role) candidate. */
    static class CandidateRecord {
        String candidateKey
        String functionKey
        String proteinId
        AnnotationType functionType
        String functionId
        CandidateRole role
        /** Pre-calibration logit of the frozen model. */
        Double margin
        /** Calibrated marginal probability that this candidate is correct (independent of all others). */
        Double probability
        String scoreStatus = 'scored'
        String unsupportedReason
        double displayThreshold
        Boolean aboveDisplayThreshold
        /** Exact-bits SHA-256 of the feature vector that was scored. */
        String featureDigest
        /** The full vector in contract order; null unless {@code retainFeatureVectors}. */
        double[] features
        /** Extractor notes (what was counted / ignored). */
        List<String> featureNotes
        List<String> sequenceSources
        List<Map> qualityEvidence = []
        /** The original raw claims behind the sequence features (empty for context-only candidates). */
        List<EvidenceClaim> supportingClaims
        /** True if this role's probability is the one written to the legacy role-free annotation. */
        boolean projected
        /** Conditional probability of being THE locus within its function/role group; null if no locus layer. */
        Double locusProbability
        /** Conditional probability that none of the group's supplied candidates is correct; null if no locus layer. */
        Double locusNullProbability
    }

    /** What happened to one function key in the legacy (role-free) output. */
    static class KeyOutcome {
        String functionKey
        /** rescored | added | retracted | removed_no_projection | scored_only */
        String action
        String reason
        /** Role whose probability became the legacy score; null if none was projectable. */
        CandidateRole projectedRole
        /** Probability before the stage (the replaced posterior); null if the key was not retained before. */
        Double priorProbability
        Double priorLogOdds
        Map<String, Double> priorContributions
        /** Learned probability of the projected role; null if no role was projectable. */
        Double learnedProbability
        /** Learned probability of every role of the key, retained even though legacy output holds one. */
        Map<String, Double> roleProbabilities
    }

    String modelId
    String modelSha256
    String modelKind
    String modelIntendedUse
    String contractId
    String contractFingerprint
    Map modelTraining
    LearnedStageSettings settings
    String locusUnavailableReason

    List<CandidateRecord> candidates = []
    List<KeyOutcome> outcomes = []
    /** Pinned posterior floors dropped for rescored keys so that no later stage can floor a learned score. */
    List<String> pinnedFloorsCleared = []

    int count(String action) { outcomes.count { it.action == action } as int }
    int unsupportedCount() { candidates.count { it.scoreStatus == 'unsupported' } as int }

    @Override
    String toString() {
        "LearnedStageReport(model=${modelId}, ${candidates.size()} candidates, rescored=${count('rescored')}, added=${count('added')}, " +
            "retracted=${count('retracted')}, removed_no_projection=${count('removed_no_projection')}, scored_only=${count('scored_only')})"
    }
}
