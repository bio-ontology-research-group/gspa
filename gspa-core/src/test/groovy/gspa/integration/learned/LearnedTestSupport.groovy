package gspa.integration.learned

import groovy.json.JsonOutput
import gspa.integration.EvidenceClaim
import gspa.integration.EvidenceType
import gspa.model.AnnotationType

/**
 * Test-only helpers. The models built here have HAND-SET weights and are marked
 * {@code intended_use = mechanism_test}: they exist to exercise the mechanisms (scores can fall,
 * null outcome is non-zero, ...) and are never evidence about prediction quality. The production
 * loader refuses such artifacts unless a caller explicitly opts in, as these tests do.
 */
class LearnedTestSupport {

    static EvidenceClaim claim(String protein, String function, double raw, String source = 'diamond',
                               EvidenceType type = EvidenceType.SEQUENCE_SIMILARITY,
                               AnnotationType ft = AnnotationType.EC, Map metadata = [:], String aspect = null) {
        new EvidenceClaim(proteinId: protein, functionType: ft, functionId: function, goAspect: aspect,
            evidenceType: type, source: source, rawScore: raw, calibratedProb: raw, metadata: metadata)
    }

    /** Logistic mechanism-test model over the CURRENT contract with explicit named weights. */
    static LearnedContextModel logisticModel(Map<String, Double> w, double intercept = 0.0d,
                                             double calibA = 1.0d, double calibB = 0.0d, double threshold = 0.5d,
                                             Map locus = null) {
        def c = FeatureContract.CURRENT
        def weights = c.names().collect { (w[it] ?: 0.0d) as double }
        Map art = [
            schema        : LearnedContextModel.SCHEMA,
            model_id      : 'mechanism-test-logistic',
            intended_use  : 'mechanism_test',
            contract      : [id: c.id, fingerprint: c.fingerprint(), features: c.names()],
            kind          : 'logistic',
            normalization : [mean: c.names().collect { 0.0d }, scale: c.names().collect { 1.0d }],
            logistic      : [weights: weights, intercept: intercept],
            calibration   : [method: 'platt', a: calibA, b: calibB],
            display_threshold: [value: threshold, advisory_only: true],
            locus_model   : locus ?: [trained: false, reason: 'test artifact without locus layer'],
            training      : [frozen_before_test: true, synthetic: true],
        ]
        mechanismModel(JsonOutput.toJson(art))
    }

    /**
     * The ONLY way to obtain a mechanism-test (non-production) model: lives in the test sources, so
     * no CLI/service code can reach it. Production loaders refuse these artifacts.
     */
    static LearnedContextModel mechanismModel(String json, FeatureContract c = FeatureContract.CURRENT) {
        LearnedContextModel.parseArtifact(json, c, true)
    }

    static LearnedContextModel mechanismFixture(String resourceName) {
        mechanismModel(LearnedTestSupport.getResourceAsStream("/learned/${resourceName}").getText('UTF-8'))
    }

    /** A stage that accepts mechanism-test models; the production class rejects them. */
    static class MechanismStage extends LearnedIntegrationStage {
        MechanismStage(LearnedStageSettings s) { super(s) }
        @Override
        protected void checkModelUse(LearnedContextModel model) { /* test API: mechanism models allowed */ }
    }

    static MechanismStage stage(LearnedContextModel model, Map overrides = [:]) {
        def s = new LearnedStageSettings(model: model)
        overrides.each { k, v -> s."$k" = v }
        new MechanismStage(s)
    }

    static Map locusLayer(double scale, double nullBias, List<Double> nullWeights = [0, 0, 0, 0, 0, 0]) {
        [trained: true,
         contract: [id: FeatureContract.LOCUS_GROUP.id, fingerprint: FeatureContract.LOCUS_GROUP.fingerprint(),
                    features: FeatureContract.LOCUS_GROUP.names()],
         scale: scale, null_bias: nullBias, null_weights: nullWeights]
    }

    /** Model with the weights used across the mechanism tests. */
    static LearnedContextModel standardModel(Map locus = null) {
        logisticModel([
            'seq.group.homology.max_raw': 4.0d,
            'seq.src.diamond.max_raw'   : 0.0d,
            'ctx.nbr.support_max'       : 1.5d,
            'ctx.nbr.conflict_max'      : -3.0d,      // unfavourable context must be able to LOWER the score
            'ctx.competitor.max_raw'    : -1.0d,
        ], -2.5d, 1.0d, 0.0d, 0.5d, locus)
    }
}
