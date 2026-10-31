package gspa.integration.learned

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import gspa.integration.*
import gspa.model.*
import spock.lang.Specification
import static gspa.integration.learned.LearnedTestSupport.*

/** Fabricated artifact metadata tests the production gate, never biological accuracy. */
class ModelInputSupportSpec extends Specification {
    def state = new IntegrationState(new Genome(id: 'g'))
    def extractor = new ContextFeatureExtractor()


    static Map testContract() {
        [schema: 'gspa.predictor-contract/1', tool: 'synthetic-test', tool_version: '1',
         tool_sha256: 'a'*64, reference_sha256: 'b'*64, parameters_sha256: 'c'*64,
         adapter_sha256: 'd'*64, score_model: 'synthetic-test/1']
    }

    static EvidenceClaim versionedClaim(String p, String f, double score, String source = 'tool_a') {
        def c = claim(p, f, score, source)
        c.metadata = [predictor_contract: testContract()]
        c
    }

    Map artifact(CandidateFeatures f) {
        Map a = new JsonSlurper().parseText(getClass().getResourceAsStream('/learned/mechanism_fixture_logistic.json').text) as Map
        a.intended_use = 'production'; a.training.synthetic = false
        def names = ModelInputSupport.featureNames(f.contract.names())
        a.input_support = [schema: 'gspa.model-input-support/3', available: true, features: names,
            score_features: ModelInputSupport.scoreNames(f.contract.names()),
            profiles: [[sources: f.sequenceSources, sequence_contracts: f.sequenceContracts, cohort_context_mode: 'UNIVERSAL',
                        lower: ModelInputSupport.scoreNames(f.contract.names()).collect { f.get(it) },
                        upper: ModelInputSupport.scoreNames(f.contract.names()).collect { f.get(it) }, flags: names.collect { (int) f.get(it) },
                        train_labels: [2,2], calibration_labels: [2,2]]]]
        a
    }

    def "quality observations override an otherwise supported production input profile"() {
        given:
        def raw = versionedClaim('p1', 'EC:1.1.1.1', .8, 'tool_a')
        def baseline = extractor.extractAll([raw], state)[0]
        // Fabricated coverage metadata tests the guard, not model accuracy.
        def model = LearnedContextModel.loadJson(JsonOutput.toJson(artifact(baseline)))
        def quality = new ClaimExtractor().fromAnnotation('p1', new Annotation(type: AnnotationType.PSEUDOGENE,
            value: 'ANF00021', source: 'antifam', score: 1d))
        def flagged = extractor.extractAll([raw, quality], state)[0]
        expect:
        model.score(baseline).probability != null
        model.unsupportedReason(flagged) == GeneModelQuality.UNSUPPORTED_REASON
        when:
        model.score(flagged)
        then:
        thrown(IllegalArgumentException)
    }

    def "pooled sources retain distinct normalized identities and ignore unusable claims"() {
        when:
        def fs = extractor.extractAll([versionedClaim('p1','EC:1.1.1.1',.8,' ToolA '), versionedClaim('p1','EC:1.1.1.1',.6,'toolB'),
            versionedClaim('p1','EC:1.1.1.1',Double.NaN,'bad'), versionedClaim('p1','EC:1.1.1.1',.2,'toola')], state)
        then:
        fs[0].sequenceSources == ['toola','toolb']
        fs[0].get('seq.src.other.reported') == 1
    }

    def "production model rejects unrelated pooled source and missing legacy support"() {
        given:
        def f = extractor.extractAll([versionedClaim('p1','EC:1.1.1.1',.8,'tool_a')], state)[0]
        def other = extractor.extractAll([versionedClaim('p1','EC:1.1.1.1',.8,'tool_b')], state)[0]
        Map a = artifact(f)
        def model = LearnedContextModel.loadJson(JsonOutput.toJson(a))
        expect:
        model.score(f).probability != null
        model.unsupportedReason(other.values, other.sequenceSources, 'UNIVERSAL', other.sequenceContracts).startsWith('unseen')
        when:
        model.score(other)
        then:
        thrown(IllegalArgumentException)
        when:
        a.remove('input_support')
        LearnedContextModel.loadJson(JsonOutput.toJson(a)).score(f)
        then:
        thrown(IllegalArgumentException)
    }

    def "mixed stage preserves unsupported posterior and never ranks a partial locus"() {
        given:
        def claims = [versionedClaim('p1','EC:1.1.1.1',.8,'tool_a'), versionedClaim('p2','EC:1.1.1.1',.8,'tool_b')]
        def f = extractor.extractAll(claims, state)[0]
        def model = LearnedContextModel.loadJson(JsonOutput.toJson(artifact(f)))
        def result = new IntegratedAnnotationSet()
        claims.each { c ->
            String key = "${c.proteinId}|EC|EC:1.1.1.1"
            result.put(new ClaimProvenance(functionKey:key, proteinId:c.proteinId, finalProbability:.7,
                supportingClaims:[c]), new Annotation(type:AnnotationType.EC,value:'EC:1.1.1.1',score:.7))
            state.set(key, Math.log(.7/.3))
        }
        state.pinnedFloors['p2|EC|EC:1.1.1.1'] = .5
        when:
        def report = new LearnedIntegrationStage(new LearnedStageSettings(model:model)).run(claims,state,result)
        then:
        report.count('rescored') == 1
        report.count('unchanged_unsupported') == 1
        result.provenance['p2|EC|EC:1.1.1.1'].finalProbability == .7
        state.pinnedFloors['p2|EC|EC:1.1.1.1'] == .5
        report.candidates.find { it.proteinId == 'p2' }.probability == null
        report.candidates.every { it.locusNullProbability == null }
        and:
        def w = new StringWriter()
        new FeatureExporter(model).export(claims,state,null,w)
        def rows = w.toString().readLines().collect { new JsonSlurper().parseText(it) }
        rows.find { it.protein_id == 'p2' }.model.score_status == 'unsupported'
        rows.find { it.protein_id == 'p2' }.model.probability == null
        rows.find { it.protein_id == 'p1' }.model.probability != null
    }

    def "support manifest cannot omit role indicators or invent class coverage"() {
        given:
        def f = extractor.extractAll([versionedClaim('p1','EC:1.1.1.1',.8)],state)[0]
        Map a = artifact(f)
        a.input_support.features = ['fn.type_ec']
        when:
        LearnedContextModel.loadJson(JsonOutput.toJson(a))
        then:
        thrown(ModelUnavailableException)
        when:
        a = artifact(f); a.input_support.profiles[0].calibration_labels = [0,4]
        def m = LearnedContextModel.loadJson(JsonOutput.toJson(a))
        then:
        m.unsupportedReason(f.values,f.sequenceSources, 'UNIVERSAL', f.sequenceContracts).contains('both_classes')
    }

    def 'legacy models abstain on an unseen context mode even with identical numeric features'() {
        given:
        def f = extractor.extractAll([versionedClaim('p1','EC:1.1.1.1',.8,'tool_a')], state)[0]
        def other = new CandidateFeatures(f.proteinId, f.functionType, f.functionId, f.role, f.contract,
            f.values, f.provenance, f.sequenceSources, [], gspa.integration.cohort.ContextMode.PATHWAY, f.sequenceContracts)
        def model = LearnedContextModel.loadJson(JsonOutput.toJson(artifact(f)))
        expect:
        model.unsupportedReason(f) == null
        model.unsupportedReason(other) == 'unseen_reaction_context_mode'
        when:
        model.score(other)
        then:
        thrown(IllegalArgumentException)
    }

    def 'mode-aware profiles require class coverage separately for each mode'() {
        given:
        def f = extractor.extractAll([versionedClaim('p1','EC:1.1.1.1',.8,'tool_a')], state)[0]
        Map a = artifact(f)
        a.input_support.schema = 'gspa.model-input-support/3'
        Map profile = a.input_support.profiles[0]
        profile.cohort_context_mode = 'UNIVERSAL'
        a.input_support.profiles << (profile + [cohort_context_mode:'PATHWAY', calibration_labels:[0,2]])
        def model = LearnedContextModel.loadJson(JsonOutput.toJson(a))
        expect:
        model.unsupportedReason(f.values, f.sequenceSources, 'UNIVERSAL', f.sequenceContracts) == null
        model.unsupportedReason(f.values, f.sequenceSources, 'PATHWAY', f.sequenceContracts) == 'profile_lacks_both_classes_in_train_and_calibration'
        model.unsupportedReason(f.values, f.sequenceSources, 'ORGANISM', f.sequenceContracts) == 'unseen_reaction_context_mode'
        when:
        a.input_support.profiles[1].cohort_context_mode = 'TYPO'
        LearnedContextModel.loadJson(JsonOutput.toJson(a))
        then:
        thrown(ModelUnavailableException)
    }

    def 'different predictor identities cannot share calibration despite identical numeric vectors'() {
        given:
        def original = versionedClaim('p', 'EC:1.1.1.1', .8)
        def altered = versionedClaim('p', 'EC:1.1.1.1', .8)
        altered.metadata.predictor_contract[field] = replacement
        def f = extractor.extractAll([original], state)[0]
        def changed = extractor.extractAll([altered], state)[0]
        def model = LearnedContextModel.loadJson(JsonOutput.toJson(artifact(f)))
        expect:
        f.digest() == changed.digest()
        model.unsupportedReason(changed) == 'unseen_predictor_contract'
        where:
        field | replacement
        'tool_version' | '2'
        'tool_sha256' | 'e'*64
        'reference_sha256' | 'e'*64
        'parameters_sha256' | 'e'*64
        'adapter_sha256' | 'e'*64
        'score_model' | 'different/1'
    }

    def 'missing or malformed weaker evidence cannot hide behind a stronger identified claim'() {
        given:
        def good = versionedClaim('p', 'EC:1.1.1.1', .8)
        def weak = claim('p', 'EC:1.1.1.1', .2, 'tool_a')
        weak.metadata = [predictor_contract: raw]
        def f = extractor.extractAll([good], state)[0]
        def mixed = extractor.extractAll([good, weak, good], state)[0]
        def model = LearnedContextModel.loadJson(JsonOutput.toJson(artifact(f)))
        expect:
        f.digest() == mixed.digest()
        mixed.sequenceContracts.size() == 2
        model.unsupportedReason(mixed) == 'candidate_predictor_contract_unavailable'
        where:
        raw << [null, [schema:'wrong'], testContract() + [tool_sha256:'broken']]
    }

    def 'sequence scores outside the frozen profile range abstain including a partial locus'() {
        given:
        def f = extractor.extractAll([versionedClaim('p', 'EC:1.1.1.1', .8)], state)[0]
        def changed = extractor.extractAll([versionedClaim('p', 'EC:1.1.1.1', .9)], state)[0]
        def model = LearnedContextModel.loadJson(JsonOutput.toJson(artifact(f)))
        expect:
        model.unsupportedReason(f) == null
        model.unsupportedReason(changed) == 'sequence_score_outside_training_calibration_range'
        when:
        model.score(changed)
        then:
        thrown(IllegalArgumentException)
    }

    def 'legacy support cannot certify predictor contracts'() {
        given:
        def f = extractor.extractAll([versionedClaim('p', 'EC:1.1.1.1', .8)], state)[0]
        def a = artifact(f)
        a.input_support.schema = schema
        def model = LearnedContextModel.loadJson(JsonOutput.toJson(a))
        expect:
        model.unsupportedReason(f) == 'model_predictor_contract_support_unavailable'
        where:
        schema << ['gspa.model-input-support/1', 'gspa.model-input-support/2']
    }

    def 'contract metadata is immutable deduplicated and survives feature export'() {
        given:
        def c = versionedClaim('p', 'EC:1.1.1.1', .8)
        def f = extractor.extractAll([c, c], state)[0]
        c.metadata.predictor_contract.tool_version = 'changed'
        expect:
        f.sequenceContracts.size() == 1
        f.sequenceContracts[0].contract.tool_version == '1'
        when:
        def writer = new StringWriter()
        new FeatureExporter().export([c],state,null,writer)
        def row = new JsonSlurper().parseText(writer.toString().trim())
        then:
        row.sequence_contracts[0].contract.tool_version == 'changed'
    }

    def 'corrupt model ranges cannot silently disable the gate'() {
        given:
        def f = extractor.extractAll([versionedClaim('p', 'EC:1.1.1.1', .8)], state)[0]
        def a = artifact(f)
        a.input_support.profiles[0].lower = []
        when:
        LearnedContextModel.loadJson(JsonOutput.toJson(a))
        then:
        thrown(ModelUnavailableException)
    }
}
