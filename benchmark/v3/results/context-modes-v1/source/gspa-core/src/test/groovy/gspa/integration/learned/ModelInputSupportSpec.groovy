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

    Map artifact(CandidateFeatures f) {
        Map a = new JsonSlurper().parseText(getClass().getResourceAsStream('/learned/mechanism_fixture_logistic.json').text) as Map
        a.intended_use = 'production'; a.training.synthetic = false
        def names = ModelInputSupport.featureNames(f.contract.names())
        a.input_support = [schema: 'gspa.model-input-support/1', available: true, features: names,
            profiles: [[sources: f.sequenceSources, flags: names.collect { (int) f.get(it) },
                        train_labels: [2,2], calibration_labels: [2,2]]]]
        a
    }

    def "quality observations override an otherwise supported production input profile"() {
        given:
        def raw = claim('p1', 'EC:1.1.1.1', .8, 'tool_a')
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
        def fs = extractor.extractAll([claim('p1','EC:1.1.1.1',.8,' ToolA '), claim('p1','EC:1.1.1.1',.6,'toolB'),
            claim('p1','EC:1.1.1.1',Double.NaN,'bad'), claim('p1','EC:1.1.1.1',.2,'toola')], state)
        then:
        fs[0].sequenceSources == ['toola','toolb']
        fs[0].get('seq.src.other.reported') == 1
    }

    def "production model rejects unrelated pooled source and missing legacy support"() {
        given:
        def f = extractor.extractAll([claim('p1','EC:1.1.1.1',.8,'tool_a')], state)[0]
        def other = extractor.extractAll([claim('p1','EC:1.1.1.1',.8,'tool_b')], state)[0]
        Map a = artifact(f)
        def model = LearnedContextModel.loadJson(JsonOutput.toJson(a))
        expect:
        model.score(f).probability != null
        model.unsupportedReason(other.values, other.sequenceSources).startsWith('unseen')
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
        def claims = [claim('p1','EC:1.1.1.1',.8,'tool_a'), claim('p2','EC:1.1.1.1',.8,'tool_b')]
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
        def f = extractor.extractAll([claim('p1','EC:1.1.1.1',.8)],state)[0]
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
        m.unsupportedReason(f.values,f.sequenceSources).contains('both_classes')
    }

    def 'legacy models abstain on an unseen context mode even with identical numeric features'() {
        given:
        def f = extractor.extractAll([claim('p1','EC:1.1.1.1',.8,'tool_a')], state)[0]
        def other = new CandidateFeatures(f.proteinId, f.functionType, f.functionId, f.role, f.contract,
            f.values, f.provenance, f.sequenceSources, [], gspa.integration.cohort.ContextMode.PATHWAY)
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
        def f = extractor.extractAll([claim('p1','EC:1.1.1.1',.8,'tool_a')], state)[0]
        Map a = artifact(f)
        a.input_support.schema = 'gspa.model-input-support/2'
        Map profile = a.input_support.profiles[0]
        profile.cohort_context_mode = 'UNIVERSAL'
        a.input_support.profiles << (profile + [cohort_context_mode:'PATHWAY', calibration_labels:[0,2]])
        def model = LearnedContextModel.loadJson(JsonOutput.toJson(a))
        expect:
        model.unsupportedReason(f.values, f.sequenceSources, 'UNIVERSAL') == null
        model.unsupportedReason(f.values, f.sequenceSources, 'PATHWAY') == 'profile_lacks_both_classes_in_train_and_calibration'
        model.unsupportedReason(f.values, f.sequenceSources, 'ORGANISM') == 'unseen_reaction_context_mode'
        when:
        a.input_support.profiles[1].cohort_context_mode = 'TYPO'
        LearnedContextModel.loadJson(JsonOutput.toJson(a))
        then:
        thrown(ModelUnavailableException)
    }
}
