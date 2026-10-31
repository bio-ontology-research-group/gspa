package gspa.integration.learned

import groovy.json.JsonSlurper
import gspa.integration.*
import gspa.integration.contextbridge.RawBaselineSnapshot
import gspa.model.*
import spock.lang.Specification
import static gspa.integration.learned.LearnedTestSupport.*

class QualityLearnedSpec extends Specification {
    def state = new IntegrationState(new Genome(id: 'g'))
    def extractor = new ContextFeatureExtractor()

    EvidenceClaim flag(String p = 'p1') {
        new ClaimExtractor().fromAnnotation(p, new Annotation(type: AnnotationType.PSEUDOGENE,
            value: 'ANF00021', source: 'antifam', score: 1d, metadata: [name: 'Spurious_ORF_21']))
    }

    def 'quality is retained separately without creating candidates or changing numeric functional features'() {
        given:
        def claims = [claim('p1', 'EC:1.1.1.1', .8d), claim('p2', 'EC:1.1.1.1', .6d)]
        def baseline = extractor.extractAll(claims, state)
        when:
        def features = extractor.extractAll(claims + [flag(), flag('quality_only')], state)
        then:
        features*.candidateKey() == baseline*.candidateKey()
        features*.asMap() == baseline*.asMap()
        features[0].qualityEvidence*.model_id == ['ANF00021']
        features[1].qualityEvidence == []
        features*.sequenceSources == baseline*.sequenceSources
        extractor.prepare([flag()], state, new ExternalContext(candidates: ['p1|PSEUDOGENE|ANF00021'] as Set)).functionKeys().isEmpty()
    }

    def 'learned stage and feature export abstain consistently without retracting a flagged function or ranking a partial locus'() {
        given:
        def claims = [claim('p1', 'EC:1.1.1.1', .8d), claim('p2', 'EC:1.1.1.1', .6d), flag()]
        def integrated = new IterativeRefiner(new EvidenceCombiner()).refine(claims, state)
        def original = integrated.provenance['p1|EC|EC:1.1.1.1'].finalProbability
        def model = standardModel(locusLayer(1d, 0d))
        when:
        def report = stage(model).run(claims, state, integrated)
        def writer = new StringWriter()
        new FeatureExporter(model).export(claims, state, null, writer)
        def rows = writer.toString().readLines().collect { new JsonSlurper().parseText(it) }
        then:
        report.candidates.size() == 2
        report.candidates[0].unsupportedReason == GeneModelQuality.UNSUPPORTED_REASON
        report.candidates[0].probability == null
        report.candidates[0].qualityEvidence.size() == 1
        report.candidates[1].scoreStatus == 'scored'
        report.candidates*.locusProbability.every { it == null }
        integrated.provenance['p1|EC|EC:1.1.1.1'].finalProbability == original
        rows[0].model.unsupported_reason == GeneModelQuality.UNSUPPORTED_REASON
        rows[0].model.probability == null
        rows[0].quality_evidence*.model_id == ['ANF00021']
        rows[1].model.score_status == 'scored'
        when:
        model.score(extractor.extractAll(claims, state)[0])
        then:
        thrown(IllegalArgumentException)
    }

    def 'quality cannot be disguised as an EC to enter Noisy OR or a cohort anchor'() {
        given:
        def q = flag()
        q.functionType = AnnotationType.EC
        q.functionId = '1.1.1.1'
        when:
        def snap = RawBaselineSnapshot.forGenome('g', [q])
        then:
        snap.calls.isEmpty()
        new EvidenceCombiner().combineLikelihood([q]) == -12d
        extractor.extractAll([q], state).isEmpty()
    }

    def 'scan coverage is never inferred from absent hits: #status'() {
        when:
        def report = GeneModelQuality.report([], status)
        then:
        report.antifam_scan_status == (status ?: 'UNKNOWN')
        report.antifam_scan_complete == complete
        report.observations.isEmpty()
        where:
        status         | complete
        null           | false
        'SUCCESS_ZERO' | true
        'SUCCESS'      | true
        'PARTIAL'      | false
        'FAILED'       | false
        'UNAVAILABLE'  | false
    }
}
