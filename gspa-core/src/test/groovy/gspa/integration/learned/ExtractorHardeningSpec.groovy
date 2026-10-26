package gspa.integration.learned

import gspa.integration.EvidenceClaim
import gspa.integration.EvidenceType
import gspa.integration.IntegrationState
import gspa.model.AnnotationType
import gspa.model.Genome
import spock.lang.Specification

import static gspa.integration.learned.LearnedTestSupport.claim

/**
 * Regression tests for defects found when auditing the foundations (scales, missingness flags,
 * complex classification, cost). Mechanism tests only; no accuracy claims.
 */
class ExtractorHardeningSpec extends Specification {

    def extractor = new ContextFeatureExtractor()
    def state = new IntegrationState(new Genome(id: 'g'))

    private CandidateFeatures one(List<EvidenceClaim> claims, ExternalContext ctx = null, IntegrationState st = state) {
        def all = extractor.extractAll(claims, st, ctx)
        assert all.size() == 1
        all[0]
    }

    def "ambiguous complex terms are distinct from homomeric and heteromeric ones"() {
        given:
        def go = { String id, String aspect -> claim('p1', id, 0.9, 'interproscan', EvidenceType.SEQUENCE_DOMAIN, AnnotationType.GO, [:], aspect) }
        state.complexClassification['GO:0000001'] = 'h'
        state.complexClassification['GO:0000002'] = 'n'
        state.complexClassification['GO:0000003'] = 'a'

        expect:
        def h = one([go('GO:0000001', 'CC')])
        def n = one([go('GO:0000002', 'CC')])
        def a = one([go('GO:0000003', 'CC')])
        [h, n, a].every { it.get('fn.in_complex_term') == 1 }
        [h.get('fn.complex_heteromeric'), n.get('fn.complex_heteromeric'), a.get('fn.complex_heteromeric')] == [0d, 1d, 0d]
        [h.get('fn.complex_ambiguous'), n.get('fn.complex_ambiguous'), a.get('fn.complex_ambiguous')] == [0d, 0d, 1d]
        h.values != a.values
        [h, n, a].every { it.role == CandidateRole.COMPLEX_SUBUNIT }
    }

    def "an unknown complex classification value is an error, not silently treated as a complex"() {
        given:
        state.complexClassification['GO:0000001'] = 'x'

        when:
        extractor.extractAll([claim('p1', 'GO:0000001', 0.9, 'diamond', EvidenceType.SEQUENCE_SIMILARITY, AnnotationType.GO)], state)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('GO:0000001')
    }

    def "GO aspect: GAF letters accepted, conflicting claim labels are not trusted, invalid labels are unknown"() {
        given:
        def go = { String aspect, String src = 'diamond' ->
            claim('p1', 'GO:0016301', 0.8, src, EvidenceType.SEQUENCE_SIMILARITY, AnnotationType.GO, [:], aspect)
        }

        expect:
        one([go('F')]).get('fn.go_aspect_mf') == 1
        one([go('mf')]).get('fn.go_aspect_mf') == 1
        one([go('weird')]).get('fn.go_aspect_unknown') == 1
        one([go('MF'), go('BP', 'interproscan')]).get('fn.go_aspect_unknown') == 1
        one([go('MF'), go('BP', 'interproscan')]).provenance.any { it.contains('disagree') }
        one([go('MF'), go('MF', 'interproscan')]).get('fn.go_aspect_mf') == 1
    }

    def "raw scores outside [0,1] are ignored with a note, like non-finite ones, and do not become evidence"() {
        when:
        def f = one([claim('p1', 'EC:1.1.1.1', 7.5d), claim('p1', 'EC:1.1.1.1', -0.2d), claim('p1', 'EC:1.1.1.1', 0.4d, 'interproscan', EvidenceType.SEQUENCE_DOMAIN)])

        then:
        f.get('seq.src.diamond.reported') == 0
        f.get('seq.group.homology.max_raw') == 0.4d
        f.provenance.any { it.contains('out-of-[0,1]') }

        and: "an out-of-range hit does not make a neighbour look supported either"
        def comp = extractor.extractAll([claim('p1', 'EC:1.1.1.1', 0.5d), claim('p2', 'EC:1.1.1.1', 3.0d)], state)
        comp.size() == 2
    }

    def "percent identity outside 0..100 is unavailable, not a feature value"() {
        expect:
        one([claim('p1', 'EC:1.1.1.1', 0.8, 'diamond', EvidenceType.SEQUENCE_SIMILARITY, AnnotationType.EC, [pident: 250.0d])]).get('seq.best_pident_available') == 0
        one([claim('p1', 'EC:1.1.1.1', 0.8, 'diamond', EvidenceType.SEQUENCE_SIMILARITY, AnnotationType.EC, [pident: 91.5d])]).get('seq.best_pident') == 0.915d
    }

    def "pipeline predictor names map to their own column; pooled unknown tools are named in the provenance"() {
        expect:
        ContextFeatureExtractor.normalizeSource('esm2-deepgoplus') == 'esm2-go'
        ContextFeatureExtractor.normalizeSource('eggNOG') == 'eggnog-mapper'
        ContextFeatureExtractor.normalizeSource('DIAMOND') == 'diamond'
        ContextFeatureExtractor.normalizeSource('esm2-centroid') == 'other'

        when:
        def f = one([claim('p1', 'EC:1.1.1.1', 0.6, 'esm2-centroid', EvidenceType.PROTEIN_LM_EMBEDDING),
                     claim('p1', 'EC:1.1.1.1', 0.3, 'newtool', EvidenceType.SEQUENCE_DEEPLEARNING)])

        then:
        f.get('seq.src.other.max_raw') == 0.6d
        f.get('seq.n_sources_reporting') == 1
        f.provenance.any { it.contains('esm2-centroid') && it.contains('newtool') }
    }

    def "an empty predictor-status map is 'not supplied', never 'supplied and clean'"() {
        expect:
        one([claim('p1', 'EC:1.1.1.1', 0.5)], new ExternalContext(predictorStatus: [:])).get('seq.cov.status_available') == 0
        one([claim('p1', 'EC:1.1.1.1', 0.5)], new ExternalContext(predictorStatus: [diamond: 'SUCCESS'])).get('seq.cov.status_available') == 1
    }

    def "FBA with zero tested conditions is unavailable, not a measured value"() {
        when:
        def f = one([claim('p1', 'EC:1.1.1.1', 0.5)], new ExternalContext(fba: ['p1|EC|EC:1.1.1.1': new ExternalContext.FbaConditional(nConditions: 0, requiredFraction: 0.0d)]))

        then:
        f.get('fba.available') == 0
        f.get('fba.suppressed_objective_forced') == 0
    }

    def "inconsistent cohort counts are rejected rather than producing fractions above one"() {
        when:
        extractor.extractAll([claim('p1', 'EC:1.1.1.1', 0.5)], state, new ExternalContext(
            cohort: ['p1|EC|EC:1.1.1.1': new ExternalContext.Cohort(nGenomes: 5, nWithOrthogroup: 9)]))

        then:
        thrown(IllegalArgumentException)
    }

    def "claim order does not change a candidate's vector or roles"() {
        given:
        def claims = [claim('p1', 'EC:1.1.1.1', 0.5), claim('p1', 'EC:1.1.1.1', 0.8, 'interproscan', EvidenceType.SEQUENCE_DOMAIN),
                      claim('p1', 'EC:1.1.1.1', 0.2, 'foldseek', EvidenceType.STRUCTURE_SIMILARITY)]
        def ctx = new ExternalContext(roles: ['p1|EC|EC:1.1.1.1': [CandidateRole.CATALYTIC, CandidateRole.ACCESSORY]])

        expect:
        extractor.extractAll(claims, state, ctx)*.values == extractor.extractAll(claims.reverse(), state, ctx)*.values
        extractor.extractAll(claims, state, ctx)*.role == [CandidateRole.CATALYTIC, CandidateRole.ACCESSORY]
    }

    def "extraction is linear: a prepared index serves many candidates and extract() agrees with extractAll()"() {
        given:
        int n = 20_000
        def claims = (1..n).collect { claim("p$it".toString(), "EC:1.1.${it % 500}.${it}".toString(), 0.1d + (it % 9) * 0.1d) }

        when:
        long t0 = System.nanoTime()
        int seen = 0
        extractor.prepare(claims, state).forEach { CandidateFeatures f -> seen++ }
        double seconds = (System.nanoTime() - t0) / 1e9

        then: "O(claims), not O(candidates x claims): 20k candidates would need 4e8 claim visits otherwise"
        seen == n
        seconds < 20

        and:
        def all = extractor.extractAll(claims.take(50), state)
        def single = extractor.extract('p7', AnnotationType.EC, claims[6].functionId, claims.take(50), state)
        single.size() == 1
        single[0].values == all.find { it.proteinId == 'p7' }.values
    }

    def "context-only candidates carry no sequence evidence and never fake it"() {
        when:
        def prepared = extractor.prepare([], state, null, ['p9|EC|EC:3.3.3.3'])
        def cands = prepared.candidatesOf('p9|EC|EC:3.3.3.3')

        then:
        cands.size() == 1
        cands[0].get('seq.n_groups_reporting') == 0
        cands[0].get('seq.max_raw_any') == 0
        cands[0].role == CandidateRole.CATALYTIC

        when:
        extractor.prepare([], state, null, ['not-a-key'])

        then:
        thrown(IllegalArgumentException)
    }
}
