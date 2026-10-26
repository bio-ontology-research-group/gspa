package gspa.integration.learned

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import gspa.integration.EvidenceType
import gspa.integration.IntegrationState
import gspa.model.AnnotationType
import gspa.model.Genome
import spock.lang.Specification

import static gspa.integration.learned.LearnedTestSupport.*

/**
 * Inference mechanics of {@link LearnedContextModel}. Two kinds of evidence are kept apart:
 * <ul>
 *   <li>mechanism tests with hand-set weights (marked mechanism_test): scores can fall as well as
 *       rise, null outcomes are non-zero, roles stay distinct, etc.</li>
 *   <li>JVM/Python numerical agreement against a fixture exported by the Python trainer from
 *       PLANTED synthetic data. Neither says anything about accuracy on real genomes.</li>
 * </ul>
 */
class LearnedContextModelSpec extends Specification {

    def extractor = new ContextFeatureExtractor()
    def state = new IntegrationState(new Genome(id: 'g'))

    private List<CandidateFeatures> feats(List claims, ExternalContext ctx = null) { extractor.extractAll(claims, state, ctx) }

    private static InputStream res(String name) { LearnedContextModelSpec.getResourceAsStream("/learned/${name}") }

    // ------------------------------------------------------------------ learned scores move both ways

    def "unfavourable context lowers the learned probability; favourable raises it; no fixed floor"() {
        given:
        def model = standardModel()
        def claims = [claim('p1', 'EC:1.1.1.1', 0.8)]
        def ctxOf = { double support, double conflict ->
            new ExternalContext(neighborhood: ['p1|EC|EC:1.1.1.1': new ExternalContext.Neighborhood(
                supportMax: support, nSupporting: 2, conflictMax: conflict)])
        }

        when:
        double sequenceOnly = model.score(feats(claims)[0]).probability
        double favourable = model.score(feats(claims, ctxOf(0.9d, 0.0d))[0]).probability
        double unfavourable = model.score(feats(claims, ctxOf(0.0d, 0.95d))[0]).probability
        double severe = model.score(feats(claims, ctxOf(0.0d, 5.0d))[0]).probability

        then:
        favourable > sequenceOnly
        unfavourable < sequenceOnly
        severe < unfavourable                    // keeps falling with more conflict: no fixed floor
        severe < 1e-3
    }

    def "a weak true candidate is scored and returned, flagged below the advisory threshold, not dropped"() {
        given:
        def model = standardModel()
        def weak = feats([claim('p1', 'EC:1.1.1.1', 0.05)])

        when:
        def scores = model.scoreAll(weak)

        then:
        scores.size() == 1
        scores[0].probability > 0
        !scores[0].aboveDisplayThreshold
        scores[0].displayThreshold == 0.5d
    }

    def "an unrelated candidate pool does not change an existing candidate's score"() {
        given:
        def model = standardModel()
        def base = [claim('p1', 'EC:1.1.1.1', 0.7)]
        def noise = (1..100).collect { claim("n$it".toString(), "EC:8.8.8.$it".toString(), 0.9) }

        expect:
        model.score(feats(base)[0]).probability == model.score(feats(base + noise).find { it.proteinId == 'p1' }).probability
    }

    def "repeated identical hits give the same probability as one hit"() {
        given:
        def model = standardModel()

        expect:
        model.score(feats([claim('p1', 'EC:1.1.1.1', 0.7)])[0]).probability ==
            model.score(feats((1..20).collect { claim('p1', 'EC:1.1.1.1', 0.7) })[0]).probability
    }

    def "complex membership and catalytic function are scored as separate candidates of the same protein"() {
        given:
        def st = new IntegrationState(new Genome(id: 'g'))
        st.complexClassification['GO:0005840'] = 'n'
        def all = extractor.extractAll([
            claim('p1', 'GO:0005840', 0.9, 'interproscan', EvidenceType.SEQUENCE_DOMAIN, AnnotationType.GO, [:], 'CC'),
            claim('p1', 'EC:2.7.7.6', 0.4)], st)
        def model = logisticModel(['seq.group.homology.max_raw': 4.0d, 'role.complex_subunit': -1.0d, 'role.catalytic': 0.5d], -2.0d)

        when:
        def scores = model.scoreAll(all).collectEntries { [it.features.role.name(), it] }

        then: "high complex-subunit probability does not transfer to the catalytic candidate"
        scores.keySet() == ['COMPLEX_SUBUNIT', 'CATALYTIC'] as Set
        scores['CATALYTIC'].candidateKey() != scores['COMPLEX_SUBUNIT'].candidateKey()
        scores['CATALYTIC'].probability != scores['COMPLEX_SUBUNIT'].probability
    }

    def "per-candidate probabilities are independent and may all be high (no cross-candidate normalisation)"() {
        given:
        def model = standardModel()
        def scores = model.scoreAll(feats([claim('p1', 'EC:1.1.1.1', 0.99), claim('p2', 'EC:1.1.1.1', 0.99),
                                           claim('p1', 'EC:2.2.2.2', 0.99)]))

        expect:
        scores.size() == 3
        scores.every { it.probability > 0.75d }
        scores*.probability.sum() > 1.0d           // would be impossible under a forced mutually exclusive softmax
    }

    // ------------------------------------------------------------------ null / locus outcome

    def "the null outcome is non-zero, a singleton is never normalised to 1, and group sums to 1"() {
        given:
        def model = standardModel(locusLayer(1.0d, 0.0d))
        def single = model.scoreAll(feats([claim('p1', 'EC:1.1.1.1', 0.99)]))
        def trio = model.scoreAll(feats([claim('p1', 'EC:1.1.1.1', 0.9), claim('p2', 'EC:1.1.1.1', 0.5), claim('p3', 'EC:1.1.1.1', 0.1)]))

        when:
        def a1 = model.assignLocus(single)
        def a3 = model.assignLocus(trio)

        then:
        a1.nullProbability > 0
        a1.candidateProbabilities.values().first() < 1.0d
        Math.abs(a1.nullProbability + a1.candidateProbabilities.values().sum() - 1.0d) < 1e-12
        a3.nullProbability > 0
        Math.abs(a3.nullProbability + a3.candidateProbabilities.values().sum() - 1.0d) < 1e-12
        a3.candidateProbabilities['p1|EC|EC:1.1.1.1|CATALYTIC'] > a3.candidateProbabilities['p3|EC|EC:1.1.1.1|CATALYTIC']
    }

    def "when the supplied candidates are all weak the learned null outcome takes the mass"() {
        given:
        def model = standardModel(locusLayer(1.0d, 0.0d))
        def weak = model.assignLocus(model.scoreAll(feats([claim('p1', 'EC:1.1.1.1', 0.02), claim('p2', 'EC:1.1.1.1', 0.01)])))
        def strong = model.assignLocus(model.scoreAll(feats([claim('p1', 'EC:1.1.1.1', 0.99), claim('p2', 'EC:1.1.1.1', 0.01)])))

        expect:
        weak.nullProbability > 0.5d
        strong.nullProbability < weak.nullProbability
    }

    def "no trained null layer: locus assignment fails instead of inventing a default null weight"() {
        given:
        def model = standardModel()           // locus_model.trained == false

        when:
        model.assignLocus(model.scoreAll(feats([claim('p1', 'EC:1.1.1.1', 0.9)])))

        then:
        def e = thrown(ModelUnavailableException)
        e.message.contains('no trained null-outcome layer')
        !model.hasLocusLayer()
    }

    def "a locus group must share one function and one role, with distinct proteins"() {
        given:
        def model = standardModel(locusLayer(1.0d, 0.0d))
        def mixed = model.scoreAll(feats([claim('p1', 'EC:1.1.1.1', 0.9), claim('p2', 'EC:2.2.2.2', 0.9)]))
        def dup = model.scoreAll(feats([claim('p1', 'EC:1.1.1.1', 0.9)])) * 2

        when:
        model.assignLocus(mixed)

        then:
        thrown(IllegalArgumentException)

        when:
        model.assignLocus(dup)

        then:
        thrown(IllegalArgumentException)
    }

    // ------------------------------------------------------------------ loading rules

    def "loading refuses: missing artifact, mechanism fixtures by default, contract drift, unfrozen artifacts"() {
        given:
        def good = new JsonSlurper().parse(res('mechanism_fixture_logistic.json')) as Map

        when: "no file"
        LearnedContextModel.load(new File('/nonexistent/model.json'))

        then:
        thrown(ModelUnavailableException)

        when: "a mechanism-test fixture in a production load"
        LearnedContextModel.load(res('mechanism_fixture_logistic.json'))

        then:
        def e = thrown(ModelUnavailableException)
        e.message.contains('mechanism-test fixture')

        when: "feature list drifted"
        def drift = new LinkedHashMap(good); drift.contract = good.contract + [fingerprint: 'deadbeef']
        LearnedContextModel.loadJson(JsonOutput.toJson(drift), FeatureContract.CURRENT, true)

        then:
        thrown(ModelUnavailableException)

        when: "not declared frozen"
        def thawed = new LinkedHashMap(good); thawed.training = good.training + [frozen_before_test: false]
        LearnedContextModel.loadJson(JsonOutput.toJson(thawed), FeatureContract.CURRENT, true)

        then:
        thrown(ModelUnavailableException)

        when: "wrong schema"
        LearnedContextModel.loadJson('{"schema": "other/1"}')

        then:
        thrown(ModelUnavailableException)

        when: "malformed tree topology"
        def gb = new JsonSlurper().parse(res('mechanism_fixture_gbdt.json')) as Map
        gb.gbdt.trees[0].left[0] = 0
        LearnedContextModel.loadJson(JsonOutput.toJson(gb), FeatureContract.CURRENT, true)

        then:
        thrown(ModelUnavailableException)
    }

    def "scoring rejects wrong-width or non-finite vectors"() {
        given:
        def model = standardModel()

        when:
        model.margin(new double[3])

        then:
        thrown(IllegalArgumentException)

        when:
        def x = new double[FeatureContract.CURRENT.size()]
        x[0] = Double.NaN
        model.margin(x)

        then:
        thrown(IllegalArgumentException)
    }

    // ------------------------------------------------------------------ JVM <-> Python agreement

    def "JVM inference reproduces the Python reference to 1e-12 (logistic and trees, incl. locus layer)"() {
        given:
        def vectors = new JsonSlurper().parse(res('mechanism_fixture_vectors.json')) as Map
        def names = FeatureContract.CURRENT.names()
        def vec = { Map sparse ->
            double[] x = new double[names.size()]
            sparse.each { k, v -> x[names.indexOf(k as String)] = (v as Number).doubleValue() }
            x
        }
        def model = LearnedContextModel.load(res("mechanism_fixture_${kind}.json"), FeatureContract.CURRENT, true)

        expect: "margins and calibrated probabilities"
        vectors.expected[kind].every { name, exp ->
            Math.abs(model.margin(vec(vectors.vectors[name])) - (exp.margin as double)) < 1e-12 &&
                Math.abs(model.probability(vec(vectors.vectors[name])) - (exp.probability as double)) < 1e-12
        }

        and: "locus assignments with the explicit null outcome"
        def cases = vectors.locus_cases.findAll { it.model == kind }
        !cases.isEmpty()
        cases.every { c ->
            def group = c.members.withIndex().collect { String m, int i ->
                def f = new CandidateFeatures("p${i}".toString(), AnnotationType.EC, 'EC:1.1.1.1', CandidateRole.CATALYTIC,
                    FeatureContract.CURRENT, vec(vectors.vectors[m]))
                model.score(f)
            }
            def got = model.assignLocus(group)
            Math.abs(got.nullProbability - (c.null_probability as double)) < 1e-12 &&
                (0..<group.size()).every { i ->
                    Math.abs(got.candidateProbabilities["p${i}|EC|EC:1.1.1.1|CATALYTIC".toString()] - (c.candidate_probabilities[i] as double)) < 1e-12
                }
        }

        where:
        kind << ['logistic', 'gbdt']
    }

    def "the exported synthetic fixtures show the mechanisms (planted structure), without claiming accuracy"() {
        given:
        def vectors = new JsonSlurper().parse(res('mechanism_fixture_vectors.json')) as Map
        def exp = vectors.expected[kind]

        expect: "unfavourable context lowers, favourable raises, relative to the same sequence evidence"
        exp.strong_sequence_unfavourable_context.probability < exp.strong_sequence.probability
        exp.strong_sequence_favourable_context.probability > exp.strong_sequence.probability

        and: "null outcome is learned and non-zero; a singleton is not forced to 1"
        def single = vectors.locus_cases.find { it.model == kind && it.case == 'singleton_strong' }
        single.null_probability > 0
        single.candidate_probabilities[0] < 1.0d

        and: "fixture is labelled synthetic and not loadable as production"
        def art = new JsonSlurper().parse(res("mechanism_fixture_${kind}.json")) as Map
        art.intended_use == 'mechanism_test'
        art.training.synthetic == true

        where:
        kind << ['logistic', 'gbdt']
    }
}
