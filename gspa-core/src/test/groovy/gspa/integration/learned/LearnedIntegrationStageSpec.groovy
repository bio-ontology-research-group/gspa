package gspa.integration.learned

import gspa.integration.ClaimProvenance
import gspa.integration.EvidenceClaim
import gspa.integration.EvidenceType
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Genome
import spock.lang.Specification

import static gspa.integration.learned.LearnedTestSupport.*

/**
 * Mechanism tests of {@link LearnedIntegrationStage} with HAND-SET test models (mechanism_test,
 * accepted only through the test API). They show the stage does what its contract says; they
 * are not evidence about prediction quality.
 */
class LearnedIntegrationStageSpec extends Specification {

    static final String EC1 = 'p1|EC|EC:1.1.1.1'

    IntegrationState state = new IntegrationState(new Genome(id: 'g'))
    IntegratedAnnotationSet result = new IntegratedAnnotationSet()

    private static double logit(double p) { Math.log(p / (1d - p)) }

    /** Retain a candidate the way the integrator does: provenance + annotation + posterior. */
    private void retain(String key, double p, List<EvidenceClaim> claims = []) {
        def parts = IntegrationState.splitFunctionKey(key)
        result.put(new ClaimProvenance(functionKey: key, proteinId: parts[0], supportingClaims: claims,
            priorContributions: ['arbitrary_boost': 3.0d], likelihoodLogOdds: 1.0d, finalProbability: p),
            new Annotation(type: AnnotationType.valueOf(parts[1]), value: parts[2], source: 'integrated', evidence: 'IEA'))
        state.set(key, logit(p))
    }

    private static ExternalContext nbr(String key, double support, double conflict) {
        new ExternalContext(neighborhood: [(key): new ExternalContext.Neighborhood(supportMax: support, nSupporting: 2, conflictMax: conflict)])
    }

    private Map snapshot() {
        [prov : result.provenance.collectEntries { k, v -> [k, [v.finalProbability, v.finalLogOdds, new LinkedHashMap(v.priorContributions)]] },
         ann  : result.annotationsByFunction.collectEntries { k, v -> [k, [v.score, new LinkedHashMap(v.metadata)]] },
         annList: result.annotations.annotations.collect { [it.value, it.score] },
         post : new LinkedHashMap(state.posteriorLogOdds),
         floors: new LinkedHashMap(state.pinnedFloors),
         processed: new LinkedHashSet(result.processedProteinIds)]
    }

    // ------------------------------------------------------------------ production gate

    def "no model: explicit failure, no defaults"() {
        when:
        new LearnedIntegrationStage(new LearnedStageSettings())

        then:
        thrown(ModelUnavailableException)
    }

    def "a mechanism-test model is refused by the production stage"() {
        when:
        new LearnedIntegrationStage(new LearnedStageSettings(model: standardModel()))

        then:
        def e = thrown(ModelUnavailableException)
        e.message.contains('not production')
    }

    def "invalid typed settings are rejected"() {
        when:
        stage(standardModel(), [retractBelow: 1.5d])

        then:
        thrown(IllegalArgumentException)
    }

    // ------------------------------------------------------------------ replacement, not addition

    def "the learned probability replaces the old posterior; the old value is audit only"() {
        given:
        def model = standardModel()
        def claims = [claim('p1', 'EC:1.1.1.1', 0.8)]
        retain(EC1, 0.97d, claims)           // an inflated legacy posterior (e.g. prior-engine boosts)
        def expected = model.score(new ContextFeatureExtractor().extractAll(claims, state)[0]).probability

        when:
        def report = stage(model).run(claims, state, result)

        then:
        result.annotationsByFunction[EC1].score == expected
        result.provenanceFor('p1', AnnotationType.EC, 'EC:1.1.1.1').finalProbability == expected
        Math.abs(state.getProbability(EC1) - expected) < 1e-9
        expected != 0.97d

        and: "the replaced posterior and its boost survive only in the report"
        report.outcomes[0].action == 'rescored'
        report.outcomes[0].priorProbability == 0.97d
        report.outcomes[0].priorContributions == ['arbitrary_boost': 3.0d]
        !result.provenance[EC1].priorContributions.containsKey('arbitrary_boost')
    }

    def "rescoring depends only on raw claims, not on the previous posterior, likelihood or boosts"() {
        given:
        def claims = [claim('p1', 'EC:1.1.1.1', 0.8)]
        def other = new IntegratedAnnotationSet()
        def otherState = new IntegrationState(new Genome(id: 'g'))
        retain(EC1, 0.99d, claims)
        other.put(new ClaimProvenance(functionKey: EC1, proteinId: 'p1', likelihoodLogOdds: -5.0d, finalProbability: 0.02d),
            new Annotation(type: AnnotationType.EC, value: 'EC:1.1.1.1'))
        otherState.set(EC1, logit(0.02d))
        otherState.setPinnedFloor(new gspa.integration.ClaimKey('p1', AnnotationType.EC, 'EC:1.1.1.1'), 4.0d)

        when:
        stage(standardModel()).run(claims, state, result)
        stage(standardModel()).run(claims, otherState, other)

        then:
        result.annotationsByFunction[EC1].score == other.annotationsByFunction[EC1].score
    }

    // ------------------------------------------------------------------ rise / fall / retract

    def "an uncertain weak candidate is recovered by favourable context and decreased by unfavourable context"() {
        given:
        def claims = [claim('p1', 'EC:1.1.1.1', 0.30d)]
        def model = standardModel()

        when: "the same retained candidate under three contexts"
        def score = { ExternalContext c ->
            def r = new IntegratedAnnotationSet()
            def st = new IntegrationState(new Genome(id: 'g'))
            r.put(new ClaimProvenance(functionKey: EC1, proteinId: 'p1', finalProbability: 0.40d), new Annotation(type: AnnotationType.EC, value: 'EC:1.1.1.1'))
            stage(model).run(claims, st, r, c)
            r.annotationsByFunction[EC1].score
        }
        double none = score(null)
        double favourable = score(nbr(EC1, 0.95d, 0.0d))
        double unfavourable = score(nbr(EC1, 0.0d, 0.95d))

        then: "the weak candidate is neither dropped nor pinned: it rises or falls with the context"
        none > 0 && none < 0.5d
        favourable > none
        unfavourable < none
    }

    def "retractBelow removes the annotation consistently; without it a falling score never removes"() {
        given:
        def claims = [claim('p1', 'EC:1.1.1.1', 0.8)]
        def severe = nbr(EC1, 0.0d, 5.0d)
        retain(EC1, 0.9d, claims)

        when: "no retraction threshold"
        def kept = stage(standardModel()).run(claims, state, result, severe)

        then:
        kept.outcomes[0].action == 'rescored'
        result.annotationsByFunction.containsKey(EC1)
        result.annotationsByFunction[EC1].score < 1e-3

        when: "with a threshold"
        retain(EC1, 0.9d, claims)
        def gone = stage(standardModel(), [retractBelow: 0.05d]).run(claims, state, result, severe)

        then:
        gone.outcomes[0].action == 'retracted'
        gone.outcomes[0].learnedProbability < 0.05d
        gone.outcomes[0].priorProbability == 0.9d
        !result.provenance.containsKey(EC1)
        !result.annotationsByFunction.containsKey(EC1)
        result.annotations.annotations.isEmpty()
        !state.posteriorLogOdds.containsKey(EC1)
        result.processedProteinIds.contains('p1')
    }

    def "pinned floors of decided keys are dropped so no later stage floors the learned probability"() {
        given:
        def claims = [claim('p1', 'EC:1.1.1.1', 0.8)]
        retain(EC1, 0.9d, claims)
        state.setPinnedFloor(new gspa.integration.ClaimKey('p1', AnnotationType.EC, 'EC:1.1.1.1'), 3.0d)
        state.pinnedFloors['unrelated|EC|EC:9.9.9.9'] = 1.0d

        when:
        def report = stage(standardModel()).run(claims, state, result, nbr(EC1, 0.0d, 5.0d))

        then:
        report.pinnedFloorsCleared == [EC1]
        !state.pinnedFloors.containsKey(EC1)
        state.pinnedFloors.containsKey('unrelated|EC|EC:9.9.9.9')
        state.getProbability(EC1) < 1e-3
    }

    // ------------------------------------------------------------------ roles

    def "an accessory competitor never raises the catalytic score; role scores are all retained"() {
        given:
        def model = logisticModel(['seq.group.homology.max_raw': 4.0d, 'role.accessory': 3.0d, 'role.catalytic': -1.0d], -2.0d)
        def claims = [claim('p1', 'EC:1.1.1.1', 0.8)]
        retain(EC1, 0.5d, claims)
        def ctx = new ExternalContext(roles: [(EC1): [CandidateRole.ACCESSORY, CandidateRole.CATALYTIC]])

        when:
        def report = stage(model).run(claims, state, result, ctx)
        def rec = report.candidates.collectEntries { [it.role.name(), it] }

        then: "two candidates, scored separately"
        rec.keySet() == ['ACCESSORY', 'CATALYTIC'] as Set
        rec['ACCESSORY'].probability > rec['CATALYTIC'].probability

        and: "the EC score is the catalytic candidate's, not the maximum or a union over roles"
        result.annotationsByFunction[EC1].score == rec['CATALYTIC'].probability
        rec['CATALYTIC'].projected && !rec['ACCESSORY'].projected
        report.outcomes[0].projectedRole == CandidateRole.CATALYTIC

        and: "role-specific scores stay auditable in report, provenance and annotation"
        report.outcomes[0].roleProbabilities['ACCESSORY'] == rec['ACCESSORY'].probability
        result.provenance[EC1].priorContributions['learned.probability.ACCESSORY'] == rec['ACCESSORY'].probability
        result.provenance[EC1].priorContributions['learned.probability.CATALYTIC'] == rec['CATALYTIC'].probability
        result.annotationsByFunction[EC1].metadata.learned_role == 'CATALYTIC'
        result.annotationsByFunction[EC1].metadata.learned_role_scores.contains('ACCESSORY=')
    }

    def "an EC with only an accessory candidate gets no legacy score at all"() {
        given:
        def claims = [claim('p1', 'EC:1.1.1.1', 0.99d)]
        retain(EC1, 0.8d, claims)
        def ctx = new ExternalContext(roles: [(EC1): [CandidateRole.ACCESSORY]])

        when:
        def report = stage(standardModel()).run(claims, state, result, ctx)

        then:
        report.outcomes[0].action == 'removed_no_projection'
        report.candidates.size() == 1 && report.candidates[0].role == CandidateRole.ACCESSORY
        !result.annotationsByFunction.containsKey(EC1)
    }

    def "a GO complex term carries the complex-subunit probability, kept apart from the catalytic EC of the same protein"() {
        given:
        state.complexClassification['GO:0005840'] = 'n'
        def claims = [
            claim('p1', 'GO:0005840', 0.9d, 'interproscan', EvidenceType.SEQUENCE_DOMAIN, AnnotationType.GO, [:], 'CC'),
            claim('p1', 'EC:2.7.7.6', 0.4d)]
        retain('p1|GO|GO:0005840', 0.5d, claims)
        retain('p1|EC|EC:2.7.7.6', 0.5d, claims)
        def model = logisticModel(['seq.group.homology.max_raw': 4.0d, 'role.complex_subunit': 1.0d, 'role.catalytic': -1.0d], -2.0d)

        when:
        def report = stage(model).run(claims, state, result)
        def byKey = report.outcomes.collectEntries { [it.functionKey, it] }

        then:
        byKey['p1|GO|GO:0005840'].projectedRole == CandidateRole.COMPLEX_SUBUNIT
        byKey['p1|EC|EC:2.7.7.6'].projectedRole == CandidateRole.CATALYTIC
        byKey['p1|GO|GO:0005840'].learnedProbability != byKey['p1|EC|EC:2.7.7.6'].learnedProbability
    }

    // ------------------------------------------------------------------ candidates

    def "weak candidates named only by structured context are added and scored with no sequence evidence"() {
        given:
        def claims = [claim('p1', 'EC:1.1.1.1', 0.8)]
        retain(EC1, 0.8d, claims)
        def newKey = 'p9|EC|EC:3.3.3.3'
        def ctx = new ExternalContext(neighborhood: [(newKey): new ExternalContext.Neighborhood(supportMax: 0.2d, nSupporting: 1, conflictMax: 0.0d)])

        when:
        def report = stage(standardModel()).run(claims, state, result, ctx)

        then:
        report.outcomes.find { it.functionKey == newKey }.action == 'added'
        def added = result.annotationsByFunction[newKey]
        added.score > 0 && added.score < 0.5d
        report.candidates.find { it.functionKey == newKey }.supportingClaims.isEmpty()
        report.candidates.find { it.functionKey == newKey }.featureNotes != null
        state.posteriorLogOdds.containsKey(newKey)
        result.annotations.size() == 2

        when: "context candidates switched off"
        def r2 = new IntegratedAnnotationSet()
        def st2 = new IntegrationState(new Genome(id: 'g'))
        def rep2 = stage(standardModel(), [addContextCandidates: false]).run(claims, st2, r2, ctx)

        then:
        r2.annotations.size() == 0
        rep2.candidates*.functionKey == [EC1]        // only the claimed function; the context-only key is not a candidate
    }

    def "retained annotations without any raw claim are still rescored, never left at their old posterior"() {
        given:
        retain('p5|EC|EC:5.5.5.5', 0.9d)

        when:
        def report = stage(standardModel()).run([], state, result)

        then:
        report.outcomes[0].action == 'rescored'
        result.annotationsByFunction['p5|EC|EC:5.5.5.5'].score < 0.1d       // no sequence evidence: the learned score is low
    }

    def "claimed functions the earlier stages dropped are reported but only added on request"() {
        given:
        def claims = [claim('p1', 'EC:1.1.1.1', 0.8)]

        when:
        def report = stage(standardModel()).run(claims, state, result)

        then:
        report.outcomes[0].action == 'scored_only'
        result.annotations.size() == 0
        report.candidates.size() == 1

        when:
        def report2 = stage(standardModel(), [addUnretainedClaims: true]).run(claims, state, result)

        then:
        report2.outcomes[0].action == 'added'
        result.annotationsByFunction[EC1].score == report2.candidates[0].probability
    }

    def "type scope leaves other keys untouched"() {
        given:
        def claims = [claim('p1', 'EC:1.1.1.1', 0.8),
                      claim('p1', 'GO:0016301', 0.8, 'diamond', EvidenceType.SEQUENCE_SIMILARITY, AnnotationType.GO)]
        retain(EC1, 0.77d, claims)
        retain('p1|GO|GO:0016301', 0.66d, claims)

        when:
        stage(standardModel(), [types: [AnnotationType.EC] as Set]).run(claims, state, result)

        then:
        result.annotationsByFunction['p1|GO|GO:0016301'].score == 0.66d
        result.annotationsByFunction[EC1].score != 0.77d
    }

    // ------------------------------------------------------------------ identity coherence

    def "state, provenance and annotations share one identity; isozymes with the same EC stay distinct"() {
        given:
        def claims = [claim('p1', 'EC:1.1.1.1', 0.9d), claim('p2', 'EC:1.1.1.1', 0.2d), claim('p2', 'EC:2.2.2.2', 0.6d)]
        claims.each { retain(it.functionKey(), 0.5d, [it]) }

        when:
        stage(standardModel()).run(claims, state, result)

        then:
        result.provenance.keySet() == result.annotationsByFunction.keySet()
        result.provenance.keySet() == state.posteriorLogOdds.keySet()
        result.annotations.annotations.collect { it.metadata.function_key } as Set == result.provenance.keySet()
        result.annotations.size() == 3
        result.annotationsByFunction.every { k, a -> a.score == result.provenance[k].finalProbability && Math.abs(state.getProbability(k) - a.score) < 1e-9 }
        result.annotationsByFunction['p1|EC|EC:1.1.1.1'].score != result.annotationsByFunction['p2|EC|EC:1.1.1.1'].score
    }

    def "an annotation without provenance is an incoherent set and is refused before any change"() {
        given:
        result.annotationsByFunction[EC1] = new Annotation(type: AnnotationType.EC, value: 'EC:1.1.1.1', score: 0.5d)
        def before = snapshot()

        when:
        stage(standardModel()).run([claim('p1', 'EC:1.1.1.1', 0.8)], state, result)

        then:
        thrown(LearnedStageException)
        snapshot() == before
    }

    // ------------------------------------------------------------------ atomicity

    def "any invalid input fails the whole stage with nothing mutated"() {
        given:
        def claims = [claim('p1', 'EC:1.1.1.1', 0.8), claim('p2', 'EC:2.2.2.2', 0.8)]
        claims.each { retain(it.functionKey(), 0.6d, [it]) }
        state.setPinnedFloor(new gspa.integration.ClaimKey('p1', AnnotationType.EC, 'EC:1.1.1.1'), 2.0d)
        def before = snapshot()

        when:
        stage(standardModel()).run(claims, state, result, ctxOf)

        then:
        thrown(LearnedStageException)
        snapshot() == before

        where:
        ctxOf << [
            nbr('p2|EC|EC:2.2.2.2', Double.NaN, 0.0d),                                           // non-finite neighbourhood
            nbr('p2|EC|EC:2.2.2.2', Double.POSITIVE_INFINITY, 0.0d),
            new ExternalContext(assemblyCompleteness: 250.0d),                                     // impossible percentage
            new ExternalContext(neighborhood: ['not a function key': new ExternalContext.Neighborhood()]),
            new ExternalContext(predictorStatus: ['diamond': 'MAYBE']),
        ]
    }

    def "a candidate whose margin overflows fails the stage even if it is the last of many"() {
        given:
        def j = FeatureContract.CURRENT.indexOf('ctx.nbr.conflict_max')
        def model = logisticModel(['ctx.nbr.conflict_max': 1.0e308d], 0.0d)
        def claims = [claim('p1', 'EC:1.1.1.1', 0.8), claim('p2', 'EC:2.2.2.2', 0.8)]
        claims.each { retain(it.functionKey(), 0.6d, [it]) }
        def ctx = new ExternalContext(neighborhood: ['p2|EC|EC:2.2.2.2': new ExternalContext.Neighborhood(supportMax: 0d, nSupporting: 0, conflictMax: 1.0e10d)])
        def before = snapshot()

        when:
        stage(model).run(claims, state, result, ctx)

        then:
        j >= 0
        def e = thrown(LearnedStageException)
        e.message.contains('p2|EC|EC:2.2.2.2')
        snapshot() == before
    }

    def "a malformed claim fails before scoring"() {
        given:
        retain(EC1, 0.6d)
        def before = snapshot()

        when:
        stage(standardModel()).run([new EvidenceClaim(proteinId: null, functionType: AnnotationType.EC, functionId: 'EC:1.1.1.1')], state, result)

        then:
        thrown(LearnedStageException)
        snapshot() == before
    }

    def "a failure while applying restores state and annotations"() {
        given: "an out-of-scope key carries an impossible probability that only the synchronise step notices"
        def claims = [claim('p1', 'EC:1.1.1.1', 0.8)]
        retain(EC1, 0.6d, claims)
        result.provenance['p7|GO|GO:0000001'] = new ClaimProvenance(functionKey: 'p7|GO|GO:0000001', proteinId: 'p7', finalProbability: 1.5d)
        def before = snapshot()

        when:
        stage(standardModel(), [types: [AnnotationType.EC] as Set]).run(claims, state, result)

        then:
        def e = thrown(LearnedStageException)
        e.message.contains('restored')
        snapshot().findAll { k, v -> k != 'ann' } == before.findAll { k, v -> k != 'ann' }
        result.annotationsByFunction[EC1].score == 0.6d
    }

    // ------------------------------------------------------------------ marginal vs null

    def "the conditional locus ranking and explicit null are reported next to, never inside, the marginal score"() {
        given:
        def model = standardModel(locusLayer(1.0d, 0.0d))
        def claims = [claim('p1', 'EC:1.1.1.1', 0.9d), claim('p2', 'EC:1.1.1.1', 0.3d)]
        claims.each { retain(it.functionKey(), 0.5d, [it]) }

        when:
        def report = stage(model).run(claims, state, result)
        def rec = report.candidates.collectEntries { [it.proteinId, it] }

        then:
        rec['p1'].locusNullProbability > 0 && rec['p1'].locusNullProbability == rec['p2'].locusNullProbability
        Math.abs(rec.values().sum { it.locusProbability } + rec['p1'].locusNullProbability - 1.0d) < 1e-12
        result.annotationsByFunction[EC1].score == rec['p1'].probability
        result.annotationsByFunction[EC1].score != rec['p1'].locusProbability

        and: "a model without a locus layer simply reports none"
        def r2 = new IntegratedAnnotationSet()
        def rep2 = stage(standardModel(), [addUnretainedClaims: true]).run(claims, new IntegrationState(new Genome(id: 'g')), r2)
        rep2.candidates.every { it.locusProbability == null && it.locusNullProbability == null }
        rep2.locusUnavailableReason != null
    }

    // ------------------------------------------------------------------ audit

    def "model hash, contract, features digest, raw claims and prior vs learned score are auditable"() {
        given:
        def model = standardModel()
        def claims = [claim('p1', 'EC:1.1.1.1', 0.8, 'diamond', EvidenceType.SEQUENCE_SIMILARITY, AnnotationType.EC, [pident: 88.0d])]
        retain(EC1, 0.97d, claims)

        when:
        def report = stage(model, [retainFeatureVectors: true]).run(claims, state, result)
        def rec = report.candidates[0]
        def ann = result.annotationsByFunction[EC1]

        then:
        report.modelSha256 == model.sha256 && report.modelSha256.size() == 64
        report.contractFingerprint == FeatureContract.CURRENT.fingerprint()
        ann.metadata.learned_model_sha256 == model.sha256
        ann.metadata.learned_contract_fingerprint == FeatureContract.CURRENT.fingerprint()
        ann.metadata.learned_replaced_posterior == 0.97d
        rec.featureDigest == new ContextFeatureExtractor().extractAll(claims, state)[0].digest()
        rec.features.length == FeatureContract.CURRENT.size()
        rec.supportingClaims.is(claims[0]) || rec.supportingClaims == claims
        result.provenance[EC1].supportingClaims == claims
        report.outcomes[0].priorProbability == 0.97d
        report.outcomes[0].learnedProbability == ann.score
    }
}
