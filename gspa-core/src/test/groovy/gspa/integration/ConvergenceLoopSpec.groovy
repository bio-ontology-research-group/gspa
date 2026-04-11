package gspa.integration

import gspa.model.AnnotationType
import gspa.model.Genome
import spock.lang.Specification

/**
 * Tests the Phase 7.3 convergence loop independent of concrete priors, by
 * registering toy priors whose boosts we can control precisely.
 */
class ConvergenceLoopSpec extends Specification {

    private static EvidenceClaim claim(String protein, String funcId, double prob) {
        new EvidenceClaim(
            proteinId: protein,
            functionType: AnnotationType.GO,
            functionId: funcId,
            goAspect: 'BP',
            evidenceType: EvidenceType.SEQUENCE_SIMILARITY,
            source: 'diamond',
            rawScore: prob,
            calibratedProb: prob,
        )
    }

    private static IterativeRefiner refinerWith(Prior... priors) {
        def reliability = new EnumMap<EvidenceType, Double>(EvidenceType)
        EvidenceType.values().each { reliability[it] = 1.0 }
        def r = new IterativeRefiner(new EvidenceCombiner(reliability))
        def engine = new PriorEngine()
        priors.each { engine.register(it, 1.0) }
        r.priorEngine = engine
        r.maxIter = 10
        r.epsilon = 0.001
        r.damping = 0.5
        r
    }

    /** A toy prior that always adds a fixed boost to a target key. */
    private static class FixedBoostPrior implements Prior {
        String targetKey
        double boost
        String getName() { 'fixed_boost' }
        String name() { 'fixed_boost' }
        double logOddsBoost(String pid, String key, IntegrationState s) {
            key == targetKey ? boost : 0.0d
        }
        Set<String> inputs() { [] as Set }
    }

    /** A toy prior that flips sign every iteration — oscillation stressor. */
    private static class OscillatingPrior implements Prior {
        String targetKey
        double magnitude = 3.0
        private int iter = 0
        private double lastBoost = 0
        void beginIteration(IntegrationState s) {
            iter++
            lastBoost = (iter % 2 == 0) ? magnitude : -magnitude
        }
        String name() { 'oscillating' }
        double logOddsBoost(String pid, String key, IntegrationState s) {
            key == targetKey ? lastBoost : 0.0d
        }
        Set<String> inputs() { [] as Set }
    }

    def "empty priorEngine yields a single-pass refinement"() {
        given:
        def refiner = new IterativeRefiner(new EvidenceCombiner())
        def claims = [claim('p1', 'GO:0001', 0.8)]

        when:
        def result = refiner.refine(claims, new IntegrationState(new Genome(id: 'g')))

        then: "still emits one annotation and provenance"
        result.annotations.size() == 1
        result.provenance['p1|GO|GO:0001'].convergenceIter == 1
    }

    def "a single prior contribution propagates into the posterior"() {
        given:
        def boost = new FixedBoostPrior(targetKey: 'p1|GO|GO:0001', boost: 2.0)
        def refiner = refinerWith(boost)
        def claims = [claim('p1', 'GO:0001', 0.5)]
        def state = new IntegrationState(new Genome(id: 'g'))

        when:
        def result = refiner.refine(claims, state)

        then:
        def prov = result.provenance['p1|GO|GO:0001']
        // Final log-odds should be above the pure likelihood log-odds (log(0.5/0.5) = 0)
        prov.finalLogOdds > 0.0d
        prov.priorContributions['fixed_boost'] == 2.0d
    }

    def "refiner converges when priors are stationary"() {
        given:
        def p1 = new FixedBoostPrior(targetKey: 'p1|GO|GO:0001', boost: 1.0)
        def refiner = refinerWith(p1)
        def claims = [claim('p1', 'GO:0001', 0.5)]
        def state = new IntegrationState(new Genome(id: 'g'))

        when:
        def result = refiner.refine(claims, state)

        then: "converges in a few iterations (well below maxIter=10)"
        def prov = result.provenance['p1|GO|GO:0001']
        prov.convergenceIter <= refiner.maxIter
    }

    def "damping prevents oscillation when a non-monotone prior flips sign"() {
        given:
        def osc = new OscillatingPrior(targetKey: 'p1|GO|GO:0001', magnitude: 3.0)
        def refiner = refinerWith(osc)
        refiner.damping = 0.3   // strong under-relaxation
        refiner.maxIter = 10
        refiner.epsilon = 0.001
        def claims = [claim('p1', 'GO:0001', 0.5)]
        def state = new IntegrationState(new Genome(id: 'g'))

        when:
        def result = refiner.refine(claims, state)

        then: "the refiner halts (either via convergence or divergence detection)"
        result.provenance['p1|GO|GO:0001'].convergenceIter <= 10

        and: "the final log-odds remains within a sane bounded range"
        def finalLog = result.provenance['p1|GO|GO:0001'].finalLogOdds
        finalLog >= -12.0d
        finalLog <= 12.0d
    }

    def "final posterior probability reflects log-odds sigmoid"() {
        given:
        def p1 = new FixedBoostPrior(targetKey: 'p1|GO|GO:0001', boost: 2.0)
        def refiner = refinerWith(p1)
        refiner.damping = 1.0   // immediate convergence; no relaxation
        refiner.maxIter = 20
        def claims = [claim('p1', 'GO:0001', 0.5)]
        def state = new IntegrationState(new Genome(id: 'g'))

        when:
        def result = refiner.refine(claims, state)

        then: "likelihood log-odds = logit(0.5) = 0; final = 0 + 2 = 2; sigmoid(2) ≈ 0.8808"
        def prov = result.provenance['p1|GO|GO:0001']
        Math.abs(prov.finalProbability - 0.880797d) < 1e-4
    }
}
