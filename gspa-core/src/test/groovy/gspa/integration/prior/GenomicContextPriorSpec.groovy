package gspa.integration.prior

import gspa.integration.IntegrationState
import gspa.integration.MetabolicGap
import gspa.model.Genome
import gspa.ontology.PathwayDatabase
import gspa.ontology.PathwayGraph
import spock.lang.Specification

class GenomicContextPriorSpec extends Specification {

    private static PathwayDatabase buildPathwayDb() {
        def db = new PathwayDatabase()
        def pwy1 = new PathwayGraph(pathwayId: 'PWY-1', pathwayName: 'Ribosome biogenesis')
        // Mock getRequiredGoTerms() to return a fixed set.
        // Since PathwayGraph is a real class, we use metaprogramming to stub.
        pwy1.metaClass.getRequiredGoTerms = { -> ['GO:0006412', 'GO:0003735', 'GO:0019843'] as Set }
        db.pathways['PWY-1'] = pwy1

        def pwy2 = new PathwayGraph(pathwayId: 'PWY-2', pathwayName: 'Glycolysis')
        pwy2.metaClass.getRequiredGoTerms = { -> ['GO:0006096', 'GO:0004340'] as Set }
        db.pathways['PWY-2'] = pwy2

        db
    }

    def "no-ops with no operons"() {
        given:
        def state = new IntegrationState(new Genome(id: 'g'))
        state.pathwayDatabase = buildPathwayDb()
        def prior = new GenomicContextPrior()

        when:
        prior.beginIteration(state)

        then:
        prior.logOddsBoost('p1', 'p1|GO|GO:0006412', state) == 0.0d
    }

    def "boosts a weak candidate claim whose pathway is the operon's consensus"() {
        given: "an operon of three genes: p1 and p2 are strongly annotated for ribosome GO terms, p3 has a weak claim for another ribosome term"
        def state = new IntegrationState(new Genome(id: 'g'))
        state.pathwayDatabase = buildPathwayDb()
        state.operons = [['p1', 'p2', 'p3']]

        // p1 strong for GO:0006412 (translation)
        state.set('p1|GO|GO:0006412', 5.0d)
        // p2 strong for GO:0003735 (structural ribosome)
        state.set('p2|GO|GO:0003735', 5.0d)
        // p3 weak for GO:0019843 (rRNA binding) — same pathway, should be boosted
        state.set('p3|GO|GO:0019843', -2.0d)
        // p3 weak for GO:0004340 (glycolysis) — different pathway, should NOT be boosted
        state.set('p3|GO|GO:0004340', -2.0d)

        def prior = new GenomicContextPrior()
        prior.alphaCtx = 1.0d

        when:
        prior.beginIteration(state)

        then: "the weak ribosome claim is boosted by alphaCtx * consensusStrength = 1.0 * 2/3"
        double boost = prior.logOddsBoost('p3', 'p3|GO|GO:0019843', state)
        Math.abs(boost - (2.0d / 3.0d)) < 1e-6

        and: "the weak glycolysis claim is not boosted"
        prior.logOddsBoost('p3', 'p3|GO|GO:0004340', state) == 0.0d
    }

    def "gap-filling multiplier applies when the boosted function matches a metabolic gap"() {
        given:
        def state = new IntegrationState(new Genome(id: 'g'))
        state.pathwayDatabase = buildPathwayDb()
        state.operons = [['p1', 'p2']]
        state.set('p1|GO|GO:0006412', 5.0d)   // strong ribosome annotation
        state.set('p2|GO|GO:0003735', -1.0d)  // weak candidate to fill a gap
        state.metabolicGaps = [new MetabolicGap(goTerm: 'GO:0003735')]

        def prior = new GenomicContextPrior()
        prior.alphaCtx = 1.0d
        prior.alphaGapCtx = 2.0d

        when:
        prior.beginIteration(state)

        then: "consensus = 1/2 = 0.5; boost = 1.0 * 0.5 * 2.0 = 1.0"
        double boost = prior.logOddsBoost('p2', 'p2|GO|GO:0003735', state)
        Math.abs(boost - 1.0d) < 1e-6
    }

    def "operons smaller than 2 are ignored"() {
        given:
        def state = new IntegrationState(new Genome(id: 'g'))
        state.pathwayDatabase = buildPathwayDb()
        state.operons = [['p1']]
        state.set('p1|GO|GO:0006412', 5.0d)
        def prior = new GenomicContextPrior()

        when:
        prior.beginIteration(state)

        then:
        prior.logOddsBoost('p1', 'p1|GO|GO:0006412', state) == 0.0d
    }
}
