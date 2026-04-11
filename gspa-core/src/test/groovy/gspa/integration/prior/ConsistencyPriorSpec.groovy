package gspa.integration.prior

import gspa.integration.IntegrationState
import gspa.model.ConsistencyViolation
import gspa.model.Genome
import gspa.ontology.ConsistencyResult
import gspa.ontology.SatConsistencyChecker
import spock.lang.Specification

class ConsistencyPriorSpec extends Specification {

    def "no-ops when the SAT checker is not wired"() {
        given:
        def state = new IntegrationState(new Genome(id: 'g'))
        def prior = new ConsistencyPrior()

        when:
        prior.beginIteration(state)

        then:
        prior.logOddsBoost('p1', 'p1|GO|GO:0015979', state) == 0.0d
    }

    def "no-ops when the SAT check is consistent"() {
        given:
        def checker = Mock(SatConsistencyChecker)
        checker.check(_) >> new ConsistencyResult(consistent: true, violations: [])
        def state = new IntegrationState(new Genome(id: 'g'))
        state.satConsistencyChecker = checker
        // Seed state with a strong posterior so currentlyAnnotatedGoTerms is non-empty.
        state.set('p1|GO|GO:0015979', 5.0d)
        def prior = new ConsistencyPrior()

        when:
        prior.beginIteration(state)

        then:
        prior.logOddsBoost('p1', 'p1|GO|GO:0015979', state) == 0.0d
    }

    def "soft-downweights claims whose GO term is in the UNSAT core"() {
        given:
        def checker = Mock(SatConsistencyChecker)
        def violation = new ConsistencyViolation(
            involvedGoTerms: ['GO:0015979', 'GO:0009523']
        )
        checker.check(_) >> new ConsistencyResult(
            consistent: false,
            violations: [violation],
        )
        def state = new IntegrationState(new Genome(id: 'g'))
        state.satConsistencyChecker = checker
        state.set('p1|GO|GO:0015979', 5.0d)    // strong, will enter MAP-annotated set
        def prior = new ConsistencyPrior()
        prior.alphaCons = 3.0d

        when:
        prior.beginIteration(state)

        then:
        prior.logOddsBoost('p1', 'p1|GO|GO:0015979', state) == -3.0d
        prior.logOddsBoost('p1', 'p1|GO|GO:0009523', state) == -3.0d
        prior.logOddsBoost('p1', 'p1|GO|GO:0006412', state) == 0.0d
    }

    def "hard-filter mode returns a massive negative penalty"() {
        given:
        def checker = Mock(SatConsistencyChecker)
        checker.check(_) >> new ConsistencyResult(
            consistent: false,
            violations: [new ConsistencyViolation(involvedGoTerms: ['GO:0015979'])],
        )
        def state = new IntegrationState(new Genome(id: 'g'))
        state.satConsistencyChecker = checker
        state.set('p1|GO|GO:0015979', 5.0d)
        def prior = new ConsistencyPrior()
        prior.hardFilter = true

        when:
        prior.beginIteration(state)

        then:
        prior.logOddsBoost('p1', 'p1|GO|GO:0015979', state) == -1000.0d
    }

    def "non-GO claims are never penalized"() {
        given:
        def checker = Mock(SatConsistencyChecker)
        checker.check(_) >> new ConsistencyResult(
            consistent: false,
            violations: [new ConsistencyViolation(involvedGoTerms: ['GO:0015979'])],
        )
        def state = new IntegrationState(new Genome(id: 'g'))
        state.satConsistencyChecker = checker
        state.set('p1|GO|GO:0015979', 5.0d)
        def prior = new ConsistencyPrior()

        when:
        prior.beginIteration(state)

        then:
        prior.logOddsBoost('p1', 'p1|EC|EC:1.1.1.1', state) == 0.0d
    }
}
