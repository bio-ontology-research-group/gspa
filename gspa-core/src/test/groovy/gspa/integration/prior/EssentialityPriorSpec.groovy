package gspa.integration.prior

import gspa.config.EssentialFunctions
import gspa.integration.IntegrationState
import gspa.model.Genome
import gspa.ontology.GoOntology
import gspa.ontology.GoReasoner
import spock.lang.Specification

class EssentialityPriorSpec extends Specification {

    def "no-ops when no essential function profile is wired"() {
        given:
        def state = new IntegrationState(new Genome(id: 'g'))
        def prior = new EssentialityPrior()

        when:
        prior.beginIteration(state)

        then:
        prior.logOddsBoost('p1', 'p1|GO|GO:0006412', state) == 0.0
    }

    def "no-ops when goReasoner is missing even if essential functions present"() {
        given:
        def state = new IntegrationState(new Genome(id: 'g'))
        state.essentialFunctions = Mock(EssentialFunctions) {
            getGoTerms() >> (['GO:0006412'] as Set)
        }
        def prior = new EssentialityPrior()

        when:
        prior.beginIteration(state)

        then:
        prior.logOddsBoost('p1', 'p1|GO|GO:0006412', state) == 0.0
    }

    def "boosts candidate claim whose function is a descendant of an uncovered essential"() {
        given: "an essential GO:0006412 translation, uncovered by current annotations"
        def goOntology = Mock(GoOntology)
        // No current annotations → nothing propagated.
        goOntology.propagateAnnotations(_) >> ([] as Set)

        def goReasoner = Mock(GoReasoner)
        // Descendants of GO:0006412 include GO:0001234 (a ribosomal subclass)
        goReasoner.getSubClasses('GO:0006412', false) >> (['GO:0001234', 'GO:0005678'] as Set)

        def ef = Mock(EssentialFunctions)
        ef.getGoTerms() >> (['GO:0006412'] as Set)

        def state = new IntegrationState(new Genome(id: 'g'))
        state.goOntology = goOntology
        state.goReasoner = goReasoner
        state.essentialFunctions = ef

        def prior = new EssentialityPrior()
        prior.alphaEss = 1.5d

        when:
        prior.beginIteration(state)

        then:
        // A descendant claim gets the boost
        prior.logOddsBoost('p1', 'p1|GO|GO:0001234', state) == 1.5d
        prior.logOddsBoost('p1', 'p1|GO|GO:0005678', state) == 1.5d
        // An unrelated claim does not
        prior.logOddsBoost('p1', 'p1|GO|GO:0009999', state) == 0.0d
        // The essential term itself is also considered a "descendant of itself"
        prior.logOddsBoost('p1', 'p1|GO|GO:0006412', state) == 1.5d
    }

    def "does not boost when the essential is already covered"() {
        given:
        def goOntology = Mock(GoOntology)
        // Simulate the essential being present in the propagated annotation set.
        goOntology.propagateAnnotations(_) >> (['GO:0006412'] as Set)
        def goReasoner = Mock(GoReasoner)
        def ef = Mock(EssentialFunctions)
        ef.getGoTerms() >> (['GO:0006412'] as Set)

        def state = new IntegrationState(new Genome(id: 'g'))
        state.goOntology = goOntology
        state.goReasoner = goReasoner
        state.essentialFunctions = ef

        def prior = new EssentialityPrior()

        when:
        prior.beginIteration(state)

        then:
        prior.logOddsBoost('p1', 'p1|GO|GO:0001234', state) == 0.0d
        0 * goReasoner.getSubClasses(_, _)   // short-circuited: nothing uncovered
    }

    def "non-GO claims are never boosted"() {
        given:
        def goOntology = Mock(GoOntology)
        goOntology.propagateAnnotations(_) >> ([] as Set)
        def goReasoner = Mock(GoReasoner)
        goReasoner.getSubClasses('GO:0006412', false) >> (['GO:0001234'] as Set)
        def ef = Mock(EssentialFunctions)
        ef.getGoTerms() >> (['GO:0006412'] as Set)
        def state = new IntegrationState(new Genome(id: 'g'))
        state.goOntology = goOntology
        state.goReasoner = goReasoner
        state.essentialFunctions = ef

        def prior = new EssentialityPrior()

        when:
        prior.beginIteration(state)

        then:
        prior.logOddsBoost('p1', 'p1|EC|EC:2.7.1.1', state) == 0.0d
        prior.logOddsBoost('p1', 'p1|PFAM|PF00001', state) == 0.0d
    }
}
