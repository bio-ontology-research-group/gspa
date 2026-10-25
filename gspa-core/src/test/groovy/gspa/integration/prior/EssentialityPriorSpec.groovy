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

    def "without goReasoner, boosts exact essential matches but cannot expand to descendants"() {
        given: "essential GO:0006412 but no GoReasoner to expand to descendants"
        def state = new IntegrationState(new Genome(id: 'g'))
        state.essentialFunctions = Mock(EssentialFunctions) {
            getGoTerms() >> (['GO:0006412'] as Set)
        }
        // goReasoner and goOntology intentionally null
        def prior = new EssentialityPrior()
        prior.alphaEss = 1.5d

        when:
        prior.beginIteration(state)

        then: "exact match to the essential itself is boosted"
        prior.logOddsBoost('p1', 'p1|GO|GO:0006412', state) == 1.5d

        and: "a non-essential claim is not boosted (no descendants expanded)"
        prior.logOddsBoost('p1', 'p1|GO|GO:0001234', state) == 0.0d
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

class EssentialityPriorTargetedSpec extends Specification {

    /** Two uncovered-category candidates: p1 at 0.45 but confidently annotated elsewhere (same aspect),
     *  p2 at 0.30, otherwise unannotated, in an operon with an essential-domain gene. */
    private IntegrationState stateWithCandidates() {
        def go = Mock(GoOntology)
        go.getDescendants('GO:E') >> (['GO:D1'] as Set)
        go.getDescendants('GO:E2') >> ([] as Set)
        go.getAncestors(_) >> ([] as Set)
        go.getAspect(_) >> 'BP'
        go.propagateAnnotations(_) >> { Set<String> t -> new LinkedHashSet<String>(t) }
        def ef = Mock(EssentialFunctions)
        ef.getGoTerms() >> (['GO:E', 'GO:E2'] as Set)
        def state = new IntegrationState(new Genome(id: 'g'))
        state.goOntology = go
        state.essentialFunctions = ef
        Closure<Double> l = { double p -> Math.log(p / (1 - p)) }
        state.likelihoodLogOdds = [
            'p1|GO|GO:D1': l(0.45), 'p1|GO|GO:X': l(0.90),   // p1: candidate, but already a confident BP annotation outside the category
            'p2|GO|GO:E' : l(0.30),                          // p2: candidate, otherwise unannotated
            'p3|GO|GO:E2': l(0.80),                          // p3: covers the other essential category (operon partner of p2)
        ] as Map<String, Double>
        state.posteriorLogOdds = new LinkedHashMap<>(state.likelihoodLogOdds)
        state.claimCounts = state.likelihoodLogOdds.collectEntries { k, v -> [k, 1] } as Map<String, Integer>
        state.operons = [['p2', 'p3']]
        state
    }

    def "targeted mode rescues the free, operon-supported candidate and pins it to the target"() {
        given:
        def state = stateWithCandidates()
        def prior = new EssentialityPrior(rescueMode: 'targeted', rescueTarget: 0.60d)

        when:
        prior.beginIteration(state)

        then: "p2 is the rescue for GO:E, not the higher-likelihood p1"
        state.rescueLog*.functionKey == ['p2|GO|GO:E']
        state.rescueLog[0].category == 'GO:E'
        state.rescueLog[0].score > 0.0d

        and: "its pinned floor is the target posterior"
        Math.abs(state.getPinnedFloor(gspa.integration.ClaimKey.parse('p2|GO|GO:E')) - Math.log(0.6d / 0.4d)) < 1e-9
        state.getPinnedFloor(gspa.integration.ClaimKey.parse('p1|GO|GO:D1')) == null

        and: "the boost reaches at least the legacy magnitude and p1 gets only the domain boost"
        prior.logOddsBoost('p2', 'p2|GO|GO:E', state) >= 1.30d * 1.5d
        prior.logOddsBoost('p1', 'p1|GO|GO:D1', state) < 1.0d
    }

    def "legacy mode keeps the v2.1.0 tier choice and sets no floor"() {
        given:
        def state = stateWithCandidates()
        def prior = new EssentialityPrior(rescueMode: 'legacy')

        when:
        prior.beginIteration(state)

        then:
        state.rescueLog*.functionKey == ['p1|GO|GO:D1']
        state.pinnedFloors.isEmpty()
        prior.logOddsBoost('p1', 'p1|GO|GO:D1', state) == 1.30d * 1.5d
    }

    def "a rescued claim is logged once across iterations"() {
        given:
        def state = stateWithCandidates()
        def prior = new EssentialityPrior()

        when:
        prior.beginIteration(state)
        prior.beginIteration(state)

        then:
        state.rescueLog.size() == 1
    }
}
