package gspa.integration.prior

import gspa.integration.ClaimKey
import gspa.integration.IntegrationState
import gspa.model.Genome
import gspa.ontology.GoOntology
import gspa.ontology.ModuleAssessorSpec
import gspa.ontology.ModuleDefinitions
import spock.lang.Specification

class CoherencePriorModuleSpec extends Specification {

    private IntegrationState stateWith(Map<String, Double> like, List<List<String>> operons = []) {
        def go = Mock(GoOntology)
        go.getDescendants(_) >> ([] as Set)
        go.getAncestors(_) >> ([] as Set)
        go.getAspect(_) >> 'MF'
        go.propagateAnnotations(_) >> { Set<String> t -> new LinkedHashSet<String>(t) }
        def state = new IntegrationState(new Genome(id: 'g'))
        state.goOntology = go
        state.moduleDefinitions = ModuleDefinitions.load(ModuleAssessorSpec.defsFile(), ModuleAssessorSpec.ec2go)
        Closure<Double> l = { double p -> Math.log(p / (1 - p)) }
        state.likelihoodLogOdds = like.collectEntries { k, v -> [k, l(v)] } as Map<String, Double>
        state.posteriorLogOdds = new LinkedHashMap<>(state.likelihoodLogOdds)
        state.claimCounts = like.collectEntries { k, v -> [k, 1] } as Map<String, Integer>
        state.operons = operons
        state
    }

    def "an activated, half-complete module gets its evidenced step rescued and pinned; the operon partner wins"() {
        given: "M00001: enzyme A present, complex present, kinase step only below threshold on two proteins"
        def state = stateWith([
            'pA|GO|GO:0000001': 0.9, 'pC|GO|GO:0000004': 0.8,
            'pK1|GO|GO:0000003': 0.40, 'pK2|GO|GO:0000003': 0.35,
        ], [['pA', 'pK2']])
        def prior = new CoherencePrior()

        when:
        prior.beginIteration(state)

        then: "pK2 (lower likelihood, but in an operon with a present module enzyme) is the rescue"
        state.rescueLog*.functionKey == ['pK2|GO|GO:0000003']
        state.rescueLog[0].prior == 'coherence'
        state.rescueLog[0].evidence.module == 'M00001'
        state.rescueLog[0].evidence.operon_partner == 'pA'
        Math.abs(state.getPinnedFloor(ClaimKey.parse('pK2|GO|GO:0000003')) - Math.log(0.6d / 0.4d)) < 1e-9
        prior.logOddsBoost('pK2', 'pK2|GO|GO:0000003', state) >= Math.log(0.6d / 0.4d) - Math.log(0.35d / 0.65d) - 1e-9
        state.getPinnedFloor(ClaimKey.parse('pK1|GO|GO:0000003')) == null
    }

    def "a module activated by a shared enzyme alone, or below half complete, is not refined"() {
        given: "only the shared kinase present (no activation) and a kinase candidate elsewhere"
        def state = stateWith(['pK|GO|GO:0000003': 0.9, 'pB|GO|GO:0000008': 0.45])
        def prior = new CoherencePrior()

        when:
        prior.beginIteration(state)

        then:
        state.rescueLog.isEmpty()
        state.pinnedFloors.isEmpty()
    }

    def "steps without any candidate stay open: no claim is invented"() {
        given: "M00001 with enzyme A and the complex present; nothing at all for the kinase"
        def state = stateWith(['pA|GO|GO:0000001': 0.9, 'pC|GO|GO:0000004': 0.8])
        def prior = new CoherencePrior()

        when:
        prior.beginIteration(state)

        then:
        state.rescueLog.isEmpty()
    }

    def "legacy pathway mode ignores the module definitions"() {
        given:
        def state = stateWith(['pA|GO|GO:0000001': 0.9, 'pC|GO|GO:0000004': 0.8, 'pK|GO|GO:0000003': 0.40])
        def prior = new CoherencePrior(pathwayMode: 'legacy')

        when:
        prior.beginIteration(state)

        then:
        state.rescueLog.isEmpty()
    }
}
