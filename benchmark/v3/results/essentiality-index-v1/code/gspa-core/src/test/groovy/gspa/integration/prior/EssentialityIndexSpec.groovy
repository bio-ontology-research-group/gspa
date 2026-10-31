package gspa.integration.prior

import gspa.config.EssentialFunctions
import gspa.integration.IntegrationState
import gspa.model.Genome
import gspa.ontology.GoOntology
import spock.lang.Specification

class EssentialityIndexSpec extends Specification {
    private static class CountedKeys extends LinkedHashMap<String, Double> {
        private long keyVisits = 0
        @Override Set<String> keySet() {
            this.@keyVisits += size()
            super.keySet()
        }
        long visits() { this.@keyVisits }
    }

    def 'irrelevant EC candidates are not rescanned for every essential GO category'() {
        given:
        def state = new IntegrationState(new Genome(id: 'synthetic'))
        state.essentialFunctions = new EssentialFunctions(functions:
            (0..<32).collectEntries { ["GO:E${it}".toString(), 'test'] })
        def values = new CountedKeys()
        (0..<2000).each { values.put("p${it}|EC|1.1.1.1".toString(), -0.5d) }
        state.posteriorLogOdds = values
        state.claimCounts = ['p0|EC|1.1.1.1': 1]
        def prior = new EssentialityPrior(rescueMode: mode)

        when:
        prior.beginIteration(state)
        long visited = values.visits()
        int candidateCount = values.size()

        then:
        visited <= candidateCount
        state.rescueLog.isEmpty()
        state.pinnedFloors.isEmpty()
        prior.logOddsBoost('p0', 'p0|EC|1.1.1.1', state) == 0.0d

        where:
        mode << ['legacy', 'targeted']
    }

    def 'overlapping closures preserve mode-specific tie order and rebuild after same-size replacement'() {
        given:
        def state = new IntegrationState(new Genome(id: 'synthetic'))
        state.essentialFunctions = new EssentialFunctions(functions: ['GO:E1': 'test', 'GO:E2': 'test'])
        state.goOntology = Mock(GoOntology) {
            getDescendants('GO:E1') >> (['GO:A', 'GO:Z'] as Set)
            getDescendants('GO:E2') >> (['GO:Z', 'GO:A'] as Set)
            getAncestors(_) >> ([] as Set)
            getAspect(_) >> 'BP'
            propagateAnnotations(_) >> { Set<String> terms -> terms }
        }
        double l0 = Math.log(0.45d / 0.55d)
        state.posteriorLogOdds = ['z|GO|GO:Z': l0, 'ec|EC|1.1.1.1': l0, 'a|GO|GO:A': l0]
        state.likelihoodLogOdds = state.posteriorLogOdds
        state.claimCounts = state.posteriorLogOdds.collectEntries { k, v -> [k, 1] }
        def prior = new EssentialityPrior(rescueMode: mode)

        when:
        prior.beginIteration(state)

        then:
        state.rescueLog*.category == ['GO:E1', 'GO:E2']
        state.rescueLog*.functionKey == [first, first]
        prior.logOddsBoost('', first, state) == 1.30d * 1.5d
        prior.logOddsBoost('ec', 'ec|EC|1.1.1.1', state) == 0.0d
        state.pinnedFloors.isEmpty() == (mode == 'legacy')

        when:
        state.posteriorLogOdds.remove('z|GO|GO:Z')
        state.posteriorLogOdds['0|GO|GO:Z'] = l0
        state.likelihoodLogOdds = state.posteriorLogOdds
        prior.beginIteration(state)

        then:
        state.rescueLog*.functionKey == [first, first, second, second]
        prior.logOddsBoost('', second, state) == 1.30d * 1.5d
        prior.logOddsBoost('z', 'z|GO|GO:Z', state) == 0.0d

        where:
        mode       | first        | second
        'legacy'   | 'z|GO|GO:Z'  | 'a|GO|GO:A'
        'targeted' | 'a|GO|GO:A'  | '0|GO|GO:Z'
    }
}
