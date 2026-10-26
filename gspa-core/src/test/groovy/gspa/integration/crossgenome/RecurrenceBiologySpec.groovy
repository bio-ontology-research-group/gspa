package gspa.integration.crossgenome

import gspa.integration.*
import gspa.integration.suggester.*
import gspa.model.AnnotationType
import spock.lang.Specification

class RecurrenceBiologySpec extends Specification {
    private static ReactionLocusCatalog negativeCatalog() {
        def catalog = new ReactionLocusCatalog(panelSize: 20)
        catalog.put(new ReactionLocusCatalog.Entry(orthogroupId: 'common_bystander', reactionId: 'missing_reaction',
            nSigWith: 1, nSigTotal: 20, nBaseWith: 20, nBaseTotal: 20))
        catalog
    }

    def "a common bystander rarely near the reaction loses support and may remain unassigned"() {
        given:
        def state = new IntegrationState()
        state.orthogroupMap = [bystander: 'common_bystander']
        def result = new IntegratedAnnotationSet(suggestions: [new SingletonSuggestion(
            source: 'reaction_local_context', proteinId: 'bystander', q: 0.6d,
            reactionId: 'missing_reaction', functionId: '1.1.1.1', functionType: AnnotationType.EC)])

        when:
        new CrossGenomeReScorer(requireCredible: false).rescore(state, result, negativeCatalog())

        then:
        def suggestion = result.suggestions.first()
        suggestion instanceof DisjunctiveSuggestion
        suggestion.candidateProbabilities.bystander < 0.6d
        suggestion.unassignedProbability > 0.5d
        Math.abs(suggestion.candidateProbabilities.values().sum() + suggestion.unassignedProbability - 1d) < 1e-12
    }

    def "recurrence can recover an alternative enzyme excluded from the displayed singleton"() {
        given:
        def state = new IntegrationState()
        state.orthogroupMap = [bystander: 'common_bystander', enzyme: 'alternative_family']
        def result = new IntegratedAnnotationSet(suggestions: [new SingletonSuggestion(
            source: 'reaction_local_context', proteinId: 'bystander', q: 0.6d,
            reactionId: 'missing_reaction', functionId: '1.1.1.1', functionType: AnnotationType.EC,
            proteinScores: [bystander: new PerProteinDecomposition(proteinId: 'bystander', q: 0.6d),
                            enzyme: new PerProteinDecomposition(proteinId: 'enzyme', q: 0.4d)])])

        when:
        new CrossGenomeReScorer(requireCredible: false).rescore(state, result, negativeCatalog())

        then:
        result.suggestions.first() instanceof SingletonSuggestion
        result.suggestions.first().proteinId == 'enzyme'
        result.suggestions.first().q < 1d
        result.suggestions.first().candidateProbabilities.keySet() == ['enzyme', 'bystander'] as Set
    }

    def "one conserved observation is uncertain even when every observed genome agrees"() {
        given:
        def small = new ReactionLocusCatalog.Entry(nSigWith: 1, nSigTotal: 1, nBaseWith: 1, nBaseTotal: 1)
        def large = new ReactionLocusCatalog.Entry(nSigWith: 100, nSigTotal: 100, nBaseWith: 100, nBaseTotal: 100)

        expect:
        small.logLRCiHalfWidth() > large.logLRCiHalfWidth()
        large.logLRCiHalfWidth() > 0d
        Double.isFinite(small.logLRCiHalfWidth())
    }

    def "a pathway suggester cannot erase reaction-network evidence when pathway inputs are absent"() {
        given:
        def reactionSuggestion = new SingletonSuggestion(source: 'reaction_local_context',
            proteinId: 'enzyme', q: 0.7d, reactionId: 'missing_reaction')
        def result = new IntegratedAnnotationSet(suggestions: [reactionSuggestion])

        when:
        new DarkMatterSuggester().suggest(new IntegrationState(), result)

        then:
        result.suggestions == [reactionSuggestion]
    }
}
