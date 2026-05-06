package gspa.integration.crossgenome

import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.integration.suggester.DisjunctiveSuggestion
import gspa.integration.suggester.SingletonSuggestion
import gspa.model.AnnotationType
import gspa.model.Genome
import spock.lang.Specification

class ConditionalLRScorerSpec extends Specification {

    private IntegrationState state(Map<String, String> orthogroups) {
        def s = new IntegrationState(new Genome(id: 'test'))
        s.orthogroupMap = orthogroups
        s
    }

    private ReactionLocusCatalog catalog(List<ReactionLocusCatalog.Entry> entries, int panelSize) {
        def cat = new ReactionLocusCatalog(panelSize: panelSize)
        entries.each { cat.put(it) }
        cat
    }

    def "pass-through when no orthogroup map or empty catalog"() {
        given:
        def s = new IntegrationState(new Genome(id: 't'))
        def integ = new IntegratedAnnotationSet(suggestions: [
            new DisjunctiveSuggestion(
                reactionId: 'r1', proteinIds: ['pA', 'pB'],
                qValues: [0.5d, 0.5d],
                functionId: 'GO:0001', functionType: AnnotationType.GO),
        ])
        def scorer = new CrossGenomeReScorer()
        def cat = new ReactionLocusCatalog()

        when:
        scorer.rescore(s, integ, cat)

        then: 'suggestions list unchanged'
        integ.suggestions.size() == 1
    }

    def "boosts orthogroup C that appears in 5/5 R-signature windows over bystander in 1/5"() {
        given: 'two candidates tied within-genome; C is panel-conserved, B is not'
        def s = state([pA: 'C', pB: 'B'])
        def integ = new IntegratedAnnotationSet(suggestions: [
            new DisjunctiveSuggestion(
                reactionId: 'r1', proteinIds: ['pA', 'pB'],
                qValues: [0.5d, 0.5d],
                functionId: 'GO:0001', functionType: AnnotationType.GO),
        ])
        def cat = catalog([
            // C present in 5/5 R-windows out of 5 R-signature genomes; base present in all 5
            new ReactionLocusCatalog.Entry(orthogroupId: 'C', reactionId: 'r1',
                    nSigWith: 5, nSigTotal: 5, nBaseWith: 5, nBaseTotal: 5),
            // B present in 1/5 R-windows; base present in 5
            new ReactionLocusCatalog.Entry(orthogroupId: 'B', reactionId: 'r1',
                    nSigWith: 1, nSigTotal: 5, nBaseWith: 5, nBaseTotal: 5),
        ], 5)
        def scorer = new CrossGenomeReScorer()
        scorer.minSupport = 3
        scorer.requireCredible = false   // low N would otherwise fail the CI test

        when:
        scorer.rescore(s, integ, cat)

        then: 'a singleton or disjunctive remains; pA beats pB'
        def s0 = integ.suggestions[0]
        if (s0 instanceof SingletonSuggestion) {
            s0.proteinId == 'pA'
        } else {
            s0.proteinIds[0] == 'pA'
            s0.qValues[0] > s0.qValues[1]
        }
    }

    def "drops LR if cg-min-support not met"() {
        given:
        def s = state([pA: 'C', pB: 'B'])
        def integ = new IntegratedAnnotationSet(suggestions: [
            new DisjunctiveSuggestion(
                reactionId: 'r1', proteinIds: ['pA', 'pB'],
                qValues: [0.5d, 0.5d],
                functionId: 'GO:0001', functionType: AnnotationType.GO),
        ])
        def cat = catalog([
            // would-be strong LR but only 2 panel genomes
            new ReactionLocusCatalog.Entry(orthogroupId: 'C', reactionId: 'r1',
                    nSigWith: 2, nSigTotal: 2, nBaseWith: 5, nBaseTotal: 5),
        ], 5)
        def scorer = new CrossGenomeReScorer()
        scorer.minSupport = 5
        scorer.requireCredible = false

        when:
        scorer.rescore(s, integ, cat)

        then: 'unchanged — min-support gate excludes the LR'
        integ.suggestions.size() == 1
        integ.suggestions[0] instanceof DisjunctiveSuggestion
        integ.suggestions[0].qValues[0] == integ.suggestions[0].qValues[1]
    }

    def "writes and reads ReactionLocusCatalog round-trip"() {
        given:
        def cat = new ReactionLocusCatalog(panelSize: 5)
        cat.put(new ReactionLocusCatalog.Entry(
            orthogroupId: 'C1', reactionId: 'rA',
            nSigWith: 4, nSigTotal: 5, nBaseWith: 5, nBaseTotal: 5))
        cat.put(new ReactionLocusCatalog.Entry(
            orthogroupId: 'C2', reactionId: 'rB',
            nSigWith: 1, nSigTotal: 4, nBaseWith: 3, nBaseTotal: 5))
        def tmp = File.createTempFile('catalog', '.tsv')
        tmp.deleteOnExit()

        when:
        cat.writeTo(tmp)
        def read = ReactionLocusCatalog.readFrom(tmp)

        then:
        read.panelSize == 5
        read.size() == 2
        read.get('C1', 'rA').nSigWith == 4
        read.get('C2', 'rB').nSigTotal == 4
    }
}
