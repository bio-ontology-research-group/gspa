package gspa.integration.suggester

import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.integration.MetabolicGap
import gspa.model.Genome
import gspa.ontology.PathwayDatabase
import gspa.ontology.PathwayGraph
import spock.lang.Specification

class DarkMatterSuggesterSpec extends Specification {

    /** Build a PathwayDatabase with two pathways whose required-term sets we control. */
    private static PathwayDatabase pathwayDb() {
        def db = new PathwayDatabase()

        def ribosome = new PathwayGraph(pathwayId: 'PWY-RIBO', pathwayName: 'Ribosome')
        ribosome.metaClass.getRequiredGoTerms = { ->
            ['GO:0006412', 'GO:0003735', 'GO:0019843', 'GO:0000028'] as Set
        }
        db.pathways['PWY-RIBO'] = ribosome

        def glycolysis = new PathwayGraph(pathwayId: 'PWY-GLYCO', pathwayName: 'Glycolysis')
        glycolysis.metaClass.getRequiredGoTerms = { ->
            ['GO:0006096', 'GO:0004340', 'GO:0004618'] as Set
        }
        db.pathways['PWY-GLYCO'] = glycolysis

        db
    }

    /** Make a state with the given operon + gap, plus a function to seed posteriors. */
    private static IntegrationState basicState(List<List<String>> operons, List<MetabolicGap> gaps) {
        def state = new IntegrationState(new Genome(id: 'test'))
        state.pathwayDatabase = pathwayDb()
        state.operons = operons
        state.metabolicGaps = gaps
        state
    }

    // ------------------------------------------------------------------
    // Smoke / no-op cases
    // ------------------------------------------------------------------

    def "no-op when no metabolic gaps"() {
        given:
        def state = basicState([['p1', 'p2']], [])
        def integrated = new IntegratedAnnotationSet()
        def suggester = new DarkMatterSuggester()

        when:
        suggester.suggest(state, integrated)

        then:
        integrated.suggestions.isEmpty()
    }

    def "no-op when no operons wired"() {
        given:
        def state = basicState([], [new MetabolicGap(pathwayId: 'PWY-RIBO', goTerm: 'GO:0003735')])
        def integrated = new IntegratedAnnotationSet()

        when:
        new DarkMatterSuggester().suggest(state, integrated)

        then:
        integrated.suggestions.isEmpty()
    }

    def "no-op when no pathway database wired"() {
        given:
        def state = new IntegrationState(new Genome(id: 't'))
        state.operons = [['p1', 'p2']]
        state.metabolicGaps = [new MetabolicGap(pathwayId: 'PWY-X', goTerm: 'GO:0000001')]
        // pathwayDatabase = null
        def integrated = new IntegratedAnnotationSet()

        when:
        new DarkMatterSuggester().suggest(state, integrated)

        then:
        integrated.suggestions.isEmpty()
    }

    // ------------------------------------------------------------------
    // Singleton case
    // ------------------------------------------------------------------

    def "operon with strong pathway consensus + single dark protein emits singleton suggestion"() {
        given: "an operon {p1,p2,p3,p4} where p1/p2/p3 have strong ribosome annotations and p4 is dark"
        def state = basicState([['p1', 'p2', 'p3', 'p4']], [
            new MetabolicGap(
                pathwayId: 'PWY-RIBO',
                reactionId: 'RXN-28S',
                ecNumber: null,
                goTerm: 'GO:0000028',   // the gap: 28S rRNA something
            )
        ])
        // 3 strong ribosome annotations push BF very high.
        state.set('p1|GO|GO:0006412', 5.0d)   // translation
        state.set('p2|GO|GO:0003735', 5.0d)   // structural ribosome
        state.set('p3|GO|GO:0019843', 5.0d)   // rRNA binding
        // p4 has no annotations at all → "dark".

        def integrated = new IntegratedAnnotationSet()
        def suggester = new DarkMatterSuggester()

        when:
        suggester.suggest(state, integrated)

        then: "exactly one suggestion is emitted"
        integrated.suggestions.size() == 1

        and: "it's a singleton pointing at p4 (the dark protein)"
        integrated.suggestions[0] instanceof SingletonSuggestion
        def ss = integrated.suggestions[0] as SingletonSuggestion
        ss.proteinId == 'p4'
        ss.functionId == 'GO:0000028'
        ss.pathwayId == 'PWY-RIBO'
        ss.bayesFactor >= new DarkMatterSuggester().bfMin
        ss.q > 0.5d
        ss.suggestionScore > 0
        ss.suggestionScore <= 0.85d
    }

    // ------------------------------------------------------------------
    // Non-dark protein still eligible
    // ------------------------------------------------------------------

    def "a protein with unrelated annotations is still a valid suggestion target"() {
        given: "operon with p1/p2 strongly ribosomal, p3 already annotated for an unrelated function"
        def state = basicState([['p1', 'p2', 'p3']], [
            new MetabolicGap(pathwayId: 'PWY-RIBO', goTerm: 'GO:0003735')
        ])
        state.set('p1|GO|GO:0006412', 5.0d)
        state.set('p2|GO|GO:0019843', 5.0d)
        // p3 is annotated for an UNRELATED GO term — not dark, just not ribosomal.
        state.set('p3|GO|GO:9999999', 5.0d)

        def integrated = new IntegratedAnnotationSet()

        when:
        new DarkMatterSuggester().suggest(state, integrated)

        then: "p3 still receives a suggestion because the operon evidence dominates"
        integrated.suggestions.size() == 1
        def s = integrated.suggestions[0] as SingletonSuggestion
        s.proteinId == 'p3'
        s.functionId == 'GO:0003735'
    }

    // ------------------------------------------------------------------
    // Disjunctive case
    // ------------------------------------------------------------------

    def "two dark proteins in a strong operon produce a disjunctive suggestion"() {
        given:
        def state = basicState([['p1', 'p2', 'p3', 'p4']], [
            new MetabolicGap(pathwayId: 'PWY-RIBO', goTerm: 'GO:0000028')
        ])
        // Two strong ribosome annotations on p1 and p2, p3 and p4 both dark.
        state.set('p1|GO|GO:0006412', 5.0d)
        state.set('p2|GO|GO:0003735', 5.0d)

        def integrated = new IntegratedAnnotationSet()

        when:
        new DarkMatterSuggester().suggest(state, integrated)

        then:
        integrated.suggestions.size() == 1

        and: "it's a disjunctive suggestion over the two dark proteins"
        integrated.suggestions[0] instanceof DisjunctiveSuggestion
        def ds = integrated.suggestions[0] as DisjunctiveSuggestion
        ds.proteinIds.toSet() == ['p3', 'p4'] as Set
        // Both dark proteins have identical L_R so q is 0.5 each (tie-broken by sort).
        ds.qValues.every { Math.abs(it - 0.5d) < 0.1d }
    }

    // ------------------------------------------------------------------
    // BF_min filter
    // ------------------------------------------------------------------

    def "operon with no pathway-matching annotations yields BF=1 and is rejected"() {
        given:
        def state = basicState([['p1', 'p2']], [
            new MetabolicGap(pathwayId: 'PWY-RIBO', goTerm: 'GO:0003735')
        ])
        // Neither protein has any ribosome annotation; both are annotated for unrelated stuff.
        state.set('p1|GO|GO:0009999', 5.0d)
        state.set('p2|GO|GO:0008888', 5.0d)

        def integrated = new IntegratedAnnotationSet()
        def suggester = new DarkMatterSuggester()

        when:
        suggester.suggest(state, integrated)

        then: "no suggestions because BF = exp(0) = 1 < bfMin"
        integrated.suggestions.isEmpty()
    }

    // ------------------------------------------------------------------
    // Suggestion score caps and shape
    // ------------------------------------------------------------------

    def "suggestionScore never exceeds the configured cap"() {
        given:
        def state = basicState([['p1', 'p2', 'p3', 'p4']], [
            new MetabolicGap(pathwayId: 'PWY-RIBO', goTerm: 'GO:0000028')
        ])
        // All four proteins strongly ribosome-annotated → huge BF.
        state.set('p1|GO|GO:0006412', 10.0d)
        state.set('p2|GO|GO:0003735', 10.0d)
        state.set('p3|GO|GO:0019843', 10.0d)

        def integrated = new IntegratedAnnotationSet()

        when:
        new DarkMatterSuggester().suggest(state, integrated)

        then:
        integrated.suggestions.every { it.suggestionScore <= 0.85d }
    }

    // ------------------------------------------------------------------
    // Provenance
    // ------------------------------------------------------------------

    def "provenance records bayes factor + operon + gap identifiers"() {
        given:
        def state = basicState([['p1', 'p2']], [
            new MetabolicGap(pathwayId: 'PWY-GLYCO', reactionId: 'RXN-HK', goTerm: 'GO:0004340')
        ])
        state.set('p1|GO|GO:0006096', 5.0d)

        def integrated = new IntegratedAnnotationSet()

        when:
        new DarkMatterSuggester().suggest(state, integrated)

        then:
        integrated.suggestions.size() == 1
        def s = integrated.suggestions[0]
        s.provenance.contains('PWY-GLYCO')
        s.provenance.contains('RXN-HK')
        s.proteinScores.size() == 2   // every operon member's decomposition
        s.proteinScores.values().every { it.operonLogOdds != 0 }
    }
}
