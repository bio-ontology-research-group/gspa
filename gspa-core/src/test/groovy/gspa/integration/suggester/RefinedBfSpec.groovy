package gspa.integration.suggester

import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.integration.MetabolicGap
import gspa.model.Contig
import gspa.model.Genome
import gspa.model.Protein
import gspa.ontology.PathwayDatabase
import gspa.ontology.PathwayGraph
import spock.lang.Specification

class RefinedBfSpec extends Specification {

    private static PathwayDatabase oneBigPathway() {
        def db = new PathwayDatabase()
        def pw = new PathwayGraph(pathwayId: 'PWY1', pathwayName: 'P1')
        pw.metaClass.getRequiredGoTerms = { ->
            ['GO:A', 'GO:B', 'GO:C', 'GO:D', 'GO:E'] as Set
        }
        db.pathways['PWY1'] = pw
        db
    }

    private static Genome genomeWithProteins(List<String> ids) {
        new Genome(id: 'g', contigs: [new Contig(id: 'c1', proteins: ids.collect { new Protein(id: it) })])
    }

    private static IntegrationState stateWith(Map<String, Double> logOdds, List<String> proteinIds) {
        def st = new IntegrationState(genomeWithProteins(proteinIds))
        st.pathwayDatabase = oneBigPathway()
        logOdds.each { k, v -> st.set(k, v) }
        st
    }

    def "refined BF rewards diversity: 4 distinct GO hits > 4 hits on same GO"() {
        given:
        def proteinsA = ['pA1', 'pA2', 'pA3', 'pA4']
        def proteinsB = ['pB1', 'pB2', 'pB3', 'pB4']

        def state = stateWith([
            // Operon A: each member hits a distinct pathway function.
            'pA1|GO|GO:A': 3.0d, 'pA2|GO|GO:B': 3.0d, 'pA3|GO|GO:C': 3.0d, 'pA4|GO|GO:D': 3.0d,
            // Operon B: all four members hit the SAME pathway function (paralog expansion).
            'pB1|GO|GO:A': 3.0d, 'pB2|GO|GO:A': 3.0d, 'pB3|GO|GO:A': 3.0d, 'pB4|GO|GO:A': 3.0d,
        ], proteinsA + proteinsB)

        def suggester = new DarkMatterSuggester()
        suggester.useRefinedBayesFactor = true

        // Build a gap; we need a pathway terms set matching functions(P).
        def pfTerms = ['GO:A', 'GO:B', 'GO:C', 'GO:D', 'GO:E'] as Set
        // Compute base rates once.
        def m = suggester.class.getDeclaredMethod('computeBaseRates', IntegrationState); m.accessible = true
        def baseRates = m.invoke(suggester, state)

        def bfMethod = suggester.class.getDeclaredMethod(
            'computeRefinedBayesFactor', List, Set, IntegrationState, Map)
        bfMethod.accessible = true

        when:
        double bfA = bfMethod.invoke(suggester, proteinsA, pfTerms, state, baseRates)
        double bfB = bfMethod.invoke(suggester, proteinsB, pfTerms, state, baseRates)

        then: "diverse operon A scores strictly higher than paralog-clump operon B"
        bfA > bfB
    }

    def "refined BF penalizes off-pathway strong annotations (purity)"() {
        given:
        def proteinsA = ['pA1', 'pA2', 'pA3']
        def proteinsB = ['pB1', 'pB2', 'pB3']

        def state = stateWith([
            // Operon A: three on-pathway annotations.
            'pA1|GO|GO:A': 3.0d, 'pA2|GO|GO:B': 3.0d, 'pA3|GO|GO:C': 3.0d,
            // Operon B: same three on-pathway annotations PLUS strong off-pathway noise.
            'pB1|GO|GO:A': 3.0d, 'pB2|GO|GO:B': 3.0d, 'pB3|GO|GO:C': 3.0d,
            'pB1|GO|GO:Z1': 5.0d, 'pB2|GO|GO:Z2': 5.0d, 'pB3|GO|GO:Z3': 5.0d,
        ], proteinsA + proteinsB)

        def suggester = new DarkMatterSuggester()
        suggester.useRefinedBayesFactor = true

        def pfTerms = ['GO:A', 'GO:B', 'GO:C', 'GO:D', 'GO:E'] as Set
        def m = suggester.class.getDeclaredMethod('computeBaseRates', IntegrationState); m.accessible = true
        def baseRates = m.invoke(suggester, state)
        def bfMethod = suggester.class.getDeclaredMethod(
            'computeRefinedBayesFactor', List, Set, IntegrationState, Map); bfMethod.accessible = true

        when:
        double bfA = bfMethod.invoke(suggester, proteinsA, pfTerms, state, baseRates)
        double bfB = bfMethod.invoke(suggester, proteinsB, pfTerms, state, baseRates)

        then: "operon A (pure pathway) scores strictly higher than B (contaminated)"
        bfA > bfB
    }

    def "refined BF weights rare pathway functions more heavily (IC)"() {
        given: "two operons each hit one pathway function; one is rare, one common"
        // 10 proteins total; GO:A annotated on 5 of them (common), GO:B on 1 (rare).
        def proteins = (1..10).collect { "p${it}".toString() }
        Map<String, Double> odds = [:]
        // GO:A on 5 proteins (base rate 0.5)
        (1..5).each { odds["p${it}|GO|GO:A".toString()] = 3.0d }
        // GO:B on 1 protein (base rate 0.1)
        odds['p10|GO|GO:B'] = 3.0d

        def state = stateWith(odds, proteins)

        def suggester = new DarkMatterSuggester()
        suggester.useRefinedBayesFactor = true
        suggester.purityWeight = 0.0d     // isolate the IC effect from purity

        def pfTerms = ['GO:A', 'GO:B'] as Set
        def m = suggester.class.getDeclaredMethod('computeBaseRates', IntegrationState); m.accessible = true
        def baseRates = m.invoke(suggester, state)
        def bfMethod = suggester.class.getDeclaredMethod(
            'computeRefinedBayesFactor', List, Set, IntegrationState, Map); bfMethod.accessible = true

        // Operon covering only the COMMON function A (via p1).
        // Operon covering only the RARE function B (via p10).
        when:
        double bfCommon = bfMethod.invoke(suggester, ['p1'], pfTerms, state, baseRates)
        double bfRare   = bfMethod.invoke(suggester, ['p10'], pfTerms, state, baseRates)

        then: "hitting the rare function yields a larger BF than hitting the common one"
        bfRare > bfCommon
    }
}
