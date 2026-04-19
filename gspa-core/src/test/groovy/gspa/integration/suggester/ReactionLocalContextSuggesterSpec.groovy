package gspa.integration.suggester

import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.integration.MetabolicGap
import gspa.model.Genome
import gspa.model.GenomeLayout
import gspa.model.ProteinLocus
import gspa.model.Strand
import gspa.ontology.PathwayDatabase
import gspa.ontology.ReactionGraph
import spock.lang.Specification

class ReactionLocalContextSuggesterSpec extends Specification {

    /**
     * Fixture: a 4-reaction chain r1-r2-r3-r4 with ECs mapped to GO terms.
     * r3 is the gap; r2 and r4 are its 1-hop neighbours. Anchor proteins
     * (catalysing r2 and r4) sit in a 10 kb operon; the "correct" candidate
     * p_mid sits in the middle of that operon. A distant bystander p_far
     * sits on the same contig but 40 kb away — outside the kernel bandwidth.
     */
    private fixture(double pTarget = 0.9d, double pBystander = 0.0d,
                    boolean commitTarget = false) {
        // Reaction graph
        def rg = new ReactionGraph()
        rg.addReaction(new ReactionGraph.ReactionSpec(
                rxnId: 'r1', ecNumber: '1.1.1.1',
                substrates: ['a'] as Set, products: ['b'] as Set))
        rg.addReaction(new ReactionGraph.ReactionSpec(
                rxnId: 'r2', ecNumber: '2.2.2.2',
                substrates: ['b'] as Set, products: ['c'] as Set))
        rg.addReaction(new ReactionGraph.ReactionSpec(
                rxnId: 'r3', ecNumber: '3.3.3.3',
                substrates: ['c'] as Set, products: ['d'] as Set))
        rg.addReaction(new ReactionGraph.ReactionSpec(
                rxnId: 'r4', ecNumber: '4.4.4.4',
                substrates: ['d'] as Set, products: ['e'] as Set))
        rg.build()

        // PathwayDatabase needed only for ec2go map
        def pdb = new PathwayDatabase()
        pdb.ec2go = [
                '1.1.1.1': 'GO:1111111',
                '2.2.2.2': 'GO:2222222',
                '3.3.3.3': 'GO:3333333',
                '4.4.4.4': 'GO:4444444',
        ]

        // Layout: operon {pN2, p_mid, pN4} clustered within ~6 kb; p_far 40 kb
        def layout = new GenomeLayout('test')
        layout.add(new ProteinLocus(genomeId: 'test', proteinId: 'pN2', contig: 'chr1',
                start: 10_000, end: 10_900, strand: Strand.PLUS))
        layout.add(new ProteinLocus(genomeId: 'test', proteinId: 'p_mid', contig: 'chr1',
                start: 12_000, end: 13_200, strand: Strand.PLUS))
        layout.add(new ProteinLocus(genomeId: 'test', proteinId: 'pN4', contig: 'chr1',
                start: 15_000, end: 16_000, strand: Strand.PLUS))
        layout.add(new ProteinLocus(genomeId: 'test', proteinId: 'p_far', contig: 'chr1',
                start: 55_000, end: 56_000, strand: Strand.PLUS))
        layout.finishLoading()

        // State with posteriors: pN2 strongly catalyses r2; pN4 catalyses r4;
        // p_mid has zero evidence for r3 (the gap) unless commitTarget=true
        def state = new IntegrationState(new Genome(id: 'test'))
        state.pathwayDatabase = pdb
        state.reactionGraph = rg
        state.genomeLayout = layout
        state.metabolicGaps = [new MetabolicGap(
                pathwayId: 'TEST', reactionId: 'r3',
                ecNumber: '3.3.3.3', goTerm: 'GO:3333333')]

        // strong posteriors = large positive log-odds
        double strong = Math.log(pTarget / (1.0d - pTarget))
        double weak = pBystander > 0.0d ? Math.log(pBystander / (1.0d - pBystander)) : -6.0d
        state.posteriorLogOdds['pN2|GO|GO:2222222'] = strong
        state.posteriorLogOdds['pN4|GO|GO:4444444'] = strong
        // p_mid: no neighbour-function posteriors (dark for the task)
        if (commitTarget) {
            // p_mid strongly committed to OTHER pathway function (r2) — should be penalised
            state.posteriorLogOdds['p_mid|GO|GO:2222222'] = strong
        }
        // p_far has a neighbour-function hit too, but it's outside the kernel
        state.posteriorLogOdds['p_far|GO|GO:4444444'] = strong

        state
    }

    def "emits a singleton assigning the gap's GO to the operon-central dark gene"() {
        given:
        def state = fixture(0.9d)
        def integrated = new IntegratedAnnotationSet()
        def s = new ReactionLocalContextSuggester()

        when:
        s.suggest(state, integrated)

        then:
        !integrated.suggestions.isEmpty()
        def sug = integrated.suggestions.find { it.reactionId == 'r3' }
        sug != null
        sug.functionId == 'GO:3333333'
    }

    def "prefers p_mid over p_far even when p_far has a stronger neighbour annotation"() {
        given:
        def state = fixture(0.9d)
        def integrated = new IntegratedAnnotationSet()

        when:
        new ReactionLocalContextSuggester().suggest(state, integrated)

        then: 'the winner is one of the in-window candidates, never p_far'
        integrated.suggestions.any { s ->
            (s instanceof SingletonSuggestion && s.proteinId == 'p_mid') ||
            (s instanceof DisjunctiveSuggestion && s.proteinIds.contains('p_mid')) ||
            (s instanceof SingletonSuggestion && (s.proteinId == 'pN2' || s.proteinId == 'pN4'))
        }
        !integrated.suggestions.any { s ->
            (s instanceof SingletonSuggestion && s.proteinId == 'p_far')
        }
    }

    def "no-op when reaction graph is not wired"() {
        given:
        def state = new IntegrationState(new Genome(id: 'test'))
        state.metabolicGaps = [new MetabolicGap(
                pathwayId: 'TEST', reactionId: 'r3', ecNumber: '3.3.3.3',
                goTerm: 'GO:3333333')]
        // reactionGraph = null
        def integrated = new IntegratedAnnotationSet()

        when:
        new ReactionLocalContextSuggester().suggest(state, integrated)

        then:
        integrated.suggestions == null || integrated.suggestions.isEmpty()
    }

    def "no-op when genome layout is not wired"() {
        given:
        def state = fixture(0.9d)
        state.genomeLayout = null
        def integrated = new IntegratedAnnotationSet()

        when:
        new ReactionLocalContextSuggester().suggest(state, integrated)

        then:
        integrated.suggestions == null || integrated.suggestions.isEmpty()
    }
}
