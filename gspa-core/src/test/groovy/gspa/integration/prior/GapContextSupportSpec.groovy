package gspa.integration.prior

import gspa.integration.IntegrationState
import gspa.integration.MetabolicGap
import gspa.integration.cohort.CohortInput
import gspa.integration.cohort.ContextMode
import gspa.model.Genome
import gspa.model.GenomeLayout
import gspa.model.ProteinLocus
import gspa.model.Strand
import gspa.ontology.ReactionGraph
import spock.lang.Specification

class GapContextSupportSpec extends Specification {

    private static double logOdds(double p) {
        Math.log(p / (1.0d - p))
    }

    /**
     * Synthetic fixture:
     *   Graph: r0 (1.1.1.1) -- r1 (2.2.2.2) -- r2 (3.3.3.3), plus isolated r9 (9.9.9.9)
     *   Layout on chr1:
     *     a1   at  1,000..1,900  (midpoint  1,450)
     *     cNear at 2,500..3,400  (midpoint  2,950; 1,500 bp from a1)
     *     cFar  at 52,000..52,900 (midpoint 52,450; >50 kb from a1, 1,500 bp from a2)
     *     a2   at 53,500..54,400 (midpoint 53,950)
     */
    private IntegrationState buildState() {
        def rg = new ReactionGraph()
        rg.addReaction(new ReactionGraph.ReactionSpec(
            rxnId: 'r0', ecNumber: '1.1.1.1',
            substrates: ['m0'] as Set, products: ['m1'] as Set))
        rg.addReaction(new ReactionGraph.ReactionSpec(
            rxnId: 'r1', ecNumber: '2.2.2.2',
            substrates: ['m1'] as Set, products: ['m2'] as Set))
        rg.addReaction(new ReactionGraph.ReactionSpec(
            rxnId: 'r2', ecNumber: '3.3.3.3',
            substrates: ['m2'] as Set, products: ['m3'] as Set))
        rg.addReaction(new ReactionGraph.ReactionSpec(
            rxnId: 'r9', ecNumber: '9.9.9.9',
            substrates: ['m8'] as Set, products: ['m9'] as Set))
        rg.bindEc('r0', '1.1.1.1')
        rg.bindEc('r1', '2.2.2.2')
        rg.bindEc('r2', '3.3.3.3')
        rg.bindEc('r9', '9.9.9.9')
        rg.build()

        def layout = new GenomeLayout('g1')
        layout.add(new ProteinLocus(genomeId: 'g1', proteinId: 'a1', contig: 'chr1', start: 1000, end: 1900, strand: Strand.PLUS))
        layout.add(new ProteinLocus(genomeId: 'g1', proteinId: 'cNear', contig: 'chr1', start: 2500, end: 3400, strand: Strand.PLUS))
        layout.add(new ProteinLocus(genomeId: 'g1', proteinId: 'cFar', contig: 'chr1', start: 52000, end: 52900, strand: Strand.PLUS))
        layout.add(new ProteinLocus(genomeId: 'g1', proteinId: 'a2', contig: 'chr1', start: 53500, end: 54400, strand: Strand.PLUS))
        layout.finishLoading()

        def state = new IntegrationState(new Genome(id: 'g1'))
        state.reactionGraph = rg
        state.genomeLayout = layout
        state.ecToGoAlternatives = [
            '1.1.1.1': ['GO:0001111'],
            '2.2.2.2': ['GO:0002222'],
            '3.3.3.3': ['GO:0003333'],
            '9.9.9.9': ['GO:0009999'],
        ]
        state
    }

    def "(a) UNIVERSAL: candidate next to an anchor of r1 gets c = 1 when r1 is the only anchored context reaction; candidate 50 kb away gets 0"() {
        given:
        def state = buildState()
        state.posteriorLogOdds['a1|GO|GO:0002222'] = logOdds(0.9d)
        def gap = new MetabolicGap(pathwayId: 'P', reactionId: 'r0', ecNumber: 'EC:1.1.1.1', goTerm: 'GO:0001111')
        def support = new GapContextSupport(mode: ContextMode.UNIVERSAL, windowBp: 10000)

        when:
        Map<String, Double> c = support.contextSupport(gap, ['cNear', 'cFar'], state)

        then:
        c['cNear'] == 1.0d
        c['cFar'] == 0.0d
        support.lastReport.gaps_evaluated == 1
        support.lastReport.gaps_with_context == 1
        support.lastReport.candidates_with_c_gt_0 == 1
    }

    def "(b) weights: with r1 (alpha) and r2 (alpha^2) both anchored and only r2's anchor near the candidate, c = alpha^2 / (alpha + alpha^2)"() {
        given:
        def state = buildState()
        state.posteriorLogOdds['a1|GO|GO:0002222'] = logOdds(0.9d)
        state.posteriorLogOdds['a2|GO|GO:0003333'] = logOdds(0.9d)
        def gap = new MetabolicGap(pathwayId: 'P', reactionId: 'r0', ecNumber: 'EC:1.1.1.1', goTerm: 'GO:0001111')
        double alpha = 0.5d
        def support = new GapContextSupport(mode: ContextMode.UNIVERSAL, windowBp: 10000, alpha: alpha, maxHops: 2)

        when:
        Map<String, Double> c = support.contextSupport(gap, ['cFar'], state)
        double expected = (alpha * alpha) / (alpha + alpha * alpha)

        then:
        Math.abs(c['cFar'] - expected) < 1e-12
    }

    def "(c) ORGANISM: with r1 excluded from the organism set, r2 is unreachable (no traversal through excluded intermediates) -> empty map when r2 was the only anchored one"() {
        given:
        def state = buildState()
        state.posteriorLogOdds['a2|GO|GO:0003333'] = logOdds(0.9d)
        def gap = new MetabolicGap(pathwayId: 'P', reactionId: 'r0', ecNumber: 'EC:1.1.1.1', goTerm: 'GO:0001111')
        def supportExcluded = new GapContextSupport(
            mode: ContextMode.ORGANISM,
            windowBp: 10000,
            organismReactions: ['r2'] as Set
        )
        def supportIncluded = new GapContextSupport(
            mode: ContextMode.ORGANISM,
            windowBp: 10000,
            organismReactions: ['r1', 'r2'] as Set
        )

        expect:
        supportExcluded.contextSupport(gap, ['cNear', 'cFar'], state).isEmpty()
        supportIncluded.contextSupport(gap, ['cFar'], state)['cFar'] == 1.0d
    }

    def "(d) PATHWAY: reactions of the same pathway id in other groups count, own-group reactions do not"() {
        given:
        def state = buildState()
        state.posteriorLogOdds['a1|GO|GO:0002222'] = logOdds(0.9d)
        state.posteriorLogOdds['a2|GO|GO:0009999'] = logOdds(0.9d)
        def gap = new MetabolicGap(pathwayId: 'P', reactionId: 'r0', ecNumber: 'EC:1.1.1.1', goTerm: 'GO:0001111')

        and: 'r1 is in r0 own group G1; r9 is in another group G2 of the same pathway P'
        List<CohortInput.PathwayGroup> groups = [
            new CohortInput.PathwayGroup(pathwayId: 'P', groupId: 'G1', source: 'test', reactionIds: ['r0', 'r1'] as Set),
            new CohortInput.PathwayGroup(pathwayId: 'P', groupId: 'G2', source: 'test', reactionIds: ['r9'] as Set),
        ]
        def support = new GapContextSupport(mode: ContextMode.PATHWAY, windowBp: 10000, pathwayGroups: groups)

        when:
        Map<String, Double> c = support.contextSupport(gap, ['cNear', 'cFar'], state)

        then: 'only r9 (anchored by a2 near cFar) counts; own-group r1 (anchored by a1 near cNear) does not'
        c['cNear'] == 0.0d
        c['cFar'] == 1.0d

        and: 'when only own-group reaction r1 is anchored, no anchored context reaction remains -> empty map'
        state.posteriorLogOdds.remove('a2|GO|GO:0009999')
        support.contextSupport(gap, ['cNear', 'cFar'], state).isEmpty()
    }

    def "(e) a gap whose EC maps to no graph reaction -> empty map"() {
        given:
        def state = buildState()
        state.posteriorLogOdds['a1|GO|GO:0002222'] = logOdds(0.9d)
        def gap = new MetabolicGap(pathwayId: 'P', reactionId: 'unknown_rxn', ecNumber: 'EC:5.5.5.5', goTerm: 'GO:0005555')
        def support = new GapContextSupport(mode: ContextMode.UNIVERSAL, windowBp: 10000)

        expect:
        support.contextSupport(gap, ['cNear', 'cFar'], state).isEmpty()
    }

    def "(f) the gap's own GO term never anchors"() {
        given:
        def state = buildState()
        // Map r1's EC (2.2.2.2) to the gap's own GO term (GO:0001111) only
        state.ecToGoAlternatives['2.2.2.2'] = ['GO:0001111']
        state.posteriorLogOdds['a1|GO|GO:0001111'] = logOdds(0.95d)
        def gap = new MetabolicGap(pathwayId: 'P', reactionId: 'r0', ecNumber: 'EC:1.1.1.1', goTerm: 'GO:0001111')
        def support = new GapContextSupport(mode: ContextMode.UNIVERSAL, windowBp: 10000)

        expect:
        support.contextSupport(gap, ['cNear', 'cFar'], state).isEmpty()
    }
}
