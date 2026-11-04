package gspa.integration.prior

import gspa.integration.IntegrationState
import gspa.integration.MetabolicGap
import gspa.integration.cohort.ContextMode
import gspa.model.AnnotationType
import gspa.model.Genome
import gspa.model.GenomeLayout
import gspa.model.ProteinLocus
import gspa.model.Strand
import gspa.ontology.ReactionGraph
import spock.lang.Specification

class GapFillingPriorSpec extends Specification {

    def "no-ops when no metabolic gaps are present"() {
        given:
        def state = new IntegrationState(new Genome(id: 'g'))
        def prior = new GapFillingPrior()

        when:
        prior.beginIteration(state)

        then:
        prior.logOddsBoost('p1', 'p1|EC|EC:2.7.1.1', state) == 0.0d
    }

    def "boosts EC claims that match a metabolic gap's EC number"() {
        given:
        def state = new IntegrationState(new Genome(id: 'g'))
        state.metabolicGaps = [
            new MetabolicGap(
                pathwayId: 'PWY-12345',
                reactionId: 'RXN-1',
                ecNumber: 'EC:2.7.1.1',
                goTerm: 'GO:0004396',
            )
        ]
        def prior = new GapFillingPrior()
        prior.alphaGap = 1.5d

        when:
        prior.beginIteration(state)

        then:
        prior.logOddsBoost('p1', 'p1|EC|EC:2.7.1.1', state) == 1.5d
        prior.logOddsBoost('p1', 'p1|GO|GO:0004396', state) == 1.5d
        prior.logOddsBoost('p1', 'p1|EC|EC:9.9.9.9', state) == 0.0d
    }

    def "gapseq-guessed gaps get a reduced boost"() {
        given:
        def state = new IntegrationState(new Genome(id: 'g'))
        state.metabolicGaps = [
            new MetabolicGap(
                pathwayId: 'PWY-1', reactionId: 'RXN-A',
                ecNumber: 'EC:1.1.1.1',
                goTerm: 'GO:0000001',
                gapseqGuessed: true,
            )
        ]
        def prior = new GapFillingPrior()
        prior.alphaGap = 2.0d
        prior.gapseqGuessedFactor = 0.5d

        when:
        prior.beginIteration(state)

        then:
        prior.logOddsBoost('p', 'p|EC|EC:1.1.1.1', state) == 1.0d   // 2.0 * 0.5
    }

    def "takes the max boost across duplicate gaps for the same function"() {
        given:
        def state = new IntegrationState(new Genome(id: 'g'))
        state.metabolicGaps = [
            new MetabolicGap(ecNumber: 'EC:2.7.1.1', gapseqGuessed: true),
            new MetabolicGap(ecNumber: 'EC:2.7.1.1', gapseqGuessed: false),
        ]
        def prior = new GapFillingPrior()
        prior.alphaGap = 2.0d
        prior.gapseqGuessedFactor = 0.5d

        when:
        prior.beginIteration(state)

        then:
        // Max of (2.0 gap, 1.0 gapseq-guessed) = 2.0
        prior.logOddsBoost('p', 'p|EC|EC:2.7.1.1', state) == 2.0d
    }

    def "non-EC, non-GO claims are not boosted"() {
        given:
        def state = new IntegrationState(new Genome(id: 'g'))
        state.metabolicGaps = [new MetabolicGap(ecNumber: 'EC:2.7.1.1')]
        def prior = new GapFillingPrior()

        when:
        prior.beginIteration(state)

        then:
        prior.logOddsBoost('p', 'p|PFAM|PF00001', state) == 0.0d
    }

    private static double logOdds(double p) {
        Math.log(p / (1.0d - p))
    }

    private IntegrationState buildContextState(List<String> candProteins, Map<String, Double> candProbs, boolean withAnchor = true) {
        def rg = new ReactionGraph()
        rg.addReaction(new ReactionGraph.ReactionSpec(
            rxnId: 'r0', ecNumber: '2.7.1.1',
            substrates: ['m0'] as Set, products: ['m1'] as Set))
        rg.addReaction(new ReactionGraph.ReactionSpec(
            rxnId: 'r1', ecNumber: '2.7.1.2',
            substrates: ['m1'] as Set, products: ['m2'] as Set))
        rg.bindEc('r0', '2.7.1.1')
        rg.bindEc('r1', '2.7.1.2')
        rg.build()

        def layout = new GenomeLayout('g')
        layout.add(new ProteinLocus(genomeId: 'g', proteinId: 'anchor', contig: 'chr1', start: 1000, end: 1900, strand: Strand.PLUS))
        layout.add(new ProteinLocus(genomeId: 'g', proteinId: 'pNear', contig: 'chr1', start: 2500, end: 3400, strand: Strand.PLUS))
        layout.add(new ProteinLocus(genomeId: 'g', proteinId: 'pFar', contig: 'chr1', start: 60000, end: 60900, strand: Strand.PLUS))
        layout.finishLoading()

        def state = new IntegrationState(new Genome(id: 'g'))
        state.reactionGraph = rg
        state.genomeLayout = layout
        state.ecToGoAlternatives = [
            '2.7.1.1': ['GO:0004396'],
            '2.7.1.2': ['GO:0004397'],
        ]
        state.metabolicGaps = [
            new MetabolicGap(pathwayId: 'PWY-1', reactionId: 'r0', ecNumber: 'EC:2.7.1.1', goTerm: 'GO:0004396')
        ]
        Map<String, Double> likes = new LinkedHashMap<>()
        Map<String, Integer> counts = new LinkedHashMap<>()
        if (withAnchor) {
            likes['anchor|GO|GO:0004397'] = logOdds(0.9d)
            counts['anchor|GO|GO:0004397'] = 2
        }
        for (String p : candProteins) {
            String key = "${p}|EC|EC:2.7.1.1".toString()
            likes[key] = logOdds(candProbs[p])
            counts[key] = 2
        }
        state.likelihoodLogOdds = new LinkedHashMap<>(likes)
        state.posteriorLogOdds = new LinkedHashMap<>(likes)
        state.claimCounts = counts
        state
    }

    def "(g) with contextSupport set but returning an empty map, boosts equal the no-support boosts for the same state"() {
        given:
        def stateNoSupport = buildContextState(['pNear', 'pFar'], [pNear: 0.41d, pFar: 0.44d], false)
        def stateEmptySupport = buildContextState(['pNear', 'pFar'], [pNear: 0.41d, pFar: 0.44d], false)
        def priorNoSupport = new GapFillingPrior()
        def priorEmptySupport = new GapFillingPrior(contextSupport: new GapContextSupport(mode: ContextMode.UNIVERSAL, windowBp: 10000))

        when:
        priorNoSupport.beginIteration(stateNoSupport)
        priorEmptySupport.beginIteration(stateEmptySupport)

        then:
        priorEmptySupport.contextSupport.lastReport.gaps_with_context == 0
        for (String key : ['pNear|EC|EC:2.7.1.1', 'pFar|EC|EC:2.7.1.1']) {
            String pid = IntegrationState.splitFunctionKey(key)[0]
            assert priorEmptySupport.logOddsBoost(pid, key, stateEmptySupport) ==
                priorNoSupport.logOddsBoost(pid, key, stateNoSupport)
        }
        priorEmptySupport.logOddsBoost('pFar', 'pFar|EC|EC:2.7.1.1', stateEmptySupport) > 0.0d
    }

    def "(h) with two candidates in the window at equal tier where the lower-likelihood one has c = 1 and the other c = 0, the c = 1 candidate receives max(1.15*mag, 0.14 - lo) * 1.5"() {
        given: 'both candidates in tier 3 ([0.40, 0.45)), pNear has lower likelihood (0.41) and c = 1, pFar has higher likelihood (0.44) and c = 0'
        def state = buildContextState(['pNear', 'pFar'], [pNear: 0.41d, pFar: 0.44d], true)
        def prior = new GapFillingPrior(
            alphaGap: 1.5d,
            contextSupport: new GapContextSupport(mode: ContextMode.UNIVERSAL, windowBp: 10000)
        )
        double loNear = state.getLikelihoodLogOdds('pNear|EC|EC:2.7.1.1')
        double expected = Math.max(1.15d * 1.5d, 0.14d - loNear) * 1.5d

        when:
        prior.beginIteration(state)

        then:
        Math.abs(prior.logOddsBoost('pNear', 'pNear|EC|EC:2.7.1.1', state) - expected) < 1e-12
        prior.logOddsBoost('pFar', 'pFar|EC|EC:2.7.1.1', state) == 0.0d
    }

    def "(i) a candidate with c = 0 that is the only candidate receives max(1.15*mag, 0.14 - lo) * 0.5"() {
        given: 'pFar is the only candidate in the window and is 50 kb from the anchor (c = 0)'
        def state = buildContextState(['pFar'], [pFar: 0.42d], true)
        def prior = new GapFillingPrior(
            alphaGap: 1.5d,
            contextSupport: new GapContextSupport(mode: ContextMode.UNIVERSAL, windowBp: 10000)
        )
        double loFar = state.getLikelihoodLogOdds('pFar|EC|EC:2.7.1.1')
        double expected = Math.max(1.15d * 1.5d, 0.14d - loFar) * 0.5d

        when:
        prior.beginIteration(state)

        then:
        Math.abs(prior.logOddsBoost('pFar', 'pFar|EC|EC:2.7.1.1', state) - expected) < 1e-12
    }

    def "(j) EC namespace repair: a prefixed EC candidate is reached by the targeted rescue only when the fix is on"() {
        given: 'integrated keys carry EC: while the gap record carries a bare EC number (production layout)'
        def gap = new MetabolicGap(pathwayId: 'PWY-9', reactionId: 'RXN-9', ecNumber: '3.2.2.3', goTerm: null)
        def mk = {
            def st = new IntegrationState(new Genome(id: 'g'))
            st.metabolicGaps = [gap]
            double lo = Math.log(0.42d / 0.58d)
            st.likelihoodLogOdds = new LinkedHashMap<String, Double>(['p1|EC|EC:3.2.2.3': lo])
            st.posteriorLogOdds = new LinkedHashMap<String, Double>(['p1|EC|EC:3.2.2.3': lo])
            st.claimCounts = new LinkedHashMap<String, Integer>(['p1|EC|EC:3.2.2.3': 2])
            st
        }
        def stateOff = mk()
        def stateOn = mk()
        def priorOff = new GapFillingPrior(alphaGap: 1.5d)
        def priorOn = new GapFillingPrior(alphaGap: 1.5d, ecNamespaceFix: true)

        when:
        priorOff.beginIteration(stateOff)
        priorOn.beginIteration(stateOn)

        then: 'today: no boost (the mismatch); with the fix: the targeted rescue boost'
        priorOff.logOddsBoost('p1', 'p1|EC|EC:3.2.2.3', stateOff) == 0.0d
        Math.abs(priorOn.logOddsBoost('p1', 'p1|EC|EC:3.2.2.3', stateOn) -
            Math.max(1.15d * 1.5d, 0.14d - stateOn.getLikelihoodLogOdds('p1|EC|EC:3.2.2.3'))) < 1e-12
    }
}
