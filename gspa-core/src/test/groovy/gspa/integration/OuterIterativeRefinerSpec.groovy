package gspa.integration

import gspa.integration.suggester.DarkMatterSuggester
import gspa.integration.suggester.SingletonSuggestion
import gspa.integration.suggester.Suggestion
import gspa.model.AnnotationType
import gspa.model.Genome
import spock.lang.Specification

/**
 * Unit tests for OuterIterativeRefiner. Uses a stubbed inner refiner
 * (returns a fixed IntegratedAnnotationSet) and a stubbed suggester
 * (scripted per-iteration promotions) to isolate the outer-loop
 * convergence / cascade-rollback / pinning logic from Phase 7's
 * full likelihood math.
 */
class OuterIterativeRefinerSpec extends Specification {

    /** Inner refiner stub: returns an empty integrated set. */
    private static IterativeRefiner trivialInner() {
        def combiner = new EvidenceCombiner()
        def r = new IterativeRefiner(combiner) {
            @Override
            IntegratedAnnotationSet refine(List<EvidenceClaim> claims, IntegrationState state) {
                new IntegratedAnnotationSet()
            }
        }
        r
    }

    /**
     * Suggester stub: on each suggest() call, pulls the next element
     * of scriptedSuggestions and writes it into integrated.suggestions.
     */
    private static DarkMatterSuggester scriptedSuggester(List<List<Suggestion>> script) {
        int[] idx = [0] as int[]
        new DarkMatterSuggester() {
            @Override
            IntegratedAnnotationSet suggest(IntegrationState state, IntegratedAnnotationSet integrated) {
                int i = idx[0]++
                integrated.suggestions = (i < script.size()) ? script[i] : []
                integrated
            }
        }
    }

    private static SingletonSuggestion singleton(String p, String pw, String rxn, String go, double q = 0.9) {
        new SingletonSuggestion(
            proteinId: p,
            pathwayId: pw,
            reactionId: rxn,
            functionId: go,
            functionType: AnnotationType.GO,
            q: q,
            suggestionScore: 0.7d,
            bayesFactor: 50.0d,
        )
    }

    private IntegrationState emptyState() {
        def s = new IntegrationState(new Genome(id: 'g'))
        s.metabolicGaps = []
        s.operons = []
        s
    }

    def "reaches fixed point immediately when suggester emits nothing"() {
        given:
        def outer = new OuterIterativeRefiner(trivialInner())
        outer.suggester = scriptedSuggester([])
        // No-op gap source so currentGaps stays empty.
        outer.gapSource = [recompute: { st, cs -> new LinkedHashSet<GapKey>() }] as OuterIterativeRefiner.GapSource

        when:
        def result = outer.refine([], emptyState())

        then:
        result.fixedPointReached
        result.outerIterationsRun == 0
        result.promotedPerIter.size() == 1   // iter 1 ran the empty-suggestion check
        result.promotedPerIter[0] == 0
    }

    def "promotes strong singleton, pins its floor, marks gap closed, terminates"() {
        given:
        def outer = new OuterIterativeRefiner(trivialInner())
        outer.suggester = scriptedSuggester([
            [singleton('p1', 'PWY-1', 'R1', 'GO:0001', 0.9d)],    // iter 0: emit 1 singleton
            [],                                                   // iter 1: nothing new → converge
        ])
        outer.gapSource = [recompute: { st, cs -> new LinkedHashSet<GapKey>() }] as OuterIterativeRefiner.GapSource
        def state = emptyState()

        when:
        def result = outer.refine([], state)

        then:
        result.fixedPointReached
        result.outerIterationsRun == 1
        result.promotedPerIter == [1, 0]
        // Pin floor was set for the promoted claim.
        def ck = new ClaimKey(proteinId: 'p1', functionType: AnnotationType.GO, functionId: 'GO:0001')
        state.getPinnedFloor(ck) != null
        state.getPinnedFloor(ck) > 0.0d        // positive log-odds for q=0.9
        // Gap is marked closed.
        state.isGapClosed(new GapKey(pathwayId: 'PWY-1', reactionId: 'R1', goTerm: null))
    }

    def "respects rising q threshold — weak singleton not promoted at late iter"() {
        given: "iteration 2's threshold is qBase + qStep*2 = 0.5 + 0.1 = 0.6; a 0.55 suggestion is rejected"
        def outer = new OuterIterativeRefiner(trivialInner())
        outer.qBase = 0.5d
        outer.qStep = 0.05d
        // k=1: threshold 0.55; q=0.6 passes. k=2: threshold 0.60; q=0.55 fails.
        outer.suggester = scriptedSuggester([
            [singleton('p1', 'PWY-1', 'R1', 'GO:0001', 0.6d)],
            [singleton('p2', 'PWY-1', 'R2', 'GO:0002', 0.55d)],
            [],
        ])
        outer.gapSource = [recompute: { st, cs -> new LinkedHashSet<GapKey>() }] as OuterIterativeRefiner.GapSource
        def state = emptyState()

        when:
        def result = outer.refine([], state)

        then:
        result.promotedPerIter[0] == 1  // iter 1: p1 promoted
        result.promotedPerIter[1] == 0  // iter 2: p2 below rising threshold
        state.isGapClosed(new GapKey(pathwayId: 'PWY-1', reactionId: 'R1', goTerm: null))
        !state.isGapClosed(new GapKey(pathwayId: 'PWY-1', reactionId: 'R2', goTerm: null))
    }

    def "rolls back on cascade (two consecutive rising iterations)"() {
        given: "script: 1 promo, 2 promo, 3 promo — triggers cascade guard after iter 3"
        def outer = new OuterIterativeRefiner(trivialInner())
        outer.maxIter = 10
        outer.suggester = scriptedSuggester([
            [singleton('p1', 'P1', 'R1', 'GO:1', 0.9d)],
            [singleton('p2', 'P2', 'R2', 'GO:2', 0.9d),
             singleton('p3', 'P2', 'R3', 'GO:3', 0.9d)],
            [singleton('p4', 'P3', 'R4', 'GO:4', 0.9d),
             singleton('p5', 'P3', 'R5', 'GO:5', 0.9d),
             singleton('p6', 'P3', 'R6', 'GO:6', 0.9d)],
        ])
        outer.gapSource = [recompute: { st, cs -> new LinkedHashSet<GapKey>() }] as OuterIterativeRefiner.GapSource
        def state = emptyState()

        when:
        def result = outer.refine([], state)

        then:
        result.cascadeRolledBack
        !result.fixedPointReached
    }

    def "pinPromotions=false leaves pinnedFloors empty even when gaps are closed"() {
        given:
        def outer = new OuterIterativeRefiner(trivialInner())
        outer.pinPromotions = false
        outer.suggester = scriptedSuggester([
            [singleton('p1', 'PWY-1', 'R1', 'GO:0001', 0.9d)],
            [],
        ])
        outer.gapSource = [recompute: { st, cs -> new LinkedHashSet<GapKey>() }] as OuterIterativeRefiner.GapSource
        def state = emptyState()

        when:
        outer.refine([], state)

        then:
        state.pinnedFloors.isEmpty()
        state.isGapClosed(new GapKey(pathwayId: 'PWY-1', reactionId: 'R1', goTerm: null))
    }
}
