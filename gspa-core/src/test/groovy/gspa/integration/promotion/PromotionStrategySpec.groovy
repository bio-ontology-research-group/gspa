package gspa.integration.promotion

import gspa.integration.IntegrationState
import gspa.integration.suggester.PerProteinDecomposition
import gspa.integration.suggester.SingletonSuggestion
import gspa.model.AnnotationType
import gspa.model.Genome
import spock.lang.Specification

class PromotionStrategySpec extends Specification {

    /** Build a minimal SingletonSuggestion with a known log-posterior. */
    private static SingletonSuggestion ss(String proteinId, String pathwayId, String rxn,
                                          String go, double logPost, double q = 0.8) {
        def dec = new PerProteinDecomposition(
            proteinId: proteinId,
            likelihoodLogOdds: logPost - 1.0d,   // arbitrary decomposition
            operonLogOdds: 0.5d,
            lmLogOdds: 0.0d,
            totalLogOdds: logPost,
            piR: 0.7d,
            q: q,
        )
        new SingletonSuggestion(
            proteinId: proteinId,
            pathwayId: pathwayId,
            reactionId: rxn,
            functionId: go,
            functionType: AnnotationType.GO,
            q: q,
            suggestionScore: 0.6d,
            bayesFactor: 50.0d,
            proteinScores: [(proteinId): dec],
        )
    }

    private IntegrationState emptyState() {
        new IntegrationState(new Genome(id: 'g'))
    }

    def "AllAboveThresholdStrategy returns every input unchanged"() {
        given:
        def strat = new AllAboveThresholdStrategy()
        def candidates = [
            ss('p1', 'PWY1', 'R1', 'GO:1', 3.0),
            ss('p2', 'PWY1', 'R1', 'GO:1', 1.0),
            ss('p3', 'PWY2', 'R2', 'GO:2', 2.0),
        ]

        when:
        def out = strat.select(candidates, emptyState(), 1)

        then:
        out.size() == 3
    }

    def "GreedyStrategy keeps best of competing suggestions for the same gap"() {
        given:
        def strat = new GreedyStrategy()
        // Three competitors on the same (PWY1, R1, GO:1) gap with descending scores.
        def candidates = [
            ss('p2', 'PWY1', 'R1', 'GO:1', 1.0),
            ss('p1', 'PWY1', 'R1', 'GO:1', 3.0),    // winner
            ss('p3', 'PWY1', 'R1', 'GO:1', 2.0),
            ss('p4', 'PWY2', 'R2', 'GO:2', 2.5),    // independent gap
        ]

        when:
        def out = strat.select(candidates, emptyState(), 1)

        then:
        out.size() == 2
        out.any { it.proteinId == 'p1' && it.pathwayId == 'PWY1' }
        out.any { it.proteinId == 'p4' && it.pathwayId == 'PWY2' }
        !out.any { it.proteinId == 'p2' }
        !out.any { it.proteinId == 'p3' }
    }

    def "GreedyStrategy respects at-most-one per protein across gaps"() {
        given:
        def strat = new GreedyStrategy()
        // Same protein wins two different gaps — the lower-score one should drop.
        def candidates = [
            ss('p1', 'PWY1', 'R1', 'GO:1', 3.0),
            ss('p1', 'PWY2', 'R2', 'GO:2', 2.0),    // same protein as above → conflict
            ss('p2', 'PWY3', 'R3', 'GO:3', 1.0),
        ]

        when:
        def out = strat.select(candidates, emptyState(), 1)

        then:
        out.size() == 2
        out.any { it.proteinId == 'p1' && it.pathwayId == 'PWY1' }
        out.any { it.proteinId == 'p2' }
        !out.any { it.pathwayId == 'PWY2' }
    }

    def "GreedyStrategy with maxPerIteration=1 degrades to strict best-first"() {
        given:
        def strat = new GreedyStrategy(maxPerIteration: 1)
        def candidates = [
            ss('p1', 'PWY1', 'R1', 'GO:1', 3.0),
            ss('p2', 'PWY2', 'R2', 'GO:2', 2.0),
            ss('p3', 'PWY3', 'R3', 'GO:3', 1.0),
        ]

        when:
        def out = strat.select(candidates, emptyState(), 1)

        then:
        out.size() == 1
        out[0].proteinId == 'p1'      // highest log-posterior
    }

    def "MaxSatStrategy picks the highest-weighted non-conflicting set"() {
        given:
        def strat = new MaxSatStrategy()
        def candidates = [
            ss('p1', 'PWY1', 'R1', 'GO:1', 3.0),
            ss('p2', 'PWY1', 'R1', 'GO:1', 1.0),    // same gap as p1 → at most one
            ss('p3', 'PWY2', 'R2', 'GO:2', 2.0),
            ss('p4', 'PWY3', 'R3', 'GO:3', 0.5),
        ]

        when:
        def out = strat.select(candidates, emptyState(), 1)

        then:
        out.size() == 3
        out.any { it.proteinId == 'p1' }
        out.any { it.proteinId == 'p3' }
        out.any { it.proteinId == 'p4' }
        !out.any { it.proteinId == 'p2' }   // lost the gap to p1
    }

    def "MaxSatStrategy respects at-most-one per protein"() {
        given:
        def strat = new MaxSatStrategy(maxPerProtein: 1)
        def candidates = [
            ss('p1', 'PWY1', 'R1', 'GO:1', 3.0),
            ss('p1', 'PWY2', 'R2', 'GO:2', 2.0),    // same protein → at most 1
            ss('p2', 'PWY3', 'R3', 'GO:3', 0.5),
        ]

        when:
        def out = strat.select(candidates, emptyState(), 1)

        then:
        out.size() == 2
        out.find { it.pathwayId == 'PWY1' }?.proteinId == 'p1'
        out.any { it.proteinId == 'p2' }
        !out.any { it.pathwayId == 'PWY2' }
    }

    def "empty candidates produce empty output regardless of strategy"() {
        expect:
        new AllAboveThresholdStrategy().select([], emptyState(), 1).isEmpty()
        new GreedyStrategy().select([], emptyState(), 1).isEmpty()
        new MaxSatStrategy().select([], emptyState(), 1).isEmpty()
        new BeamSearchStrategy().select([], emptyState(), 1).isEmpty()
    }

    def "BeamSearch picks top candidate per gap when no cross-gap conflicts"() {
        given:
        def strat = new BeamSearchStrategy(beamWidth: 5, candidatesPerGap: 3)
        def candidates = [
            ss('p1', 'PWY1', 'R1', 'GO:1', 3.0),
            ss('p2', 'PWY1', 'R1', 'GO:1', 1.0),
            ss('p3', 'PWY2', 'R2', 'GO:2', 2.5),
            ss('p4', 'PWY2', 'R2', 'GO:2', 2.0),
        ]

        when:
        def out = strat.select(candidates, emptyState(), 1)

        then:
        out.size() == 2
        out.any { it.proteinId == 'p1' && it.pathwayId == 'PWY1' }
        out.any { it.proteinId == 'p3' && it.pathwayId == 'PWY2' }
    }

    def "BeamSearch with width > 1 can beat greedy when top-1 blocks a better downstream"() {
        given: "gap G1 top-1 is p1 (score 3.0); gap G2 top-1 is ALSO p1 (score 10.0); gap G2 top-2 is p2 (score 1.0)"
        // Greedy processes G2 first (higher top score 10.0), commits p1 → G2.
        // Then for G1, p1 is taken, so falls back to top-2 or skips → suboptimal.
        // Beam with width ≥ 2 can explore "give p1 to G1, p3 to G2" variant.
        def candidates = [
            ss('p1', 'PWY1', 'R1', 'GO:1', 3.0),        // G1 top-1: p1
            ss('p2', 'PWY1', 'R1', 'GO:1', 0.5),        // G1 top-2: p2
            ss('p1', 'PWY2', 'R2', 'GO:2', 10.0),       // G2 top-1: p1 (same as G1!)
            ss('p3', 'PWY2', 'R2', 'GO:2', 9.0),        // G2 top-2: p3
        ]

        when:
        def greedyOut = new GreedyStrategy().select(candidates, emptyState(), 1)
        def beam2Out = new BeamSearchStrategy(beamWidth: 5, candidatesPerGap: 3).select(candidates, emptyState(), 1)

        then: "beam finds the joint-max assignment (p1→G1 3.0, p3→G2 9.0 = 12.0) strictly better than"
        and:  "greedy's (p1→G2 10.0, p2→G1 0.5 = 10.5)"
        def greedySum = greedyOut.collect { it.proteinScores[it.proteinId].totalLogOdds }.sum() ?: 0.0
        def beamSum = beam2Out.collect { it.proteinScores[it.proteinId].totalLogOdds }.sum() ?: 0.0
        beamSum > greedySum
    }

    def "BeamSearch degrades to single-path when width = 1"() {
        given:
        def strat = new BeamSearchStrategy(beamWidth: 1, candidatesPerGap: 3)
        def candidates = [
            ss('p1', 'PWY1', 'R1', 'GO:1', 3.0),
            ss('p2', 'PWY1', 'R1', 'GO:1', 2.0),
            ss('p3', 'PWY2', 'R2', 'GO:2', 2.5),
        ]

        when:
        def out = strat.select(candidates, emptyState(), 1)

        then:
        out.size() == 2
        out.any { it.proteinId == 'p1' }
        out.any { it.proteinId == 'p3' }
    }
}
