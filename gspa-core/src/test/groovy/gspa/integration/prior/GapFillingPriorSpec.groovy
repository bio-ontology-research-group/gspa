package gspa.integration.prior

import gspa.integration.IntegrationState
import gspa.integration.MetabolicGap
import gspa.model.AnnotationType
import gspa.model.Genome
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
}
