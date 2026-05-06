package gspa.integration.prior

import gspa.integration.IntegrationState
import gspa.model.Genome
import spock.lang.Specification

class HomologyTransferPriorSpec extends Specification {

    def "no-op when orthogroupMap or consensus is null"() {
        given:
        def state = new IntegrationState(new Genome(id: 'test'))
        def prior = new HomologyTransferPrior()

        expect:
        prior.logOddsBoost('p1', 'p1|GO|GO:0006412', state) == 0.0d
    }

    def "no-op when protein has no orthogroup assignment"() {
        given:
        def state = new IntegrationState(new Genome(id: 'test'))
        state.orthogroupMap = ['p1': 'cl1']
        state.orthogroupConsensus = ['cl1|GO|GO:0006412': 0.9d]
        def prior = new HomologyTransferPrior()

        expect:
        prior.logOddsBoost('p99_unmapped', 'p99_unmapped|GO|GO:0006412', state) == 0.0d
    }

    def "no-op when cluster has no consensus for this function"() {
        given:
        def state = new IntegrationState(new Genome(id: 'test'))
        state.orthogroupMap = ['p1': 'cl1']
        state.orthogroupConsensus = ['cl1|GO|GO:0001111': 0.9d]   // different function
        def prior = new HomologyTransferPrior()

        expect:
        prior.logOddsBoost('p1', 'p1|GO|GO:0006412', state) == 0.0d
    }

    def "no-op when cluster consensus is below minConsensus"() {
        given:
        def state = new IntegrationState(new Genome(id: 'test'))
        state.orthogroupMap = ['p1': 'cl1']
        state.orthogroupConsensus = ['cl1|GO|GO:0006412': 0.3d]   // weak consensus
        def prior = new HomologyTransferPrior()

        expect:
        prior.logOddsBoost('p1', 'p1|GO|GO:0006412', state) == 0.0d
    }

    def "boosts weak claim toward strong cluster consensus"() {
        given: "cluster consensus 0.9 and current posterior ~0.2"
        def state = new IntegrationState(new Genome(id: 'test'))
        state.orthogroupMap = ['p1': 'cl1']
        state.orthogroupConsensus = ['cl1|GO|GO:0006412': 0.9d]
        // seed current posterior at log-odds ~-1.39 (p=0.2)
        state.set('p1|GO|GO:0006412', -1.39d)

        def prior = new HomologyTransferPrior()
        prior.alpha = 1.0d

        when:
        double boost = prior.logOddsBoost('p1', 'p1|GO|GO:0006412', state)

        then: "log-odds pushed toward logit(0.9)=2.197 by alpha times the diff"
        boost > 0.0d
        boost <= prior.maxBoost
        // boost = alpha * (logit(0.9) - current) = 1.0 * (2.197 - (-1.39)) = 3.587 → clipped to maxBoost=3.0
        Math.abs(boost - 3.0d) < 1e-6
    }

    def "zero boost when current claim is already as strong as consensus"() {
        given:
        def state = new IntegrationState(new Genome(id: 'test'))
        state.orthogroupMap = ['p1': 'cl1']
        state.orthogroupConsensus = ['cl1|GO|GO:0006412': 0.7d]
        state.set('p1|GO|GO:0006412', 1.5d)   // current prob ~0.82 > consensus 0.7

        def prior = new HomologyTransferPrior()

        expect:
        prior.logOddsBoost('p1', 'p1|GO|GO:0006412', state) == 0.0d
    }

    def "respects minDelta — doesn't fire for tiny improvements"() {
        given: "consensus only marginally above current"
        def state = new IntegrationState(new Genome(id: 'test'))
        state.orthogroupMap = ['p1': 'cl1']
        state.orthogroupConsensus = ['cl1|GO|GO:0006412': 0.52d]
        state.set('p1|GO|GO:0006412', Math.log(0.50d / 0.50d))  // logit(0.5) = 0; prob=0.5

        def prior = new HomologyTransferPrior()
        prior.minDelta = 0.05d   // require 5-point prob gap

        expect:
        prior.logOddsBoost('p1', 'p1|GO|GO:0006412', state) == 0.0d
    }
}
