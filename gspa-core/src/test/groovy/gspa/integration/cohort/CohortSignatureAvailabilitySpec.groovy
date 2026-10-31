package gspa.integration.cohort

import spock.lang.Specification

import static gspa.integration.cohort.CohortFixtures.*

/** An untestable reaction signature cannot supply negative recurrence evidence. */
class CohortSignatureAvailabilitySpec extends Specification {
    private CohortInput cohort(String reaction = 'LONE', Map peer = [:]) {
        def input = newInput()
        addGenome(input, [id: 'G0', clade: 'target', cluster: 'target', status: 'COMPLETE',
                          closed: true, family: 'OG1', candidate: reaction])
        addGenome(input, [id: 'P0', clade: 'peer', cluster: 'peer', status: 'COMPLETE',
                          closed: true, family: 'OG1'] + peer)
        input
    }

    def 'a present family without a testable signature is unknown rather than negative'() {
        given:
        def input = cohort()
        if (edge) input.addEdge(new CohortInput.InteractionEdge(reactionA: 'LONE', reactionB: 'N1', source: 'independent'))

        when:
        def result = new CohortInference(input, graph(), new CohortConfig(minAnchors: minimum)).run()
        def f = result.record('G0', 'G0_20', 'LONE').features

        then:
        f.cohort_raw_genomes_unknown_context == 1d
        f.cohort_raw_genomes_negative == 0d
        f.cohort_effective_n == 0d
        f.cohort_unknown_eff == 1d
        f.cohort_adjusted_available == 0d
        f.cohort_frac_unknown_context == 1d
        f.cohort_frac_negative == 0d
        !f.containsKey('cohort_perm_p')
        !f.containsKey('cohort_enrichment_log')
        f.signature_testable == 0d
        !f.containsKey('local_context_positive')
        !f.containsKey('local_enrichment_log')

        where:
        edge  | minimum
        false | 1
        true  | 2
    }

    def 'independent interaction evidence makes an isolated reaction testable'() {
        given:
        def input = cohort('LONE', [anchors: anchors])
        input.addEdge(new CohortInput.InteractionEdge(reactionA: 'LONE', reactionB: 'N1', source: 'independent'))

        when:
        def f = features(input, null, 'G0', 20, 'LONE')

        then:
        f.cohort_effective_n == 1d
        f.cohort_raw_genomes_supported == supported
        f.cohort_raw_genomes_negative == 1d - supported
        f.cohort_raw_genomes_unknown_context == 0d

        where:
        anchors | supported
        ['N1']  | 1d
        []      | 0d
    }

    def 'an untestable signature preserves family presence and certified copy counts'() {
        when:
        def f = features(cohort('LONE', [paralogs: [30], weight: 0.7d]), null, 'G0', 20, 'LONE')

        then:
        f.cohort_unknown_eff == 0.7d
        f.cohort_multicopy_fraction == 1d
        f.cohort_raw_genomes_true_absent == 0d
        f.cohort_raw_genomes_unobserved_absent == 0d
        f.cohort_support_eff == 0d
        f.cohort_negative_eff == 0d
    }

    def 'family absence is still distinguished from untestable context'() {
        when:
        def f = features(cohort('LONE', [family: null, search: certified]), null, 'G0', 20, 'LONE')

        then:
        f.cohort_raw_genomes_true_absent == (certified ? 1d : 0d)
        f.cohort_raw_genomes_unobserved_absent == (certified ? 0d : 1d)
        f.cohort_raw_genomes_unknown_context == 0d
        f.cohort_effective_n == 0d

        where:
        certified << [true, false]
    }
}
