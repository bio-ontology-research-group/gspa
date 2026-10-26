package gspa.integration.cohort

import spock.lang.Specification

import static gspa.integration.cohort.CohortFixtures.*

/**
 * Proteins that really exist but have no usable coordinate (missing, split, multi-locus) are declared
 * unplaced: the family member is PRESENT, its context is UNKNOWN, and nothing is invented or called absent.
 * SYNTHETIC genomes; these check mechanisms, not accuracy.
 */
class CohortUnplacedSpec extends Specification {

    private CohortInput target() {
        addGenome(newInput(), [id: 'G0', clade: 'cT', cluster: 'sT', status: 'COMPLETE', closed: true,
                               family: 'OG1', candidate: 'T', anchors: ['N1', 'N2']])
    }

    private CohortInput withPeers(Map peer, int n = 3) {
        def in0 = target()
        n.times { addGenome(in0, [id: "P${it}", clade: "c${it}", cluster: "s${it}", status: 'COMPLETE', closed: true,
                                  family: 'OG1', anchors: ['N1', 'N2']] + peer) }
        in0
    }

    def 'a family member that is unplaced counts as present with unknown context, not as absent or negative'() {
        when:
        def f = features(withPeers([unplaced: [20]]))

        then: 'every peer holds the family: no absence of any kind, no negative, context unknown'
        f.cohort_frac_unknown_context == 1d
        f.cohort_frac_true_absent == 0d
        f.cohort_frac_unobserved_absent == 0d
        f.cohort_frac_negative == 0d
        f.cohort_frac_supported == 0d
    }

    def 'complete absence is a different state from an unobserved context'() {
        when:
        def absent = features(withPeers([family: null, anchors: []]))
        def unplaced = features(withPeers([unplaced: [20]]))

        then:
        absent.cohort_frac_true_absent == 1d
        absent.cohort_frac_unknown_context == 0d
        unplaced.cohort_frac_true_absent == 0d
        unplaced.cohort_frac_unknown_context == 1d
    }

    def 'one placed copy without anchors plus one unplaced paralog is not called negative'() {
        when: 'the placed copy sits in a fully observed window without anchors'
        def placedOnly = features(withPeers([family: 'OG1', anchors: [], famIdx: 20]))
        def withUnplaced = features(withPeers([family: 'OG1', anchors: [], famIdx: 20, unplaced: [30], paralogs: [30]]))

        then: 'alone it is negative evidence'
        placedOnly.cohort_frac_negative == 1d

        and: 'the unresolved paralog may be the supported copy, so the genome stays unknown'
        withUnplaced.cohort_frac_negative == 0d
        withUnplaced.cohort_frac_unknown_context == 1d
        withUnplaced.cohort_frac_true_absent == 0d
    }

    def 'one placed supported copy plus one unplaced paralog stays supported'() {
        when:
        def f = features(withPeers([unplaced: [30], paralogs: [30]]))

        then:
        f.cohort_frac_supported == 1d
        f.cohort_frac_negative == 0d
    }

    def 'an unplaced focal protein has no local context and no anchor, but keeps its family and raw candidate'() {
        given:
        def in0 = addGenome(newInput(), [id: 'G0', clade: 'cT', cluster: 'sT', status: 'COMPLETE', closed: true,
                                         family: 'OG1', candidate: 'T', anchors: ['N1', 'N2'], unplaced: [20]])
        3.times { addGenome(in0, [id: "P${it}", clade: "c${it}", cluster: "s${it}", status: 'COMPLETE', closed: true,
                                  family: 'OG1', anchors: ['N1', 'N2']]) }

        when:
        def f = features(in0)

        then: 'local context is missing, never an invented zero-anchor window'
        f.local_layout_missing == 1d
        !f.containsKey('local_anchor_count')
        !f.containsKey('local_context_positive')
        and: 'family-level recurrence and the raw candidate are retained'
        f.family_known == 1d
        f.baseline_confidence == 0.3d
        f.cohort_frac_supported == 1d
    }

    def 'ids that are neither placed nor declared are still rejected; declared ids are accepted'() {
        given:
        def in0 = target()
        in0.addClaim(new CohortInput.CandidateClaim(genomeId: 'G0', proteinId: 'G0_typo', reactionId: 'U1',
                baselineConfidence: 0.9d, provenance: 'seq@snap'))

        when:
        in0.validate()

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('G0_typo')
        e.message.contains('not declared unplaced')

        when: 'the same id declared unplaced'
        in0.genomes.G0.unplacedProteins.add('G0_typo')
        in0.validate()

        then:
        noExceptionThrown()
    }

    def 'a family assignment or role evidence on an undeclared missing protein is rejected too'() {
        given:
        def in0 = target()
        in0.assign('G0', new CohortInput.FamilyAssignment(proteinId: 'G0_nope', orthogroupId: 'OG1'))

        when:
        in0.validate()

        then:
        thrown(IllegalArgumentException)
    }

    def 'declaring a placed protein unplaced, or a blank id, is a contradiction'() {
        given:
        def in0 = target()
        in0.genomes.G0.unplacedProteins.add(id)

        when:
        in0.validate()

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains(fragment)

        where:
        id      | fragment
        'G0_5'  | 'has a coordinate'
        ' '     | 'blank'
    }
}
