package gspa.integration.cohort

import spock.lang.Specification

import static gspa.integration.cohort.CohortFixtures.*

/** Regressions for the pre-merge review of the cohort package. */
class CohortReviewSpec extends Specification {

    private CohortInput targetOnly(Map extra = [:]) {
        addGenome(newInput(), [id: 'G0', clade: 'cT', cluster: 'sT', status: 'COMPLETE', closed: true,
                               family: 'OG1', candidate: 'T', anchors: ['N1', 'N2']] + extra)
    }

    private CohortInput withUnknownPeers(int n, Map peer = [:]) {
        def in0 = targetOnly()
        n.times { addGenome(in0, [id: "U${it}", status: 'COMPLETE', closed: true, family: 'OG1', anchors: ['N1', 'N2']] + peer) }
        in0
    }

    // --- 2. unknown phylogeny --------------------------------------------------

    def "duplicated genomes with unknown metadata cannot increase the independent effective N"() {
        when:
        def one = features(withUnknownPeers(1))
        def many = features(withUnknownPeers(12))

        then: 'all unknown clades form ONE unresolved block, whatever its size'
        one.cohort_effective_n == 1d
        many.cohort_effective_n == 1d
        many.cohort_n_peer_clades == 1d
        many.cohort_informative_clades == 1d
        many.cohort_raw_genomes_supported == 12d
        and: 'adjusted statistics are withheld; descriptive counts remain'
        many.cohort_adjusted_available == 0d
        !many.containsKey('cohort_support_fraction_lo90')
        !many.containsKey('cohort_support_fraction_hi90')
        !many.containsKey('cohort_perm_p')
        many.cohort_phylogeny_known_fraction == 0d
        many.cohort_support_fraction_raw == 1d
    }

    def "unknown clades together with a resolved clade count the unresolved block once"() {
        given:
        def in0 = withUnknownPeers(5)
        addGenome(in0, [id: 'K0', clade: 'cK', cluster: 'sK', status: 'COMPLETE', closed: true, family: 'OG1', anchors: ['N1', 'N2']])

        when:
        def f = features(in0)

        then:
        f.cohort_effective_n == 2d
        f.cohort_n_peer_clades == 2d
        f.cohort_adjusted_available == 1d
        f.containsKey('cohort_support_fraction_lo90')
    }

    def "clade or cluster holdout of an unknown-metadata target does not leak its unknown relatives"() {
        given:
        def in0 = newInput()
        addGenome(in0, [id: 'G0', status: 'COMPLETE', closed: true, family: 'OG1', candidate: 'T', anchors: ['N1']])
        3.times { addGenome(in0, [id: "U${it}", status: 'COMPLETE', closed: true, family: 'OG1', anchors: ['N1']]) }
        addGenome(in0, [id: 'K0', clade: 'cK', cluster: 'sK', status: 'COMPLETE', closed: true, family: 'OG1', anchors: ['N1']])

        when:
        def clade = features(in0) { it.holdout = HoldoutLevel.CLADE }
        def cluster = features(in0) { it.holdout = HoldoutLevel.STRAIN_CLUSTER }
        def genome = features(in0) { it.holdout = HoldoutLevel.GENOME }

        then: 'the whole unresolved block goes with the target; only the resolved genome remains'
        clade.cohort_n_peer_genomes == 1d
        cluster.cohort_n_peer_genomes == 1d
        genome.cohort_n_peer_genomes == 4d
    }

    def "a cluster listed under several clades is rejected, a cluster without clade inherits its clade"() {
        when:
        def bad = newInput()
        addGenome(bad, [id: 'A', clade: 'c1', cluster: 's', family: null])
        addGenome(bad, [id: 'B', clade: 'c2', cluster: 's', family: null])
        new CohortInference(bad, graph(), new CohortConfig())

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('several clades')

        when:
        def ok = newInput()
        addGenome(ok, [id: 'A', clade: 'c1', cluster: 's', family: null])
        addGenome(ok, [id: 'B', cluster: 's', family: null])
        def model = new CohortInference(ok, graph(), new CohortConfig()).cladeModel()

        then:
        model.clade['B'] == 'c1'
        model.cladeKnown('B')
        model.units(['A', 'B'])*.genomes == [['A', 'B']]
    }

    // --- 3. search certification and unannotated labelling --------------------

    def "a COMPLETE assembly without a certified family search gives unobserved, not true, absence"() {
        given:
        def uncertified = targetOnly()
        addGenome(uncertified, [id: 'A', clade: 'c1', cluster: 's1', status: 'COMPLETE', closed: true, search: false, family: null])
        def certified = targetOnly()
        addGenome(certified, [id: 'A', clade: 'c1', cluster: 's1', status: 'COMPLETE', closed: true, search: true, family: null])
        def draftCertified = targetOnly()
        addGenome(draftCertified, [id: 'A', clade: 'c1', cluster: 's1', status: 'DRAFT', closed: false, search: true, family: null])

        when:
        def u = features(uncertified)
        def c = features(certified)
        def d = features(draftCertified)

        then:
        u.cohort_frac_unobserved_absent == 1d
        u.cohort_frac_true_absent == 0d
        c.cohort_frac_true_absent == 1d
        c.cohort_frac_unobserved_absent == 0d
        d.cohort_frac_unobserved_absent == 1d
    }

    def "copy number is unknown without a certified search"() {
        given:
        def in0 = targetOnly()
        addGenome(in0, [id: 'A', clade: 'c1', cluster: 's1', status: 'COMPLETE', closed: true, search: false, family: 'OG1', anchors: ['N1']])
        addGenome(in0, [id: 'B', clade: 'c2', cluster: 's2', status: 'COMPLETE', closed: true, search: true, family: 'OG1', anchors: ['N1']])

        when:
        def f = features(in0)

        then: 'only the certified genome counts as a known single copy'
        f.cohort_single_copy_effective_n == 1d
        f.cohort_multicopy_fraction == 0d
    }

    def "contradictory search certification is rejected"() {
        given:
        def in0 = targetOnly()
        in0.genomes['G0'].familySearchCoverage = 0.8d

        when:
        new CohortInference(in0, graph(), new CohortConfig())

        then:
        thrown(IllegalArgumentException)
    }

    // --- 6. input validation ---------------------------------------------------

    def "malformed membership, ids, evidence and duplicates are rejected before computing"() {
        when:
        def in0 = mutate.call(targetOnly())
        new CohortInference(in0, graph(), new CohortConfig())

        then:
        thrown(IllegalArgumentException)

        where:
        label                | mutate
        'NaN weight'         | { CohortInput i -> i.assign('G0', new CohortInput.FamilyAssignment(proteinId: 'G0_3', orthogroupId: 'X', weight: Double.NaN)); i }
        'negative weight'    | { CohortInput i -> i.assign('G0', new CohortInput.FamilyAssignment(proteinId: 'G0_3', orthogroupId: 'X', weight: -0.1d)); i }
        'weight above one'   | { CohortInput i -> i.assign('G0', new CohortInput.FamilyAssignment(proteinId: 'G0_3', orthogroupId: 'X', weight: 1.5d)); i }
        'sum above one'      | { CohortInput i -> i.assign('G0', new CohortInput.FamilyAssignment(proteinId: 'G0_20', orthogroupId: 'Y', weight: 0.4d)); i }
        'unknown protein'    | { CohortInput i -> i.assign('G0', new CohortInput.FamilyAssignment(proteinId: 'nope', orthogroupId: 'X')); i }
        'unknown genome'     | { CohortInput i -> i.assign('ZZ', new CohortInput.FamilyAssignment(proteinId: 'G0_3', orthogroupId: 'X')); i }
        'claim protein'      | { CohortInput i -> i.addClaim(new CohortInput.CandidateClaim(genomeId: 'G0', proteinId: 'nope', reactionId: 'N1', baselineConfidence: 0.5d)); i }
        'claim NaN'          | { CohortInput i -> i.addClaim(new CohortInput.CandidateClaim(genomeId: 'G0', proteinId: 'G0_3', reactionId: 'Z', baselineConfidence: Double.NaN)); i }
        'duplicate claim'    | { CohortInput i -> i.addClaim(new CohortInput.CandidateClaim(genomeId: 'G0', proteinId: 'G0_20', reactionId: 'T', baselineConfidence: 0.5d)); i }
        'role out of range'  | { CohortInput i -> i.addRoleEvidence(new CohortInput.RoleEvidence(genomeId: 'G0', proteinId: 'G0_3', scores: [(Role.ACCESSORY): 1.2d] as Map)); i }
        'role NaN'           | { CohortInput i -> i.addRoleEvidence(new CohortInput.RoleEvidence(genomeId: 'G0', proteinId: 'G0_3', scores: [(Role.REGULATOR): Double.NaN] as Map)); i }
        'role protein'       | { CohortInput i -> i.addRoleEvidence(new CohortInput.RoleEvidence(genomeId: 'G0', proteinId: 'nope', scores: [(Role.REGULATOR): 0.5d] as Map)); i }
        'organism genome'    | { CohortInput i -> i.setOrganismReactions('ZZ', ['T']); i }
        'half-open edge'     | { CohortInput i -> i.addEdge(new CohortInput.InteractionEdge(reactionA: 'T')); i }
    }

    def "duplicate genomes and layout genome mismatches are rejected"() {
        when:
        def in0 = targetOnly()
        addGenome(in0, [id: 'G0', family: null])

        then:
        thrown(IllegalArgumentException)

        when:
        def in1 = targetOnly()
        in1.genomes['G0'].layout.genomeId = 'OTHER'
        new CohortInference(in1, graph(), new CohortConfig())

        then:
        thrown(IllegalArgumentException)
    }

    def "malformed trees are rejected instead of falling back to independence"() {
        given:
        def ids = ['A', 'B', 'C'] as Set

        when:
        CladeModel.cutTree(newick, 0.05d, ids)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains(fragment)

        where:
        newick                              | fragment
        '((A:0.1,B:0.1):0.1,A:0.1);'        | 'Duplicate'
        '((A:0.1,B:0.1):0.1,D:0.1);'        | 'not a genome'
        '(A:0.1,B:0.1);'                    | 'missing'
        '((A:0.1,B:-0.1):0.1,C:0.1);'       | 'non-negative'
        '((A:0.1,B:NaN):0.1,C:0.1);'        | 'non-negative'
        '((A:0.1,B:0.1):0.1,C:0.1); junk'   | 'Trailing'
        '((A:0.1,B:0.1):0.1,C:0.1)'         | 'end with'
        '((A:0.1,B:0.1):0.1,C);'            | 'no branch length'
        '((A:0.1,B:0.1):0.1,:0.1);'         | 'unnamed'
        '((A:0.1,B:0.1):0.1,C:0.1;'         | 'Malformed'
    }

    def "a tree that does not cover the cohort genomes is rejected by the engine"() {
        given:
        def in0 = targetOnly()
        in0.treeNewick = '(G0:0.1,X:0.1);'

        when:
        new CohortInference(in0, graph(), new CohortConfig())

        then:
        thrown(IllegalArgumentException)
    }

    // --- 7. scale --------------------------------------------------------------

    def "dense candidate sets are not capped: every input claim yields a record"() {
        given:
        def in0 = newInput()
        12.times { i ->
            addGenome(in0, [id: "S${i}", clade: "c${i % 4}", cluster: "s${i}", status: 'COMPLETE', closed: true,
                            nGenes: 300, famIdx: 150, family: 'OG1', anchors: ['N1', 'N2'], candidate: 'T', candConf: 0.05d])
            // extra weak EC-style candidates on every gene
            (0..299).each { k ->
                if (k != 150) in0.addClaim(new CohortInput.CandidateClaim(genomeId: "S${i}", proteinId: "S${i}_${k}",
                        reactionId: 'LONE', baselineConfidence: 0.01d, provenance: 'seq@snap'))
            }
        }

        when:
        def res = new CohortInference(in0, graph(), new CohortConfig(permutations: 20)).run()

        then:
        res.records.size() == in0.claims.size()
        res.record('S0', 'S0_150', 'T').baselineConfidence == 0.05d
        res.records.every { r -> r.features.values().every { !it.isNaN() } }
    }
}
